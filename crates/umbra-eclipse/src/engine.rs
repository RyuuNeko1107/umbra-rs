//! 日食エンジンの外殻（`docs/api-draft.md` §3.2・ISSUE-043 S5b, 確定A1）。
//!
//! [`EclipseEngine`]`<E, D, O>` は暦 [`Ephemeris`]・ΔT [`DeltaTModel`]・地球姿勢
//! [`EarthOrientation`] にジェネリックなまま日食計算を駆動する外殻（**`Box<dyn>` 不使用**,
//! 確定A1）。本スライス（S5b）では構築 [`EclipseEngine::new`]・標準構築
//! [`standard_engine`]／型エイリアス [`StandardEngine`]・瞬時要素
//! [`EclipseEngine::instantaneous_elements`]・未実装 [`EclipseEngine::path`]
//! （`Err(NotImplemented)`）のみを提供する。`search`/`local_circumstances`/
//! `next_visible_eclipse` は後続スライス（S6/S7/S8）で追加する。
//!
//! 瞬時要素は供給源 [`InstantaneousEvaluator`](crate::source::InstantaneousEvaluator)
//! （ISSUE-043 S3・暦ジェネリックな直接評価）を退化区間 `[time, time]` で構築して
//! `.at(time)` を呼ぶ。`AnalyticalEphemeris` ＋ `AstrometryOptions::standard()` では
//! `besselian_elements_at`（ISSUE-021）と同等の瞬時要素になる（S2 回帰ブリッジ済）。

use std::time::{SystemTime, UNIX_EPOCH};

use umbra_core::constants::EARTH_EQUATORIAL_RADIUS_M;
use umbra_core::deltat::{decimal_year, DeltaTModel};
use umbra_core::ellipsoid::{observer_geocentric, Ellipsoid, GeocentricObserver};
use umbra_core::eop::EarthOrientation;
use umbra_core::{
    EspenakMeeusDeltaT, IersEopData, JulianDate2, Observer, Radians, SolverError, TimeData,
    TimeInterval, TimeRange, TimeScales, TtInstant, UtcInstant,
};
use umbra_ephemeris::{AnalyticalEphemeris, AstrometryOptions, Ephemeris};
use umbra_geo::{GeoLine, GeoPoint};

use crate::axis_intercept::surface_point_for_fundamental;
use crate::bessel_poly::{BesselFitError, BesselianPolynomial};
use crate::besselian::InstantaneousBesselianElements;
use crate::calc_metadata::CalculationMetadata;
use crate::candidates::new_moon_candidates;
use crate::config::EngineConfig;
use crate::conjunction::{solve_conjunction, ConjunctionKind, RootConfig};
use crate::eclipse_filter::assess_eclipse_possibility;
use crate::error::EclipseError;
use crate::global::{classify_global_kind, solve_greatest_eclipse};
use crate::global_contacts::solve_global_contact_set;
use crate::horizontal::{classify_visibility, sun_horizontal, RefractionModel, Visibility};
use crate::local_contacts::{solve_local_contacts, ContactInstant};
use crate::local_maximum::solve_local_maximum;
use crate::magnitude::{eclipse_magnitude, eclipse_obscuration, EclipseMagnitude, Obscuration};
use crate::path::{EclipsePath, PathOptions};
use crate::position_angle::contact_position_angle;
use crate::projection::{project_observer_to_fundamental, ObserverFundamental};
use crate::results::{
    GlobalCircumstances, LocalCircumstances, LocalContact, LocalContactSet, SolarEclipse,
    VisibleSolarEclipse,
};
use crate::source::{BesselianSource, InstantaneousEvaluator};

use std::collections::HashSet;

/// 1 日 = 86400 SI 秒。
const SECONDS_PER_DAY: f64 = 86_400.0;
/// 探索段の Brent 反復上限。
const SEARCH_ROOT_MAX_ITER: usize = 200;
/// ベッセル多項式 fit の開始次数（NASA 慣習 3。残差未達なら内部で自動昇次, ≤6）。
const BESSEL_FIT_START_DEGREE: usize = 3;
/// ベッセル多項式 fit の残差ゲート（Re。x/y/l1/l2 の最大残差上限。v0.1 標準・要精緻化）。
const BESSEL_FIT_TOLERANCE: f64 = 1.0e-4;
/// UNIX エポックの UTC ユリウス日（generated_at の壁時計変換に使用）。
const UNIX_EPOCH_JD: f64 = 2_440_587.5;
/// `next_visible_eclipse` の探索窓幅（日）。≈半年＝日食季間隔オーダー。窓ごとに `search`（重い全球解）
/// を回し最初の可視で打ち切るため、1 窓に含む日食を ~1 件に抑えて遅延評価する（無駄な全球解を減らす）。
const NEXT_VISIBLE_WINDOW_DAYS: f64 = 183.0;
/// `next_visible_eclipse` の探索 horizon（日, ≈11 年）。これを超えて可視日食が無ければ `Ok(None)`。
/// 任意地点は通常数年内に部分食を見られるため実地点では到達しない安全弁（accuracy/可視性の §7）。
const NEXT_VISIBLE_MAX_HORIZON_DAYS: f64 = 4000.0;

/// 探索範囲（UTC）。
pub type UtcRange = TimeRange<UtcInstant>;

/// 日食エンジン（ジェネリック E/D/O 維持・`Box<dyn>` 不使用, 確定A1）。
///
/// `earth_orientation`（極運動・UT1）と `time_scales`（UTC↔TT↔UT1 facade）は後続スライス
/// （S6 search の UTC 範囲変換・S7 局地の観測者 ITRS/極運動）が消費する。本 S5b では
/// `instantaneous_elements`（ephemeris + delta_t）と `path`（未実装）のみが使う。
#[derive(Debug)]
#[allow(dead_code)] // earth_orientation / time_scales は S6/S7 で消費（結線され次第この許容を外す）。
pub struct EclipseEngine<E: Ephemeris, D: DeltaTModel, O: EarthOrientation> {
    /// 天体暦バックエンド。
    ephemeris: E,
    /// ΔT モデル（μ の UT1 変換）。
    delta_t: D,
    /// 地球姿勢（UT1−UTC・極運動）。
    earth_orientation: O,
    /// 時刻系変換 facade（`TimeData` から構築, 確定B3）。
    time_scales: TimeScales,
    /// エンジン設定。
    config: EngineConfig,
}

impl<E: Ephemeris, D: DeltaTModel, O: EarthOrientation> EclipseEngine<E, D, O> {
    /// 引数の `TimeData` から `TimeScales::new(time)` を構築して保持する（確定B3）。
    pub fn new(
        ephemeris: E,
        delta_t: D,
        earth_orientation: O,
        time: TimeData,
        config: EngineConfig,
    ) -> Self {
        Self {
            ephemeris,
            delta_t,
            earth_orientation,
            time_scales: TimeScales::new(time),
            config,
        }
    }

    /// 1 時刻の瞬時ベッセル要素（検証/CLI inspect 用）。
    ///
    /// 供給源 [`InstantaneousEvaluator`] を半径（config の太陽/月モデル）・見かけ補正
    /// [`AstrometryOptions::standard`]（標準/参照とも全補正 ON）・退化区間 `[time, time]` で構築し
    /// `.at(time)` を評価する。μ は `delta_t`（`DeltaTModel`）由来。
    /// 注（精度・後続）: EOP coverage 内では `delta_t`（Espenak–Meeus 外挿）由来 UT1 が EOP 由来
    /// UT1 と数秒差を持ちうる（μ→局地接触）。EOP 由来 UT1 を μ に使う精緻化は S7/精度工程で扱う。
    pub fn instantaneous_elements(
        &self,
        time: TtInstant,
    ) -> Result<InstantaneousBesselianElements, EclipseError> {
        let re_km = EARTH_EQUATORIAL_RADIUS_M / 1000.0;
        let r_sun_km = self.config.solar_radius_model.radius_km();
        let r_moon_km = self.config.lunar_radius_model.k() * re_km;
        let evaluator = InstantaneousEvaluator::new(
            &self.ephemeris,
            &self.delta_t,
            r_sun_km,
            r_moon_km,
            AstrometryOptions::standard(),
            TimeInterval {
                start: time,
                end: time,
            },
        );
        evaluator.at(time)
    }

    /// UTC 範囲の日食を探索して `SolarEclipse` のリストを返す（ISSUE-043 S6c）。
    ///
    /// パイプライン: 新月候補（ISSUE-016）→ 合（ISSUE-017）→ 早期棄却（ISSUE-018, 偽陰性不可）→
    /// 残候補ごとに供給源（暦ジェネリック `InstantaneousEvaluator`）を構築し、最大食（S6a）・
    /// 全球接触（S6b-i/ii）・種別（S6b-iii）・ベッセル多項式 fit（ISSUE-022）を計算して
    /// `SolarEclipse`（event_key・GlobalCircumstances・bessel・metadata）を組み立てる。
    pub fn search(&self, range: UtcRange) -> Result<Vec<SolarEclipse>, EclipseError> {
        let root_config = RootConfig {
            x_tolerance_days: self.config.root_tolerance_seconds / SECONDS_PER_DAY,
            max_iterations: SEARCH_ROOT_MAX_ITER,
        };
        let re_km = EARTH_EQUATORIAL_RADIUS_M / 1000.0;
        let r_sun_km = self.config.solar_radius_model.radius_km();
        let r_moon_km = self.config.lunar_radius_model.k() * re_km;
        // earth_model は現状 WGS84 のみ（projection/global と同様に定数で扱う）。
        let _ = self.config.earth_model;
        let ellipsoid = Ellipsoid::WGS84;
        let generated_at = utc_now();

        let mut eclipses = Vec::new();
        for candidate in new_moon_candidates(range)? {
            // 合 → 早期棄却（偽陰性不可・偽陽性可）。棄却された朔は日食でない。
            let conjunction =
                solve_conjunction(&candidate, ConjunctionKind::EclipticLongitude, root_config)?;
            if !assess_eclipse_possibility(&conjunction).possible {
                continue;
            }

            // 候補窓の供給源（暦ジェネリック直接評価, ISSUE-037）。全球 solver を直接源で駆動する。
            let source = InstantaneousEvaluator::new(
                &self.ephemeris,
                &self.delta_t,
                r_sun_km,
                r_moon_km,
                AstrometryOptions::standard(),
                candidate.search_window,
            );

            // 種別。`None` は全球で日食なし（早期棄却の偽陽性）→ 採用しない。
            let Some(kind) = classify_global_kind(&source, root_config)? else {
                continue;
            };
            let solution = solve_greatest_eclipse(&source, &self.delta_t, &self.config)?;
            let contacts = solve_global_contact_set(&source, &ellipsoid, root_config)?;

            // ベッセル多項式（結果の携帯表現）。fit 区間は全球部分食 [P1,P4]（無ければ候補窓）、
            // エポックは最大食 TT。次数は開始 3 から自動昇次（≤6）で残差ゲートを満たす。
            let fit_interval = match (contacts.p1, contacts.p4) {
                (Some(p1), Some(p4)) => TimeInterval {
                    start: p1.time_tt,
                    end: p4.time_tt,
                },
                _ => candidate.search_window,
            };
            let tolerance = BesselFitError {
                max_x: BESSEL_FIT_TOLERANCE,
                max_y: BESSEL_FIT_TOLERANCE,
                max_l1: BESSEL_FIT_TOLERANCE,
                max_l2: BESSEL_FIT_TOLERANCE,
            };
            let bessel = BesselianPolynomial::fit(
                &source,
                solution.greatest.time_tt,
                fit_interval,
                BESSEL_FIT_START_DEGREE,
                tolerance,
            )?;

            let global = GlobalCircumstances {
                kind,
                partial_begin: contacts.p1,
                central_begin: contacts.u1,
                greatest: solution.greatest,
                central_end: contacts.u4,
                partial_end: contacts.p4,
                gamma: solution.gamma,
            };
            let event_key = format_event_key(solution.greatest.time_utc, candidate.lunation_number);
            let metadata = self.build_metadata(solution.greatest.time_tt, generated_at);
            eclipses.push(SolarEclipse {
                event_key,
                kind,
                global,
                bessel,
                metadata,
            });
        }
        Ok(eclipses)
    }

    /// 観測地点の局地条件を計算する（ISSUE-043 S7b）。
    ///
    /// 既定 Standard は直接瞬時計算（`InstantaneousEvaluator`, ISSUE-037, fit 誤差ゼロ, B2）。
    /// パイプライン: 観測者射影（ISSUE-024）→ 局地最大食（ISSUE-026, dm/dt=0）→ 食分/食面積
    /// （ISSUE-027）→ 太陽高度方位/位置角/可視性（ISSUE-028・S7a PA）→ C1-C4 接触（ISSUE-025）→
    /// `LocalCircumstances`（接触は UTC/TT 両方, A3 LocalContactSet）。
    ///
    /// 接触は C1-C4＋最大食（A3 LocalContactSet）。部分食地点は内接 C2/C3 が `None`、食域外は
    /// 全接触 `None`＋`Visibility::NotVisible`。可視性は実 C1/C4 高度で 6 値判定する（ISSUE-028）。
    /// 遠方観測者（探索窓内で局地最小をブラケット不能）は全球最大食時刻に錨を打ち NotVisible・食分0。
    pub fn local_circumstances(
        &self,
        eclipse: &SolarEclipse,
        observer: Observer,
    ) -> Result<LocalCircumstances, EclipseError> {
        let root_config = RootConfig {
            x_tolerance_days: self.config.root_tolerance_seconds / SECONDS_PER_DAY,
            max_iterations: SEARCH_ROOT_MAX_ITER,
        };
        let re_km = EARTH_EQUATORIAL_RADIUS_M / 1000.0;
        let r_sun_km = self.config.solar_radius_model.radius_km();
        let r_moon_km = self.config.lunar_radius_model.k() * re_km;
        let _ = self.config.earth_model;
        let ellipsoid = Ellipsoid::WGS84;
        let generated_at = utc_now();

        // 観測者 → 地心動径成分 ρsinφ′/ρcosφ′（WGS84 扁平・標高込み, ISSUE-024）＋ 東経・測地緯度。
        let geodetic_latitude = observer.latitude.radians();
        let east_longitude = observer.longitude.radians();
        let geo_obs = observer_geocentric(&ellipsoid, geodetic_latitude.0, observer.elevation.0);

        // 探索窓: 全球部分食 [P1,P4]（無ければベッセル fit 区間）。局地接触/最大食はこの中。
        let window = match (eclipse.global.partial_begin, eclipse.global.partial_end) {
            (Some(p1), Some(p4)) => TimeInterval {
                start: p1.time_tt,
                end: p4.time_tt,
            },
            _ => eclipse.bessel.fit_interval,
        };

        // 供給源（search と同一の暦ジェネリック直接評価, ISSUE-037）。
        let source = InstantaneousEvaluator::new(
            &self.ephemeris,
            &self.delta_t,
            r_sun_km,
            r_moon_km,
            AstrometryOptions::standard(),
            window,
        );

        // 局地最大食（dm/dt=0, ISSUE-026）。窓内に内部極小をブラケットできない遠方観測者は
        // 全球最大食時刻に錨を打ち NotVisible・食分0 を返す（S7b 確定）。
        match solve_local_maximum(&source, &geo_obs, east_longitude, window, root_config) {
            Ok(max) => {
                let elements = source.at(max.time_tt)?;
                let of = project_observer_to_fundamental(&geo_obs, east_longitude, &elements);
                // ζ 補正半径 L1'=l1−ζ·tanf1 / L2'=l2−ζ·tanf2（符号付き, global solve_greatest_eclipse と同方式）。
                let l1p = elements.l1 - of.zeta * elements.tan_f1;
                let l2p = elements.l2 - of.zeta * elements.tan_f2;
                let magnitude = eclipse_magnitude(max.min_separation, l1p, l2p);
                // 視半径比 ρ=(L1'−L2')/(L1'+L2')、視半径平面の中心離隔 separation=(1+ρ)·m/L1'。
                let radius_ratio = (l1p - l2p) / (l1p + l2p);
                let separation = (1.0 + radius_ratio) * max.min_separation / l1p;
                let obscuration = eclipse_obscuration(separation, 1.0, radius_ratio);
                // 最大食 LocalContact（最大食は接触点が月中心方向 ⇒ PA は σ=+1）。
                let maximum = self.build_local_contact(
                    &elements,
                    &of,
                    max.time_utc,
                    max.time_tt,
                    geodetic_latitude,
                    east_longitude,
                    false,
                );
                // 局地接触 C1-C4（ISSUE-025）。部分食地点は内接 C2/C3 が None。各接触は
                // `build_contact_at` で時刻＋高度方位＋PA＋可視を付与（C2/C3 は内接 ⇒ 皆既で σ=−1）。
                let contacts =
                    solve_local_contacts(&source, &geo_obs, east_longitude, window, root_config)?;
                let c1 = self.contact_local(
                    &source,
                    &geo_obs,
                    geodetic_latitude,
                    east_longitude,
                    contacts.c1,
                    false,
                )?;
                let c2 = self.contact_local(
                    &source,
                    &geo_obs,
                    geodetic_latitude,
                    east_longitude,
                    contacts.c2,
                    true,
                )?;
                let c3 = self.contact_local(
                    &source,
                    &geo_obs,
                    geodetic_latitude,
                    east_longitude,
                    contacts.c3,
                    true,
                )?;
                let c4 = self.contact_local(
                    &source,
                    &geo_obs,
                    geodetic_latitude,
                    east_longitude,
                    contacts.c4,
                    false,
                )?;

                // 可視性（実 C1/C4 高度で精緻化）。in_eclipse は食分>0 で判定。
                let in_eclipse = magnitude.0 > 0.0;
                let visibility = classify_visibility(
                    in_eclipse,
                    c1.map(|c| c.sun_altitude),
                    maximum.sun_altitude,
                    c4.map(|c| c.sun_altitude),
                );
                Ok(LocalCircumstances {
                    contacts: LocalContactSet {
                        c1,
                        c2,
                        maximum,
                        c3,
                        c4,
                    },
                    magnitude,
                    obscuration,
                    maximum_altitude: maximum.sun_altitude,
                    visibility,
                    metadata: self.build_metadata(max.time_tt, generated_at),
                })
            }
            Err(EclipseError::Solver(SolverError::RootNotBracketed)) => {
                // 遠方観測者: 局地最小が窓内に無い → 非可視。全球最大食時刻に錨（食分0・接触なし）。
                let time_tt = eclipse.global.greatest.time_tt;
                let time_utc = eclipse.global.greatest.time_utc;
                let elements = source.at(time_tt)?;
                let of = project_observer_to_fundamental(&geo_obs, east_longitude, &elements);
                let maximum = self.build_local_contact(
                    &elements,
                    &of,
                    time_utc,
                    time_tt,
                    geodetic_latitude,
                    east_longitude,
                    false,
                );
                Ok(LocalCircumstances {
                    contacts: LocalContactSet {
                        c1: None,
                        c2: None,
                        maximum,
                        c3: None,
                        c4: None,
                    },
                    magnitude: EclipseMagnitude(0.0),
                    obscuration: Obscuration(0.0),
                    maximum_altitude: maximum.sun_altitude,
                    visibility: Visibility::NotVisible,
                    metadata: self.build_metadata(time_tt, generated_at),
                })
            }
            Err(other) => Err(other),
        }
    }

    /// 1 時点の `LocalContact` を組み立てる（時刻 ＋ 太陽高度方位 ＋ 位置角 PA ＋ 可視）。
    ///
    /// 太陽高度は幾何学的高度（大気差なし, conventions §7 既定）、方位は北0東回り（ISSUE-028）、
    /// 位置角 PA は天の北0東回り（ISSUE-043 S7a, `umbral_interior` は皆既内接 C2/C3 のみ true）、
    /// `visible` は太陽が地平上（幾何高度 ≥ 0）か。S7b-i は最大食のみ、S7b-ii で C1-C4 も使う。
    #[allow(clippy::too_many_arguments)]
    fn build_local_contact(
        &self,
        elements: &InstantaneousBesselianElements,
        observer_fundamental: &ObserverFundamental,
        time_utc: UtcInstant,
        time_tt: TtInstant,
        geodetic_latitude: Radians,
        east_longitude: Radians,
        umbral_interior: bool,
    ) -> LocalContact {
        let horizontal = sun_horizontal(
            geodetic_latitude,
            east_longitude,
            time_tt,
            RefractionModel::None,
            &self.delta_t,
        );
        let position_angle =
            contact_position_angle(elements, observer_fundamental, umbral_interior);
        LocalContact {
            time_utc,
            time_tt,
            sun_altitude: horizontal.altitude_geometric,
            sun_azimuth: horizontal.azimuth,
            position_angle,
            visible: horizontal.altitude_geometric.0 >= 0.0,
        }
    }

    /// 接触時刻 [`ContactInstant`]（存在すれば）を rich [`LocalContact`] へ昇格する（ISSUE-043 S7b-ii）。
    ///
    /// `None`（その地点に当該接触なし。部分食地点の内接 C2/C3 等）はそのまま `None`。`interior` は
    /// 内接接触（C2/C3）か外接（C1/C4）かで、内接かつ本影 `l2<0`（皆既）のときのみ PA を σ=−1
    /// （接触点が月中心の反対側, S7a）にする。供給源 `source` で接触時刻の瞬時要素・観測者射影を取る。
    fn contact_local<B: BesselianSource>(
        &self,
        source: &B,
        geo_obs: &GeocentricObserver,
        geodetic_latitude: Radians,
        east_longitude: Radians,
        contact: Option<ContactInstant>,
        interior: bool,
    ) -> Result<Option<LocalContact>, EclipseError> {
        let Some(contact) = contact else {
            return Ok(None);
        };
        let elements = source.at(contact.time_tt)?;
        let of = project_observer_to_fundamental(geo_obs, east_longitude, &elements);
        // 内接 C2/C3 かつ皆既（l2<0）のみ σ=−1（接触点が月中心の反対側）。外接・金環内接は σ=+1。
        let umbral_interior = interior && elements.l2 < 0.0;
        Ok(Some(self.build_local_contact(
            &elements,
            &of,
            contact.time_utc,
            contact.time_tt,
            geodetic_latitude,
            east_longitude,
            umbral_interior,
        )))
    }

    /// 計算メタデータ（レシピ＋生成時刻印）を組み立てる（accuracy.md §0）。暦名は `ephemeris.metadata()`、
    /// ΔT 名は `delta_t.model_name()`、ΔT 不確かさは最大食 TT の十進年で評価、月半径名は Debug
    /// （バリアント名）、地球モデルは WGS84 固定（現状単一）、ライブラリ版は crate version。
    fn build_metadata(&self, time_tt: TtInstant, generated_at: UtcInstant) -> CalculationMetadata {
        let em = self.ephemeris.metadata();
        let year = decimal_year(time_tt.jd2());
        CalculationMetadata {
            library_version: env!("CARGO_PKG_VERSION").to_string(),
            ephemeris_model: em.model,
            ephemeris_version: em.version,
            delta_t_model: self.delta_t.model_name().to_string(),
            delta_t_uncertainty_seconds: self.delta_t.uncertainty_seconds(year),
            earth_model: "WGS84".to_string(),
            lunar_radius_model: self.config.lunar_radius_model.name().to_string(),
            accuracy_profile: self.config.accuracy,
            generated_at,
        }
    }

    /// `after` 以降で観測者が最初に「見える」日食を返す（ISSUE-043 S8）。
    ///
    /// `after` 以降を窓刻みで [`search`](Self::search) 走査し、各日食の
    /// [`local_circumstances`](Self::local_circumstances) を評価して、可視性が「見える」種別
    /// （[`next_visible_is_observable`]）になる最初を `Some(VisibleSolarEclipse)` で返す。探索 horizon
    /// 内に見える日食が無ければ `Ok(None)`（「該当なし」はエラーにしない, api-draft §0）。
    ///
    /// 注（性能）: `search`/`local_circumstances` は直接瞬時計算（ISSUE-037）で日食 1 件あたり重い。
    /// 窓刻みで遅延評価し最初の可視で打ち切るが、可視日食が遠い/無い場合は horizon まで走査する。
    pub fn next_visible_eclipse(
        &self,
        after: UtcInstant,
        observer: Observer,
    ) -> Result<Option<VisibleSolarEclipse>, EclipseError> {
        let after_jd = after.jd2().jd();
        let horizon_end_jd = after_jd + NEXT_VISIBLE_MAX_HORIZON_DAYS;
        // event_key 重複排除（窓境界・候補オーバーハングで同一日食が隣接窓に現れるのを 1 回に）。
        let mut seen: HashSet<String> = HashSet::new();
        let mut start_jd = after_jd;
        while start_jd < horizon_end_jd {
            let end_jd = (start_jd + NEXT_VISIBLE_WINDOW_DAYS).min(horizon_end_jd);
            let range = UtcRange {
                start: UtcInstant::from_jd2(JulianDate2::from_jd(start_jd)),
                end: UtcInstant::from_jd2(JulianDate2::from_jd(end_jd)),
            };
            // search は新月候補昇順ゆえ日食も昇順。最初に「見える」ものが時系列最初の可視日食。
            for eclipse in self.search(range)? {
                // `after` より前（first 窓の取りこぼし）・既評価（境界重複）はスキップ。
                if eclipse.global.greatest.time_utc.jd2().jd() < after_jd {
                    continue;
                }
                if !seen.insert(eclipse.event_key.clone()) {
                    continue;
                }
                let local = self.local_circumstances(&eclipse, observer)?;
                if next_visible_is_observable(local.visibility) {
                    return Ok(Some(VisibleSolarEclipse { eclipse, local }));
                }
            }
            start_jd = end_jd;
        }
        Ok(None)
    }

    /// 日食経路を生成する（M9 第1スライス＝中心線トラック）。
    ///
    /// **中心食**（全球 U1/U4 接触＝`central_begin`/`central_end` が両方 `Some`）では、`[U1.time_tt,
    /// U4.time_tt]` を `options.sample_interval_seconds` 刻みでサンプルし、各時刻のベッセル要素
    /// （`bessel.at`）から影軸の地表貫通点（[`shadow_axis_surface_point`]）を結んだ中心線 [`GeoLine`] を
    /// `center_line` に返す。軸が地球を外す端の時刻（`RootNotBracketed`）はスキップする。**非中心**
    /// （部分/非中心食）では `center_line = None`。`greatest_point` は常に `global.greatest.position`。
    ///
    /// 本スライスでは北/南限界線・部分食域・経路サンプル（`samples`）は未生成（`None` / 空）。これらは
    /// 帯幅・中心食継続式（算法 §8.11/8.12・要一次資料確認）や限界線追跡を要するため後続スライス（M9）で
    /// 実装する。`bessel.at` / 影軸貫通の `RootNotBracketed` 以外の `Err` は伝播する。
    pub fn path(
        &self,
        eclipse: &SolarEclipse,
        options: PathOptions,
    ) -> Result<EclipsePath, EclipseError> {
        let greatest_point = eclipse.global.greatest.position;
        // 中心食（U1/U4 両方 Some）でのみ中心線・限界線を追跡。片方でも None なら経路なし。
        let (center_line, northern_limit, southern_limit) =
            match (&eclipse.global.central_begin, &eclipse.global.central_end) {
                (Some(u1), Some(u4)) => {
                    let (center, limits) = trace_central(
                        &eclipse.bessel,
                        u1.time_tt,
                        u4.time_tt,
                        options.sample_interval_seconds,
                        options.include_limits,
                    )?;
                    match limits {
                        Some((north, south)) => (Some(center), Some(north), Some(south)),
                        None => (Some(center), None, None),
                    }
                }
                _ => (None, None, None),
            };
        Ok(EclipsePath {
            center_line,
            northern_limit,
            southern_limit,
            partial_limit: None,
            greatest_point,
            samples: Vec::new(),
        })
    }
}

/// 可視性が「見える」種別か（`next_visible_eclipse` の採否判定, ISSUE-043 S8）。
///
/// 地平上で日食を観測できる `FullyVisible`/`PartialVisible`/`SunriseEclipse`/`SunsetEclipse` を
/// `true`、観測不能な `NotVisible`（食域外）/`BelowHorizon`（最大食も地平下）を `false` とする。
///
/// 注: `Visibility` は `#[non_exhaustive]`。将来バリアントが追加されると `matches!` の暗黙の既定で
/// `false`（不可視＝安全側）になる。新たに「見える」種別を足す場合は本関数のアーム追加が必要。
fn next_visible_is_observable(visibility: Visibility) -> bool {
    matches!(
        visibility,
        Visibility::FullyVisible
            | Visibility::PartialVisible
            | Visibility::SunriseEclipse
            | Visibility::SunsetEclipse
    )
}

/// 中心食の中心線（と任意で南北限界線）を `[start_tt, end_tt]` を `interval_seconds` 刻みでサンプルして
/// 追跡する（M9.1 中心線 / M9.3 限界線・[`EclipseEngine::path`] から呼ぶ）。
///
/// 各サンプル時刻で `bessel.at` の瞬時要素から影軸∩WGS84 地表点（中心線）を求める。`include_limits` の
/// とき、ζ補正本影半径 `L2' = l2 − ζ·tan f2` を影の運動方向 (x',y') に垂直へ ±|L2'| オフセットした 2 点も
/// 地表へ射影し、**高緯度側を北限・低緯度側を南限**とする（geometric 近似・厳密な錐接線解は後続。
/// accuracy.md / ISSUE-045 に明記）。中心線と南北限界線が**同じサンプル列**になるよう、軸または縁が
/// 地表を外す（`RootNotBracketed`）/ 影速度ゼロのサンプルは**3 本とも**スキップする（lockstep）。
/// `RootNotBracketed` 以外の `Err` は伝播。始点・終点を必ず含む（端は span にクランプ）。`interval_seconds`
/// 非正は始点のみ（無限ループ回避）。前提 `start_tt ≤ end_tt`（U1≤U4・逆順は始点のみの無害な縮退）。
#[allow(clippy::type_complexity)]
fn trace_central(
    bessel: &BesselianPolynomial,
    start_tt: TtInstant,
    end_tt: TtInstant,
    interval_seconds: f64,
    include_limits: bool,
) -> Result<(GeoLine, Option<(GeoLine, GeoLine)>), EclipseError> {
    let ellipsoid = Ellipsoid::WGS84;
    // 影速度 (x', y')（基本面・hour 単位。垂直方向の比のみ使うので単位は相殺）。
    let x_deriv = bessel.x.derivative();
    let y_deriv = bessel.y.derivative();
    let epoch = bessel.epoch_tt;
    let span_seconds = end_tt.jd2().days_since(start_tt.jd2()) * SECONDS_PER_DAY;

    let mut center_points = Vec::new();
    let mut north_points = Vec::new();
    let mut south_points = Vec::new();
    let mut t_sec = 0.0_f64;
    loop {
        let t = TtInstant::from_jd2(start_tt.jd2().add_days(t_sec / SECONDS_PER_DAY));
        if let Some((center, limits)) = sample_central_point(
            bessel,
            &x_deriv,
            &y_deriv,
            epoch,
            t,
            include_limits,
            &ellipsoid,
        )? {
            center_points.push(center);
            if let Some((north, south)) = limits {
                north_points.push(north);
                south_points.push(south);
            }
        }
        if t_sec >= span_seconds || interval_seconds <= 0.0 {
            break;
        }
        t_sec = (t_sec + interval_seconds).min(span_seconds);
    }

    let limits = include_limits.then(|| (GeoLine::new(north_points), GeoLine::new(south_points)));
    Ok((GeoLine::new(center_points), limits))
}

/// 1 サンプル時刻の中心線点（と任意で南北限界点）を求める。軸/縁が地表を外す（`RootNotBracketed`）/
/// 影速度ゼロなら `Ok(None)`（lockstep スキップ）。他の `Err` は伝播。
#[allow(clippy::type_complexity)]
fn sample_central_point(
    bessel: &BesselianPolynomial,
    x_deriv: &crate::polynomial::Polynomial,
    y_deriv: &crate::polynomial::Polynomial,
    epoch: TtInstant,
    t: TtInstant,
    include_limits: bool,
    ellipsoid: &Ellipsoid,
) -> Result<Option<(GeoPoint, Option<(GeoPoint, GeoPoint)>)>, EclipseError> {
    let elements = bessel.at(t)?;
    let (center, zeta0) = match surface_point_for_fundamental(
        elements.x,
        elements.y,
        elements.declination,
        elements.mu,
        ellipsoid,
    ) {
        Ok(v) => v,
        Err(EclipseError::Solver(SolverError::RootNotBracketed)) => return Ok(None),
        Err(e) => return Err(e),
    };
    if !include_limits {
        return Ok(Some((center, None)));
    }

    // ζ補正本影半径 |L2'| = |l2 − ζ·tan f2|（本影 l2<0 なら |l2|+ζ·tan f2、反本影 l2>0 も同式）。
    let umbral_radius = (elements.l2 - zeta0 * elements.tan_f2).abs();
    // 影の運動方向 (x', y') に垂直な単位ベクトル n = (−y', x')/|v|（基本面）。
    let t_hours = t.jd2().days_since(epoch.jd2()) * 24.0;
    let vx = x_deriv.eval(t_hours);
    let vy = y_deriv.eval(t_hours);
    let speed = vx.hypot(vy);
    if speed == 0.0 {
        return Ok(None); // 影速度ゼロ（退行）はスキップ（lockstep 維持）。
    }
    let (nx, ny) = (-vy / speed, vx / speed);
    // 帯の北縁/南縁＝軸から ±|L2'| だけ垂直オフセットした基本面点を地表へ射影。
    let edge = |sign: f64| -> Result<Option<GeoPoint>, EclipseError> {
        let xi = elements.x + sign * umbral_radius * nx;
        let eta = elements.y + sign * umbral_radius * ny;
        match surface_point_for_fundamental(xi, eta, elements.declination, elements.mu, ellipsoid) {
            Ok((p, _)) => Ok(Some(p)),
            Err(EclipseError::Solver(SolverError::RootNotBracketed)) => Ok(None),
            Err(e) => Err(e),
        }
    };
    let (Some(edge_a), Some(edge_b)) = (edge(1.0)?, edge(-1.0)?) else {
        return Ok(None); // どちらかの縁が地表を外すサンプルは 3 本ともスキップ。
    };
    // 高緯度側＝北限・低緯度側＝南限（geometric 近似・経路がほぼ東西走行の前提, ISSUE-045）。
    let (north, south) = if edge_a.lat.degrees().0 >= edge_b.lat.degrees().0 {
        (edge_a, edge_b)
    } else {
        (edge_b, edge_a)
    };
    Ok(Some((center, Some((north, south)))))
}

/// 現在の壁時計 UTC を `UtcInstant` で返す（`generated_at` 用）。std 時計 → UNIX 秒 → UTC-JD。
/// fingerprint からは除外される時刻印なので再現性に影響しない（確定: std 時計で生成時に刻む）。
/// 時計取得失敗時は UNIX エポックにフォールバック（生成時刻印のみゆえ無害）。
fn utc_now() -> UtcInstant {
    let secs = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_secs_f64())
        .unwrap_or(0.0);
    UtcInstant::from_jd2(JulianDate2::from_jd(UNIX_EPOCH_JD + secs / SECONDS_PER_DAY))
}

/// `event_key` = `"YYYY-MM-DD#K"`（最大食の UTC 暦日 ＋ `#` ＋ lunation 番号, ISSUE-043）。
/// lunation は **符号付き 10 進**（Meeus lunation index, `k=0`↔2000-01-06 朔。1900–2100 では
/// 約 −1240…+1240）。負（2000 年朔より前）も素直に `-K` で表す（固定幅ゼロ詰めは符号で崩れるため
/// 使わない。日付は固定幅 YYYY-MM-DD）。
fn format_event_key(greatest_utc: UtcInstant, lunation_number: i64) -> String {
    let (year, month, day, _, _, _) = greatest_utc.to_gregorian();
    format!("{year:04}-{month:02}-{day:02}#{lunation_number}")
}

/// 標準エンジン型エイリアス（確定A1, dyn 不使用）。
pub type StandardEngine = EclipseEngine<AnalyticalEphemeris, EspenakMeeusDeltaT, IersEopData>;

/// 同梱データで標準エンジンを構築（`EngineConfig::standard()`）。`O` は `time` の EOP を複製して
/// 保持する（`time.eop()` と整合）。例: `standard_engine(bundled_time_data())`。
pub fn standard_engine(time: TimeData) -> StandardEngine {
    let earth_orientation = time.eop().clone();
    EclipseEngine::new(
        AnalyticalEphemeris::new(),
        EspenakMeeusDeltaT,
        earth_orientation,
        time,
        EngineConfig::standard(),
    )
}

#[cfg(test)]
mod tests {
    //! ISSUE-043 S5b 受け入れテスト（strict・EclipseEngine 外殻）。
    //!
    //! ## オラクル戦略（実装方針に立ち入らず、確定仕様の公開 IF だけを縛る）
    //! - **瞬時要素 = 既知の独立オラクル**: `instantaneous_elements(2017 最大食 TT)` の `gamma()` が
    //!   NASA 公表 gamma≈0.4367 と 4 桁一致する（既存 `besselian.rs`/`source.rs` の実日食ゲートと
    //!   同値域 `[0.43, 0.44]`）。これは時刻系/EOP を使わず ephemeris+delta_t のみで決まるため、
    //!   合成 EOP/TimeData でも 2017 gamma が出る（確定仕様）。
    //! - **path() = `Err(NotImplemented)`**: 引数を使わず常にこの variant を返す（`matches!`）。
    //! - **StandardEngine / standard_engine がコンパイル・動作**（受け入れ §78, `Box<dyn>` 不使用＝
    //!   E/D/O の単相化がコンパイルできることで担保）。
    //!
    //! ## red 設計（本体未実装）
    //! `EclipseEngine::new`/`instantaneous_elements`/`path`/`standard_engine`/`StandardEngine`/
    //! `UtcRange` は本体未実装（外殻のみ・各メソッドは `unimplemented!`）。テストは未存在の
    //! 振る舞い（戻り値・gamma 値）を要求するため、`unimplemented!` の panic／未存在シンボルの
    //! コンパイルエラーで red になる。実装は本体側で追加する。

    use super::*;

    use umbra_core::constants::EARTH_EQUATORIAL_RADIUS_M;
    use umbra_core::{
        DataSetMetadata, EopRecord, EspenakMeeusDeltaT, IersEopData, JulianDate2, LeapSecondTable,
        TimeData, TtInstant, UtcInstant,
    };
    use umbra_ephemeris::AnalyticalEphemeris;

    use crate::config::{AccuracyProfile, EngineConfig, LunarRadiusModel, SolarRadiusModel};
    use crate::global::SolarEclipseKind;
    use crate::horizontal::Visibility;
    use crate::magnitude::{EclipseMagnitude, Obscuration};
    use crate::path::PathOptions;
    use crate::results::{GlobalCircumstances, GreatestEclipse, SolarEclipse};

    // ------------------------------------------------------------------
    // 時刻ヘルパ（既存 source.rs/besselian.rs と同一エポック）
    // ------------------------------------------------------------------

    /// TT 時刻を 2 要素 JD から構築するヘルパ。
    fn tt(jd1: f64, jd2: f64) -> TtInstant {
        TtInstant::from_jd2(JulianDate2::new(jd1, jd2))
    }

    /// 2017-08-21 最大食付近の TT（besselian.rs/source.rs テストと同一エポック）。
    /// NASA 公表 gamma≈0.4367 の独立オラクルが効くエポック。
    fn tt_2017_max() -> TtInstant {
        tt(2_457_986.5, 7.685_322_222_222_222e-1)
    }

    /// J2000.0（TT）。別エポックでのサニティ用。
    fn tt_j2000() -> TtInstant {
        tt(2_451_545.0, 0.0)
    }

    // ------------------------------------------------------------------
    // 合成 TimeData / IersEopData（eop.rs/timescales.rs テストの 2020 レコード流用）
    //
    // instantaneous_elements は時刻系/EOP を使わず ephemeris+delta_t のみ使うため、
    // 合成 EOP でも 2017 gamma が出る（確定仕様）。`standard_engine` の §78 経路のみ
    // 同梱データ（feature-gated）を使う。
    // ------------------------------------------------------------------

    /// provenance 完全な代表 EOP metadata（全フィールド非空, timescales.rs テスト流用）。
    fn eop_metadata() -> DataSetMetadata {
        DataSetMetadata {
            name: "iers-eop-c04".to_string(),
            version: "EOP 14 C04".to_string(),
            source: "IERS Earth Orientation Center, datacenter.iers.org".to_string(),
            license: "public-domain".to_string(),
            valid_from: "2020-01-01".to_string(),
            valid_to: "2020-01-02".to_string(),
            checksum: "0000000000000000000000000000000000000000000000000000000000000000"
                .to_string(),
        }
    }

    /// 2 点 {58849, 58850}（2020-01-01/02）のみの合成 EOP（eop.rs/timescales.rs と同一値）。
    fn synthetic_eop() -> IersEopData {
        IersEopData::from_records(
            vec![
                EopRecord::new(58849, -0.177_122_2, 0.076_609, 0.282_358),
                EopRecord::new(58850, -0.177_580_6, 0.074_635, 0.282_666),
            ],
            "EOP 14 C04".to_string(),
            eop_metadata(),
        )
        .expect("two adjacent ascending 2020 records build")
    }

    /// 同梱閏秒 + 合成 EOP の合成 TimeData。
    fn synthetic_time_data() -> TimeData {
        TimeData::new(LeapSecondTable::bundled(), synthetic_eop())
    }

    /// 合成データから StandardEngine を構築する（feature 非依存・pipeline テスト用）。
    /// `standard_engine` は §78 で同梱データ経由を別途検証するため、ここでは `new` を直接使う。
    fn standard_engine_from_synthetic() -> StandardEngine {
        let time = synthetic_time_data();
        // O は time の EOP の複製（standard_engine と同じ整合・確定仕様）。
        let eop = time.eop().clone();
        EclipseEngine::new(
            AnalyticalEphemeris::new(),
            EspenakMeeusDeltaT,
            eop,
            time,
            EngineConfig::standard(),
        )
    }

    // ------------------------------------------------------------------
    // 最小 SolarEclipse（path() テスト用・results.rs テストの構築パターン）
    //
    // path は引数を使わず常に Err(NotImplemented) を返すので最小値でよい。
    // BesselianPolynomial は results.rs と同様 pub フィールドで最小構築する。
    // ------------------------------------------------------------------

    /// path() テスト用の UTC 瞬時。
    fn utc(year: i32, month: u8, day: u8, hour: u8, minute: u8, second: f64) -> UtcInstant {
        UtcInstant::from_gregorian(year, month, day, hour, minute, second).expect("有効な UTC 日時")
    }

    /// 最小 GeoPoint。
    fn geo(lat: f64, lon: f64) -> umbra_geo::GeoPoint {
        umbra_geo::GeoPoint::from_degrees(lat, lon).expect("有効な地表点")
    }

    /// 最小 BesselianPolynomial（results.rs の minimal_bessel パターン）。
    fn minimal_bessel() -> crate::bessel_poly::BesselianPolynomial {
        use crate::polynomial::Polynomial;
        let c = |v: f64| Polynomial {
            coefficients: vec![v],
        };
        crate::bessel_poly::BesselianPolynomial {
            epoch_tt: tt(2_451_545.0, 0.0),
            x: c(0.20),
            y: c(-0.30),
            d: c(0.2070),
            mu: c(1.2),
            l1: c(0.5400),
            l2: c(-0.0090),
            tan_f1: 0.004_65,
            tan_f2: 0.004_63,
            fit_interval: umbra_core::TimeInterval {
                start: tt(2_451_544.9, 0.0),
                end: tt(2_451_545.1, 0.0),
            },
            fit_error: crate::bessel_poly::BesselFitError {
                max_x: 1.0e-7,
                max_y: 2.0e-7,
                max_l1: 3.0e-7,
                max_l2: 4.0e-7,
            },
        }
    }

    /// 最小 CalculationMetadata（results.rs の metadata パターン）。
    fn metadata() -> crate::calc_metadata::CalculationMetadata {
        crate::calc_metadata::CalculationMetadata {
            library_version: "0.1.0".to_string(),
            ephemeris_model: "ELP/MPP02+VSOP87D".to_string(),
            ephemeris_version: "2024a".to_string(),
            delta_t_model: "EspenakMeeus".to_string(),
            delta_t_uncertainty_seconds: 0.5,
            earth_model: "WGS84".to_string(),
            lunar_radius_model: "IauMean".to_string(),
            accuracy_profile: crate::config::AccuracyProfile::Standard,
            generated_at: utc(2026, 6, 18, 0, 0, 0.0),
        }
    }

    /// 最小の SolarEclipse（path の引数用・全フィールドは results.rs と同パターン）。
    fn minimal_eclipse() -> SolarEclipse {
        let greatest = GreatestEclipse {
            time_utc: utc(2024, 4, 8, 18, 17, 0.0),
            time_tt: tt(2_460_409.0, 0.123),
            position: geo(25.0, -104.0),
            magnitude: EclipseMagnitude(1.0566),
            obscuration: Obscuration(1.0),
            path_width: None,
            central_duration: None,
            sun_altitude: umbra_core::Degrees(70.3),
        };
        let global = GlobalCircumstances {
            kind: SolarEclipseKind::Total,
            partial_begin: None,
            central_begin: None,
            greatest,
            central_end: None,
            partial_end: None,
            gamma: 0.3431,
        };
        SolarEclipse {
            event_key: "2024-04-08#1252".to_string(),
            kind: SolarEclipseKind::Total,
            global,
            bessel: minimal_bessel(),
            metadata: metadata(),
        }
    }

    // ==================================================================
    // 1. instantaneous_elements: 2017 実日食ゲート（NASA gamma≈0.4367）
    // ==================================================================

    /// 実日食オラクル（最重要）: 合成 TimeData で構築した StandardEngine（new 経由）の
    /// `instantaneous_elements(2017 最大食 TT)` の gamma が NASA 公表 gamma≈0.4367 と 4 桁一致
    /// （`[0.43, 0.44]`）し、`time_tt` ラベルが入力 TT に一致する。
    ///
    /// 殺す変異: 内部評価器の暦/ΔT/半径引数の取り違え・options 既定値改変・退化区間 `[time,time]`
    /// の取り違え・`time_tt` ラベルのずれ・`at(time)` 呼び忘れ。
    #[test]
    fn instantaneous_elements_2017_gamma_matches_nasa() {
        let engine = standard_engine_from_synthetic();
        let e = engine
            .instantaneous_elements(tt_2017_max())
            .expect("2017 最大食での瞬時要素評価は成功する");
        // NASA gamma=0.4367 を [0.43,0.44] で締める（既存実日食ゲートと同値域）。
        assert!(
            (0.43..0.44).contains(&e.gamma()),
            "gamma = {} (NASA 0.4367)",
            e.gamma()
        );
        // time_tt ラベルが入力 TT を保持する（確定仕様: 入力に一致）。
        assert_eq!(e.time_tt, tt_2017_max(), "time_tt label preserved");
    }

    /// 半径配線の絶対値検証: 標準 config の engine の l1/l2 が、独立に算出した標準半径
    /// （IauMean k=0.2725076 × Re[km], Iau2015 太陽 696000 km）で評価した `besselian_elements_at`
    /// と一致する。gamma は半径非依存ゆえ上のテストでは捕捉できない **re_km の m→km 変換
    /// （`/1000`）・`k * re_km` の積**の取り違えを l1/l2 の絶対値で撃破する。
    #[test]
    fn instantaneous_elements_radii_wiring_matches_independent_besselian() {
        use crate::besselian::besselian_elements_at;
        use umbra_core::constants::SOLAR_RADIUS_KM;

        let engine = standard_engine_from_synthetic();
        let t = tt_2017_max();
        let got = engine
            .instantaneous_elements(t)
            .expect("2017 最大食での瞬時要素評価は成功する");

        // 独立に組み立てた標準半径（config の k/radius モデルを介さず直接）。
        let re_km = EARTH_EQUATORIAL_RADIUS_M / 1000.0;
        let r_moon_km = 0.272_507_6 * re_km;
        let want = besselian_elements_at(t, SOLAR_RADIUS_KM, r_moon_km, &EspenakMeeusDeltaT)
            .expect("独立評価は成功する");

        // l1/l2 は半径依存。re_km の % 化（137km）や k+re の和（≈6378km）では桁違いにずれる。
        assert!(
            (got.l1 - want.l1).abs() < 1e-9,
            "l1 = {}, want {} (太陽半径配線)",
            got.l1,
            want.l1
        );
        assert!(
            (got.l2 - want.l2).abs() < 1e-9,
            "l2 = {}, want {} (月半径配線 k×Re)",
            got.l2,
            want.l2
        );
    }

    // ==================================================================
    // 2. instantaneous_elements: 別エポックの Ok・有限・サニティ
    // ==================================================================

    /// 別エポック（J2000）の `instantaneous_elements` が `Ok` を返し、要素が有限で
    /// 半影半径 l1 が正であるサニティを縛る（実太陽・月位置で評価可能）。
    ///
    /// 殺す変異: 特定エポックへのハードコード・非有限/NaN の素通し・l1 符号反転や 0 固定。
    #[test]
    fn instantaneous_elements_other_epoch_is_ok_and_sane() {
        let engine = standard_engine_from_synthetic();
        let e = engine
            .instantaneous_elements(tt_j2000())
            .expect("J2000 での瞬時要素評価は成功する");
        // 主要量が有限。
        assert!(e.x.is_finite(), "x = {}", e.x);
        assert!(e.y.is_finite(), "y = {}", e.y);
        assert!(e.l1.is_finite(), "l1 = {}", e.l1);
        assert!(e.l2.is_finite(), "l2 = {}", e.l2);
        assert!(e.gamma().is_finite(), "gamma = {}", e.gamma());
        // 半影半径 l1 は正（半影は常に正の円錐半径）。
        assert!(e.l1 > 0.0, "l1 = {} (半影は正)", e.l1);
        // time_tt ラベルは入力 TT に一致。
        assert_eq!(e.time_tt, tt_j2000(), "time_tt label preserved (J2000)");
    }

    // ==================================================================
    // 3. path(): 実装済み（M9.1+）— 非中心食ゲート（中心線は中心食のみ）
    // ==================================================================

    /// `path()` は実装済み（M9.1）。中心食でない最小 SolarEclipse（`central_begin`/`central_end`
    /// が `None`）では `Ok(EclipsePath)` を返し、`center_line` は `None`（中心線は中心食のみ）、
    /// `greatest_point` は `global.greatest.position` の passthrough、限界線・部分食域・samples は空。
    ///
    /// 殺す変異: 非中心でも center_line を Some にする・greatest_point を別地点にする・Err を返す。
    #[test]
    fn path_non_central_returns_ok_without_center_line() {
        let engine = standard_engine_from_synthetic();
        let eclipse = minimal_eclipse();
        let path = engine
            .path(&eclipse, PathOptions::default())
            .expect("path は実装済み・Ok を返す");
        assert!(
            path.center_line.is_none(),
            "中心食でない（central_begin/end None）ので center_line は None"
        );
        assert!(
            path.northern_limit.is_none() && path.southern_limit.is_none(),
            "非中心食は限界線を持たない"
        );
        assert!(path.partial_limit.is_none(), "部分食域は未実装で None");
        assert!(path.samples.is_empty(), "samples は未実装で空");
        // greatest_point は global.greatest.position の passthrough。
        assert_eq!(
            path.greatest_point, eclipse.global.greatest.position,
            "greatest_point は global.greatest.position を passthrough"
        );
    }

    /// 中心食でない eclipse では PathOptions に依らず `center_line` は `None`（中心線は中心食ゲートのみ）。
    /// `include_limits=true` を含む非既定 options でも結果が変わらないことで「非中心は options で中心線・
    /// 限界線を生まない」ことを縛る。
    ///
    /// 殺す変異: 非中心でも options 次第で center_line/限界線を生成する。
    #[test]
    fn path_non_central_is_options_invariant() {
        let engine = standard_engine_from_synthetic();
        let eclipse = minimal_eclipse();
        let custom = PathOptions {
            sample_interval_seconds: 5.0,
            include_limits: true,
            split_antimeridian: false,
        };
        let path = engine.path(&eclipse, custom).expect("path は Ok");
        assert!(
            path.center_line.is_none(),
            "非中心食は options に依らず center_line None"
        );
        assert!(
            path.northern_limit.is_none() && path.southern_limit.is_none(),
            "非中心食は include_limits=true でも限界線を生成しない"
        );
    }

    // ==================================================================
    // 4. §78: standard_engine(bundled_time_data()) — 同梱データ（feature-gated）
    // ==================================================================

    /// §78 標準構築（feature-gated）: `standard_engine(bundled_time_data())` が StandardEngine を
    /// 返し、その `instantaneous_elements(2017)` が動作して NASA gamma と 4 桁一致する。
    /// 同梱データ経路（O = time の EOP 複製）でも瞬時要素が正しく出ることを縛る。
    ///
    /// 殺す変異: standard_engine が new/標準 config/EOP 複製を取り違える・同梱経路で評価が壊れる。
    #[cfg(feature = "bundled-data")]
    #[test]
    fn standard_engine_bundled_instantaneous_elements_works() {
        let engine = standard_engine(umbra_ephemeris::bundled_time_data());
        let e = engine
            .instantaneous_elements(tt_2017_max())
            .expect("同梱データの StandardEngine で 2017 瞬時要素は成功する");
        assert!(
            (0.43..0.44).contains(&e.gamma()),
            "bundled standard_engine gamma = {} (NASA 0.4367)",
            e.gamma()
        );
        assert_eq!(
            e.time_tt,
            tt_2017_max(),
            "time_tt label preserved (bundled)"
        );
    }

    // ==================================================================
    // 5. ジェネリック性 / dyn 不使用（確定A1）
    // ==================================================================

    /// StandardEngine の E/D/O 単相化（`AnalyticalEphemeris`/`EspenakMeeusDeltaT`/`IersEopData`）が
    /// コンパイル・構築できる＝`Box<dyn>` 不使用のジェネリック維持（確定A1, 受け入れ §78）。
    /// 型エイリアス `StandardEngine` と `new` の単相化がコンパイルできることで担保する。
    ///
    /// 殺す変異: ジェネリックを `Box<dyn>` 化する・StandardEngine エイリアスの型パラメタを取り違える。
    #[test]
    fn standard_engine_type_is_monomorphic_and_constructs() {
        // 型注釈で StandardEngine への単相化を強制（dyn 化していればここで型不一致になる）。
        let engine: StandardEngine = standard_engine_from_synthetic();
        // 構築した単相エンジンが実際に動く（瞬時要素を 1 回評価できる）。
        let e = engine
            .instantaneous_elements(tt_2017_max())
            .expect("単相 StandardEngine は動作する");
        assert!(e.gamma().is_finite(), "gamma = {}", e.gamma());
    }

    /// `UtcRange` 型エイリアスが `TimeRange<UtcInstant>` として使える（再エクスポート確認）。
    /// 殺す変異: UtcRange エイリアスの欠落・別型への差し替え。
    #[test]
    fn utc_range_alias_is_time_range_of_utc() {
        let start = utc(2024, 1, 1, 0, 0, 0.0);
        let end = utc(2024, 12, 31, 0, 0, 0.0);
        let range: UtcRange = UtcRange { start, end };
        assert_eq!(range.start, start, "UtcRange.start");
        assert_eq!(range.end, end, "UtcRange.end");
    }

    // ==================================================================
    // 6. config の半径モデルが instantaneous_elements に効く（軽い検証）
    // ==================================================================

    /// config の月半径モデルが瞬時要素評価に反映される: 既定（IauMean, k=0.2725076）と
    /// EspenakUmbral（k=0.272281, より小）の config で構築した 2 エンジンの
    /// `instantaneous_elements(2017)` の本影半径 l2 が異なる（月半径が小さいほど本影が変わる）。
    ///
    /// 半径取り違えの精密検証は besselian 側で済むため、ここでは「config の半径モデルが評価に
    /// 伝わっている（無視されていない）」ことだけを軽く縛る。
    /// 殺す変異: config.lunar_radius_model を無視して半径をハードコードする。
    #[test]
    fn config_lunar_radius_model_affects_instantaneous_elements() {
        // 既定（IauMean）。
        let default_engine = standard_engine_from_synthetic();
        let default_e = default_engine
            .instantaneous_elements(tt_2017_max())
            .expect("既定 config で評価成功");

        // 月半径モデルだけ EspenakUmbral（k が小さい）に変えた config。
        let mut umbral_config = EngineConfig::standard();
        umbral_config.lunar_radius_model = LunarRadiusModel::EspenakUmbral;
        let time = synthetic_time_data();
        let eop = time.eop().clone();
        let umbral_engine = EclipseEngine::new(
            AnalyticalEphemeris::new(),
            EspenakMeeusDeltaT,
            eop,
            time,
            umbral_config,
        );
        let umbral_e = umbral_engine
            .instantaneous_elements(tt_2017_max())
            .expect("EspenakUmbral config で評価成功");

        // k が異なる（0.2725076 vs 0.272281）ので本影半径 l2 は変わるはず。
        assert!(
            (default_e.l2 - umbral_e.l2).abs() > 1e-9,
            "月半径モデル差が l2 に効く: IauMean l2={}, EspenakUmbral l2={}",
            default_e.l2,
            umbral_e.l2
        );
    }

    /// config の太陽半径が瞬時要素評価で消費されることを、定数オラクルと併せて軽く縛る。
    /// `EngineConfig::standard().solar_radius_model.radius_km()` が
    /// `EARTH_EQUATORIAL_RADIUS_M/1000` を介した月半径計算とは別物（太陽半径 696000 km）であり、
    /// 既定 config の太陽半径モデルが Iau2015（696000 km）であることを固定する。
    /// （半径そのものの取り違え検出は besselian 側が担うため、ここは config 配線の存在のみ。）
    /// 殺す変異: 既定 config の solar_radius_model を別モデルに差し替える。
    #[test]
    fn standard_config_solar_radius_is_iau2015() {
        let c = EngineConfig::standard();
        assert_eq!(
            c.solar_radius_model,
            SolarRadiusModel::Iau2015,
            "既定 config の太陽半径モデルは Iau2015"
        );
        assert_eq!(
            c.solar_radius_model.radius_km(),
            696_000.0,
            "Iau2015 太陽半径は 696000 km"
        );
        // 月半径計算に使う地球赤道半径定数が想定値（km 換算の足場）。
        assert_eq!(
            EARTH_EQUATORIAL_RADIUS_M / 1000.0,
            6378.137,
            "地球赤道半径 [km]（月半径 = k·Re の足場）"
        );
    }

    // ==================================================================
    // 7. search(): UTC 範囲の日食探索（ISSUE-043 S6c）
    //
    // ## オラクル戦略（full pipeline・実装方針に立ち入らない）
    // - **kind / gamma = NASA 公表値（ballpark）**: 2017-08-21 は皆既（Total）・gamma≈0.4367。
    //   kind は厳密一致（独立事実）、gamma は `[0.40, 0.47]` の ballpark で締める（要素レベルの
    //   既存ゲートより緩い帯。search は最大食を別途解くため要素 gamma と微差しうる）。
    // - **ephemeris_model = 暦自身の metadata().model**（独立オラクル）: ハードコード値ではなく
    //   `AnalyticalEphemeris::new().metadata().model` と一致することで、暦メタの転記を縛る。
    // - **event_key = "YYYY-MM-DD#KKKK" 形式**（独立・形式オラクル）: `#` で 2 分割し、左が
    //   最大食 UTC 日付プレフィクス（`2017-08-21`）、右が非空の数字列（lunation）であることを縛る。
    // - **bessel が最大食を bracket**: `fit_interval.start <= greatest.time_tt <= fit_interval.end`
    //   を JD（`jd2().jd()`）の不等式で縛る（多項式が最大食を含む区間で fit される確定仕様）。
    // - **空結果**: 新月のない短窓 ⇒ 候補ゼロ ⇒ `Vec::new()`（独立・構造オラクル）。
    //
    // ## red 設計（本体未実装）
    // `search` は現状 `Err(EclipseError::NotImplemented)` を返すため、`Ok` を要求する本群は
    // `expect`／`assert!` で失敗する（想定どおりの理由で red）。
    // ==================================================================

    /// 緯度・経度が有限かつ範囲内の妥当な `GeoPoint` であることを縛るヘルパ。
    fn assert_valid_geo_point(p: umbra_geo::GeoPoint) {
        let lat = p.lat.degrees().0;
        let lon = p.lon.degrees().0;
        assert!(lat.is_finite(), "緯度が有限: {lat}");
        assert!(lon.is_finite(), "経度が有限: {lon}");
        assert!((-90.0..=90.0).contains(&lat), "緯度が範囲内: {lat}");
        assert!((-180.0..=180.0).contains(&lon), "経度が範囲内: {lon}");
    }

    /// event_key が `"YYYY-MM-DD#KKKK"` 形式（指定日付プレフィクス ＋ `#` ＋ 非空数字列）か検証する。
    /// 形式オラクル: `#` で 2 分割し、左が `date_prefix`、右がすべて ASCII 数字で非空であること。
    fn assert_event_key_format(event_key: &str, date_prefix: &str) {
        let (date, lunation) = event_key
            .split_once('#')
            .unwrap_or_else(|| panic!("event_key は '#' を含む: {event_key:?}"));
        assert_eq!(
            date, date_prefix,
            "event_key の日付部は最大食 UTC 日付（{date_prefix}）: {event_key:?}"
        );
        assert!(
            !lunation.is_empty() && lunation.chars().all(|c| c.is_ascii_digit()),
            "event_key の lunation 部は非空の数字列: {event_key:?}"
        );
    }

    /// **主要テスト**: 2017-08-21 皆既日食を含む 1 か月窓を探索すると、その日食が見つかる。
    /// kind=Total・gamma≈0.4367（ballpark `[0.40,0.47]`）・全 4 接触 Some（中心皆既）・
    /// event_key 形式・bessel が最大食を bracket・metadata レシピ（暦 metadata と一致・ΔT/地球/月
    /// モデル・標準プロファイル・非空 library_version・正の ΔT 不確かさ）を縛る。
    ///
    /// 殺す変異: pipeline 段の取りこぼし（候補→合→可能性→ソルバ）、kind 誤判定、gamma 計算ミス、
    /// 接触の取りこぼし（中心食で None）、event_key 組み立て誤り、bessel 区間ずれ、metadata 転記漏れ。
    #[test]
    fn search_finds_2017_08_21_total_eclipse() {
        let engine = standard_engine_from_synthetic();
        let range = UtcRange {
            start: utc(2017, 8, 1, 0, 0, 0.0),
            end: utc(2017, 9, 1, 0, 0, 0.0),
        };
        let results = engine.search(range).expect("2017 年 8 月の探索は成功する");

        // 2017-08-21 を最大食日付に持つ日食を 1 つ取り出す（8 月の日食は 1 つ）。
        let eclipse = results
            .iter()
            .find(|e| e.event_key.starts_with("2017-08-21"))
            .expect("2017-08-21 の皆既日食が結果に含まれる");

        // kind は皆既（NASA 事実・厳密一致）。
        assert_eq!(
            eclipse.kind,
            SolarEclipseKind::Total,
            "2017-08-21 は皆既日食（NASA 事実）"
        );

        // gamma は NASA≈0.4367 を ballpark `[0.40, 0.47]` で締める。
        assert!(
            (0.40..=0.47).contains(&eclipse.global.gamma),
            "gamma = {} (NASA≈0.4367)",
            eclipse.global.gamma
        );

        // 最大食: 食分 > 1（皆既）・位置は妥当な地表点。
        assert!(
            eclipse.global.greatest.magnitude.0 > 1.0,
            "皆既なので食分 > 1: {}",
            eclipse.global.greatest.magnitude.0
        );
        assert_valid_geo_point(eclipse.global.greatest.position);

        // 中心皆既なので 4 接触（P1/U1/U4/P4）がすべて Some。
        assert!(
            eclipse.global.partial_begin.is_some(),
            "P1(partial_begin)=Some（中心皆既）"
        );
        assert!(
            eclipse.global.central_begin.is_some(),
            "U1(central_begin)=Some（中心皆既）"
        );
        assert!(
            eclipse.global.central_end.is_some(),
            "U4(central_end)=Some（中心皆既）"
        );
        assert!(
            eclipse.global.partial_end.is_some(),
            "P4(partial_end)=Some（中心皆既）"
        );

        // 4 接触の時系列順 P1 < U1 < U4 < P4（GlobalCircumstances のフィールド swap を撃破）。
        // 全 Some だけでは partial_begin↔central_begin や central_end↔partial_end の取り違えを
        // 見逃すため、TT の単調増加を縛る。
        let p1 = eclipse.global.partial_begin.unwrap().time_tt.jd2().jd();
        let u1 = eclipse.global.central_begin.unwrap().time_tt.jd2().jd();
        let u4 = eclipse.global.central_end.unwrap().time_tt.jd2().jd();
        let p4 = eclipse.global.partial_end.unwrap().time_tt.jd2().jd();
        assert!(
            p1 < u1 && u1 < u4 && u4 < p4,
            "接触の時系列順 P1<U1<U4<P4: P1={p1} U1={u1} U4={u4} P4={p4}"
        );

        // event_key 形式: "2017-08-21#<digits>"。
        assert_event_key_format(&eclipse.event_key, "2017-08-21");

        // bessel が最大食 TT を bracket（多項式 fit 区間が最大食を含む確定仕様）。
        let start_jd = eclipse.bessel.fit_interval.start.jd2().jd();
        let greatest_jd = eclipse.global.greatest.time_tt.jd2().jd();
        let end_jd = eclipse.bessel.fit_interval.end.jd2().jd();
        assert!(
            start_jd <= greatest_jd && greatest_jd <= end_jd,
            "bessel.fit_interval が最大食を bracket: start={start_jd} \
             greatest={greatest_jd} end={end_jd}"
        );

        // fit_interval は全球部分食 [P1,P4]（候補窓ではない）。`(Some,Some)` match arm の削除＝
        // 候補窓フォールバックを撃破する（候補窓も最大食を bracket するため上の bracket 検証だけでは
        // 見逃す）。
        assert_eq!(
            eclipse.bessel.fit_interval.start,
            eclipse.global.partial_begin.unwrap().time_tt,
            "bessel.fit_interval.start は P1"
        );
        assert_eq!(
            eclipse.bessel.fit_interval.end,
            eclipse.global.partial_end.unwrap().time_tt,
            "bessel.fit_interval.end は P4"
        );

        // 半径配線ピン（S5b と同方針）: search の `r_moon = k()*Re`（line 135 の `k()+Re` 取り違えを
        // 撃破）。bessel を最大食 TT で評価した l1/l2 が、**未変異**の `instantaneous_elements`
        // （同一 config・別経路の半径計算）と一致する。gamma は半径非依存ゆえ上の判定では捕捉
        // できない月/太陽半径配線を l1/l2 の絶対値（fit 残差 1e-4 ≪ 1e-3 ≪ 半径取り違え誤差）で縛る。
        let want = engine
            .instantaneous_elements(eclipse.global.greatest.time_tt)
            .expect("最大食 TT の瞬時要素（未変異の半径配線）");
        let got = eclipse
            .bessel
            .at(eclipse.global.greatest.time_tt)
            .expect("bessel 多項式の最大食評価");
        assert!(
            (got.l1 - want.l1).abs() < 1.0e-3,
            "bessel l1={} は独立 instantaneous_elements l1={}（半径配線）に一致",
            got.l1,
            want.l1
        );
        assert!(
            (got.l2 - want.l2).abs() < 1.0e-3,
            "bessel l2={} は独立 instantaneous_elements l2={}（月半径 k×Re 配線）に一致",
            got.l2,
            want.l2
        );

        // metadata レシピ: 暦モデルは暦自身の metadata().model（独立オラクル）と一致。
        let m = &eclipse.metadata;
        assert_eq!(
            m.ephemeris_model,
            AnalyticalEphemeris::new().metadata().model,
            "ephemeris_model は暦の metadata().model を転記"
        );
        assert_eq!(
            m.ephemeris_version,
            AnalyticalEphemeris::new().metadata().version,
            "ephemeris_version は暦の metadata().version を転記"
        );
        assert_eq!(m.delta_t_model, "Espenak-Meeus", "ΔT モデル名");
        assert_eq!(m.earth_model, "WGS84", "地球モデル名");
        assert_eq!(
            m.lunar_radius_model, "IauMean",
            "月半径モデル名（Standard）"
        );
        assert_eq!(
            m.accuracy_profile,
            AccuracyProfile::Standard,
            "精度プロファイルは Standard"
        );
        assert!(!m.library_version.is_empty(), "library_version は非空");
        assert!(
            m.delta_t_uncertainty_seconds > 0.0,
            "ΔT 不確かさは正: {}",
            m.delta_t_uncertainty_seconds
        );
    }

    /// **空結果**: 新月を含まない短窓（中旬・約 3 日）を探索すると候補ゼロ ⇒ `Ok(vec![])`。
    /// 2017-08-21 が新月（皆既）なので、そこから十分離れた 8 月初旬の 3 日窓には新月がない。
    /// （新月を含むが日食でない窓ではなく、確実に候補ゼロな短窓を選ぶ＝偽陽性に強い構造オラクル。）
    ///
    /// 殺す変異: 候補ゼロでも非空を返す・常に Some を返す・空判定の反転。
    #[test]
    fn search_empty_when_no_new_moon_in_window() {
        let engine = standard_engine_from_synthetic();
        // 2017-08-21 の新月から 2 週間ほど前の 3 日窓（新月なし＝候補ゼロ）。
        let range = UtcRange {
            start: utc(2017, 8, 5, 0, 0, 0.0),
            end: utc(2017, 8, 8, 0, 0, 0.0),
        };
        let results = engine
            .search(range)
            .expect("新月のない短窓でも探索は成功する");
        assert!(
            results.is_empty(),
            "新月を含まない短窓は日食ゼロ（空 Vec）: {} 件",
            results.len()
        );
    }

    // ==================================================================
    // 8. 純ヘルパの独立ピン（FAST・search() パイプラインを呼ばない）
    //
    // format_event_key / utc_now / build_metadata は search() の構成部品だが、
    // 直接呼べる純ヘルパ・私有メソッドなので、月単位の遅い探索を回さずに高速に縛る。
    // ==================================================================

    /// `format_event_key` の文字列整形を 3 ケースで exact 検証（FAST・パイプライン非依存）。
    /// 日付は固定幅 `YYYY-MM-DD`、lunation は符号付き 10 進（ゼロ詰めなし）。
    ///
    /// 手計算オラクル:
    /// - `utc(2017,8,21,...)`, k=211 → `"2017-08-21#211"`。
    /// - `utc(1950,9,12,...)`, k=-589 → `"1950-09-12#-589"`（先頭 `-` を確認。`:04` 幅指定なら
    ///   `-589` が `-589`→`0-589` 等に崩れる＝符号付き 10 進であることを縛る）。
    /// - `utc(2000,1,6,...)`, k=0 → `"2000-01-06#0"`（小さい正値にゼロ詰めが入らない）。
    ///
    /// 殺す変異: 日付フィールド幅（`:02`/`:04`）の改変・lunation のゼロ詰め化・区切り `#` の変更。
    #[test]
    fn format_event_key_positive_and_negative_lunation() {
        // 2017-08-21 最大食・lunation 211 → "2017-08-21#211"。
        assert_eq!(
            format_event_key(utc(2017, 8, 21, 18, 25, 0.0), 211),
            "2017-08-21#211"
        );
        // 負（2000 年朔より前）・lunation -589 → "1950-09-12#-589"（先頭 '-' を保持）。
        let neg = format_event_key(utc(1950, 9, 12, 3, 0, 0.0), -589);
        assert_eq!(neg, "1950-09-12#-589");
        assert!(
            neg.starts_with("1950-09-12#-"),
            "負 lunation は符号付き 10 進（先頭 '-'）: {neg:?}"
        );
        // k=0（小さい正値）→ ゼロ詰めなしの "2000-01-06#0"。
        assert_eq!(
            format_event_key(utc(2000, 1, 6, 0, 0, 0.0), 0),
            "2000-01-06#0"
        );
    }

    /// `build_metadata` がレシピ各フィールドを正しく充填する（FAST・search() を呼ばない）。
    /// 暦モデル/版は暦自身の `metadata()`（独立オラクル）、ΔT モデル名・地球/月モデル・
    /// 精度プロファイル・非空 library_version・正の ΔT 不確かさ・生成時刻印の逐語保持を縛る。
    ///
    /// 殺す変異: フィールド転記漏れ・取り違え・generated_at の上書き・モデル名のハードコード誤り。
    #[test]
    fn build_metadata_populates_recipe_fields() {
        let engine = standard_engine_from_synthetic();
        let stamp = utc(2026, 6, 18, 0, 0, 0.0);
        let md = engine.build_metadata(tt_2017_max(), stamp);

        // 暦モデル/版は暦自身の metadata()（独立オラクル）。
        assert_eq!(
            md.ephemeris_model,
            AnalyticalEphemeris::new().metadata().model,
            "ephemeris_model は暦の metadata().model"
        );
        assert_eq!(
            md.ephemeris_version,
            AnalyticalEphemeris::new().metadata().version,
            "ephemeris_version は暦の metadata().version"
        );
        assert_eq!(md.delta_t_model, "Espenak-Meeus", "ΔT モデル名");
        assert_eq!(md.earth_model, "WGS84", "地球モデル名");
        assert_eq!(
            md.lunar_radius_model, "IauMean",
            "月半径モデル名（Standard config）"
        );
        assert_eq!(
            md.accuracy_profile,
            crate::config::AccuracyProfile::Standard,
            "精度プロファイルは Standard"
        );
        assert!(!md.library_version.is_empty(), "library_version は非空");
        assert!(
            md.delta_t_uncertainty_seconds > 0.0,
            "ΔT 不確かさは正: {}",
            md.delta_t_uncertainty_seconds
        );
        // 渡した時刻印は逐語保持（壁時計で上書きしない）。
        assert_eq!(md.generated_at, stamp, "generated_at は渡した値を逐語保持");
    }

    /// `utc_now` が最近の壁時計 UTC を返す（FAST・パイプライン非依存）。
    /// 年が `[2024, 2100)` に収まることで UNIX_EPOCH_JD 定数・秒→JD 換算の足場を縛る
    /// （誤ったエポック JD や単位は現在から大きく外れる）。緩い範囲で時計差を許容する。
    ///
    /// 殺す変異: UNIX_EPOCH_JD の値ずれ・`secs / SECONDS_PER_DAY` の単位取り違え。
    #[test]
    fn utc_now_returns_plausible_recent_instant() {
        let now = utc_now();
        let (y, ..) = now.to_gregorian();
        assert!(
            (2024..2100).contains(&y),
            "utc_now の年が最近（2024..2100）: {y}"
        );
    }

    // ==================================================================
    // 9. local_circumstances(): 観測地点の局地条件（ISSUE-043 S7b-ii）
    //
    // ## オラクル戦略（追認回避・物理事実＋構造契約を主軸）
    // エンジンの内部結線コードを写経しない。取得経路は既存 `search_finds_2017_08_21_total_eclipse`
    // と同一（`standard_engine_from_synthetic()` → `search(2017-08 窓)` →
    // `event_key.starts_with("2017-08-21")` の SolarEclipse 抽出）。その実 SolarEclipse に対し:
    //
    // - **中心食地点**（皆既帯中心付近 37.5°N, 西経89.2°, 標高200m。local_maximum.rs の
    //   central_observer と同緯度経度）: NASA 事実として 2017-08-21 は皆既。よって食分 > 1・
    //   食面積 ≈ 1（>0.99）・最大食時 太陽地平上（北米昼）・visible==true。**中心食地点ゆえ
    //   C1/C2/C3/C4 すべて Some**（皆既帯内に内接 C2/C3 が存在する）。
    // - **部分食地点**（高緯度 60°N, 西経100°。local_maximum.rs の partial_observer 相当）:
    //   部分食（0 < magnitude < 1, obscuration < 1）。中心食地点と magnitude が明確に異なる
    //   （observer 配線が効いている＝lat/lon 取り違え変異を撃破）。**本影外ゆえ C1/C4=Some・
    //   C2/C3=None**（内接なし）。
    // - **見えない観測者**（南半球 −40°S, 東経140°。2017 北米日食帯から十分離れ物理的に非可視）:
    //   visibility==NotVisible・magnitude.0==0.0・obscuration.0==0.0。C1-C4 すべて None（不変）。
    //
    // ## S7b-ii 確定契約（C1-C4 充填 ＋ 可視性精緻化）
    // - **C1-C4 接触集合を充填**: 中心食地点で c1/c2/c3/c4 すべて Some、部分食地点で c1/c4=Some・
    //   c2/c3=None、食なし地点で c1-c4 すべて None。
    // - **接触は時系列順**: 中心食で c1<c2<max<c3<c4（time_tt の JD）、部分食で c1<max<c4。
    // - **可視性精緻化**: 食ありで C1/C4 とも地平上 → FullyVisible（2017 北米中緯度・高緯度地点は
    //   日中ゆえ FullyVisible 想定）。食なし → NotVisible（不変）。
    // - **maximum・magnitude・obscuration・metadata は S7b-i と同じ（不変）**。
    //
    // ## 追認回避（独立オラクル）
    // - 物理事実: 2017-08-21 は皆既（NASA）。中心食地点で 4 接触・順序・FullyVisible、
    //   部分食地点で c2/c3=None という幾何的事実。接触順序 C1<C2<最大<C3<C4 は日食の幾何の独立事実。
    // - PA は値域/象限のみ（厳密値は使わない）。内接 C2/C3（皆既 σ=−1）と外接 C1/C4（σ=+1）で
    //   接触点の向きが反転する構造（PA が一定以上異なる）を縛る（式は写経しない）。
    // - solve_local_contacts（ISSUE-025・独立検証済）を test 内で直接呼んで突合するのは同一
    //   プリミティブの追認になりうるので避け、物理的順序・Some/None 構造・可視性で縛る。
    //
    // ## red 設計（S7b-i 実装に対する想定 red）
    // 現状 S7b-i は c1-c4=None・可視日食でも PartialVisible を返すため、本群の更新後アサーション
    // （c1-c4 Some・FullyVisible）は assertion 失敗で red になる。invisible/anchor テスト
    // （c1-c4=None・NotVisible を縛る）は緑のまま。
    // ==================================================================

    /// 2017-08-21 皆既日食を search で取得し、その SolarEclipse を返すヘルパ（取得経路は
    /// `search_finds_2017_08_21_total_eclipse` と同一・追認なしの独立取得）。
    fn search_2017_total(engine: &StandardEngine) -> SolarEclipse {
        let range = UtcRange {
            start: utc(2017, 8, 1, 0, 0, 0.0),
            end: utc(2017, 9, 1, 0, 0, 0.0),
        };
        let results = engine.search(range).expect("2017 年 8 月の探索は成功する");
        results
            .into_iter()
            .find(|e| e.event_key.starts_with("2017-08-21"))
            .expect("2017-08-21 の皆既日食が結果に含まれる")
    }

    /// 中心食地点の観測者（皆既帯中心付近 37.5°N, 西経89.2°, 標高200m）。
    /// local_maximum.rs の central_observer と同緯度経度（西経は負）。
    fn central_observer() -> Observer {
        Observer::from_degrees(37.5, -89.2, 200.0).expect("有効な中心食地点観測者")
    }

    /// 部分食のみの観測者（高緯度 60°N, 西経100°）。local_maximum.rs の partial_observer 相当。
    fn partial_observer() -> Observer {
        Observer::from_degrees(60.0, -100.0, 200.0).expect("有効な部分食地点観測者")
    }

    /// 2017 北米日食帯から十分離れた南半球の観測者（−40°S, 東経140°）＝物理的に非可視。
    /// 2017-08-21 の食域（北米・北大西洋）から地球の反対側に近く、月影は到達しない。
    fn invisible_observer() -> Observer {
        Observer::from_degrees(-40.0, 140.0, 0.0).expect("有効な非可視地点観測者")
    }

    /// **主要テスト（中心食地点・皆既）**: 2017-08-21 を中心食地点で見ると皆既条件になる。
    /// 食分 > 1（皆既=NASA 事実）・食面積 ≈ 1（>0.99）・最大食 太陽地平上（北米昼・alt>0）・
    /// visible==true・最大食 UTC は 18 時台 ballpark（17.5〜19.0 時）・
    /// maximum_altitude==maximum.sun_altitude。
    ///
    /// **S7b-ii 確定契約**: 中心食地点ゆえ C1/C2/C3/C4 すべて Some（皆既帯内に内接 C2/C3 が存在）。
    /// 4 接触は日中ゆえ各々 visible==true。接触の時系列順 c1<c2<max<c3<c4（time_tt の JD 単調増加）。
    /// 食あり・全接触地平上ゆえ visibility==FullyVisible（日中の中心食）。
    ///
    /// 殺す変異: observer の lat/lon 取り違え（部分食/非可視と magnitude 差で別途撃破）、
    /// ζ 補正欠落（皆既で magnitude>1・obscuration≈1 を縛る）、maximum_altitude と sun_altitude の
    /// 不一致、可視性の取り違え（FullyVisible を縛る）、C1-C4 を誤って None のまま放置する
    /// （S7b-i の取りこぼし）、c1↔c4 / c2↔c3 の時系列取り違え。
    #[test]
    fn local_circumstances_central_site_is_total() {
        let engine = standard_engine_from_synthetic();
        let eclipse = search_2017_total(&engine);
        let lc = engine
            .local_circumstances(&eclipse, central_observer())
            .expect("中心食地点の局地条件は成功する");

        // 皆既（NASA 事実）: 食分 > 1・食面積 ≈ 1。ζ 補正欠落だとここが崩れる。
        assert!(
            lc.magnitude.0 > 1.0,
            "中心食地点は皆既なので食分 > 1: {}",
            lc.magnitude.0
        );
        assert!(
            lc.obscuration.0 > 0.99,
            "皆既なので食面積 ≈ 1（>0.99）: {}",
            lc.obscuration.0
        );

        // 最大食時 太陽は地平上（北米昼）・visible==true。
        let mx = lc.contacts.maximum;
        assert!(
            mx.sun_altitude.0 > 0.0,
            "中心食地点の最大食時 太陽は地平上: alt = {}",
            mx.sun_altitude.0
        );
        assert!(mx.visible, "太陽地平上ゆえ visible==true");

        // maximum_altitude は maximum.sun_altitude と一致（別フィールドからの転記漏れを撃破）。
        assert_eq!(
            lc.maximum_altitude, mx.sun_altitude,
            "maximum_altitude == contacts.maximum.sun_altitude"
        );

        // S7b-ii: 食あり・全接触地平上（日中の中心食）ゆえ FullyVisible。
        assert_eq!(
            lc.visibility,
            Visibility::FullyVisible,
            "S7b-ii: 食あり・C1/C4 とも地平上（日中）ゆえ FullyVisible"
        );

        // 最大食 UTC は 2017-08-21 18 時台 ballpark（17.5〜19.0 時）。窓/時刻取り違えを撃破。
        let (y, mo, d, h, mi, _) = mx.time_utc.to_gregorian();
        assert_eq!((y, mo, d), (2017, 8, 21), "最大食 UTC 日付は 2017-08-21");
        let hour_frac = f64::from(h) + f64::from(mi) / 60.0;
        assert!(
            (17.5..=19.0).contains(&hour_frac),
            "最大食 UTC は 18 時台 ballpark（17.5〜19.0h）: {hour_frac}h"
        );

        // S7b-ii: 中心食地点（皆既帯内）ゆえ C1/C2/C3/C4 すべて Some。
        // None のまま放置（S7b-i の取りこぼし）を撃破する。
        let c1 = lc.contacts.c1.expect("中心食地点: C1=Some");
        let c2 = lc.contacts.c2.expect("中心食地点: C2=Some（皆既内接あり）");
        let c3 = lc.contacts.c3.expect("中心食地点: C3=Some（皆既内接あり）");
        let c4 = lc.contacts.c4.expect("中心食地点: C4=Some");

        // 4 接触は日中ゆえ各々 visible==true（太陽地平上）。
        assert!(c1.visible, "C1 は日中ゆえ visible");
        assert!(c2.visible, "C2 は日中ゆえ visible");
        assert!(c3.visible, "C3 は日中ゆえ visible");
        assert!(c4.visible, "C4 は日中ゆえ visible");

        // 接触の時系列順 c1 < c2 < max < c3 < c4（time_tt の JD 単調増加）。
        // 日食の幾何的事実（独立オラクル）。c1↔c4 / c2↔c3 の取り違えを撃破する。
        let j_c1 = c1.time_tt.jd2().jd();
        let j_c2 = c2.time_tt.jd2().jd();
        let j_mx = mx.time_tt.jd2().jd();
        let j_c3 = c3.time_tt.jd2().jd();
        let j_c4 = c4.time_tt.jd2().jd();
        assert!(
            j_c1 < j_c2 && j_c2 < j_mx && j_mx < j_c3 && j_c3 < j_c4,
            "接触の時系列順 c1<c2<max<c3<c4: c1={j_c1} c2={j_c2} max={j_mx} c3={j_c3} c4={j_c4}"
        );
    }

    /// **部分食地点**: 高緯度 60°N/西経100° は本影外＝部分食。0 < magnitude < 1・obscuration < 1。
    /// 中心食地点と magnitude が明確に異なる（observer 配線が効いている）。
    ///
    /// **S7b-ii 確定契約**: 部分食地点は本影外ゆえ C1/C4=Some・**C2/C3=None**（内接なし）。
    /// 接触の時系列順 c1<max<c4（部分食地点）。可視性は「見える種別」（食あり・地平上ゆえ
    /// NotVisible でない。北米日中の高緯度ゆえ FullyVisible 想定）。
    ///
    /// 殺す変異: observer の lat/lon を無視して常に中心食地点の値を返す（magnitude 差で撃破）、
    /// 部分食地点を皆既扱い（magnitude<1 を縛る）、食なし扱い（magnitude>0 を縛る）、
    /// 部分食地点で C2/C3 を誤って Some にする（内接ありと取り違え）、C1/C4 を None のまま放置。
    #[test]
    fn local_circumstances_partial_site_is_partial() {
        let engine = standard_engine_from_synthetic();
        let eclipse = search_2017_total(&engine);

        let central = engine
            .local_circumstances(&eclipse, central_observer())
            .expect("中心食地点の局地条件は成功する");
        let partial = engine
            .local_circumstances(&eclipse, partial_observer())
            .expect("部分食地点の局地条件は成功する");

        // 部分食: 0 < magnitude < 1・obscuration < 1。
        assert!(
            partial.magnitude.0 > 0.0,
            "部分食地点でも食はある: magnitude = {}",
            partial.magnitude.0
        );
        assert!(
            partial.magnitude.0 < 1.0,
            "部分食地点は皆既でない: magnitude = {}",
            partial.magnitude.0
        );
        assert!(
            partial.obscuration.0 < 1.0,
            "部分食地点の食面積 < 1: {}",
            partial.obscuration.0
        );

        // observer 配線が効いている: 中心食地点と部分食地点で magnitude が明確に異なる
        // （lat/lon 取り違え＝両地点を同一視する変異を撃破）。
        assert!(
            (central.magnitude.0 - partial.magnitude.0).abs() > 0.05,
            "観測者で結果が変わる（中心食 {} vs 部分食 {}）",
            central.magnitude.0,
            partial.magnitude.0
        );

        // S7b-ii: 部分食地点は外接 C1/C4=Some・内接 C2/C3=None（本影外ゆえ内接なし）。
        let c1 = partial.contacts.c1.expect("部分食地点: C1=Some（外接）");
        let c4 = partial.contacts.c4.expect("部分食地点: C4=Some（外接）");
        assert!(
            partial.contacts.c2.is_none(),
            "部分食地点: C2=None（本影外ゆえ内接なし）"
        );
        assert!(
            partial.contacts.c3.is_none(),
            "部分食地点: C3=None（本影外ゆえ内接なし）"
        );

        // 接触の時系列順 c1 < max < c4（部分食地点・独立な幾何的事実）。
        let j_c1 = c1.time_tt.jd2().jd();
        let j_mx = partial.contacts.maximum.time_tt.jd2().jd();
        let j_c4 = c4.time_tt.jd2().jd();
        assert!(
            j_c1 < j_mx && j_mx < j_c4,
            "部分食の接触順 c1<max<c4: c1={j_c1} max={j_mx} c4={j_c4}"
        );

        // 可視性は「見える種別」: 食あり・地平上ゆえ NotVisible でない（緩い構造契約）。
        assert_ne!(
            partial.visibility,
            Visibility::NotVisible,
            "部分食地点（食あり・日中）は NotVisible でない"
        );
    }

    /// **見えない観測者**: 2017 北米日食帯から離れた南半球（−40°S, 東経140°）は非可視。
    /// 観測可能な契約（内部分岐に依らず保証されるべき値）を縛る: visibility==NotVisible・
    /// magnitude.0==0.0・obscuration.0==0.0・C1-C4 None。
    ///
    /// 注（分岐）: この地点は探索窓 [P1,P4] 内に局地最接近の極小を持つため `solve_local_maximum`
    /// は成功し（`Ok` 分岐）、`min_sep ≥ L1` で食分 0＝NotVisible になる（maximum は局地最接近時刻
    /// ＝窓内）。全球最大食時刻への錨打ちは「窓内に極小をブラケットできない」遠方観測者
    /// （`RootNotBracketed` 分岐）専用であり、本地点はそちらを通らない。よって maximum.time_tt は
    /// 局地最接近時刻（[P1,P4] 内）であって全球最大食 TT と一致するとは限らない。ここでは
    /// **非可視の観測可能契約**（NotVisible・食分/食面積 0）と maximum が窓内であることを縛る。
    ///
    /// 殺す変異: 非可視地点で食ありを返す（magnitude/obscuration 0 を縛る）、可視性分類の取り違え
    /// （NotVisible を縛る）、maximum を窓外の時刻にする。
    #[test]
    fn local_circumstances_invisible_site_is_not_visible() {
        let engine = standard_engine_from_synthetic();
        let eclipse = search_2017_total(&engine);
        let lc = engine
            .local_circumstances(&eclipse, invisible_observer())
            .expect("非可視地点でも局地条件は成功する（NotVisible を返す）");

        // 食なし（観測可能契約・内部分岐に依らない）。
        assert_eq!(
            lc.magnitude.0, 0.0,
            "非可視地点は食なし: magnitude = {}",
            lc.magnitude.0
        );
        assert_eq!(
            lc.obscuration.0, 0.0,
            "非可視地点は食面積 0: {}",
            lc.obscuration.0
        );
        assert_eq!(
            lc.visibility,
            Visibility::NotVisible,
            "非可視地点は NotVisible"
        );

        // maximum は全球部分食窓 [P1,P4] 内（局地最接近 or 全球最大食錨。窓外時刻への取り違えを撃破）。
        let p1_jd = eclipse
            .global
            .partial_begin
            .expect("2017 は全球 P1 を持つ")
            .time_tt
            .jd2()
            .jd();
        let p4_jd = eclipse
            .global
            .partial_end
            .expect("2017 は全球 P4 を持つ")
            .time_tt
            .jd2()
            .jd();
        let max_jd = lc.contacts.maximum.time_tt.jd2().jd();
        assert!(
            p1_jd <= max_jd && max_jd <= p4_jd,
            "非可視地点の maximum は全球部分食窓 [P1,P4] 内: P1={p1_jd} max={max_jd} P4={p4_jd}"
        );

        // S7b-i: C1-C4 接触集合は常に None。
        assert!(lc.contacts.c1.is_none(), "S7b-i: c1=None");
        assert!(lc.contacts.c2.is_none(), "S7b-i: c2=None");
        assert!(lc.contacts.c3.is_none(), "S7b-i: c3=None");
        assert!(lc.contacts.c4.is_none(), "S7b-i: c4=None");
    }

    /// **錨分岐（窓内でブラケット不能な観測者）**: `solve_local_maximum` が `RootNotBracketed`
    /// を返す場合、maximum を全球最大食時刻に錨打ちし NotVisible・食分0 を返す（S7b 確定）。
    /// **退化窓 `[t,t]`（fit_interval の start==end）で本分岐を決定的に励起する**: 退化窓では
    /// 粗走査の全サンプルが同値→最小が窓端 (min_i==0)→機構的にブラケット不成立（local_maximum.rs
    /// の `constant_m2_flat_bottom_is_graceful_root_not_bracketed` と同経路）。partial_begin/end は
    /// `minimal_eclipse` で None ゆえ探索窓は bessel.fit_interval へフォールバックする。
    ///
    /// 殺す変異: 錨時刻を greatest 以外にする（time_tt/time_utc 参照差し替え）、錨分岐の magnitude/
    /// obscuration を 0 以外にする、visibility を NotVisible 以外にする、錨分岐そのものの削除。
    #[test]
    fn local_circumstances_unbracketable_window_anchors_at_global_greatest() {
        let engine = standard_engine_from_synthetic();
        // 退化窓 [greatest, greatest] の合成日食（partial_begin/end=None → 窓=bessel.fit_interval）。
        let mut eclipse = minimal_eclipse();
        let anchor_tt = eclipse.global.greatest.time_tt;
        eclipse.bessel.fit_interval = umbra_core::TimeInterval {
            start: anchor_tt,
            end: anchor_tt,
        };
        // 前提: partial_begin/end が None（窓は bessel.fit_interval へフォールバックする）。
        assert!(
            eclipse.global.partial_begin.is_none() && eclipse.global.partial_end.is_none(),
            "前提: 合成日食は全球接触 P1/P4 を持たない（窓は fit_interval）"
        );

        let lc = engine
            .local_circumstances(&eclipse, central_observer())
            .expect("退化窓でも局地条件は成功する（NotVisible の錨を返す）");

        // 錨: maximum は全球最大食時刻（TT/UTC 両方）に一致。
        assert_eq!(
            lc.contacts.maximum.time_tt, eclipse.global.greatest.time_tt,
            "錨分岐: maximum.time_tt は全球最大食 TT"
        );
        assert_eq!(
            lc.contacts.maximum.time_utc, eclipse.global.greatest.time_utc,
            "錨分岐: maximum.time_utc は全球最大食 UTC"
        );
        // 食なし・非可視。
        assert_eq!(lc.magnitude.0, 0.0, "錨分岐は食分 0");
        assert_eq!(lc.obscuration.0, 0.0, "錨分岐は食面積 0");
        assert_eq!(lc.visibility, Visibility::NotVisible, "錨分岐は NotVisible");
        // C1-C4 None。
        assert!(
            lc.contacts.c1.is_none()
                && lc.contacts.c2.is_none()
                && lc.contacts.c3.is_none()
                && lc.contacts.c4.is_none(),
            "S7b-i: C1-C4 は None"
        );
    }

    /// **UTC/TT 整合**: maximum.time_utc == tt_to_utc(maximum.time_tt)（post-1972）。
    /// 中心食地点で確認（time_utc/time_tt の取り違えを撃破）。
    ///
    /// 殺す変異: time_utc に time_tt を入れる（または逆）、UTC 変換忘れ。
    #[test]
    fn local_circumstances_maximum_utc_tt_consistent() {
        let engine = standard_engine_from_synthetic();
        let eclipse = search_2017_total(&engine);
        let lc = engine
            .local_circumstances(&eclipse, central_observer())
            .expect("中心食地点の局地条件は成功する");

        let mx = lc.contacts.maximum;
        let expected_utc = umbra_core::time::tt_to_utc(mx.time_tt)
            .expect("最大食 TT は post-1972 で UTC 変換可能");
        // 同一瞬時（JD 差 < 1ms 相当）。
        let got_jd = mx.time_utc.jd2().jd();
        let want_jd = expected_utc.jd2().jd();
        assert!(
            (got_jd - want_jd).abs() < 1.0 / SECONDS_PER_DAY,
            "maximum.time_utc == tt_to_utc(maximum.time_tt): got_jd={got_jd} want_jd={want_jd}"
        );
    }

    /// **alt/az/PA の健全性**: 中心食地点の最大食 ＋ C1-C4 各 LocalContact について、
    /// sun_altitude ∈ [-90,90]・sun_azimuth ∈ [0,360)・position_angle ∈ [0,360) かつすべて有限。
    /// 厳密値は flaky ゆえ範囲のみ縛る（最大食 σ=+1・内接 C2/C3 σ=−1 でも値域は同じ）。
    ///
    /// 殺す変異: 角度の値域違反（ラジアン/度取り違え等）、PA を非有限/未計算で素通し、
    /// C1-C4 接触の角度値が値域外（接触ごとの座標計算崩れ）。
    #[test]
    fn local_circumstances_angles_in_valid_ranges() {
        let engine = standard_engine_from_synthetic();
        let eclipse = search_2017_total(&engine);
        let lc = engine
            .local_circumstances(&eclipse, central_observer())
            .expect("中心食地点の局地条件は成功する");

        // 角度値域を 1 つの LocalContact について縛るヘルパ。
        let assert_angles = |c: &LocalContact, label: &str| {
            assert!(c.sun_altitude.0.is_finite(), "{label}: sun_altitude 有限");
            assert!(c.sun_azimuth.0.is_finite(), "{label}: sun_azimuth 有限");
            assert!(
                c.position_angle.0.is_finite(),
                "{label}: position_angle 有限"
            );
            assert!(
                (-90.0..=90.0).contains(&c.sun_altitude.0),
                "{label}: sun_altitude ∈ [-90,90]: {}",
                c.sun_altitude.0
            );
            assert!(
                (0.0..360.0).contains(&c.sun_azimuth.0),
                "{label}: sun_azimuth ∈ [0,360): {}",
                c.sun_azimuth.0
            );
            assert!(
                (0.0..360.0).contains(&c.position_angle.0),
                "{label}: position_angle ∈ [0,360): {}",
                c.position_angle.0
            );
        };

        // 最大食。
        assert_angles(&lc.contacts.maximum, "maximum");
        // C1-C4（中心食地点ゆえすべて Some）。各接触の角度も値域内。
        assert_angles(&lc.contacts.c1.expect("中心食 C1=Some"), "c1");
        assert_angles(&lc.contacts.c2.expect("中心食 C2=Some"), "c2");
        assert_angles(&lc.contacts.c3.expect("中心食 C3=Some"), "c3");
        assert_angles(&lc.contacts.c4.expect("中心食 C4=Some"), "c4");
    }

    /// **中心食 PA の内接/外接差（S7b-ii 固有）**: 中心食地点（皆既）の内接 C2/C3（皆既 σ=−1・
    /// 接触点が月中心の反対側）と外接 C1/C4（σ=+1・接触点が月中心方向）で接触点の向きが反転し、
    /// PA が一定以上異なる。厳密値は flaky ゆえ「内接と外接の PA が十分離れている」構造で縛る
    /// （PA の式は写経しない＝追認回避）。σ 反転（内接で −1 を掛ける）は PA を約 180° 回すため、
    /// 内接 C2 と外接 C1 の PA 差は周期 [0,180] に畳んで 90° 超になることを期待する。
    ///
    /// 殺す変異: 内接 C2/C3 に外接と同じ σ=+1 を使う（umbral_interior フラグの取り違え）、
    /// 全接触で PA を同一値にする（接触ごとの座標を使わない）。
    #[test]
    fn local_circumstances_central_inner_outer_pa_differ() {
        let engine = standard_engine_from_synthetic();
        let eclipse = search_2017_total(&engine);
        let lc = engine
            .local_circumstances(&eclipse, central_observer())
            .expect("中心食地点の局地条件は成功する");

        let c1 = lc.contacts.c1.expect("中心食 C1=Some");
        let c2 = lc.contacts.c2.expect("中心食 C2=Some（内接）");

        // PA 周期差を [0,180] に畳む（北0近傍の巻き戻り・180°反転を頑健に測る）。
        let folded_diff = |a: f64, b: f64| -> f64 {
            let mut d = (a - b).rem_euclid(360.0);
            if d > 180.0 {
                d = 360.0 - d;
            }
            d
        };

        // 内接 C2（皆既 σ=−1）と外接 C1（σ=+1）は接触点の向きが反転 → PA が大きく異なる。
        // σ 反転は約 180° 回すため、畳んだ差は 90° 超を期待する（内接/外接で同 σ を使う変異を撃破）。
        let diff_c2_c1 = folded_diff(c2.position_angle.0, c1.position_angle.0);
        assert!(
            diff_c2_c1 > 90.0,
            "内接 C2(σ=−1) と外接 C1(σ=+1) の PA は十分異なる（畳んだ差 {diff_c2_c1}° > 90°）: \
             C2={} C1={}",
            c2.position_angle.0,
            c1.position_angle.0
        );
    }

    /// **接触の UTC/TT 整合（C1）**: ある接触（C1）で time_utc == tt_to_utc(time_tt)。
    /// 最大食以外の接触でも時刻系の組が整合していることを縛る。
    ///
    /// 殺す変異: 接触の time_utc に time_tt を入れる（または逆）、接触の UTC 変換忘れ。
    #[test]
    fn local_circumstances_contact_utc_tt_consistent() {
        let engine = standard_engine_from_synthetic();
        let eclipse = search_2017_total(&engine);
        let lc = engine
            .local_circumstances(&eclipse, central_observer())
            .expect("中心食地点の局地条件は成功する");

        let c1 = lc.contacts.c1.expect("中心食 C1=Some");
        let expected_utc =
            umbra_core::time::tt_to_utc(c1.time_tt).expect("C1 の TT は post-1972 で UTC 変換可能");
        let got_jd = c1.time_utc.jd2().jd();
        let want_jd = expected_utc.jd2().jd();
        assert!(
            (got_jd - want_jd).abs() < 1.0 / SECONDS_PER_DAY,
            "C1.time_utc == tt_to_utc(C1.time_tt): got_jd={got_jd} want_jd={want_jd}"
        );
    }

    /// **中心食 vs 部分食で内接 C2/C3 の Some/None が分かれる（観測者で接触構造が変わる）**:
    /// 同一日食に対し中心食地点では C2/C3 が Some、部分食地点では C2/C3 が None になる。
    /// 観測者によって接触集合の構造（内接の有無）が切り替わることを 1 テストで対比して縛る。
    ///
    /// 殺す変異: 観測者に依らず常に C2/C3 を Some（または常に None）にする、内接判定が
    /// observer を見ていない。
    #[test]
    fn local_circumstances_inner_contacts_depend_on_observer() {
        let engine = standard_engine_from_synthetic();
        let eclipse = search_2017_total(&engine);

        let central = engine
            .local_circumstances(&eclipse, central_observer())
            .expect("中心食地点の局地条件は成功する");
        let partial = engine
            .local_circumstances(&eclipse, partial_observer())
            .expect("部分食地点の局地条件は成功する");

        // 中心食地点: 内接 C2/C3 あり。
        assert!(
            central.contacts.c2.is_some() && central.contacts.c3.is_some(),
            "中心食地点は内接 C2/C3=Some"
        );
        // 部分食地点: 内接 C2/C3 なし。
        assert!(
            partial.contacts.c2.is_none() && partial.contacts.c3.is_none(),
            "部分食地点は内接 C2/C3=None"
        );
        // 外接 C1/C4 はどちらの地点でも存在（食ありゆえ部分食の外接は共通）。
        assert!(
            central.contacts.c1.is_some() && central.contacts.c4.is_some(),
            "中心食地点も外接 C1/C4=Some"
        );
        assert!(
            partial.contacts.c1.is_some() && partial.contacts.c4.is_some(),
            "部分食地点も外接 C1/C4=Some"
        );
    }

    /// **metadata レシピ**: local_circumstances の metadata が search と同じレシピで充填される。
    /// 暦モデル/版は暦自身の metadata()（独立オラクル）、ΔT モデル名="Espenak-Meeus"・地球="WGS84"・
    /// 月半径="IauMean"・精度=Standard・非空 library_version・正の ΔT 不確かさ。
    ///
    /// 殺す変異: metadata 転記漏れ・モデル名のハードコード誤り・ΔT 不確かさ 0 固定。
    #[test]
    fn local_circumstances_metadata_recipe() {
        let engine = standard_engine_from_synthetic();
        let eclipse = search_2017_total(&engine);
        let lc = engine
            .local_circumstances(&eclipse, central_observer())
            .expect("中心食地点の局地条件は成功する");

        let m = &lc.metadata;
        assert_eq!(
            m.ephemeris_model,
            AnalyticalEphemeris::new().metadata().model,
            "ephemeris_model は暦の metadata().model を転記"
        );
        assert_eq!(
            m.ephemeris_version,
            AnalyticalEphemeris::new().metadata().version,
            "ephemeris_version は暦の metadata().version を転記"
        );
        assert_eq!(m.delta_t_model, "Espenak-Meeus", "ΔT モデル名");
        assert_eq!(m.earth_model, "WGS84", "地球モデル名");
        assert_eq!(
            m.lunar_radius_model, "IauMean",
            "月半径モデル名（Standard）"
        );
        assert_eq!(
            m.accuracy_profile,
            AccuracyProfile::Standard,
            "精度プロファイルは Standard"
        );
        assert!(!m.library_version.is_empty(), "library_version は非空");
        assert!(
            m.delta_t_uncertainty_seconds > 0.0,
            "ΔT 不確かさは正: {}",
            m.delta_t_uncertainty_seconds
        );
    }

    // ==================================================================
    // 8. next_visible_eclipse / next_visible_is_observable（ISSUE-043 S8）
    //
    // ## オラクル戦略（実装方針に立ち入らず、確定仕様の公開 IF だけを縛る）
    // - **純ヘルパ `next_visible_is_observable`（主力・FAST）**: 「地平上で日食を観測できる
    //   高度状態か」の定義そのものを 6 値で直接表化する独立オラクル。FullyVisible/PartialVisible/
    //   SunriseEclipse/SunsetEclipse は太陽が（一部でも）地平上で食が観測できる → true、
    //   NotVisible（食域外）/BelowHorizon（最大食も地平下）は観測不能 → false。
    //   実装の `matches!` を写経せず、各バリアントの意味（地平上で観測可能か）から true/false を決める。
    // - **統合 happy-path（SLOW・1件）**: 物理事実（2017-08-21 は北米中緯度=central_observer で
    //   皆既可視, NASA）と構造（Option/event_key/可視種別/local 整合）で縛る。`after` を日食の直前
    //   （2017-08-01）に置き、central_observer で皆既可視日食が最初の探索窓で見つかるようにして
    //   解く日食を 1 件に抑える（コスト最小化）。
    // - **統合 skip（SLOW・1件）**: invisible_observer（−40°S,140°E。2017-08-21 は NotVisible=
    //   既存 local_circumstances_invisible_site_is_not_visible で確認済）で呼ぶと 2017-08-21 を
    //   **スキップ**して後続の可視日食を返す。これは「可視性を見ずに最初の日食を返す」バグの
    //   唯一のガード。skip 先の具体日付は geography 依存で flaky ゆえハードコードせず、
    //   「2017-08-21 でない可視日食」という構造で縛る。
    // - **None（horizon 内に可視日食なし）はテストしない**: 全 horizon 走査が極めて遅いため。
    //
    // ## red 設計（本体未実装）
    // `next_visible_eclipse` / `next_visible_is_observable` は現状 `unimplemented!` ゆえ panic で
    // red になる。純ヘルパテストは `super::next_visible_is_observable(...)` 呼び出しで即 panic。
    // 統合テストも `next_visible_eclipse` の入口（unimplemented!）で即 panic するため、search の
    // 実走（数百秒）は **red 段階では発生しない**（unimplemented で即落ちる）。
    // ==================================================================

    // ------------------------------------------------------------------
    // 8a. next_visible_is_observable（純関数・FAST・search を呼ばない・主力）
    // ------------------------------------------------------------------

    /// **見える種別 → true（FullyVisible）**: 全経過が地平上 → 観測可能。
    /// 殺す変異: FullyVisible アームを false にする・反転する。
    #[test]
    fn next_visible_is_observable_fully_visible_is_true() {
        assert!(
            super::next_visible_is_observable(Visibility::FullyVisible),
            "FullyVisible は地平上で全経過観測可能 → true"
        );
    }

    /// **見える種別 → true（PartialVisible）**: 一部の接触のみ地平上でも食は観測可能。
    /// 殺す変異: PartialVisible アームを false にする・反転する。
    #[test]
    fn next_visible_is_observable_partial_visible_is_true() {
        assert!(
            super::next_visible_is_observable(Visibility::PartialVisible),
            "PartialVisible は一部地平上で観測可能 → true"
        );
    }

    /// **見える種別 → true（SunriseEclipse）**: 日の出中に食が進行＝地平上で観測可能。
    /// 殺す変異: SunriseEclipse アームを false にする・反転する。
    #[test]
    fn next_visible_is_observable_sunrise_eclipse_is_true() {
        assert!(
            super::next_visible_is_observable(Visibility::SunriseEclipse),
            "SunriseEclipse は日の出中に食を観測可能 → true"
        );
    }

    /// **見える種別 → true（SunsetEclipse）**: 日没中に食が終了＝地平上で観測可能。
    /// 殺す変異: SunsetEclipse アームを false にする・反転する。
    #[test]
    fn next_visible_is_observable_sunset_eclipse_is_true() {
        assert!(
            super::next_visible_is_observable(Visibility::SunsetEclipse),
            "SunsetEclipse は日没中に食を観測可能 → true"
        );
    }

    /// **見えない種別 → false（NotVisible）**: 食域外＝そもそも食がない → 観測不能。
    /// 殺す変異: NotVisible アームを true にする・反転する（最重要・skip ロジックの根拠）。
    #[test]
    fn next_visible_is_observable_not_visible_is_false() {
        assert!(
            !super::next_visible_is_observable(Visibility::NotVisible),
            "NotVisible は食域外 → 観測不能 → false"
        );
    }

    /// **見えない種別 → false（BelowHorizon）**: 最大食も含め全接触が地平下 → 観測不能。
    /// 殺す変異: BelowHorizon アームを true にする・反転する。
    #[test]
    fn next_visible_is_observable_below_horizon_is_false() {
        assert!(
            !super::next_visible_is_observable(Visibility::BelowHorizon),
            "BelowHorizon は最大食も地平下 → 観測不能 → false"
        );
    }

    /// **網羅メタ確認（全 6 値の table 一括）**: 6 バリアントを 1 つの真理値表で一括検証する。
    /// 「観測可能な高度状態か」の定義を 6 値で直接表化した独立オラクル（実装の matches! を写経せず
    /// 意味＝地平上で観測できるか から true/false を決める）。個別テストの取りこぼし防止と、
    /// 「全部 true / 全部 false に潰す」変異の撃破を兼ねる。
    /// 殺す変異: いずれかのアームの true/false 取り違え・反転、定数 true/定数 false への退化。
    #[test]
    fn next_visible_is_observable_truth_table_is_exhaustive() {
        // (variant, 地平上で日食を観測できるか) を意味から直接列挙（実装非参照の独立表）。
        let cases = [
            (Visibility::FullyVisible, true),
            (Visibility::PartialVisible, true),
            (Visibility::SunriseEclipse, true),
            (Visibility::SunsetEclipse, true),
            (Visibility::NotVisible, false),
            (Visibility::BelowHorizon, false),
        ];
        for (v, expected) in cases {
            assert_eq!(
                super::next_visible_is_observable(v),
                expected,
                "{v:?} の観測可能判定は {expected} であるべき"
            );
        }
    }

    // ------------------------------------------------------------------
    // 8b. next_visible_eclipse 統合 happy-path（SLOW・1件・central_observer）
    //
    // コスト: search（日食 1 件の全球解 ≈ 300s）を 1 回想定（after を日食直前に置き
    //   最初の探索窓で皆既可視日食が見つかる）。red 段階では unimplemented! で即 panic ゆえ
    //   実走しない。将来 green で過度に遅い/不安定なら #[ignore] 付与を検討する。
    // ------------------------------------------------------------------

    /// **happy-path（central_observer・2017-08-21 皆既可視）**: `after = 2017-08-01`、中心食地点
    /// （37.5°N,−89.2°E）で呼ぶと 2017-08-21 の皆既日食を `Some(VisibleSolarEclipse)` で返す。
    /// - `Ok(Some(vse))`（可視日食が見つかる）。
    /// - `vse.eclipse.event_key` が `"2017-08-21"` で始まる（最初の可視日食＝2017 北米皆既, NASA 事実）。
    /// - `vse.local.visibility` は「見える種別」（`next_visible_is_observable == true`）。
    /// - `vse.local.magnitude.0 > 1.0`（中心食地点の皆既＝食分>1, NASA 事実）。
    /// - `vse.local` は `local_circumstances(&vse.eclipse, central_observer())` と整合
    ///   （同じ可視性・食分。再計算は追加で遅いので visibility と magnitude の一致のみで縛る）。
    ///
    /// 殺す変異: 可視性を見ずに最初の日食をそのまま返す（happy では event_key 一致で漏れるが
    /// 8c の skip で撃破）、`None` を返す（Some を要求）、`local` を別観測者/別日食で計算する
    /// （visibility/magnitude 整合で撃破）、返す local.visibility が「見える種別」でない。
    #[test]
    fn next_visible_eclipse_central_site_returns_2017_total() {
        let engine = standard_engine_from_synthetic();
        let after = utc(2017, 8, 1, 0, 0, 0.0);
        let vse = engine
            .next_visible_eclipse(after, central_observer())
            .expect("next_visible_eclipse はエラーにならない")
            .expect("中心食地点では 2017-08-21 皆既が可視日食として見つかる");

        // 最初の可視日食は 2017-08-21（北米皆既・NASA 事実）。
        assert!(
            vse.eclipse.event_key.starts_with("2017-08-21"),
            "最初の可視日食は 2017-08-21: event_key = {}",
            vse.eclipse.event_key
        );

        // 返る local.visibility は「見える種別」（採否判定と整合）。
        assert!(
            super::next_visible_is_observable(vse.local.visibility),
            "返る日食の可視性は見える種別: {:?}",
            vse.local.visibility
        );

        // 中心食地点ゆえ皆既（食分 > 1, NASA 事実）。
        assert!(
            vse.local.magnitude.0 > 1.0,
            "中心食地点は皆既なので食分 > 1: {}",
            vse.local.magnitude.0
        );

        // local は local_circumstances(&eclipse, central_observer()) と整合（同一観測者・日食の局地条件）。
        let recomputed = engine
            .local_circumstances(&vse.eclipse, central_observer())
            .expect("同一日食・同一観測者の局地条件は成功する");
        assert_eq!(
            vse.local.visibility, recomputed.visibility,
            "vse.local は local_circumstances と同じ可視性"
        );
        assert_eq!(
            vse.local.magnitude, recomputed.magnitude,
            "vse.local は local_circumstances と同じ食分"
        );
    }

    // ------------------------------------------------------------------
    // 8c. next_visible_eclipse 統合 skip（SLOW・1件・invisible_observer）
    //
    // コスト警告: invisible_observer（−40°S,140°E）の次の可視日食まで複数件の全球解を要する
    //   可能性があり非常に遅い（数百〜千秒超）。red 段階では unimplemented! で即 panic ゆえ
    //   実走しない。将来 green で過度に遅い/不安定なら #[ignore] 付与を検討する（要コメント）。
    //   ただし skip ロジックの唯一のガードゆえ設計としては入れておく。
    // ------------------------------------------------------------------

    /// **skip（invisible_observer・2017-08-21 は不可視ゆえスキップ）**: `after = 2017-08-01`、
    /// 非可視地点（−40°S,140°E。2017-08-21 は NotVisible=既存テストで確認済）で呼ぶと、
    /// 2017-08-21 を**スキップ**して後続の可視日食を返す。
    /// - `Ok(Some(vse))`（後続のどこかに可視日食がある）。
    /// - `vse.eclipse.event_key` が `"2017-08-21"` で**始まらない**（不可視日食をスキップした証拠）。
    /// - `vse.local.visibility` は「見える種別」（`next_visible_is_observable == true`）。
    ///
    /// skip 先の具体日付は geography 依存で flaky ゆえハードコードしない（「2017-08-21 でない可視
    /// 日食」という構造で縛る）。
    ///
    /// 殺す変異（このテスト固有の主目的）: 可視性を評価せず単に最初の日食（2017-08-21）を返す
    /// （event_key が 2017-08-21 で始まらない＝スキップを縛ることで撃破）、不可視日食を返す
    /// （見える種別を縛ることで撃破）。
    #[test]
    fn next_visible_eclipse_skips_invisible_eclipse() {
        let engine = standard_engine_from_synthetic();
        let after = utc(2017, 8, 1, 0, 0, 0.0);
        let vse = engine
            .next_visible_eclipse(after, invisible_observer())
            .expect("next_visible_eclipse はエラーにならない")
            .expect("非可視地点でも後続に可視日食が存在する");

        // 2017-08-21 はこの観測者では不可視ゆえスキップされる（最初の日食をそのまま返さない）。
        assert!(
            !vse.eclipse.event_key.starts_with("2017-08-21"),
            "不可視の 2017-08-21 はスキップされる: event_key = {}",
            vse.eclipse.event_key
        );

        // 返る日食はこの観測者で「見える種別」。
        assert!(
            super::next_visible_is_observable(vse.local.visibility),
            "スキップ後に返る日食は見える種別: {:?}",
            vse.local.visibility
        );
    }

    // ------------------------------------------------------------------
    // 8d. None ケースは意図的にテストしない
    //
    // horizon 内に可視日食が無いケースは全 horizon 走査が必要で極めて遅いため（数百秒×多数件）、
    // ここではテストしない。`Ok(None)`（該当なしはエラーにしない）の契約は型シグネチャ
    // （`Result<Option<..>, ..>`）と docstring で表現し、happy/skip の `Some` 経路で
    // 「見つかれば Some」を縛るに留める。
    // ------------------------------------------------------------------

    // ==================================================================
    // 9. sample_central_point: 独立オラクル（中心点 + 南北限界）
    //
    // オラクル戦略（strict・追認回避）: `sample_central_point` の各算術
    //   （umbral_radius=|l2−ζ₀·tanf2|, t_hours=days_since(epoch)·24, 速度ベクトル (vx,vy)・
    //    速さ・単位法線 (−vy,vx)/|v|, 縁基本面点 (x±|L2'|·n), 南北＝高緯度側）
    // を、bessel.at(t) が返した瞬時要素から**テスト内で独立に再計算**して突き合わせる。
    // 地表射影だけは別テストで縛り済みの `surface_point_for_fundamental` を再利用してよい
    // （test A で独立にピン済み）。x,y を**二次**多項式にすることで影速度が t_hours に依存し、
    // 時間スケール（·24 や days_since の方向）の取り違えが (vx,vy) 経由で南北点へ伝播し検出される。
    // ------------------------------------------------------------------

    /// 二次 x,y を持つ合成ベッセル多項式。epoch_tt をフィット中心に置き、t_hours が綺麗に出るようにする。
    /// l2<0（皆既）・tan_f2>0 で umbral_radius=|l2−ζ₀·tanf2| は綺麗な正値。x,y は単位円内に十分収まる
    /// 小さな値（サンプル時 x≈0.1, y≈0.05 付近）で中心・両縁とも地表に当たる。
    fn quadratic_bessel(epoch: TtInstant) -> crate::bessel_poly::BesselianPolynomial {
        use crate::polynomial::Polynomial;
        let p = |coeffs: Vec<f64>| Polynomial {
            coefficients: coeffs,
        };
        crate::bessel_poly::BesselianPolynomial {
            epoch_tt: epoch,
            // x(t)=0.08 + 0.03 t + 0.01 t²（t は epoch からの hour）— 速度 x'(t)=0.03+0.02 t は t 依存。
            x: p(vec![0.08, 0.03, 0.01]),
            // y(t)=0.04 + 0.02 t − 0.005 t²  — 速度 y'(t)=0.02−0.01 t は t 依存。
            y: p(vec![0.04, 0.02, -0.005]),
            d: p(vec![0.2070]),
            mu: p(vec![1.2]),
            l1: p(vec![0.5400]),
            l2: p(vec![-0.0090]),
            tan_f1: 0.004_65,
            tan_f2: 0.004_63,
            // fit_interval はサンプル時刻（epoch + 0.5h ≈ epoch + 0.0208 d）を十分含む幅。
            fit_interval: umbra_core::TimeInterval {
                start: tt(2_451_544.5, 0.0),
                end: tt(2_451_545.5, 0.0),
            },
            fit_error: crate::bessel_poly::BesselFitError {
                max_x: 1.0e-7,
                max_y: 2.0e-7,
                max_l1: 3.0e-7,
                max_l2: 4.0e-7,
            },
        }
    }

    /// include_limits=true: 中心点・北限・南限の 3 点が、bessel.at(t) の要素から独立に組んだ
    /// umbral_radius・単位法線・縁基本面点（surface_point_for_fundamental 射影）と一致する。
    /// t_hours≠0（epoch+0.5h）で速度が t 依存。各座標値を絶対値で突き合わせ、
    /// umbral_radius・法線・時間スケール・縁オフセットの符号/演算子変異を撃破する。
    #[test]
    fn sample_central_point_independent_oracle_center_and_limits() {
        // surface_point_for_fundamental は `use super::*` 経由でスコープ内（engine.rs L29 で import 済み）。
        let epoch = tt(2_451_545.0, 0.0);
        let bessel = quadratic_bessel(epoch);
        let x_deriv = bessel.x.derivative();
        let y_deriv = bessel.y.derivative();
        // サンプル時刻 = epoch + 0.5 hour（t_hours≠0 を保証）。
        let t = TtInstant::from_jd2(epoch.jd2().add_days(0.5 / 24.0));
        let ellipsoid = Ellipsoid::WGS84;

        // 関数戻り（被テスト）。
        let (center, limits) =
            sample_central_point(&bessel, &x_deriv, &y_deriv, epoch, t, true, &ellipsoid)
                .expect("中心・両縁とも地表に当たる構成では Ok")
                .expect("speed≠0 かつ縁が外れないので Some");
        let (north, south) = limits.expect("include_limits=true なら Some((north,south))");

        // ---- 独立オラクル ----
        // bessel.at(t) は別途 bessel_poly.rs で検証済み。ここでは要素を真値として受ける。
        let elements = bessel.at(t).expect("区間内サンプルは評価成功");

        // 中心点: x,y を直接射影（test A で縛った関数を再利用）。
        let (exp_center, zeta0) = surface_point_for_fundamental(
            elements.x,
            elements.y,
            elements.declination,
            elements.mu,
            &ellipsoid,
        )
        .expect("中心軸は地表に当たる");
        assert!(
            (center.lat.degrees().0 - exp_center.lat.degrees().0).abs() < 1e-7
                && (center.lon.degrees().0 - exp_center.lon.degrees().0).abs() < 1e-7,
            "center {:?} expected {:?}（中心軸射影）",
            center,
            exp_center
        );

        // umbral_radius=|l2−ζ₀·tan_f2|（独立算術）。
        let umbral_radius = (elements.l2 - zeta0 * elements.tan_f2).abs();
        // t_hours と速度ベクトル（独立: epoch から days_since·24）。
        let t_hours = t.jd2().days_since(epoch.jd2()) * 24.0;
        assert!(
            t_hours.abs() > 1e-9,
            "t_hours={t_hours} は非零（速度の t 依存を観測）"
        );
        let vx = x_deriv.eval(t_hours);
        let vy = y_deriv.eval(t_hours);
        let speed = vx.hypot(vy);
        assert!(speed > 0.0, "速度ゼロでない構成");
        // 単位法線 n=(−vy, vx)/|v|（独立算術）。
        let nx = -vy / speed;
        let ny = vx / speed;

        // 両縁の基本面点 (x ± |L2'|·n) を独立に組み、射影する。
        let project_edge = |sign: f64| {
            let xi = elements.x + sign * umbral_radius * nx;
            let eta = elements.y + sign * umbral_radius * ny;
            surface_point_for_fundamental(xi, eta, elements.declination, elements.mu, &ellipsoid)
                .expect("縁は地表に当たる構成")
                .0
        };
        let edge_plus = project_edge(1.0);
        let edge_minus = project_edge(-1.0);
        // 高緯度側が north。
        let (exp_north, exp_south) = if edge_plus.lat.degrees().0 >= edge_minus.lat.degrees().0 {
            (edge_plus, edge_minus)
        } else {
            (edge_minus, edge_plus)
        };

        // 北限・南限の絶対値一致（umbral_radius・法線・時間スケール・縁オフセットを縛る）。
        assert!(
            (north.lat.degrees().0 - exp_north.lat.degrees().0).abs() < 1e-7
                && (north.lon.degrees().0 - exp_north.lon.degrees().0).abs() < 1e-7,
            "north {:?} expected {:?}",
            north,
            exp_north
        );
        assert!(
            (south.lat.degrees().0 - exp_south.lat.degrees().0).abs() < 1e-7
                && (south.lon.degrees().0 - exp_south.lon.degrees().0).abs() < 1e-7,
            "south {:?} expected {:?}",
            south,
            exp_south
        );
        // north が確かに高緯度側（南北判定の取り違えを縛る補助確認）。
        assert!(
            north.lat.degrees().0 >= south.lat.degrees().0,
            "north.lat={} ≥ south.lat={} のはず",
            north.lat.degrees().0,
            south.lat.degrees().0
        );
    }

    /// include_limits=false: 中心点のみ・limits は None（早期 return 経路を縛る）。
    #[test]
    fn sample_central_point_no_limits_returns_none() {
        let epoch = tt(2_451_545.0, 0.0);
        let bessel = quadratic_bessel(epoch);
        let x_deriv = bessel.x.derivative();
        let y_deriv = bessel.y.derivative();
        let t = TtInstant::from_jd2(epoch.jd2().add_days(0.5 / 24.0));
        let ellipsoid = Ellipsoid::WGS84;

        let (center, limits) =
            sample_central_point(&bessel, &x_deriv, &y_deriv, epoch, t, false, &ellipsoid)
                .expect("中心軸は地表に当たる")
                .expect("speed/縁を評価せず Some を返す");
        assert!(limits.is_none(), "include_limits=false なら limits は None");
        // 中心点は include_limits の有無に依らないので、true ケースと同一であることまでは
        // ここで再検証しない（中心の絶対値は上のテストで縛り済み）。サニティのみ。
        assert!(
            center.lat.degrees().0.is_finite() && center.lon.degrees().0.is_finite(),
            "中心点は有限"
        );
    }
}
