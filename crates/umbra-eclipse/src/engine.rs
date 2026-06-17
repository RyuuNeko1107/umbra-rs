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

use umbra_core::constants::EARTH_EQUATORIAL_RADIUS_M;
use umbra_core::deltat::DeltaTModel;
use umbra_core::eop::EarthOrientation;
use umbra_core::{
    EspenakMeeusDeltaT, IersEopData, TimeData, TimeInterval, TimeRange, TimeScales, TtInstant,
    UtcInstant,
};
use umbra_ephemeris::{AnalyticalEphemeris, AstrometryOptions, Ephemeris};

use crate::besselian::InstantaneousBesselianElements;
use crate::config::EngineConfig;
use crate::error::EclipseError;
use crate::path::{EclipsePath, PathOptions};
use crate::results::SolarEclipse;
use crate::source::{BesselianSource, InstantaneousEvaluator};

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

    /// v0.1 未実装（経路は umbra-geo・M9）。panic でなく `Err(NotImplemented)`（PATH/045）。
    pub fn path(
        &self,
        _eclipse: &SolarEclipse,
        _options: PathOptions,
    ) -> Result<EclipsePath, EclipseError> {
        Err(EclipseError::NotImplemented)
    }
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

    use crate::config::{EngineConfig, LunarRadiusModel, SolarRadiusModel};
    use crate::error::EclipseError;
    use crate::global::SolarEclipseKind;
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
    // 3. path(): 常に Err(NotImplemented)
    // ==================================================================

    /// `path()` は引数を使わず常に `Err(EclipseError::NotImplemented)` を返す（panic でない）。
    /// 最小 SolarEclipse と既定 PathOptions で呼ぶ。
    ///
    /// 殺す変異: NotImplemented 以外の variant への差し替え・Ok を返す・panic 化（unimplemented! 等）。
    #[test]
    fn path_returns_not_implemented() {
        let engine = standard_engine_from_synthetic();
        let eclipse = minimal_eclipse();
        let r = engine.path(&eclipse, PathOptions::default());
        assert!(
            matches!(r, Err(EclipseError::NotImplemented)),
            "expected Err(EclipseError::NotImplemented), got {r:?}"
        );
    }

    /// `path()` は PathOptions の値に依らず常に `Err(NotImplemented)`（引数を使わない確定仕様）。
    /// 既定とは異なる PathOptions でも同じ結果になることで「引数を見て分岐しない」ことを縛る。
    ///
    /// 殺す変異: PathOptions/eclipse の中身で戻り値を分岐させる。
    #[test]
    fn path_ignores_options_and_always_not_implemented() {
        let engine = standard_engine_from_synthetic();
        let eclipse = minimal_eclipse();
        let custom = PathOptions {
            sample_interval_seconds: 5.0,
            include_limits: false,
            split_antimeridian: false,
        };
        let r = engine.path(&eclipse, custom);
        assert!(
            matches!(r, Err(EclipseError::NotImplemented)),
            "非既定 options でも NotImplemented, got {r:?}"
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
}
