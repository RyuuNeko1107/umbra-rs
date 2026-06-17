//! 全球の日食種別判定（`docs/issues/ISSUE-023`、`docs/physical-models.md` §C3）。
//!
//! 影軸の地心最小距離 `gamma`（Re）と本影半径 `l2`（符号付き）で分類する（Meeus Ch.54 基準）:
//! ```text
//! |gamma| < 0.9972                    → 中心食（l2<0 皆既 / l2>0 金環）
//! 0.9972 ≤ |gamma| < 0.9972 + |l2|    → 非中心 皆既/金環
//! 0.9972 + |l2| ≤ |gamma| < 1.5433+l2 → 部分食
//! |gamma| ≥ 1.5433 + l2               → 日食なし
//! ```
//! 注: ハイブリッド（中心線上で l2 が符号反転）は単一時刻では判定不能。全球パス（時系列）で
//! 判定し本関数は瞬時の Total/Annular を返す（要確認: 0.9972/1.5433 の式番号・有効桁。§C3）。

use crate::axis_intercept::shadow_axis_surface_point;
use crate::besselian::BesselianElements;
use crate::config::EngineConfig;
use crate::conjunction::RootConfig;
use crate::error::EclipseError;
use crate::horizontal::{sun_horizontal, RefractionModel};
use crate::local_maximum::solve_local_maximum;
use crate::magnitude::{eclipse_magnitude, eclipse_obscuration};
use crate::projection::project_observer_to_fundamental;
use crate::results::GreatestEclipse;
use crate::source::BesselianSource;
use umbra_core::deltat::DeltaTModel;
use umbra_core::ellipsoid::{observer_geocentric, Ellipsoid, GeocentricObserver};
use umbra_core::Radians;

/// 1 日 = 86400 SI 秒（root_tolerance を日へ換算）。
const SECONDS_PER_DAY: f64 = 86_400.0;
/// 最大食時刻 Brent 求根の反復上限。
const GREATEST_ROOT_MAX_ITER: usize = 200;

/// 中心食境界（≈1 − 扁平縮約。Meeus）。
const CENTRAL_LIMIT: f64 = 0.9972;
/// 半影限界（Meeus）。
const PENUMBRA_LIMIT: f64 = 1.5433;

/// 太陽食の種別。
#[non_exhaustive]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SolarEclipseKind {
    /// 部分食。
    Partial,
    /// 金環食。
    Annular,
    /// 皆既食。
    Total,
    /// ハイブリッド（金環↔皆既。全球パスで判定）。
    Hybrid,
    /// 非中心の皆既。
    NonCentralTotal,
    /// 非中心の金環。
    NonCentralAnnular,
}

/// 瞬時ベッセル要素から日食種別を判定する。`None` は（その時刻に全球で）日食なし。
///
/// ハイブリッドは返さない（時系列が必要。上記注）。中心/非中心は l2 符号で皆既/金環を分ける。
pub fn classify(elements: &BesselianElements) -> Option<SolarEclipseKind> {
    let g = elements.gamma(); // ≥ 0
    let l2 = elements.l2;
    let total = l2 < 0.0; // l2<0 皆既 / l2>0 金環（正本 B1）

    if g < CENTRAL_LIMIT {
        Some(if total {
            SolarEclipseKind::Total
        } else {
            SolarEclipseKind::Annular
        })
    } else if g < CENTRAL_LIMIT + l2.abs() {
        Some(if total {
            SolarEclipseKind::NonCentralTotal
        } else {
            SolarEclipseKind::NonCentralAnnular
        })
    } else if g < PENUMBRA_LIMIT + l2 {
        Some(SolarEclipseKind::Partial)
    } else {
        None
    }
}

/// 最大食（global greatest eclipse）の解（時刻・地表点・食分/食面積/太陽高度）と gamma。
///
/// `kind` と全球接触 P1/U1/U4/P4・帯幅・中心食継続は S6b、`SolarEclipse` 組立は S6c の責務。
/// `path_width`/`central_duration` は本スライスでは常に `None`（S6b で充足）。
#[derive(Clone, Copy, Debug, PartialEq)]
#[allow(dead_code)] // S6c（classify_global / search 結線）が消費するまで未使用。
pub(crate) struct GreatestEclipseSolution {
    /// 最大食の地表点・食分・食面積・太陽高度・時刻（path/duration は None）。
    pub greatest: GreatestEclipse,
    /// 影軸の地心最小距離 gamma（Re, 符号なし）。
    pub gamma: f64,
}

/// 中心食の最大食状況（時刻・gamma・地表点・食分/食面積・太陽高度）を解く（ISSUE-043 S6a-ii）。
///
/// 検証済みプリミティブの合成: 最大食時刻と gamma は地心観測者（ρ=0）に対する
/// [`solve_local_maximum`] で得る（m²=(0−x)²+(0−y)²=x²+y²=gamma², ISSUE-026）。地表点は
/// [`shadow_axis_surface_point`]（S6a-i）、太陽高度は [`sun_horizontal`]（ISSUE-028）。食分・食面積は
/// 観測者 ζ で補正した本影/半影半径 L1'=l1−ζ·tanf1, L2'=l2−ζ·tanf2 から
/// [`eclipse_magnitude`]/[`eclipse_obscuration`] で評価する（中心点では m=0・視半径比 ρ=(L1'−L2')/(L1'+L2')）。
///
/// 中心食でない（軸が地表を外す）と地表点が無く `Err(Solver(RootNotBracketed))`（S6b で部分/非中心を扱う）。
#[allow(dead_code)] // S6c が消費するまで未使用。
pub(crate) fn solve_greatest_eclipse<B, M>(
    source: &B,
    delta_t: &M,
    config: &EngineConfig,
) -> Result<GreatestEclipseSolution, EclipseError>
where
    B: BesselianSource,
    M: DeltaTModel,
{
    // earth_model は現状 WGS84 のみ（projection/horizontal/axis_intercept と同様に定数で扱う）。
    let _ = config.earth_model;
    let ellipsoid = Ellipsoid::WGS84;

    // 1. 最大食時刻・gamma: 地心観測者（ρsin=ρcos=0）の局地最大食で得る。射影は ξ=η=ζ=0 となり
    //    m² = (0−x)² + (0−y)² = x² + y² = gamma²（影軸の地心距離²）。経度は ρcos=0 ゆえ無関係。
    let geocenter = GeocentricObserver {
        rho_sin_phi_prime: 0.0,
        rho_cos_phi_prime: 0.0,
    };
    let root_config = RootConfig {
        x_tolerance_days: config.root_tolerance_seconds / SECONDS_PER_DAY,
        max_iterations: GREATEST_ROOT_MAX_ITER,
    };
    let max = solve_local_maximum(
        source,
        &geocenter,
        Radians::new(0.0),
        source.fit_interval(),
        root_config,
    )?;
    let gamma = max.min_separation; // = √(x²+y²) at t_max

    // 2. 最大食時刻の瞬時ベッセル要素。
    let elements = source.at(max.time_tt)?;

    // 3. 影軸の地表貫通点（中心食でなければ Err(Solver(RootNotBracketed))・S6b で部分/非中心）。
    let position = shadow_axis_surface_point(&elements, &ellipsoid)?;

    // 4. 地表点の観測者 ζ（検証済み前方射影し直して取得。中心点なので ξ=x, η=y, m=0）。
    let phi = position.lat.radians();
    let lambda = position.lon.radians();
    let obs = observer_geocentric(&ellipsoid, phi.0, 0.0);
    let zeta = project_observer_to_fundamental(&obs, lambda, &elements).zeta;

    // 5. 食分・食面積（中心点 m=0・観測者 ζ で補正した半径）。
    //    L1'=l1−ζ·tanf1（半影）, L2'=l2−ζ·tanf2（本影, 符号付き）。視半径比 ρ=(L1'−L2')/(L1'+L2')。
    let l1p = elements.l1 - zeta * elements.tan_f1;
    let l2p = elements.l2 - zeta * elements.tan_f2;
    let magnitude = eclipse_magnitude(0.0, l1p, l2p);
    let radius_ratio = (l1p - l2p) / (l1p + l2p);
    let obscuration = eclipse_obscuration(0.0, 1.0, radius_ratio);

    // 6. 太陽の幾何学的高度（大気差なし, conventions §7 既定）。
    let sun_altitude =
        sun_horizontal(phi, lambda, max.time_tt, RefractionModel::None, delta_t).altitude_geometric;

    let greatest = GreatestEclipse {
        time_utc: max.time_utc,
        time_tt: max.time_tt,
        position,
        magnitude,
        obscuration,
        path_width: None,
        central_duration: None,
        sun_altitude,
    };
    Ok(GreatestEclipseSolution { greatest, gamma })
}

#[cfg(test)]
mod tests {
    use super::SolarEclipseKind::{Annular, NonCentralAnnular, NonCentralTotal, Partial, Total};
    use super::*;
    use umbra_core::constants::{EARTH_EQUATORIAL_RADIUS_M, SOLAR_RADIUS_KM};
    use umbra_core::{JulianDate2, Radians, TdbInstant};
    use umbra_ephemeris::{Body, Ephemeris, EphemerisFrame, MockEphemeris, Origin};

    /// gamma=`g`（x=g, y=0）, 本影半径 `l2` のベッセル要素を作る。
    fn elem(g: f64, l2: f64) -> BesselianElements {
        BesselianElements {
            x: g,
            y: 0.0,
            declination: Radians(0.0),
            l1: 0.53,
            l2,
            tan_f1: 0.0047,
            tan_f2: 0.0046,
        }
    }

    #[test]
    fn central_total_and_annular_by_l2_sign() {
        assert_eq!(classify(&elem(0.5, -0.02)), Some(Total));
        assert_eq!(classify(&elem(0.5, 0.02)), Some(Annular));
    }

    #[test]
    fn non_central_band_by_l2_sign() {
        // 0.9972 ≤ g < 0.9972+|l2|（|l2|=0.02 → 上限 1.0172）。
        assert_eq!(classify(&elem(1.0, -0.02)), Some(NonCentralTotal));
        assert_eq!(classify(&elem(1.0, 0.02)), Some(NonCentralAnnular));
    }

    #[test]
    fn partial_band() {
        assert_eq!(classify(&elem(1.2, 0.01)), Some(Partial));
    }

    #[test]
    fn no_eclipse_when_gamma_too_large() {
        assert_eq!(classify(&elem(2.0, 0.01)), None);
    }

    #[test]
    fn central_to_noncentral_boundary() {
        // g=0.9972 ちょうどは中心食でない（< 厳密）→ 非中心。直下は中心。
        assert_eq!(classify(&elem(0.9972, -0.02)), Some(NonCentralTotal));
        assert_eq!(classify(&elem(0.9971, -0.02)), Some(Total));
    }

    #[test]
    fn noncentral_to_partial_boundary() {
        // 境界は実装と同じ計算 CENTRAL_LIMIT+|l2| で踏む（リテラル 1.0172 では f64 が
        // ぴったり一致せず < / <= を区別できない）。境界ちょうどは非中心でない → 部分。
        let b = CENTRAL_LIMIT + 0.02;
        assert_eq!(classify(&elem(b, -0.02)), Some(Partial));
        assert_eq!(classify(&elem(b - 1e-6, -0.02)), Some(NonCentralTotal));
    }

    #[test]
    fn l2_exactly_zero_is_annular_not_total() {
        // l2==0（皆既/金環の連続境界）は total=(l2<0)=false → 金環側に倒す（< 厳密）。
        assert_eq!(classify(&elem(0.5, 0.0)), Some(Annular));
    }

    #[test]
    fn partial_to_none_boundary() {
        // g=1.5433+l2 ちょうどは日食なし → 直下は部分。
        assert_eq!(classify(&elem(1.5433 + 0.01, 0.01)), None);
        assert_eq!(classify(&elem(1.5433 + 0.01 - 1e-6, 0.01)), Some(Partial));
    }

    #[test]
    fn partial_to_none_boundary_negative_l2() {
        // 皆既側(l2<0)では上限が 1.5433+l2 < 1.5433 に縮む（符号付き）。
        // `+l2` を `+|l2|` や `+0.01` に取り違えるとこのテストで露見する（H1）。
        let l2 = -0.02;
        assert_eq!(classify(&elem(1.5433 + l2 - 1e-6, l2)), Some(Partial));
        assert_eq!(classify(&elem(1.5433 + l2 + 1e-6, l2)), None);
    }

    #[test]
    fn partial_and_none_bands_with_negative_l2() {
        assert_eq!(classify(&elem(1.2, -0.02)), Some(Partial));
        assert_eq!(classify(&elem(2.0, -0.02)), None);
    }

    #[test]
    fn matches_mock_configurations() {
        let t = TdbInstant::from_jd2(JulianDate2::from_jd(2_451_545.0));
        let r_sun = SOLAR_RADIUS_KM;
        let r_moon = 0.2725076 * (EARTH_EQUATORIAL_RADIUS_M / 1000.0);
        let kind = |m: &MockEphemeris| {
            let pos = |b| {
                m.state(b, t, Origin::Geocenter, EphemerisFrame::Icrs)
                    .unwrap()
                    .position
            };
            let e = crate::besselian::besselian_elements(
                pos(Body::Sun),
                pos(Body::Moon),
                r_sun,
                r_moon,
            )
            .unwrap();
            classify(&e)
        };
        assert_eq!(kind(&MockEphemeris::central_total()), Some(Total));
        assert_eq!(kind(&MockEphemeris::clear_annular()), Some(Annular));
        assert_eq!(kind(&MockEphemeris::clear_partial()), Some(Partial));
        assert_eq!(kind(&MockEphemeris::shadow_misses_earth()), None);
        // 非中心皆既（暦→ベッセル→分類の貫通で NonCentralTotal バンドを踏む, M2）。
        assert_eq!(
            kind(&MockEphemeris::non_central_total()),
            Some(NonCentralTotal)
        );
    }

    // ====================================================================
    // solve_greatest_eclipse（ISSUE-043 S6a-ii・全球最大食組立）
    // ====================================================================
    //
    // ## オラクル戦略（追認回避）
    // 主オラクルは **独立内部再計算**。solver の内部手法（local_maximum・逆射影・magnitude/
    // obscuration の合成）には依存せず、外部観測可能な振る舞いを **検証済プリミティブ**
    // `source.at(t)`（ISSUE-037）＋ `project_observer_to_fundamental`（ISSUE-024）＋
    // `observer_geocentric`（ISSUE-010/011）から **別経路** で縛る:
    //   - gamma 再計算: 返った time_tt で √(x²+y²) を `source.at` から直接組む（1e-6 一致）。
    //   - 局地最小性: time_tt±δ の gamma が time_tt 以上（局地最小）。300s 両側で狭義増加。
    //   - 地表点往復: 返った position を前方射影し直すと (ξ,η)=(e.x,e.y)・ζ>0（逆射影の独立検証）。
    //   - 食面積==1（皆既の強い縛り）: 中心点では太陽が完全に隠れる。
    // NASA 公表値（gamma≈0.4367・greatest≈18:25UTC・37.0N/87.7W・magnitude≈1.031・alt 61–64°）は
    // **ballpark のみ**（k/ΔT 慣習差で秒値は再現しない）。range check に限定して用いる。
    //
    // 探索窓は `source.fit_interval()`。2017-08-21 最大食（TT-JD≈2457986.768）を内部に括る
    // 区間 [2457986, 2457988] を渡す（時刻 solver がブラケットできるよう最小が窓内部に来る）。
    //
    // 注（追補・実装レビューの結線網羅）: 当初省略した 2 経路を後段で追加した（実装が確定し
    // 金環分岐・非中心分岐が end-to-end で踏めるようになったため）:
    //   - 金環（2023-10-14）: l2>0 ⇒ magnitude<1 / obscuration<1 / radius_ratio<1 の分岐を踏む
    //     （`greatest_annular_2023_*`）。NASA 値は ballpark のみ・solver が窓内の真の最小を探す。
    //   - 軸が地表を外す→RootNotBracketed: 時不変 `ConstantSource` は時不変ゆえ「軸ミス」ではなく
    //     `solve_local_maximum` のブラケット不成立で別理由 Err になる。代わりに **時変** 合成供給源
    //     （gamma に内部極小を持つが極小値が >1）を使い、`solve_local_maximum` は成功し
    //     `shadow_axis_surface_point` 側で RootNotBracketed が起きることを縛る（`greatest_*_axis_miss_*`）。

    use umbra_core::ellipsoid::{observer_geocentric, Ellipsoid};
    use umbra_core::{EspenakMeeusDeltaT, SolverError, TimeInterval, TtInstant};

    use crate::besselian::InstantaneousBesselianElements;
    use crate::projection::project_observer_to_fundamental;
    use crate::source::{BesselianSource, DirectBesselianSource};

    /// 1 日 = 86400 SI 秒。
    const SECONDS_PER_DAY: f64 = 86_400.0;
    /// 太陽物理半径[km]（local_maximum.rs / axis_intercept.rs と同一）。
    const G_R_SUN: f64 = SOLAR_RADIUS_KM;
    /// 月半径[km]（k·Re, IAU 慣習 k=0.2725076・同上）。
    const G_R_MOON: f64 = 0.2725076 * (EARTH_EQUATORIAL_RADIUS_M / 1000.0);

    /// 2017-08-21 最大食を内部に括る探索窓（TT-JD 2457986〜2457988, axis_intercept.rs と同形）。
    /// 最大食 TT-JD≈2457986.768 が区間内部にあり、時刻 solver が最小をブラケットできる。
    fn solve_window_2017() -> TimeInterval<TtInstant> {
        TimeInterval {
            start: TtInstant::from_jd2(JulianDate2::new(2_457_986.0, 0.0)),
            end: TtInstant::from_jd2(JulianDate2::new(2_457_988.0, 0.0)),
        }
    }

    /// この供給源の TT-JD（単一 f64）から TtInstant。
    fn g_tt_jd(jd: f64) -> TtInstant {
        TtInstant::from_jd2(JulianDate2::from_jd(jd))
    }

    /// 2017-08-21 中心皆既を括る `DirectBesselianSource`（fit_interval=探索窓）を解く。
    fn solve_2017<'d>(
        dt: &'d EspenakMeeusDeltaT,
    ) -> (
        DirectBesselianSource<'d, EspenakMeeusDeltaT>,
        GreatestEclipseSolution,
    ) {
        let src = DirectBesselianSource::new(G_R_SUN, G_R_MOON, dt, solve_window_2017());
        let config = crate::config::EngineConfig::standard();
        let sol = solve_greatest_eclipse(&src, dt, &config)
            .expect("2017 central total eclipse should yield a greatest-eclipse solution");
        (src, sol)
    }

    /// 影軸交点 (x,y) から geocentric な gamma=√(x²+y²) を `source.at` 由来で独立再計算する。
    fn gamma_at<B: BesselianSource>(src: &B, t: TtInstant) -> f64 {
        let e = src
            .at(t)
            .expect("source.at should succeed near 2017 eclipse");
        (e.x * e.x + e.y * e.y).sqrt()
    }

    /// 観点1（gamma ballpark, 緩め・NASA range check）: gamma ∈ [0.40, 0.47]。
    /// NASA 公表 gamma≈0.4367 を k/ΔT 慣習差の余裕込みで括る（絶対基準にしない）。
    #[test]
    fn greatest_gamma_in_nasa_ballpark() {
        let dt = EspenakMeeusDeltaT;
        let (_src, sol) = solve_2017(&dt);
        assert!(
            (0.40..=0.47).contains(&sol.gamma),
            "gamma={} not in NASA ballpark [0.40,0.47] (NASA≈0.4367)",
            sol.gamma
        );
    }

    /// 観点2（gamma 独立再計算, tight）: 返った gamma が、返った time_tt で `source.at` から
    /// 別経路に組んだ √(x²+y²) と 1e-6 Re 一致する。gamma が「その時刻の geocentric 軸距離」で
    /// あることを直接縛る（追認回避の主オラクル）。
    #[test]
    fn greatest_gamma_matches_independent_recomputation() {
        let dt = EspenakMeeusDeltaT;
        let (src, sol) = solve_2017(&dt);
        let g = gamma_at(&src, sol.greatest.time_tt);
        assert!(
            (sol.gamma - g).abs() < 1e-6,
            "gamma={} must match independent √(x²+y²)={g} at time_tt (tol 1e-6 Re)",
            sol.gamma
        );
    }

    /// 観点3（局地最小性, tight・主オラクル）: 返った time_tt で gamma が局地最小である。
    /// time_tt±60s / ±300s の gamma（`source.at` から独立再計算）が time_tt の gamma 以上で、
    /// ±300s では両側で狭義増加する（平坦/最大取り違えを弾く）。
    #[test]
    fn greatest_time_is_local_minimum_of_geocentric_gamma() {
        let dt = EspenakMeeusDeltaT;
        let (src, sol) = solve_2017(&dt);
        let t0 = sol.greatest.time_tt.jd2().jd();
        let g0 = gamma_at(&src, sol.greatest.time_tt);

        for &delta in &[60.0_f64, 300.0] {
            let plus = gamma_at(&src, g_tt_jd(t0 + delta / SECONDS_PER_DAY));
            let minus = gamma_at(&src, g_tt_jd(t0 - delta / SECONDS_PER_DAY));
            assert!(
                plus >= g0,
                "gamma(t+{delta}s)={plus} must be ≥ gamma(t_max)={g0} (local min)"
            );
            assert!(
                minus >= g0,
                "gamma(t−{delta}s)={minus} must be ≥ gamma(t_max)={g0} (local min)"
            );
        }
        // 300s 両側で狭義増加（谷であること＝最小/最大取り違えを弾く）。
        let p300 = gamma_at(&src, g_tt_jd(t0 + 300.0 / SECONDS_PER_DAY));
        let m300 = gamma_at(&src, g_tt_jd(t0 - 300.0 / SECONDS_PER_DAY));
        assert!(
            p300 > g0 && m300 > g0,
            "300s away on both sides must strictly increase: +300s={p300}, −300s={m300}, center={g0}"
        );
    }

    /// 観点4（TT/UTC 整合, tight）: time_utc == tt_to_utc(time_tt)（2017 は post-1972 で変換可能）。
    #[test]
    fn greatest_utc_matches_tt_to_utc() {
        let dt = EspenakMeeusDeltaT;
        let (_src, sol) = solve_2017(&dt);
        let want = umbra_core::time::tt_to_utc(sol.greatest.time_tt)
            .expect("2017 is post-1972, tt_to_utc must succeed");
        assert_eq!(
            sol.greatest.time_utc, want,
            "time_utc must equal tt_to_utc(time_tt)"
        );
    }

    /// 観点5（UTC 時刻 ballpark, 緩め）: time_utc の時刻が 17.5〜18.8 時（北米中緯度の最大食は
    /// 18時台 UTC）。NASA 公表 ≈18:25 UTC を range check で括る。
    #[test]
    fn greatest_utc_hour_in_ballpark() {
        let dt = EspenakMeeusDeltaT;
        let (_src, sol) = solve_2017(&dt);
        let (_y, _mo, _d, hh, mm, _ss) = sol.greatest.time_utc.to_gregorian();
        let hour = f64::from(hh) + f64::from(mm) / 60.0;
        assert!(
            (17.5..=18.8).contains(&hour),
            "greatest UTC hour {hour} not ~18:xx (NASA≈18:25 UTC)"
        );
    }

    /// 観点6（地表点 ballpark, 緩め）: 最大食地点が lat∈[30,42]N・lon∈[−95,−80]E。
    /// NASA 公表 ≈37.0N/87.7W を range check で括る。
    #[test]
    fn greatest_position_in_nasa_ballpark() {
        let dt = EspenakMeeusDeltaT;
        let (_src, sol) = solve_2017(&dt);
        let lat = sol.greatest.position.lat.degrees().0;
        let lon = sol.greatest.position.lon.degrees().0;
        assert!(
            (30.0..=42.0).contains(&lat),
            "greatest lat {lat}° not in ballpark [30,42]N (NASA≈37.0N)"
        );
        assert!(
            (-95.0..=-80.0).contains(&lon),
            "greatest lon {lon}°E not in ballpark [−95,−80] (NASA≈−87.7)"
        );
    }

    /// 観点7（地表点 往復, tight・主オラクル）: 返った position を **検証済み前方射影** へ通すと
    /// (ξ,η)=(e.x,e.y)・ζ>0 に往復一致する（影軸が太陽側で地表を貫く点であることの独立検証）。
    /// 逆射影の内部式は再実装せず、`project_observer_to_fundamental`（ISSUE-024）＋
    /// `observer_geocentric`（ISSUE-010/011）の独立経路で縛る（axis_intercept.rs と同方針）。
    #[test]
    fn greatest_position_roundtrips_through_forward_projection() {
        let dt = EspenakMeeusDeltaT;
        let (src, sol) = solve_2017(&dt);
        let e = src
            .at(sol.greatest.time_tt)
            .expect("source.at should succeed at greatest-eclipse time");

        let phi = sol.greatest.position.lat.radians().0;
        let lam = sol.greatest.position.lon.radians().0;
        let obs = observer_geocentric(&Ellipsoid::WGS84, phi, 0.0);
        let r = project_observer_to_fundamental(&obs, Radians::new(lam), &e);
        assert!(
            (r.xi - e.x).abs() < 1e-9,
            "ξ={} must round-trip to e.x={} (tol 1e-9)",
            r.xi,
            e.x
        );
        assert!(
            (r.eta - e.y).abs() < 1e-9,
            "η={} must round-trip to e.y={} (tol 1e-9)",
            r.eta,
            e.y
        );
        assert!(r.zeta > 0.0, "ζ={} must be sunward (>0)", r.zeta);
    }

    /// 観点8（magnitude, total）: 中心点の食分 > 1（皆既）かつ ballpark < 1.08。
    /// NASA 公表 magnitude≈1.031 を range check で括る（皆既は厳密に 1 超でなければならない）。
    #[test]
    fn greatest_magnitude_exceeds_one_for_total() {
        let dt = EspenakMeeusDeltaT;
        let (_src, sol) = solve_2017(&dt);
        let mag = sol.greatest.magnitude.0;
        assert!(mag > 1.0, "total eclipse magnitude {mag} must exceed 1");
        assert!(
            mag < 1.08,
            "total eclipse magnitude {mag} above ballpark <1.08 (NASA≈1.031)"
        );
    }

    /// 観点9（obscuration==1, total・強い縛り）: 皆既の中心点では太陽面が完全に隠れる。
    /// obscuration ≈ 1.0 を 1e-9 で縛る（中心点で太陽が月に完全内包される幾何の直接帰結）。
    #[test]
    fn greatest_obscuration_is_one_for_total() {
        let dt = EspenakMeeusDeltaT;
        let (_src, sol) = solve_2017(&dt);
        let obsc = sol.greatest.obscuration.0;
        assert!(
            (obsc - 1.0).abs() < 1e-9,
            "total eclipse central obscuration {obsc} must be exactly 1.0 (tol 1e-9)"
        );
    }

    /// 観点9b（ζ補正半径の結線ピン, tight）: 返った magnitude/obscuration が、**検証済み前方射影**
    /// 由来の ζ で **正しい符号** の補正半径 L1'=l1−ζ·tanf1・L2'=l2−ζ·tanf2 から組んだ値に厳密一致する。
    /// ballpark（mag<1.08）・obscuration==1（皆既で ζ 補正に鈍感）では捕捉できない、ζ 補正項の
    /// 符号・演算の取り違え（`l1−ζ·tanf1` → `+`/`/` 等）を撃破する結線ピン（S5b の半径結線ピンと同方針）。
    /// minus 符号は錐頂点が観測者より太陽側（apex ζ > 観測者 ζ）で半径が観測者へ向け縮む幾何に由来
    /// （独立レビューで besselian.rs apex 定義から確認済）。ζ は impl と同じ前方射影法だが本テストは
    /// **正符号を明示** するため、impl が符号を誤れば不一致で落ちる。
    #[test]
    fn greatest_magnitude_obscuration_match_zeta_corrected_radii() {
        let dt = EspenakMeeusDeltaT;
        let (src, sol) = solve_2017(&dt);
        let e = src
            .at(sol.greatest.time_tt)
            .expect("source.at should succeed at greatest-eclipse time");
        // ζ を検証済み前方射影から独立に取得（地表点を射影し直す）。
        let phi = sol.greatest.position.lat.radians().0;
        let lam = sol.greatest.position.lon.radians().0;
        let obs = observer_geocentric(&Ellipsoid::WGS84, phi, 0.0);
        let zeta = project_observer_to_fundamental(&obs, Radians::new(lam), &e).zeta;
        // 正しい符号（minus）の ζ 補正半径。
        let l1p = e.l1 - zeta * e.tan_f1;
        let l2p = e.l2 - zeta * e.tan_f2;
        let want_mag = crate::magnitude::eclipse_magnitude(0.0, l1p, l2p);
        let ratio = (l1p - l2p) / (l1p + l2p);
        let want_obsc = crate::magnitude::eclipse_obscuration(0.0, 1.0, ratio);
        assert_eq!(
            sol.greatest.magnitude, want_mag,
            "magnitude は ζ補正半径(minus)由来でなければならない"
        );
        assert_eq!(
            sol.greatest.obscuration, want_obsc,
            "obscuration は ζ補正半径(minus)由来でなければならない"
        );
    }

    /// 観点10（太陽高度 ballpark, 緩め）: 最大食地点の太陽幾何高度 ∈ [55,70]°。
    /// NASA 公表 ≈61–64° を range check で括る。
    #[test]
    fn greatest_sun_altitude_in_ballpark() {
        let dt = EspenakMeeusDeltaT;
        let (_src, sol) = solve_2017(&dt);
        let alt = sol.greatest.sun_altitude.0;
        assert!(
            (55.0..=70.0).contains(&alt),
            "greatest sun altitude {alt}° not in ballpark [55,70] (NASA≈61–64)"
        );
    }

    /// 観点11（path/duration は本スライス非責務）: path_width・central_duration はともに None。
    /// S6b（帯幅・中心食継続）が充足するまで本関数では常に None を返す契約を縛る。
    #[test]
    fn greatest_path_and_duration_are_none() {
        let dt = EspenakMeeusDeltaT;
        let (_src, sol) = solve_2017(&dt);
        assert_eq!(
            sol.greatest.path_width, None,
            "path_width must be None in this slice (S6b territory)"
        );
        assert_eq!(
            sol.greatest.central_duration, None,
            "central_duration must be None in this slice (S6b territory)"
        );
    }

    // ====================================================================
    // 追補A: 金環食 end-to-end（l2>0 ⇒ magnitude<1 / obscuration<1 / radius_ratio<1 分岐）
    // ====================================================================
    //
    // 2017 は皆既（l2<0・magnitude>1・obscuration==1）のみで、金環分岐
    // （L2'>0 ⇒ magnitude<1, 視半径比 ρ<1）が end-to-end で踏まれていなかった。
    // 2023-10-14 金環日食（greatest≈17:59:28 UTC, NASA ballpark: gamma≈0.375・
    // magnitude≈0.952・obscuration≈0.91・≈11.4°N/83.1°W）を実 epoch で踏む。
    // NASA 公表値は **ballpark のみ**（k/ΔT 慣習差で秒値・小数下位は再現しない）。

    /// 2023-10-14 金環食 greatest（TT-JD≈2_460_232.25）を内部に括る 1.5 日探索窓。
    /// gamma の極小が区間内部にあり、時刻 solver がブラケットできる。
    fn solve_window_2023() -> TimeInterval<TtInstant> {
        TimeInterval {
            start: TtInstant::from_jd2(JulianDate2::new(2_460_231.5, 0.0)),
            end: TtInstant::from_jd2(JulianDate2::new(2_460_233.0, 0.0)),
        }
    }

    /// 2023-10-14 金環食を括る `DirectBesselianSource`（半径・config は 2017 と同一）を解く。
    fn solve_2023<'d>(
        dt: &'d EspenakMeeusDeltaT,
    ) -> (
        DirectBesselianSource<'d, EspenakMeeusDeltaT>,
        GreatestEclipseSolution,
    ) {
        let src = DirectBesselianSource::new(G_R_SUN, G_R_MOON, dt, solve_window_2023());
        let config = crate::config::EngineConfig::standard();
        let sol = solve_greatest_eclipse(&src, dt, &config)
            .expect("2023 annular eclipse should yield a greatest-eclipse solution");
        (src, sol)
    }

    /// 追補A（金環 end-to-end）: 2023-10-14 金環食を解き、金環の契約と独立再計算を縛る。
    ///
    /// **本テストの主眼**は金環分岐（l2>0 ⇒ L2'>0）の網羅: 中心点でも太陽は完全には隠れず
    /// `magnitude<1`・`obscuration<1`・視半径比 ρ=(L1'−L2')/(L1'+L2')<1 を踏む（皆既と区別）。
    /// NASA 値（gamma≈0.375・≈11.4°N/83.1°W・≈17:59 UTC・mag≈0.952・obsc≈0.91）は **ballpark のみ**。
    /// 主オラクルは 2017 と同じ独立再計算（gamma=√(x²+y²)・前方射影往復・TT/UTC 整合）。
    #[test]
    fn greatest_annular_2023_contract_and_independent_checks() {
        let dt = EspenakMeeusDeltaT;
        let (src, sol) = solve_2023(&dt);

        // --- 金環の契約（本テストの主眼・l2>0 / radius_ratio<1 分岐を踏む）---
        let mag = sol.greatest.magnitude.0;
        let obsc = sol.greatest.obscuration.0;
        assert!(
            mag < 1.0,
            "annular central magnitude {mag} must be < 1 (金環は太陽を覆い切らない)"
        );
        assert!(
            mag > 0.85,
            "annular central magnitude {mag} below ballpark >0.85 (NASA≈0.952)"
        );
        assert!(
            obsc < 1.0,
            "annular central obscuration {obsc} must be < 1 (環が残る)"
        );
        assert!(
            obsc > 0.7,
            "annular central obscuration {obsc} below ballpark >0.7 (NASA≈0.91)"
        );

        // --- ballpark（NASA range check のみ・絶対基準にしない）---
        assert!(
            (0.30..=0.45).contains(&sol.gamma),
            "gamma={} not in NASA ballpark [0.30,0.45] (NASA≈0.375)",
            sol.gamma
        );
        let lat = sol.greatest.position.lat.degrees().0;
        let lon = sol.greatest.position.lon.degrees().0;
        assert!(
            (0.0..=22.0).contains(&lat),
            "annular lat {lat}° not in ballpark [0,22]N (NASA≈11.4N)"
        );
        assert!(
            (-92.0..=-74.0).contains(&lon),
            "annular lon {lon}°E not in ballpark [−92,−74] (NASA≈−83.1)"
        );
        let (_y, _mo, _d, hh, mm, _ss) = sol.greatest.time_utc.to_gregorian();
        let hour = f64::from(hh) + f64::from(mm) / 60.0;
        assert!(
            (17.0..=19.0).contains(&hour),
            "annular greatest UTC hour {hour} not ~18:xx (NASA≈17:59 UTC)"
        );

        // --- 独立再計算（tight・2017 と同じ主オラクル）---
        // gamma == √(x²+y²) at time_tt（source.at 由来の別経路, tol 1e-6 Re）。
        let g = gamma_at(&src, sol.greatest.time_tt);
        assert!(
            (sol.gamma - g).abs() < 1e-6,
            "gamma={} must match independent √(x²+y²)={g} at time_tt (tol 1e-6)",
            sol.gamma
        );

        // 地表点往復: 前方射影し直すと (ξ,η)=(e.x,e.y)・ζ>0。
        let e = src
            .at(sol.greatest.time_tt)
            .expect("source.at should succeed at 2023 greatest-eclipse time");
        let phi = sol.greatest.position.lat.radians().0;
        let lam = sol.greatest.position.lon.radians().0;
        let obs = observer_geocentric(&Ellipsoid::WGS84, phi, 0.0);
        let r = project_observer_to_fundamental(&obs, Radians::new(lam), &e);
        assert!(
            (r.xi - e.x).abs() < 1e-9,
            "ξ={} must round-trip to e.x={}",
            r.xi,
            e.x
        );
        assert!(
            (r.eta - e.y).abs() < 1e-9,
            "η={} must round-trip to e.y={}",
            r.eta,
            e.y
        );
        assert!(r.zeta > 0.0, "ζ={} must be sunward (>0)", r.zeta);

        // TT/UTC 整合。
        let want_utc = umbra_core::time::tt_to_utc(sol.greatest.time_tt)
            .expect("2023 is post-1972, tt_to_utc must succeed");
        assert_eq!(
            sol.greatest.time_utc, want_utc,
            "time_utc must equal tt_to_utc(time_tt)"
        );

        // path/duration は本スライス非責務（常に None）。
        assert_eq!(
            sol.greatest.path_width, None,
            "path_width must be None (S6b)"
        );
        assert_eq!(
            sol.greatest.central_duration, None,
            "central_duration must be None (S6b)"
        );
    }

    // ====================================================================
    // 追補B: 非中心（軸が地表を外す）→ Err(Solver(RootNotBracketed)) の結線
    // ====================================================================
    //
    // 「中心食でない＝影軸が地表を外す」経路（`shadow_axis_surface_point` 失敗）は
    // `solve_greatest_eclipse` レベルで未テストだった。静的 `ConstantSource` は時不変ゆえ
    // `solve_local_maximum` がブラケットできず **別理由** で Err になる（軸ミス検証にならない）。
    // 代わりに **時変** 合成供給源を使う: gamma=√(x²+y²)=|x| が中心で内部極小を持ち、かつ
    // その極小値が >1（X_MIN=1.2）になるよう x=X_MIN+K·(jd−center)²（y=0）を返す。
    // これにより `solve_local_maximum` は内部放物線最小を **成功裏に** ブラケットし、
    // 続く `shadow_axis_surface_point` が gamma=1.2>1 で地表交点を見つけられず
    // `Err(Solver(RootNotBracketed))` を返す ⇒ 結線（伝播）を縛る。

    /// 中心極小値 >1 を持つ時変合成供給源。`solve_local_maximum` を成功させた上で
    /// `shadow_axis_surface_point` を軸ミスで失敗させるための、本テスト専用の最小実装。
    struct AxisMissSource {
        /// gamma の極小がここ（区間内部）に来る。
        center_jd: f64,
        /// x = X_MIN + K·(jd−center)²。X_MIN>1 ゆえ全域で gamma=|x|>1（軸は地表を外す）。
        window: TimeInterval<TtInstant>,
    }

    /// 中心での x（=gamma の極小値）。1 超ゆえ軸は地表に届かない。
    const AXIS_MISS_X_MIN: f64 = 1.2;
    /// 放物線の曲率（小さく正で、窓内でも x>1 を保ちつつ明瞭な内部極小を作る）。
    const AXIS_MISS_K: f64 = 50.0;

    impl BesselianSource for AxisMissSource {
        /// x=X_MIN+K·(jd−center)²（>1）, y=0, 他は有限固定値。gamma=√(x²+0)=x が center で極小。
        fn at(&self, t: TtInstant) -> Result<InstantaneousBesselianElements, EclipseError> {
            let dj = t.jd2().jd() - self.center_jd;
            Ok(InstantaneousBesselianElements {
                x: AXIS_MISS_X_MIN + AXIS_MISS_K * dj * dj,
                y: 0.0,
                declination: Radians(0.0),
                mu: Radians(0.0),
                l1: 0.5,
                l2: -0.01,
                tan_f1: 0.0047,
                tan_f2: 0.0046,
                time_tt: t,
            })
        }

        /// center_jd を内部（center ± 0.05 day）に持つ窓。gamma 極小（=X_MIN）が厳密に内部に来る。
        fn fit_interval(&self) -> TimeInterval<TtInstant> {
            self.window
        }
    }

    /// 追補B（軸ミス結線）: 時変合成供給源で `solve_local_maximum` が **成功** し、その後
    /// `shadow_axis_surface_point` が gamma>1 で失敗して `Err(Solver(RootNotBracketed))` が
    /// `solve_greatest_eclipse` から伝播することを縛る。
    ///
    /// Err の出所が `solve_local_maximum` ではなく `shadow_axis_surface_point` である根拠:
    /// 供給源の gamma=|x|=X_MIN+K·(jd−center)² は center で **内部放物線極小** を持つので
    /// `solve_local_maximum` はブラケットに成功する（≠ 時不変 ConstantSource の RootNotBracketed）。
    /// その極小値（=center の gamma）は 1.2>1（下で独立確認）。よって地表交点が無く、
    /// RootNotBracketed は影軸貫通段で発生し、上位へ透過する。
    #[test]
    fn greatest_noncentral_axis_miss_propagates_root_not_bracketed() {
        let dt = EspenakMeeusDeltaT;
        let center_jd = 2_457_986.768;
        let src = AxisMissSource {
            center_jd,
            window: TimeInterval {
                start: TtInstant::from_jd2(JulianDate2::new(center_jd - 0.05, 0.0)),
                end: TtInstant::from_jd2(JulianDate2::new(center_jd + 0.05, 0.0)),
            },
        };

        // center の gamma が 1.2 > 1（軸ミスの前提）であることを独立に確認する。
        let g_center = gamma_at(&src, g_tt_jd(center_jd));
        assert!(
            (g_center - AXIS_MISS_X_MIN).abs() < 1e-12,
            "center gamma {g_center} should equal X_MIN={AXIS_MISS_X_MIN}"
        );
        assert!(
            g_center > 1.0,
            "center gamma {g_center} must be > 1 (axis misses Earth)"
        );

        let config = crate::config::EngineConfig::standard();
        let r = solve_greatest_eclipse(&src, &dt, &config);
        assert_eq!(
            r,
            Err(EclipseError::Solver(SolverError::RootNotBracketed)),
            "axis-miss (interior gamma min >1) must propagate Solver(RootNotBracketed), got {r:?}"
        );
    }
}
