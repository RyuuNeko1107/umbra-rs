//! 瞬時ベッセル要素（`docs/issues/ISSUE-021`、`docs/algorithms/06-besselian.md`）。
//!
//! 太陽・月の地心位置と両半径から、基本面（z=影軸）における幾何量を求める:
//! - x, y: 影軸が基本面を貫く点の座標（基本面基底 x̂,ŷ 上、単位 Re）。`gamma = √(x²+y²)`。
//! - d: 影軸（太陽方向）の赤緯。
//! - l1, l2: 半影/本影錐の基本面での半径（Re, 符号付き）。**l2<0 ⇒ 皆既 / l2>0 ⇒ 金環**（正本 B1）。
//! - tan f1, tan f2: 半影/本影錐の半頂角の正接。
//!
//! `besselian_elements`（μ を除く幾何）に加え、`besselian_mu`（μ = θ_ERA(UT1) − α_axis, CIO ベース,
//! ISSUE-039）と `besselian_elements_at`（実 apparent 位置 CIRS から μ 込みのフル瞬時要素を構成,
//! ISSUE-021）を提供する。μ のみ UT1（→ΔT）依存、他は TT 基準。

use umbra_core::constants::EARTH_EQUATORIAL_RADIUS_M;
use umbra_core::deltat::{tt_to_ut1, DeltaTModel};
use umbra_core::{Radians, TtInstant, UnitVector3, Ut1Instant, Vector3};
use umbra_ephemeris::apparent::{moon_apparent_cirs, sun_apparent_cirs};
use umbra_ephemeris::frames::earth_rotation_angle;

use crate::error::EclipseError;
use crate::fundamental::fundamental_plane_basis;
use crate::shadow::shadow_cone;

/// 瞬時ベッセル要素（μ を除く幾何由来分）。長さは Re 単位。
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BesselianElements {
    /// 影軸交点の x 座標（Re, 東）。
    pub x: f64,
    /// 影軸交点の y 座標（Re, 北）。
    pub y: f64,
    /// 影軸（太陽方向）の赤緯 d。
    pub declination: Radians,
    /// 半影錐半径 l1（Re, 正）。
    pub l1: f64,
    /// 本影錐半径 l2（Re, 符号付き。l2<0=皆既 / l2>0=金環）。
    pub l2: f64,
    /// 半影半頂角の正接 tan f1。
    pub tan_f1: f64,
    /// 本影半頂角の正接 tan f2。
    pub tan_f2: f64,
}

impl BesselianElements {
    /// 影軸の地心最小距離 `gamma = √(x²+y²)`（Re）。
    pub fn gamma(&self) -> f64 {
        (self.x * self.x + self.y * self.y).sqrt()
    }
}

/// 太陽中心・月中心（地心 km）と太陽半径・月半径（km）から瞬時ベッセル要素を構成する。
pub fn besselian_elements(
    sun: Vector3,
    moon: Vector3,
    r_sun_km: f64,
    r_moon_km: f64,
) -> Result<BesselianElements, EclipseError> {
    let re_km = EARTH_EQUATORIAL_RADIUS_M / 1000.0;
    let cone = shadow_cone(sun, moon, r_sun_km, r_moon_km)?;

    // ẑ = 太陽方向（影軸の逆向き）。基本面基底を張る。
    let z_axis = (sun - moon)
        .normalized()
        .ok_or(EclipseError::DegenerateGeometry)?;
    let basis = fundamental_plane_basis(z_axis)?;
    let z = z_axis.get();

    // 影軸交点 = 月位置を基本面（x̂,ŷ）へ射影（軸成分は基本面で 0）。
    let x = moon.dot(basis.x_axis.get()) / re_km;
    let y = moon.dot(basis.y_axis.get()) / re_km;

    // d = 影軸の赤緯（ẑ の z 成分 = sin d）。
    let declination = Radians(z.z.asin());

    let tan_f1 = cone.penumbra_half_angle.0.tan();
    let tan_f2 = cone.umbra_half_angle.0.tan();

    // 錐半径 = (頂点の ζ 座標) × tan f。ζ = 頂点·ẑ（太陽向き正）。
    // 半影頂点は太陽側 ζ1>0 → l1>0。本影頂点 ζ2 が反太陽側(地球を越える)なら ζ2<0 → l2<0=皆既。
    let z1 = cone.penumbra_apex.dot(z) / re_km;
    let z2 = cone.umbra_apex.dot(z) / re_km;
    let l1 = z1 * tan_f1;
    let l2 = z2 * tan_f2;

    Ok(BesselianElements {
        x,
        y,
        declination,
        l1,
        l2,
        tan_f1,
        tan_f2,
    })
}

/// 瞬時ベッセル要素（μ 含む・TT ラベル付き）。長さ Re 無次元、角 rad。
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct InstantaneousBesselianElements {
    /// 影軸交点の x 座標（Re, 東）。
    pub x: f64,
    /// 影軸交点の y 座標（Re, 北）。
    pub y: f64,
    /// 影軸（太陽方向）の赤緯 d。
    pub declination: Radians,
    /// 影軸の見かけグリニッジ時角 μ = θ_ERA(UT1) − α_axis（CIO ベース, [0,2π)）。
    pub mu: Radians,
    /// 半影錐半径 l1（Re, 正）。
    pub l1: f64,
    /// 本影錐半径 l2（Re, 符号付き。l2<0=皆既 / l2>0=金環）。
    pub l2: f64,
    /// 半影半頂角の正接 tan f1。
    pub tan_f1: f64,
    /// 本影半頂角の正接 tan f2。
    pub tan_f2: f64,
    /// この要素の TT 時刻ラベル（μ のみ UT1 由来）。
    pub time_tt: TtInstant,
}

impl InstantaneousBesselianElements {
    /// 影軸の地心最小距離 `gamma = √(x²+y²)`（Re）。
    pub fn gamma(&self) -> f64 {
        (self.x * self.x + self.y * self.y).sqrt()
    }
}

/// 影軸（CIRS, 太陽向き単位ベクトル）と UT1 からベッセル時角 μ を構成する。
///
/// `μ = θ_ERA(UT1) − α_axis`、`α_axis = atan2(axis.y, axis.x)`（CIRS 赤経）、`[0,2π)` 正規化。
/// θ_ERA と α_axis はともに CIO 起点で測るため μ は CIO で完全に閉じる（GAST も EO も登場しない。
/// `docs/algorithms/06-besselian.md` B5）。μ は θ_ERA 経由で **UT1（→ΔT）依存**。
pub fn besselian_mu(axis_cirs: UnitVector3, time_ut1: Ut1Instant) -> Radians {
    let a = axis_cirs.get();
    let alpha_axis = a.y.atan2(a.x);
    let era = earth_rotation_angle(time_ut1).0;
    Radians::new(era - alpha_axis).normalized_two_pi()
}

/// 時刻 TT における瞬時ベッセル要素を、実 apparent 位置（CIRS）から構成する。
///
/// 太陽・月の見かけ地心位置（光行時間＋光行差＋歳差章動, CIRS）を `*_apparent_cirs(time_tt)` で得て、
/// 幾何（x,y,d,l1,l2,tan f）を [`besselian_elements`] で、時角 μ を [`besselian_mu`] で構成する。
/// μ は `tt_to_ut1(time_tt, delta_t)` 由来の UT1 を使う（要素は TT ラベル、μ のみ UT1 依存）。
/// `r_sun_km` = 太陽物理半径、`r_moon_km` = 月半径（= k·Re, k は LunarRadiusModel）。
pub fn besselian_elements_at<M: DeltaTModel>(
    time_tt: TtInstant,
    r_sun_km: f64,
    r_moon_km: f64,
    delta_t: &M,
) -> Result<InstantaneousBesselianElements, EclipseError> {
    let sun = sun_apparent_cirs(time_tt);
    let moon = moon_apparent_cirs(time_tt);
    let ut1 = tt_to_ut1(time_tt, delta_t);
    instantaneous_from_cirs(sun, moon, r_sun_km, r_moon_km, time_tt, ut1)
}

/// 見かけ CIRS 位置（太陽・月, km）と UT1 から瞬時ベッセル要素を組み立てる共有ルーチン。
///
/// 幾何（x,y,d,l1,l2,tan f）を [`besselian_elements`] で、時角 μ を [`besselian_mu`]（axis=(sun−moon)
/// 正規化）で構成する。[`besselian_elements_at`]（具象 VSOP/ELP apparent）と
/// `InstantaneousEvaluator`（ISSUE-043 S2 の `apparent_cirs<E>` 経由・暦ジェネリック）が共用し、
/// apparent 取得元のみが異なる（組立は同一）。
pub(crate) fn instantaneous_from_cirs(
    sun_cirs: Vector3,
    moon_cirs: Vector3,
    r_sun_km: f64,
    r_moon_km: f64,
    time_tt: TtInstant,
    time_ut1: Ut1Instant,
) -> Result<InstantaneousBesselianElements, EclipseError> {
    let geom = besselian_elements(sun_cirs, moon_cirs, r_sun_km, r_moon_km)?;
    let axis = (sun_cirs - moon_cirs)
        .normalized()
        .ok_or(EclipseError::DegenerateGeometry)?;
    let mu = besselian_mu(axis, time_ut1);
    Ok(InstantaneousBesselianElements {
        x: geom.x,
        y: geom.y,
        declination: geom.declination,
        mu,
        l1: geom.l1,
        l2: geom.l2,
        tan_f1: geom.tan_f1,
        tan_f2: geom.tan_f2,
        time_tt,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use umbra_core::constants::SOLAR_RADIUS_KM;
    use umbra_core::{JulianDate2, TdbInstant};
    use umbra_ephemeris::{Body, Ephemeris, EphemerisFrame, MockEphemeris, Origin};

    const R_SUN: f64 = SOLAR_RADIUS_KM;
    const R_MOON: f64 = 0.2725076 * (EARTH_EQUATORIAL_RADIUS_M / 1000.0);

    fn elems(m: &MockEphemeris) -> BesselianElements {
        let t = TdbInstant::from_jd2(JulianDate2::from_jd(2_451_545.0));
        let pos = |b| {
            m.state(b, t, Origin::Geocenter, EphemerisFrame::Icrs)
                .unwrap()
                .position
        };
        besselian_elements(pos(Body::Sun), pos(Body::Moon), R_SUN, R_MOON).unwrap()
    }

    #[test]
    fn central_total_has_zero_gamma_and_negative_l2() {
        let e = elems(&MockEphemeris::central_total());
        assert!(e.gamma() < 1e-6, "gamma = {}", e.gamma());
        assert!(e.l2 < 0.0, "l2 = {} (皆既は負)", e.l2);
        assert!(e.l1 > 0.0);
        // |l2| は実日食域（~0.005–0.05 Re）。符号だけでなく大きさも縛る（H3）。
        assert!((0.005..0.05).contains(&e.l2.abs()), "|l2| = {}", e.l2.abs());
    }

    #[test]
    fn clear_annular_has_zero_gamma_and_positive_l2() {
        let e = elems(&MockEphemeris::clear_annular());
        assert!(e.gamma() < 1e-6, "gamma = {}", e.gamma());
        assert!(e.l2 > 0.0, "l2 = {} (金環は正)", e.l2);
        assert!((0.005..0.05).contains(&e.l2.abs()), "|l2| = {}", e.l2.abs());
    }

    #[test]
    fn partial_axis_misses_earth_disc() {
        // |gamma| が 1 を超え 1.55 未満（軸は地球を外すが半影は届く）。
        let g = elems(&MockEphemeris::clear_partial()).gamma();
        assert!((1.0..1.55).contains(&g), "gamma = {g}");
    }

    #[test]
    fn shadow_miss_has_large_gamma() {
        let g = elems(&MockEphemeris::shadow_misses_earth()).gamma();
        assert!(g > 1.55, "gamma = {g}");
    }

    #[test]
    fn penumbra_radius_is_about_half_earth_radius() {
        // 実際の日食の l1 は ~0.53–0.57 Re。
        let e = elems(&MockEphemeris::central_total());
        assert!((0.5..0.6).contains(&e.l1), "l1 = {}", e.l1);
        assert!(e.l1 > e.l2.abs(), "l1 should exceed |l2|");
    }

    #[test]
    fn declination_is_near_zero_for_equatorial_mock() {
        // Mock の太陽は赤道上（+x）。
        let e = elems(&MockEphemeris::central_total());
        assert!(e.declination.0.abs() < 1e-6);
    }

    #[test]
    fn tan_half_angles_match_cone() {
        let e = elems(&MockEphemeris::central_total());
        assert!(e.tan_f1 > e.tan_f2 && e.tan_f2 > 0.0);
    }

    #[test]
    fn gamma_combines_both_components() {
        // 両成分が非ゼロ: x²+y² の各項・符号を区別する（Mock 構成は片成分のみ非ゼロ）。
        let e = BesselianElements {
            x: 3.0,
            y: 4.0,
            declination: Radians(0.0),
            l1: 0.5,
            l2: -0.01,
            tan_f1: 0.0047,
            tan_f2: 0.0046,
        };
        assert!((e.gamma() - 5.0).abs() < 1e-12);
    }

    #[test]
    fn out_of_plane_moon_sets_y_and_negative_declination() {
        // 月を +z へオフセット → 軸 ẑ=(sun−moon) は −z 成分 → d<0、射影は ŷ 方向（y≠0, x=0）。
        let re_km = EARTH_EQUATORIAL_RADIUS_M / 1000.0;
        let sun = Vector3::new(149_597_870.7, 0.0, 0.0);
        let moon = Vector3::new(384_400.0, 0.0, 8_000.0);
        let e = besselian_elements(sun, moon, R_SUN, R_MOON).unwrap();
        assert!(e.declination.0 < 0.0, "d = {}", e.declination.0);
        assert!(e.x.abs() < 1e-6, "x = {}", e.x);
        // 独立オラクル: ŷ≈天の北(z軸)なので y ≈ moon.z/Re（マジック値でなく入力から導出）。
        assert!((e.y - 8_000.0 / re_km).abs() < 0.01, "y = {}", e.y);
    }

    #[test]
    fn minus_z_moon_gives_positive_declination() {
        // 月を −z へ → ẑ は +z 成分 → d>0（d 符号の両方向を踏む, L3）。
        let sun = Vector3::new(149_597_870.7, 0.0, 0.0);
        let moon = Vector3::new(384_400.0, 0.0, -8_000.0);
        let e = besselian_elements(sun, moon, R_SUN, R_MOON).unwrap();
        assert!(e.declination.0 > 0.0, "d = {}", e.declination.0);
    }

    #[test]
    fn both_x_and_y_nonzero_match_independent_oracle() {
        // 月を +y(赤道東) と +z(北) の両方へオフセット → x,y とも非ゼロ。
        // 独立オラクル: x̂≈+y_eq, ŷ≈+z_eq なので x≈moon.y/Re, y≈moon.z/Re。
        // これで x̂/ŷ の取り違え・Re 割り忘れを検出する（H2）。
        let re_km = EARTH_EQUATORIAL_RADIUS_M / 1000.0;
        let sun = Vector3::new(149_597_870.7, 0.0, 0.0);
        let moon = Vector3::new(384_400.0, 5_000.0, 8_000.0);
        let e = besselian_elements(sun, moon, R_SUN, R_MOON).unwrap();
        assert!((e.x - 5_000.0 / re_km).abs() < 0.01, "x = {}", e.x);
        assert!((e.y - 8_000.0 / re_km).abs() < 0.01, "y = {}", e.y);
        assert!(e.x < e.y, "x={} should differ from y={}", e.x, e.y);
    }
}

// μ（ベッセル時角）＋ 実 apparent 位置での瞬時ベッセル要素（ISSUE-021/039）。
#[cfg(test)]
mod mu_tests {
    #![allow(clippy::excessive_precision)]

    use super::*;
    use core::f64::consts::TAU;
    use umbra_core::constants::SOLAR_RADIUS_KM;
    use umbra_core::{EspenakMeeusDeltaT, JulianDate2, TtInstant};

    const R_SUN: f64 = SOLAR_RADIUS_KM;
    const R_MOON: f64 = 0.2725076 * (EARTH_EQUATORIAL_RADIUS_M / 1000.0);

    fn unit(x: f64, y: f64, z: f64) -> UnitVector3 {
        Vector3::new(x, y, z).normalized().expect("non-zero axis")
    }
    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    // 一次オラクル: pyerfa 2.0.1.5 `erfa.era00(jd1,jd2)`（liberfa = SOFA, 独立実装）。
    //   μ_expected = (era00(jd1,jd2) − atan2(axis.y, axis.x)) mod 2π。
    //   各 case の (jd1,jd2,axis) は下記。era00 を pyerfa で取得し python 側で μ を算出・転記。

    /// T1: μ = θ_ERA(UT1) − α_axis を pyerfa 厳密一致で検証（4 象限・南北 d, tol 1e-12）。
    #[test]
    fn besselian_mu_matches_erfa_era_minus_alpha() {
        struct Case {
            jd1: f64,
            jd2: f64,
            ax: f64,
            ay: f64,
            az: f64,
            mu: f64,
        }
        let cases = [
            Case {
                jd1: 2_451_545.0,
                jd2: 0.0,
                ax: 1.0,
                ay: 0.5,
                az: 0.2,
                mu: 4.431_313_603_822_950_14e0,
            },
            Case {
                jd1: 2_458_000.5,
                jd2: 0.25,
                ax: -1.0,
                ay: 0.7,
                az: -0.3,
                mu: 5.032_118_767_655_029_86e0,
            },
            Case {
                jd1: 2_400_000.5,
                jd2: 53_736.0,
                ax: -0.8,
                ay: -0.6,
                az: 0.1,
                mu: 4.250_924_800_099_582_51e0,
            },
            Case {
                jd1: 2_459_731.75,
                jd2: 0.0,
                ax: 1.0,
                ay: -0.9,
                az: 0.95,
                mu: 3.736_615_603_114_827_88e-1,
            },
        ];
        for (i, c) in cases.iter().enumerate() {
            let ut1 = Ut1Instant::from_jd2(JulianDate2::new(c.jd1, c.jd2));
            let mu = besselian_mu(unit(c.ax, c.ay, c.az), ut1).0;
            assert!(
                close(mu, c.mu, 1e-12),
                "case {i}: mu={mu} expected {}",
                c.mu
            );
        }
    }

    /// T2: μ は常に [0,2π) へ正規化される（θ_ERA<α で差が負になる axis を含む）。
    #[test]
    fn besselian_mu_is_normalized_into_zero_two_pi() {
        let axes = [
            unit(1.0, 0.0, 0.0),
            unit(-1.0, 0.001, 0.0),
            unit(-1.0, -0.001, 0.0),
            unit(0.0, 1.0, 0.5),
            unit(0.0, -1.0, -0.5),
        ];
        for a in axes {
            for &(p1, p2) in &[(2_451_545.0, 0.0), (2_458_000.5, 0.5)] {
                let mu = besselian_mu(a, Ut1Instant::from_jd2(JulianDate2::new(p1, p2))).0;
                assert!((0.0..TAU).contains(&mu), "mu out of [0,2π): {mu}");
            }
        }
    }

    /// T3: α を引かない恒等実装(μ=θ_ERA)では落ちる。α≠0 で μ≠θ_ERA を独立に固定。
    #[test]
    fn besselian_mu_subtracts_alpha_axis() {
        let ut1 = Ut1Instant::from_jd2(JulianDate2::new(2_458_000.5, 0.0));
        let axis = unit(1.0, 1.0, 0.0); // α_axis = +π/4
        let era = earth_rotation_angle(ut1).0;
        let alpha = axis.get().y.atan2(axis.get().x);
        let expected = Radians::new(era - alpha).normalized_two_pi().0;
        let mu = besselian_mu(axis, ut1).0;
        assert!(close(mu, expected, 1e-12), "mu={mu} expected={expected}");
        let era_norm = Radians::new(era).normalized_two_pi().0;
        assert!(
            (mu - era_norm).abs() > 1e-6,
            "α not subtracted: mu={mu} era={era_norm}"
        );
    }

    /// T4: α_axis は (x,y) のみ依存・z 無関係（atan2 引数取り違え検出）。
    #[test]
    fn besselian_mu_ignores_axis_z_component() {
        let ut1 = Ut1Instant::from_jd2(JulianDate2::new(2_457_000.0, 0.3));
        let mu_a = besselian_mu(unit(0.6, 0.4, 0.1), ut1).0;
        let mu_b = besselian_mu(unit(0.6, 0.4, 0.9), ut1).0;
        assert!(
            close(mu_a, mu_b, 1e-12),
            "mu depends on z: {mu_a} vs {mu_b}"
        );
    }

    /// T11: besselian_mu は有限。
    #[test]
    fn besselian_mu_is_finite() {
        let ut1 = Ut1Instant::from_jd2(JulianDate2::new(2_451_545.0, 0.0));
        assert!(besselian_mu(unit(1.0, 0.2, -0.3), ut1).0.is_finite());
    }

    /// T5: 瞬時要素 gamma() = √(x²+y²)。
    #[test]
    fn instantaneous_gamma_combines_both_components() {
        let e = InstantaneousBesselianElements {
            x: 3.0,
            y: 4.0,
            declination: Radians(0.0),
            mu: Radians(1.0),
            l1: 0.5,
            l2: -0.01,
            tan_f1: 0.0047,
            tan_f2: 0.0046,
            time_tt: TtInstant::from_jd2(JulianDate2::new(2_451_545.0, 0.0)),
        };
        assert!(close(e.gamma(), 5.0, 1e-12), "gamma = {}", e.gamma());
    }

    /// 2017-08-21 最大食（greatest eclipse 18:25:32 UTC + ΔT≈69.184s → 18:26:41.184 TT）。
    fn tt_2017_max() -> TtInstant {
        TtInstant::from_jd2(JulianDate2::new(2_457_986.5, 7.685_322_222_222_221_72e-1))
    }

    /// T6: geom 部 == besselian_elements(apparent..)（合成同一性）。
    #[test]
    fn elements_at_geom_matches_besselian_elements() {
        let tt = tt_2017_max();
        let inst = besselian_elements_at(tt, R_SUN, R_MOON, &EspenakMeeusDeltaT).unwrap();
        let geom = besselian_elements(sun_apparent_cirs(tt), moon_apparent_cirs(tt), R_SUN, R_MOON)
            .unwrap();
        assert!(close(inst.x, geom.x, 1e-12), "x");
        assert!(close(inst.y, geom.y, 1e-12), "y");
        assert!(close(inst.declination.0, geom.declination.0, 1e-12), "d");
        assert!(close(inst.l1, geom.l1, 1e-12), "l1");
        assert!(close(inst.l2, geom.l2, 1e-12), "l2");
        assert!(close(inst.tan_f1, geom.tan_f1, 1e-12), "tan_f1");
        assert!(close(inst.tan_f2, geom.tan_f2, 1e-12), "tan_f2");
    }

    /// T7: μ == besselian_mu((sun−moon)正規化, tt_to_ut1(tt,&Δ))。
    #[test]
    fn elements_at_mu_matches_besselian_mu() {
        let tt = tt_2017_max();
        let inst = besselian_elements_at(tt, R_SUN, R_MOON, &EspenakMeeusDeltaT).unwrap();
        let axis = (sun_apparent_cirs(tt) - moon_apparent_cirs(tt))
            .normalized()
            .unwrap();
        let ut1 = tt_to_ut1(tt, &EspenakMeeusDeltaT);
        assert!(
            close(inst.mu.0, besselian_mu(axis, ut1).0, 1e-12),
            "mu mismatch"
        );
    }

    /// T8: time_tt ラベルは入力 tt を保持。
    #[test]
    fn elements_at_preserves_input_tt_label() {
        let tt = tt_2017_max();
        let inst = besselian_elements_at(tt, R_SUN, R_MOON, &EspenakMeeusDeltaT).unwrap();
        assert_eq!(inst.time_tt, tt);
    }

    /// T9: 2017-08-21 最大食 end-to-end サニティ。実装値 gamma=0.43671 は **NASA 公表 gamma=0.4367
    /// と 4 桁一致**（apparent→ベッセルの全チェーン検証）。d=11.86°≈太陽見かけ赤緯。
    #[test]
    fn elements_at_2017_total_eclipse_sanity() {
        let e = besselian_elements_at(tt_2017_max(), R_SUN, R_MOON, &EspenakMeeusDeltaT).unwrap();
        // NASA gamma=0.4367 を [0.43,0.44] で締め（モデル差 ΔT/k/平均月縁の余裕）。
        assert!(
            (0.43..0.44).contains(&e.gamma()),
            "gamma = {} (NASA 0.4367)",
            e.gamma()
        );
        // 太陽見かけ赤緯 ≈ +11.86° = 0.2070 rad（±0.3°）。
        assert!(
            close(e.declination.0, 0.2070, 5.0e-3),
            "d = {} rad",
            e.declination.0
        );
        assert!((0.0..TAU).contains(&e.mu.0), "mu = {}", e.mu.0);
        assert!((0.53..0.57).contains(&e.l1), "l1 = {}", e.l1);
        assert!(e.l2 < 0.0, "l2 = {} (皆既は負)", e.l2);
        assert!((0.003..0.05).contains(&e.l2.abs()), "|l2| = {}", e.l2.abs());
    }

    /// T10: 実位置で全フィールド有限。
    #[test]
    fn elements_at_all_fields_finite() {
        let e = besselian_elements_at(tt_2017_max(), R_SUN, R_MOON, &EspenakMeeusDeltaT).unwrap();
        for (n, v) in [
            ("x", e.x),
            ("y", e.y),
            ("d", e.declination.0),
            ("mu", e.mu.0),
            ("l1", e.l1),
            ("l2", e.l2),
            ("tan_f1", e.tan_f1),
            ("tan_f2", e.tan_f2),
        ] {
            assert!(v.is_finite(), "{n} not finite: {v}");
        }
    }
}
