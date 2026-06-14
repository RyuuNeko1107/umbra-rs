//! 瞬時ベッセル要素（`docs/issues/ISSUE-021`、`docs/algorithms/06-besselian.md`）。
//!
//! 太陽・月の地心位置と両半径から、基本面（z=影軸）における幾何量を求める:
//! - x, y: 影軸が基本面を貫く点の座標（基本面基底 x̂,ŷ 上、単位 Re）。`gamma = √(x²+y²)`。
//! - d: 影軸（太陽方向）の赤緯。
//! - l1, l2: 半影/本影錐の基本面での半径（Re, 符号付き）。**l2<0 ⇒ 皆既 / l2>0 ⇒ 金環**（正本 B1）。
//! - tan f1, tan f2: 半影/本影錐の半頂角の正接。
//!
//! μ（影軸のグリニッジ時角）は GAST=ERA(UT1) を要するため本段では算出しない（ISSUE-039）。

use umbra_core::constants::EARTH_EQUATORIAL_RADIUS_M;
use umbra_core::{Radians, Vector3};

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
    }

    #[test]
    fn clear_annular_has_zero_gamma_and_positive_l2() {
        let e = elems(&MockEphemeris::clear_annular());
        assert!(e.gamma() < 1e-6, "gamma = {}", e.gamma());
        assert!(e.l2 > 0.0, "l2 = {} (金環は正)", e.l2);
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
        // 月を +z へオフセット → 軸 ẑ=(sun−moon) は −z 成分 → d<0、射影は ŷ 方向（y≠0）。
        let sun = Vector3::new(149_597_870.7, 0.0, 0.0);
        let moon = Vector3::new(384_400.0, 0.0, 8_000.0);
        let e = besselian_elements(sun, moon, R_SUN, R_MOON).unwrap();
        assert!(e.declination.0 < 0.0, "d = {}", e.declination.0);
        assert!((1.2..1.3).contains(&e.y), "y = {}", e.y);
        assert!(e.x.abs() < 1e-6, "x = {}", e.x);
    }
}
