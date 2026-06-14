//! 月影円錐の幾何（`docs/issues/ISSUE-019`、`docs/algorithms/05-shadow-cone.md`）。
//!
//! 太陽中心・月中心・両半径から、本影/半影円錐（軸・頂点・半頂角）を構成する。
//! 半頂角は `sin f1 = (R_sun + R_moon)/D`（半影）、`sin f2 = (R_sun − R_moon)/D`（本影）、
//! D は太陽-月中心間距離。頂点は軸上で月から `R_moon/sin f` の距離。
//! 本影頂点が基本面（地球側）より遠ければ皆既、手前なら金環（反本影）— 符号判定は §6 ベッセルで。

use umbra_core::{Radians, UnitVector3, Vector3};

use crate::error::EclipseError;

/// 月影円錐（半影・本影）。長さは km、角度は rad。
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ShadowCone {
    /// 軸上の基準点（月中心）。
    pub axis_origin: Vector3,
    /// 軸方向（太陽 → 月、すなわち影の伸びる向き）。
    pub axis_direction: UnitVector3,
    /// 本影錐の頂点。
    pub umbra_apex: Vector3,
    /// 半影錐の頂点（太陽側）。
    pub penumbra_apex: Vector3,
    /// 本影半頂角 f2。
    pub umbra_half_angle: Radians,
    /// 半影半頂角 f1。
    pub penumbra_half_angle: Radians,
}

/// 太陽中心 `sun`・月中心 `moon`（地心 km）と太陽半径・月半径（km）から影円錐を構成する。
pub fn shadow_cone(
    sun: Vector3,
    moon: Vector3,
    r_sun_km: f64,
    r_moon_km: f64,
) -> Result<ShadowCone, EclipseError> {
    let axis_direction = (moon - sun)
        .normalized()
        .ok_or(EclipseError::DegenerateGeometry)?;
    let dist = (moon - sun).norm();

    // 半頂角（asin 引数は丸め誤差を考慮しクランプ。numerical-policy §A5）。
    let sin_f1 = ((r_sun_km + r_moon_km) / dist).clamp(-1.0, 1.0); // 半影
    let sin_f2 = ((r_sun_km - r_moon_km) / dist).clamp(-1.0, 1.0); // 本影
    let f1 = sin_f1.asin();
    let f2 = sin_f2.asin();
    // sin_f1 は正の半径・正の距離で常に >0。sin_f2 ≤ 0 は R_sun ≤ R_moon（非物理含む）で
    // 本影頂点が無限遠/裏返る → 退化（仕様 05-shadow-cone §境界）。
    if sin_f2 <= 0.0 {
        return Err(EclipseError::DegenerateGeometry);
    }

    let axis = axis_direction.get();
    // 半影頂点: 太陽側（−軸方向）に月から R_moon/sin f1。
    let penumbra_apex = moon - axis.scale(r_moon_km / sin_f1);
    // 本影頂点: 反太陽側（+軸方向）に月から R_moon/sin f2。
    let umbra_apex = moon + axis.scale(r_moon_km / sin_f2);

    Ok(ShadowCone {
        axis_origin: moon,
        axis_direction,
        umbra_apex,
        penumbra_apex,
        umbra_half_angle: Radians(f2),
        penumbra_half_angle: Radians(f1),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use umbra_core::constants::{EARTH_EQUATORIAL_RADIUS_M, SOLAR_RADIUS_KM};
    use umbra_core::{JulianDate2, TdbInstant};
    use umbra_ephemeris::{Body, Ephemeris, EphemerisFrame, MockEphemeris, Origin};

    const R_SUN: f64 = SOLAR_RADIUS_KM;
    // 月半径 = k·Re（IauMean k=0.2725076, conventions §9）。
    const R_MOON: f64 = 0.2725076 * (EARTH_EQUATORIAL_RADIUS_M / 1000.0);

    fn pos(m: &MockEphemeris, body: Body) -> Vector3 {
        let t = TdbInstant::from_jd2(JulianDate2::from_jd(2_451_545.0));
        m.state(body, t, Origin::Geocenter, EphemerisFrame::Icrs)
            .unwrap()
            .position
    }

    #[test]
    fn axis_points_from_sun_toward_moon() {
        let m = MockEphemeris::central_total();
        let c = shadow_cone(pos(&m, Body::Sun), pos(&m, Body::Moon), R_SUN, R_MOON).unwrap();
        // 月(+x 近距離)は太陽(+x 1AU)より地球側 → 軸は −x（影は地球向き）。
        assert!((c.axis_direction.get().x + 1.0).abs() < 1e-9);
    }

    #[test]
    fn penumbra_is_wider_than_umbra() {
        let m = MockEphemeris::central_total();
        let c = shadow_cone(pos(&m, Body::Sun), pos(&m, Body::Moon), R_SUN, R_MOON).unwrap();
        assert!(c.penumbra_half_angle.0 > c.umbra_half_angle.0);
    }

    #[test]
    fn half_angles_match_sine_definition() {
        let m = MockEphemeris::central_total();
        let (s, mo) = (pos(&m, Body::Sun), pos(&m, Body::Moon));
        let dist = (mo - s).norm();
        let c = shadow_cone(s, mo, R_SUN, R_MOON).unwrap();
        assert!((c.penumbra_half_angle.0.sin() - (R_SUN + R_MOON) / dist).abs() < 1e-12);
        assert!((c.umbra_half_angle.0.sin() - (R_SUN - R_MOON) / dist).abs() < 1e-12);
    }

    #[test]
    fn total_umbra_apex_passes_earth_center() {
        // 近地点（大きい月）: 本影頂点が地心(x=0)より遠側（x<0）→ 本影が地球に届く＝皆既。
        let m = MockEphemeris::central_total();
        let c = shadow_cone(pos(&m, Body::Sun), pos(&m, Body::Moon), R_SUN, R_MOON).unwrap();
        assert!(c.umbra_apex.x < 0.0, "umbra_apex.x = {}", c.umbra_apex.x);
    }

    #[test]
    fn annular_umbra_apex_falls_short_of_earth_center() {
        // 遠地点（小さい月）: 本影頂点が地心の手前（x>0）→ 反本影が地球に届く＝金環。
        let m = MockEphemeris::clear_annular();
        let c = shadow_cone(pos(&m, Body::Sun), pos(&m, Body::Moon), R_SUN, R_MOON).unwrap();
        assert!(c.umbra_apex.x > 0.0, "umbra_apex.x = {}", c.umbra_apex.x);
    }

    #[test]
    fn penumbra_apex_is_on_sun_side_at_correct_distance() {
        let m = MockEphemeris::central_total();
        let moon = pos(&m, Body::Moon);
        let c = shadow_cone(pos(&m, Body::Sun), moon, R_SUN, R_MOON).unwrap();
        // 太陽側 = +x 方向（月より大きい x）。
        assert!(c.penumbra_apex.x > moon.x);
        // 月からの距離は厳密に R_moon / sin f1（除算の取り違えを検出）。
        let expected = R_MOON / c.penumbra_half_angle.0.sin();
        assert!(((c.penumbra_apex - moon).norm() - expected).abs() < 1e-6);
    }

    #[test]
    fn umbra_apex_distance_matches_definition() {
        let m = MockEphemeris::central_total();
        let moon = pos(&m, Body::Moon);
        let c = shadow_cone(pos(&m, Body::Sun), moon, R_SUN, R_MOON).unwrap();
        let expected = R_MOON / c.umbra_half_angle.0.sin();
        assert!(((c.umbra_apex - moon).norm() - expected).abs() < 1e-6);
    }

    #[test]
    fn degenerate_when_sun_coincides_with_moon() {
        let p = Vector3::new(1.0, 2.0, 3.0);
        assert_eq!(
            shadow_cone(p, p, R_SUN, R_MOON).unwrap_err(),
            EclipseError::DegenerateGeometry
        );
    }

    #[test]
    fn degenerate_when_radii_equal() {
        // R_sun == R_moon → sin f2 = 0（本影頂点が無限遠）→ DegenerateGeometry。
        let m = MockEphemeris::central_total();
        assert_eq!(
            shadow_cone(pos(&m, Body::Sun), pos(&m, Body::Moon), R_SUN, R_SUN).unwrap_err(),
            EclipseError::DegenerateGeometry
        );
    }

    #[test]
    fn degenerate_when_moon_larger_than_sun() {
        // 非物理 R_moon > R_sun → sin f2 < 0 → DegenerateGeometry（不等号 <= の検証）。
        let m = MockEphemeris::central_total();
        assert_eq!(
            shadow_cone(pos(&m, Body::Sun), pos(&m, Body::Moon), R_MOON, R_SUN).unwrap_err(),
            EclipseError::DegenerateGeometry
        );
    }

    #[test]
    fn tilted_axis_apex_offsets_match_definition() {
        // 斜め配置（一般軸）でも頂点 = 月 ± (R_moon/sin f)·軸 が各成分で一致。
        let sun = Vector3::new(120_000_000.0, 80_000_000.0, 20_000_000.0);
        let moon = Vector3::new(300_000.0, 150_000.0, 40_000.0);
        let c = shadow_cone(sun, moon, R_SUN, R_MOON).unwrap();
        let axis = c.axis_direction.get();
        let u_expected = moon + axis.scale(R_MOON / c.umbra_half_angle.0.sin());
        let p_expected = moon - axis.scale(R_MOON / c.penumbra_half_angle.0.sin());
        assert!((c.umbra_apex - u_expected).norm() < 1e-6);
        assert!((c.penumbra_apex - p_expected).norm() < 1e-6);
    }
}
