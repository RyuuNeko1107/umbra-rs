//! 地球楕円体と観測者の地心座標（`docs/issues/ISSUE-010`/`ISSUE-011`、`docs/conventions.md` §4）。
//!
//! 既定は WGS84。測地緯度 → 地心緯度、観測者の扁平補正済み地心動径成分 ρsinφ′/ρcosφ′
//! （視差・ベッセル観測者射影で使用）、および地球固定直交座標（ITRS/ECEF）への変換を提供する。

use crate::constants::{EARTH_EQUATORIAL_RADIUS_M, WGS84_FLATTENING};
use crate::vector::Vector3;

/// 回転楕円体（長半径 a と扁平率 f）。
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Ellipsoid {
    /// 長半径（赤道半径）a \[m\]。
    pub a_m: f64,
    /// 扁平率 f。
    pub f: f64,
}

impl Ellipsoid {
    /// WGS84（conventions §4.1）。
    pub const WGS84: Ellipsoid = Ellipsoid {
        a_m: EARTH_EQUATORIAL_RADIUS_M,
        f: WGS84_FLATTENING,
    };

    /// 短半径 b = a(1 − f) \[m\]。
    pub fn b_m(&self) -> f64 {
        self.a_m * (1.0 - self.f)
    }

    /// 第一離心率の二乗 e² = f(2 − f)。
    pub fn e2(&self) -> f64 {
        self.f * (2.0 - self.f)
    }
}

/// 観測者の扁平補正済み地心動径成分（単位: 地球赤道半径 Re）。
///
/// ρsinφ′・ρcosφ′（Meeus *Astronomical Algorithms* Ch.11）。視差・ベッセル観測者射影で使う。
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GeocentricObserver {
    /// ρ·sin(φ′)（Re 単位）。
    pub rho_sin_phi_prime: f64,
    /// ρ·cos(φ′)（Re 単位）。
    pub rho_cos_phi_prime: f64,
}

/// 測地緯度 `geodetic_lat`（rad）・楕円体高 `height_m`（m）から ρsinφ′/ρcosφ′ を求める。
///
/// Meeus Ch.11。簡約緯度 u = atan((b/a)·tanφ) を用いる（要確認: 式番号）。
pub fn observer_geocentric(
    ellipsoid: &Ellipsoid,
    geodetic_lat: f64,
    height_m: f64,
) -> GeocentricObserver {
    let b_over_a = 1.0 - ellipsoid.f;
    let u = (b_over_a * geodetic_lat.tan()).atan();
    let h_over_a = height_m / ellipsoid.a_m;
    GeocentricObserver {
        rho_sin_phi_prime: b_over_a * u.sin() + h_over_a * geodetic_lat.sin(),
        rho_cos_phi_prime: u.cos() + h_over_a * geodetic_lat.cos(),
    }
}

/// 測地緯度 → 地心緯度（rad）。`tan φ′ = (1 − e²) tan φ`。
pub fn geodetic_to_geocentric_latitude(ellipsoid: &Ellipsoid, geodetic_lat: f64) -> f64 {
    ((1.0 - ellipsoid.e2()) * geodetic_lat.tan()).atan()
}

/// 測地座標（緯度 φ・東経 λ・楕円体高 h）→ 地球固定直交座標 ITRS/ECEF（km）。
///
/// 標準式: N = a / √(1 − e² sin²φ)、X=(N+h)cosφcosλ、Y=(N+h)cosφsinλ、Z=(N(1−e²)+h)sinφ。
pub fn geodetic_to_ecef_km(
    ellipsoid: &Ellipsoid,
    geodetic_lat: f64,
    east_longitude: f64,
    height_m: f64,
) -> Vector3 {
    let e2 = ellipsoid.e2();
    let sin_lat = geodetic_lat.sin();
    let cos_lat = geodetic_lat.cos();
    let n = ellipsoid.a_m / (1.0 - e2 * sin_lat * sin_lat).sqrt();
    let x_m = (n + height_m) * cos_lat * east_longitude.cos();
    let y_m = (n + height_m) * cos_lat * east_longitude.sin();
    let z_m = (n * (1.0 - e2) + height_m) * sin_lat;
    Vector3::new(x_m / 1000.0, y_m / 1000.0, z_m / 1000.0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use core::f64::consts::PI;

    const WGS84: Ellipsoid = Ellipsoid::WGS84;

    #[test]
    fn semi_minor_and_e2() {
        assert!(
            (WGS84.b_m() - 6_356_752.314_245).abs() < 1e-3,
            "b = {}",
            WGS84.b_m()
        );
        // e² ≈ 0.00669437999014
        assert!((WGS84.e2() - 0.006_694_379_990_14).abs() < 1e-12);
    }

    #[test]
    fn observer_at_equator_sea_level() {
        let o = observer_geocentric(&WGS84, 0.0, 0.0);
        assert!((o.rho_cos_phi_prime - 1.0).abs() < 1e-12);
        assert!(o.rho_sin_phi_prime.abs() < 1e-12);
    }

    #[test]
    fn observer_at_pole_uses_b_over_a() {
        let o = observer_geocentric(&WGS84, PI / 2.0, 0.0);
        assert!(o.rho_cos_phi_prime.abs() < 1e-9);
        assert!((o.rho_sin_phi_prime - (1.0 - WGS84.f)).abs() < 1e-9);
    }

    #[test]
    fn geocentric_latitude_is_smaller_at_mid_latitude() {
        let phi = 45.0_f64.to_radians();
        let phi_prime = geodetic_to_geocentric_latitude(&WGS84, phi);
        assert!(phi_prime < phi);
        // 45°では差は約 11.5′
        assert!((phi - phi_prime).to_degrees() * 60.0 > 11.0);
    }

    #[test]
    fn ecef_equator_and_pole() {
        let eq = geodetic_to_ecef_km(&WGS84, 0.0, 0.0, 0.0);
        assert!((eq.x - 6378.137).abs() < 1e-3);
        assert!(eq.y.abs() < 1e-9 && eq.z.abs() < 1e-9);

        let pole = geodetic_to_ecef_km(&WGS84, PI / 2.0, 0.0, 0.0);
        assert!((pole.z - WGS84.b_m() / 1000.0).abs() < 1e-6);
        assert!(pole.x.abs() < 1e-6);
    }
}
