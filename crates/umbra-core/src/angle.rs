//! 角度型（`docs/conventions.md` §2）。
//!
//! 内部表現はラジアン、公開入出力では度も提供する。正規化は用途別に分け、混在させない。
//! - [`Radians::normalized_signed`] … `[-π, π)`（経度・時角など循環量）
//! - [`Radians::normalized_two_pi`] … `[0, 2π)`（赤経・恒星時など）

use core::f64::consts::PI;

const TWO_PI: f64 = 2.0 * PI;

/// ラジアン角。
#[derive(Clone, Copy, Debug, PartialEq, PartialOrd)]
#[cfg_attr(feature = "serde", derive(serde::Serialize))]
#[cfg_attr(feature = "serde", serde(transparent))]
pub struct Radians(pub f64);

/// 度角。
#[derive(Clone, Copy, Debug, PartialEq, PartialOrd)]
#[cfg_attr(feature = "serde", derive(serde::Serialize))]
#[cfg_attr(feature = "serde", serde(transparent))]
pub struct Degrees(pub f64);

impl Radians {
    /// 値を包む。
    pub const fn new(value: f64) -> Self {
        Radians(value)
    }

    /// 度へ変換。
    pub fn to_degrees(self) -> Degrees {
        Degrees(self.0 * 180.0 / PI)
    }

    /// `[-π, π)` へ正規化（循環量。conventions §2）。
    pub fn normalized_signed(self) -> Self {
        let x = self.0.rem_euclid(TWO_PI); // [0, 2π)
        Radians(if x >= PI { x - TWO_PI } else { x }) // [-π, π)
    }

    /// `[0, 2π)` へ正規化。
    pub fn normalized_two_pi(self) -> Self {
        Radians(self.0.rem_euclid(TWO_PI))
    }
}

impl Degrees {
    /// 値を包む。
    pub const fn new(value: f64) -> Self {
        Degrees(value)
    }

    /// ラジアンへ変換。
    pub fn to_radians(self) -> Radians {
        Radians(self.0 * PI / 180.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-12;

    fn close(a: f64, b: f64) -> bool {
        (a - b).abs() < EPS
    }

    #[test]
    fn degrees_radians_round_trip() {
        let d = Degrees::new(133.508);
        assert!(close(d.to_radians().to_degrees().0, 133.508));
    }

    #[test]
    fn signed_normalization_maps_into_half_open_interval() {
        // 3π → [-π,π) では -π
        assert!(close(Radians::new(3.0 * PI).normalized_signed().0, -PI));
        // -π/2 はそのまま
        assert!(close(
            Radians::new(-PI / 2.0).normalized_signed().0,
            -PI / 2.0
        ));
        // 2π → 0
        assert!(close(Radians::new(TWO_PI).normalized_signed().0, 0.0));
        // ちょうど π は半開区間の下端 -π へ
        assert!(close(Radians::new(PI).normalized_signed().0, -PI));
    }

    #[test]
    fn two_pi_normalization_is_non_negative() {
        let r = Radians::new(-0.1).normalized_two_pi().0;
        assert!((0.0..TWO_PI).contains(&r));
        assert!(close(r, TWO_PI - 0.1));
    }
}
