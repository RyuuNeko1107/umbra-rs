//! 物理・規約定数（`docs/conventions.md` §4.1）。
//!
//! いわゆる *magic number* はすべて本モジュールへ集約し、コードからは本定数を参照する
//! （`docs/conventions.md` §11）。値はいずれも定義値または採用固定値。
//! EOP・ΔT・閏秒のような時変・データ駆動の量はここに含めない。

/// 真空中の光速 \[km/s\]。IAU の SI 定義値（無誤差）。光行時間・光行差で使用。
pub const SPEED_OF_LIGHT_KM_S: f64 = 299_792.458;

/// 天文単位 \[km\]。IAU 2012 (Resolution B2) 定義値。
pub const ASTRONOMICAL_UNIT_KM: f64 = 149_597_870.7;

/// TT − TAI \[s\]。IAU 1991 による定数オフセット（閏秒は UTC↔TAI 側で別途扱う）。
pub const TT_MINUS_TAI_SECONDS: f64 = 32.184;

/// 地球赤道半径（WGS84 長半径 a）\[m\]。ベッセル無次元化の基準 Re。
pub const EARTH_EQUATORIAL_RADIUS_M: f64 = 6_378_137.0;

/// 太陽公称半径 \[km\]。IAU 2015 Resolution B3。影円錐・視半径に使用（conventions §9）。
pub const SOLAR_RADIUS_KM: f64 = 696_000.0;

/// WGS84 扁平率 f。
pub const WGS84_FLATTENING: f64 = 1.0 / 298.257_223_563;

/// J2000.0 のユリウス日（TT, 2000-01-01 12:00:00 TT）。
pub const J2000_JD: f64 = 2_451_545.0;

/// 1 ユリウス世紀の日数。
pub const JULIAN_CENTURY_DAYS: f64 = 36_525.0;

/// 1 ユリウス千年の日数（VSOP87 の引数 T で使用）。
pub const JULIAN_MILLENNIUM_DAYS: f64 = 365_250.0;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn light_time_one_au_is_about_499_seconds() {
        // 太陽光が 1 AU を進む時間 ≈ 499 s（光行時間の桁感の妥当性確認）。
        let light_time = ASTRONOMICAL_UNIT_KM / SPEED_OF_LIGHT_KM_S;
        assert!(
            (light_time - 499.0).abs() < 1.0,
            "light_time = {light_time}"
        );
    }

    #[test]
    fn julian_century_and_millennium_are_consistent() {
        assert_eq!(JULIAN_MILLENNIUM_DAYS, JULIAN_CENTURY_DAYS * 10.0);
    }

    #[test]
    fn wgs84_flattening_matches_inverse_298() {
        // f = 1/298.257223563 ≈ 0.00335281（runtime 比較で定数 assert を避ける）。
        let inv_f = 1.0 / WGS84_FLATTENING;
        assert!((inv_f - 298.257_223_563).abs() < 1e-6, "1/f = {inv_f}");
    }
}
