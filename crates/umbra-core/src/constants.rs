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

/// 太陽の Schwarzschild 半径相当（無次元）。SOFA/ERFA `ab.c` の定義値 `SRS = 2·G·M_⊙/(c²·au)`。
/// 相対論的光行差 `iauAb` の微項 `w2 = SRS/s`（`s` = 太陽-観測者距離 \[au\]）で使用する。
/// これは `iauAb` 内に常時含まれる微小項（~SRS/au ≈ 2e-8）で、角度依存の太陽光偏向 `iauLd`
/// （Standard 既定 OFF）とは別物（`docs/algorithms/03-ephemeris.md` E8/E9）。
pub const SRS: f64 = 1.974_125_743_36e-8;

/// 秒角 → ラジアン変換係数（= π / 648000）。SOFA の `DAS2R` 相当。
/// 歳差・章動・黄道傾斜など、秒角で表される級数係数のラジアン化に用いる
/// （IAU2006/2000A の係数は秒角単位。`docs/algorithms/02-frames.md`）。
pub const ARCSEC_TO_RAD: f64 = core::f64::consts::PI / 648_000.0;

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
    fn arcsec_to_rad_is_pi_over_648000() {
        // 180° = 648000″ = π rad。半周分の秒角がちょうど π になる。
        assert!((648_000.0 * ARCSEC_TO_RAD - core::f64::consts::PI).abs() < 1e-15);
        // 1″ ≈ 4.8481368e-6 rad（桁感）。
        assert!((ARCSEC_TO_RAD - 4.848_136_811_095_36e-6).abs() < 1e-18);
    }

    #[test]
    fn wgs84_flattening_matches_inverse_298() {
        // f = 1/298.257223563 ≈ 0.00335281（runtime 比較で定数 assert を避ける）。
        let inv_f = 1.0 / WGS84_FLATTENING;
        assert!((inv_f - 298.257_223_563).abs() < 1e-6, "1/f = {inv_f}");
    }
}
