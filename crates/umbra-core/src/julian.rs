//! 2要素ユリウス日（`docs/conventions.md` §6 / `docs/numerical-policy.md` §A1）。
//!
//! 巨大な JD（≈2.45e6）と微小な日数差を 1 つの `f64` に押し込むと、エポック差を取る際に
//! 約 4.6e-5 s の桁落ちが生じ、±1s の精度目標を直接侵食する。これを避けるため、JD を
//! `part1`（整数日側）と `part2`（小数日側）の 2 要素で保持し、エポック減算は整数部側で行う。

use crate::constants::{J2000_JD, JULIAN_CENTURY_DAYS, JULIAN_MILLENNIUM_DAYS};

/// 2要素ユリウス日。`jd = part1 + part2`。
#[derive(Clone, Copy, Debug, PartialEq, PartialOrd)]
pub struct JulianDate2 {
    /// 大きい側（通常は整数日に正規化）。
    pub part1: f64,
    /// 小さい側（小数日）。
    pub part2: f64,
}

impl JulianDate2 {
    /// 2 要素から構築（正規化はしない）。
    pub const fn new(part1: f64, part2: f64) -> Self {
        JulianDate2 { part1, part2 }
    }

    /// 単一の JD から構築し正規化する。
    pub fn from_jd(jd: f64) -> Self {
        JulianDate2 {
            part1: jd,
            part2: 0.0,
        }
        .normalized()
    }

    /// `part2` を `[-0.5, 0.5)` に寄せ、整数分を `part1` へ移す。
    pub fn normalized(self) -> Self {
        let extra = self.part2.round();
        JulianDate2 {
            part1: self.part1 + extra,
            part2: self.part2 - extra,
        }
    }

    /// 合計 JD（表示・粗い比較用。精度クリティカルな差分には使わない）。
    pub fn jd(self) -> f64 {
        self.part1 + self.part2
    }

    /// `self − earlier` を日数で返す。2要素のまま差を取り桁落ちを避ける
    /// （`jd()` 同士の差は巨大 JD の合算で精度を失うため、間隔計算にはこちらを使う）。
    pub fn days_since(self, earlier: JulianDate2) -> f64 {
        (self.part1 - earlier.part1) + (self.part2 - earlier.part2)
    }

    /// 日数オフセットを加算（光行時間など）。`part2` へ足して再正規化する。
    pub fn add_days(self, days: f64) -> Self {
        JulianDate2 {
            part1: self.part1,
            part2: self.part2 + days,
        }
        .normalized()
    }

    /// J2000.0 からの経過ユリウス世紀。エポック減算を整数部側で厳密に行う。
    pub fn julian_centuries_since_j2000(self) -> f64 {
        ((self.part1 - J2000_JD) + self.part2) / JULIAN_CENTURY_DAYS
    }

    /// J2000.0 からの経過ユリウス千年（VSOP87 の引数 T）。
    pub fn julian_millennia_since_j2000(self) -> f64 {
        ((self.part1 - J2000_JD) + self.part2) / JULIAN_MILLENNIUM_DAYS
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn centuries_at_j2000_is_zero() {
        let t = JulianDate2::new(J2000_JD, 0.0);
        assert_eq!(t.julian_centuries_since_j2000(), 0.0);
    }

    #[test]
    fn half_day_offset_is_exact() {
        let t = JulianDate2::new(J2000_JD, 0.5);
        assert!((t.julian_centuries_since_j2000() - 0.5 / JULIAN_CENTURY_DAYS).abs() < 1e-18);
    }

    /// 2要素表現は、巨大 JD 近傍の微小オフセットを失わない（桁落ち対策の回帰）。
    #[test]
    fn tiny_offset_is_preserved_via_two_part() {
        let tiny_days = 1e-9; // ≈ 8.6e-5 s
        let t = JulianDate2::new(J2000_JD, tiny_days);
        let centuries = t.julian_centuries_since_j2000();
        let expected = tiny_days / JULIAN_CENTURY_DAYS;
        // 相対誤差で厳密一致に近いこと（part2 が独立保持されるため）。
        assert!((centuries - expected).abs() < expected * 1e-9);
        assert!(centuries > 0.0);
    }

    #[test]
    fn days_since_preserves_subsecond_offset() {
        // 巨大 JD 近傍の 1 秒差を桁落ちなく取り出せる（jd() 同士の差では失われる）。
        let base = JulianDate2::new(2_460_000.5, 0.0);
        let plus_one_sec = base.add_days(1.0 / 86_400.0);
        let d = plus_one_sec.days_since(base) * 86_400.0;
        assert!((d - 1.0).abs() < 1e-9, "d = {d}");
    }

    #[test]
    fn add_days_normalizes_part2() {
        let t = JulianDate2::new(J2000_JD, 0.4).add_days(0.4);
        assert!(t.part2.abs() < 0.5);
        assert!((t.jd() - (J2000_JD + 0.8)).abs() < 1e-12);
    }

    #[test]
    fn from_jd_round_trips() {
        let t = JulianDate2::from_jd(2_460_000.25);
        assert!((t.jd() - 2_460_000.25).abs() < 1e-9);
    }

    #[test]
    fn millennia_at_j2000_is_zero() {
        assert_eq!(
            JulianDate2::new(J2000_JD, 0.0).julian_millennia_since_j2000(),
            0.0
        );
    }

    #[test]
    fn one_millennium_after_j2000() {
        let t = JulianDate2::new(J2000_JD + JULIAN_MILLENNIUM_DAYS, 0.0);
        assert!((t.julian_millennia_since_j2000() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn millennia_is_one_tenth_of_centuries() {
        let t = JulianDate2::new(J2000_JD + 12_345.0, 0.25);
        let m = t.julian_millennia_since_j2000();
        let c = t.julian_centuries_since_j2000();
        assert!(m != 0.0 && (m * 10.0 - c).abs() < 1e-15);
    }
}
