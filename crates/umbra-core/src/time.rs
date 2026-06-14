//! 時刻系と変換 UTC ↔ TAI ↔ TT（`docs/issues/ISSUE-006`、`docs/algorithms/01-time-scales.md`）。
//!
//! 各時刻は対応する時刻系の [`JulianDate2`] として保持する。変換は一様時系 TAI を経由する:
//! `TAI = UTC + ΔAT(UTC)`（ΔAT=閏秒）、`TT = TAI + 32.184 s`（定数, conventions §4.1）。
//!
//! ΔAT は IERS の閏秒テーブル（1972– の公開事実データ）を組み込む。**1972 年より前は
//! 本テーブルでは未定義**で `TimeError::MissingLeapSecondData` を返す（その領域は ΔT 経由
//! = ISSUE-007 で扱う）。UT1/TDB は EOP/ΔT が必要なため本モジュールでは扱わない。

use crate::calendar::{gregorian_to_jd2, jd2_to_gregorian};
use crate::constants::TT_MINUS_TAI_SECONDS;
use crate::error::{DomainError, TimeError};
use crate::julian::JulianDate2;

const SECONDS_PER_DAY: f64 = 86_400.0;

/// 閏秒テーブル `(year, month, day, TAI−UTC[s])`（各 0h UTC で発効）。IERS, 1972–。
const LEAP_SECONDS: &[(i32, u8, u8, f64)] = &[
    (1972, 1, 1, 10.0),
    (1972, 7, 1, 11.0),
    (1973, 1, 1, 12.0),
    (1974, 1, 1, 13.0),
    (1975, 1, 1, 14.0),
    (1976, 1, 1, 15.0),
    (1977, 1, 1, 16.0),
    (1978, 1, 1, 17.0),
    (1979, 1, 1, 18.0),
    (1980, 1, 1, 19.0),
    (1981, 7, 1, 20.0),
    (1982, 7, 1, 21.0),
    (1983, 7, 1, 22.0),
    (1985, 7, 1, 23.0),
    (1988, 1, 1, 24.0),
    (1990, 1, 1, 25.0),
    (1991, 1, 1, 26.0),
    (1992, 7, 1, 27.0),
    (1993, 7, 1, 28.0),
    (1994, 7, 1, 29.0),
    (1996, 1, 1, 30.0),
    (1997, 7, 1, 31.0),
    (1999, 1, 1, 32.0),
    (2006, 1, 1, 33.0),
    (2009, 1, 1, 34.0),
    (2012, 7, 1, 35.0),
    (2015, 7, 1, 36.0),
    (2017, 1, 1, 37.0),
];

fn leap_threshold_jd(entry: &(i32, u8, u8, f64)) -> f64 {
    gregorian_to_jd2(entry.0, entry.1, entry.2, 0, 0, 0.0)
        .expect("leap-second table dates are valid")
        .jd()
}

/// 協定世界時 UTC の瞬時。
#[derive(Clone, Copy, Debug, PartialEq, PartialOrd)]
pub struct UtcInstant(JulianDate2);
/// 国際原子時 TAI の瞬時。
#[derive(Clone, Copy, Debug, PartialEq, PartialOrd)]
pub struct TaiInstant(JulianDate2);
/// 地球時 TT の瞬時。
#[derive(Clone, Copy, Debug, PartialEq, PartialOrd)]
pub struct TtInstant(JulianDate2);

impl UtcInstant {
    /// UTC スケールの JD から構築。
    pub fn from_jd2(jd: JulianDate2) -> Self {
        UtcInstant(jd)
    }
    /// グレゴリオ暦（UTC）から構築。
    pub fn from_gregorian(
        year: i32,
        month: u8,
        day: u8,
        hour: u8,
        minute: u8,
        second: f64,
    ) -> Result<Self, DomainError> {
        Ok(UtcInstant(gregorian_to_jd2(
            year, month, day, hour, minute, second,
        )?))
    }
    /// UTC スケールの JD。
    pub fn jd2(self) -> JulianDate2 {
        self.0
    }
    /// グレゴリオ暦（UTC）へ。
    pub fn to_gregorian(self) -> (i32, u8, u8, u8, u8, f64) {
        jd2_to_gregorian(self.0)
    }
}

impl TaiInstant {
    /// TAI スケールの JD から構築。
    pub fn from_jd2(jd: JulianDate2) -> Self {
        TaiInstant(jd)
    }
    /// TAI スケールの JD。
    pub fn jd2(self) -> JulianDate2 {
        self.0
    }
}

impl TtInstant {
    /// TT スケールの JD から構築。
    pub fn from_jd2(jd: JulianDate2) -> Self {
        TtInstant(jd)
    }
    /// TT スケールの JD。
    pub fn jd2(self) -> JulianDate2 {
        self.0
    }
}

/// 世界時 UT1 の瞬時（地球回転。TT − ΔT。ΔT は `crate::deltat` 参照）。
#[derive(Clone, Copy, Debug, PartialEq, PartialOrd)]
pub struct Ut1Instant(JulianDate2);

impl Ut1Instant {
    /// UT1 スケールの JD から構築。
    pub fn from_jd2(jd: JulianDate2) -> Self {
        Ut1Instant(jd)
    }
    /// UT1 スケールの JD。
    pub fn jd2(self) -> JulianDate2 {
        self.0
    }
}

/// 太陽系力学時 TDB の瞬時（Reference 暦用。TT との差は周期項で最大 ~1.7 ms。
/// TT↔TDB 変換は精度が要る段階で別途実装する）。
#[derive(Clone, Copy, Debug, PartialEq, PartialOrd)]
pub struct TdbInstant(JulianDate2);

impl TdbInstant {
    /// TDB スケールの JD から構築。
    pub fn from_jd2(jd: JulianDate2) -> Self {
        TdbInstant(jd)
    }
    /// TDB スケールの JD。
    pub fn jd2(self) -> JulianDate2 {
        self.0
    }
}

/// 時刻範囲 `[start, end]`（任意の時刻型に対する区間）。
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct TimeRange<T> {
    /// 開始。
    pub start: T,
    /// 終了。
    pub end: T,
}

/// 時間区間（フィット区間など。`TimeRange` と別用途で使い分ける）。
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct TimeInterval<T> {
    /// 開始。
    pub start: T,
    /// 終了。
    pub end: T,
}

/// ΔAT = TAI − UTC（秒）。1972 年より前は `MissingLeapSecondData`。
///
/// 最終エントリ以降は最後の値を据え置く（閏秒は約半年前に告知され、2017 以降は増えていない）。
pub fn tai_minus_utc(utc: UtcInstant) -> Result<f64, TimeError> {
    let jd = utc.0.jd();
    if jd < leap_threshold_jd(&LEAP_SECONDS[0]) {
        return Err(TimeError::MissingLeapSecondData);
    }
    let mut dat = LEAP_SECONDS[0].3;
    for entry in LEAP_SECONDS {
        if jd >= leap_threshold_jd(entry) {
            dat = entry.3;
        } else {
            break;
        }
    }
    Ok(dat)
}

/// UTC → TAI。
pub fn utc_to_tai(utc: UtcInstant) -> Result<TaiInstant, TimeError> {
    let dat = tai_minus_utc(utc)?;
    Ok(TaiInstant(utc.0.add_days(dat / SECONDS_PER_DAY)))
}

/// TAI → TT（定数 +32.184 s）。
pub fn tai_to_tt(tai: TaiInstant) -> TtInstant {
    TtInstant(tai.0.add_days(TT_MINUS_TAI_SECONDS / SECONDS_PER_DAY))
}

/// TT → TAI（定数 −32.184 s）。
pub fn tt_to_tai(tt: TtInstant) -> TaiInstant {
    TaiInstant(tt.0.add_days(-TT_MINUS_TAI_SECONDS / SECONDS_PER_DAY))
}

/// UTC → TT。
pub fn utc_to_tt(utc: UtcInstant) -> Result<TtInstant, TimeError> {
    Ok(tai_to_tt(utc_to_tai(utc)?))
}

/// TAI → UTC。ΔAT は UTC 依存だが、tai を UTC とみなして ΔAT を引く単純法を用いる
/// （閏秒境界の ±ΔAT 秒以内でのみ 1 秒ずれうる。出力用途では十分。conventions / 要確認）。
pub fn tai_to_utc(tai: TaiInstant) -> Result<UtcInstant, TimeError> {
    let dat = tai_minus_utc(UtcInstant(tai.0))?;
    Ok(UtcInstant(tai.0.add_days(-dat / SECONDS_PER_DAY)))
}

/// TT → UTC。
pub fn tt_to_utc(tt: TtInstant) -> Result<UtcInstant, TimeError> {
    tai_to_utc(tt_to_tai(tt))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn utc(y: i32, mo: u8, d: u8) -> UtcInstant {
        UtcInstant::from_gregorian(y, mo, d, 0, 0, 0.0).unwrap()
    }

    #[test]
    fn delta_at_known_values() {
        assert_eq!(tai_minus_utc(utc(2020, 6, 1)).unwrap(), 37.0);
        assert_eq!(tai_minus_utc(utc(2017, 1, 1)).unwrap(), 37.0);
        assert_eq!(tai_minus_utc(utc(2016, 12, 31)).unwrap(), 36.0);
        assert_eq!(tai_minus_utc(utc(2000, 1, 1)).unwrap(), 32.0);
        assert_eq!(tai_minus_utc(utc(1985, 7, 1)).unwrap(), 23.0);
        assert_eq!(tai_minus_utc(utc(1972, 1, 1)).unwrap(), 10.0);
    }

    #[test]
    fn delta_at_steps_exactly_on_boundary() {
        // 2017-01-01 0h で 36 → 37 へ跳ぶ。境界当日は新しい値。
        assert_eq!(tai_minus_utc(utc(2016, 12, 31)).unwrap(), 36.0);
        assert_eq!(tai_minus_utc(utc(2017, 1, 1)).unwrap(), 37.0);
    }

    #[test]
    fn before_1972_is_missing() {
        let pre = UtcInstant::from_gregorian(1971, 12, 31, 0, 0, 0.0).unwrap();
        assert_eq!(
            tai_minus_utc(pre).unwrap_err(),
            TimeError::MissingLeapSecondData
        );
    }

    #[test]
    fn utc_to_tt_offset_is_dat_plus_32_184() {
        // 2020: ΔAT=37 → TT−UTC = 69.184 s。
        let u = utc(2020, 1, 1);
        let tt = utc_to_tt(u).unwrap();
        let diff_s = tt.jd2().days_since(u.jd2()) * SECONDS_PER_DAY;
        assert!((diff_s - 69.184).abs() < 1e-6, "diff = {diff_s}");
    }

    #[test]
    fn j2000_utc_to_tt_offset_is_64_184() {
        // 2000-01-01 12:00 UTC: ΔAT=32 → TT−UTC = 64.184 s。
        let u = UtcInstant::from_gregorian(2000, 1, 1, 12, 0, 0.0).unwrap();
        let tt = utc_to_tt(u).unwrap();
        let diff_s = tt.jd2().days_since(u.jd2()) * SECONDS_PER_DAY;
        assert!((diff_s - 64.184).abs() < 1e-6, "diff = {diff_s}");
    }

    #[test]
    fn tai_tt_offset_is_exactly_32_184() {
        let tai = TaiInstant::from_jd2(JulianDate2::from_jd(2_460_000.0));
        let tt = tai_to_tt(tai);
        assert!((tt.jd2().days_since(tai.jd2()) * SECONDS_PER_DAY - 32.184).abs() < 1e-9);
    }

    #[test]
    fn utc_instant_gregorian_round_trip() {
        // UtcInstant::from_gregorian → to_gregorian の往復（委譲先 calendar とは別に本型を検証）。
        let u = UtcInstant::from_gregorian(2035, 9, 2, 1, 30, 15.5).unwrap();
        let (y, mo, d, h, mi, s) = u.to_gregorian();
        assert_eq!((y, mo, d, h, mi), (2035, 9, 2, 1, 30));
        assert!((s - 15.5).abs() < 1e-4, "s = {s}");
    }

    #[test]
    fn round_trip_utc_tt_utc() {
        let u = UtcInstant::from_gregorian(2035, 9, 2, 1, 30, 15.0).unwrap();
        let back = tt_to_utc(utc_to_tt(u).unwrap()).unwrap();
        assert!(back.jd2().days_since(u.jd2()).abs() * SECONDS_PER_DAY < 1e-6);
    }
}
