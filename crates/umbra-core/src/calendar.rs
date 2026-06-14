//! グレゴリオ暦 ⇔ ユリウス日（`docs/issues/ISSUE-005`、`docs/algorithms/01-time-scales.md`）。
//!
//! 時刻系に依存しない暦変換。全期間を**プロレプティック・グレゴリオ暦**で一貫させる
//! （1582年以前もユリウス暦へ切替えない。conventions §3 / 要確認: NASA カタログのユリウス暦慣習との差）。
//! アルゴリズムは Meeus *Astronomical Algorithms* 2nd ed. Ch.7（Gregorian 分岐を常用）。

use crate::error::DomainError;
use crate::julian::JulianDate2;

/// グレゴリオ暦（年・月・日・時・分・秒）→ [`JulianDate2`]。
///
/// 年は負も可（プロレプティック）。月 1–12・日 1–31・時 0–23・分 0–59・秒 \[0, 60) を検証する。
pub fn gregorian_to_jd2(
    year: i32,
    month: u8,
    day: u8,
    hour: u8,
    minute: u8,
    second: f64,
) -> Result<JulianDate2, DomainError> {
    if !(1..=12).contains(&month) {
        return Err(DomainError::OutOfRange { what: "month" });
    }
    if !(1..=31).contains(&day) {
        return Err(DomainError::OutOfRange { what: "day" });
    }
    if hour >= 24 {
        return Err(DomainError::OutOfRange { what: "hour" });
    }
    if minute >= 60 {
        return Err(DomainError::OutOfRange { what: "minute" });
    }
    if !(0.0..60.0).contains(&second) {
        return Err(DomainError::OutOfRange { what: "second" });
    }

    // 1〜2月は前年の13〜14月として扱う（Meeus 7）。
    let (y, m) = if month <= 2 {
        (year - 1, i32::from(month) + 12)
    } else {
        (year, i32::from(month))
    };
    let a = (f64::from(y) / 100.0).floor();
    let b = 2.0 - a + (a / 4.0).floor(); // 常に Gregorian 補正（プロレプティック）
                                         // 0h（深夜）の JD。day は整数。
    let jd0 = (365.25 * (f64::from(y) + 4716.0)).floor()
        + (30.6001 * (f64::from(m) + 1.0)).floor()
        + f64::from(day)
        + b
        - 1524.5;
    let day_fraction = f64::from(hour) / 24.0 + f64::from(minute) / 1440.0 + second / 86400.0;
    Ok(JulianDate2::new(jd0, day_fraction).normalized())
}

/// [`JulianDate2`] → グレゴリオ暦 `(year, month, day, hour, minute, second)`。
///
/// Meeus Ch.7（Gregorian 分岐を常用）。
// 出力は有界な暦フィールド（年は i32、月日時分は小範囲）。値域は floor 演算で保証されるため、
// f64→整数の切り詰め/符号落ちは意図的であり安全。
#[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
pub fn jd2_to_gregorian(jd2: JulianDate2) -> (i32, u8, u8, u8, u8, f64) {
    let jd = jd2.jd() + 0.5;
    let z = jd.floor();
    let f = jd - z;

    // プロレプティック・グレゴリオ: z の大小に関わらず Gregorian 補正を適用。
    let alpha = ((z - 1_867_216.25) / 36_524.25).floor();
    let a = z + 1.0 + alpha - (alpha / 4.0).floor();
    let b = a + 1524.0;
    let c = ((b - 122.1) / 365.25).floor();
    let d = (365.25 * c).floor();
    let e = ((b - d) / 30.6001).floor();

    let day_with_frac = b - d - (30.6001 * e).floor() + f;
    let day = day_with_frac.floor();
    let month_f = if e < 14.0 { e - 1.0 } else { e - 13.0 };
    let year_f = if month_f > 2.0 {
        c - 4716.0
    } else {
        c - 4715.0
    };

    let mut rem_seconds = (day_with_frac - day) * 86400.0;
    let hour = (rem_seconds / 3600.0).floor();
    rem_seconds -= hour * 3600.0;
    let minute = (rem_seconds / 60.0).floor();
    let second = rem_seconds - minute * 60.0;

    (
        year_f as i32,
        month_f as u8,
        day as u8,
        hour as u8,
        minute as u8,
        second,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn j2000_noon() {
        let jd = gregorian_to_jd2(2000, 1, 1, 12, 0, 0.0).unwrap();
        assert!((jd.jd() - 2_451_545.0).abs() < 1e-9);
    }

    #[test]
    fn midnight_is_half_less() {
        let jd = gregorian_to_jd2(2000, 1, 1, 0, 0, 0.0).unwrap();
        assert!((jd.jd() - 2_451_544.5).abs() < 1e-9);
    }

    #[test]
    fn meeus_example_1987() {
        // Meeus Ch.7: 1987-01-27.0 → JD 2446822.5
        let jd = gregorian_to_jd2(1987, 1, 27, 0, 0, 0.0).unwrap();
        assert!((jd.jd() - 2_446_822.5).abs() < 1e-9);
    }

    #[test]
    fn round_trips_modern_and_future() {
        for (y, mo, d, h, mi, s) in [
            (2000, 1, 1, 12, 0, 0.0),
            (2024, 2, 29, 12, 0, 0.0), // 閏日（e≥14 分岐 / 月=2 境界）
            (2000, 2, 15, 8, 20, 5.0), // 2月（month_f>2 の境界を踏む）
            (2035, 9, 2, 1, 30, 15.5),
            (1900, 12, 31, 23, 59, 59.0),
            (-1000, 6, 15, 6, 0, 0.0), // プロレプティック（負の年）
        ] {
            let jd = gregorian_to_jd2(y, mo, d, h, mi, s).unwrap();
            let (ry, rmo, rd, rh, rmi, rs) = jd2_to_gregorian(jd);
            assert_eq!(
                (ry, rmo, rd, rh, rmi),
                (y, mo, d, h, mi),
                "date {y}-{mo}-{d}"
            );
            assert!((rs - s).abs() < 1e-4, "seconds {rs} vs {s}");
        }
    }

    #[test]
    fn rejects_invalid_fields() {
        assert!(gregorian_to_jd2(2000, 13, 1, 0, 0, 0.0).is_err());
        assert!(gregorian_to_jd2(2000, 1, 32, 0, 0, 0.0).is_err());
        assert!(gregorian_to_jd2(2000, 1, 1, 24, 0, 0.0).is_err());
        assert!(gregorian_to_jd2(2000, 1, 1, 0, 0, 60.0).is_err());
    }
}
