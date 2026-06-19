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

/// [`JulianDate2`] → グレゴリオ暦で、**秒を 0.1 秒へ丸めた**表示用の暦成分を返す。
///
/// [`jd2_to_gregorian`] は JD 往復の ±eps により、暦上ちょうど `hh:mm:00` の瞬間を
/// `hh:(mm-1):59.9995` のように返すことがある（例: `16:00:00` → `15:59:59.9995`）。これを
/// `{:.1}` で素通しすると `15:59:60.0` のような不正表記になる。本関数は秒を 0.1 秒へ丸め、
/// 60.0 到達分を 分→時→日 へ繰り上げる（日跨ぎは整数日 +1 で年月日を再導出）。返す秒は必ず
/// `[0.0, 60.0)` の 0.1 刻み（0.0, 0.1, …, 59.9）で、`60.0` 以上にはならない。表示専用で、
/// 完全精度が要る用途は [`JulianDate2`] を直接使うこと（本関数は丸めにより情報を落とす）。
// 秒→tenths（i64）と tenths→f64 は値域が小さく（0..=600）安全。暦フィールドは jd2_to_gregorian
// が有界化済み。意図的な丸め用変換のため lint を許容する。
#[allow(clippy::cast_possible_truncation, clippy::cast_precision_loss)]
pub fn jd2_to_gregorian_deciseconds(jd: JulianDate2) -> (i32, u8, u8, u8, u8, f64) {
    let (mut year, mut month, mut day, mut hour, mut minute, second) = jd2_to_gregorian(jd);
    // 秒を 1/10 秒（整数 tenths, 0..=600）へ丸める。second ∈ [0,60) なので tenths ∈ 0..=600。
    let mut tenths = (second * 10.0).round() as i64;
    if tenths >= 600 {
        tenths -= 600;
        minute += 1;
    }
    if minute >= 60 {
        minute -= 60;
        hour += 1;
    }
    if hour >= 24 {
        hour -= 24;
        // 日跨ぎ: 整数日を 1 進めて年月日のみ採用（時刻は繰り上げ済みの 00:00:00.x）。
        let (next_year, next_month, next_day, ..) =
            jd2_to_gregorian(JulianDate2::new(jd.part1 + 1.0, jd.part2));
        year = next_year;
        month = next_month;
        day = next_day;
    }
    let second = (tenths as f64) / 10.0;
    (year, month, day, hour, minute, second)
}

#[cfg(test)]
mod tests {
    // deciseconds テストは秒を 0.1 秒グリッド整数（tenths）へ `round() as i64` して比較する
    // （`round()` の f64 等値比較は float_cmp を踏むため整数化が適切）。値域は 0..=600 で安全。
    #![allow(clippy::cast_possible_truncation)]

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

    // ==================================================================
    // jd2_to_gregorian_deciseconds（0.1 秒丸め＋桁上げ。:60.0 回帰防止）
    // ==================================================================
    //
    // ## オラクル戦略（丸め方向に依存しない）
    // 秒が 0.1 秒グリッド上にある civil 時刻を gregorian_to_jd2 で JD に往復させ、
    // jd2_to_gregorian_deciseconds が**元の civil 成分**を厳密復元することを縛る
    // （往復誤差 < 0.05 秒なので最近接 0.1 秒丸めで元成分に戻る）。素通し／無丸めだと
    // 16:00:00 / 深夜 が 59:59.9995 として現れて落ちる。

    /// 【THE 回帰】既知ドリフト分境界 2024-04-08 16:00:00.0。
    /// gregorian_to_jd2 で往復すると jd2_to_gregorian は 15:59:59.9995 を返す（実測）。
    /// jd2_to_gregorian_deciseconds は最近接 0.1 秒へ丸め桁上げし (2024,4,8,16,0,0.0) を復元する。
    /// 殺す変異: 丸めなし／jd2_to_gregorian 素通し（15,59,59.9995 を返す）、桁上げを分へ伝播しない。
    #[test]
    fn deciseconds_known_drift_whole_minute_recovers_16_00_00() {
        let jd = gregorian_to_jd2(2024, 4, 8, 16, 0, 0.0).unwrap();
        let (y, mo, d, h, mi, s) = jd2_to_gregorian_deciseconds(jd);
        assert_eq!(
            (y, mo, d, h, mi),
            (2024, 4, 8, 16, 0),
            "暦日時分が復元される"
        );
        assert_eq!(
            (s * 10.0).round() as i64,
            0,
            "秒は 0.0（59.9995 でない）: {s}"
        );
    }

    /// 通常の分境界 2024-04-08 18:17:00.0 も丸めで分・秒が正しく復元される。
    /// 殺す変異: 丸めなし（往復 eps で 16〜17 秒や 59.999 秒が漏れる）、秒のみ丸め桁上げ非伝播。
    #[test]
    fn deciseconds_plain_whole_minute_recovers_18_17_00() {
        let jd = gregorian_to_jd2(2024, 4, 8, 18, 17, 0.0).unwrap();
        let (y, mo, d, h, mi, s) = jd2_to_gregorian_deciseconds(jd);
        assert_eq!(
            (y, mo, d, h, mi),
            (2024, 4, 8, 18, 17),
            "暦日時分が復元される"
        );
        assert_eq!((s * 10.0).round() as i64, 0, "秒は 0.0: {s}");
    }

    /// グリッド上の端数秒 2024-04-08 10:20:30.7 が 0.1 秒グリッド上で復元される（307 = 30.7 秒）。
    /// 殺す変異: 端数秒を 0 へ潰す（整数秒丸め）、最近接でなく floor/ceil する、秒成分を誤配線。
    #[test]
    fn deciseconds_fractional_on_grid_recovers_30_7() {
        let jd = gregorian_to_jd2(2024, 4, 8, 10, 20, 30.7).unwrap();
        let (y, mo, d, h, mi, s) = jd2_to_gregorian_deciseconds(jd);
        assert_eq!(
            (y, mo, d, h, mi),
            (2024, 4, 8, 10, 20),
            "暦日時分が復元される"
        );
        assert_eq!(
            (s * 10.0).round() as i64,
            307,
            "秒は 30.7（0.1 秒グリッド）: {s}"
        );
    }

    /// 深夜 2024-04-09 00:00:00.0（往復で 04-08 23:59:59.9995 になりうる）。
    /// 丸め桁上げが 秒→分→時→日 を貫通し (2024,4,9,0,0,0.0) を復元する。
    /// 殺す変異: 時→日 桁上げ破壊（04-08 23:59:60 や 04-09 24:00 を返す）、丸めなし。
    #[test]
    fn deciseconds_midnight_carries_hour_to_day() {
        let jd = gregorian_to_jd2(2024, 4, 9, 0, 0, 0.0).unwrap();
        let (y, mo, d, h, mi, s) = jd2_to_gregorian_deciseconds(jd);
        assert_eq!(
            (y, mo, d, h, mi),
            (2024, 4, 9, 0, 0),
            "深夜が翌日 0:00 へ桁上げ"
        );
        assert_eq!((s * 10.0).round() as i64, 0, "秒は 0.0: {s}");
    }

    /// 年境界 2025-01-01 00:00:00.0（往復で 2024-12-31 23:59:59.9995 になりうる）。
    /// 丸め桁上げが 日→月→年 のロールオーバーを貫通し (2025,1,1,0,0,0.0) を復元する。
    /// 殺す変異: 日→月→年 ロールオーバー破壊（2024-12-31 23:59:60 や 2024-13-01 を返す）、丸めなし。
    #[test]
    fn deciseconds_year_boundary_rolls_day_month_year() {
        let jd = gregorian_to_jd2(2025, 1, 1, 0, 0, 0.0).unwrap();
        let (y, mo, d, h, mi, s) = jd2_to_gregorian_deciseconds(jd);
        assert_eq!(
            (y, mo, d, h, mi),
            (2025, 1, 1, 0, 0),
            "年末が翌年 1/1 0:00 へ桁上げ"
        );
        assert_eq!((s * 10.0).round() as i64, 0, "秒は 0.0: {s}");
    }

    /// 不変条件: 返り値の秒は常に [0.0, 60.0) かつ 0.1 秒グリッド上（0.1 の倍数）。NEVER ≥ 60.0。
    /// 16:00:00 / 深夜 などドリフト境界を含む複数入力で確認する。
    /// 殺す変異: 60.0 を返す（桁上げ漏れ）、生の秒（59.9995 等・非グリッド）を返す。
    #[test]
    fn deciseconds_second_always_in_grid_below_sixty() {
        for (y, mo, d, h, mi, sec) in [
            (2024, 4, 8, 16, 0, 0.0),   // ドリフト境界
            (2024, 4, 9, 0, 0, 0.0),    // 深夜
            (2025, 1, 1, 0, 0, 0.0),    // 年境界
            (2024, 4, 8, 18, 17, 0.0),  // 通常分境界
            (2024, 4, 8, 10, 20, 30.7), // グリッド端数
            (2024, 4, 8, 23, 59, 59.9), // 分末端のグリッド秒
        ] {
            let jd = gregorian_to_jd2(y, mo, d, h, mi, sec).unwrap();
            let (_, _, _, _, _, s) = jd2_to_gregorian_deciseconds(jd);
            assert!(
                (0.0..60.0).contains(&s),
                "秒は [0.0, 60.0)（60.0 を返さない）: got {s} for {y}-{mo}-{d} {h}:{mi}:{sec}"
            );
            assert!(
                ((s * 10.0).round() - s * 10.0).abs() < 1e-6,
                "秒は 0.1 秒グリッド上（生の秒でない）: got {s} for {y}-{mo}-{d} {h}:{mi}:{sec}"
            );
        }
    }
}
