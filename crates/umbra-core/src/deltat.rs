//! ΔT (= TT − UT1) モデルと UT1 変換（`docs/issues/ISSUE-007`、`docs/algorithms/01-time-scales.md`）。
//!
//! [`EspenakMeeusDeltaT`] は Espenak & Meeus（NASA TP-2006-214141, `eclipse.gsfc.nasa.gov`）の
//! **区分多項式**。公開された数式（バルクデータではない）なのでライセンス問題なく組み込める。
//! 厳密な近代 ΔT は EOP（UT1−UTC 実測, 将来の ISSUE-007 完全版）から得るが、本モデルは
//! その予測値・外挿値と不確実性帯（accuracy.md §0）を与える。

use crate::calendar::jd2_to_gregorian;
use crate::julian::JulianDate2;
use crate::time::{TtInstant, Ut1Instant};

const SECONDS_PER_DAY: f64 = 86_400.0;

/// ΔT = TT − UT1 モデル。
pub trait DeltaTModel {
    /// 十進年に対する ΔT（秒）。
    fn delta_t_seconds(&self, decimal_year: f64) -> f64;
    /// ΔT の不確実性（秒, 1σ 目安。accuracy.md §0 の不確実性帯）。
    fn uncertainty_seconds(&self, decimal_year: f64) -> f64;
}

/// Espenak & Meeus の ΔT 区分多項式（NASA TP-2006-214141）。
#[derive(Debug, Clone, Copy, Default)]
pub struct EspenakMeeusDeltaT;

/// 1820 を基準とした長期放物線 ΔT = −20 + 32·u²（u = (y−1820)/100）。範囲外の既定。
fn long_term(y: f64) -> f64 {
    let u = (y - 1820.0) / 100.0;
    -20.0 + 32.0 * u * u
}

impl DeltaTModel for EspenakMeeusDeltaT {
    fn delta_t_seconds(&self, y: f64) -> f64 {
        if y < 1900.0 {
            long_term(y)
        } else if y < 1920.0 {
            let t = y - 1900.0;
            -2.79 + 1.494119 * t - 0.0598939 * t * t + 0.0061966 * t * t * t
                - 0.000197 * t * t * t * t
        } else if y < 1941.0 {
            let t = y - 1920.0;
            21.20 + 0.84493 * t - 0.076100 * t * t + 0.0020936 * t * t * t
        } else if y < 1961.0 {
            let t = y - 1950.0;
            29.07 + 0.407 * t - t * t / 233.0 + t * t * t / 2547.0
        } else if y < 1986.0 {
            let t = y - 1975.0;
            45.45 + 1.067 * t - t * t / 260.0 - t * t * t / 718.0
        } else if y < 2005.0 {
            let t = y - 2000.0;
            63.86 + 0.3345 * t - 0.060374 * t * t
                + 0.0017275 * t * t * t
                + 0.000651814 * t * t * t * t
                + 0.00002373599 * t * t * t * t * t
        } else if y < 2050.0 {
            let t = y - 2000.0;
            62.92 + 0.32217 * t + 0.005589 * t * t
        } else if y < 2150.0 {
            -20.0 + 32.0 * ((y - 1820.0) / 100.0) * ((y - 1820.0) / 100.0) - 0.5628 * (2150.0 - y)
        } else {
            long_term(y)
        }
    }

    fn uncertainty_seconds(&self, y: f64) -> f64 {
        // 粗い目安（要確認）。実測 EOP 域では UT1−UTC から精密に得るべきで、本値は
        // モデル予測の帯。順序しきい値（不連続）にして各境界・式を検証可能にする。
        if y < 1900.0 {
            5.0 + 0.1 * (1900.0 - y) // 古い年代ほど増大
        } else if y < 1955.0 {
            2.0
        } else if y < 2006.0 {
            0.5 // 近代（おおむね実測 EOP 域）
        } else {
            1.0 + 0.5 * (y - 2006.0) // 2006 以降の外挿は年々増大
        }
    }
}

/// 十進年 `year + (month − 0.5)/12`（Espenak 慣習）を JD から求める。
pub fn decimal_year(jd: JulianDate2) -> f64 {
    let (year, month, ..) = jd2_to_gregorian(jd);
    f64::from(year) + (f64::from(month) - 0.5) / 12.0
}

/// TT → UT1（`UT1 = TT − ΔT`）。
pub fn tt_to_ut1<M: DeltaTModel>(tt: TtInstant, model: &M) -> Ut1Instant {
    let dt = model.delta_t_seconds(decimal_year(tt.jd2()));
    Ut1Instant::from_jd2(tt.jd2().add_days(-dt / SECONDS_PER_DAY))
}

/// UT1 → TT（`TT = UT1 + ΔT`）。ΔT は UT1 の年から評価する（TT 年との差は ΔT に無視可能）。
pub fn ut1_to_tt<M: DeltaTModel>(ut1: Ut1Instant, model: &M) -> TtInstant {
    let dt = model.delta_t_seconds(decimal_year(ut1.jd2()));
    TtInstant::from_jd2(ut1.jd2().add_days(dt / SECONDS_PER_DAY))
}

#[cfg(test)]
mod tests {
    use super::*;

    const DT: EspenakMeeusDeltaT = EspenakMeeusDeltaT;

    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    #[test]
    fn piece_values_match_formulas_at_interior_points() {
        // 各区分の内部点（非ゼロ t）で手計算オラクルと一致 → 区分内の演算子取り違えを検出。
        assert!(
            close(DT.delta_t_seconds(1910.0), 10.3884, 1e-3),
            "{}",
            DT.delta_t_seconds(1910.0)
        );
        assert!(
            close(DT.delta_t_seconds(1930.0), 24.1329, 1e-3),
            "{}",
            DT.delta_t_seconds(1930.0)
        );
        assert!(
            close(DT.delta_t_seconds(1955.0), 31.0468, 1e-3),
            "{}",
            DT.delta_t_seconds(1955.0)
        );
        assert!(
            close(DT.delta_t_seconds(1980.0), 50.5148, 1e-3),
            "{}",
            DT.delta_t_seconds(1980.0)
        );
        assert!(
            close(DT.delta_t_seconds(1995.0), 60.7954, 1e-3),
            "{}",
            DT.delta_t_seconds(1995.0)
        );
        assert!(
            close(DT.delta_t_seconds(2030.0), 77.6151, 1e-3),
            "{}",
            DT.delta_t_seconds(2030.0)
        );
        assert!(
            close(DT.delta_t_seconds(2100.0), 202.74, 1e-2),
            "{}",
            DT.delta_t_seconds(2100.0)
        );
        assert!(
            close(DT.delta_t_seconds(2200.0), 442.08, 1e-2),
            "{}",
            DT.delta_t_seconds(2200.0)
        );
        assert!(close(DT.delta_t_seconds(1850.0), long_term(1850.0), 1e-9));
    }

    #[test]
    fn piece_starts_match_constants_at_boundaries() {
        // 各区分の開始年（境界）で、その区分の定数項に一致 → 境界比較の向き（< vs <=）を検出。
        assert!(close(DT.delta_t_seconds(1900.0), -2.79, 1e-9));
        assert!(close(DT.delta_t_seconds(1920.0), 21.20, 1e-9));
        assert!(close(DT.delta_t_seconds(1950.0), 29.07, 1e-9)); // 1941–1961 区分 t=0
        assert!(close(DT.delta_t_seconds(1975.0), 45.45, 1e-9)); // 1961–1986 区分 t=0
        assert!(close(DT.delta_t_seconds(2000.0), 63.86, 1e-9)); // 1986–2005 区分 t=0
    }

    #[test]
    fn delta_t_2000_matches_known_value() {
        // ΔT(2000.0) ≈ 63.8 s（観測既知）に近い。
        assert!(close(DT.delta_t_seconds(2000.0), 63.8, 0.1));
    }

    #[test]
    fn uncertainty_branches_boundaries_and_growth() {
        // 各分岐の内部値（式の演算子取り違えを検出）。
        assert!(close(DT.uncertainty_seconds(1980.0), 0.5, 1e-12));
        assert!(close(DT.uncertainty_seconds(1930.0), 2.0, 1e-12));
        assert!(close(
            DT.uncertainty_seconds(1800.0),
            5.0 + 0.1 * 100.0,
            1e-9
        )); // 15.0
        assert!(close(
            DT.uncertainty_seconds(2106.0),
            1.0 + 0.5 * 100.0,
            1e-9
        )); // 51.0
            // 境界（不連続なので比較の向き < / <= を検出）。
        assert!(close(DT.uncertainty_seconds(1900.0), 2.0, 1e-12)); // y<1900 false
        assert!(close(DT.uncertainty_seconds(1955.0), 0.5, 1e-12)); // y<1955 false
        assert!(close(DT.uncertainty_seconds(2006.0), 1.0, 1e-12)); // y<2006 false
                                                                    // 将来へ増大。
        assert!(DT.uncertainty_seconds(2100.0) > DT.uncertainty_seconds(2030.0));
    }

    #[test]
    fn decimal_year_uses_mid_month_convention() {
        let jd = crate::calendar::gregorian_to_jd2(2000, 7, 2, 0, 0, 0.0).unwrap();
        // 7月 → 2000 + (7-0.5)/12 = 2000.5417
        assert!(close(decimal_year(jd), 2000.0 + 6.5 / 12.0, 1e-9));
    }

    #[test]
    fn tt_to_ut1_subtracts_delta_t() {
        let tt =
            TtInstant::from_jd2(crate::calendar::gregorian_to_jd2(2010, 1, 1, 0, 0, 0.0).unwrap());
        let ut1 = tt_to_ut1(tt, &DT);
        let dt = DT.delta_t_seconds(decimal_year(tt.jd2()));
        let diff_s = ut1.jd2().days_since(tt.jd2()) * SECONDS_PER_DAY;
        assert!(close(diff_s, -dt, 1e-6), "diff = {diff_s}, dt = {dt}");
        assert!(dt > 60.0 && dt < 75.0, "dt(2010) = {dt}");
    }

    #[test]
    fn ut1_tt_round_trip() {
        let tt =
            TtInstant::from_jd2(crate::calendar::gregorian_to_jd2(2035, 9, 2, 1, 30, 0.0).unwrap());
        let back = ut1_to_tt(tt_to_ut1(tt, &DT), &DT);
        assert!(back.jd2().days_since(tt.jd2()).abs() * SECONDS_PER_DAY < 1e-3);
    }
}
