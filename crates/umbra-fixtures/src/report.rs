//! 誤差統計と許容プロファイル（`docs/issues/ISSUE-030` S30a）。
//!
//! 検証レポータの**純粋な基盤**: 誤差列の記述統計 [`ErrorStats`] と、合否判定に使う
//! 許容値 [`ToleranceProfile`]（`docs/accuracy.md` §2）を提供する。ゴールデン比較・層分解・
//! JPL DE 差分・1900〜2100 全走査・JSON/CLI 出力は後続スライス（本モジュールは比較の前段＝
//! 統計と許容の定義に徹する）。
//!
//! 設計規律（conventions §11 / accuracy.md §4）: 統計は**誤差を隠さない**。pass 判定が通っても
//! 数値（max/mean/p95）は必ず保持する。許容を pass のために拡大しない。

use umbra_eclipse::SolarEclipse;

use crate::types::GoldenEclipse;

/// 1 日の秒数（JD 差 → 秒の換算）。
const SECONDS_PER_DAY: f64 = 86_400.0;

/// 1 項目の誤差記述統計（**絶対誤差**ベース, accuracy.md §3.4）。
///
/// 各 metric（接触秒・最大食秒・食分・食面積・高度…）の誤差列から、最大絶対誤差・平均絶対誤差・
/// 95 パーセンタイルを保持する。`units` は表示・取り違え防止のための単位識別子（例 `"s"`/`"deg"`）。
///
/// percentile は **R-7（線形補間, NumPy 既定）** に固定する（ISSUE-030 §65 の「規約を 1 つに固定」）:
/// 昇順絶対誤差 `a` と分位 `p` に対し `h = (n-1)·p`、`lo = ⌊h⌋`、
/// `a[lo] + (h-lo)·(a[lo+1]-a[lo])`（端は `a[n-1]`）。`n==1` は単一値、`n==0` は 0.0。
#[derive(Clone, Debug, PartialEq)]
pub struct ErrorStats {
    /// 標本数（誤差列の長さ）。
    pub n: usize,
    /// 最大絶対誤差 `max|e|`。
    pub max_abs: f64,
    /// 平均絶対誤差 `mean|e|`。
    pub mean_abs: f64,
    /// 絶対誤差の 95 パーセンタイル（R-7 線形補間）。
    pub p95: f64,
    /// 単位識別子（`"s"` / `"deg"` / 無次元は `""` 等）。
    pub units: &'static str,
}

impl ErrorStats {
    /// 誤差列（**符号付き可**）から絶対誤差の記述統計を作る。
    ///
    /// 各要素を絶対値化してから `max|e|` / `mean|e|` / `p95(|e|)` を求める（符号は相殺しない）。
    /// 空列は `n=0`・`max_abs=mean_abs=p95=0.0`（`units` は保持・NaN を出さない）。
    #[allow(clippy::cast_precision_loss)]
    pub fn from_errors(errors: &[f64], units: &'static str) -> Self {
        let n = errors.len();
        if n == 0 {
            return Self {
                n: 0,
                max_abs: 0.0,
                mean_abs: 0.0,
                p95: 0.0,
                units,
            };
        }
        let mut abs: Vec<f64> = errors.iter().map(|e| e.abs()).collect();
        let sum: f64 = abs.iter().sum();
        let mean_abs = sum / n as f64;
        // 全 [0,∞) の有限値（入力は有限前提）。total_cmp で NaN 非依存・unwrap なしに昇順化。
        abs.sort_by(f64::total_cmp);
        let max_abs = abs[n - 1];
        let p95 = percentile_r7_sorted(&abs, 0.95);
        Self {
            n,
            max_abs,
            mean_abs,
            p95,
            units,
        }
    }

    /// 最大絶対誤差が許容以下か（pass 判定）。`max_abs <= tolerance`（境界は pass・inclusive）。
    pub fn within(&self, tolerance: f64) -> bool {
        self.max_abs <= tolerance
    }
}

/// 昇順済みスライスの分位（R-7 線形補間）。`p ∈ [0,1]`。空は 0.0、単一要素はその値。
#[allow(
    clippy::cast_precision_loss,
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss
)]
fn percentile_r7_sorted(sorted_abs: &[f64], p: f64) -> f64 {
    let n = sorted_abs.len();
    if n == 0 {
        return 0.0;
    }
    if n == 1 {
        return sorted_abs[0];
    }
    let h = (n - 1) as f64 * p;
    let lo = h.floor() as usize;
    if lo + 1 >= n {
        sorted_abs[n - 1]
    } else {
        sorted_abs[lo] + (h - lo as f64) * (sorted_abs[lo + 1] - sorted_abs[lo])
    }
}

/// 合否判定に使う許容プロファイル（plan §18・モデル別, accuracy.md §2）。
///
/// `note_utc_is_delta_t_limited` は「**UTC 絶対は ΔT/UT1 予測律速**（accuracy.md §0(b)/§2.3）」を
/// 表すフラグ。UTC 列の判定は不確実性帯を考慮する必要がある（計算誤差と混同しない）。値は TT 基準の
/// 幾何許容（接触秒・最大食秒）と無次元/度の許容。
#[derive(Clone, Debug, PartialEq)]
pub struct ToleranceProfile {
    /// 局地接触時刻の許容（秒・TT 基準幾何。accuracy.md §2.1L 目標 ±2 s）。
    pub contact_seconds: f64,
    /// 最大食時刻の許容（秒・TT 基準幾何。accuracy.md §2.1 目標 ±1.5 s）。
    pub maximum_seconds: f64,
    /// 食分の許容（無次元。accuracy.md §2.2 ±0.0005）。
    pub magnitude: f64,
    /// 食面積の許容（無次元。accuracy.md §2.2 相当）。
    pub obscuration: f64,
    /// 太陽高度の許容（度・表示精度）。
    pub altitude_degrees: f64,
    /// UTC 絶対が ΔT/UT1 予測律速であることの注記フラグ（accuracy.md §0(b)/§2.3）。
    pub note_utc_is_delta_t_limited: bool,
}

impl ToleranceProfile {
    /// 本番標準プロファイル（accuracy.md §2 の設計目標）。
    pub fn standard() -> Self {
        Self {
            contact_seconds: 2.0,
            maximum_seconds: 1.5,
            magnitude: 0.0005,
            obscuration: 0.0005,
            altitude_degrees: 0.1,
            note_utc_is_delta_t_limited: true,
        }
    }

    /// 高精度参照プロファイル（回帰・差分テストの第一義オラクル。standard を厳格化）。
    pub fn reference() -> Self {
        Self {
            contact_seconds: 1.0,
            maximum_seconds: 1.0,
            magnitude: 0.0002,
            obscuration: 0.0002,
            altitude_degrees: 0.05,
            note_utc_is_delta_t_limited: true,
        }
    }
}

/// 1 日食の全球条件の誤差（**符号付き = computed − golden**, accuracy.md §3.4）。
///
/// 最大食時刻は **TT 基準の幾何誤差**（golden が TT を持てば TT 差、無ければ UTC 差で代替）の秒。
/// γ は Re、食分は無次元。符号は computed が後/大なら正。集計時に [`ErrorStats`] が絶対値化する。
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GlobalErrors {
    /// 最大食時刻誤差 \[s\]（computed − golden）。
    pub greatest_seconds: f64,
    /// γ 誤差 \[Re\]（computed − golden）。
    pub gamma: f64,
    /// 食分誤差（無次元, computed − golden）。
    pub magnitude: f64,
}

/// computed の全球条件を golden と比較し、符号付き誤差（computed − golden）を返す（純粋）。
///
/// 最大食時刻は golden が TT(=TD) を持てば **TT 差**（純幾何・ΔT 非依存, accuracy.md §0(a)）、
/// 持たなければ UTC 差で代替する（その場合は ΔT 律速の注記対象, §0(b)）。いずれも秒。
pub fn compare_global(computed: &SolarEclipse, golden: &GoldenEclipse) -> GlobalErrors {
    let greatest = &computed.global.greatest;
    // 差分は JulianDate2 の 2 要素を保った days_since で取る（単一 f64 の jd() 同士を引くと
    // JD≈2.45e6 でエポック減算の桁落ち（~4.6e-5 s）が出るため。julian §2要素表現の理由）。
    let greatest_seconds = match golden.greatest_time_tt {
        Some(golden_tt) => greatest.time_tt.jd2().days_since(golden_tt.jd2()) * SECONDS_PER_DAY,
        None => {
            greatest
                .time_utc
                .jd2()
                .days_since(golden.greatest_time_utc.jd2())
                * SECONDS_PER_DAY
        }
    };
    GlobalErrors {
        greatest_seconds,
        gamma: computed.global.gamma - golden.gamma,
        magnitude: greatest.magnitude.0 - golden.magnitude,
    }
}

/// 全球比較の metric 別統計＋合否（accuracy.md §3.4: pass でも統計を必ず出す）。
#[derive(Clone, Debug, PartialEq)]
pub struct GlobalReport {
    /// 最大食時刻誤差の統計（単位 `"s"`）。
    pub greatest: ErrorStats,
    /// γ 誤差の統計（単位 `"Re"`）。許容未設定のため合否ゲートには使わない（統計のみ）。
    pub gamma: ErrorStats,
    /// 食分誤差の統計（無次元・単位 `""`）。
    pub magnitude: ErrorStats,
    /// 合否（`greatest` が `maximum_seconds` 以内 **かつ** `magnitude` が許容以内）。γ は非ゲート。
    pub pass: bool,
}

/// 複数日食の [`GlobalErrors`] を metric 別に集計し、[`ToleranceProfile`] で合否判定する。
///
/// 各 metric を [`ErrorStats::from_errors`]（絶対値化）で統計化する。合否は **最大食時刻**
/// （`maximum_seconds`）と **食分**（`magnitude`）のみゲート（γ は許容未設定＝統計のみ・非ゲート）。
/// 空入力は全 metric が空統計（n=0・全 0.0）で `pass = true`（vacuous, accuracy.md §3.4）。
pub fn aggregate_global(errors: &[GlobalErrors], profile: &ToleranceProfile) -> GlobalReport {
    let greatest = ErrorStats::from_errors(
        &errors
            .iter()
            .map(|e| e.greatest_seconds)
            .collect::<Vec<_>>(),
        "s",
    );
    let gamma = ErrorStats::from_errors(&errors.iter().map(|e| e.gamma).collect::<Vec<_>>(), "Re");
    let magnitude =
        ErrorStats::from_errors(&errors.iter().map(|e| e.magnitude).collect::<Vec<_>>(), "");
    let pass = greatest.within(profile.maximum_seconds) && magnitude.within(profile.magnitude);
    GlobalReport {
        greatest,
        gamma,
        magnitude,
        pass,
    }
}

#[cfg(test)]
mod tests {
    use super::percentile_r7_sorted;

    /// 非公開 `percentile_r7_sorted` の**上端境界**を直接縛る。`from_errors` は p=0.95 固定で
    /// `lo+1 < n` が恒真のため、防御分岐 `if lo+1 >= n { sorted_abs[n-1] }` は公開 API 経由では
    /// 到達不能（cargo-mutants で 88-89 行が生存）。一般 p（=1.0 で最大、0.0 で最小、0.5 で中央）を
    /// 直接呼んでガードと true 分岐を実効化する。
    /// 殺す変異: ガード `lo+1>=n` の `+`→`*`（到達時 else で OOB）、true 分岐 `sorted_abs[n-1]` の
    ///   `n-1` 改変（`n+1`/`n/1` で OOB）。
    #[test]
    fn percentile_r7_boundaries_direct() {
        let sorted = [1.0_f64, 2.0, 3.0];
        // p=1.0: h=(3-1)*1.0=2.0, lo=2=n-1 → ガード true 分岐 sorted_abs[n-1]=3.0。
        assert!(
            (percentile_r7_sorted(&sorted, 1.0) - 3.0).abs() < 1e-12,
            "p=1.0 は最大値 3.0"
        );
        // p=0.0: h=0, lo=0 → else 下端 sorted_abs[0]=1.0。
        assert!(
            (percentile_r7_sorted(&sorted, 0.0) - 1.0).abs() < 1e-12,
            "p=0.0 は最小値 1.0"
        );
        // p=0.5: h=(3-1)*0.5=1.0, lo=1, lo+1<n → 通常分岐 sorted_abs[1]=2.0。
        assert!(
            (percentile_r7_sorted(&sorted, 0.5) - 2.0).abs() < 1e-12,
            "p=0.5 は中央値 2.0"
        );
    }
}
