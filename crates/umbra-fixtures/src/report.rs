//! 誤差統計と許容プロファイル（`docs/issues/ISSUE-030` S30a）。
//!
//! 検証レポータの**純粋な基盤**: 誤差列の記述統計 [`ErrorStats`] と、合否判定に使う
//! 許容値 [`ToleranceProfile`]（`docs/accuracy.md` §2）を提供する。ゴールデン比較・層分解・
//! JPL DE 差分・1900〜2100 全走査・JSON/CLI 出力は後続スライス（本モジュールは比較の前段＝
//! 統計と許容の定義に徹する）。
//!
//! 設計規律（conventions §11 / accuracy.md §4）: 統計は**誤差を隠さない**。pass 判定が通っても
//! 数値（max/mean/p95）は必ず保持する。許容を pass のために拡大しない。

use umbra_core::{TtInstant, UtcInstant};
use umbra_eclipse::{EclipseError, LocalCircumstances, SolarEclipse, SolarEclipseKind};

use crate::types::{GoldenEclipse, GoldenLocation};

/// 1 日の秒数（JD 差 → 秒の換算）。
const SECONDS_PER_DAY: f64 = 86_400.0;

/// 接触 1 点の時刻誤差（秒, computed − golden）。golden が TT(=TD) を持てば **TT 差**
/// （純幾何・ΔT 非依存, accuracy.md §0(a)）、無ければ UTC 差で代替する（§0(b)）。
/// 差分は [`umbra_core::JulianDate2::days_since`] で 2 要素を保ち、単一 f64 の jd() 同士の減算で
/// 生じる JD≈2.45e6 の桁落ち（~4.6e-5 s）を回避する（julian §2要素表現の理由）。
fn contact_time_error_seconds(
    computed_utc: UtcInstant,
    computed_tt: TtInstant,
    golden_utc: UtcInstant,
    golden_tt: Option<TtInstant>,
) -> f64 {
    match golden_tt {
        Some(g_tt) => computed_tt.jd2().days_since(g_tt.jd2()) * SECONDS_PER_DAY,
        None => computed_utc.jd2().days_since(golden_utc.jd2()) * SECONDS_PER_DAY,
    }
}

/// 1 項目の誤差記述統計（**絶対誤差**ベース, accuracy.md §3.4）。
///
/// 各 metric（接触秒・最大食秒・食分・食面積・高度…）の誤差列から、最大絶対誤差・平均絶対誤差・
/// 95 パーセンタイルを保持する。`units` は表示・取り違え防止のための単位識別子（例 `"s"`/`"deg"`）。
///
/// percentile は **R-7（線形補間, NumPy 既定）** に固定する（ISSUE-030 §65 の「規約を 1 つに固定」）:
/// 昇順絶対誤差 `a` と分位 `p` に対し `h = (n-1)·p`、`lo = ⌊h⌋`、
/// `a[lo] + (h-lo)·(a[lo+1]-a[lo])`（端は `a[n-1]`）。`n==1` は単一値、`n==0` は 0.0。
#[derive(Clone, Debug, PartialEq, serde::Serialize)]
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
    let greatest_seconds = contact_time_error_seconds(
        greatest.time_utc,
        greatest.time_tt,
        golden.greatest_time_utc,
        golden.greatest_time_tt,
    );
    GlobalErrors {
        greatest_seconds,
        gamma: computed.global.gamma - golden.gamma,
        magnitude: greatest.magnitude.0 - golden_magnitude_engine_convention(golden),
    }
}

/// golden（NASA 5MCSE）の全球食分をエンジン規約（太陽直径の隠蔽率）へ換算する。
///
/// NASA「Eclipse Magnitude」は**月/太陽の見かけ直径比**（glossary: "strictly a ratio of diameters"）。
/// 一方エンジンの `greatest.magnitude` は標準の隠蔽率 `(l1'−m)/(l1'+l2')` で、**中心食の最大食点
/// （m=0）では `l1'/(l1'+l2') = (1+直径比)/2`**（l1'∝rs+rm, l2'∝rs−rm）。よって中心食
/// （Total/Annular/Hybrid）では直径比を `(1+ratio)/2` へ換算して apples-to-apples 比較する。
/// 部分食（および非中心・将来種別）はエンジンの最大食点が縁端（m>0）で NASA も隠蔽率を報告する
/// ため換算しない（同義）。残差は ΔT/暦由来の真の精度誤差のみとなる（accuracy.md §3.1）。
fn golden_magnitude_engine_convention(golden: &GoldenEclipse) -> f64 {
    match golden.kind_expected {
        SolarEclipseKind::Total | SolarEclipseKind::Annular | SolarEclipseKind::Hybrid => {
            (1.0 + golden.magnitude) / 2.0
        }
        _ => golden.magnitude,
    }
}

/// 全球比較の metric 別統計＋合否（accuracy.md §3.4: pass でも統計を必ず出す）。
#[derive(Clone, Debug, PartialEq, serde::Serialize)]
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

/// 1 地点の局地条件の誤差（**符号付き = computed − golden**, accuracy.md §3.4）。
///
/// 時刻は秒（TT 基準・golden が TT を持てば TT 差）。接触 C1〜C4 は両方存在するもののみ
/// `contact_seconds`（時系列順）へ集め、Some/None が食い違う接触は `contact_presence_mismatches`
/// で数える（時刻誤差には混ぜない）。
#[derive(Clone, Debug, PartialEq)]
pub struct LocalErrors {
    /// 最大食接触の時刻誤差 \[s\]（computed − golden）。
    pub maximum_seconds: f64,
    /// c1,c2,c3,c4 のうち**両方 Some**の接触の時刻誤差 \[s\]（時系列順）。
    pub contact_seconds: Vec<f64>,
    /// c1..c4 で computed/golden の Some/None が食い違った数（0..=4）。
    pub contact_presence_mismatches: usize,
    /// 食分誤差（無次元, computed − golden）。
    pub magnitude: f64,
    /// 食面積誤差（無次元, computed − golden）。
    pub obscuration: f64,
    /// 最大食での太陽高度誤差 \[deg\]（computed − golden）。
    pub max_altitude_deg: f64,
    /// 可視性が一致するか（`computed.visibility == golden.visibility_expected`）。
    pub visibility_matches: bool,
}

/// computed の局地条件を golden 地点と比較し、符号付き誤差（computed − golden）を返す（純粋）。
///
/// 各接触時刻は [`contact_time_error_seconds`]（TT 優先・days_since・2 要素保持）。最大食は常に存在。
/// C1〜C4 は両方 Some のみ時刻誤差化し、Some/None 食い違いは `contact_presence_mismatches` に計上。
pub fn compare_local(computed: &LocalCircumstances, golden: &GoldenLocation) -> LocalErrors {
    let contacts = &computed.contacts;
    let maximum_seconds = contact_time_error_seconds(
        contacts.maximum.time_utc,
        contacts.maximum.time_tt,
        golden.maximum.time_utc,
        golden.maximum.time_tt,
    );
    let mut contact_seconds = Vec::new();
    let mut contact_presence_mismatches = 0usize;
    // 時系列順 c1,c2,c3,c4 で computed(Option<LocalContact>) と golden(Option<GoldenContact>) を対にする。
    let pairs = [
        (contacts.c1.as_ref(), golden.c1.as_ref()),
        (contacts.c2.as_ref(), golden.c2.as_ref()),
        (contacts.c3.as_ref(), golden.c3.as_ref()),
        (contacts.c4.as_ref(), golden.c4.as_ref()),
    ];
    for (computed_contact, golden_contact) in pairs {
        match (computed_contact, golden_contact) {
            (Some(local), Some(gold)) => contact_seconds.push(contact_time_error_seconds(
                local.time_utc,
                local.time_tt,
                gold.time_utc,
                gold.time_tt,
            )),
            // 地平下の computed 接触（visible=false）は USNO（golden）が省略する慣習と一致＝
            // 不一致に数えない（両表現とも「観測不能」で合致。日没/日の出食の C1/C4 等）。
            // 地平上の接触を golden が持たない場合は真の不一致。
            (Some(local), None) => {
                if local.visible {
                    contact_presence_mismatches += 1;
                }
            }
            (None, Some(_)) => contact_presence_mismatches += 1,
            (None, None) => {}
        }
    }
    LocalErrors {
        maximum_seconds,
        contact_seconds,
        contact_presence_mismatches,
        magnitude: computed.magnitude.0 - golden.magnitude,
        obscuration: computed.obscuration.0 - golden.obscuration,
        max_altitude_deg: computed.maximum_altitude.0 - golden.max_altitude_deg,
        visibility_matches: computed.visibility == golden.visibility_expected,
    }
}

/// 地点別比較の metric 別統計＋合否（accuracy.md §3.4: pass でも統計を必ず出す）。
#[derive(Clone, Debug, PartialEq, serde::Serialize)]
pub struct LocalReport {
    /// 最大食接触時刻誤差の統計（単位 `"s"`）。
    pub maximum: ErrorStats,
    /// C1〜C4 接触時刻誤差の統計（全地点・全接触をフラット集計, 単位 `"s"`）。
    pub contacts: ErrorStats,
    /// 食分誤差の統計（無次元・単位 `""`）。
    pub magnitude: ErrorStats,
    /// 食面積誤差の統計（無次元・単位 `""`）。
    pub obscuration: ErrorStats,
    /// 最大食高度誤差の統計（単位 `"deg"`）。
    pub max_altitude: ErrorStats,
    /// 可視性が不一致だった地点数（全地点合計）。
    pub visibility_mismatches: usize,
    /// 接触の Some/None 食い違い数（全地点合計）。
    pub contact_presence_mismatches: usize,
    /// 合否（全 metric が許容以内、かつ可視性不一致・接触食い違いがいずれも 0）。
    pub pass: bool,
}

/// 複数地点の [`LocalErrors`] を metric 別に集計し、[`ToleranceProfile`] で合否判定する。
///
/// 接触は全地点・全接触をフラット化して 1 つの [`ErrorStats`] にする。合否は最大食/接触/食分/食面積/
/// 高度の許容内、かつ可視性不一致 0・接触食い違い 0。空入力は全空統計・mismatch 0・`pass = true`（vacuous）。
pub fn aggregate_local(errors: &[LocalErrors], profile: &ToleranceProfile) -> LocalReport {
    let maximum = ErrorStats::from_errors(
        &errors.iter().map(|e| e.maximum_seconds).collect::<Vec<_>>(),
        "s",
    );
    let contacts = ErrorStats::from_errors(
        &errors
            .iter()
            .flat_map(|e| e.contact_seconds.iter().copied())
            .collect::<Vec<_>>(),
        "s",
    );
    let magnitude =
        ErrorStats::from_errors(&errors.iter().map(|e| e.magnitude).collect::<Vec<_>>(), "");
    let obscuration = ErrorStats::from_errors(
        &errors.iter().map(|e| e.obscuration).collect::<Vec<_>>(),
        "",
    );
    let max_altitude = ErrorStats::from_errors(
        &errors
            .iter()
            .map(|e| e.max_altitude_deg)
            .collect::<Vec<_>>(),
        "deg",
    );
    let visibility_mismatches = errors.iter().filter(|e| !e.visibility_matches).count();
    let contact_presence_mismatches = errors.iter().map(|e| e.contact_presence_mismatches).sum();
    let pass = maximum.within(profile.maximum_seconds)
        && contacts.within(profile.contact_seconds)
        && magnitude.within(profile.magnitude)
        && obscuration.within(profile.obscuration)
        && max_altitude.within(profile.altitude_degrees)
        && visibility_mismatches == 0
        && contact_presence_mismatches == 0;
    LocalReport {
        maximum,
        contacts,
        magnitude,
        obscuration,
        max_altitude,
        visibility_mismatches,
        contact_presence_mismatches,
        pass,
    }
}

/// ゴールデンから「当日の食」と「指定地点の局地条件」を計算する注入インターフェース。
///
/// 実装は実エンジン（`EclipseEngine` 経由・SLOW）でもモック（テスト）でもよい。これにより
/// [`report_against_golden`] のオーケストレーション（ループ・取りこぼし計数・集計）を、実エンジンを
/// 走らせずに検証できる（エンジン結線は SLOW 統合テストで担保, ISSUE-030 S30d）。
pub trait GoldenComputer {
    /// `golden` の食日に起こる食を返す（見つからなければ `None`＝取りこぼし）。
    fn eclipse_on(&self, golden: &GoldenEclipse) -> Result<Option<SolarEclipse>, EclipseError>;
    /// 指定食・指定地点の局地条件を返す。
    fn local_at(
        &self,
        eclipse: &SolarEclipse,
        location: &GoldenLocation,
    ) -> Result<LocalCircumstances, EclipseError>;
}

/// ゴールデン照合の総合レポート（全球＋地点別の集計と被覆カウント）。
#[derive(Clone, Debug, PartialEq, serde::Serialize)]
pub struct GoldenReport {
    /// 全球条件の集計レポート。
    pub global: GlobalReport,
    /// 地点別条件の集計レポート。
    pub local: LocalReport,
    /// computer が食を返した（=照合できた）golden 数。
    pub eclipses_found: usize,
    /// computer が `None` を返した golden 数（取りこぼし＝異常シグナル）。
    pub eclipses_missing: usize,
    /// 比較した地点の総数（照合できた食の地点合計）。
    pub locations_compared: usize,
}

/// golden 群に `computer` を適用し、全球＋地点別比較を集計したレポートを返す（純粋オーケストレーション）。
///
/// 各 golden について `eclipse_on` で食を取得（`None` は `eclipses_missing` に計上し以降スキップ）、
/// 見つかれば [`compare_global`] と、各地点の `local_at` → [`compare_local`] を収集する。最後に
/// [`aggregate_global`]/[`aggregate_local`] で集計する。`eclipse_on`/`local_at` のエラーは伝播する。
pub fn report_against_golden<C: GoldenComputer>(
    computer: &C,
    golden: &[GoldenEclipse],
    profile: &ToleranceProfile,
) -> Result<GoldenReport, EclipseError> {
    let mut global_errors = Vec::new();
    let mut local_errors = Vec::new();
    let mut eclipses_missing = 0usize;
    for g in golden {
        match computer.eclipse_on(g)? {
            Some(eclipse) => {
                global_errors.push(compare_global(&eclipse, g));
                for location in &g.locations {
                    let local = computer.local_at(&eclipse, location)?;
                    local_errors.push(compare_local(&local, location));
                }
            }
            None => eclipses_missing += 1,
        }
    }
    let eclipses_found = global_errors.len();
    let locations_compared = local_errors.len();
    Ok(GoldenReport {
        global: aggregate_global(&global_errors, profile),
        local: aggregate_local(&local_errors, profile),
        eclipses_found,
        eclipses_missing,
        locations_compared,
    })
}

/// 1 metric の [`ErrorStats`] を 1 行に整形する（`render_text` 用）。`n / max / mean / p95`＋単位。
fn stats_line(label: &str, stats: &ErrorStats) -> String {
    format!(
        "  {label}: n={n} max={max:.4}{u} mean={mean:.4}{u} p95={p95:.4}{u}\n",
        n = stats.n,
        max = stats.max_abs,
        mean = stats.mean_abs,
        p95 = stats.p95,
        u = stats.units,
    )
}

/// [`GoldenReport`] を人間可読サマリ（複数行テキスト）に整形する（ISSUE-030 §82・誤差を隠さない）。
///
/// 被覆カウント（発見/取りこぼし/比較地点）、全球・地点別の各 metric 統計（n/max/mean/p95＋単位）、
/// 可視性/接触 presence 不一致、global/local の合否を出す。pass でも統計を必ず表示する。
pub fn render_text(report: &GoldenReport) -> String {
    let mut out = String::new();
    out.push_str("=== Golden comparison report ===\n");
    out.push_str(&format!(
        "eclipses: {} found, {} missing | locations compared: {}\n",
        report.eclipses_found, report.eclipses_missing, report.locations_compared,
    ));
    out.push_str(&format!("GLOBAL  pass: {}\n", report.global.pass));
    out.push_str(&stats_line("greatest", &report.global.greatest));
    out.push_str(&stats_line("gamma", &report.global.gamma));
    out.push_str(&stats_line("magnitude", &report.global.magnitude));
    out.push_str(&format!("LOCAL   pass: {}\n", report.local.pass));
    out.push_str(&stats_line("maximum", &report.local.maximum));
    out.push_str(&stats_line("contacts", &report.local.contacts));
    out.push_str(&stats_line("magnitude", &report.local.magnitude));
    out.push_str(&stats_line("obscuration", &report.local.obscuration));
    out.push_str(&stats_line("max_altitude", &report.local.max_altitude));
    out.push_str(&format!(
        "  visibility mismatches: {}  contact presence mismatches: {}\n",
        report.local.visibility_mismatches, report.local.contact_presence_mismatches,
    ));
    out
}

/// [`GoldenReport`] を機械可読 JSON（pretty・末尾改行）に整形する（ISSUE-030 §82・CI/履歴比較用）。
pub fn render_json(report: &GoldenReport) -> Result<String, serde_json::Error> {
    let mut out = serde_json::to_string_pretty(report)?;
    out.push('\n');
    Ok(out)
}

// ============================================================
// DE 差分・誤差層分解（accuracy.md §4 — 同一パイプライン 2 エンジン差分法）
// ============================================================

/// 1 metric の DE 差分・誤差層分解（accuracy.md §4 / §3.1-1）。
///
/// 同一のベッセル/接触パイプラインに**解析暦（analytical）**と **JPL DE（de）**を通し、metric の
/// 誤差を暦層と幾何/数値層へ帰属する。各層は computed−reference の符号付き誤差列を [`ErrorStats`]
/// で絶対値統計化する:
/// - `ephemeris`: analytical − DE（同一パイプライン → 差は**暦のみ**）。
/// - `geometry`: DE − golden オラクル（DE 入力 → 残差は**幾何/数値＋慣習差**）。
/// - `total`: analytical − golden（実測総誤差）。各サンプルで **符号付き total ≡ ephemeris + geometry**
///   （同じ 3 実数 a,d,g から (a−d)+(d−g)=a−g）。
///
/// **設計（ISSUE-030 §38-45 スケッチからの確定逸脱・accuracy.md §4 に記録）**: issue が描く 6 物理層
/// （time/sun/moon/shadow/poly/contact の個別）は**エンジン内部の計装が必要で 030 のスコープ外**。
/// 同一パイプライン 2 エンジン差分で**清浄に切り出せる粒度は暦層 vs 幾何/数値層の 2 バケット**に限る
/// ため、ここへ集約する。誤差を隠さない原則（conventions §11）は「内部層を個別測定せず統合表示する」
/// と明記して担保する。sun/月の**生の方向残差**（0.040″/0.268″）は `umbra-ephemeris` の DE 差分
/// （`tests/de_diff.rs`・accuracy.md §3.3）が暦層で別途担保しており、本レポートはそれを**日食 metric
/// への寄与**として再表現する位置づけ。
#[derive(Clone, Debug, PartialEq, serde::Serialize)]
pub struct LayeredError {
    /// 暦層: analytical − DE（同一パイプライン → 暦差のみ）。
    pub ephemeris: ErrorStats,
    /// 幾何/数値層: DE − golden オラクル（DE 入力 → 幾何/数値＋慣習差）。
    pub geometry: ErrorStats,
    /// 総合: analytical − golden（実測総誤差・恒等性 ephemeris⊕geometry の確認）。
    pub total: ErrorStats,
}

/// 1 metric の 3 層誤差列を貯める累積器（符号付き = computed − reference）。
#[derive(Default)]
struct LayerAccumulator {
    eph: Vec<f64>,
    geo: Vec<f64>,
    tot: Vec<f64>,
}

impl LayerAccumulator {
    /// スカラ metric（食分・食面積・高度）: a,d,g の 3 値から 3 層を導出して push。
    fn push_scalar(&mut self, a: f64, d: f64, g: f64) {
        self.eph.push(a - d);
        self.geo.push(d - g);
        self.tot.push(a - g);
    }

    /// 時刻 metric: 各層の差（秒）を**直接** push する（2 要素 days_since で桁落ち回避済みの値を渡す）。
    fn push_diffs(&mut self, eph: f64, geo: f64, tot: f64) {
        self.eph.push(eph);
        self.geo.push(geo);
        self.tot.push(tot);
    }

    /// 3 層を [`ErrorStats`]（絶対値統計）へ確定する。空列でも `units` を保持する。
    fn finish(&self, units: &'static str) -> LayeredError {
        LayeredError {
            ephemeris: ErrorStats::from_errors(&self.eph, units),
            geometry: ErrorStats::from_errors(&self.geo, units),
            total: ErrorStats::from_errors(&self.tot, units),
        }
    }
}

/// computed−computed の時刻差（秒）。golden が TT を持てば **TT 差**、無ければ **UTC 差**で測る。
///
/// golden 比較（[`contact_time_error_seconds`]）と**同じ時刻表現**を使うことで、各サンプルの
/// 符号付き恒等性 `(a−d)+(d−g)=a−g` を厳密に保つ（TT/UTC が層間で食い違うと恒等性が崩れる）。
/// 2 要素 [`umbra_core::JulianDate2::days_since`] で JD≈2.45e6 の桁落ちを回避する。
fn computed_pair_time_error_seconds(
    a_utc: UtcInstant,
    a_tt: TtInstant,
    d_utc: UtcInstant,
    d_tt: TtInstant,
    golden_has_tt: bool,
) -> f64 {
    if golden_has_tt {
        a_tt.jd2().days_since(d_tt.jd2()) * SECONDS_PER_DAY
    } else {
        a_utc.jd2().days_since(d_utc.jd2()) * SECONDS_PER_DAY
    }
}

/// DE 差分・層分解レポート（全球＋地点別 metric 毎の [`LayeredError`] ＋被覆カウント）。
///
/// 各 metric は ephemeris（暦層）/ geometry（幾何/数値層）/ total（総合）の 3 層に分解される
/// （[`LayeredError`]）。合否（pass/fail）は持たない**純粋な誤差層分解**（誤差を隠さない・原因層の
/// 特定が目的）。合否判定は [`report_against_golden`] / [`ToleranceProfile`] の領分。
#[derive(Clone, Debug, PartialEq, serde::Serialize)]
pub struct DifferentialReport {
    /// 全球最大食時刻の層分解（単位 `"s"`）。
    pub global_greatest_seconds: LayeredError,
    /// 全球食分の層分解（無次元・単位 `""`）。
    pub global_magnitude: LayeredError,
    /// 地点最大食接触時刻の層分解（単位 `"s"`）。
    pub local_maximum_seconds: LayeredError,
    /// 地点 C1〜C4 接触時刻の層分解（3 者 Some の接触のみ・単位 `"s"`）。
    pub local_contact_seconds: LayeredError,
    /// 地点食分の層分解（無次元・単位 `""`）。
    pub local_magnitude: LayeredError,
    /// 地点食面積の層分解（無次元・単位 `""`）。
    pub local_obscuration: LayeredError,
    /// 地点最大食高度の層分解（単位 `"deg"`）。
    pub local_max_altitude_deg: LayeredError,
    /// analytical も DE も食を返した（=層分解できた）golden 数。
    pub eclipses_compared: usize,
    /// analytical / DE のどちらか（または両方）が `None` を返した golden 数。
    pub eclipses_missing: usize,
    /// 比較した地点の総数（層分解できた食の地点合計）。
    pub locations_compared: usize,
}

/// 解析暦 computer と DE computer を同一 golden へ適用し、metric 毎に誤差を 3 層へ分解する（純）。
///
/// 各 golden について `analytical.eclipse_on` と `de.eclipse_on` を**両方**呼ぶ（どちらの `Err` も
/// 即伝播）。**両方 `Some`** の時のみ層分解対象とし、地点ごとに `local_at` を解いて metric を集める。
/// どちらかが `None` なら `eclipses_missing` に計上し、その食の地点は評価しない。各 metric は
/// ephemeris=analytical−DE / geometry=DE−golden / total=analytical−golden（[`LayerAccumulator`]）。
/// 接触秒は c1..c4 のうち **analytical・DE・golden の 3 者すべてが `Some`** の接触のみ寄与（時系列順）。
/// 空 golden は全層空統計・カウント 0・`Ok`（vacuous）。
pub fn report_differential<A: GoldenComputer, D: GoldenComputer>(
    analytical: &A,
    de: &D,
    golden: &[GoldenEclipse],
) -> Result<DifferentialReport, EclipseError> {
    let mut global_greatest = LayerAccumulator::default();
    let mut global_magnitude = LayerAccumulator::default();
    let mut local_maximum = LayerAccumulator::default();
    let mut local_contacts = LayerAccumulator::default();
    let mut local_magnitude = LayerAccumulator::default();
    let mut local_obscuration = LayerAccumulator::default();
    let mut local_altitude = LayerAccumulator::default();
    let mut eclipses_compared = 0usize;
    let mut eclipses_missing = 0usize;
    let mut locations_compared = 0usize;

    for g in golden {
        // 両エンジンを必ず呼ぶ（どちらの Err も伝播）。両方 Some の時のみ層分解する。
        let a_eclipse = analytical.eclipse_on(g)?;
        let d_eclipse = de.eclipse_on(g)?;
        let (Some(a_eclipse), Some(d_eclipse)) = (a_eclipse, d_eclipse) else {
            eclipses_missing += 1;
            continue;
        };
        eclipses_compared += 1;

        // --- 全球 ---
        let total_global = compare_global(&a_eclipse, g);
        let geo_global = compare_global(&d_eclipse, g);
        let a_greatest = &a_eclipse.global.greatest;
        let d_greatest = &d_eclipse.global.greatest;
        global_greatest.push_diffs(
            computed_pair_time_error_seconds(
                a_greatest.time_utc,
                a_greatest.time_tt,
                d_greatest.time_utc,
                d_greatest.time_tt,
                g.greatest_time_tt.is_some(),
            ),
            geo_global.greatest_seconds,
            total_global.greatest_seconds,
        );
        // 食分は ephemeris = a−d（両エンジン出力・換算不要）、geometry/total は compare_global
        // （golden をエンジン規約へ換算済み）を流用＝恒等性が保たれる。
        global_magnitude.push_diffs(
            a_greatest.magnitude.0 - d_greatest.magnitude.0,
            geo_global.magnitude,
            total_global.magnitude,
        );

        // --- 地点別 ---
        for location in &g.locations {
            locations_compared += 1;
            let a_local = analytical.local_at(&a_eclipse, location)?;
            let d_local = de.local_at(&d_eclipse, location)?;
            let total_local = compare_local(&a_local, location);
            let geo_local = compare_local(&d_local, location);

            // 最大食接触時刻。
            local_maximum.push_diffs(
                computed_pair_time_error_seconds(
                    a_local.contacts.maximum.time_utc,
                    a_local.contacts.maximum.time_tt,
                    d_local.contacts.maximum.time_utc,
                    d_local.contacts.maximum.time_tt,
                    location.maximum.time_tt.is_some(),
                ),
                geo_local.maximum_seconds,
                total_local.maximum_seconds,
            );
            // スカラ metric（食分/食面積/高度）。golden は compare_local と同じ raw 値で素通し。
            local_magnitude.push_scalar(
                a_local.magnitude.0,
                d_local.magnitude.0,
                location.magnitude,
            );
            local_obscuration.push_scalar(
                a_local.obscuration.0,
                d_local.obscuration.0,
                location.obscuration,
            );
            local_altitude.push_scalar(
                a_local.maximum_altitude.0,
                d_local.maximum_altitude.0,
                location.max_altitude_deg,
            );

            // C1〜C4 接触: analytical・DE・golden の 3 者すべて Some の接触のみ層分解に寄与（時系列順）。
            let triples = [
                (
                    a_local.contacts.c1.as_ref(),
                    d_local.contacts.c1.as_ref(),
                    location.c1.as_ref(),
                ),
                (
                    a_local.contacts.c2.as_ref(),
                    d_local.contacts.c2.as_ref(),
                    location.c2.as_ref(),
                ),
                (
                    a_local.contacts.c3.as_ref(),
                    d_local.contacts.c3.as_ref(),
                    location.c3.as_ref(),
                ),
                (
                    a_local.contacts.c4.as_ref(),
                    d_local.contacts.c4.as_ref(),
                    location.c4.as_ref(),
                ),
            ];
            for (a_contact, d_contact, g_contact) in triples {
                let (Some(a_c), Some(d_c), Some(g_c)) = (a_contact, d_contact, g_contact) else {
                    continue;
                };
                local_contacts.push_diffs(
                    computed_pair_time_error_seconds(
                        a_c.time_utc,
                        a_c.time_tt,
                        d_c.time_utc,
                        d_c.time_tt,
                        g_c.time_tt.is_some(),
                    ),
                    contact_time_error_seconds(
                        d_c.time_utc,
                        d_c.time_tt,
                        g_c.time_utc,
                        g_c.time_tt,
                    ),
                    contact_time_error_seconds(
                        a_c.time_utc,
                        a_c.time_tt,
                        g_c.time_utc,
                        g_c.time_tt,
                    ),
                );
            }
        }
    }

    Ok(DifferentialReport {
        global_greatest_seconds: global_greatest.finish("s"),
        global_magnitude: global_magnitude.finish(""),
        local_maximum_seconds: local_maximum.finish("s"),
        local_contact_seconds: local_contacts.finish("s"),
        local_magnitude: local_magnitude.finish(""),
        local_obscuration: local_obscuration.finish(""),
        local_max_altitude_deg: local_altitude.finish("deg"),
        eclipses_compared,
        eclipses_missing,
        locations_compared,
    })
}

/// 1 metric の [`LayeredError`] を 3 層ブロック（ephemeris/geometry/total）に整形する（`render_text` 用）。
fn layered_block(label: &str, layered: &LayeredError) -> String {
    let mut out = format!("{label}:\n");
    out.push_str(&stats_line("ephemeris", &layered.ephemeris));
    out.push_str(&stats_line("geometry", &layered.geometry));
    out.push_str(&stats_line("total", &layered.total));
    out
}

/// [`DifferentialReport`] を人間可読サマリ（複数行テキスト）に整形する（誤差を隠さない: 全層を表示）。
///
/// 被覆カウント（層分解できた食/取りこぼし/比較地点）と、全 metric の 3 層（ephemeris/geometry/total）
/// の n/max/mean/p95＋単位を出す。pass/fail は持たない（純粋な誤差層分解）。
pub fn render_differential_text(report: &DifferentialReport) -> String {
    let mut out = String::new();
    out.push_str("=== DE differential layered report ===\n");
    out.push_str(&format!(
        "eclipses: {} compared, {} missing | locations compared: {}\n",
        report.eclipses_compared, report.eclipses_missing, report.locations_compared,
    ));
    out.push_str(&layered_block(
        "global_greatest_seconds",
        &report.global_greatest_seconds,
    ));
    out.push_str(&layered_block("global_magnitude", &report.global_magnitude));
    out.push_str(&layered_block(
        "local_maximum_seconds",
        &report.local_maximum_seconds,
    ));
    out.push_str(&layered_block(
        "local_contact_seconds",
        &report.local_contact_seconds,
    ));
    out.push_str(&layered_block("local_magnitude", &report.local_magnitude));
    out.push_str(&layered_block(
        "local_obscuration",
        &report.local_obscuration,
    ));
    out.push_str(&layered_block(
        "local_max_altitude_deg",
        &report.local_max_altitude_deg,
    ));
    out
}

/// [`DifferentialReport`] を機械可読 JSON（pretty・末尾改行）に整形する（CI/履歴比較用）。
pub fn render_differential_json(report: &DifferentialReport) -> Result<String, serde_json::Error> {
    let mut out = serde_json::to_string_pretty(report)?;
    out.push('\n');
    Ok(out)
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
