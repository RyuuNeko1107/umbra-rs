//! D6: 日食フィルタの偽陰性ゼロ・マージン実余裕統計（accuracy.md §3.4・ISSUE-018/030）。
//!
//! 日食候補フィルタ（[`crate::eclipse_filter`]）は `possible = (separation < bare_limit + margin)`
//! で朔を早期棄却する（`bare_limit = s_sun + s_moon + π_moon`、`margin = SAFETY_MARGIN_RAD`）。
//! 既知の実日食/grazing 部分食が全て採用されることは [`crate::eclipse_filter`] の単体テストが**サンプルで**
//! 担保するが、本モジュールは**全期間スキャンでマージンがどれだけ余裕を持って効いたか**（実消費・実余裕）を
//! 統計化する（D6・「採用マージンに対し実際にどれだけ余裕があったか」）。
//!
//! - **実消費** = 採用候補のうち `max(0, separation − bare_limit)` の最大（grazing でマージンに頼って
//!   拾われた量）。`bare_limit` 内（separation < bare_limit）の候補は消費 0。
//! - **実余裕** = `margin − 実消費`（最も棄却に近い採用候補がマージンに残した余裕）。`>0` ならマージンは
//!   枯渇しておらず、採用側で偽陰性は起きていない（採用候補 ⊇ 真の日食ゆえ保守的な下限余裕）。
//!
//! 注: 本統計は採用候補（偽陽性を含みうる）に対する保守的な下限余裕で、真の日食の余裕 ≥ 本値。真偽分類
//! （ベッセル/全球）は行わない＝候補→合→フィルタの軽量前段のみ（engine.search の前段と同一経路）。
//! coarse-scan のさらに内側（合 solver の刻み等）の余裕は対象外。

use crate::candidates::new_moon_candidates;
use crate::conjunction::{solve_conjunction, ConjunctionKind, RootConfig};
use crate::eclipse_filter::{assess_eclipse_possibility, SAFETY_MARGIN_RAD};
use crate::engine::UtcRange;
use crate::error::EclipseError;

/// 日食フィルタの偽陰性ゼロ保証マージン \[rad\]（[`crate::eclipse_filter`] の `SAFETY_MARGIN_RAD`・D6）。
pub const ECLIPSE_FILTER_SAFETY_MARGIN_RAD: f64 = SAFETY_MARGIN_RAD;

/// `scan_filter_margins` の合 solver 刻み（search 既定と同じ厳しさ・x_tolerance_days=1e-7≈8.6 ms）。
const SCAN_ROOT_CONFIG: RootConfig = RootConfig {
    x_tolerance_days: 1e-7,
    max_iterations: 100,
};

/// 1 候補のマージン標本（[`aggregate_filter_margins`] の入力）。
///
/// `separation_rad` は合の月-太陽角距離、`bare_limit_rad` は食限（マージン抜き）= s_sun+s_moon+π_moon、
/// `possible` はフィルタの採用判定（`separation < bare_limit + margin`）。
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MarginSample {
    /// 合の月-太陽角距離 \[rad\]。
    pub separation_rad: f64,
    /// 食限（マージン抜き）= s_sun + s_moon + π_moon \[rad\]。
    pub bare_limit_rad: f64,
    /// フィルタが採用（日食を起こしうる）と判定したか。
    pub possible: bool,
}

/// D6 マージン実余裕統計（偽陰性ゼロ・マージンの全期間での消費/余裕）。
#[derive(Clone, Debug, PartialEq)]
pub struct FilterMarginStats {
    /// 走査した候補（朔）の総数。
    pub candidates_scanned: usize,
    /// 採用（possible=true）された候補数。
    pub accepted: usize,
    /// 棄却（possible=false）された候補数。
    pub rejected: usize,
    /// 採用したマージン \[rad\]（`ECLIPSE_FILTER_SAFETY_MARGIN_RAD`）。
    pub safety_margin_rad: f64,
    /// 採用候補が実際に消費したマージン最大 = max over accepted of `max(0, separation − bare_limit)`
    /// （採用が無ければ 0.0）。grazing でマージンに頼って拾われた最大量。
    pub max_margin_consumed_rad: f64,
    /// 採用側の実余裕 = `safety_margin − max_margin_consumed`（最も棄却に近い採用候補の余裕・通常 0..=safety）。
    pub min_accepted_slack_rad: f64,
}

/// 標本群から D6 統計を集計する（純）。
///
/// `candidates_scanned`=標本数、`accepted`=possible 数、`rejected`=差。`max_margin_consumed_rad` は
/// **採用標本**の `max(0, separation − bare_limit)` の最大（採用無し/全て bare_limit 内なら 0.0）。
/// `min_accepted_slack_rad = safety_margin_rad − max_margin_consumed_rad`。棄却標本は消費に寄与しない。
pub fn aggregate_filter_margins(
    samples: &[MarginSample],
    safety_margin_rad: f64,
) -> FilterMarginStats {
    let candidates_scanned = samples.len();
    let accepted = samples.iter().filter(|s| s.possible).count();
    let rejected = candidates_scanned - accepted;
    // 採用候補の消費 max(0, sep−bare) の最大。0.0 始点ゆえ採用無し/全 bare_limit 内は 0.0（NaN 不出）。
    let max_margin_consumed_rad = samples
        .iter()
        .filter(|s| s.possible)
        .map(|s| (s.separation_rad - s.bare_limit_rad).max(0.0))
        .fold(0.0_f64, f64::max);
    let min_accepted_slack_rad = safety_margin_rad - max_margin_consumed_rad;
    FilterMarginStats {
        candidates_scanned,
        accepted,
        rejected,
        safety_margin_rad,
        max_margin_consumed_rad,
        min_accepted_slack_rad,
    }
}

/// 範囲を candidate→合→フィルタで実走し D6 統計を返す（実パイプライン・**やや低速**）。
///
/// `engine.search` の前段（[`new_moon_candidates`]→[`solve_conjunction`]→[`assess_eclipse_possibility`]）と
/// **同一経路**で各朔のマージン標本を集め、[`aggregate_filter_margins`] で集計する。Besselian/全球分類は
/// 行わない（軽量）。1972 年より前を含む範囲は [`new_moon_candidates`] が `Err(EclipseError::Time)`。
///
/// **エラー伝播（意図的）**: `solve_conjunction` の失敗（`RootNotBracketed` 等）は `?` で即伝播し統計
/// 全体を `Err` にする（`engine.search` と同一方針）。[`new_moon_candidates`] が候補窓内に合の存在を
/// 偽陰性ゼロで保証するため合の失敗は本来起こらず、起きた場合は**偽陰性ゼロ保証の破れを示す異常**で、
/// 黙ってスキップせず表面化させるべき（本ツール自体が偽陰性ゼロの検証＝失敗を隠すと本末転倒）。
pub fn scan_filter_margins(range: UtcRange) -> Result<FilterMarginStats, EclipseError> {
    let mut samples = Vec::new();
    for candidate in new_moon_candidates(range)? {
        let conjunction = solve_conjunction(
            &candidate,
            ConjunctionKind::EclipticLongitude,
            SCAN_ROOT_CONFIG,
        )?;
        let assessment = assess_eclipse_possibility(&conjunction);
        samples.push(MarginSample {
            separation_rad: assessment.min_separation.0,
            bare_limit_rad: assessment.bare_limit.0,
            possible: assessment.possible,
        });
    }
    Ok(aggregate_filter_margins(
        &samples,
        ECLIPSE_FILTER_SAFETY_MARGIN_RAD,
    ))
}
