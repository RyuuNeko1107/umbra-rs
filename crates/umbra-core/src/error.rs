//! 基盤エラー型（`docs/api-draft.md` §1.6 / `docs/conventions.md` §11）。
//!
//! 下位の具体的な失敗を表す。上位（`umbra-eclipse`）の `EclipseError` はこれらを
//! `#[from]` で透過ラップする方針（`docs/issues/ISSUE-044`）。

use thiserror::Error;

/// 値域違反（入力の正規化前チェックなど）。
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq, Eq, Error)]
pub enum DomainError {
    /// `what` が許容範囲外。
    #[error("value out of range: {what}")]
    OutOfRange {
        /// 範囲外だった量の名前。
        what: &'static str,
    },
}

/// 数値解法（求根・最小化）の失敗。
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq, Eq, Error)]
pub enum SolverError {
    /// 区間端で符号が同じ＝根がブラケットされていない。
    #[error("root is not bracketed: f(a) and f(b) have the same sign")]
    RootNotBracketed,
    /// 反復上限内に収束しなかった。
    #[error("solver did not converge within the iteration limit")]
    DidNotConverge,
    /// 数値的不安定（NaN/Inf など）を検出。
    #[error("numerical instability encountered")]
    NumericalInstability,
}
