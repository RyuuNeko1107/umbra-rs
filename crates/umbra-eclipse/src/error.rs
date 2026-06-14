//! 日食エンジンのエラー（`docs/issues/ISSUE-044` の最小版。今後 From 変換等で拡張）。

use thiserror::Error;

/// 日食計算のエラー。
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq, Eq, Error)]
pub enum EclipseError {
    /// 影幾何が退化している（太陽・月が同一点など）。
    #[error("degenerate shadow geometry")]
    DegenerateGeometry,
}
