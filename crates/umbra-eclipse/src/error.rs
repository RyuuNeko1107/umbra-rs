//! 日食エンジンのエラー（`docs/issues/ISSUE-044` の最小版＋ISSUE-022 fit 残差超過）。
//!
//! `BesselFitExceededTolerance` は実測残差 `achieved` を保持する（誤差を隠さない, conventions §11）。
//! `BesselFitError`（f64 を含む）を載せるため列挙体は `Eq` を持たない（`PartialEq` のみ）。

use thiserror::Error;
use umbra_core::TimeError;

use crate::bessel_poly::BesselFitError;

/// 日食計算のエラー。
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq, Error)]
pub enum EclipseError {
    /// 影幾何が退化している（太陽・月が同一点など）。
    #[error("degenerate shadow geometry")]
    DegenerateGeometry,

    /// 時刻系変換に失敗した（例: 1972 年より前で UTC↔TT が定義できない, ISSUE-016）。
    #[error("time scale conversion failed: {0}")]
    Time(#[from] TimeError),

    /// ベッセル多項式 fit の残差が許容を超えた（ISSUE-022）。実測残差と要求許容を保持する。
    #[error("Besselian polynomial fit residual exceeded tolerance (achieved {achieved:?}, tolerance {tolerance:?})")]
    BesselFitExceededTolerance {
        /// 最良次数で実測した残差（誤差を隠さない）。
        achieved: BesselFitError,
        /// 要求された許容。
        tolerance: BesselFitError,
    },

    /// fit 区間が不正（経過時間で start ≥ end、または非有限, ISSUE-022）。
    #[error("invalid fit interval")]
    InvalidFitInterval,

    /// 多項式を fit 区間外で評価しようとした（多項式は区間内のみ妥当, ISSUE-022/037）。
    #[error("evaluation outside fit interval")]
    EvaluationOutsideFitInterval,
}
