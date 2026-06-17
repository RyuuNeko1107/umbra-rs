//! 日食エンジンの集約エラー（`docs/issues/ISSUE-044`・確定A6 / api-draft §3.5）。
//!
//! 各層のエラー（時刻 `TimeError` / 暦 `EphemerisError` / 求根 `SolverError` / 定義域 `DomainError`）を
//! `#[from]` + `#[error("{0}")]` で**透過ラップ**し、`?` で [`EclipseError`] に集約する
//! （同義の二重 variant は作らない・確定ERR）。直 variant は**日食固有**（下位層に無い失敗）のみ。
//! - 透過ラップ: [`EclipseError::Time`] / [`EclipseError::Ephemeris`] / [`EclipseError::Solver`] /
//!   [`EclipseError::Domain`]。`Display` は内側のメッセージへ委譲し、`#[from]` 由来の `#[source]` に
//!   より**内側エラーを `source()` として公開**する（原因連鎖を辿れる・§受け入れテスト3）。
//!   完全な `#[error(transparent)]` は `Display` と `source()` の**両方**を内側へ委譲し、leaf な基盤
//!   エラー（`TimeError` 等は `source()`＝None）では原因連鎖が途切れるため、`Display` のみ委譲する
//!   `#[error("{0}")]` を採用する（内側を 1 段の cause として露出。確定ERR の「透過」意図＝メッセージ
//!   非二重化を満たしつつ source 連鎖も提供）。
//! - 日食固有: [`EclipseError::DegenerateGeometry`] / [`EclipseError::BesselFitExceededTolerance`] /
//!   [`EclipseError::InvalidFitInterval`] / [`EclipseError::EvaluationOutsideFitInterval`] /
//!   [`EclipseError::NotImplemented`]（PATH: 未実装機能は panic でなくこれを返す, ISSUE-045）。
//!
//! 求根のブラケット失敗は下位 `SolverError::RootNotBracketed` を `Solver` で透過ラップして表すため、
//! 同義の直 `RootNotBracketed` variant は設けない（確定ERR §要確認の解決・既存呼出と整合）。
//! `BesselFitExceededTolerance` は実測残差 `achieved` を保持する（誤差を隠さない, conventions §11）。
//! `BesselFitError`（f64 を含む）を載せるため列挙体は `Eq` を持たない（`PartialEq` のみ）。

use thiserror::Error;
use umbra_core::{DomainError, SolverError, TimeError};
use umbra_ephemeris::EphemerisError;

use crate::bessel_poly::BesselFitError;

/// 日食計算の集約エラー（前方互換のため `#[non_exhaustive]`, api-draft §0）。
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq, Error)]
pub enum EclipseError {
    /// 時刻系変換の失敗（閏秒/EOP/不正暦, ISSUE-006/007）。透過ラップ。
    #[error("{0}")]
    Time(#[from] TimeError),

    /// 天体暦バックエンドの失敗（対応範囲外/データ欠落, ISSUE-012）。透過ラップ。
    #[error("{0}")]
    Ephemeris(#[from] EphemerisError),

    /// 数値求根/最小化の失敗（ブラケット不成立・収束せず・不安定, ISSUE-008/017）。透過ラップ。
    #[error("{0}")]
    Solver(#[from] SolverError),

    /// 定義域違反（観測者・範囲などの入力不正, ISSUE-001）。透過ラップ。
    #[error("{0}")]
    Domain(#[from] DomainError),

    /// 影幾何が退化している（太陽・月が同一点など）。日食固有。
    #[error("degenerate shadow geometry")]
    DegenerateGeometry,

    /// ベッセル多項式 fit の残差が許容を超えた（ISSUE-022）。実測残差と要求許容を保持する。
    #[error("Besselian polynomial fit residual exceeded tolerance (achieved {achieved:?}, tolerance {tolerance:?})")]
    BesselFitExceededTolerance {
        /// 最良次数で実測した残差（誤差を隠さない）。
        achieved: BesselFitError,
        /// 要求された許容。
        tolerance: BesselFitError,
    },

    /// fit 区間が不正（経過時間で start ≥ end、または非有限, ISSUE-022）。日食固有。
    #[error("invalid fit interval")]
    InvalidFitInterval,

    /// 多項式を fit 区間外で評価しようとした（多項式は区間内のみ妥当, ISSUE-022/037）。日食固有。
    #[error("evaluation outside fit interval")]
    EvaluationOutsideFitInterval,

    /// 未実装機能（PATH: `path()` 等。panic でなくこれを返す, ISSUE-043/045）。
    /// 「対応年代外」専用語義の流用はしない（未実装は本 variant）。
    #[error("not implemented")]
    NotImplemented,
}

#[cfg(test)]
mod tests {
    use super::EclipseError;
    use crate::bessel_poly::BesselFitError;
    use std::error::Error;
    use umbra_core::{DomainError, SolverError, TimeError};
    use umbra_ephemeris::EphemerisError;

    // ============================================================
    // 受け入れテスト1: From/`?` による各層（4 層）の透過ラップ変換
    // ============================================================

    // 各層のエラーを返す関数を ? で集約し、期待するラップ variant に一致するか確認する。

    fn time_via_question() -> Result<(), EclipseError> {
        fn inner() -> Result<(), TimeError> {
            Err(TimeError::MissingLeapSecondData)
        }
        inner()?;
        Ok(())
    }

    fn ephemeris_via_question() -> Result<(), EclipseError> {
        fn inner() -> Result<(), EphemerisError> {
            Err(EphemerisError::OutOfSupportedRange)
        }
        inner()?;
        Ok(())
    }

    fn solver_via_question() -> Result<(), EclipseError> {
        fn inner() -> Result<(), SolverError> {
            Err(SolverError::RootNotBracketed)
        }
        inner()?;
        Ok(())
    }

    fn domain_via_question() -> Result<(), EclipseError> {
        fn inner() -> Result<(), DomainError> {
            Err(DomainError::OutOfRange { what: "lat" })
        }
        inner()?;
        Ok(())
    }

    /// `TimeError` が `?` で `EclipseError::Time(_)` に透過ラップされる。
    #[test]
    fn time_error_wraps_via_question_mark() {
        let r = time_via_question();
        assert!(
            matches!(r, Err(EclipseError::Time(TimeError::MissingLeapSecondData))),
            "expected Err(EclipseError::Time(MissingLeapSecondData)), got {r:?}"
        );
    }

    /// `EphemerisError` が `?` で `EclipseError::Ephemeris(_)` に透過ラップされる（NEW）。
    #[test]
    fn ephemeris_error_wraps_via_question_mark() {
        let r = ephemeris_via_question();
        assert!(
            matches!(
                r,
                Err(EclipseError::Ephemeris(EphemerisError::OutOfSupportedRange))
            ),
            "expected Err(EclipseError::Ephemeris(OutOfSupportedRange)), got {r:?}"
        );
    }

    /// `SolverError` が `?` で `EclipseError::Solver(_)` に透過ラップされる。
    #[test]
    fn solver_error_wraps_via_question_mark() {
        let r = solver_via_question();
        assert!(
            matches!(r, Err(EclipseError::Solver(SolverError::RootNotBracketed))),
            "expected Err(EclipseError::Solver(RootNotBracketed)), got {r:?}"
        );
    }

    /// `DomainError` が `?` で `EclipseError::Domain(_)` に透過ラップされる（NEW）。
    #[test]
    fn domain_error_wraps_via_question_mark() {
        let r = domain_via_question();
        assert!(
            matches!(
                r,
                Err(EclipseError::Domain(DomainError::OutOfRange {
                    what: "lat"
                }))
            ),
            "expected Err(EclipseError::Domain(OutOfRange{{lat}})), got {r:?}"
        );
    }

    // ============================================================
    // 受け入れテスト2: 透過 Display は内側の Display と一致（#[error(transparent)] 契約）
    // ============================================================

    /// `EclipseError::Time` の Display は内側 `TimeError` の Display と一致する。
    #[test]
    fn time_display_is_transparent() {
        let inner = TimeError::MissingLeapSecondData;
        let wrapped = EclipseError::from(inner.clone());
        assert_eq!(wrapped.to_string(), inner.to_string());
    }

    /// `EclipseError::Ephemeris` の Display は内側 `EphemerisError` の Display と一致する（NEW）。
    #[test]
    fn ephemeris_display_is_transparent() {
        let inner = EphemerisError::DataUnavailable;
        let wrapped = EclipseError::from(inner.clone());
        assert_eq!(wrapped.to_string(), inner.to_string());
    }

    /// `EclipseError::Solver` の Display は内側 `SolverError` の Display と一致する。
    #[test]
    fn solver_display_is_transparent() {
        let inner = SolverError::DidNotConverge;
        let wrapped = EclipseError::from(inner.clone());
        assert_eq!(wrapped.to_string(), inner.to_string());
    }

    /// `EclipseError::Domain` の Display は内側 `DomainError` の Display と一致する（NEW）。
    #[test]
    fn domain_display_is_transparent() {
        let inner = DomainError::OutOfRange { what: "lon" };
        let wrapped = EclipseError::from(inner.clone());
        assert_eq!(wrapped.to_string(), inner.to_string());
    }

    // ============================================================
    // 受け入れテスト3: source() 連鎖（#[from] が source() を自動実装）
    // ============================================================

    /// `Time` ラップは source() を返し、その Display は内側と一致する。
    #[test]
    fn time_wrap_has_source() {
        let inner = TimeError::MissingEarthOrientationData;
        let wrapped = EclipseError::from(inner.clone());
        let src = Error::source(&wrapped);
        assert!(src.is_some(), "expected source() to be Some, got None");
        assert_eq!(src.unwrap().to_string(), inner.to_string());
    }

    /// `Ephemeris` ラップは source() を返し、その Display は内側と一致する（NEW）。
    #[test]
    fn ephemeris_wrap_has_source() {
        let inner = EphemerisError::OutOfSupportedRange;
        let wrapped = EclipseError::from(inner.clone());
        let src = Error::source(&wrapped);
        assert!(src.is_some(), "expected source() to be Some, got None");
        assert_eq!(src.unwrap().to_string(), inner.to_string());
    }

    /// `Solver` ラップは source() を返し、その Display は内側と一致する。
    #[test]
    fn solver_wrap_has_source() {
        let inner = SolverError::NumericalInstability;
        let wrapped = EclipseError::from(inner.clone());
        let src = Error::source(&wrapped);
        assert!(src.is_some(), "expected source() to be Some, got None");
        assert_eq!(src.unwrap().to_string(), inner.to_string());
    }

    /// `Domain` ラップは source() を返し、その Display は内側と一致する（NEW）。
    #[test]
    fn domain_wrap_has_source() {
        let inner = DomainError::OutOfRange { what: "alt" };
        let wrapped = EclipseError::from(inner.clone());
        let src = Error::source(&wrapped);
        assert!(src.is_some(), "expected source() to be Some, got None");
        assert_eq!(src.unwrap().to_string(), inner.to_string());
    }

    // ============================================================
    // 受け入れテスト4: NotImplemented の存在・構築・返却（PATH, NEW）
    // ============================================================

    fn returns_not_implemented() -> Result<(), EclipseError> {
        Err(EclipseError::NotImplemented)
    }

    /// `EclipseError::NotImplemented` が存在し、関数から返せる。
    #[test]
    fn not_implemented_is_constructible_and_returnable() {
        let r = returns_not_implemented();
        assert!(
            matches!(r, Err(EclipseError::NotImplemented)),
            "expected Err(EclipseError::NotImplemented), got {r:?}"
        );
    }

    /// `NotImplemented` の Display は非空。
    #[test]
    fn not_implemented_display_is_non_empty() {
        let msg = EclipseError::NotImplemented.to_string();
        assert!(!msg.is_empty(), "NotImplemented Display must be non-empty");
    }

    // ============================================================
    // 受け入れテスト5: 日食固有 variant が維持される
    // ============================================================

    /// `DegenerateGeometry`/`InvalidFitInterval`/`EvaluationOutsideFitInterval` が
    /// 引き続き存在し構築可能（コンパイルレベル + matches!）。
    #[test]
    fn eclipse_specific_variants_preserved() {
        let degenerate = EclipseError::DegenerateGeometry;
        let invalid = EclipseError::InvalidFitInterval;
        let outside = EclipseError::EvaluationOutsideFitInterval;
        assert!(matches!(degenerate, EclipseError::DegenerateGeometry));
        assert!(matches!(invalid, EclipseError::InvalidFitInterval));
        assert!(matches!(
            outside,
            EclipseError::EvaluationOutsideFitInterval
        ));
    }

    // ============================================================
    // 受け入れテスト5b: BesselFitExceededTolerance 構築 + PartialEq
    // ============================================================

    /// 最も複雑な variant `BesselFitExceededTolerance { achieved, tolerance }` を
    /// 実 `BesselFitError` 値で構築でき、`matches!` で当該 variant に一致し、
    /// `PartialEq` が反射的（値はそのクローンと等しい）であることを固定する。
    /// この variant（f64 を含む `BesselFitError` を 2 つ載せる）が `Eq` を妨げるため、
    /// `PartialEq` 反射性をここで明示的に縛る。
    #[test]
    fn bessel_fit_exceeded_tolerance_constructs_and_is_reflexively_eq() {
        let achieved = BesselFitError {
            max_x: 1.0e-3,
            max_y: 2.0e-3,
            max_l1: 3.0e-3,
            max_l2: 4.0e-3,
        };
        let tolerance = BesselFitError {
            max_x: 1.0e-6,
            max_y: 1.0e-6,
            max_l1: 1.0e-6,
            max_l2: 1.0e-6,
        };
        let e = EclipseError::BesselFitExceededTolerance {
            achieved,
            tolerance,
        };
        assert!(
            matches!(e, EclipseError::BesselFitExceededTolerance { .. }),
            "expected BesselFitExceededTolerance variant, got {e:?}"
        );
        // PartialEq 反射性: 値はそのクローンと等しい（最も複雑な variant を pin）。
        assert_eq!(
            e,
            e.clone(),
            "BesselFitExceededTolerance must equal its clone"
        );
    }

    // ============================================================
    // 受け入れテスト6: non_exhaustive 前方互換（_ アーム要求）
    // ============================================================

    /// `EclipseError` は `#[non_exhaustive]` のため、（特に外部 crate からの）match では
    /// ワイルドカード `_` アームが必須。この match は現状コンパイルするが、契約として
    /// `_` アームを置くことで前方互換の意図を固定する（variant 追加で壊れないことを示す）。
    #[test]
    fn non_exhaustive_requires_wildcard_arm() {
        let e = EclipseError::DegenerateGeometry;
        let classified = match e {
            EclipseError::DegenerateGeometry => "degenerate",
            // `_` は #[non_exhaustive] のため REQUIRED（前方互換: 将来 variant 追加でも壊れない）。
            _ => "other",
        };
        assert_eq!(classified, "degenerate");
    }
}
