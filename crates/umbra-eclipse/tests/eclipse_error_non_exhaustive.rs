//! 外部 crate からの `EclipseError` match で `#[non_exhaustive]` 契約を実検証する。
//!
//! 統合テストは本体クレートとは**別クレート**としてコンパイルされる。この境界では
//! `#[non_exhaustive]` 列挙体に対する `match` はワイルドカード `_` アームが**コンパイル必須**
//! となる（既知 variant を網羅しても `_` が無ければ E0004 で落ちる）。in-crate のユニットテスト
//! では `#[non_exhaustive]` は自クレートに対して効かないため `_` を省いてもコンパイルでき、
//! 「`_` が要求される」契約を真には強制できない。よってこのファイルこそが前方互換契約
//! （将来 variant を足しても外部の網羅 match が壊れない）の実検証である。

use umbra_eclipse::EclipseError;

/// 外部 crate から `EclipseError` を match する。既知 variant を 2 つ明示し、末尾に
/// `_ =>` アームを置く。この `_` は `#[non_exhaustive]` かつ外部クレートであるため REQUIRED
/// （削ると E0004 でコンパイルが通らない）。これが前方互換契約の実コンパイル検証。
#[test]
fn external_crate_match_requires_wildcard_arm() {
    let e: EclipseError = EclipseError::DegenerateGeometry;
    let classified = match e {
        EclipseError::DegenerateGeometry => "degenerate",
        EclipseError::NotImplemented => "not-implemented",
        // `_` は外部クレート + #[non_exhaustive] のため REQUIRED（前方互換契約の核）。
        _ => "other",
    };
    assert_eq!(classified, "degenerate");
}
