# ISSUE-044: EclipseError 集約・From 変換（thiserror・各層エラーの統合）

- crate: umbra-eclipse（`EclipseError`）／ 基盤エラー（`DomainError`/`TimeError`/`SolverError`）は umbra-core（ISSUE-001/007/008/012 由来）
- 依存: ISSUE-001（`DomainError`）, ISSUE-006/007（`TimeError`）, ISSUE-008（`SolverError`）, ISSUE-012（`EphemerisError`）, ISSUE-041（定数 crate との依存方針整合）, api-draft §1.6/§3.5
- milestone: M9（中核公開 API のエラー集約。エンジン結線 043 と同時期に確定）
- モード(tdd-workflow): **strict**（**確定A6**: thiserror 採用・公開エラー列挙・From 変換は SemVer 公開境界。誤ると下流の `?` 連鎖が壊れるため strict）

## 目的
**確定A6**（milestone0-review.md）に従い、各層のエラーを集約する公開エラー型 `EclipseError` を **thiserror** で実装する（api-draft §3.5/§6, architecture §10）。
- `EclipseError`（非網羅 enum）を thiserror で定義し、`Display` / `std::error::Error` を実装。
- `From<TimeError>` / `From<EphemerisError>` / `From<SolverError>` / `From<DomainError>` を `#[from]` + `#[error(transparent)]` で実装し、各層から `?` で集約（透過ラップ版 ERR、api-draft §3.5）。
- 基盤エラー（`DomainError`/`TimeError`/`SolverError`/`EphemerisError`）も `Display`/`Error` を満たす（thiserror、api-draft §1.6/§2）。
- `cargo-deny` allow-list に thiserror を追加（確定A6、data-sources §6）。

## 非目的
- エラーを発生させる各層のロジック（時刻 006/007・solver 008・暦 012・幾何 019-028）。本 Issue は**型集約と変換**のみ。
- エラーメッセージの国際化（v0.1 は英語 Display）。
- `anyhow` 等のアプリ側エラー（ライブラリは型付きエラーを返す。CLI 031/032 が必要なら別途）。
- 各層のエラー発生ロジック。本 Issue は型集約のみ。なお `path()` の v0.1 未実装は **panic でなく `Err(EclipseError::NotImplemented)`**（PATH/ISSUE-043/045）であり、本 Issue はその `NotImplemented` variant を提供する。

## 公開インターフェース
api-draft §3.5 / §1.6 を転記・**確定A6**（thiserror）で具体化:
```rust
// umbra-core 基盤エラー（thiserror, Display/Error 実装）
#[non_exhaustive] #[derive(Debug, thiserror::Error)]
pub enum DomainError { #[error("out of range: {what}")] OutOfRange { what: &'static str } }

#[non_exhaustive] #[derive(Debug, thiserror::Error)]
pub enum TimeError {
    #[error("invalid date")] InvalidDate,
    #[error("missing leap second data")] MissingLeapSecondData,
    #[error("missing earth orientation data")] MissingEarthOrientationData,
}
#[non_exhaustive] #[derive(Debug, thiserror::Error)]
pub enum SolverError {
    #[error("root not bracketed")] RootNotBracketed,
    #[error("did not converge")] DidNotConverge,
    #[error("numerical instability")] NumericalInstability,
}
#[non_exhaustive] #[derive(Debug, thiserror::Error)]
pub enum EphemerisError {
    #[error("out of supported range")] OutOfSupportedRange,
    #[error("data unavailable")] DataUnavailable,
    #[error("io error")] Io(/* 軽量化: feature jpl 無効時 std::io を引かない, ISSUE-012 */),
}

// umbra-eclipse 集約エラー（thiserror, 透過ラップ版に一本化 ERR）
// 下位エラー（Time/Ephemeris/Solver/Domain）は #[from] + #[error(transparent)] でラップ。
// 同義の二重 variant は作らない。直 variant は日食固有（下位と重複しないもの）のみ。
#[non_exhaustive] #[derive(Debug, thiserror::Error)]
pub enum EclipseError {
    // 下位エラーは透過ラップ（#[from]）。時刻/EOP/閏秒は Time、求根は Solver、定義域は Domain、暦は Ephemeris に集約。
    #[error(transparent)] Time(#[from] TimeError),            // InvalidDate / MissingLeapSecondData / MissingEarthOrientationData
    #[error(transparent)] Ephemeris(#[from] EphemerisError),  // OutOfSupportedRange / DataUnavailable / Io
    #[error(transparent)] Solver(#[from] SolverError),        // RootNotBracketed / DidNotConverge / NumericalInstability
    #[error(transparent)] Domain(#[from] DomainError),        // 観測者/範囲などの定義域違反

    // 日食固有（下位に無い失敗のみ直 variant）
    #[error("degenerate eclipse geometry")] DegenerateGeometry,
    #[error("besselian fit exceeded tolerance")] BesselFitExceededTolerance,
    #[error("solver root not bracketed")] RootNotBracketed,
    /// 未実装機能（PATH: path() 等。panic でなくこれを返す。UnsupportedTimeRange を流用しない）。
    #[error("not implemented")] NotImplemented,
}
```
- `#[non_exhaustive]`（前方互換、api-draft §0）。`From<TimeError/EphemerisError/SolverError/DomainError>` を `#[from]` で実装し透過ラップ（ERR）。
- **確定（ERR）**: `MissingLeapSecondData`/`MissingEarthOrientationData`/`InvalidDate` は `From<TimeError>` 経由に一本化（直 variant にしない）。`DidNotConverge`/`NumericalInstability` は `From<SolverError>` 経由。`InvalidObserver`/定義域違反は `From<DomainError>` 経由。`UnsupportedTimeRange` は「対応年代外」専用語義に保ち、**未実装には流用しない**（未実装は `NotImplemented`、PATH）。`RootNotBracketed` は確定 ERR の「日食固有・重複しない直 variant」列挙に含まれるため直 variant として残置するが、下位 `SolverError::RootNotBracketed`（#[from] 透過ラップ）と語義が重なりうる。**要確認**: 日食固有のブラケット失敗を別契機として残すか、下位由来に一本化するか（同義二重 variant 回避の原則と整合させ実装レビューで確定）。

## 数式・アルゴリズムの出典
- 数式なし（型設計）。出典は API 規約のみ:
  - thiserror 採用: 確定A6（milestone0-review.md）、api-draft §6 未確定「thiserror で実装するか」の確定。
  - `std::error::Error` / `Display` 実装: Rust API ガイドライン（エラー型の慣習）。
  - `non_exhaustive`・From 変換: api-draft §0/§3.5（前方互換・`?` 集約）。

## 単位 / 時刻系 / 座標系
- 該当なし（エラー型）。エラーは値・単位を持たない（メッセージは英語 Display）。
- 時刻系/座標系の不整合自体はエラーの**契機**だが、本 Issue は型のみ。

## アルゴリズム概要
1. umbra-core に `DomainError`/`TimeError`/`SolverError`（と ISSUE-012 の `EphemerisError`）を thiserror で定義（Display/Error）。
2. umbra-eclipse に `EclipseError` を thiserror で定義、`#[from]` + `#[error(transparent)]` で `TimeError`/`EphemerisError`/`SolverError`/`DomainError` を透過ラップ（同義二重 variant は作らない）。直 variant は日食固有のみ（`DegenerateGeometry`/`BesselFitExceededTolerance`/`RootNotBracketed`/`NotImplemented`）。
3. 各層（042 時刻・043 エンジン・solver・暦）は `?` で `EclipseError` へ自動変換。
4. `Display` メッセージは原因が分かる英文。`source()`（Error トレイト）で連鎖を辿れる（`#[from]` が自動実装）。
5. `cargo-deny` allow-list に thiserror を追加（A6、data-sources §6）。
- 設計原則（ERR）: 透過ラップ版に一本化。同義二重 variant を作らず、下位由来は `#[from]`（transparent）で伝播、直 variant は日食固有のみ（下記）。

## 受け入れテスト
accuracy.md テストレベル **L1（型・契約）**。
- **From 変換（ERR 必須）**: `TimeError`/`EphemerisError`/`SolverError`/`DomainError` が `?` で `EclipseError` に透過ラップ変換される（コンパイル＋実行テスト）。
- **NotImplemented**: `EclipseError::NotImplemented` が存在し、未実装機能（path() 等、PATH/ISSUE-045）がこれを返せる（`UnsupportedTimeRange` を流用しない）。
- **Display/Error**: 各エラーが `Display` を実装し、`std::error::Error::source()` で原因連鎖を返す（`#[from]` variant）。
- **non_exhaustive**: 外部 crate からの match で `_` アームが要求される（前方互換、コンパイルテスト）。
- **下流統合**: ISSUE-042（valid_to 超過 → `MissingLeapSecondData`/`MissingEarthOrientationData`）・ISSUE-043（solver/暦エラー）が `EclipseError` に集約されて伝播する。
- **依存方針**: `cargo-deny` allow-list に thiserror が載り、ライセンス機械チェックを通る（data-sources §6）。
- **軽量化**: `EphemerisError::Io` が feature `jpl` 無効時に std::io 重依存を引かない契約（ISSUE-012 と整合、要確認項目をテストで固定）。

## 許容誤差
- 該当なし（型集約。数値許容誤差なし）。
- 「エラーを握り潰さない・誤差を隠さない」原則（conventions §11, accuracy.md §0）と整合: valid_to 超過や solver 不収束を沈黙させず型で表面化する。

## 実装メモ
- **確定（ERR/A6）厳守**: thiserror 採用、**透過ラップ版に一本化**（`#[from]` + `#[error(transparent)]` で TimeError/EphemerisError/SolverError/DomainError をラップ）、`Display`/`std::error::Error`。`cargo-deny` allow-list に追加（data-sources §6 / ISSUE-046 CI）。
- **variant 重複の撤廃（確定 ERR）**: `MissingLeapSecondData`/`MissingEarthOrientationData`/`InvalidDate` は `From<TimeError>` 経由に一本化し直 variant を作らない（旧二重定義を撤廃）。`DidNotConverge`/`NumericalInstability` は `From<SolverError>` 経由、`InvalidObserver`/定義域違反は `From<DomainError>` 経由（ISSUE-042 の valid_to 超過挙動＝`Missing*Data` と整合）。直 variant は日食固有（`DegenerateGeometry`/`BesselFitExceededTolerance`/`RootNotBracketed`/`NotImplemented`）のみ。`UnsupportedTimeRange` は「対応年代外」専用に保ち未実装に流用しない。
  - **要確認**: 直 `RootNotBracketed` と下位 `SolverError::RootNotBracketed`（透過ラップ）の語義重複。日食固有の別契機として残すか下位由来に一本化するか、実装レビューで確定。
  - **実装時確定（実装レビュー済み）**: 直 `RootNotBracketed` は**設けない**（下位 `SolverError::RootNotBracketed` を `Solver` で透過ラップして表現。既存呼出 `EclipseError::Solver(SolverError::RootNotBracketed)` と整合・同義二重 variant 回避）。
- **実装時確定（透過ラップの実体）**: 4 層ラップは `#[error(transparent)]` ではなく **`#[error("{0}")]`** を採用。`#[error(transparent)]` は `Display` と `source()` の**両方**を内側へ委譲するため、leaf な基盤エラー（`TimeError` 等は `source()`＝None）では §受け入れテスト「source() で原因連鎖を返す」が満たせない。`#[error("{0}")]` は `Display` を内側に委譲（メッセージ非二重化＝「透過」意図）しつつ `#[from]`（＝`#[source]`）で内側を 1 段の cause として `source()` 公開し、Display 一致と source 連鎖の**両アクセプタンスを充足**する。
- **実装時確定（直 variant の集合）**: 日食固有の直 variant は `DegenerateGeometry` / `BesselFitExceededTolerance{achieved,tolerance}`（ISSUE-022・実測残差保持で誤差を隠さない §11）/ `InvalidFitInterval` / `EvaluationOutsideFitInterval`（ISSUE-022/037）/ `NotImplemented`。`#[non_exhaustive]` のため後続追加は前方互換。
- **実装時確定（Io 繰延）**: `EphemerisError::Io` は feature `jpl` 本実装（ISSUE-012/036・M10）まで未追加。現 `EphemerisError` は `OutOfSupportedRange`/`DataUnavailable` のみで、`Ephemeris(#[from])` ラップは現状の variant で From 要件を満たす。
- `EphemerisError::Io` は ISSUE-012 の「feature `jpl` 無効時 std::io 非依存」方針と整合（軽量表現、要確認）。
- 基盤エラー（DomainError/TimeError/SolverError）は umbra-core、集約 `EclipseError` は umbra-eclipse（crate 帰属、architecture §10）。
- レビュー重点: thiserror、透過ラップ版一本化（From の網羅: Time/Ephemeris/Solver/Domain）、二重 variant 撤廃、`NotImplemented` 追加（PATH）、non_exhaustive、source() 連鎖、cargo-deny allow-list、Io 軽量化。
