# ISSUE-001: プロジェクト規約と単位系（workspace 雛形・lint・CI 骨子）

- crate: umbra-core（ただし workspace ルートも対象）
- 依存: なし
- モード(tdd-workflow): standard（コードは雛形中心で trivial 寄りだが、規約・lint・deny.toml は公開仕様と全 crate の前提を固定するため standard とする）

## 目的
Cargo workspace の雛形と、全 crate が従う規約基盤を確立する。具体的には:
- `architecture.md §1` のクレート構成（umbra-core / -ephemeris / -eclipse / -geo / -cli / -fixtures + xtask）の空クレートと workspace `Cargo.toml`。
- `rustfmt.toml` / `clippy` 設定（lint レベルの統一）。
- `deny.toml`（cargo-deny: ライセンス allow-list と advisory）。
- CI 骨子（fmt / clippy / test / cargo-deny / cargo-audit / verify-generated の各ジョブ枠）。
- 各 crate の lib ルートに `conventions.md` / `accuracy.md` / `data-sources.md` への参照コメント（module doc）を置き、実装コメントから規約を参照する運用（conventions §10）を仕込む。

## 非目的
- 実際の数式・型の実装（002 以降）。
- 係数生成パイプラインの実装（xtask の中身。data-sources §5。本 issue は空コマンド枠のみ）。
- ライセンスの法的最終確定（data-sources §6 のチェックリストは別工程。ここでは deny.toml の allow-list 雛形まで）。

## 公開インターフェース
本 issue は型を公開しない。確立する構造物:
- `Cargo.toml`（workspace、`[workspace.package]` で version/edition/license を共有）。
- 各 crate の `src/lib.rs` に module doc コメント（規約 doc 参照）。
- `rustfmt.toml`, `clippy` 設定（`[lints]` テーブル or `lib.rs` の `#![deny(...)]`）。
- `deny.toml`。
- `.github/workflows/*.yml`（or 同等 CI 定義）の骨子。
- `cargo xtask` のサブコマンド枠: `generate-coefficients` / `verify-generated`（architecture §11）。

## 数式・アルゴリズムの出典
数式なし。準拠する社内文書:
- conventions.md §1（単位表）, §10（実装コメント運用）, §11（禁止事項）。
- architecture.md §1（クレート構成）, §11（データ管理・xtask）。
- data-sources.md §0/§6（ライセンスリスクと cargo-deny 要件）。
- 外部ツール: `cargo-deny`（EmbarkStudios）, `cargo-audit`（RustSec）。

## 単位 / 時刻系 / 座標系
- 該当なし（雛形）。ただし `conventions.md §1` の単位表を「正本」として lib doc に転記参照する。生 f64 を単位付き量として公開 API に出さない方針（conventions §1 / api-draft §0）を `#![warn]` 運用ポリシーとして明記。

## アルゴリズム概要
1. workspace `Cargo.toml` を作成。`resolver = "2"`、`[workspace.package]` に共有メタ（license は MIT OR Apache-2.0 を仮置き＝data-sources §6 で確定するまで「要確認」コメント）。
2. 6 crate + xtask の空 lib/bin を生成。各 `lib.rs` の module doc に conventions/accuracy/data-sources への相対パス参照を記載。
3. `rustfmt.toml`（max_width 等は既定踏襲、明示）。clippy は workspace `[lints.clippy]` で `all` + `pedantic` の必要分。`unwrap_used` / `panic` を計算経路で warn。
4. `deny.toml`: licenses allow-list（MIT/Apache-2.0/BSD 系）、GPL を明示 deny（data-sources §0 の GPL 非混入要件）。advisories で RustSec を有効化。
5. CI 骨子: fmt → clippy → test → deny → audit → `cargo xtask verify-generated` の各ジョブを雛形定義（中身は後続 issue で埋める）。
- 禁止事項（規約レベル・conventions §11）を CI/lint で機械的に効かせる足場を置くこと。magic number / 無条件 Newton は後続 issue で個別にガード。

## 受け入れテスト
accuracy.md のテストレベル該当: なし（インフラ）。代わりに CI ゲートで検証:
- `cargo build --workspace` が通る（空クレートでも）。
- `cargo fmt --check` / `cargo clippy --workspace -- -D warnings` が通る。
- `cargo deny check` が通り、GPL ライセンスのダミー依存を追加すると **fail する**こと（deny の有効性を 1 ケースで確認。テスト後に削除）。
- `cargo xtask --help` が `generate-coefficients` / `verify-generated` を表示する。
- 異常系: deny.toml の allow-list に無いライセンスの crate を入れると CI が落ちる（境界確認）。

## 許容誤差
数値誤差なし。accuracy.md にバジェット項目なし（インフラ層）。
- 「厳密一致」基準: lint 違反ゼロ、deny check ゼロ違反を合格条件とする（根拠: conventions §11 の禁止事項を CI で漏れなく弾くため、警告許容ゼロが妥当）。

## 実装メモ
- license フィールドは data-sources §6 が未解決のため `# TODO(要確認): data-sources §6` コメントを必ず残す（推測で確定しない）。
- `umbra-fixtures` は通常依存に含めない（architecture §1）。workspace member には入れるが、各 crate の `[dependencies]` には入れず `[dev-dependencies]` のみ許容する構成にする。
- xtask はライブラリ本体と分離（実行時ネットワーク禁止の境界。accuracy.md §5 / data-sources §4 冒頭）。
- レビュー重点: deny.toml の GPL deny が ELP/MPP02 GPL 再実装（data-sources §0/§2.2）混入を実際に防げる設定か。
