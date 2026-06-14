# ISSUE-046: xtask 骨子・CI・cargo-deny/audit 運用（生成物 checksum 差分検査の機械化）

- crate: xtask（バイナリ crate・workspace ローカル）＋ CI（GitHub Actions 等のワークフロー定義）
- 依存: ISSUE-001（workspace 雛形・lint・deny.toml 骨子・CI 骨子）, ISSUE-033（VSOP87 係数生成＝generate/verify の利用側）, ISSUE-034（ELP/MPP02 係数生成＝同）, ISSUE-040（章動係数 IAU2000A データ管理＝同）, ISSUE-007（EOP/閏秒/ΔT データの versioned/checksum 管理）, umbra-core（**`DataSetMetadata` は umbra-core の純粋型・確定 B3**。xtask とライブラリ本体で共有）
- milestone: M2 暦（係数パイプライン 033/034/040 と同時に必要。骨子は M0/M1 から、checksum 差分検査は M2 で実効化）
- モード(tdd-workflow): strict（生成物の**再現性**・provenance/checksum・**ライセンス機械チェック**（cargo-deny/audit）は data-sources §0/§5/§6 の法的整合と全データ層の信頼性に直結。CI ゲートは公開仕様級の運用契約。strict）

## 目的
`xtask` の**サブコマンド分岐の骨子**と、**CI** による生成物 checksum 差分検査・ライセンス機械チェックを定義する（architecture §11, data-sources §5/§6, milestone0-review §Minor「046 xtask骨子・CI実体・cargo-deny/audit運用」）。
- **xtask 骨子**: 033/034/040 が使うサブコマンド分岐（`generate-coefficients` / `verify-generated` など）と、各データセット共通の **`DataSetMetadata` 共有型**（name/version/source/license/valid_from/valid_to/checksum、architecture §11）を一箇所に集約する。
- **generate/verify の契約**: `generate` は一次データ→packed 生成物を**決定的**に作る。`verify` は generated/ と再生成物の **checksum 差分**を検出する（architecture §11, data-sources §5）。
- **CI**: 生成物 checksum 差分検査（コミット済み generated と xtask 再生成の不一致を fail）、`cargo-deny`（ライセンス allow-list・GPL 非混入・thiserror 等の許可）、`cargo-audit`（脆弱性）。
- **data-sources §6 チェックリストの CI 化**: ライセンス allow-list 機械チェック・GPL 派生物非混入・依存とデータのライセンス整合を CI ゲートにする。

## 非目的
- 各係数の**生成ロジック本体**（VSOP87=ISSUE-033、ELP/MPP02=ISSUE-034、章動 IAU2000A=ISSUE-040）。本 Issue は**サブコマンド分岐・共有型・CI ゲート**の骨子で、各 generate の中身は各 Issue。
- ライブラリ本体の機能テスト・差分テスト（ISSUE-029/030/036）。本 Issue は生成物の再現性とライセンス運用。
- データの一次取得（ネットワーク）をライブラリ本体に入れること（**禁止**。取得は xtask 隔離・accuracy.md §5, data-sources §4 冒頭）。
- 暫定オラクル戦略・二段ゲートの定義（ISSUE-047）。本 Issue はその CI 実行基盤を提供しうるが、方針は 047。

## 公開インターフェース
xtask は内部ツール（公開 API ではない）。CLI サブコマンド契約と共有型を定義:
```text
cargo xtask generate-coefficients [--dataset vsop87|elp-mpp02|nutation-iau2000a|all]
cargo xtask verify-generated      [--dataset ...]   # generated/ と再生成の checksum 差分を検査（差分あれば非0終了）
cargo xtask verify-data           # EOP/閏秒/ΔT の valid_to 期限・checksum 検査（accuracy.md §5）
cargo xtask check-licenses        # cargo-deny + データ NOTICE/provenance 整合（data-sources §6）
```
```rust
/// 全データセット共通の出所・完全性メタデータ（architecture §11, data-sources §5）。
/// 置き場: **umbra-core の純粋型（確定 B3）**。xtask・ライブラリ本体（TimeData 等, ISSUE-042）で共有。
#[derive(Clone, Debug, /* serde gate */)]
pub struct DataSetMetadata {
    pub name: String, pub version: String,
    pub source: String,            // 一次配布元・引用（data-sources §2/§3）
    pub license: String,           // ライセンス区分（GPL 派生物は不可・§0）
    pub valid_from: String, pub valid_to: String,
    pub checksum: String,          // 生成物の決定的ハッシュ
}
```
- CI ワークフロー（リポジトリ `.github/workflows/` 等）:
  - `fmt` / `clippy`（-D warnings）/ `test`（ISSUE-001 の lint 骨子を実効化）。
  - `verify-generated`: 生成物が再現可能か（差分 fail）。
  - `cargo-deny check`（advisories/bans/licenses/sources）+ `cargo-audit`。
  - `check-licenses`: data-sources §6 チェックリストの機械化部分。

## 数式・アルゴリズムの出典
- 本 Issue は**運用・ツール骨子**で数式を持たない。
- 規約・運用の出典:
  - architecture §11（コードと係数分離、`cargo xtask generate-coefficients` / `verify-generated`、CI で生成済みとの差分確認、`DataSetMetadata`）。
  - data-sources §0（GPL 派生物非混入）/ §5（取り込みパイプラインと完全性）/ §6（OSS 公開前チェックリスト＝CI 化対象）。
  - accuracy.md §5（実行時ネットワーク禁止・更新は xtask 隔離・versioned+checksum）。
  - conventions §11（出典不明データ／magic number 禁止）。

## 単位 / 時刻系 / 座標系
- 本 Issue は数値計算を持たない（単位・時刻系・座標系の直接の対象外）。
- データ完全性: checksum はバイト列に対する決定的ハッシュ（アルゴリズムは固定し `DataSetMetadata.version` に紐付け）。
- `valid_from`/`valid_to`: 各データの対応年代（data-sources §1）。期限切れ検知は `verify-data` で。

## アルゴリズム概要
1. `xtask` バイナリに `clap` 等でサブコマンド分岐（generate-coefficients / verify-generated / verify-data / check-licenses）。033/034/040 はこの分岐の `--dataset` 配下に generate 実体を登録する。
2. `DataSetMetadata` 共有型を **umbra-core の純粋型として定義（確定 B3）**。各 generate（xtask）と TimeData（ISSUE-042）が共有し provenance＋checksum を付与。
3. `verify-generated`: 一次データ→再生成し、コミット済み generated/ と checksum を比較。不一致なら非0終了（CI fail）。決定的生成が前提（data-sources §5）。
4. `cargo-deny`: `deny.toml`（ISSUE-001 骨子）に **licenses allow-list**（MIT / Apache-2.0 等）と **bans**（GPL 系・取り込み不可）を定義。`thiserror` 等の採用クレートを allow（milestone0-review A6）。`cargo-audit` で脆弱性。
5. data-sources §6 チェックリストのうち機械化可能な項（ライセンス allow-list 検査・GPL 非混入・依存/データ整合）を CI ジョブ化。手動確認項（一次配布元の利用条件）はチェックリストとして残す。
6. CI で全ジョブをゲート化（main へのマージ条件）。巨大 DE データを要する差分テストは CI 必須にしない（nightly/手動・accuracy.md §3.1）。

## 受け入れテスト
- サブコマンド分岐: `cargo xtask --help` が generate-coefficients/verify-generated/verify-data/check-licenses を列挙。未知 dataset はエラー。
- 決定性: `generate-coefficients` を2回実行して生成物の checksum が一致（決定的、data-sources §5）。
- 差分検出: generated/ を意図的に1バイト改変すると `verify-generated` が**非0終了**（CI fail を再現）。
- `DataSetMetadata` 充足: 033/034/040 の生成物が name/version/source/license/valid_from/valid_to/checksum を**非空**で持つ（provenance 欠落を fail）。
- ライセンスゲート: `cargo-deny check licenses` が allow-list 外を fail。GPL 系データ/コードの混入を bans で fail（data-sources §0）。`cargo-audit` が既知脆弱性で fail。
- 期限検知: `verify-data` が EOP/閏秒の `valid_to` 超過を検出（accuracy.md §5・期限切れ挙動）。
- チェックリスト機械化: data-sources §6 のうち CI 化した項目が実際にジョブとして走る（メタテスト）。
- ネットワーク禁止確認: ライブラリ本体ビルド/テストがネットワークアクセスを要しない（取得は xtask のみ）。
- 二段オラクルゲート（ISSUE-047 連動）: **本 Issue は運用基盤のため数値ゲート対象外**。ただし 047 の M2 暫定/M10 最終ゲートを **CI で実行する基盤**（feature `jpl` ジョブの分離・nightly 化）を提供する。

## 許容誤差
- 本 Issue は数値計算を持たないため**許容誤差なし**。生成物は **bit 完全一致**（checksum 不一致＝fail。近似許容なし）。

## 実装メモ
- **milestone0-review §Minor 確定事項**の反映: xtask 骨子・CI 実体・cargo-deny/audit 運用。033/034/040 が乗るサブコマンド分岐と `DataSetMetadata` 共有型を本 Issue で固定する。
- `DataSetMetadata` の置き場は **umbra-core の純粋型（確定 B3）**。xtask とライブラリ本体（ISSUE-042 TimeData）で共有する。ISSUE-033 の依存行はこれに合わせる。ライブラリ本体に xtask 専用ロジックは持ち込まない（型のみ core 共有）。
- `cargo-deny` allow-list には milestone0-review A6 の `thiserror` を含める。bans に GPL 系（ytliu0/ElpMpp02 派生等・data-sources §0/参照）を明示。
- 巨大データ（JPL DE・ISSUE-036）の差分テストジョブは CI 必須にせず nightly/手動（accuracy.md §3.1, data-sources §2.3）。feature `jpl` ジョブを分離。
- data-sources §6 の手動確認項（一次配布元の再配布可否）は機械化できないため、チェックリストを CI 出力に残し、OSS 公開前ゲート（plan §26）に接続する。
- レビュー重点: 生成の決定性、checksum 差分の確実な fail、ライセンス bans/allow-list の網羅、ネットワーク隔離、`DataSetMetadata` 共有の一貫性。
