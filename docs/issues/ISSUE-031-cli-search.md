# ISSUE-031: CLI search（umbra search --from --to --accuracy）

- crate: umbra-cli
- 依存: ISSUE-016〜019（新月・合・候補・影円錐）, ISSUE-020〜023（基本面・ベッセル要素・多項式・全球分類＝`EclipseEngine::search` 経路）, ISSUE-007（時刻データ同梱）, ISSUE-001
- モード(tdd-workflow): standard（薄い CLI ラッパ。コア計算は umbra-eclipse が担保済みで、本 issue は引数解釈・出力整形が主。境界の堅牢性は要るが数式は持たないため standard）

## 目的
期間内の太陽食を列挙する CLI サブコマンド `umbra search` を実装する（architecture §1 umbra-cli, api-draft §3.2 `EclipseEngine::search`）。
- 引数: `--from <date> --to <date> --accuracy <fast|standard|reference>`。
- 出力: 各日食の event_key・種別・最大食時刻（UTC+TT）・gamma・最大食分・食面積（人間可読＋ `--format json`）。
- 「該当なし」はエラーにせず空リスト（api-draft §0/§3.2 `Result<Option<_>>` 方針）。

## 非目的
- 局地条件（ISSUE-032 `umbra local`）。
- 経路/GeoJSON（`umbra path`・別 issue）、inspect/bessel/validate/bench（architecture §1、別 issue）。
- search アルゴリズム本体（ISSUE-020〜023）。本 issue は CLI 境界。

## 公開インターフェース
architecture §8 / api-draft §3.2、CLI 仕様:
```
umbra search --from <DATE> --to <DATE> [--accuracy <fast|standard|reference>]
             [--format <text|json>] [--kind <all|total|annular|partial|hybrid>]
```
- `--from/--to`: 日付（RFC3339 or `YYYY-MM-DD`。UTC 既定、`UtcInstant::from_rfc3339`/`from_gregorian`・api-draft §1.3）。`UtcRange { start, end }`。
- `--accuracy`: `AccuracyProfile`（api-draft §3.1。既定 standard）。`EngineConfig::{fast,standard,reference}`。
- 内部: `standard_engine(TimeData::bundled())`（api-draft §3.2/§5）→ `engine.search(range)`。
- 出力（text）: 1 日食 1 行 + 詳細。出力（json）: `SolarEclipse` の serde（feature `serde`・api-draft §0）＋ `CalculationMetadata`（accuracy.md §0）。
```rust
// clap 由来の引数構造体（例）
pub struct SearchArgs { pub from: String, pub to: String, pub accuracy: AccuracyArg, pub format: FormatArg, pub kind: KindFilter }
pub fn run_search(args: SearchArgs) -> Result<(), CliError>;
```

## 数式・アルゴリズムの出典
- 本 issue は計算を持たない（CLI 境界）。列挙は `EclipseEngine::search`（architecture §3 データフロー、偽陰性なし方針）。
- 日付パース: RFC3339（`from_rfc3339`）/ グレゴリオ暦（`from_gregorian`・api-draft §1.3）。閏秒境界は TimeScales（ISSUE-006/007）。
- 出典明記対象は主にコア側。CLI は出力に `CalculationMetadata`（ephemeris/ΔT モデル・不確実性帯、accuracy.md §0）を必ず添える。

## 単位 / 時刻系 / 座標系
- 入力日付: UTC（既定・公開入出力は UTC・conventions §6）。
- 出力時刻: 最大食は **UTC+TT 両方**（accuracy.md §0）。将来日食は ΔT 不確実性帯を併記（metadata、accuracy.md §0/§2.3）。
- 角度: gamma 無次元、食分/食面積 無次元。表示は度（必要時）。
- 座標系: 最大食点は GeoPoint（測地緯度・東経・conventions §3）。

## アルゴリズム概要
1. 引数パース（clap）。日付→`UtcRange`、`--accuracy`→`EngineConfig`。
2. `standard_engine(TimeData::bundled())` 構築（同梱 EOP/閏秒/ΔT・accuracy.md §5、実行時ネットワーク禁止）。
3. `engine.search(range)` 実行。`--kind` でフィルタ。
4. 整形出力（text/json）。各日食に metadata（モデル・ΔT 不確実性帯）。
5. 該当なし → 空出力（エラーにしない）。
- 境界/堅牢性: 不正日付 → 明確なエラーメッセージ＋非0終了コード。`--from > --to` → エラー。範囲が暦対応外（1900〜2100 外・accuracy.md §6）→ `UnsupportedTimeRange` を分かりやすく提示。閏秒/EOP データ期限切れ → `MissingLeapSecondData`/`MissingEarthOrientationData` を案内（api-draft §6 期限切れ挙動・要確認）。

## 受け入れテスト
テストレベル **L5（全球日食）＋ CLI 統合**:
- 既知期間: 2024-01-01〜2024-12-31 で 2024 の日食件数・種別が NASA カタログ（data-sources §4、ISSUE-029 フィクスチャ）と整合（慣習を揃えて・accuracy.md §3.1）。
- 出力に最大食 UTC+TT 両方・gamma・食分・metadata（ΔT 不確実性帯）が含まれる（accuracy.md §0）。
- `--accuracy fast/standard/reference` で profile が metadata に反映。
- `--format json` が serde で妥当（パース可能・列挙タグ安定・api-draft §0）。
- 該当なし期間 → 空出力・終了コード0。
- 異常系: 不正日付/`from>to`/対応外範囲 → 非0終了＋明確メッセージ。
- 再現性: 同入力で決定的出力（実行時ネットワークなし・data-sources §0）。

## 許容誤差
CLI は計算許容を新設しない。出力値の許容はコア（ISSUE-020〜023）＝accuracy.md §2:
- 最大食時刻 ±1〜2s（TT 基準）、食分 ±0.0005。UTC は ΔT 律速（accuracy.md §0/§2.3、不確実性帯を出力）。
- CLI 自体の要件: 値を**改変・丸めで精度を捏造しない**（有効桁を metadata の精度に合わせる）。日付パースは閏秒境界で正しい（ISSUE-006）。
- 根拠: CLI は表示層。コア精度を素通しし、UTC/TT 分離と不確実性帯を必ず提示（accuracy.md §0）。

## 実装メモ
- clap で引数定義。`--accuracy` 既定 standard、`--format` 既定 text。
- 同梱データ（`TimeData::bundled()`）の取得 API は api-draft §6 未確定。bundled 前提で実装し、期限切れ挙動を案内（要確認）。
- 出力に必ず `CalculationMetadata`（モデル・ΔT 不確実性帯・accuracy.md §0）。将来日食は UTC 律速を明示。
- 実行時ネットワーク禁止（data-sources §0、accuracy.md §5）。
- レビュー重点: 日付パース堅牢性、UTC+TT 両方表示、ΔT 不確実性帯、終了コード/エラーメッセージ、json 安定性、対応年代外の扱い。
