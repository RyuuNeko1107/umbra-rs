---
milestone: M11
---

# ISSUE-048: `umbra path` CLI サブコマンド（経路・GeoJSON 出力）

`EclipseEngine::path`（M9 完了・中心線／南北限界線／部分食域／`to_geojson`）を CLI へ露出する。
ISSUE-031 が「経路/GeoJSON（`umbra path`・別 issue）」として先送りした分で、architecture §1 の
`umbra-cli: search/local/path/...` を満たす。

## 目的

指定日の日食について経路を計算し、**text**（人間可読サマリ）または **geojson**（`EclipsePath::to_geojson`
の生出力）で表示する。計算は既存 API の呼び出しのみで、**新しい天文計算は行わない**。

## 非目的

- 経路計算そのものの変更（ISSUE-045 / M9 で完了済み）。
- `bessel` / `inspect` / `validate` / `bench` サブコマンド（別 issue）。
- ファイル出力（`--output`）。標準出力のみ（リダイレクトで足りる）。

## 公開インターフェース

```rust
pub struct PathArgs {
    pub date: String,                 // YYYY-MM-DD（UTC）。当日に起こる日食の経路
    pub accuracy: AccuracyArg,        // 既定 standard（search/local と同一）
    pub format: PathFormatArg,        // 既定 text
    pub interval: f64,                // サンプル間隔 [s]（既定 PathOptions::default）
    pub no_limits: bool,              // 指定時 include_limits=false
}

pub enum PathFormatArg { Text, Geojson }

pub fn run_path(args: &PathArgs) -> Result<String, CliError>;
```

`Command::Path(PathArgs)` を追加する。

## 確定仕様

1. **日付解決**: `--date` を既存 `parse_date` で `UtcInstant` に。その日 00:00 から**翌日 00:00 まで**を
   `UtcRange` として `search` し、**最初に見つかった日食**を対象とする（`run_local` の日付解決と同型）。
2. **日食が無い日は「成功」**（`run_local` と同一方針＝エラーにしない・架空の値を作らない）:
   - `text`: `No solar eclipse on {date}.` の 1 行。
   - `geojson`: **空の `FeatureCollection`**（`{"type":"FeatureCollection","features":[]}`）。
     `null` は GeoJSON として妥当でないので返さない。空 FeatureCollection は「feature が無い」ことの
     正確な表現であり、誤差・未提供を隠していない（conventions §11）。
   **`CliError::NoEclipse` は追加しない**（当初案を撤回）。
3. **エンジン構築**: `search`/`local` と同一（`AccuracyArg::Standard` は `standard_engine`、
   `Reference` は `EngineConfig::reference()` で構築）。同梱データのみ・実行時ネットワークなし。
4. **`PathOptions`**: `sample_interval_seconds = --interval`、`include_limits = !--no-limits`、
   `split_antimeridian` は `PathOptions::default()` の値を維持（経路計算側の規約・(3g) は GeoJSON 出力層）。
   `--interval` が**非正・非有限**なら `CliError::InvalidInterval`（`search` を呼ぶ前に fast-fail）。

   **サンプル数の上限（実装レビュー指摘・2026-09-14）**: 正の有限値でも、極端に小さい `--interval`
   （例 `1e-300`）は `path` の走査回数を事実上無限にし、**CLI が無反応のままハングする**。よって
   日食を特定した後・`path` を呼ぶ前に、`P1`〜`P4`（`partial_begin`/`partial_end`）の秒数を
   `--interval` で割った**推定サンプル数**が `MAX_PATH_SAMPLES = 100_000` を超えるなら
   `CliError::IntervalTooSmall { interval, estimated_samples }` を返す。
   既定 60 s では実日食で数百点なので、上限は通常利用を制限しない。
   `P1`/`P4` が `None`（起こり得ない想定）のときは検査を行わない（捏造した span で判定しない）。
   **本質的な露出はエンジン側（`EclipseEngine::path` は任意の呼び出し元から同じ入力を受ける）**にあり、
   CLI 境界での防御は暫定。エンジン側の入力検証は後続（本 issue の非目的）。
5. **text 出力**（決定的・1 行 1 項目）。**ラベルは下表で固定**する（テストが縛る契約）:

   | ラベル | 内容 | 要素が `None` のとき |
   |---|---|---|
   | `種別` | `SolarEclipseKind` | （常に有る） |
   | `最大食` | 最大食の UTC 時刻と地点（緯度経度） | （常に有る） |
   | `中心線` | 点数と始点・終点 | `なし` |
   | `北限` | 点数 | `なし` |
   | `南限` | 点数 | `なし` |
   | `部分食域` | 環数と外環頂点数 | `なし` |
   | `samples` | 件数（`None` ではなく 0 件で表現） | （`0`） |

   **`None` の要素は行を省略せず「なし」と明示**し、空の成功出力を作らない
   （conventions §11「誤差・未提供を隠さない」）。`--no-limits` は `include_limits=false` を通じて
   北限・南限・**部分食域**を `None` にし、`samples` を空にする（M9.7 lockstep 契約）。
6. **geojson 出力**: `EclipsePath::to_geojson()` の結果をそのまま印字（末尾に改行）。**再整形しない**
   （二重の整形規約を持たない）。
7. **エラー**: 不正日付・非正/非有限 interval は `search` 前に fast-fail、過小 interval は `path` 前に fast-fail。エンジンエラーは `CliError::Eclipse` で透過。
   部分食で中心線が無い場合も**成功**（text は「中心線: なし」）。

## 受け入れテスト戦略

FAST（合成・実エンジンだが 1 日範囲）:
- 不正日付・非正 interval・非有限 interval が **fast-fail**（エンジンを呼ばない）。
- 極端に小さい interval が `IntervalTooSmall` で弾かれ、**`path` を呼ばない**（ハング防止）。
- 日食の無い日が**成功**し、text は 1 行の告知・`geojson` は**空 FeatureCollection**。
- `--format geojson` が**妥当な JSON**で `FeatureCollection`、`role` の集合が `EclipsePath` の
  Some 要素と一致。
- `--no-limits` で北限/南限/部分食域が「なし」になり、`samples` が空になる（M9.7 の lockstep 契約）。
- text が **`None` 要素を「なし」と明示**する（空欄にしない）。

SLOW（実日食 1 件）: 2024-04-08 で text に皆既・中心線点数 > 0 が出る／geojson が MultiPolygon を含みうる。

## 依存

ISSUE-043（Engine 結線・完了）、ISSUE-045（path 本体・M9 完了）、ISSUE-031/032（CLI 基盤・完了）。

## 完了状況（2026-09-14・strict）

**完了**。`umbra path`（text / geojson）を実装し、architecture §1 の `umbra-cli: search/local/path/...` のうち
`path` を満たした。ISSUE-031 が「別 issue」として先送りした分の回収。

### 確定した仕様からの追補

- **サンプル数上限 `MAX_PATH_SAMPLES = 100_000`**（§4）: 実装レビューが「正の有限値でも極端に小さい
  `--interval`（例 `1e-300`）は走査回数が事実上無限になり **CLI が無反応のままハングする**」ことを指摘。
  日食特定後・`path` 呼び出し前に P1〜P4 の秒数から推定サンプル数を検査し、超過なら
  `CliError::IntervalTooSmall`。検査は純関数 `check_path_sample_count(Option<f64>, f64)` に分離し、
  境界（厳密 `>`）・`None` 分岐・エラー payload を単体テストで固定した。
  **本質的な露出はエンジン側**（`EclipseEngine::path` は任意の呼び出し元から同じ入力を受ける）にあり、
  CLI 境界の防御は暫定。**エンジン側の入力検証は後続**（本 issue の非目的）。
- **`CliError::NoEclipse` は追加せず**、日食の無い日は成功（text は 1 行の告知・geojson は空
  `FeatureCollection`）。`run_local` と同一方針。

### 検証

`crates/umbra-cli/src/lib.rs` の ISSUE-048 節に **20 テスト**（crate 合計 65・全通過）:

| 区分 | 本数 | 内容 |
|---|---|---|
| fast-fail | 5 | 不正日付／`interval` が 0・負・NaN・±∞（いずれも**エンジンを呼ばない**） |
| 日食なし | 3 | text の告知行・geojson の空 FeatureCollection・過小 interval でも**告知が優先**（検査順序） |
| 出力契約 | 5 | geojson の `role` 集合が `EclipsePath` の `Some` 要素と一致／`--no-limits` で北限・南限・部分食域が「なし」かつ `samples` 0／部分食で `None` 要素が行ごと「なし」と明示／実 2024 の種別・最大食・中心線点数 |
| ハング防止 | 5 | 過小 interval の `IntervalTooSmall`／閾値上は成功／既定 interval は不採択にならない／span が **P1〜P4** であること（bit-exact） |
| 純関数単体 | 4 | 境界が厳密 `>`（2 進厳密な span/interval で固定）／1 ULP 上は拒否／`None` span は `Ok`／payload が `interval` と推定値 |

fmt / clippy `-D warnings` / 全テスト通過（Docker 内）。
