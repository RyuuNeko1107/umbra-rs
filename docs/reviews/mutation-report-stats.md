# ミューテーション生存変異の許容判断 — `umbra-fixtures::report`（ISSUE-030 S30a/S30b）

`cargo mutants --package umbra-fixtures --file report.rs`
（cargo-mutants 27.1.0, Docker 内）。

- **S30a**: 35 mutants → 32 caught / 3 unviable / 0 missed（境界単体テスト追加後）。
- **S30b**（compare_global / aggregate_global 追加後の全体）: **46 mutants: 41 caught / 5 unviable / 0 missed**。
- **S30c**（compare_local / aggregate_local / contact_time_error_seconds 追加後の全体）:
  **69 mutants: 62 caught / 7 unviable / 0 missed**。可視性カウント検査は match/mismatch 取り違えを
  確実に撃破するため非対称（false 3 / true 1）に強化。
- **S30d**（report_against_golden / GoldenComputer 追加後の全体）:
  **72 mutants: 64 caught / 8 unviable / 0 missed**。実エンジン SLOW 統合テスト
  （`report_against_golden_real_engine_one_golden`）は各変異で再実行すると非現実的なため
  `cargo mutants … -- -- --skip real_engine` で除外（オーケストレーションはモック注入の高速テストが縛る）。

対象（純粋・検証基盤）:
- `ErrorStats::from_errors`（絶対誤差の記述統計 max|e| / mean|e| / p95）。
- `ErrorStats::within`（`max_abs <= tolerance` の合否境界）。
- 非公開 `percentile_r7_sorted`（R-7 線形補間パーセンタイル）。
- `ToleranceProfile::standard()` / `reference()`（accuracy.md §2 の許容定数）。
- `compare_global`（全球条件の符号付き誤差 computed−golden。時刻は days_since×86400）。
- `aggregate_global`（metric 別 ErrorStats 集計＋ greatest/magnitude ゲートの pass。γ 非ゲート）。
- `contact_time_error_seconds`（TT 優先・days_since・2要素保持の共有時刻誤差ヘルパ）。
- `compare_local`（地点別。接触 C1–C4 の Some/None 対・presence 不一致計上・可視性一致）。
- `aggregate_local`（接触フラット集計＋ 7 ゲート pass・可視性/presence 不一致カウント）。
- `report_against_golden`（GoldenComputer 注入のオーケストレーション。found/missing/locations 計数・集計）。

## caught（統計・契約のコア）

- **abs 適用**（符号相殺なし）・**mean = sum/n**・**max = 昇順末尾**・**ソート有無**は
  `from_errors_*`（基本/負値/未ソート 11 点）が捕捉。
- **R-7 補間の本体**（`a[lo] + (h-lo)·(a[lo+1]-a[lo])` の `+`/`-`/添字、`h=(n-1)·p`）は
  既知 p95（2.9 / 9.5）と n=1 分岐テストが捕捉。
- **`within` の境界 `<=`**（`<`/`>=`、max↔mean 取り違え）は inclusive 境界テストが捕捉。
- **許容定数 12 値**（standard/reference 各 6）は exact フィールド検証＋
  「reference は standard より厳格」テストが捕捉。

## 初回生存 3 件 → 追加テストで撃破（0 missed）

初回 run で `percentile_r7_sorted` の**上端ガード分岐**に 3 件生存:
- `report.rs:88` `lo + 1 >= n` の `+`→`*`
- `report.rs:89` true 分岐 `sorted_abs[n - 1]` の `n-1`→`n+1` / `n/1`

原因: 公開経路 `from_errors` は **p=0.95 固定**で、`lo = ⌊(n-1)·0.95⌋` は常に `lo+1 < n`
（`0.95 ≥ 1` が偽ゆえ `lo == n-1` に到達しない）。よってガード true 分岐
`if lo+1 >= n { sorted_abs[n-1] }` は公開 API 経由では**到達不能**で、変異が生存した。

許容ではなく**撃破**を選択（防御分岐は一般 p のため残す）: 非公開関数を直接呼ぶ単体テスト
`report::tests::percentile_r7_boundaries_direct`（p=1.0→最大 / p=0.0→最小 / p=0.5→中央）を追加し、
ガードと true 分岐を実効化。再 run で **0 missed**（到達時に変異は OOB/誤値となり caught）。

## unviable 3 件

コンパイル不能な変異（型不整合等）。生存ではない。

## 結論

`umbra-fixtures::report` は**生存変異ゼロ**（等価/許容扱いの残置なし）。tdd-workflow 工程7 充足。
