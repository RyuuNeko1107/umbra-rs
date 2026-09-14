# mutation レビュー: `search` の範囲フィルタ（ISSUE-050）

対象: `crates/umbra-eclipse/src/engine.rs` の `EclipseEngine::search` 末尾の `retain`
（最大食時刻が閉区間 `[start, end]` に入る日食だけを残す）。正本は
`docs/issues/ISSUE-050-search-range-filter.md`。

## 実行

```
docker compose -p umbra-rs run --rm rust bash -c \
  "cargo install cargo-mutants --locked && \
   cargo mutants --package umbra-eclipse --re 'EclipseEngine<.*>::search' \
     -j 2 --timeout 700 -- --test search_range_filter"
```

killer を `--test search_range_filter`（SLOW 5 本）に限定した。`umbra-eclipse` 全スイートは約 530 秒で
baseline が timeout に掛かるため（`mutation-path-sample-limit.md` と同じ制約）。

## 結果（2026-09-15）

| 総数 | caught | missed | unviable | 所要 |
|---|---|---|---|---|
| 13 | 8 | **4** | 1 | 44 分 |

### 本 issue が導入した行（`retain` の判定・247 行）は **全て caught**

| 変異 | 結果 |
|---|---|
| `&&` → `\|\|`（両端判定の結合） | caught |
| 開始側 `>= 0.0` → `< 0.0` | caught |
| 終了側 `>= 0.0` → `< 0.0` | caught |

`search_range_filter.rs` の 5 本（再現・閉区間の両端・1 秒外・P1〜P4 重なりのみ・広範囲の不変条件）が
判定式を厳密に縛れていることの確認になっている。

### missed 4 件は**本 issue の変更範囲外**（既存コード）

| 変異 | 位置 | 判定 |
|---|---|---|
| `/` → `%` / `*`（161 行） | 候補窓の走査に関わる既存の算術 | 本走の killer（範囲フィルタ 5 本）は候補生成を縛っていない。**全スイート（453 本）側で担保**されるべき箇所で、本 issue で新たに生じた穴ではない |
| `*` → `+`（166 行） | 同上 | 同上 |
| match arm `(Some(p1), Some(p4))` の削除（201 行） | 全球接触の組み立て | 同上 |

**誤解を避けるための明示**: これらは「ISSUE-050 のテストが弱い」のではなく、
**killer を 5 本に絞った走なので、絞った範囲の外は当然撃てない**というだけである。
`search` 全体の mutation を全スイート killer で回すには baseline の timeout 対策
（テストの FAST/SLOW 分離）が要る＝`mutation-cli-path.md` と同じ残作業。

## 実回帰ガード

通常 CI の `cargo test -p umbra-eclipse --test search_range_filter`（5 本）が、範囲判定の
criterion（最大食時刻）・境界（閉区間）・順序・範囲外除外を縛る。
