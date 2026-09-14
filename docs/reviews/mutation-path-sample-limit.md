# mutation レビュー: `path()` のサンプル数上限（ISSUE-049・`check_path_sample_limit`）

対象: `crates/umbra-eclipse/src/engine.rs` の `check_path_sample_limit`。
正本は `docs/issues/ISSUE-049-path-sample-limit.md`。CLI 側（ISSUE-048）は
`docs/reviews/mutation-cli-path.md`。

## 実行

```
docker compose -p umbra-rs run --rm rust bash -c \
  "cargo install cargo-mutants --locked && \
   cargo mutants --package umbra-eclipse --re 'check_path_sample_limit' \
     -j 2 --timeout 500 -- --test path_limits"
```

killer は `crates/umbra-eclipse/tests/path_limits.rs` の ISSUE-049 節 9 本（合成 fixture・FAST）。
**テスト対象を `--test path_limits` に限定**した（`umbra-eclipse` 全スイートは約 530 秒で
baseline が 400 秒 timeout に掛かり「unmutated tree で失敗」になるため。path_limits 単体は約 103 秒）。

## 結果（2026-09-14）

| 総数 | caught | missed | timeout | unviable |
|---|---|---|---|---|
| 12 | 10 | **0** | 2 | 0 |

**生存 0**。timeout 2 件はいずれも**ガードを壊すと実際にハングする**ことの検出で、caught 扱い
（本プロジェクトの既存方針＝ループ制御の変異は timeout で検出。`mutation-limb-bulge.md` /
`mutation-polygon-union.md` と同じ）。

| timeout した変異 | 壊れ方 |
|---|---|
| `span_of` の `days_since(...) * SECONDS_PER_DAY` → `/` | span が 86400² 分の 1 に潰れ、推定サンプル数が上限を下回って**素通り**→ 走査が終わらない |
| 部分食域の match guard `options.include_limits` → `false` | P1〜P4 が検査対象から外れ、`include_limits=true` の極小 interval が**素通り**→ 走査が終わらない |

この 2 件が timeout になること自体が、**検査対象 span の選び方（§確定仕様 2）が実際に効いている**
ことの証明になっている（U1〜U4 だけを見ていたら部分食域の走査を止められない）。

## 実回帰ガード

通常 CI の `cargo test -p umbra-eclipse --test path_limits`（ISSUE-049 節 9 本）が、
上限の境界（厳密 `>`）・NaN・非正・走査区間なし・P1〜P4 の採用を縛る。
