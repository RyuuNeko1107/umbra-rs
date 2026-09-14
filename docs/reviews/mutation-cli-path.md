# mutation レビュー: `umbra path` CLI（ISSUE-048）

対象: `crates/umbra-cli/src/lib.rs` の `run_path` / `format_path_text` / `format_limit_line` /
`check_path_sample_count`。正本は `docs/issues/ISSUE-048-cli-path.md`。

## 実行

```
docker compose -p umbra-rs run --rm rust bash -c \
  "cargo install cargo-mutants --locked && \
   cargo mutants --package umbra-cli \
     --re 'run_path|format_path_text|format_limit_line|check_path_sample_count' -j 2 --timeout 600"
```

## 結果（部分実行・2026-09-14）

**38 変異中 6 件を実行した時点で中断**（下記「中断の理由」）。実行分は **全て caught**（うち 5 件は TIMEOUT）。

| 変異 | 結果 |
|---|---|
| `check_path_sample_count -> Ok(())`（検査を無効化） | **TIMEOUT** |
| `check_path_sample_count` の `span / interval` → `%` / `*` | **TIMEOUT** |
| `check_path_sample_count` の `>` → `==` / `<` | **TIMEOUT** |

`check_path_sample_count` が生む変異は上記 5 件で**全て**であり、**いずれも TIMEOUT＝ハング検出で caught**。
これは本ガードの存在意義そのものの証明になっている: **ガードを壊すと `run_path` が実際に停止しなくなる**
（`run_path_tiny_interval_is_interval_too_small` が返らなくなる）。本プロジェクトで timeout を caught として
扱う既存方針（ループ制御の変異＝`mutation-limb-bulge.md` / `mutation-polygon-union.md`）と同じ扱い。

## 中断の理由（誤魔化さずに記録）

`umbra-cli` のテストスイートは**実エンジンの SLOW テスト**（1 日範囲の `search` を複数回）を含むため
baseline だけで約 236 秒かかる。cargo-mutants は変異ごとに全スイートを回すので、38 変異の完走には
**3〜4 時間**かかる見積りだった（timeout 到達分は 1 件 600 秒）。マシン負荷を優先して中断した
（`docker-light` の運用方針）。

**未実行の 32 件**は `run_path` / `format_path_text` / `format_limit_line` の変異で、
**text 整形のラベル・件数・分岐**が中心。これらは ISSUE-048 節の 20 テスト
（`--no-limits` で「なし」・`None` 要素の明示・geojson の `role` 集合一致・実 2024 の種別/点数）が
**通常の `cargo test` で**縛っている。ただし **mutation で撃てたことは確認できていない**ので、
「テストは存在するが判別力は未証明」の状態である。

**残作業**: CLI のテストを FAST（合成 fixture）と SLOW（実エンジン）に分離し、mutation は FAST のみを
killer にして完走させる。現状は実エンジン依存のため mutation の費用が実用域を超えている。
