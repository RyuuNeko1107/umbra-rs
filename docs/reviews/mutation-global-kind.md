# Mutation review — `umbra-eclipse::global::classify_global_kind` / `l2_changes_sign`（ISSUE-043 S6b-iii）

`cargo mutants` の生存変異の列挙と許容判断。退行ガード `mutation.yml` は eclipse を含む。

## 結果

- `classify_global_kind`（最大食基本種別＋中心食区間 [U1,U4] の l2 符号反転による Hybrid 上書き）と
  `l2_changes_sign`（[U1,U4] の l2 符号走査）の**挙動ロジック**は全 caught。
  - 走査区間演算 `t0 + (t1 − t0)·frac` の `− → /`（`t1/t0 ≈ 1` で走査が ≈[U1, U1+1日] にずれる）は、
    **純皆既で [U1,U4] 内は l2<0・区間外で l2>0** になる合成源
    （`classify_global_total_not_hybrid_when_l2_negative_in_central_interval`）が
    「正しい [U1,U4] でのみ Total、ずれた区間では Hybrid 誤判定」を縛って**撃破**。
  - 内部符号反転（金環-皆既-金環・両端同符号）の合成ハイブリッド
    （`classify_global_hybrid_interior_crossing`）が、端点だけでは検出できない＝サンプル密度を
    load-bearing 化し、SAMPLES 削減系の変異を撃破。
- 生存（許容）= **`l2_changes_sign` の符号判定境界 `l2 > 0` / `l2 < 0` の `>→>=` / `<→<=`** のみ。

## 生存（許容）= l2==0 ちょうどの境界等価変異

- `if l2 > 0.0` の `> → >=`、`if l2 < 0.0` の `< → <=`（global.rs `l2_changes_sign`）。
- 差が出るのは**サンプル点で l2 がちょうど 0.0** のときだけ。l2(t) は時刻の連続関数で、走査格子点が
  厳密に零点に当たるのは測度ゼロ（実 ephemeris・合成源いずれでも起きない）。よって挙動不変＝等価。
- `delta_t_seconds` / `BesselianPolynomial::at` / `assess_eclipse_possibility` の `< → <=`
  境界等価変異と同カテゴリ（接合点・境界一点のみで差・到達不能）。`mutation.yml` で
  `--exclude-re 'with >= in.*l2_changes_sign'` / `'with <= in.*l2_changes_sign'` 除外する。
- l2==0（皆既↔金環の連続境界）でどちらのフラグも立てない現挙動は意図どおり（厳密 0 は反転と見なさない）。

種別の実契約（Total/Annular/Hybrid/Partial/None・中心食のみ Hybrid 判定・純皆既/金環を Hybrid 誤判定
しない・内部反転の検出）は通常 CI の `cargo test`（`global::tests::classify_global_*` 群）と本体の
全 caught 変異が担保する。
