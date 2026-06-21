# mutation レビュー: 南北限界線 厳密錐接線解（`solve_limit_edge` / `sample_central_point`）

対象: ISSUE-045 残(5)・M9.4（南北限界線の厳密錐接線解）。`crates/umbra-eclipse/src/engine.rs` の
`solve_limit_edge`（相対速度包絡の不動点反復）と `sample_central_point`（中心線＋限界点サンプル）。

## 実行
```
cargo mutants -p umbra-eclipse --file crates/umbra-eclipse/src/engine.rs \
  --re 'solve_limit_edge|sample_central_point' --no-shuffle -- --lib -- limits
```
（限界線核は `path()` の限界生成からのみ到達＝lib の `limits*` テスト〔錐exact 1e-7・包絡⊥ 1e-9・南北割当〕が
killer。非到達の他 lib テストを省いて高速化。）

## 結果（2026-06-21）
**65 mutants: 57 caught・5 missed・3 unviable・0 timeout。**

- caught 57: 相対速度 rel の各項（x'/y'/μ' 係数・cos d/sin d・ζ/η/ξ）、半径 `|l2−ζ·tan f2|`、オフセット
  `±radius·n̂`、法線 `(−rel_y,rel_x)/|rel|`、南北割当、`t_hours` 時間スケール（`rigorous_bessel` の x 二次化で
  load-bearing 化）等の算術・比較を全捕捉。`< → >`（収束判定の反転＝発散）も caught。
- unviable 3: 戻り値を `Ok(Some((Default::default(), …)))` 等へ置換する変異。`GeoPoint` は `Default` 非実装で
  コンパイル不能＝unviable（偽の生存ではない）。

## 生存 5 件＝等価変異（許容・収束判定の微細構造）
すべて不動点反復の収束判定 `step_converged`（engine.rs 773–775）に対する変異:

| # | 変異 | 等価の理由 |
|---|---|---|
| 1,2,3 | `< → <=`（ξ/η/ζ の各収束比較） | 反復は機械精度（残差 ~1e-15）まで過収束し、差分が許容 `LIMIT_FIXED_POINT_TOL=1e-12` に**ちょうど一致する一点**は浮動小数では到達不能。境界 `<`/`<=` の差は出ない。 |
| 4,5 | `&& → ||`（3 成分 AND→OR） | ξ,η,ζ は同率で収束する。反復1回目は 3 成分とも ≫TOL（オフセット ~0.01 Re）、2回目で 3 成分とも <TOL に同時到達するため、AND と OR は**同じ反復で break**。仮に 1 反復早く止めても、その時点で点は cone/⊥ 条件を test 許容（1e-7/1e-9）の遥か内側で満たすため検出不能。 |

これは既存の粗走査/求根機構（`descending_sign_change_bracket` / `solve_zero_in_window` / `scan_point_count`）
と同カテゴリ＝**解像度・収束機構の微細構造は、過収束する真解に影響しない等価変異**。`< → >`（発散）のような
**意味を変える**変異は caught で、収束の*検出*境界のみが等価。よって**許容**し、CI mutation gate からは
`--exclude-re 'with <= in.*solve_limit_edge'` / `--exclude-re 'with || in.*solve_limit_edge'` で除外する
（`solve_limit_edge` 内の `<`/`&&` は収束判定にしか現れないため、この 2 パターンは当該等価変異のみを的確に
捕捉し、算術・`<→>` 等の可殺変異は除外しない）。

実回帰ガード: 通常 CI の `cargo test -p umbra-eclipse`（FAST 合成 μ'≠0 の前方射影 2 条件＋SLOW 実 2024-04-08
皆既の NASA 帯幅 [185,215]km・2 条件）が限界線の正しさを縛る。

---

## 追記: M9 残(3) 3a 錐半径引数化の再走（2026-06-21）

`solve_limit_edge` は M9.6 で `engine.rs` から `axis_intercept.rs` へ移設済み。M9.7 残(3) 3a で**錐半径を
引数化**（`cone_l`/`cone_tan_f` を追加し、本影 `(l2,tan f2)`／半影 `(l1,tan f1)` の両縁を解けるよう一般化）。
半径以外は不変。新たに axis_intercept.rs 内 FAST 単体テスト（半影・本影縁の錐exact＋包絡⊥を前方射影で機械精度
検証・半影/本影の半径差・sign 両側・地表外し/相対速度ゼロの `Ok(None)`）を killer に追加。

### 実行
```
cargo mutants -p umbra-eclipse --re 'solve_limit_edge' --no-shuffle \
  --exclude-re 'with <= in.*solve_limit_edge' --exclude-re 'with \|\| in.*solve_limit_edge' \
  -- --lib -- penumbral umbral_edge sign_gives edge_off zero_relative
```

### 結果
**52 mutants: 51 caught・1 unviable・0 missed。** 半径 `cone_l − ζ·cone_tan_f` の引数化後も全算術・比較・
sign・rel 項を新テストが捕捉（生存 0）。収束判定の等価 2 件（`<→<=`・`&&→||`）のみ除外（上表と同根拠）。

### ⚠ CI 除外パターンの正規表現バグ修正（重要）
旧 mutation.yml は `--exclude-re 'with || in.*solve_limit_edge'`（および `…central_width_and_duration`）と
**`||` を未エスケープ**で書いていた。`||` は正規表現で**空 alternation＝全文字列マッチ**になるため、cargo-mutants
の `--exclude-re` が**全 mutant を除外** → CI mutation ジョブが「Found 0 mutants／No mutants found」で**事実上
no-op**（生存を検出できない）だった。リテラル `||`（mutant 名 "replace && with ||"）を除外するには
**`\|\|` とエスケープ必須**。本スライスで mutation.yml の該当 2 パターンを `\|\|` へ修正（`solve_limit_edge` /
`central_width_and_duration`）。M9.7 で追加した `sample_central_point` 用パターンは当初から `\|\|` でエスケープ済み。
この修正で CI mutation が実際に変異をテストするようになる（本再走の 52 mutants がその実証）。
