# mutation レビュー: 部分食域 3 領域ユニオン（umbra-eclipse 側・`trace_terminator_lunes` / `cone_terminator_intersections_detailed` / `build_partial_limit`）

対象: ISSUE-045 M9 残(3) サブスライス **(3f)** の umbra-eclipse 側。umbra-geo 側（`union_rings`）は
`docs/reviews/mutation-polygon-union.md`。正本は `docs/algorithms/11-path-partial-domain.md` §11.6。

- `axis_intercept.rs` `cone_terminator_intersections_detailed`: terminator 交点を基本面 ξ 付きで返す（morning/evening 分類の供給源）。
- `engine.rs` `trace_terminator_lunes`: ξ 符号で morning/evening に分け、緯度 hi/lo の lockstep リボンで 2 lune を組む。
- `engine.rs` `build_partial_limit`: 帯 ∪ morning ∪ evening を `union_rings` で合成し最大面積成分を返す。

## 実行

```
# (A) 統合: path_limits（FAST 合成 bessel ＋ SLOW 実 2024）を killer にする
docker compose -p umbra-rs run --rm rust bash -c \
  "cargo install cargo-mutants --locked && \
   cargo mutants --package umbra-eclipse \
     --re 'trace_terminator_lunes|cone_terminator_intersections_detailed|build_partial_limit' \
     -j 2 --timeout 600 -- --test path_limits"
# (B) 単体: detailed の算術は axis_intercept の --lib 単体テストで縛る
docker compose -p umbra-rs run --rm rust bash -c \
  "cargo mutants --package umbra-eclipse --re 'cone_terminator_intersections' \
     --exclude-re 'in scan_periodic_sign_change_roots' -j 1 --timeout 300 \
     -- --lib -- two_intersections d_zero cone_not_reaching cone_radius deterministic terminator_ellipse"
```

## 結果（2026-09-14）

| 走 | 総数 | caught | missed | timeout | unviable |
|---|---|---|---|---|---|
| (A) 統合（path_limits） | 58 | 39 | **9** | 3 | 7 |
| (B) 単体（`--lib`・detailed 40 変異） | 40 | 36 | **0** | 0 | 4 |

(A) の missed 9 のうち **5 は `cone_terminator_intersections_detailed` の算術**（`299: - → +`、`301: / → %`・`/ → *`・`* → /`、
`310: / → %`＝楕円係数 k・√k・η=sinθ/√k）で、統合テストは実 2024・合成 bessel とも扁平率の効果（~21 km）に鈍感なため
撃てないが、**(B) の `--lib` 単体（往復オラクル・d=π/2 二円閉形式・楕円条件）が全て撃つ**。通常 CI（全テスト）では caught。
→ **実質 missed は `trace_terminator_lunes` の 4 件**（下記）。

timeout 3 は `trace_terminator_lunes` のループ制御（`839: || → &&`、`<= → >`、`+ → *` 系）＝break 条件が消えて無限ループ
＝ハング検出で caught 扱い（(3c-i)/(3c-iii) と同カテゴリ）。

`build_partial_limit` は 3 変異: `Ok(None)` 固定・`.len() >= 3` → `<` は caught、`Ok(Some(Default))` は unviable。

## 生存 4 件（`engine.rs:829` `(*xi < 0.0) == is_morning` の分類式）— **除外**

| 変異 | 判定 | 理由 |
|---|---|---|
| `< → >`、`== → !=` | **等価（公開 IF 上）** | morning/evening の**ラベルが入れ替わる**だけ。2 lune は `build_partial_limit` で**両方とも**ユニオンに入り、ユニオンは対称なので `partial_limit` は不変。ラベルは関数外に出ない |
| `< → <=` | **等価（測度ゼロ）** | `ξ = 0` ちょうどの交点だけ帰属が変わる。§11.6(a) が「ξ=0 は evening」と定めた strict 不等号の境界で、連続量ゆえ実現しない |
| `< → ==` | **非等価・受入不感** | `ξ == 0.0` が常に偽になり全交点が evening 側へ入る（morning lune 空・evening lune が両 limb を跨ぐ 1 本のリボン）。両 limb の hi/lo を結ぶリボンは各時刻の半影 ∩ terminator 弦を覆い、頂点集合は同一（閉半影内・昼面側の頂点契約は保たれる）。中心線全点内包・帯>皆既帯・ζ≈0 頂点出現のいずれでも差が出ず、**実 2024 SLOW 含む全受入で不可分** |

`< → ==` を撃つには「lune の面積・形状」を独立オラクルで縛る必要があるが、部分食域の許容は ballpark（数十 km・
§許容誤差）で、外周は 3 領域の**和**として消費されるため、個々の lune の形状を固定するのは過仕様になる
（(3c-iii) の `&& → ||` 除外と同じ判断・`docs/reviews/mutation-limb-bulge.md`）。以上 4 件は
`mutation.yml` で `--exclude-re 'in trace_terminator_lunes'` の**関数内比較演算子**に限定して除外する
（`replace (<|==) with .* in.*trace_terminator_lunes`。ループ制御の `||`/`<=`/`+` は timeout で caught のまま残す）。

## 実回帰ガード

通常 CI の `cargo test -p umbra-eclipse`（`--lib` の `cone_terminator_intersections*` 単体＋`path_limits` 27 本・
実 2024 SLOW の中心線 194 点全点内包）が、分類の供給源・lune 組立・ユニオン消費の実挙動を縛る。
