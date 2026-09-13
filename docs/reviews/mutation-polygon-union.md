# mutation レビュー: 多角形ユニオン（`umbra_geo::union_rings` / `clip.rs`）

対象: ISSUE-045 M9 残(3) サブスライス **(3f)**（部分食域＝帯 ∪ morning lune ∪ evening lune の多角形ユニオン）。
`crates/umbra-geo/src/clip.rs` 全体。正本は `docs/algorithms/11-path-partial-domain.md` §11.6 / §11.7。

## 実行

```
docker compose -p umbra-rs run --rm rust bash -c \
  "cargo install cargo-mutants --locked && \
   cargo mutants --package umbra-geo --file crates/umbra-geo/src/clip.rs -j 2"
```

killer は `crates/umbra-geo/tests/union.rs`（FAST・29 本）。合成多角形の**厳密な頂点列・面積・環数・向き**で縛る
（座標表に依存しない・機械精度）。

## 結果

| 走 | 総数 | caught | missed | timeout | unviable |
|---|---|---|---|---|---|
| 初回（21 テスト） | 280 | 230 | **43** | 6 | 1 |
| 判別力追加後（29 テスト）＋ dead code 除去 | 275 | 234 | **35** | 5 | 1 |

追加で撃てた 8 件は、三角形（頂点ちょうど 3）のリングを入力・出力・**穴**に持ち込むテスト群
（`union_regression_single_triangle_*` / `*_two_overlapping_triangles_*` / `*_triangle_and_rectangle_*` /
`*_three_strips_form_triangular_hole`）が `pip` と `normalize` の `< 3` 判定を固定したもの。
あわせて、到達不能だった `while diff > TAU` の正規化ループ（`back`・`atan2` がともに `(-π, π]` ゆえ
第 1 ループ後に必ず `(0, TAU]` へ収まる＝上側折返しは起こり得ない）を**死コードとして削除**した
（missed 3 + timeout 1 が消滅）。

## timeout（＝ハング検出で caught 扱い）

| 変異 | 挙動 |
|---|---|
| `drop_collinear` の `delete !`（`if !removed`） | 反復停止条件が消えて無限ループ |
| `union_rings` の `back - …` → `/`、`while diff <= 0` → `>`、`diff += TAU` → `-=` / `*=` | 角度正規化が収束せず無限ループ |

いずれも「壊れたら止まらなくなる」ことをテスト実行時間が検出している。

## 生存変異の分類（許容可否）

### A. 等価変異 — 許容（31 件）

| 系統 | 変異 | 等価の理由 |
|---|---|---|
| ray-casting の向き | `pip:101` `>` → `<` | even-odd は**右向き半直線と左向き半直線で交差数の偶奇が一致**する（全直線の交差数が偶数）。判定結果は不変 |
| 厳密等値の境界 | `pip:101` `>` → `>=`、`intersection_points:152/161` `>` → `>=`、`probe_sides:239` `>` → `>=`、`drop_collinear:208` `>` → `>=` | 閾値と**厳密に等しい**ときだけ差が出る測度ゼロ境界。しかも境界判定点は法線オフセットで境界から外してある |
| 到達不能な退化ガード | `point_segment_distance:177` / `union_rings:298` の `rr > 0.0` → `>=`、`drop_collinear:205/206` の `l1/l2 > 0.0` → `>=` | `normalize` が連続重複点を潰すため**長さ 0 の辺は存在しない**。ガードは防御的で到達しない |
| 死んだ分岐 | `union_rings:420` `a2 > 0.0` → `>=` | 直前で `a2.abs() <= AREA2_EPS` を除外済みゆえ `a2 == 0.0` に到達しない |
| 正準化の向き | `union_rings:321` `ka <= kb` → `>` | 無向端点対の**代表の選び方**が反転するだけ。写像は依然として一意（重複除去の結果は不変） |
| 自己ペアの追加 | `union_rings:283` `(i + 1)..` → `(i * 1)..` | `intersection_points(e, e)` は共線分岐で自分の端点を返すのみ＝既にある分割点 0/1 と同じ |
| 末尾重複点の除去 | `normalize:117` `v.len() > 1` → `==` / `>=` | 閉リング入力で末尾の重複点が残るが、長さ 0 の辺は断片を生まない（`key` 一致で捨てられる）。`union_is_invariant_to_closed_or_open_input_rings` が示すとおり結果は不変 |
| 許容幅のスケール | `drop_collinear:207` `*` → `/`（×2）、`probe_sides:236` `*` → `+` / `/`、`probe_sides:240` `*` → `/` | **解像度・許容幅の係数**を変える変異。共線判定は厳密共線で `cross ≈ 0`（~1e-16）・非共線で桁違いに大きいので採否が変わらない。境界プローブ幅は「中点から他の辺までの最短距離」項が支配するため、辺長側の係数を変えても左右の membership は変わらない。既存の `descending_sign_change_bracket` / `scan_periodic_sign_change_roots` と**同カテゴリ**（`docs/reviews/mutation-axis-intercept.md` / `mutation-rise-set.md`） |
| 同値の tie-break | `union_rings:393` `diff < d` → `<=` | `diff` が完全一致するのは**同一方向の出辺が 2 本ある**場合だが、無向重複除去でそれは残らない |

### B. 生存・**未解決**（4 件・許容するがリスクを明示）

| 変異 | 位置 | 判定 |
|---|---|---|
| `back - atan2(…)` → `back + atan2(…)` | `clip.rs:387`（環再結合の後継半辺選択） | **等価ではない**（分岐頂点での選択が変わりうる）。判別力を狙って追加した非対称分岐テスト（`union_regression_three_regions_meeting_at_one_point_stay_separate` / `*_asymmetric_point_contact_rectangle_and_triangle`）でも、**実 2024-04-08 の SLOW（`real_2024_eclipse_partial_limit_is_plausible`・ガード付き手動変異で確認）でも撃てなかった** |
| `signed_area2(o).abs() < …` → `>` / `<=` / `==` | `clip.rs:445`（穴を含む最小面積の外環へ割り当てる比較） | **等価ではない**。ただし差が出るのは**穴が 2 つ以上の外環に含まれる二重入れ子**（穴の中の島が、さらに穴を持つ）配置のみ。現行テスト・実 2024 のいずれもこの構造を作らない |
| `param_on_segment -> false` ほか同関数内の算術 7 件 | `clip.rs:129–135`（共線重なりの分割点登録） | **等価ではない**（共線部分重なりの分割が失われる）。狙って追加した B 系テスト（`union_regression_collinear_partial_overlap_split_in_its_interior` ほか）でも撃てなかった。現行の全配置では、必要な分割点が**非平行交差**または**リング頂点**として別経路で供給されているため差が出ない |

**許容判断**: いずれも「誤った結果を返しうる経路が存在するが、現行の受け入れオラクル（合成 29 本＋実 2024 SLOW）が
到達しない」もの。部分食域の許容誤差は ballpark（数十 km・algorithms §許容誤差）で、`partial_limit` は
最大面積成分のみを使うため実害は限定的だが、**`union_rings` は `pub` API なので契約としては未証明**である。

**残作業（要確認へ登録）**:
1. 分岐頂点の後継選択を撃つ判別テスト（4 本以上の境界半辺が非対称な角度で集まり、選択を誤ると成分分割が
   変わる配置の構成）。
2. 二重入れ子（穴の中の島がさらに穴を持つ）での穴割当テスト。
3. 共線部分重なりの分割が**唯一の**供給源になる配置の構成。

以上が構成できない場合は、該当ロジックを「到達不能ゆえ削除」か「不変条件で防御」のどちらかに倒すべきで、
現状の「実装はあるが証明されていない」状態を放置しない。

## 実回帰ガード

通常 CI の `cargo test -p umbra-geo`（union 29 本・geojson 21 本）と
`cargo test -p umbra-eclipse --test path_limits`（27 本・実 2024 SLOW の中心線 194 点全点内包を含む）が、
ユニオンの実挙動を縛る。
