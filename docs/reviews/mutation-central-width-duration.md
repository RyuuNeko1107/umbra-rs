# mutation レビュー: 中心食 帯幅/継続・純関数（`central_width_and_duration` / `wrap_to_pi` / `great_circle_distance_km`）

対象: ISSUE-045 M9.6（中心食の帯幅 `path_width` と中心線継続 `central_duration`）の核心関数
`central_width_and_duration`（`crates/umbra-eclipse/src/global.rs`）と、それが使う 2 純関数
`wrap_to_pi`（同 global.rs・μ' 数値差分の ±2π 折返し）・`great_circle_distance_km`
（`crates/umbra-eclipse/src/axis_intercept.rs`・haversine 帯幅換算）。

## 実行
```
cargo mutants -p umbra-eclipse \
  --re 'great_circle_distance_km|central_width_and_duration|wrap_to_pi' --no-shuffle \
  -- --lib -- greatest_ wrap great_circle
```
（核心は `solve_greatest_eclipse` 経由の `greatest_*` 中心食オラクル＋純関数の高速 in-module
ユニットテスト〔`wrap_*` / `great_circle_*`〕が killer。テストフィルタに `wrap`/`great_circle` を
含め新ユニットテストを走らせる。）

## 結果（2026-06-21・工程7 再走）
**93 mutants: 90 caught・3 missed・4 unviable・0 timeout。**

### 生存 3 件
| # | 変異 | 区分 | 理由 |
|---|---|---|---|
| 1 | `global.rs:249 > → >=`（`rel_speed > 0.0`） | 等価（許容） | `rel_speed = rel_x.hypot(rel_y)` が**厳密 0.0** になるのは rel_x・rel_y がともに厳密 0 のときのみ。実中心食では rel≈0.2 Re/h で `> 0` と `>= 0` は同じ分岐を取る＝到達不能。 |
| 2 | `global.rs:249 && → ||`（`rel_speed > 0.0 && rel_speed.is_finite()`） | 等価（許容） | `&&`→`||` が結果を変えるのは「`rel_speed==0.0`（有限）」または「`rel_speed=+∞`」のときのみ。前者は #1 と同じく到達不能、後者は rel 成分が無限大＝瞬時ベッセル要素が非有限となる構成で実中心食・合成中心食いずれでも生成不能。よって全到達入力で `&&` と `||` は同値。 |

注: `central_width_and_duration` 関数内に `>` は 249 行の 1 箇所のみ、`&&` も 249 行の 1 箇所のみ
（他比較は別関数 `solve_limit_edge` に属す）。よって除外パターンは当該等価変異だけを的確に捕捉する。

### 当初 missed だった `great_circle_distance_km` の積項を追撃して caught 化
工程7 初回走で `axis_intercept.rs` の haversine `h` 式中 `lat1.cos() * lat2.cos()` の `* → /`
（`cos φ1 / cos φ2`）が missed だった。既存の純関数テストが
**赤道上（cos φ=1 で `*`=`/`）・同経度（dlon=0 で当該項が ×0）** のケースしか持たず、この積が
値を持つ斜め 2 点を欠いていたため。緯度・経度がともに非自明に異なる
`great_circle_distance_diagonal_pair_known_value`（10°N,0° ↔ 50°N,40°E ⇒ 5763.650 km・
haversine 既知式の手計算）を追加して caught 化（`/` 化で 7332.8 km へずれ tol 0.5 km を大きく外す）。

## 退化していたオラクルの非退化化（mutation 強化の本体）
工程7 で `central_width_and_duration` の rel_x/rel_y・数値中心差分の算術変異
（240/245/246 行系）が多数 missed だった根因は、独立オラクル
`greatest_central_duration_matches_two_l2p_over_rel_speed` の合成源 `CentralDurationSource` が
**退化**していたこと:
- `declination = 0`（sin d = 0）⇒ rel_x の `η·sin d` 項・rel_y の `μ'·ξ·sin d` 全項が ×0 で消える。
- `y = y0`（定数・y' = 0）⇒ rel_y の `vy` が 0・y の数値中心差分が 0。

これらが当該行の算術変異を**等価化**していた。`CentralDurationSource` に `y1_per_hour`（y を一次・
y'≠0）と `declination`（非零定数 D=0.2）を追加し、オラクルの期待 vx/vy/μ' を**ソースの既知傾き**
（X1=0.40・Y1=0.02・MU1=0.26 per hour）から組んで実装の数値微分とは独立に
`rel=(vx−μ'(ζcosd−η sind), vy−μ'ξsind)`・`duration=2|L2'|/|rel|·3600` を tol 1e-2s で照合。
gamma≪1 を保つ小傾きで中心食（軸が地表貫通）を維持。これで 240（微分スケール）・
245（rel_x 全項）・246（rel_y 全項）の変異がすべて caught に転じた。

## CI 除外提案（mutation.yml）
等価 2 件のみを的確に除外（killable な算術・分岐は除外しない）:
```
--exclude-re 'with >= in.*central_width_and_duration'
--exclude-re 'with || in.*central_width_and_duration'
```
実回帰ガード: 通常 CI の `cargo test -p umbra-eclipse`（合成中心食の独立オラクル＋純関数の
既知値ユニットテスト＋実 2017/2023/2024 の NASA ballpark）が帯幅・継続・折返し・大圏距離の
正しさを縛る。
