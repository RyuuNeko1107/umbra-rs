# mutation レビュー: rise/set limb 点（`cone_terminator_intersections` / `scan_periodic_sign_change_roots`）

対象: ISSUE-045 M9 残(3) サブスライス 3b（rise/set limb 点＝錐縁 ∩ WGS84 terminator 楕円）。
`crates/umbra-eclipse/src/axis_intercept.rs` の `cone_terminator_intersections`（k・残差・媒介・測地座標化の
数値核）と、そこから分離した resolution 機構 `scan_periodic_sign_change_roots`（周期残差の符号反転走査＋Brent）。

## 実行
```
cargo mutants -p umbra-eclipse \
  --re 'cone_terminator_intersections|scan_periodic_sign_change_roots' --no-shuffle \
  --exclude-re 'in scan_periodic_sign_change_roots' \
  -- --lib -- two_intersections d_zero cone_not_reaching cone_radius deterministic terminator_ellipse
```
killer は axis_intercept.rs の FAST 単体テスト（前方射影往復で ζ≈0・面内距離=cone_l、d=π/2 二円閉形式一致、
空 Vec、半径 load-bearing、決定性、扁平 k>1）。

## 結果（2026-06-21）
**40 mutants: 39 caught・1 unviable・0 missed**（`scan_periodic_sign_change_roots` 内は除外）。
- caught 39: `k = sin²d + cos²d/(1−f)²` の各項、媒介 `ξ=cos θ, η=sin θ/√k`、円残差 `(ξ−x)²+(η−y)²−cone_l²`、
  中心 (x,y)、`fundamental_to_geodetic` への ζ=0 引数等の算術を全捕捉。
- unviable 1: 戻り値を `Default` 等へ置換する型不成立。

## 機構分離の経緯（初回 2 missed → 分離で解消）
初回走（分離前・45 mutants）で `prev_r * r < 0.0`（符号反転判定）の **2 件が missed**:
| 変異 | 等価の理由 |
|---|---|
| `* → /`（`prev_r * r` → `prev_r / r`） | 積と商は**符号判定が同値**（sign(a·b)=sign(a/b), b≠0）。r=0（グリッドが根に一致）は測度ゼロ。よって同じ根をブラケットし真根を変えない。 |
| `< → <=`（`< 0.0` → `<= 0.0`） | 厳密ゼロ（`prev_r·r==0`＝グリッド点が根）でのみ差。浮動小数では測度ゼロで到達不能。 |

両者は `descending_sign_change_bracket` / `solve_zero_in_window` / `scan_point_count` と**同カテゴリ＝解像度・
符号ブラケット機構の等価変異**（真根は Brent が決める）。ただし `cone_terminator_intersections` には殺せる
算術（残差・k・媒介）も混在するため、関数レベルの粗除外（`with / in.*cone_terminator_intersections` は
残差の `*→/` も巻き込む）は不可。そこで**機構を `scan_periodic_sign_change_roots` へ分離**し、
`descending_sign_change_bracket` と同様に `--exclude-re 'in scan_periodic_sign_change_roots'` で wholesale 除外
（数値核は除外対象外＝全 caught を維持）。これで数値核 0 missed・機構は精度に影響しない解像度要素として除外、
という清明な境界になる。

実回帰ガード: 通常 CI の `cargo test -p umbra-eclipse`（FAST 前方射影往復・d=π/2 閉形式・空 Vec・半径
load-bearing・扁平 k>1）が rise/set 交点の正しさを縛る。

## 注記（後続 (3c)）
接する端（P1/P4 の 1 点接触）は符号反転が無く拾わない（交点数 0↔2 の遷移点）。外周端の扱いは (3c) 外周組立の
責務（§11.3/11.4）。
