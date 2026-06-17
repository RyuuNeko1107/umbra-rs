# Mutation review — `umbra-eclipse::axis_intercept`（ISSUE-043 S6a-i 逆射影）

`cargo mutants` の生存変異の列挙と許容判断。退行ガード `mutation.yml` は eclipse を含む。

## 結果

- `shadow_axis_surface_point`（影軸の地表貫通点＝前方射影の逆）の**実数学**
  （逆回転 `px=ζcosd−ysind / pz=ζsind+ycosd`・海面子午線楕円拘束 `ρcos²+(ρsin/(1−f))²−1`・
  局地時角 `H=atan2(x,px)`・東経 `λ=H−μ`・測地緯度 `φ=atan2(ρsin, ρcos·(1−f)²)`）は**全 caught**。
  前方射影 `project_observer_to_fundamental`（ISSUE-024）＋ `observer_geocentric`（ISSUE-010/011）を
  独立オラクルとする往復一致テストが、符号・係数・(1−f)² 因子・太陽側選択（ζ>0）を実効的に縛る。
- 生存（許容）= **`descending_sign_change_bracket`（粗走査ブラケット）内の算術変異のみ**。
  `conjunction::solve_zero_in_window` / `local_maximum::scan_point_count` と同じ「粗走査機構」カテゴリで、
  `mutation.yml` の `--exclude-re 'in descending_sign_change_bracket'` で除外する。

## 生存（許容）= 粗走査機構の等価変異（`descending_sign_change_bracket`）

初回試走（抽出前・本体内インライン）で 2 件が missed:

1. `step = ZETA_SCAN_MAX / ZETA_SCAN_STEPS`（刻み計算）の `/ → %`。
2. `r_lo * r_hi <= 0.0`（符号反転判定）の `* → /`。

いずれも**等価変異**である（独立に検証）:

- **`* → /`（符号判定）**: `a·b ≤ 0 ⟺ a/b ≤ 0`（到達しうる `r_hi ≠ 0` で恒真。`r_hi` は
  `residual(ZETA_SCAN_MAX) > 0` から走査点値しか取らず厳密 0 にならない）。乗算/除算で符号判定は不変。
- **`/ → %`（刻み）**: 物理的に妥当な中心食（軸が地表に当たる・ζ>0）では、残差
  `r(ζ)=Aζ²+Bζ+C`（`A=cos²d+sin²d/(1−f)²>0`）は **区間 `[0, ZETA_SCAN_MAX]` に太陽側の単一根 ζ₊** を持ち
  （負根 ζ₋<0 は区間外）。`C=r(0)>0` は gamma≳1 を要すが、そこでは扁平由来の `B`（係数 ≈0.0067）が
  小さすぎて頂点を負へ押せず**実根なし**＝軸が地表を外す。よって「2 正根」領域は**到達不能**。
  単一根なら刻み・点配置によらず Brent が同じ ζ₊ へ収束するため、刻み計算は結果に影響しない。
  20,240 ケース（φ∈[−88°,88°], λ, d∈[−50°,50°], μ）の精走査 vs 粗走査(刻み=1.05)スイープで
  **発散 0 件**を確認済（テスト設計時の独立解析）。

抽出により上記 2 件は `descending_sign_change_bracket` 内に局在化する。加えて cargo-mutants は
**関数全体置換**（`descending_sign_change_bracket` を定数 `Some((0.0, 1.0))` / `Some((1.0, 0.0))` で
置換）も missed にする — 定数ブラケット `[0, 1]` は単一根 ζ₊（∈(0,1]）を含むため Brent が同じ ζ₊ を
精解し、軸が外れる場合は `residual` が同符号で Brent 自身が `RootNotBracketed` を返すため、いずれも
end-to-end の挙動が不変だからである（ブラケットの「具体値」ではなく「根を括る」ことだけが load-bearing で、
それは Brent が再検証する）。これらも等価。よって関数名そのもの `--exclude-re 'descending_sign_change_bracket'`
で（演算子変異・関数全体置換の両形を）除外する（`solver.rs` / `solve_zero_in_window` / `scan_point_count`
除外と同方針）。逆射影の実契約
（往復一致・太陽側選択・半球符号・経度正規化・軸が外す→RootNotBracketed・実 2017 大局妥当性）は
通常 CI の `cargo test`（`axis_intercept::tests` 群）と本体の全 caught 変異が担保する。
