# ミューテーション生存変異の許容判断 — `umbra-ephemeris::apparent`（ISSUE-015 S2/S3）

`cargo mutants --package umbra-ephemeris --file crates/umbra-ephemeris/src/apparent.rs`
（cargo-mutants 27.1.0, Docker 内）。**78 mutants: 64 caught / 10 unviable / 4 missed**。

対象:
- S2 `light_time_correct`（`sun_/moon_light_time_corrected_gcrs` 本体）= 光行時間補正。
  出力 `s = B_geo(t−τ) + (E(t−τ)−E(t))`、一次近似 `−v_E·τ`、SOFA `iauAtciq` の light-time ステップ。
- S3 `apply_iau_ab` / `aberrated_gcrs`（`sun_/moon_aberrated_gcrs` 本体）= 恒星光行差（SOFA `iauAb` 逐語）。

## caught（物理・契約のコア）

- **S2 第2項の符号・合成**（`B_geo(t−τ) + v_E.scale(−τ)`）、放射時刻 `add_days(−τ/86400)`、`next=|s|/c`。
- **S3 `iauAb` の主項**: `pnat[i]·bm1`、`w1 = 1 + pdv/(1+bm1)`、`w1·v[i]`、距離スケール `.scale(dist)`。
  恒星光行差 ~20.5″ を変える演算は全 caught。erfa.ab 厳密一致テスト（tol 1e-12）と距離不変・apex 向き・
  20.5″角テストが捕捉。
- **S3 w2 微項（SRS/s 相対論補正）の内部演算**: `w2 = SRS/s_au` の除算、`w2·(v−pdv·pnat)` の積・
  `(v−pdv·pnat)` の差・符号。実エポックでは w2 項の単位ベクトル寄与が ~1e-12 で tol に埋もれ生存するため、
  **純関数 `apply_iau_ab` を増幅速度（|v|≈0.06, s=0.5）で erfa.ab と直接突合するテスト**
  （`apply_iau_ab_matches_erfa_ab_amplified`）を追加し、w2 項を ~1e-9 に励起して 7 件すべて caught
  （s=0.5 で `SRS/s` と `SRS*s`/`SRS%s` も区別）。
- **S2 τ0 初期値 `|B_geo(t)|/c` の `/` と収束判定 `<`→`>`（1反復化）**: S2 単独のゆるいテストでは生存したが、
  S3 の **erfa.ab 厳密一致テスト（tol 1e-12）が pnat=unit(s_S2) 経由で S2 出力の最終ビットに依存する**ため、
  これらが招く ~1e-9〜1e-11 km の差を検出して caught（S3 追加の副次的な締め）。

## 生存 4 件（S2 不動点反復の収束打切り簿記・等価変異・許容）

`light_time_correct` の収束判定 `(next − tau).abs() < 1e-6`（行 83）に対する変異で、**いずれも「常に上限5反復まで回す」結果になり、5反復で完全収束するため戻り値がビット同一**:

| 変異 | 効果 | なぜ等価か |
|---|---|---|
| 83:44 `<`→`==` | float 厳密一致はほぼ不成立 → 早期 break せず常に5反復 | 5反復で完全収束（2反復で実質収束）、戻り値同一 |
| 83:44 `<`→`<=` | 境界 1e-6 ちょうどは到達不能 → 同上 | 同上 |
| 83:31 `−`→`+` | `next+tau`≈998s(太陽)/2.6s(月) は < 1e-6 不成立 → 同上 | 同上 |
| 83:31 `−`→`/` | `next/tau`≈1.0 は < 1e-6 不成立 → 同上 | 同上 |

収束判定は早期打切りの**最適化**にすぎず、不動点反復は τ0 によらず同じ τ へ強収縮（縮小率 ≈ v_E/c ≈ 1e-4/反復、
2反復で収束）。これら4変異は早期 break を無効化して全反復を回すだけで、収束値（位置・τ）を変えない。
**真の等価変異**で、係数符号（遠方エポックで励起可）と異なり原理的に kill 不可。numerical-policy §A3 が
「固定回数でなく相対収束判定」を要求するため固定反復化による kill も不可（仕様逸脱）。solver.rs の
Brent/golden 加速ヒューリスティクス（二分法フォールバックで根は契約どおり）と同 category。

> 注: 反転 `<`→`>`（1反復化）は S3 追加前は生存だったが、上記のとおり S3 の 1e-12 厳密一致テストが
> 1反復解と完全収束解の ~1e-11 差を検出するため、現在は caught。残るのは「より多く反復させる」4変異のみ。

## 判断
4 件すべて、不動点反復の収束打切り簿記に対する等価変異で、収束した戻り値を変えない。物理・契約のコア
（S2 第2項・S3 iauAb 全項・w2 微項含む）は 64 caught で捕捉済み。よって**全件許容**。

## ゲート注記
`light_time_correct` 内の生存（行 83 `<`/`−`）は、同関数内の捕捉済み演算と cargo-mutants の
レンダリング説明文が同一（行番号を含まない）のため `--exclude-re` で等価分のみを選択除外できない。
moon.rs（docs/reviews/mutation-moon.md）と同様、focused `--file` 実行＋本文書での列挙を許容の正本とする。
