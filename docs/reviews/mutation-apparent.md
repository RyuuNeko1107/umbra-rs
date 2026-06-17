# ミューテーション生存変異の許容判断 — `umbra-ephemeris::apparent`（ISSUE-015 S2/S3/S4）

`cargo mutants --package umbra-ephemeris --file crates/umbra-ephemeris/src/apparent.rs`
（cargo-mutants 27.1.0, Docker 内）。**80 mutants: 64 caught / 12 unviable / 4 missed**。

対象:
- S2 `light_time_correct`（`sun_/moon_light_time_corrected_gcrs` 本体）= 光行時間補正。
  出力 `s = B_geo(t−τ) + (E(t−τ)−E(t))`、一次近似 `−v_E·τ`、SOFA `iauAtciq` の light-time ステップ。
- S3 `apply_iau_ab` / `aberrated_gcrs`（`sun_/moon_aberrated_gcrs` 本体）= 恒星光行差（SOFA `iauAb` 逐語）。
- S4 `sun_/moon_apparent_cirs` = 歳差章動 GCRS→CIRS（`gcrs_to_cirs_matrix·aberrated_gcrs`）。
  合成（行列適用先・順序・transpose）は end-to-end erfa オラクル＋合成同一性＋回転適用テストが捕捉、
  **新規生存なし**。

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

---

## ISSUE-043 S2: ジェネリック `apparent_cirs<E: Ephemeris>`

`cargo mutants --package umbra-ephemeris --file crates/umbra-ephemeris/src/apparent.rs
--re 'apparent_cirs|earth_velocity_gcrs|geocentric_gcrs_of|AstrometryOptions'`
（cargo-mutants 27.1.0, Docker 内）。**30 mutants: 13 caught / 8 unviable / 9 missed**。

`apparent_cirs` は具象チェーン（`light_time_correct`/`aberrated_gcrs`）と**同一式**を任意 `Ephemeris`
で回す（幾何位置・地球速度=−地心太陽速度・太陽距離を `eph.state` から導出）。回帰ブリッジテスト
（`generic_*_standard_matches_concrete_*` = AnalyticalEphemeris + standard が具象とビット級一致）と
Mock 幾何経路・全フラグ組合せ・エラー透過が物理/契約コアを 13 caught で捕捉。生存 9 件は**具象チェーンと
同一カテゴリの等価/許容下限変異**:

| 変異 | 効果 | なぜ許容か |
|---|---|---|
| 265:33 `/`→`%`/`*`（τ0 = `\|g0\|/c`） | 光行時間反復の初期値 | 不動点反復は縮小率 ≈ v_E/c ≈ 1e-4/反復で強収縮、τ0 によらず 5 反復内に同一 τ へ収束（巨大 τ0 でも収束）。戻り値不変＝真の等価変異（上記 `light_time_correct` 83 行と同種） |
| 272:48 `<`→`==`/`<=`（収束判定） | 早期 break の打切り条件 | 5 反復で完全収束、戻り値同一（早期 break 無効化のみ） |
| 272:35 `−`→`+`/`/`（`(next−tau).abs()`） | 収束判定式 | `next+tau`/`next/tau` は < 1e-6 不成立 → 全反復化のみ、戻り値不変 |
| 293:35 `/`→`%`/`*`（s_au = `\|sun_pos\|/AU`） | `iauAb` の w2 = SRS/s_au 微項 | w2 は相対論補正 ~0.004″、単位ベクトル寄与 ~2e-12（実 v 微小）。回帰ブリッジ tol（太陽 1 km ≈ 1.4e-3″）下。`apply_iau_ab` 自体は amplified テストで w2 を捕捉済み |
| 294:24 `−`→`+`（bm1 = `√(1−\|v\|²)`） | 光行差の bm1 係数 | 実 \|v\|≈1e-4 で 1∓1e-8、単位ベクトルへの寄与が tol 下。`apply_iau_ab` の amplified テストが bm1 経路を別途捕捉 |

## 判断
9 件すべて、(a) 光行時間の不動点反復収束簿記（τ0・収束判定）＝具象 `light_time_correct` と同種の真の
等価変異、(b) 恒星光行差の w2/bm1 微項＝実 v が微小で単位ベクトル寄与が許容下限・物理予算（0.10″）下、
のいずれか。物理・契約コア（幾何位置取得・地球速度符号・光行時間第2項・CIRS 回転・velocity 要否・
エラー透過）は 13 caught で捕捉、回帰ブリッジが具象との一致を担保。よって**全件許容**。
focused `--re` 実行＋本文書での列挙を正本とする（`light_time_correct` 同様、cargo-mutants の説明文が
行番号を含まず regex で等価分のみ選択除外できないため）。besselian.rs が S3 でジェネリック経路へ
切替後は、下流の高精度テスト（2017 gamma 等）が収束簿記変異の一部を捕捉する見込み。
