# §3 天体暦・見かけ位置（VSOP87D / ELP-MPP02 + 補正チェーン）

> 正本: `algorithms.md §0`（記号表）・`conventions.md §5.1`（補正順序 = SOFA iauAtciq）・`§5.2`（μ/CIO）・`§9`（半径）。本セクションはこれらに**厳密準拠**する。
> 状態: ドラフト（Milestone 0）。一次資料未確認の式番号・符号は「**要確認**」を残す。
> 関連 Issue: ISSUE-012（trait）, 013（太陽 VSOP87D）, 014（月 ELP/MPP02）, 015（見かけ位置）, 033/034（係数生成）。

---

## 目的と入出力（単位・時刻系・座標系・フレーム）

太陽・月の**見かけ地心位置** `r_sun`, `r_moon`（algorithms.md §0）を、解析暦の幾何位置から補正チェーンを経て供給する。

| 項目 | 内容 |
|---|---|
| 入力 | 観測時刻 `time_tt: TtInstant`、`AstrometryOptions`（light_time / aberration / precession_nutation / relativistic_deflection） |
| 出力 | `Position<F>`（F = `Gcrs` の見かけ位置、または §2 連鎖で `Cirs` まで）、見かけ距離 `R_sun`/`R_moon`（km）、視半径 `s_sun`/`s_moon`（rad） |
| 単位 | 位置 = km, 速度 = km/s, 角度 = rad, c = km/s（conventions §1/§4.1）。暦内部 AU 許容 → trait 境界で km（× 149597870.7） |
| 時刻系 | 入力 `TtInstant`。暦級数評価は **TDB**（TT≈TDB, 差 ≲2ms → 太陽 ≲0.001″/月 ≲0.03″ 許容, metadata 記録）。光行時間補正後の放射時刻も TDB で再評価（ISSUE-013/014/015, conventions §6） |
| 座標系/フレーム | VSOP87D/ELP native = 黄道・平均分点 of date（球面）→ 直交。GCRS 内で補正 → §2 で CIRS（conventions §5/§5.1） |

**確定方針（厳守, conventions §5.1 / D3）**:
- 見かけ位置補正は **SOFA `iauAtciq` の順序に固定**: light-time → **GCRS 内で deflection → aberration** → その後 **frame bias + IAU2006 歳差 + IAU2000A 章動（§2）で CIRS**。
- **「章動後に aberration」を禁止**（順序逆転は数″誤差, ISSUE-015）。
- aberration の観測者速度 = 地球重心速度を **VSOP87 級数の項別解析微分**で取得（numerical-policy §A2(1)）。

---

## 記号（algorithms.md §0 参照。固有補助記号のみ）

| 記号 | 意味 | 単位 |
|---|---|---|
| `T` | J2000 TDB からの**ユリウス千年** = `(JD_TDB − 2451545.0)/365250`（VSOP87 引数, §0 millennia） | 無次元 |
| `t_c` | J2000 TDB からのユリウス**世紀** = `(JD_TDB − 2451545.0)/36525`（ELP/MPP02 引数。採用版に従う, ISSUE-014/034） | 無次元 |
| `L, B, R` | VSOP87D の地球日心黄経・黄緯・動径（of date, 球面） | rad, rad, AU |
| `V, U, r` | ELP/MPP02 の月地心黄経・黄緯・距離（of date） | rad, rad, km |
| `A_k, B_k, C_k` | 級数の振幅・位相(rad)・振動数(rad/千年 or rad/世紀) | — |
| `D, l, l′, F, Ω` | Delaunay 基本引数（月・摂動引数） | rad |
| `τ` | 光行時間（light-time） | s |
| `r_⊕`, `v_⊕` | 地球重心の日心位置・速度（aberration 用） | km, km/s |
| `v̂`, `β`, `γ` | 観測者速度方向、`β=|v|/c`、`γ=1/√(1−β²)`（iauAb 用） | — |
| `s_sun, s_moon` | 太陽・月の見かけ視半径（§0, conventions §9） | rad |

---

## 数式（番号付き・各式に出典）

### 3.1 VSOP87D 太陽暦（地球日心 → 太陽地心）

**(E1) VSOP87 級数**
出典: Bretagnon & Francou (1988), A&A 202, 309。VSOP87**D** = 黄道・平均分点 of date・球面（heliocentric）。

各座標変数 `s ∈ {L, B, R}` を、べき α=0..5 ごとの周期項総和で評価する:

```
s = Σ_{α=0..5}  T^α · Σ_k  A_k · cos( B_k + C_k · T )
```

- `T` = ユリウス千年（TDB）。`A_k`（L,B 無次元 / R は AU）, `B_k` (rad), `C_k` (rad/千年) は packed 係数から（ISSUE-033, T 単位一致）。
- 級数和は **振幅昇順ナイブ和**（numerical-policy §A1, Kahan 不要）。位相 `B_k + C_k·T` は **mod 2π 縮約しない**（numerical-policy §A1）。

**(E2) 球面 → 日心黄道直交（of date）**
出典: 標準球面→直交変換。

```
r_earth_helio = R · ( cos B · cos L,  cos B · sin L,  sin B )      (黄道 of date 直交, AU)
```

**(E3) 太陽地心 = 地球日心の反転**
出典: ISSUE-013（太陽地心方向 = −地球日心方向）。

```
r_sun_geo(geom) = − r_earth_helio
```

`frame=Icrs/CIRS` 要求時は §2（黄道 of date → ICRS は frame bias + 黄道傾斜, `iauObl06`/§2）へ委譲。

### 3.2 ELP/MPP02 月暦（地心黄道 of date）

**(E4) 基本引数（Delaunay 等）**
出典: Chapront & Francou (2002), A&A 387, 700（MPP02 DE fit 版）。基本引数 `D, l, l′, F`（および平均経度）の時間多項式は MPP02 採用版（ISSUE-034 が原データから抽出, 時間単位 `t_c` または `T` を packed メタと一致）。

**(E5) 主問題 + 摂動級数**
出典: Chapront & Francou (2002)。経度 `V`・緯度 `U`・距離 `r` を、主問題（Delaunay 引数の整数結合に対する正弦/余弦項）＋摂動級数（惑星摂動・地球扁平・潮汐・相対論・固有摂動）で評価:

```
V = V_mean + Σ_k  A_k · sin( Σ (整数係数)·(基本引数) )
U =          Σ_k  A_k · sin( Σ (整数係数)·(基本引数) )
r =          Σ_k  A_k · cos( Σ (整数係数)·(基本引数) )
```

- 振幅・整数引数係数は packed から（ISSUE-034）。級数和は §A1 同様（昇順・位相非縮約）。
- **要確認**: 主問題級数の sin/cos 区分・平均経度の加算位置は MPP02 採用版に厳密準拠（系統誤差化を防ぐ, ISSUE-014）。

**(E6) 球面 → 地心黄道直交（of date）**
```
r_moon_geo(geom) = r · ( cos U · cos V,  cos U · sin V,  sin U )      (黄道 of date 直交, km)
```

### 3.3 見かけ地心位置の補正チェーン（SOFA iauAtciq 順, conventions §5.1）

> **適用順序（D3 固定）**: ① light-time → ② GCRS 内 deflection → ③ GCRS 内 aberration → ④ frame bias + IAU2006 歳差 + IAU2000A 章動（§2, GCRS→CIRS）。

**(E7) 光行時間反復（light-time）**
出典: ISSUE-015 / numerical-policy §A3。放射時刻 `t − τ` を反復で解く（地心見かけ, 観測者 = 地心）:

```
τ_0 = |r_target(t)| / c
τ_{k+1} = |r_target(t − τ_k)| / c
収束: |τ_{k+1} − τ_k| < 1e-6 s   （上限 5 反復, 通常 2–3）
```

- 後退時刻 `t − τ` で暦（E1–E6）を**再評価**して幾何方向ベクトルを確定（速度外挿は使わない, numerical-policy §A3）。
- `c = 299792.458 km/s`（定義値, conventions §4.1）。月 τ≈1.28 s, 太陽 τ≈499 s。
- `τ` の加算は `JulianDate2` の part2 に足して再正規化（numerical-policy §A1/§A3）。

**(E8) 重力偏向（deflection, Standard は太陽のみ・既定 OFF / M1 扱い）**
出典: SOFA `iauLd`（`iauLdsun` 経由）。Klioner 系の点質量偏向式。

```
e = (u·q − (e_SO·u)(e_SO·q)) を用い、
u' = u + (2·G·M_⊙ / (c² · d_SO)) · [ q·(e_SO·u) − e_SO·(u·q) ] / (1 + e_SO·q)
```

- `u` = 観測者→天体 の単位方向（GCRS）, `q` = 太陽→天体 方向, `e_SO` = 太陽→観測者 方向, `d_SO` = 太陽-観測者距離。SOFA 定数 `SRS = 2·G·M_⊙/c²`（AU 単位の Schwarzschild 半径相当, IAU 公称 GM_⊙ から）。
- 偏向は太陽中心との角離隔が小さいとき SOFA `iauLd` の閾値で抑制（太陽縁内で 0 に落とす）。**要確認**: SOFA `iauLd` の正確な定数 SRS と抑制閾値（一次は SOFA `ld.c`）。
- **Standard は既定 OFF（M1）**。日食=太陽縁近傍のため寄与上限を一度実測（ISSUE-015 M1 受入テスト, ≪0.05″ 想定）、省略を metadata に記録（必須, ISSUE-015）。

**(E9) 年周光行差（aberration, 相対論的）**
出典: SOFA `iauAb`。式（ISSUE-015 §38 に記載の SOFA `iauAb` 形）:

```
u' = ( u/β_r + (1 + (u·v_c)/(1 + 1/γ)) · v_c ) / (1 + u·v_c)
```

ここで `v_c = v_⊕/c`（地球重心速度を c で無次元化）, `β_r = 1/γ = √(1 − |v_c|²)`, `u` = 補正前単位方向。**要確認**: SOFA `iauAb` の引数規約（`bm1 = √(1−|v|²)`, `bdv = u·v`）を一次（SOFA `ab.c`）で最終確認。

- **観測者速度 = 地球重心の日心速度 `v_⊕`**。VSOP87 級数の**項別解析微分**で取得（numerical-policy §A2(1)）:

**(E10) 速度の解析微分（級数項別）**
出典: numerical-policy §A2(1)。
```
d/dT [ A_k · cos(B_k + C_k·T) ] = − A_k · C_k · sin(B_k + C_k·T)
```
- 解析微分は厳密・ほぼ無コスト。必要精度（光行差 0.05″ 確保に δv < 73 m/s, §A2(1)）を桁で上回る。月（ELP）も同様に項別微分で速度供給。

**(E11) frame bias + 歳差 + 章動（最後に適用）**
出典: §2 (F1)（`Q(t)` = GCRS→CIRS, IAU2006/2000A, CIO ベース）。
```
r_CIRS = Q(t) · r_GCRS_apparent
```
偏向・光行差は GCRS（J2000 軸）で先に適用済み。bias/歳差/章動は**最後にまとめて**回す（conventions §5.1 / D3）。

### 3.4 見かけ視半径（conventions §9）

**(E12) 太陽視半径**
出典: conventions §9（IAU2015 公称 R_sun = 696000 km）。
```
s_sun = asin( R_sun_phys / R_sun )   ≈ R_sun_phys / R_sun      (R_sun_phys = 696000 km)
```

**(E13) 月視半径**
出典: conventions §9（月半径係数 k, Re = WGS84 a）。
```
s_moon = asin( k · Re / R_moon )   ≈ k · Re / R_moon
```
- `k`: `IauMean`=0.2725076（既定）/ `EspenakUmbral`=0.272281 / `EspenakPenumbral`=0.2725076（conventions §9, LunarRadiusModel）。選択を metadata に記録。NASA 照合時は Espenak 慣習へ切替（系統差を accuracy.md）。
- `asin` 引数は `[-1,1]` クランプ（丸め誤差考慮, numerical-policy §A5）。

---

## 手順（実装順・数値注意・補正の適用順序）

1. **時刻引数生成**: `T`（千年, VSOP）・`t_c`（世紀, ELP）を `JulianDate2`（TDB）から生成。エポック 2451545.0 減算は整数部側で厳密（numerical-policy §A1）。TT→TDB は同一視＋差を metadata（ISSUE-012/013/014）。
2. **幾何地心位置**: 太陽 (E1→E3)、月 (E4→E6)。級数和は振幅昇順・位相非縮約（§A1）。係数は packed（手書き禁止, conventions §11 / ISSUE-033/034）。
3. **補正チェーン（D3 順, conventions §5.1）**:
   1. light-time (E7): 後退時刻で暦再評価 → 幾何方向確定。
   2. deflection (E8): opts ON 時のみ（既定 OFF, 太陽のみ）。OFF は metadata 記録。
   3. aberration (E9): 地球重心速度 (E10 解析微分) で `iauAb`。
   4. frame bias + 歳差 + 章動 (E11, §2 `Q(t)`): **最後に**適用。内部粗スキャン（非公開）は IAU2000B 許容（公開出力は 2000A）。
4. **視半径** (E12/E13): 見かけ距離から。k 選択を metadata。
5. **出力**: `Position<F>` と適用補正記録。Standard は light_time/aberration/precession_nutation を**強制 ON**（型/エンジンで固定, ISSUE-015）。

**数値注意（横断）**:
- light-time 反復は相対収束判定（1e-6 s, 上限 5）、固定回数にしない（numerical-policy §A3 / ISSUE-015）。
- aberration 速度は暦と同一バックエンドの地球速度（解析微分既定, 速度 None の暦は対称差分要求, ISSUE-015）。
- 補正順序入替（章動後 aberration 等）を**禁止**し回帰テストで保護（ISSUE-015 D3）。

---

## 境界・特異・異常系

- **TDB/TT 同一視**: 差 ≲2ms の位置影響を metadata（max_residual_arcsec 帰属外, accuracy §2.4）。
- **light-time 未収束**: 5 反復超で `SolverDidNotConverge`（実運用では起きない安全弁, numerical-policy §A3）。
- **deflection の太陽縁近傍**: SOFA `iauLd` 抑制閾値で 0 へ連続的に落とす（日食配置で寄与上限を M1 実測, ISSUE-015）。
- **asin/acos クランプ**: 視半径・方向の `[-1,1]` クランプ（numerical-policy §A5）。
- **L/V 正規化**: 暦層は `[0,2π)`、連続性が要る求解（合）は呼び出し側で連続化（conventions §2 / ISSUE-013）。
- **黄道版/時間単位の不一致**（J2000 vs of date, 世紀 vs 千年）: 系統誤差化するため packed メタと厳密一致＋ラウンドトリップ＋DE 差分で検出（ISSUE-014/034）。

---

## 検証（二段ゲート: M2 暫定 = Mock+SOFA+NASA / M10 = JPL DE 差分）

accuracy §3.1/§3.3、ISSUE-013/014/015 受入テスト準拠。基準値は **DE オラクルから動的取得・ハードコード禁止**（conventions §11）。GPL 実装の数値を期待値に貼らない（基準は DE, ISSUE-014/034）。

- **M2 暫定ゲート（ISSUE-047）**: DE 取込前は MockEphemeris（幾何足場）＋ SOFA 参照（検証のみ・移植禁止）＋ NASA 公開値（第二義・慣習差明記）で打切り次数を**暫定**確定（保証値化しない, accuracy §3.3）。
- **M10 最終ゲート（DE 差分, accuracy §3.1/§4 層分解）**:
  - 太陽（補正前幾何, E1–E3）: DE440 差分で **残差 0.05″ 級**（accuracy §2.4 / ISSUE-013）。
  - 月（補正前幾何, E4–E6）: DE440 差分で **残差 0.1″ 級**（accuracy §2.4 / ISSUE-014）。距離 r の δr/r を別途ガード（視半径→食分, さらに地平視差 π≈3422″·(δr/r) 経由で局地接触へ, accuracy §2.1L / ISSUE-014）。
  - 見かけ位置（補正後, E7–E11）: Standard 補正後の方向角差を 1900–2100 で測定、**月 ≲0.1″ / 太陽 ≲0.05″**（暦残差と補正残差を層分解, ISSUE-015）。
- **補正分解テスト**（ISSUE-015）: light_time のみ / +歳差章動 / +aberration を段階適用し、各段の寄与が既知オーダー（aberration ≈20″ 級, 太陽 ≈20.5″; 月 light-time ≈数″移動）に一致。
- **適用順序固定の回帰（D3）**: 正順（iauAtciq）と入替（章動後 aberration）で数″差が出ること、正順が DE 同等パイプラインと一致することを固定。順序が変わったら fail（ISSUE-015）。
- **deflection 省略上限（M1）**: 太陽縁近傍配置で ON/OFF 差を実測し ≪0.05″ を metadata 注記（ISSUE-015）。
- **速度**: 解析微分 vs 対称差分 vs DE 速度の一致（角速度許容内, ISSUE-013/014）。
- **反転テスト**: |地球日心方向 + 太陽地心方向| ≈ 0（ISSUE-013）。
- **フレーム整合**: EclipticOfDate↔Icrs ラウンドトリップ ≲1mas（§2 連携）。

---

## 出典

- Bretagnon & Francou (1988), A&A 202, 309（VSOP87, VSOP87D 黄道 of date 球面）。
- Chapront & Francou (2002), A&A 387, 700（ELP/MPP02, DE fit 版・主問題＋摂動）。
- SOFA `iauAb`（年周光行差, 相対論）, `iauLd`/`iauLdsun`（太陽重力偏向）, `iauAtciq`（補正順序の正本）, `iauPnm06a`/`iauC2i06a`（§2, frame bias+歳差+章動）, `iauObl06`（黄道傾斜）（参照のみ・移植しない, data-sources §0 / conventions §11）。
- IERS Conventions (2010), ch.5（歳差章動・フレーム, §2 経由）。
- conventions §5/§5.1/§5.2/§9, accuracy §2.1/§2.1L/§2.4/§3, numerical-policy §A1/§A2/§A3/§A5, data-sources §2.1/§2.2/§0。
- 関連 Issue: ISSUE-012, 013, 014, 015, 033, 034。
