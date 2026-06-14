# §5 影円錐幾何（Step5）

> 正本: `algorithms.md §0`（記号表）・`conventions.md §4`（距離 km・Re）・`§9`（半径モデル k）・`§5`（フレーム）・`accuracy.md §2.2`（acos/asin クランプ）・`numerical-policy.md §A5`（クランプ）。本セクションはこれらに**厳密準拠**する。
> 状態: ドラフト（Milestone 0）。一次資料で未確認の式番号・符号は「**要確認**」を残し、推測で式番号を書かない。
> 関連 Issue: ISSUE-019（影円錐幾何・金環判定）。後続: ISSUE-020（基底）, ISSUE-021（ベッセル要素 tan f1/tan f2, l1/l2 への写像）。

---

## 目的と入出力（単位・時刻系・座標系・フレーム）

太陽（有限サイズの光源）と月（遮蔽球）から、**半影（penumbra）・本影（umbra）・反本影（antumbra）** の円錐幾何を構成する（architecture §6 `ShadowCone`, ISSUE-019）。本セクションは**地心ベクトルでの純幾何**まで（基本面基底＝§6/ISSUE-020、Re 無次元の l1/l2/tan f への写像＝§6/ISSUE-021）。

| 項目 | 内容 |
|---|---|
| 入力 | 太陽・月の見かけ地心位置 `r_sun`, `r_moon`（km, §3 ISSUE-015）、太陽半径 `R_sun_phys` = 696000 km（conventions §9）、月半径 `r_moon_phys = k·Re`（conventions §9 の k 選択） |
| 出力 | `ShadowCone`: `axis_direction`（影の伸びる向き）, `axis_origin`（月中心）, `umbra_apex`/`penumbra_apex`（頂点）, `umbra_half_angle`(=f2 の素), `penumbra_half_angle`(=f1 の素)。金環判定 `UmbraApexLocation` |
| 単位 | 距離 = km（幾何計算, conventions §4）。半角 = rad。Re = WGS84 a（conventions §4.1） |
| 時刻系 | 入力は特定 TT の見かけ位置（呼出側が `TtInstant` 供給, conventions §6）。**本節は時刻非依存の純幾何**（ISSUE-019） |
| 座標系/フレーム | 地心（GCRS/見かけ, §3 出力）。FundamentalPlane への基底変換は §6/ISSUE-020 |

**確定方針（厳守, ISSUE-019）**: 本影＝**収束**錐（太陽より月が小さい）、半影＝**発散**錐。本影頂点の地球側/反地球側判定が**金環/皆既の境界**を決める（種別最終確定は §8/ISSUE-023）。

---

## 記号（algorithms.md §0 参照。本節固有の補助記号のみ定義）

§0「ベッセル要素」表（f1, f2, l1, l2）・「天体位置」表を正本とする。本節固有:

| 記号 | 意味 | 単位 |
|---|---|---|
| `R_sun_phys` | 太陽実半径 = 696000 km（IAU2015 公称, conventions §9） | km |
| `r_moon_phys` | 月実半径 = `k · Re`（k: conventions §9 の選択） | km |
| `D` | 太陽中心-月中心 間距離 = `|r_sun − r_moon|` | km |
| `û` | 月→太陽 方向の単位ベクトル = `(r_sun − r_moon)/D` | — |
| `axis_direction` | 影の伸びる向き = `−û`（月→太陽の反対, 地心→反太陽側） | — |
| `f1`, `f2` | 半影/本影 円錐の半頂角（§0。本節で `sin f`/`tan f` を算出） | rad |
| `L_umbra` | 月中心→本影頂点 の距離（収束点までの距離） | km |
| `L_pen` | 半影頂点（後方発散錐の見かけ頂点）の月中心からの距離 | km |
| `umbra_apex` | 本影頂点（収束点）の地心位置 | km |
| `penumbra_apex` | 半影頂点（発散錐の見かけ頂点）の地心位置 | km |

> 向き規約の最重要点（ISSUE-020/021 と整合）: ベッセル慣習の**影軸 z 方向は「地心→太陽」向き**（§6/ISSUE-020）。一方、本節の `axis_direction`（影の伸びる向き）は `−û = 太陽→月→地球` 側＝反太陽向き。両者は**逆向き**であることを明記し、§6 へ渡す際の符号を取り違えない（ISSUE-019 実装メモ「NASA のベッセル軸は太陽中心と月中心を結ぶ線」, ISSUE-020 レビュー重点）。

---

## 数式（番号付き・各式に出典）

出典（共通）: Explanatory Supplement to the Astronomical Almanac (3rd ed.), **Ch.11「Eclipses of the Sun and Moon」Besselian elements 節**（影錐の半角・頂点距離の定義）。補助: Meeus, *Astronomical Algorithms* (2nd ed.), **Ch.54「Eclipses」**（l1, l2, f1, f2 の構成）。**要確認**: Explanatory Supplement Ch.11 の影錐半角・頂点距離の正確な式番号（一次で確定）。

### 5.1 太陽-月軸と影軸方向

**(C1) 軸距離と方向**
```
D              = | r_sun − r_moon |
û              = ( r_sun − r_moon ) / D          # 月 → 太陽 単位方向
axis_direction = − û                              # 影の伸びる向き（反太陽側）
axis_origin    = r_moon                           # 影軸の基準点（月中心）
```

### 5.2 半影・本影の半頂角

**(C2) 半影の半頂角 f1（発散錐）**
出典: Explanatory Supplement Ch.11 / Meeus Ch.54（**要確認**: 式番号）。
```
sin f1 = ( R_sun_phys + r_moon_phys ) / D
```
- 半影は **外接共通接線**（太陽・月の外側で交差）で作る発散錐。`f1` は大きい方の半角。

**(C3) 本影の半頂角 f2（収束錐）**
出典: Explanatory Supplement Ch.11 / Meeus Ch.54（**要確認**: 式番号）。
```
sin f2 = ( R_sun_phys − r_moon_phys ) / D
```
- 本影は **内接共通接線**で作る収束錐。太陽 > 月（`R_sun_phys > r_moon_phys`）ゆえ `f2 > 0` の小さい値。`R_sun_phys ≤ r_moon_phys`（非物理）は縮退（後述）。
- §6/ISSUE-021 が保持するのは `tan f1 = tan(f1)`, `tan f2 = tan(f2)`（NASA 表記, §0）。**符号: f2（本影/収束）は正の小さい値、f1（半影/発散）は正の大きい値**（ISSUE-021）。

### 5.3 頂点位置

**(C4) 本影頂点（収束点）**
出典: ISSUE-019（`L = r_moon_phys / sin f2`）/ Explanatory Supplement Ch.11（**要確認**: 式番号）。
```
L_umbra    = r_moon_phys / sin f2
umbra_apex = r_moon + L_umbra · axis_direction        # 月中心から反太陽側へ L_umbra
```
- `umbra_apex` は本影錐が一点に収束する点。ここを越えて軸を延長した発散錐が**反本影（antumbra）**（C6）。

**(C5) 半影頂点（後方発散錐の見かけ頂点）**
出典: ISSUE-019（半影頂点も同様に算出）。
```
L_pen        = r_moon_phys / sin f1
penumbra_apex = r_moon − L_pen · axis_direction       # 半影錐は太陽側に見かけ頂点（発散方向の逆）
```
- 半影は発散錐のため見かけ頂点は **月の太陽側（`+û` 方向）** にある。符号を本影（収束＝反太陽側）と取り違えない（ISSUE-019 実装メモ）。**要確認**: 半影頂点を §6 で直接使うか（NASA 流は l1/l2 を基本面で評価し頂点は中間量）— §6/ISSUE-021 の写像と整合確認。

### 5.4 金環/皆既の判定（本影頂点 vs 地球面）

**(C6) 本影頂点の地球側/反地球側（金環判定の核）**
出典: Explanatory Supplement Ch.11（central eclipse の umbral cone reaching Earth 条件）/ NASA Espenak の annular 定義（ISSUE-019）。
```
本影頂点が地球面より手前（地球側） ⇒ 本影が地表に届かない ⇒ 反本影が地表 ⇒ 金環（BeforeEarthSurface）⇒ l2 > 0
本影頂点が地球面以遠（反地球側）   ⇒ 本影が地表に届く           ⇒ 皆既（OnOrBeyondEarth）⇒ l2 < 0
```
- 判定: 影軸が地球面を貫く点（軸と WGS84 楕円体/球の交点）を求め、その交点までの距離と `L_umbra`（C4）を比較する。`umbra_apex` の地心距離と、軸の地球貫通点での地球面までの距離を比較してもよい（ISSUE-019）。**要確認**: 地球面との交点を球近似（Re）で取るか楕円体で取るか — 金環/皆既の閾値精度に効くため §6/§8 と整合（conventions §4, ISSUE-019）。
- §0 の符号規約との接続: この判定は §6/ISSUE-021 の **`l2 < 0 ⇒ 皆既`（金環は l2 > 0）** へ写像される。すなわち本影頂点が基本面の地球側にある（金環）と `l2` が正符号で出る（§0, ISSUE-021）。本節は頂点の幾何的位置を供給し、`l2` の符号化は §6 の責務。

**(C7) 反本影（antumbra）**
出典: ISSUE-019。
```
反本影 = 本影錐を頂点 umbra_apex の先（反太陽側）へ延長した発散錐
```
- 金環食で地表に見える円環はこの反本影内。本影半角と反本影は**同一錐の表裏**であり、`umbra_half_angle`（=f2）1つで両方を表し、`umbra_apex` の前/後で本影/反本影領域を切り分ける（ISSUE-019 実装メモ）。

---

## 手順（実装順・数値注意・補正の適用順序）

1. **軸構成 (C1)**: `r_sun − r_moon` から `D`, `û`, `axis_direction = −û`, `axis_origin = r_moon`。
2. **半角 (C2/C3)**: `f1 = asin(clamp((R_sun_phys + r_moon_phys)/D, −1, 1))`, `f2 = asin(clamp((R_sun_phys − r_moon_phys)/D, −1, 1))`。クランプ必須（accuracy.md §2.2 / numerical-policy §A5）。
3. **頂点 (C4/C5)**: `sin f2 ≈ 0`（頂点が無限遠＝皆既/金環境界）で 0 除算回避（極限処理 or 大距離扱い）。`L_umbra = r_moon_phys / sin f2`。
4. **金環判定 (C6)**: 影軸の地球貫通点を求め頂点と比較 → `UmbraApexLocation`。
5. **縮退検出**: `D ≈ 0`、`R_sun_phys ≤ r_moon_phys`（非物理）、軸が地球を全く貫かない等で `DegenerateGeometry`（api-draft §3.5, ISSUE-019）。

**k の使い分け（conventions §9, ISSUE-019）**: NASA 照合時は **本影系=EspenakUmbral(0.272281)・半影系=EspenakPenumbral(0.2725076)** を config 経由で使い分ける。既定 `IauMean(0.2725076)` は単一値。どのモデルで計算したか metadata に残す（ISSUE-019）。

**数値注意（横断, numerical-policy）**:
- asin 引数を `[-1,1]` クランプ（丸め誤差込み, accuracy.md §2.2 / §A5）。
- `sin f2 ≈ 0` の 0 除算回避（皆既/金環境界＝頂点無限遠, ISSUE-019）。極端配置（月遠・本影が地球に届かず頂点距離が発散）でも破綻しないこと。
- 純幾何は **f64 機械精度（≪ バジェット）** で計算する（誤差は入力位置律速, ISSUE-019）。

---

## 境界・特異・異常系

- **`D ≈ 0`**: `DegenerateGeometry`（軸方向不定, ISSUE-019）。
- **`R_sun_phys ≤ r_moon_phys`（非物理）**: `sin f2 ≤ 0` → 本影錐が定義できない → `DegenerateGeometry`（ISSUE-019）。
- **`sin f2 ≈ 0`（皆既/金環境界＝頂点無限遠）**: 0 除算を避け、頂点を大距離 or 極限処理。判定 (C6) はハイブリッド境界として安定に（最終分類は §8/ISSUE-023, ISSUE-019）。
- **軸が地球を全く貫かない**: 影が地球を外す配置。(C6) は判定不能のため上流（§4 Step3 が既に粗棄却済み）か `DegenerateGeometry`。エラーで止めず可能性を潰さない設計余地を残す（ISSUE-019 実装メモ）。
- **asin クランプ境界**: 引数が ±1 近傍で丸め誤差により域外化 → クランプで安定（accuracy.md §2.2）。
- **半影頂点の符号**: 発散錐ゆえ太陽側（本影と逆）。符号取り違えは l1 写像（§6）を鏡像化（ISSUE-019）。

---

## 検証（基準値の出典。実装へ値コピー禁止 = conventions §11）

accuracy.md §3.1、ISSUE-019 受入テスト準拠。基準値は数式オラクル/DE/NASA から取得し**実装へハードコードしない**（conventions §11）。

- **半角の解析検証（L1）**: 既知の `R_sun_phys, r_moon_phys, D` で `sin f1/f2` を手計算（オラクル＝数式, 実装非コピー）と照合。**本影=収束・半影=発散** を符号で確認（ISSUE-019）。
- **MockEphemeris 人工ケース（L4, accuracy.md §3.1）**:
  - 明確な皆既（月近・本影頂点が地球面以遠）→ `OnOrBeyondEarth`。
  - 明確な金環（月遠・本影頂点が地球面手前）→ `BeforeEarthSurface`。
  - 境界（頂点が地球面ちょうど）→ ハイブリッド境界の判定安定性（最終分類は §8/ISSUE-023）。
- **DE/NASA 整合（第二義, data-sources §4.1）**: 既知の皆既/金環食で頂点判定が NASA 種別と一致（**k 慣習を Espenak に揃える**, accuracy.md §3.1）。
- **k 値感度**: `EspenakUmbral(0.272281)` vs `IauMean(0.2725076)` で本影半角・頂点が変わり皆既/金環境界がずれることを定量確認（conventions §9 系統差を accuracy.md に記録, ISSUE-019）。
- **縮退系**: `D≈0`、`R_sun_phys ≤ r_moon_phys` → `DegenerateGeometry`。asin クランプ境界（ISSUE-019）。

---

## 許容誤差

- accuracy.md §2.1 幾何バジェットの一部（**影幾何誤差**, §4 層分解）。半角・頂点の誤差は最終的に l1/l2（§6/ISSUE-021）経由で食分・継続時間に効く（食分 0.001 ≈ 1.9″, accuracy.md §2.2）。
- 半角誤差は**入力位置・距離の誤差**（§3 の月 0.1″/太陽 0.05″）で律速され、本節の純幾何は **f64 機械精度（≪ バジェット）** で計算する（誤差を層ごとに分解, accuracy.md §0/§4, ISSUE-019）。
- **金環/皆既境界（ハイブリッド）付近は k 選択で結果が変わる**（conventions §9）。系統差を accuracy.md へ記録し、テストでは k を固定して比較（誤差を隠さない, accuracy.md §0, ISSUE-019）。

---

## 出典

- Explanatory Supplement to the Astronomical Almanac (3rd ed.), Ch.11「Eclipses of the Sun and Moon」Besselian elements 節（影錐半角・頂点距離・central eclipse の umbral cone 条件）。**要確認**: 式番号。
- Meeus, *Astronomical Algorithms* (2nd ed.), Ch.54「Eclipses」（l1, l2, f1, f2 構成）。
- Espenak/NASA（GSFC eclipse / EclipseWise）: annular（反本影）定義・本影頂点 vs 地球面（data-sources §4.1, 第二義照合）。
- conventions §4/§4.1/§5/§9/§11, accuracy.md §0/§2.1/§2.2/§3.1/§4, numerical-policy §A5, algorithms.md §0。
- §3（見かけ地心位置・距離 ISSUE-015）, §6（基底・ベッセル要素 ISSUE-020/021）, §8（種別確定 ISSUE-023）。
- 関連 Issue: ISSUE-019, ISSUE-020, ISSUE-021, ISSUE-023。
- **要確認**: Explanatory Supplement Ch.11 影錐式の式番号（C2–C4）。半影頂点の §6 での使用要否（C5）。地球面交点を球/楕円体どちらで取るか（C6）。
