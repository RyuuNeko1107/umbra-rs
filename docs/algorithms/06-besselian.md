# §6 基本面・瞬時ベッセル要素（Step6-7）

> 正本: `algorithms.md §0`（記号表 = x,y,d,μ,l1,l2,f1,f2 の符号・単位 Re・**l2<0=皆既**）・`conventions.md §5`（FundamentalPlane）・`§5.2`（μ=GAST−α, ERA/CIO, **分点 GST 禁止**）・`§9`（半径 k）・`accuracy.md §0(a)`（μ のみ UT1 依存）。本セクションはこれらに**厳密準拠**する。
> 状態: ドラフト（Milestone 0）。一次資料で未確認の式番号・符号は「**要確認**」を残し、推測で式番号を書かない。
> 関連 Issue: ISSUE-020（基本面基底）, ISSUE-021（瞬時ベッセル要素）, ISSUE-037（直接評価器 BesselianSource）。供給連携: ISSUE-039（μ・GAST 素材, ERA/CIO）。

---

## 目的と入出力（単位・時刻系・座標系・フレーム）

**1つの TT 時刻**における瞬時ベッセル要素 `InstantaneousBesselianElements`（x, y, d, μ, l1, l2, tan f1, tan f2）を、影円錐（§5/ISSUE-019）と基本面基底（ISSUE-020）から算出する（architecture §6, api-draft §3.3, ISSUE-021）。これを**任意 TT で連続供給する層**が `BesselianSource`（直接評価器, ISSUE-037）。

| 項目 | 内容 |
|---|---|
| 入力 | `ShadowCone`（§5）, `FundamentalPlaneBasis`（基底・d・α_axis, ISSUE-020）, 月の見かけ地心位置（km, §3）, `ERA`（影軸時角 μ 構成, ISSUE-039 供給, UT1 依存）, `Re = WGS84 a`（km） |
| 出力 | `InstantaneousBesselianElements`: `x, y`（Re 無次元）, `d`(rad), `μ`(rad), `l1, l2`（Re 無次元）, `tan_f1, tan_f2`（無次元）, `time_tt` |
| 単位 | x,y,l1,l2 = **Re 無次元**（conventions §1/§4.1, Re = WGS84 a）。tan f1/f2 無次元。d, μ = rad |
| 時刻系 | **TT 基準**（conventions §6, accuracy.md §0(a)）。**μ の恒星時(ERA)部分のみ UT1 由来**。要素自体は TT 時刻でラベル |
| 座標系/フレーム | **FundamentalPlane**（conventions §5）。**X̂=東向き, Ŷ=天の北の基本面射影, Ẑ=影軸（地心→太陽向き）**（ISSUE-020, §0） |

**確定方針（厳守, §0/ISSUE-021）**: x=東正, y=北正。`l2 < 0 ⇒ 皆既`（金環は `l2 > 0`）。μ は **ERA 経由 CIO ベース**で構成し、**分点 GST 混在禁止**（conventions §5.2 / D4）。直接評価器（ISSUE-037）が任意 TT で本カーネルを呼んで供給する（fit 誤差ゼロ）。

---

## 記号（algorithms.md §0 参照。本節固有の補助記号のみ定義）

§0「ベッセル要素」表（x, y, d, μ, l1, l2, f1, f2）を正本とする。本節固有・補助:

| 記号 | 意味 | 単位 |
|---|---|---|
| `X̂, Ŷ, Ẑ` | 基本面基底（東/北射影/影軸=地心→太陽向き, ISSUE-020） | — |
| `α_axis` | 影軸赤経（CIRS 基準, ISSUE-020 の α_z, §2 と同記号） | rad |
| `d` | 影軸赤緯（= Ẑ の赤緯, ISSUE-020, §0） | rad |
| `g` | 影軸の地心位置ベクトル（基本面交点を張る量。月中心 r_moon を基本面へ射影） | km |
| `GAST` | 見かけ恒星時（ERA 経由 CIO ベース, §2 F9, conventions §5.2） | rad |
| `ζ_apex1`, `ζ_apex2` | 半影/本影 頂点の基本面（z）座標（Re。l1/l2 構成の中間量） | Re |
| `f1, f2` | 半影/本影 半頂角（§5 C2/C3 で算出済み, §0） | rad |
| `Re` | 地球赤道半径 = WGS84 a = 6378.137 km（無次元化基準, conventions §4.1） | km |

> 記号衝突注意: 本節の `g`（影軸地心位置, 中間量）は §4 の連続化合関数 `g(t)` とは別物（衝突回避のため本節は添字なし km ベクトルとして明示）。`α_axis` は §2 の `α_axis`・ISSUE-020 の `α_z` と同一。

---

## 数式（番号付き・各式に出典）

出典（共通）: **第一義** Explanatory Supplement to the Astronomical Almanac (3rd ed.), Ch.11「Eclipses」Besselian elements 節（x,y,d,μ,l1,l2,f1,f2 の正式定義）。**補助** Meeus, *Astronomical Algorithms* (2nd ed.), Ch.54「Eclipses」（実用式）。**NASA 定義** Espenak（GSFC eclipse / EclipseWise / NASA TP-2006-214141, data-sources §4.1）。**要確認**: Explanatory Supplement Ch.11 の x,y,l1,l2,μ の正確な式番号（一次で確定）。

### 6.1 基本面基底（z=影軸, y=北の射影, x=右手系東）— ISSUE-020

**(B1) 影軸方向 → 赤経・赤緯**
出典: ISSUE-020（基本面定義, Explanatory Supplement Ch.11 / Meeus Ch.54）。**影軸単位ベクトル Ẑ は地心→太陽向き**（ベッセル慣習, §5 の `axis_direction = −û` の逆向き＝`+û` 方向に対応, §5 記号注記参照）:
```
Ẑ      = (cos d · cos α_axis,  cos d · sin α_axis,  sin d)        # 赤道直交系（CIRS, of date）
α_axis = atan2( Ẑ_y, Ẑ_x )                                       # 影軸赤経（CIRS 基準, CIO 一貫）
d      = asin( clamp( Ẑ_z, −1, 1 ) )                             # 影軸赤緯
```

**(B2) 北射影 Ŷ と 東 X̂（右手系）**
出典: ISSUE-020（標準ベッセル基底, Meeus Ch.54）。**北 = 天の北極（赤道系 (0,0,1)）, 黄道北ではない**（ISSUE-020 レビュー重点）:
```
n = (0,0,1)                                  # 天の北極
Ŷ = normalize( n − (n·Ẑ) Ẑ )                # 天の北を基本面へ射影（北向き）
X̂ = normalize( Ŷ × Ẑ )                      # 右手系・東向き
```
- 直交性・右手性を検証（`X̂·Ŷ ≈ 0`, `det[X̂ Ŷ Ẑ] ≈ +1`, ISSUE-020）。
- **極端配置の縮退（最重要, architecture §6）**: `Ẑ ≈ (0,0,±1)`（天の極）で `n − (n·Ẑ)Ẑ ≈ 0` → Ŷ 定義不能。縮退を検出し代替射影 or `DegenerateGeometry`（現実の日食は太陽が黄道上 |δ|≤23.5° ゆえ発火しないが防御必須, ISSUE-020）。

### 6.2 瞬時ベッセル要素（x, y, d, μ, l1, l2, tan f1, tan f2）— ISSUE-021

**(B3) x, y（影軸の基本面交点, Re 無次元）**
出典: Meeus Ch.54 の x,y 定義（東正/北正）/ Explanatory Supplement Ch.11（**要確認**: 式番号）。影軸の地心位置ベクトル `g`（月中心 r_moon を採用）を基本面 (X̂, Ŷ) へ射影し **Re で割って無次元化**:
```
x = ( g · X̂ ) / Re          # 東向き成分 / Re（無次元）
y = ( g · Ŷ ) / Re          # 北向き成分 / Re（無次元）
```
- 符号: x=東正, y=北正（§0, ISSUE-021）。`x = y = 0` は影軸が地心を貫く完全中心食。
- **要確認**: 影軸の基本面交点を「月中心 r_moon」で取るか「影軸と基本面 z=0 の交点」で取るか（厳密定義は Explanatory Supplement Ch.11 一次で確定。NASA は影軸と基本面の交点。ISSUE-021 §出典）。

**(B4) d（影軸赤緯）**
出典: ISSUE-020（B1）/ §0。
```
d = （B1 の影軸赤緯。基底 declination と同一）
```

**(B5) μ（影軸の見かけグリニッジ時角, ERA 経由 CIO ベース）**
出典: **§0・conventions §5.2 / D4**・§2 (F4)・ISSUE-021/039。
```
μ = θ_ERA − α_axis           ( rad, [0,2π) 正規化 )
```
- `θ_ERA = ERA(UT1)`（§2 F4, `iauEra00`）は **CIO の見かけグリニッジ時角**。`α_axis` は **CIRS 基準**（B1, CIO 起点で測った影軸赤経）。ともに CIO 起点で測るため、**時角 μ = θ_ERA − α_axis は CIO で完全に閉じる**（EO も GAST も登場しない）。
  - 補足: 分点系で書けば `μ = GAST − α_equinox`。CIO ↔ 分点は `GAST = θ_ERA − EO` かつ `α_equinox = α_axis − EO` の関係にあり、`μ = (θ_ERA − EO) − (α_axis − EO) = θ_ERA − α_axis` と **EO は相殺する**。旧記述 `μ = GAST − α_axis(CIRS)` は GAST（分点起点）と α_axis（CIO 起点）を混ぜて EO を二重計上する誤りだった（D4「分点 GST 禁止」＝CIO/分点の二重定義排除に反する）。**正本は `μ = θ_ERA − α_axis`**（ISSUE-021/039 実装で確定）。
- **μ は UT1（→δUT1）依存**（θ_ERA 経由）。全球 gamma・最大食「時刻 TT」が純 TT で閉じるのに対し、**μ を介する局地接触時刻には δUT1 が混入**する（accuracy.md §0(a) 脚注・§2.1L, ISSUE-021 I1）。要素は TT 時刻でラベルするが μ のみ UT1 依存である旨を doc/metadata に明記（ISSUE-039/021）。
- 正規化: μ は赤経基準で `[0,2π)`（conventions §2, two_pi 規約）。
- **旧記述の分点 GAST 直接構成は採らない**（CIO/分点二重定義の排除, D4）。

**(B6) tan f1, tan f2（円錐半頂角の tan）**
出典: §5 C2/C3（半角 f1, f2）/ §0/ISSUE-021。
```
tan_f1 = tan( f1 )          # 半影（発散）: 正の大きい値
tan_f2 = tan( f2 )          # 本影（収束）: 正の小さい値
```

**(B7) l1, l2（基本面での半影/本影 円錐半径, Re 無次元）**
出典: Explanatory Supplement Ch.11 の l1, l2 定義式（**要確認**: 式番号）/ Meeus Ch.54 / ISSUE-021。基本面（z=0 平面）での錐の切口半径を **Re 無次元**で:
```
l1 = ( 基本面での半影外縁半径 ) / Re        # 半影外半径
l2 = ( 基本面での本影縁半径 )   / Re        # 本影半径。金環で負符号
```
- 概念形（頂点と半角からの切口半径, ISSUE-021 §出典）: 基本面（z=0）から各頂点までの z 距離 `ζ_apex` と半角 `f` から `l = ζ_apex · tan f`（中間量 `ζ_apex1`/`ζ_apex2`）。**正確な無次元化・符号付き定義は Explanatory Supplement Ch.11 を一次で確定（要確認）**。
- **符号規約（必須, §0/ISSUE-021）**: `l2 < 0 ⇒ 皆既`（本影頂点が基本面の**反地球側**＝地球面以遠に届く, §5 C6 `OnOrBeyondEarth`）。`l2 > 0 ⇒ 金環`（本影頂点が基本面の**地球側**＝届かず反本影が地表, §5 C6 `BeforeEarthSurface`）。皆既↔金環境界で `l2` は連続に符号反転する（§8/ISSUE-023 種別境界の素, ISSUE-021）。
  > 注意（向き定義の明示）: §0 の符号規約「l2<0=皆既」を本実装の正本とする。NASA/Espenak 公開値の l2 符号と整合するよう構成し、もし符号が逆なら NASA 対応表（下記 D5）で系統差として明記する（誤差に誤帰属しない, accuracy.md §0）。**要確認**: NASA 公開 l2 の符号（total/annular）と本実装の一致を実値で確認（ISSUE-021 検証）。

### 6.3 直接評価器（任意 TT 供給）— ISSUE-037

**(B8) BesselianSource::at(time_tt)**
出典: ISSUE-037（オーケストレーション層, 式は持たない）。
```
at(time_tt):  §3(ISSUE-015 見かけ位置) → §5(ISSUE-019 影円錐) → B1/B2(ISSUE-020 基底)
              → ISSUE-039(GAST/ERA) → B3..B7(ISSUE-021 besselian_elements_at)
              ⇒ InstantaneousBesselianElements（毎回暦再評価 ⇒ fit 誤差ゼロ）
```
- 各呼出は純関数的・副作用なし（キャッシュ無し, ISSUE-037）。Standard 局地計算の**既定供給源**（多項式 §7/ISSUE-022 と `BesselianSource` 契約で差し替え可能）。

---

## NASA 表記との対応（D5・必須記載, ISSUE-021）

NASA/Espenak ベッセル要素と本実装の単位・規約対応（系統差を幾何誤差に誤帰属しないため, accuracy.md §0）:

| 量 | NASA/Espenak 表記・単位 | 本実装 | 変換・規約 |
|---|---|---|---|
| x, y | Re 無次元（基本面座標, 東正/北正） | Re 無次元 | 同一 |
| d | 度（影軸赤緯） | Radians | 度→rad |
| **μ** | **度**（基本面の経度基準 ※hour で公開されることもあり, **要確認**） | Radians | 度→rad。**D4: ERA(iauEra00)経由 CIO ベース**（分点 GST 禁止） |
| **μ′（μ の時間変化率）** | **≈15°/hour**（地球自転 ≈15.041°/hour + 影軸赤経変化分。NASA 度・hour 系, **要確認**） | Radians/SI秒（内部微分） | hour→SI秒（μ′ ≈ 15°/3600 s を基準換算） |
| l1, l2 | Re 無次元（半影/本影半径） | Re 無次元 | 同一。**正本(B1): l2<0=皆既 / l2>0=金環**（§0 と一致。NASA 値の符号が逆なら D5 対応表で系統差記録） |
| tan f1, tan f2 | 無次元 | 無次元 | 同一 |
| ΔT 規約 | NASA カタログ採用 ΔT（Espenak/Meeus） | conventions §9 / §1 | 比較時に NASA の ΔT へ揃える（系統差を記録） |
| GAST/時角 供給源 | NASA は分点系の場合あり | **CIO(ERA)基準・ISSUE-039 供給** | 分点 GST 禁止（D4）。CIO↔分点の系統差を記録 |

- **μ の NASA 単位（度 vs hour）と μ′≈15°/hour の正確な定義は一次資料（NASA TP-2006-214141 / GSFC eclipse site）で最終確認（要確認）**。系統差は誤差として隠さず accuracy.md に記録（ISSUE-021/039）。

---

## 手順（実装順・数値注意）

1. **基底 (B1/B2, ISSUE-020)**: 影軸方向（§5 から地心→太陽向き Ẑ）→ α_axis, d → Ŷ（北射影）→ X̂（東, 右手系）。直交性・det=+1 検証。極端配置の縮退防御。
2. **x, y (B3, ISSUE-021)**: 影軸地心位置 `g` を (X̂, Ŷ) へ射影 → Re で割り無次元化。
3. **d (B4)**: 基底 declination。
4. **μ (B5)**: `GAST = θ_ERA − EO`（§2 F9, ERA/CIO, ISSUE-039 供給）→ `μ = GAST − α_axis` を `[0,2π)` 正規化。**分点 GST 禁止**（D4）。
5. **tan f1/f2 (B6)**: §5 の f1, f2 から。
6. **l1, l2 (B7)**: 基本面切口半径を Re 無次元で。**金環: l2 を負符号で出す**（§5 C6 と整合）。
7. **構成**: `InstantaneousBesselianElements`（`time_tt` を TT でラベル）。
8. **直接供給 (B8, ISSUE-037)**: 上記を任意 TT で再評価（fit 誤差ゼロ）。

**数値注意（横断, numerical-policy §A5）**:
- asin/acos 引数を `[-1,1]` クランプ（accuracy.md §2.2）。
- μ の正規化規約を固定（`[0,2π)`, conventions §2）。μ→μ+2π 不変をプロパティテスト（§2/ISSUE-039）。
- l2 の符号は皆既↔金環境界で**連続に反転**すること（§8/ISSUE-023 種別境界の素, ISSUE-021）。
- 基底構成は f64 機械精度（直交性残差 < 1e-12, 中心線 sub-km へ余裕, ISSUE-020）。
- μ の恒星時部分は UT1（ERA）由来。TT 要素に UT1 が混じる点を実装メモ/metadata に明記（accuracy.md §0(a), ISSUE-021 I1/ISSUE-039）。

---

## 境界・特異・異常系

- **極端配置（Ẑ ≈ 天の極, B2）**: Ŷ 定義不能 → 代替射影 or `DegenerateGeometry`。現実の |d|≤23.5° では正常（fuzz/property L8/L9 で踏むため必ず実装, ISSUE-020）。
- **完全中心食 (x=y=0)**: 正常値（影軸が地心貫通）。MockEphemeris で検証（ISSUE-021）。
- **皆既↔金環境界 (l2=0)**: `sin f2 ≈ 0`（§5 頂点無限遠）と対応。l2 が連続に符号反転（ISSUE-021/023）。
- **μ の ±2π 折返し**: `[0,2π)` 正規化を明示、μ→μ+2π 不変（ISSUE-039）。求解側（§4/§8）で連続性が要るときは呼出側で連続化（本層は素の値, §2 と同方針）。
- **ERA/UT1 欠損**: EOP データ層の責務（§2/ISSUE-035）。欠損は上流へ（本層は値を受けて μ を組むのみ）。
- **z 向き規約の取り違え**: Ẑ が地心→太陽の逆だと x,y 符号が全反転（鏡像）。§5（axis_direction=−û）との符号橋渡しを最重要レビュー項目とする（ISSUE-019/020/021）。

---

## 検証（基準値の出典。実装へ値コピー禁止 = conventions §11）

accuracy.md §3.1、ISSUE-020/021/037 受入テスト準拠。基準値は fixtures/DE/NASA から動的取得（ハードコード禁止, conventions §11）。GPL 実装の数値を期待値に貼らない（基準は DE/数式オラクル）。

- **基底（ISSUE-020, L1/L4）**: 任意影軸で `|X̂|=|Ŷ|=|Ẑ|=1`・相互直交・`det=+1`（プロパティ L8 で多数ランダム Ẑ）。x=東（赤経増加方向）・y=北（赤緯増加方向）を既知方向で確認（オラクル＝幾何手計算, 非コピー）。極端配置で NaN/縮退せず安定 or `DegenerateGeometry`（必須, architecture §6）。既知影軸（太陽が春分点方向）で d≈0・α_axis 既知。
- **瞬時要素（ISSUE-021, L4）**: NASA 5千年カタログの既知日食（皆既/金環/部分/ハイブリッド各）の最大食付近 x,y,d,μ,l1,l2,tan f1,tan f2 を NASA 公開値と比較。**k 慣習を Espenak（EspenakUmbral/Penumbral）に揃え ΔT も合わせる**（系統差を accuracy.md に記録, 第二義, accuracy.md §3.1）。x,y は <0.001 Re 級（要 M2 実測）。
- **DE 差分（第一義, accuracy.md §3.1）**: 解析暦 vs DE440 で同一ベッセル計算を通し x,y,l1,l2 の差を層分解（暦残差 vs 幾何残差, §4）。
- **MockEphemeris 人工ケース（accuracy.md §3.1）**: 完全中心（x=y=0）, 明確な金環（**l2>0**）, 明確な皆既（**l2<0**）, 影が地球を外す（|gamma|>1+l1）。各で符号・値を解析オラクルと照合。
  > 検証時の符号統一: §0 の正本「**l2<0=皆既 / l2>0=金環**」に全テストを合わせる。NASA 値の符号が逆なら D5 対応表で系統差として記録（誤差化しない, accuracy.md §0）。
- **符号規約テスト（必須・品質基準, ISSUE-021）**: x=東正・y=北正・`l2<0`(皆既)・`l2>0`(金環) を境界をまたいで**連続に**確認。NASA 対応表（D5）と一致。
- **単位テスト**: x,y,l1,l2 が Re 無次元（km 値/Re と一致）。
- **直接評価器（ISSUE-037, L4/L7）**: `InstantaneousEvaluator::at(t)` が同一入力で `besselian_elements_at`（ISSUE-021）と**完全一致**。直接（fit 誤差0）vs 多項式（§7/ISSUE-022）の x,y,l1,l2 残差を fit 区間で実測（L7, accuracy.md §3.2）。`&dyn BesselianSource` 経由で多項式と差し替え可能。
- **CIO 統一の証明（μ, ISSUE-039 連携）**: GAST（B5）が分点 GST（`iauGst06a`）と「CIO−分点の既知量」だけ差を持つ＝内部で分点経路を使っていない証明（§2 検証, 系統差は accuracy.md 記録）。μ→μ+2π 不変。ERA は UT1・CIP/s は TT で評価（時刻系分離回帰, ISSUE-039）。

---

## 許容誤差

- accuracy.md §2.1 幾何バジェット（**影幾何誤差**, §4 層分解）。最終 ±1.5s（最大食）・食分 ±0.0005 へ寄与。
- x,y の誤差は最大食時刻・gamma に直結（感度 0.5″/s, 1″≈2s）。l1,l2 の誤差は食分（0.001食分≈1.9″, accuracy.md §2.2）。
- **直接計算（ISSUE-021/037）の fit 誤差はゼロ**（暦再評価）。多項式（§7/ISSUE-022）との残差は L7 サブテスト（accuracy.md §3.2）。
- 基底構成は純幾何で f64 機械精度（≪ バジェット, 直交性残差 < 1e-12, ISSUE-020）。誤差は入力影軸方向（§3/§5）律速。
- **μ のみ UT1（δUT1）依存**: 全球 gamma・最大食「時刻 TT」は純 TT、**局地接触時刻は μ 経由で δUT1 混入**（accuracy.md §0(a)/§2.1L, ISSUE-021 I1）。系統差（k/ΔT/CIO 慣習）は誤差として隠さず accuracy.md に記録（accuracy.md §0）。

---

## 出典

- Explanatory Supplement to the Astronomical Almanac (3rd ed.), Ch.11「Eclipses」Besselian elements 節（x,y,d,μ,l1,l2,f1,f2 正式定義・基本面定義）。**要確認**: 式番号。
- Meeus, *Astronomical Algorithms* (2nd ed.), Ch.54「Eclipses」（x,y,l1,l2,f1,f2 実用式）。
- Espenak/NASA: Besselian Elements 解説（GSFC eclipse / EclipseWise）・NASA TP-2006-214141（μ 単位・μ′・k 慣習, data-sources §4.1, 第二義照合）。
- conventions §5/§5.2/§6/§9/§11, accuracy.md §0/§2.1/§2.2/§3.1/§4, numerical-policy §A5, algorithms.md §0。
- §2（GAST=θ_ERA−EO, μ=GAST−α_axis, F4/F9/F10/F11）, §3（見かけ位置 ISSUE-015）, §5（影円錐・半角・金環判定 ISSUE-019）, §7（多項式 ISSUE-022）, §8（種別 ISSUE-023）。
- 関連 Issue: ISSUE-020, ISSUE-021（D5）, ISSUE-037, ISSUE-039（μ/GAST 供給）。
- **要確認**: Explanatory Supplement Ch.11 の x,y,l1,l2,μ 式番号（B3/B5/B7）。影軸の基本面交点を月中心 r_moon で取るか影軸×z=0 交点で取るか（B3）。NASA 公開 l2 符号（total/annular）と §0「l2<0=皆既」の実値一致（B7/D5）。μ の NASA 単位（度 vs hour）と μ′≈15°/hour 定義（D5）。
