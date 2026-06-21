# §11 部分食域（Partial eclipse domain, Milestone 9 残(3)）

本節は `EclipsePath::partial_limit`（`Option<GeoPolygon>`, api-draft §4）の**部分食が見える全球領域の外周**を、
ベッセル要素から第一原理で構成する数式・手順・サブスライス分解を定める。中心線・南北（本影）限界線・
経路サンプル列は M9.1–M9.7 で実装済み（ISSUE-045 / §8 全球条件を参照）。本節は**残り (3) 部分食域**の設計。

正本は `algorithms.md §0`（記号表）と §8（gamma・全球接触・限界線の包絡条件）。本節はそれに厳密準拠する。

> 状態: **設計ドラフト**（2026-06-21・未実装）。概念出典は Explanatory Supplement to the Astronomical
> Almanac Ch.11 / Meeus *Astronomical Algorithms* Ch.54 / NASA Espenak（部分食限界・rise/set 曲線）。
> **一次資料の式番号は未確認のため転記しない**（house rule・M9.4/M9.6 と同方針）。数式は第一原理
> （錐∩楕円体・移動半影の包絡・半影縁∩terminator）から導出し、数値オラクルは NASA 2024-04-08 公開
> limits 表を用いる。未確定箇所は「**要確認**」を残す。

---

## 目的と入出力

- **目的**: ベッセル要素（多項式 ISSUE-022 / 直接供給 ISSUE-037）から、**いずれかの時刻に部分食以上が
  見える地表点の集合**の外周を閉じた `GeoPolygon`（外環＋必要なら穴・反子午線跨ぎは MultiPolygon）として返す。
- **入力**: `BesselianSource`（x, y, d, μ, l1, l2, tan f1, tan f2）、全球接触 P1/P4（§8・部分食区間 [P1,P4]）、
  `PathOptions`（`include_limits`・`split_antimeridian`）、WGS84。
- **出力**: `EclipsePath.partial_limit: Option<GeoPolygon>`。中心食/部分食を問わず部分食域は存在しうる
  （部分食 only の日食でも `partial_limit` は Some になりうる＝中心食限定の `center_line` とは独立）。
- **単位/時刻系/フレーム**: §8 と同一（角度 rad、ベッセル長 Re、FundamentalPlane〔影軸 z・x̂=東・ŷ=天の北射影〕、
  地表点は WGS84 測地座標へ）。時刻は各境界点の TT（必要なら UTC 併記は後続）。

---

## 記号（algorithms.md §0・§8 を正本。本節固有の補助のみ）

| 記号 | 意味 | 出典/備考 |
|---|---|---|
| `L1(ζ) = l1 − ζ·tan f1` | ζ補正**半影**半径（Re） | §8 と同形。本影 `L2(ζ)=l2−ζ·tan f2` の半影版 |
| `rel` | 影の地表に対する相対速度（基本面） | §8 限界線・M9.4 と同一: `rel=(x′−μ′(ζcosd−η sind), y′−μ′ξ sind)` |
| 半影限界 | 移動する**半影**縁の包絡＝部分食域の南北端（昼面） | 式 11.1（M9.4 本影限界の l1 版） |
| terminator | 基本面 ζ=0 の地表大円（太陽が地平＝日の出入り） | 式 11.2 |
| rise/set 曲線 | 半影縁が terminator 上を通る点の時系列軌跡（日の出/日没に食接触） | 式 11.3 |

---

## 数式（第一原理・番号は本節内ローカル）

### 11.0 部分食域の定義

地表点 P=(ξ,η,ζ)（WGS84 上）は時刻 t に**部分食以上**を見る ⇔ 半影内:
```
(ξ − x)² + (η − y)² ≤ L1(ζ)²        L1(ζ) = l1 − ζ·tan f1
```
等号が半影縁（P 接触＝部分食の開始/終了）。部分食域 = `⋃_{t∈[P1,P4]} { P : 上式 }` の地表射影。
その**外周**が `partial_limit`。外周は次の 2 曲線族が交互に連なって閉じる:
**(A) 南北半影限界（昼面の包絡）** と **(B) rise/set 曲線（terminator 上の半影縁）**。

### 11.1 南北半影限界（昼面の包絡）— 式 11.1

部分食域の南北端は、移動する**半影**縁の包絡。各縁点 P は §8 限界線と同一の 2 条件を満たす（半径を
本影 L2 から半影 L1 に替えるのみ）:
- **(1) 錐exact（自己整合ζ）**: `(ξ−x)²+(η−y)² = L1(ζ)²`、`L1=l1−ζ·tan f1`（ζ は P 自身）。
- **(2) 包絡（路限界方向）**: 軸からのオフセットが `rel` に直交。

> 実装: 既存 `axis_intercept::solve_limit_edge` は半径に本影 `(l2, tan f2)` をハードコードしている。
> **半径 `(radius_l, tan_f)` を引数化**して半影 `(l1, tan f1)` でも解けるよう一般化する（本影限界 M9.4 は
> `(l2, tan f2)` を渡す呼び出しに退化）。不動点反復・WGS84 楕円体投影・南北割当は M9.4 と完全同一。

**昼面限定の注意**: 半影は大きく（l1 ≈ 0.53 Re）、|gamma| が大きい日食では半影限界が地表（昼面 ζ>0）を
外れて terminator（ζ=0）に達する。半影限界は **ζ>0 の昼面でのみ存在**し、ζ→0 で rise/set 曲線（11.3）へ
連結する。`solve_limit_edge` が `RootNotBracketed`（縁が地表を外す）を返す端が連結点の近傍。

### 11.2 terminator（日の出入り境界）— 式 11.2（**WGS84 厳密**・要確認#2 解決）

基本面 ζ=0 の地表点（太陽が地平＝日の出入り）の軌跡。`surface_point_for_fundamental` の残差
`r(ζ)=ρcos²+(ρsin/(1−f))²−1`（子午線楕円拘束）に ζ=0 を代入すると、px=−η sin d・py=ξ・pz=η cos d より
```
ξ² + k·η² = 1,      k = sin²d + cos²d/(1−f)²        （WGS84 terminator 楕円・基本面 (ξ,η)）
```
（球 f=0 で k=1 ⇒ 単位円。扁平 f>0 で k≥1 ⇒ η 方向に縮む楕円＝極扁平）。**球近似（k=1, 単位円）は
中心線・限界線の WGS84-exact パイプラインと不整合**（球面上に置いた点を WGS84 前方射影で戻すと ζ≠0、
~Re·f≈21 km の残差）になり往復オラクルが緩むため、本節は **terminator を上記楕円で厳密化**する（要確認#2 解決）。

### 11.3 rise/set 曲線（terminator 上の半影縁）— 式 11.3（**WGS84 厳密**）

時刻 t に「部分食の接触が日の出/日没ちょうどに起こる」点 = 半影縁 ∩ terminator（ζ=0）。基本面で **円**
（半影縁）∩ **楕円**（terminator・11.2）を厳密に解く:
```
楕円（terminator）:  ξ² + k·η² = 1                       （k = sin²d + cos²d/(1−f)²）
円（半影縁・ζ=0）:   (ξ−x)² + (η−y)² = l1²               （ζ=0 で半影半径は L1(0)=l1）
```
2 二次曲線の交点（幾何的に通常 ≤2）。楕円から `ξ = ±√(1−k·η²)`（|η|≤1/√k）を円へ代入した η の 1 変数
残差を、既存の粗走査＋Brent（`descending_sign_change_bracket`＋`brent_root`・降順走査機構）で根求めし、
各根 (ξ,η) を `fundamental_to_geodetic(ξ,η,0,d,μ)` で測地座標化する。**点は 11.2 楕円上にあるので WGS84
前方射影で ζ≈0 へ往復一致**＝中心線/限界線と同じ WGS84-exact オラクルで検証できる（球近似の不整合を排除）。

t を [P1,P4] で動かすと 2 根が 2 曲線（朝側・夕側 limb）を描く。P1・P4 では円が楕円に接し交点 1 点（外周の
尖点＝半影が地球に最初/最後に触れる点）。**rise/set/begin/end の 4 分類は出力ラベル（metadata）であり、
外周ポリゴン自体には不要**（要確認#1 解決）＝外環は 2 根曲線の**両方**（朝側 limb 弧・夕側 limb 弧）を使う。
朝/夕の別は経度の subsolar 点に対する東西（自転で先行/後行する側）で決め、begin/end は半影の進入/退出
（in-plane で半影が点を覆い始める/終わる）で決める（必要なら後続でラベル付け）。

### 11.4 外周ポリゴンの構成 — 式 11.4

外環 = 北半影限界（昼面・西→東）→ 日没側 rise/set 曲線（terminator を北→南）→ 南半影限界（東→西）→
日の出側 rise/set 曲線（南→北）→ 始点へ閉じる、を時刻順・地理順に連結。
- 連結点: 半影限界が ζ→0 で terminator に達する点（11.1 の `RootNotBracketed` 端）と rise/set 曲線端を接合。
- 反子午線（|Δlon|>180）跨ぎは MultiPolygon に分割（GeoLine の M9.5 ±180 補間をポリゴン環へ拡張）。
- 環の向きは RFC 7946（外環 CCW・穴 CW）。**要確認**: 高緯度/極で領域が terminator に張り付き
  単連結でなくなる退化（穴の有無）。通常の中緯度日食は単連結 1 外環。

---

## サブスライス分解（実装順・各 strict）

- **(3a) 半影限界の一般化** ✅（2026-06-21 実装済み）: `solve_limit_edge` を錐半径引数化（`cone_l`/`cone_tan_f`・
  `(l1,tan f1)`/`(l2,tan f2)` 両対応）。本影限界 M9.4 は退化呼び出しで完全回帰。mutation 52/51 caught・0 missed
  （docs/reviews/mutation-limit-line.md 追記）。**API 露出なし**（南北半影限界点列の生成・partial_limit は (3c) で組む）。
- **(3b) rise/set 曲線** ✅（2026-06-21 実装済み）: `cone_terminator_intersections`＝**円∩terminator 楕円**（WGS84 厳密）。
  楕円を θ 媒介（ξ=cosθ, η=sinθ/√k）した円残差の符号反転を粗走査＋Brent（機構は `scan_periodic_sign_change_roots`
  に分離）→`fundamental_to_geodetic(ξ,η,0,…)`。WGS84 前方射影の往復（ζ≈0・面内距離=cone_l）＋d=π/2 二円閉形式で
  検証。mutation 40/39 caught・0 missed（機構は wholesale 除外・docs/reviews/mutation-rise-set.md）。**API 露出なし**
  （t を [P1,P4] で動かした曲線化・接点端の扱いは (3c)）。
- **(3c) 外周組立**: (3a)+(3b) を 11.4 の順序で連結し `GeoPolygon` 外環を構成。反子午線 MultiPolygon 分割。
  `partial_limit` を `path()` で Some に（部分食 only の日食含む）。オラクル: 実 2024 で外環が NASA 部分食
  限界（北限/南限/rise/set）の ballpark を内包、頂点数・閉性・環向き。
- **(3d) GeoPolygon GeoJSON**: `GeoPolygon::geojson_geometry`（Polygon／反子午線 MultiPolygon・環閉包
  RFC 7946 §3.1.6）＋ `EclipsePath::to_geojson` に `partial_limit` feature（`role="partial_limit"`）を
  決定的順序で追加。オラクル: 既存 M9.2/M9.5 GeoJSON テストにポリゴン版を追加。

---

## 受け入れテスト戦略

- **FAST（合成・機械精度）**: 半影限界 2 条件（11.1）、円∩terminator 楕円の往復・閉形式（11.3）、外環の閉性・向き・
  反子午線分割（11.4）、GeoJSON 構造（3d）。値は非対称にして取り違え変異を撃つ。
- **SLOW（実 2024-04-08）**: 外周が NASA 公開部分食限界（北限・南限・eclipse begins/ends at sunrise/sunset）の
  代表座標を ballpark 域で内包。中心線・本影限界（実装済み）との整合（partial ⊃ umbral path）。
- **二段オラクル（ISSUE-047）**: M2 暫定（NASA limits 表 ballpark）／M10 最終（DE 差分）。

---

## 許容誤差

- 半影限界・外周位置: 中心線（accuracy §2.1 sub-km）より緩く、**ballpark（数十 km 規模）**で可（半影縁は
  本影縁より太陽縁の鋭さが鈍く、NASA 公開値の桁も粗い）。fit 残差は `BesselianPolynomial.fit_error` でガード。
- rise/set 曲線は terminator を WGS84 楕円（11.2）で厳密化したので中心線・限界線と同精度域。
  fit 残差は `BesselianPolynomial.fit_error` でガード。**近似を残す箇所は明記**（conventions §11）。

---

## 要確認（一次資料・設計判断の未決事項）

1. ~~rise/set 曲線の 4 分類の取捨~~ **解決（11.3）**: 外周ポリゴンは 2 根曲線（朝側・夕側 limb）の**両方**を使う。
   begin/end×sunrise/sunset の 4 分類は出力ラベル（metadata）で外周自体には不要。
2. ~~terminator の WGS84 厳密化~~ **解決（11.2/11.3）**: 球近似は WGS84-exact パイプラインと不整合（~21km・
   往復オラクルが緩む）ため、terminator を楕円 `ξ²+k·η²=1`（k=sin²d+cos²d/(1−f)²）で厳密化。円∩楕円を求根。
3. 部分食域の**単連結性**（高緯度での穴・退化、11.4）。通常の中緯度日食は単連結 1 外環の想定。
4. NASA 2024 公開 limits 表の rise/set 曲線座標の入手（SLOW オラクルの粒度。当面は北限/南限/帯と
   partial ⊃ umbral path の包含で代用可）。
5. ~~`solve_limit_edge` 引数化の M9.4 mutation 影響~~ **解決（3a 実装済み）**: 退化呼び出しで本影回帰・
   mutation 52/51 caught・0 missed を確認（docs/reviews/mutation-limit-line.md）。
