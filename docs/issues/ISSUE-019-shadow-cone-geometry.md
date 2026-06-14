# ISSUE-019: Shadow cone geometry（本影/半影/反本影円錐・金環の頂点判定）

- crate: umbra-eclipse
- 依存: ISSUE-015（見かけ地心位置・距離）, ISSUE-003（Vector3, UnitVector3）, ISSUE-002（Radians）, ISSUE-010（WGS84/Re）, umbra-core（TtInstant, EarthModel）
- モード(tdd-workflow): strict（影円錐の頂点位置と半角はベッセル l1/l2/tan f1/tan f2（ISSUE-021）の直接の素材。本影頂点の地球側/反地球側判定が皆既/金環の境界を決める。幾何の誤りが種別誤判定に直結。strict）

## 目的
太陽（光源・有限サイズ）と月（遮蔽球）から、**本影（umbra）・半影（penumbra）・反本影（antumbra, 本影の延長錐）** の円錐幾何を構成する（architecture §6 `ShadowCone`）。
- 円錐の頂点位置・軸方向・半角を算出。
- **金環食判定の核**: 本影頂点が地球の手前（地球側）か向こう（反地球側）かを明示判定（architecture §6「金環食は本影頂点の地球側/反地球側を明示判定」）。

## 非目的
- ベッセル基本面の基底構成（ISSUE-020）。本 Issue は地心ベクトルでの円錐幾何まで。
- ベッセル要素 x,y,l1,l2,tan f1,tan f2 への変換（ISSUE-021。本 Issue の半角・頂点を Re 単位の l1/l2/tan f へ写像するのは ISSUE-021）。
- 種別の最終確定（ISSUE-023）。本 Issue は幾何量のみ供給。
- 観測者位置での影半径評価（局地, 別 Issue）。

## 公開インターフェース
architecture §6 の `ShadowCone` に準拠。crate 内部中心（`pub(crate)`）、検証用に `pub` 候補。

```rust
pub struct ShadowCone {
    pub axis_origin: Vector3,         // 影軸の基準点（月中心, 地心座標 km）
    pub axis_direction: UnitVector3,  // 月→太陽 反対向き（影の伸びる向き）
    pub umbra_apex: Vector3,          // 本影頂点（収束点）
    pub penumbra_apex: Vector3,       // 半影頂点（太陽-月の後方発散錐の見かけ頂点）
    pub umbra_half_angle: Radians,    // tan f2 の素（本影/反本影）
    pub penumbra_half_angle: Radians, // tan f1 の素（半影）
}

/// 本影頂点と地球の位置関係（金環判定）。
pub enum UmbraApexLocation {
    BeforeEarthSurface,  // 頂点が地球面の手前 → 本影が地球に達しない → 金環側
    OnOrBeyondEarth,     // 頂点が地球面以遠 → 本影が地球に達する → 皆既側
}

pub(crate) fn shadow_cone(
    sun_geocentric_km: Vector3,
    moon_geocentric_km: Vector3,
    sun_radius_km: f64,    // R_sun（conventions §9）
    moon_radius_km: f64,   // k·Re（conventions §9, k 選択）
) -> Result<ShadowCone, EclipseError>;  // DegenerateGeometry を返しうる

pub(crate) fn umbra_apex_location(cone: &ShadowCone, earth: EarthModel) -> UmbraApexLocation;
```

## 数式・アルゴリズムの出典
- **影円錐の幾何（外接共通接線による2球の円錐）**: Explanatory Supplement to the Astronomical Almanac (3rd ed.), **Ch.11「Eclipses of the Sun and Moon ...」の Besselian elements 節**（影錐の半角・頂点距離の定義）。および Meeus, *Astronomical Algorithms* (2nd ed.), **Ch.54「Eclipses」**（l1, l2, f1, f2 の構成）。
- **半角（half-angle）**:
  - 本影/反本影の半角 f2: `sin f2 = (R_sun − k·Re) / D`（D = 太陽-月中心間距離）。本影は太陽より月が小さいため錐が収束。
  - 半影の半角 f1: `sin f1 = (R_sun + k·Re) / D`。半影は発散錐。
  - 出典: Explanatory Supplement Ch.11 / Meeus Ch.54。**符号・どちらが収束/発散かを実装コメントに明記**（conventions §10）。NASA 表記 tan f1/tan f2 との対応は ISSUE-021 で確定。
- **本影頂点距離**: 月中心から本影頂点までの距離 `L = k·Re / sin f2`（収束点）。頂点 = 月中心 + L·(月→太陽の反対向き)。
- **金環/皆既の判定（頂点 vs 地球面）**: 本影頂点の地心距離と、その軸が地球面を貫く点での地球面までの距離を比較。頂点が地球面より手前なら本影が届かず**反本影（antumbra）**が地表に当たる＝金環。出典: Explanatory Supplement Ch.11（central eclipse の umbral cone reaching Earth の条件）/ NASA Espenak の annular 定義。
- **反本影（antumbra）**: 本影錐を頂点の先へ延長した発散錐。金環食の見かけ円環はこの反本影内。

## 単位 / 時刻系 / 座標系
- 単位: 距離 km（幾何計算, conventions §4）。半角はラジアン。Re=WGS84 a（conventions §4）。
- 時刻系: 入力は特定 TT の見かけ位置（呼出側が TtInstant で供給, conventions §6）。本 Issue は時刻非依存の純幾何。
- 座標系: 地心（GCRS/見かけ。ISSUE-015 出力）。FundamentalPlane への基底変換は ISSUE-020。
- 月半径は `k·Re`（k は conventions §9: IauMean=0.2725076 / EspenakUmbral=0.272281 / EspenakPenumbral=0.2725076）。**本影系は EspenakUmbral、半影系は EspenakPenumbral を使い分ける選択を config 経由で受ける**（NASA 照合時）。

## アルゴリズム概要
1. 太陽-月中心間距離 D、月→太陽方向 û を算出。軸方向 = −û（影が伸びる向き）。
2. 半影半角 `f1 = asin((R_sun + r_moon)/D)`、本影半角 `f2 = asin((R_sun − r_moon)/D)`（r_moon = k·Re）。
3. 本影頂点 `umbra_apex = moon_center + (r_moon / sin f2)·axis_direction`。半影頂点（後方発散錐の見かけ頂点）も同様に算出。
4. `umbra_apex_location`: 頂点の地心距離と地球面（軸の地球貫通点での半径）を比較し金環/皆既を判定。
5. 縮退検出: `D ≈ 0`、`R_sun ≤ r_moon`（非物理）、軸が地球を全く貫かない等で `DegenerateGeometry`（api-draft §3.5）。
- 数値安定性: asin 引数を `[-1,1]` クランプ（accuracy.md §2.2）。`sin f2 ≈ 0`（頂点が無限遠＝皆既/金環境界）で 0 除算回避（極限処理 or DegenerateGeometry 手前で大距離扱い）。極端配置（月が遠く本影が地球に届かない）で頂点距離が発散しても破綻しない。

## 受け入れテスト
accuracy.md テストレベル **L4（ベッセル幾何の素）** ＋ L1（純幾何）。
- **MockEphemeris 人工ケース（accuracy.md §3.1）**:
  - 明確な皆既配置（月近・本影頂点が地球面以遠）→ `OnOrBeyondEarth`。
  - 明確な金環配置（月遠・本影頂点が地球面手前）→ `BeforeEarthSurface`。
  - 境界（頂点が地球面ちょうど）→ ハイブリッド境界として判定の安定性確認（ISSUE-023 が最終分類）。
- 半角の解析検証（L1）: 既知の R_sun, r_moon, D で `sin f1/f2` を手計算（オラクル＝数式, 実装からコピーしない）と照合。本影が収束・半影が発散することを符号で確認。
- **DE/NASA 整合（第二義, data-sources §4.1）**: 既知の皆既/金環食で頂点判定が NASA の種別と一致（k 慣習を Espenak に揃える, accuracy.md §3.1）。
- k 値感度: `EspenakUmbral`(0.272281) と `IauMean`(0.2725076) で本影半角・頂点が変わり、皆既/金環境界の判定がずれることを定量確認（conventions §9 系統差の記録）。
- 縮退系: D≈0、R_sun≤r_moon → `DegenerateGeometry`。asin クランプ境界。

## 許容誤差
- accuracy.md §2.1 幾何バジェットの一部（影幾何誤差, §4 層分解）。半角・頂点の誤差は最終的に l1/l2（ISSUE-021）経由で食分・継続時間に効く（食分 0.001 ≈ 1.9″, accuracy.md §2.2）。
- 半角誤差は入力位置・距離の誤差（ISSUE-015 の月0.1″/太陽0.05″）で律速され、本 Issue の純幾何は f64 機械精度（≪ バジェット）で計算すること。
- **金環/皆既境界（ハイブリッド）付近は k 選択で結果が変わる**（conventions §9）。系統差を accuracy.md へ記録し、テストでは k を固定して比較（誤差を隠さない, accuracy.md §0）。

## 実装メモ
- `axis_origin` を月中心にするか太陽-月軸上の別基準にするかは ISSUE-020/021 の基本面構成と整合させる（要レビュー）。NASA のベッセル軸は月影軸＝太陽中心と月中心を結ぶ線。
- 本影半角と反本影は同一錐の表裏。`umbra_half_angle` 1つで両方を表し、頂点位置で本影/反本影領域を切り分ける（実装コメントで明記）。
- k の本影/半影使い分け（EspenakUmbral/EspenakPenumbral）は config 由来。既定 IauMean は単一値（conventions §9）。どのモデルで計算したか metadata に残す。
- `DegenerateGeometry` は呼出側（ISSUE-023）で「中心食でない部分食」等として扱う余地を残す（エラーで止めず可能性を潰さない）。設計はレビュー。
- レビュー重点: 本影=収束/半影=発散の符号、頂点距離の 0 除算回避、金環判定の地球面貫通計算、k 使い分けの metadata 記録。
