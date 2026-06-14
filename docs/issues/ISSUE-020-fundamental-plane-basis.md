# ISSUE-020: Fundamental plane basis（z=月影軸・y=天球北射影・x=右手系東）

- crate: umbra-eclipse
- 依存: ISSUE-019（影軸方向）, ISSUE-015（見かけ地心位置）, ISSUE-003（Vector3/UnitVector3/Matrix3）, ISSUE-002（Radians）, umbra-core（Position<FundamentalPlane>, ReferenceFrame）
- モード(tdd-workflow): strict（基本面の基底はベッセル要素 x,y,d,μ（ISSUE-021）の座標系定義そのもの。基底の符号・向きを誤ると全要素が誤る。極端配置で基底が壊れない数値安定性が必須。strict）

## 目的
月影軸を z 軸とする**ベッセル基本面（fundamental plane）**の正規直交基底を構成する（architecture §6, conventions §5 FundamentalPlane）。
- **z = 月影軸方向**（地心→影軸, 太陽方向に一致する向き規約）。
- **y = 天球北を基本面へ射影した方向**。
- **x = 右手系で東向き**（x = y × z 相当, 右手系統一, conventions §5）。
- 極端配置（z が天の北極に極めて近い等）で基底が縮退しないようにする（architecture §6）。

## 非目的
- 影円錐の半角・頂点（ISSUE-019）。本 Issue は基底（回転行列）構成のみ。
- ベッセル要素 x,y,d,μ,l1,l2 の数値算出（ISSUE-021。本 Issue の基底＝ISSUE-021 が投影に使う座標系）。
- 観測者の基本面座標への投影（局地, 別 Issue。本 Issue は基底定義まで）。

## 公開インターフェース
conventions §5 / api-draft §1.4（`FundamentalPlane`, `Position<F>`, `Matrix3` 内部）に整合。crate 内部中心。

```rust
/// 基本面の正規直交基底（地心赤道系 → 基本面系 の回転）。
pub(crate) struct FundamentalPlaneBasis {
    pub x_axis: UnitVector3,  // 東（右手系）
    pub y_axis: UnitVector3,  // 天球北の射影
    pub z_axis: UnitVector3,  // 月影軸（地心→太陽向き）
    pub rotation: Matrix3,    // 赤道系ベクトル → 基本面系ベクトル（行 = 各軸）
    pub declination: Radians, // d = z 軸の赤緯（ISSUE-021 の d）
    pub hour_angle_ref: Radians, // μ 構成用: z 軸の赤経 → グリニッジ恒星時と組む（ISSUE-021）
}

pub(crate) fn fundamental_plane_basis(
    shadow_axis_direction: UnitVector3,  // ISSUE-019: 地心→太陽向きの影軸
) -> Result<FundamentalPlaneBasis, EclipseError>;  // DegenerateGeometry（極端配置）
```

- `Position<FundamentalPlane>`（api-draft §1.4）への変換は `rotation` 適用で行う（型でフレーム判別, conventions §5）。
- `declination` (d) と `hour_angle_ref` は ISSUE-021 が d, μ を組み立てる素。

## 数式・アルゴリズムの出典
- **基本面の定義**: Explanatory Supplement to the Astronomical Almanac (3rd ed.), **Ch.11 Besselian elements 節**。基本面 = 地球中心を通り月影軸に垂直な平面。z 軸 = 影軸（月→太陽の反対の伸長方向だが、ベッセル慣習では地心から見て太陽方向の単位ベクトル方向を z にとる）。
- **軸の構成（標準的なベッセル基底）**: Meeus, *Astronomical Algorithms* (2nd ed.), **Ch.54「Eclipses」**。z 軸の赤経 α・赤緯 d（= 影軸の方向）を求め、
  - z = (cos d cos α, cos d sin α, sin d)（赤道直交系）。
  - y = 天の北極ベクトル ẑ_eq=(0,0,1) を z に直交する成分へ射影し正規化（北向き）。
  - x = y × z（右手系・東向き）。
  - 出典の符号規約（x=東/y=北/z=太陽向き）を実装コメントに明記（conventions §10）。NASA Espenak の x,y は基本面上の影中心座標で、この x(東)/y(北) 基底と一致。
- **μ（エフェメリス時角/赤経基準）**: μ = グリニッジ見かけ恒星時 − α_z（影軸赤経）。NASA 表記の μ（基本面の x 軸が地球経度に対し回る量）に対応（ISSUE-021 で完成）。本 Issue は α_z（hour_angle_ref の素）を提供。
- **d** = 影軸の赤緯（= z 軸の赤緯）。NASA 表記の d と一致。

## 単位 / 時刻系 / 座標系
- 座標系: 入力 = 地心赤道（CIRS/of date, ISSUE-015/035 連鎖）の影軸単位ベクトル。出力 = 基本面基底（FundamentalPlane フレーム, conventions §5）。
- 時刻系: 基底自体は瞬時量で TT 基準（conventions §6, ベッセル要素は TT）。μ の恒星時部分は UT1 由来（ISSUE-021 で結合）。本 Issue の d/α_z は TT 時刻の幾何。
- 単位: 角度ラジアン。基底は無次元単位ベクトル。
- 右手系統一（conventions §5）。x=東, y=北, z=太陽向き。

## アルゴリズム概要
1. 影軸方向（地心→太陽向き単位ベクトル, ISSUE-019）から赤経 α_z = atan2(z_y, z_x)、赤緯 d = asin(z_z) を算出（赤道直交系前提）。
2. z 軸 = 影軸単位ベクトル。
3. y 軸 = 北極ベクトル (0,0,1) の z 直交成分 `n − (n·z)z` を正規化。
4. x 軸 = y × z（右手系・東向き）を正規化。直交性・右手性を検証（`x·y≈0`, `det(rotation)≈+1`）。
5. `rotation` = 各軸を行に持つ行列（赤道→基本面）。
6. d, α_z を保持（ISSUE-021 の d, μ 構成へ）。
- **極端配置の数値安定性（最重要）**: z が (0,0,±1)（天の北/南極）に極めて近いと `n − (n·z)z ≈ 0` で y が定義不能。この縮退を検出し、代替射影（例: 別の参照軸 x_eq=(1,0,0) を使う）or `DegenerateGeometry` を返す。日食の影軸が天の極に来ることは現実にはほぼ無い（太陽は黄道上 |δ|≤23.5°）ので、防御的処理＋テストで担保（architecture §6「極端配置で基底が壊れない」）。asin 引数クランプ（accuracy.md §2.2）。

## 受け入れテスト
accuracy.md テストレベル **L4（ベッセル）** ＋ L1（線形代数）。
- 直交性・右手性（L1）: 任意の影軸方向で `|x|=|y|=|z|=1`、相互直交（内積≈0）、`det(rotation)=+1`（右手系）。プロパティテスト（L8）で多数のランダム z。
- 符号規約テスト: x が東（赤経増加方向）、y が北（赤緯増加方向）を指すことを既知方向で確認（オラクル＝幾何の手計算, 実装非コピー）。
- d/α_z 検証: 既知の影軸（例: 太陽が春分点方向）で d≈0, α_z 既知。
- **極端配置テスト（必須, architecture §6）**: z を (0,0,1) に漸近させ、基底が NaN/縮退せず安定 or `DegenerateGeometry` を返すこと。現実の日食範囲（|d|≤23.5°）では常に正常。
- **MockEphemeris 整合**: 人工配置の影軸から構成した基底が、ISSUE-021 で NASA 公開 d/μ と整合する素を与えること（ISSUE-021 と連結テスト）。

## 許容誤差
- accuracy.md §2.1 幾何バジェット（影幾何誤差, §4 層分解）の一部。基底の誤差は d, x, y（ISSUE-021）経由で最大食時刻・中心線に効く（中心線 sub-km ≲0.5km, accuracy.md §2.1）。
- 基底構成は純幾何で f64 機械精度（≪ バジェット）。誤差は入力影軸方向（ISSUE-019/015）律速。
- 直交性残差は `< 1e-12`（f64 で達成可能, 中心線 sub-km へ余裕）。

## 実装メモ
- z 軸の向き規約（地心→太陽 か 太陽→地心 か）を ISSUE-019/021/NASA と厳密に揃える。**ベッセル慣習では z は地心から太陽（影軸）方向**。誤ると x,y の符号が反転し全要素が鏡像になる。最重要レビュー項目。
- y の射影に使う「北」は天の北極（赤道系 (0,0,1)）。黄道北ではない（NASA 基本面は赤道基準）。混同禁止。
- 極端配置の代替射影は現実には発火しないが、fuzz/property（L8/L9）で踏むため必ず実装。
- μ は本 Issue では未完（恒星時結合は ISSUE-021）。α_z のみ供給し責務を分離。
- レビュー重点: z 向き規約、x=東/y=北の符号、右手系 det=+1、極端配置の縮退防御、北=天の北極（黄道北でない）。
