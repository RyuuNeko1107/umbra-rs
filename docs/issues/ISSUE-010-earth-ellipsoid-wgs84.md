# ISSUE-010: Earth ellipsoid（WGS84・測地⇔地心緯度）

- crate: umbra-core
- 依存: ISSUE-001, ISSUE-002（角度）
- モード(tdd-workflow): strict（公開仕様＝`EarthModel`、観測者位置・中心線位置の基礎。測地/地心緯度の取り違えは誤差を直接生むため strict）

## 目的
WGS84 地球楕円体モデルと測地緯度⇔地心緯度の変換を提供する（conventions §4）。
- `EarthModel::Wgs84`（既定）。a = 6 378 137.0 m, 1/f = 298.257223563。
- 測地緯度 (geodetic) ⇔ 地心緯度 (geocentric) の変換（標高依存）。
- ベッセル無次元化の基準赤道半径 Re = WGS84 a を本モデルから供給。

## 非目的
- 観測者→ITRS/ECEF ベクトル化（ISSUE-011）。本 issue は楕円体パラメータと緯度変換まで。
- ジオイド高/正高の変換（conventions §4: 既定は楕円体高。正高入力変換は別途明示・対象外）。
- 大気差・重力モデル。

## 公開インターフェース
api-draft §1.5 / conventions §4 に整合:
```rust
#[non_exhaustive] #[derive(Clone, Copy, Debug)] pub enum EarthModel { Wgs84 }

impl EarthModel {
    pub fn equatorial_radius(self) -> Meters;     // a = 6378137.0
    pub fn inverse_flattening(self) -> f64;        // 1/f = 298.257223563
    pub fn flattening(self) -> f64;                // f
    pub fn polar_radius(self) -> Meters;           // b = a(1−f)
    pub fn eccentricity_squared(self) -> f64;      // e² = 2f − f²
    /// ベッセル無次元化基準（conventions §4: Re = WGS84 a）
    pub fn reference_equatorial_radius(self) -> Meters;
}

/// 測地緯度＋楕円体高 → 地心緯度（および地心距離成分）
pub fn geodetic_to_geocentric_lat(
    model: EarthModel,
    geodetic_lat: GeodeticLatitude,
    height: Meters,
) -> Radians;   // geocentric_lat（変数名で区別・conventions §3）
```
- 変数名で `geodetic_lat` / `geocentric_lat` を区別（conventions §3）。生 f64 角度を渡さない。

## 数式・アルゴリズムの出典
- **WGS84 定数**: NIMA TR8350.2 / NGA WGS84 定義。a = 6 378 137.0 m（厳密定義値）、1/f = 298.257223563（厳密定義値）。conventions §4 に固定済み。
- **測地⇔地心緯度・観測者の ρsinφ' / ρcosφ'**: **Meeus, "Astronomical Algorithms", 2nd ed., Chapter 11 "The Earth's Globe"**（式 11.1〜11.4: `u`, `ρ·sin φ'`, `ρ·cos φ'`, 標高 H を含む）。地心緯度 φ' = atan2(ρsinφ', ρcosφ')。
- e² = 2f − f²（標準楕円体関係。Meeus Ch.11）。
- 要確認: Meeus Ch.11 は b/a = 0.99664719（旧値）を例示。WGS84 の f から自前計算した値を使い、Meeus 例値は検証参照に留める（系統差をコメント）。

## 単位 / 時刻系 / 座標系
- 入力: 測地緯度 `GeodeticLatitude`（rad）、楕円体高 `Meters`。
- 出力: 地心緯度 `Radians`。半径は `Meters`、Re も `Meters`（ベッセル側で Re 無次元化に使用）。
- 時刻系: なし。
- 座標系: 楕円体（WGS84≈ITRS 軸、conventions §5）。本 issue は緯度変換まで（ベクトル化は ISSUE-011）。

## アルゴリズム概要
1. `EarthModel::Wgs84` の定数を関数で供給（a, 1/f は定義値、f/b/e² は導出）。すべて定数化＋出典コメント（magic number 禁止・conventions §11）。
2. `geodetic_to_geocentric_lat`: Meeus Ch.11 式で `u = atan((b/a)·tan φ)`、`ρsinφ' = (b/a)·sin u + (H/a)·sin φ`、`ρcosφ' = cos u + (H/a)·cos φ`。地心緯度 = `atan2(ρsinφ', ρcosφ')`。
3. 標高 H は楕円体高（conventions §4）。H=0 で純楕円体面の地心緯度。
- 数値安定性: 極（±90°）・赤道（0°）で `atan2` が安定。`tan φ` の発散を避け `atan2` で表現。禁止: 地心/測地の取り違え（変数名厳守）、b/a の magic 値直書き。

## 受け入れテスト
accuracy.md テストレベル **L1（純数学）** ＋ L6（局地条件）の前提。
- 既知値（オラクル＝Meeus AA 2nd ed. Ch.11 worked example、WGS84 定義値。実装からコピーしない）:
  - 赤道（φ=0）: 地心緯度 = 0。
  - 極（φ=±90°）: 地心緯度 = ±90°。
  - Meeus Ch.11 の Palomar 例（φ=33°21'22", H=1706 m）の `ρsinφ'`, `ρcosφ'` 近似一致（旧 b/a との系統差をコメント）。
  - 中緯度（φ=45°, H=0）: 地心緯度 < 測地緯度（差は最大 ~11.5′ 付近）。
- 定数: a, 1/f が WGS84 定義値と**厳密一致**。e²=2f−f², b=a(1−f) の導出一致。
- 境界値: φ=0, ±90°, 高標高（H=8000 m）、負の標高（−400 m, 死海等）。
- プロパティ（L8）: |地心緯度| ≤ |測地緯度|（H=0 で常に成立、符号同じ）。

## 許容誤差
accuracy.md §2.3「観測者/楕円体/標高 sub-m」「WGS84 で十分」から:
- 定数 a, 1/f: **厳密一致**（定義値）。
- 地心緯度変換: 角度で **≤ 1e-9 rad（≈0.2 mas）** 目標。根拠: §2.3 で観測者位置 sub-m を要求。地表で 1e-9 rad ≈ 6.4 mm（赤道半径×角）で sub-m に十分余裕。
- 地心距離成分（ρ）: 相対 ≤ 数 ULP。
- 根拠: accuracy.md §2.3 は観測者・楕円体・標高を sub-m とし WGS84 で十分とする。緯度変換誤差をそれより 2〜3 桁小さく保つ。

## 実装メモ
- conventions §4: Re（ベッセル無次元化基準）= WGS84 a を固定。NASA 慣習との差異（NASA は別の Re を使う場合あり）は accuracy.md に記録する旨コメント。
- 地心緯度は ISSUE-011（観測者→ITRS ベクトル）と umbra-eclipse（ベッセル面投影）で消費。`geodetic_lat`/`geocentric_lat` の命名規約（conventions §3）を厳守。
- `EarthModel` は `#[non_exhaustive]`（将来の楕円体追加に備える・api-draft §0）。
- レビュー重点: 測地/地心の取り違え無し、WGS84 定義値の厳密性、Meeus 旧 b/a 例値との系統差の扱い、極・赤道境界。
