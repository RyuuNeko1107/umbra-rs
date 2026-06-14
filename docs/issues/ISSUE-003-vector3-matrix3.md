# ISSUE-003: Vector3 / Matrix3（最小線形代数・回転行列）

- crate: umbra-core
- 依存: ISSUE-001, ISSUE-002（回転角に `Radians` を使う）
- モード(tdd-workflow): standard（数値の正確さは重要だが、初等線形代数で公開面も限定的。回転連鎖がフレーム変換の土台となるため trivial にはしない）

## 目的
日食パイプラインに必要な最小限の 3 次元線形代数を提供:
- `Vector3`（加減・スカラ倍・内積・外積・ノルム・正規化）。
- `UnitVector3`（正規化済み不変量）。
- `Matrix3`（行列積・行列×ベクトル・転置・恒等）と基本回転行列 `Rx / Ry / Rz`。
これらは歳差章動・ERA・極運動の回転連鎖（conventions §5）と影円錐/ベッセル基底（architecture §6）で使う。

## 非目的
- 一般の N 次元線形代数、固有値分解、外部 crate（nalgebra 等）への依存。
- フレーム型 `Position<F>` の変換ロジック本体（型は api-draft §1.4 にあるが、変換連鎖の実装は ephemeris/eclipse 側）。ここでは素の `Vector3`/`Matrix3` まで。
- 四元数（v1.0 では回転行列で統一。必要になれば別 issue）。

## 公開インターフェース
api-draft §1.4 を転記・具体化:
```rust
#[derive(Clone, Copy, Debug, PartialEq)] pub struct Vector3 { pub x: f64, pub y: f64, pub z: f64 }
#[derive(Clone, Copy, Debug)] pub struct UnitVector3(/* Vector3, 正規化済み */);

impl Vector3 {
    pub fn dot(self, o: Self) -> f64;
    pub fn cross(self, o: Self) -> Self;
    pub fn norm(self) -> f64;
    pub fn normalize(self) -> Result<UnitVector3, DomainError>; // 零ベクトルは Err
    // add/sub/scale（演算子 or メソッド）
}

// Matrix3 は内部利用中心。必要分のみ pub（api-draft §1.4 コメント準拠）
pub struct Matrix3 { /* [[f64;3];3] */ }
impl Matrix3 {
    pub fn identity() -> Self;
    pub fn rotation_x(theta: Radians) -> Self;
    pub fn rotation_y(theta: Radians) -> Self;
    pub fn rotation_z(theta: Radians) -> Self;
    pub fn transpose(self) -> Self;
    pub fn mul_mat(self, o: Self) -> Self;
    pub fn mul_vec(self, v: Vector3) -> Vector3;
}
```
- 距離量を持つベクトルは `Kilometers`/`Meters` ではなく素の `Vector3`（成分は呼び出し側が単位を管理。conventions §1: 単位変換は境界に集約）。**変数名で単位を明示**すること。

## 数式・アルゴリズムの出典
- 基本回転行列の定義: `IERS Conventions 2010 (IERS TN 36)` および `IAU SOFA` の `iauRx / iauRy / iauRz` と同一符号規約（右手系・能動/受動の別を実装コメントで固定）。SOFA の `Rx(φ)` は座標フレームを角度 φ 回転させる**受動回転**（フレーム回転）。本プロジェクトはフレーム連鎖（conventions §5）で使うため SOFA と同じ受動回転を採用する。
- 外積・内積・ノルム: 初等定義（`Meeus AA 2nd ed.` でも使用される標準ベクトル演算）。
- 要確認: 能動回転 vs 受動回転の符号。SOFA `iauRz` の行列形を正本としてコメントに転記し、テストで固定する。

## 単位 / 時刻系 / 座標系
- 入力: 回転角は `Radians`（ISSUE-002）。ベクトル成分は無次元 or 呼び出し側単位（km/AU/Re）。
- 出力: 同型。
- 座標系: 右手系（conventions §5）。回転は受動回転（フレーム軸の回転）で統一。

## アルゴリズム概要
1. `Vector3` 基本演算を素直に実装。`norm` は `hypot` 連鎖 or `sqrt(dot)`（オーバーフロー懸念が小さい日食スケールでは `sqrt(x²+y²+z²)` で可、ただしテストで桁落ち確認）。
2. `normalize`: ノルムが `0`（or 極小閾値）なら `DomainError`。閾値は magic number 化を避け定数化＋根拠コメント。
3. 回転行列 `Rz(θ)` 等は SOFA と同一の行列要素で定義。実装コメントに行列を明記。
4. `mul_mat` / `mul_vec` は素朴な三重/二重ループ（3×3 固定でアンロール可）。
- 数値安定性: 回転行列は直交。連鎖後の直交性ドリフトは小さいが、テストで `R·Rᵀ ≈ I` を確認。禁止: magic number の角度・正規化閾値直書き。

## 受け入れテスト
accuracy.md テストレベル **L1（純数学）**。
- 既知値（オラクル＝手計算/数学的恒等式、実装からコピーしない）:
  - `cross(x̂, ŷ) = ẑ`（右手系確認）。
  - `Rz(π/2)·x̂ = -ŷ`（受動回転の符号固定。能動だと `+ŷ`。どちらかをコメントの規約と一致させ厳密にテスト）。
  - `Rx(θ)·Rx(-θ) = I`。
  - `dot`, `norm`: ピタゴラス三つ組（3,4,0→5 など）で厳密一致。
- 直交性: 任意 θ で `R·Rᵀ` と `I` の各要素差 ≤ 許容。
- 境界/異常系: `normalize(0ベクトル)` → `Err(DomainError)`。極小ノルムの扱い。
- プロパティ（L8 候補）: `mul_vec(R, v)` のノルム保存（回転はノルム不変）。

## 許容誤差
accuracy.md に専用行なし。フレーム回転は §2.1 で「歳差章動+フレーム 0.05″（実力 ~1mas）」に含まれるため、本数値層はそれを劣化させない:
- 回転後ノルム保存 / 直交性 `R·Rᵀ=I`: 各要素 `≤ 1e-12`（根拠: 3×3 直交行列の丸めは数 ULP、十分に 1mas=4.8e-9 rad より小さく取る）。
- 内積・外積・ノルム: 相対誤差 数 ULP。ピタゴラス三つ組など整数系は厳密一致を期待。

## 実装メモ
- `UnitVector3` は構築時に正規化を保証する不変量。外部から生成分を直接いじれないようにする（フィールド非公開）。
- 距離単位はベクトルに型付けしない代わりに、呼び出し側で `// 単位: km`（or AU/Re）コメントを必須化（conventions §10）。`architecture §6` の `ShadowCone` 等が km/Re を混ぜないよう、境界で変換。
- レビュー重点: 受動/能動回転の符号が conventions §5 のフレーム連鎖方向と整合するか。SOFA 行列要素のコメント転記が正確か。
