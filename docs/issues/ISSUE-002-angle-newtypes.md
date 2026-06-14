# ISSUE-002: Angle newtypes（Radians / Degrees、用途別正規化）

- crate: umbra-core
- 依存: ISSUE-001
- モード(tdd-workflow): strict（角度正規化は全幾何計算・接触/合 solver の連続性の前提。公開仕様かつ符号規約に直結するため strict）

## 目的
角度量の薄い newtype `Radians` / `Degrees` を提供し、**用途別の正規化関数**を分離して実装する。
- `normalized_signed`: `[-π, π)`（経度・時角など循環量。度では `[-180°, 180°)`）。
- `normalized_two_pi`: `[0, 2π)`（赤経・恒星時など）。
- `Radians ⇔ Degrees` の相互変換。

## 非目的
- 角度差の「連続関数化」（±π 折返し除去）そのもの。これは合・接触 solver 側（umbra-eclipse / algorithms）の責務。ここでは正規化関数の提供までで、混在禁止の型的土台を作る。
- 三角関数ラッパの提供（生 `f64` の `sin/cos` を内部で使うのは可）。

## 公開インターフェース
api-draft §1.1 を転記・具体化:
```rust
#[derive(Clone, Copy, Debug, PartialEq, PartialOrd)] pub struct Radians(pub f64);
#[derive(Clone, Copy, Debug, PartialEq, PartialOrd)] pub struct Degrees(pub f64);

impl Radians {
    pub fn to_degrees(self) -> Degrees;
    pub fn normalized_signed(self) -> Self;   // [-π, π)
    pub fn normalized_two_pi(self) -> Self;    // [0, 2π)
}
impl Degrees {
    pub fn to_radians(self) -> Radians;
    // 必要なら Degrees 側にも同名の正規化を提供（内部は Radians 経由 or 度で直接）
}
```
- 生 `f64` を「角度」として他 API へ渡さない（conventions §1）。

## 数式・アルゴリズムの出典
- 正規化は初等的（剰余演算）だが、IEEE754 で堅牢な実装として **`rem_euclid`** を基礎にする（標準ライブラリ仕様。負値でも非負剰余を返す）。
- 度↔ラジアン換算: `Meeus, "Astronomical Algorithms", 2nd ed.` の角度規約に準拠（π の定義は `std::f64::consts::PI`）。一般式のため特定章番号は不要だが、π 定数の出典として標準ライブラリを明記。
- `[-π,π)` 化は `((x + π).rem_euclid(2π)) - π` の標準手法（要確認: 端点 +π が `-π` に落ちる扱いを区間定義 `[-π, π)` と一致させる）。

## 単位 / 時刻系 / 座標系
- 入力単位: ラジアン（`Radians`）/ 度（`Degrees`）。
- 出力単位: 同上（型保存）。
- 時刻系 / 座標系: 無関係（純粋な角度ユーティリティ）。

## アルゴリズム概要
1. `to_degrees` / `to_radians`: `* 180/π` / `* π/180`。π は `std::f64::consts::PI`（magic number 禁止・conventions §11）。
2. `normalized_two_pi`: `x.rem_euclid(2π)`。結果は `[0, 2π)`。
3. `normalized_signed`: `let t = (x + π).rem_euclid(2π) - π;` 結果は `[-π, π)`。境界 `x == π` は `-π` 側へ寄せる（区間が右半開のため）。
- 数値安定性: 大きな入力角でも `rem_euclid` は安定。ただし丸めで端点ぎりぎりが反対端へ落ちうるので、区間の半開性をテストで固定。
- 禁止: 区間が曖昧な単一 `normalize()` を作って用途混在させること（conventions §2）。signed と two_pi を**別関数**に保つ。

## 受け入れテスト
accuracy.md テストレベル **L1（純数学）**。
- 厳密/準厳密一致テスト:
  - `Degrees(180).to_radians()` ≈ `π`（オラクル: 数学定数 π。許容 §下記）。
  - `Radians(π).to_degrees()` ≈ `Degrees(180)`。
  - `normalized_two_pi`: 入力 `-0.1, 0, 2π, 2π+0.1, -10π+0.3` → すべて `[0, 2π)` に入り、元角と `2π` の整数倍だけ異なる。
  - `normalized_signed`: 入力 `π, -π, 3π, -3π+0.001` → `[-π, π)`。特に `π → -π`、`-π → -π`（端点規約）を固定。
- 境界値: `0`, `±π`, `±2π`, 非常に大きい角（例 `1e6 * 2π`）。
- 異常系: `NaN` / `±inf` 入力時の挙動を明示（伝播させる＝そのまま `NaN` を返す）。テストで仕様固定。
- プロパティ（L8 候補）: 任意 x で `normalized_two_pi(x)` と x の差が `2π` の整数倍（許容内）。

## 許容誤差
accuracy.md にこの層の専用バジェット行はない。根拠付きで設定:
- 度↔ラジアン換算: 相対誤差 ≤ 数 ULP（理想は 1〜2 ULP）。根拠: 単純積で丸めは最小。実テストは絶対許容 `1e-12 rad` 程度を上限ガードに。
- 正規化の周期一致: 元角との差が `2π·n` から `≤ 1e-9 rad` 以内（大入力での丸め累積を考慮した実用上限。conventions §2 の「混在させない」を満たすこと自体は厳密＝関数分離は型で保証）。
- accuracy.md §2.1 で「solver 収束 0.05″」「角度バジェット 0.75″」を扱うため、この基盤層の誤差はそれより 3 桁以上小さく（≪0.001″）保つのが妥当。

## 実装メモ
- `PartialOrd` は単位混在比較を生まないよう、同型同士のみ。`Radians` と `Degrees` の直接比較は型で不可能（意図通り）。
- 端点規約（`[-π, π)` の右半開）はベッセル/時角計算の不連続点に効くため、テストで明示固定し、コメントに conventions §2 を参照。
- レビュー重点: `rem_euclid` の端点丸め（`π` が含まれるか）と、NaN 伝播仕様。`Eq`/`Hash` は f64 のため付けない（意図通り api-draft は付与していない）。
