# ISSUE-009: One-dimensional minimizer（Brent 最小化 or 黄金分割）

- crate: umbra-core
- 依存: ISSUE-001
- モード(tdd-workflow): strict（最大食＝中心間距離関数の最小化に使う。最大食時刻精度を律速し、収束保証が必要なため strict。初期は黄金分割可だが公開仕様として厳格に）

## 目的
1 次元の単峰関数の最小値（と最小点）を求める最小化器を提供する（architecture §12, conventions §8 最大食）。
- 既定実装は **黄金分割探索**（初期・堅牢）。将来 **Brent 最小化**（放物線補間＋黄金分割）へ拡張可能な API。
- 3 点ブラケット `a < b < c` かつ `f(b) < f(a), f(b) < f(c)` を必須入力。

## 非目的
- 求根（ISSUE-008）。
- 多次元最適化。
- ブラケットの自動探索（粗走査で極小区間を検出するのは呼び出し側 = umbra-eclipse / architecture §12）。

## 公開インターフェース
architecture §12 / api-draft §1.6 に整合:
```rust
pub struct MinimizeConfig {
    pub x_tolerance: f64,       // 独立変数の許容（呼出側単位）
    pub max_iterations: u32,
}
pub struct Minimum { pub x: f64, pub f_x: f64 }

/// 3点ブラケット (a,b,c) 必須。単峰を仮定。
pub fn minimize_1d<F: FnMut(f64) -> f64>(
    f: F,
    a: f64, b: f64, c: f64,
    config: MinimizeConfig,
) -> Result<Minimum, SolverError>;   // RootNotBracketed（=ブラケット不正）/ DidNotConverge / NumericalInstability
```
- 時刻に使う場合、独立変数（日/秒）は呼び出し側が `JulianDate2` 差分から橋渡し（conventions §6、境界で変換）。

## 数式・アルゴリズムの出典
- **黄金分割探索**: 黄金比 `φ = (1+√5)/2`、縮小率 `1/φ ≈ 0.618`。`Numerical Recipes（Press et al.）"golden" / "Golden Section Search in One Dimension"`。Kiefer (1953) が原典。
- **Brent 最小化**（将来拡張）: **Brent (1973), "Algorithms for Minimization without Derivatives", Chapter 5**（放物線補間＋黄金分割のハイブリッド、`localmin`）。NR `brent`、GSL `gsl_min_fminimizer_brent`（参照のみ・移植しない）。
- 黄金比定数は `(1.0 + 5.0_f64.sqrt())/2.0`（magic number 化しない・conventions §11）。

## 単位 / 時刻系 / 座標系
- 入力/出力: 無次元 `f64`（呼び出し側が単位管理）。
- 時刻系/座標系: なし。最大食では「中心間距離（Re 単位など）の時刻関数」の最小化に使う（呼び出し側で単位コメント）。

## アルゴリズム概要
1. 入力検証: `a < b < c` かつ `f(b) < f(a)`, `f(b) < f(c)`（谷ブラケット）でなければ `RootNotBracketed`（=ブラケット不正）。
2. 黄金分割: 長い側に内分点を取り、関数値比較で区間を `1/φ` ずつ縮小。最小を内包し続ける。
3. 収束: `|c−a| ≤ x_tolerance·(|b|+eps)` 程度（NR の相対許容式）で停止し、最小点と値を返す。
4. `max_iterations` 超過で `DidNotConverge`。NaN/inf で `NumericalInstability`。
- 数値安定性: 最小値近傍は 2 次なので独立変数の許容は √eps オーダが理論下限（NR 注記）。これを下回る要求は警告/丸め扱い。禁止: magic number 許容、谷でない区間での暴走、無条件 Newton 的勾配法。

## 受け入れテスト
accuracy.md テストレベル **L1（純数学）**。
- 既知最小（オラクル＝解析解）:
  - `f(x)=(x−0.3)²` on (0,0.2,1) → 最小点 0.3、値 0。
  - `f(x)=−sin(x)` near π/2 → 最小点 π/2。
  - 中心間距離を模した擬似谷（放物線＋微小高次）で最小点復元。
- 収束: 区間が `1/φ` で縮小（黄金分割の理論縮小率）。
- 境界値: 最小がブラケット端近く、平坦に近い谷（緩勾配）、非常に狭い初期ブラケット。
- 異常系: `f(b)≥f(a)` 等の谷でないブラケット → エラー。`max_iterations` 過小 → `DidNotConverge`。NaN → `NumericalInstability`。
- プロパティ（L8）: 返り `x` で `f(x) ≤ f(a), f(b), f(c)`（許容内）。

## 許容誤差
accuracy.md §2.1（最大食時刻 ±1.5s、solver 収束 0.05″、root_tolerance 目標の 1/10）から:
- 最大食（中心間距離最小化）に使う場合: 独立変数の収束 ≤ 目標時刻の 1/10（≤ 0.15 s）。ただし 2 次最小のため独立変数許容は実効 √eps（≈1.5e-8 相対）が下限 — 距離関数値の最小化精度に注意（NR 注記: 最小近傍では x 精度が √eps に律速）。最大食「時刻」は距離が極小で平坦なため、必要なら最小点近傍で距離の対称性を使った別法をレビューで検討。
- 純数学テストの最小点一致: `≤ 1e-6`（独立変数。2 次最小の √eps 下限を踏まえた現実的値）。値一致は `≤ 1e-12`。
- 根拠: 黄金分割は線形収束（縮小率 0.618）。最大食の時刻精度 ±1.5s に対し十分な反復で到達可能。Brent 最小化へ移行で反復削減（accuracy.md L7 で profile 毎に採用判断）。

## 実装メモ
- 最大食の独立変数は連続時刻（折返し無し）。距離関数の極小近傍は平坦なので、x 許容を過度に小さくしても値が改善しない点を doc 明記（NR 注記の √eps 下限）。
- 初期は黄金分割で可（architecture §12）。Brent 最小化拡張時も同 API を保つ（`SolverError` 共通）。
- レビュー重点: 谷ブラケット検証、収束許容の √eps 下限の扱い、最大食での平坦谷への対処方針、無条件勾配法の不在。
