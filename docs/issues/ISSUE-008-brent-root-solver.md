# ISSUE-008: Brent root solver（導関数不要・ブラケット必須）

- crate: umbra-core
- 依存: ISSUE-001
- モード(tdd-workflow): strict（合・接触のゼロ点求解の標準解法。収束保証と数値安定性が接触時刻精度を律速。安全性（無条件 Newton 禁止・conventions §11）に関わるため strict）

## 目的
導関数不要・ブラケット必須の **Brent 法**求根器を実装する（architecture §12, conventions §11）。
- 符号変化区間 `[a, b]`（`f(a)·f(b) < 0`）を必須入力とし、ブラケットが取れない場合はエラー。
- 収束基準（区間幅 or 関数値）を呼び出し側が指定（`root_tolerance_seconds` 等、目標の 1/10 以下・accuracy.md §2.1）。

## 非目的
- 1 次元最小化（ISSUE-009）。
- ブラケットの自動探索（粗走査での符号変化検出は合 solver 側 = umbra-eclipse / architecture §12）。本 issue はブラケット既知前提。
- Newton 法・割線法単独（conventions §11 で無条件 Newton 禁止。Brent 内部の逆二次補間＋二分法フォールバックのみ）。

## 公開インターフェース
architecture §12 / api-draft §1.6（`SolverError`）に整合:
```rust
pub struct RootConfig {
    pub x_tolerance: f64,       // 独立変数の許容（例: 日 or 秒。呼出側単位）
    pub max_iterations: u32,
}
pub enum Bracket { /* a, b, f(a), f(b) を保持し符号反転を保証 */ }

/// f はクロージャ。ブラケット必須。
pub fn brent_root<F: FnMut(f64) -> f64>(
    f: F,
    a: f64, b: f64,
    config: RootConfig,
) -> Result<f64, SolverError>;     // RootNotBracketed / DidNotConverge / NumericalInstability
```
- `SolverError`（api-draft §1.6）: `RootNotBracketed`, `DidNotConverge`, `NumericalInstability`。
- 独立変数が時刻のとき、呼び出し側は `JulianDate2` の差分（日/秒）を `f64` 引数に橋渡しする（型混在を避け、境界で変換。conventions §6）。

## 数式・アルゴリズムの出典
- **Brent (1973), "Algorithms for Minimization without Derivatives", Chapter 4（zeroin）**。原典アルゴリズム。
- 等価実装: **Numerical Recipes（Press et al.）の `zbrent`**（逆二次補間＋割線＋二分法のハイブリッド、ブラケット維持）。GSL `gsl_root_fsolver_brent` も同手順（参照のみ・コード移植しない）。
- 逆二次補間（inverse quadratic interpolation）と二分法フォールバックの切替条件は Brent 1973 / NR `zbrent` の判定式に従う（実装コメントに式を転記）。

## 単位 / 時刻系 / 座標系
- 入力/出力: 無次元 `f64`（呼び出し側が単位・時刻系を管理）。本 issue は純数値求根。
- 時刻系/座標系: なし。時刻のゼロ点求解に使う場合、独立変数の単位（日 or 秒）は呼び出し側が固定しコメント。

## アルゴリズム概要
1. 入力検証: `f(a)·f(b) < 0` でなければ `RootNotBracketed`（端点が 0 ちょうどなら即返す境界規約を固定）。
2. Brent 反復: 逆二次補間 → 妥当でなければ割線 → さらに不適なら二分法。各反復でブラケット `[a,b]` を縮小し**常に符号反転を維持**。
3. 収束: 区間幅 `|b−a| ≤ x_tolerance`（or 機械精度項を含む NR の許容式）で停止。
4. `max_iterations` 超過で `DidNotConverge`。NaN/inf 発生で `NumericalInstability`。
- 数値安定性: 許容式に機械イプシロン項を含める（NR `zbrent` の `tol1 = 2·eps·|b| + 0.5·tol`）。禁止: magic number の許容/反復上限直書き（`RootConfig` 経由）、無条件 Newton、ブラケット喪失。

## 受け入れテスト
accuracy.md テストレベル **L1（純数学）**。
- 既知根（オラクル＝解析解。実装からコピーしない）:
  - `f(x)=x²−2` on [1,2] → √2（厳密値と比較）。
  - `f(x)=cos(x)−x` on [0,1] → Dottie number（高精度参照値）。
  - `f(x)=x³−x−2` 等、逆二次補間が効く非線形。
- 収束特性: 反復回数が二分法より少ない（超線形）ことをカウントで確認。
- 境界値: 根が端点 `a` or `b`、ほぼ重根に近い緩い勾配、区間が極めて狭い。
- 異常系: `f(a)·f(b) > 0` → `RootNotBracketed`。`max_iterations` を 1 等に絞り `DidNotConverge`。`f` が NaN を返す → `NumericalInstability`。
- プロパティ（L8）: 返り値 `r` で `|f(r)|` が許容内、かつ `r ∈ [a,b]`。

## 許容誤差
accuracy.md §2.1「solver 収束 0.05″」「root_tolerance を目標の 1/10 以下」から:
- 時刻ゼロ点に使う場合の既定: 最大食 ±1.5s 目標 → root 収束は ≤ 0.15 s 相当（目標の 1/10）。角度感度 0.5″/s（§2.1）より 0.15 s ≈ 0.075″ ≪ 0.05″ バジェットに収まるよう、呼び出し側はさらに厳しく設定可。
- 純数学テストの根一致: 厳密解に対し `≤ 1e-10`（独立変数）程度。根拠: f64 で達成可能かつ時刻 µs 換算より十分小さい。
- 合格条件: 「許容を通すためだけに拡大しない」（conventions §11）。`RootConfig` で明示し、テストは目標の 1/10 ルールを検証。

## 実装メモ
- 端点が厳密に 0 の扱い（即返す/許容内とみなす）を固定しテスト。
- 時刻求解での独立変数は「合・接触」では ±π 折返しを除いた**連続関数**を呼び出し側が用意する前提（conventions §2 / architecture §12）。本 solver は連続 `f` を仮定する旨を doc に明記。
- レビュー重点: ブラケット維持（符号反転を絶対に失わない）、収束許容式の機械精度項、`DidNotConverge`/`NumericalInstability` の区別、無条件 Newton が混入していないこと。
