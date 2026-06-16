# Mutation review — `umbra-eclipse::bessel_poly` / `polynomial` / `source`（ISSUE-022 / 037）

`cargo mutants` の生存変異の列挙と許容判断。退行ガード `mutation.yml` は eclipse を含む。

## `source.rs`（ISSUE-037, BesselianSource trait + DirectBesselianSource）

- mutants: 5、**全 unviable**（viable 変異ゼロ）。`at()`/`fit_interval()`/`new()` は純粋なパススルーで
  返り値置換のみが候補だが、`InstantaneousBesselianElements`/`TimeInterval` が `Default` 非実装のため
  コンパイル不能 → unviable。kill すべき survivor なし（`apparent` 純合成と同 category）。

## `polynomial.rs`（ISSUE-022 S1, Polynomial + Chebyshev fit カーネル）

- 最終: **102 mutants / 100 caught / 2 unviable / 0 missed**。
- 初回 2 missed（`v * t_prev` / `cheb[0] * t_prev[0]` の `* → /`）は被演算子が定数 `1.0`（T₀ 値/係数）で
  `×1.0 ≡ ÷1.0` の等価変異 → **冗長な `*1.0` を除去して変異自体を消滅**（survivor ゼロ化）。

## `bessel_poly.rs`（ISSUE-022 S2, BesselianPolynomial::fit / at / unwrap_mu）

- 最終: **83 mutants / 74 caught / 7 unviable / 2 missed（等価, 下記）**。
- 初回 29 missed → テスト強化で 2 へ。強化の要点:
  - **unwrap_mu を branchless（`round(diff/TAU)` 一括補正）へ書換**＋**μ 減少（逆方向 wrap）テスト追加** →
    while ループ/recompute の生存（片方向のみ励起・1 回補正で recompute 値が無関係）を一掃。
  - **非多項式（三角）源テスト追加**（対称＋非対称区間）→ `t_center` / `half_width` / 残差サンプル span の
    幾何変異を kill。多項式源はサンプリング配置に不変（自己整合写像）かつ実暦は滑らかすぎて検出できないため、
    曲率のある sin/cos 源を締めた許容（1e-3）で張り、配置を崩すと外挿で残差ゲートが落ちるよう構成。
  - **全成分要求テスト**（`within` の `&&`）／**非有限区間テスト**（区間検査 `||` の各項）／
    **degree=0 定数フィットテスト**（`m = deg + 3` が 0 ノードにならない: 0 ノードは NaN 多項式を生み
    `f64::max` が NaN を無視して残差 0 と誤報し素通りする経路）を追加。

### 生存（許容）= 2 件・等価変異

`at()` の区間端チェック `if t < t0 - EPS || t > t1 + EPS`（`INTERVAL_EPS_HOURS = 1e-6` hour）:

| 変異 | 理由（等価） |
|---|---|
| `bessel_poly.rs:243 < → <=` | 差が出るのは `t == t0 − EPS` ちょうどの一点のみ。EPS 自体が端点を含めるための数値マージン（任意のtolerance）であり、その境界点を含む/含まないは契約上無差別。FP で当該点を厳密に踏むテストは構成不能・無意味。 |
| `bessel_poly.rs:243 > → >=` | 同上（`t == t1 + EPS` の一点）。 |

`delta_t_seconds` の区分境界 `< → <=`（接合点連続で差 ≪ 不確実性, `mutation-umbra-core.md`）と同 category の
**境界等価変異**。退行ガード `mutation.yml` で `--exclude-re 'with <= in.*BesselianPolynomial>::at'` /
`'with >= in.*BesselianPolynomial>::at'` を除外し、この 2 件を除いて eclipse は **0 missed** を要求する。
