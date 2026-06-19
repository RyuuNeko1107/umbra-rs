//! 汎用多項式型と Chebyshev 最小二乗フィット数値カーネル（ISSUE-022 スライス1）。
//!
//! 日食・暦に依存しない純粋な数値計算:
//! - [`Polynomial`]: 単項式（monomial）係数による多項式（昇べき, Horner 評価, 解析微分）。
//! - [`chebyshev_nodes`]: 第一種 Chebyshev ノード。
//! - [`fit_chebyshev_monomial`]: Chebyshev 基底での最小二乗フィット → monomial(t) へ基底変換。
//!
//! フィットは Chebyshev 基底（τ∈[-1,1]）で行い、出力時に monomial(t) へ基底変換する
//! （単項式 Vandermonde 直接フィットの悪条件を避ける, numerical-policy §A4）。

// 本モジュールの整数→浮動小数変換は小さなループ添字・多項式次数のみ（degree は数次、ノード数も小）で
// f64 に厳密表現できる範囲。精度クリティカルな天文量の変換ではないため、添字算術に限り許容する。
#![allow(clippy::cast_precision_loss, clippy::cast_possible_truncation)]

use core::f64::consts::PI;

/// 単項式（monomial）係数による多項式。係数は昇べき: `p(t) = Σ coefficients[i] · tⁱ`。
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize))]
pub struct Polynomial {
    /// 昇べき順の係数（`coefficients[i]` が tⁱ の係数）。
    pub coefficients: Vec<f64>,
}

impl Polynomial {
    /// Horner 法で `p(t)` を評価する。空係数は 0。
    pub fn eval(&self, t: f64) -> f64 {
        let mut acc = 0.0;
        for &c in self.coefficients.iter().rev() {
            acc = acc * t + c;
        }
        acc
    }

    /// 導関数 `p'(t)`（変数 t に関する微分）。
    /// `p(t)=Σ aᵢ tⁱ ⇒ p'(t)=Σ_{i≥1} i·aᵢ t^{i-1}`。定数（または空）の導関数は `[0]`。
    pub fn derivative(&self) -> Polynomial {
        if self.coefficients.len() <= 1 {
            return Polynomial {
                coefficients: vec![0.0],
            };
        }
        let coefficients = self
            .coefficients
            .iter()
            .enumerate()
            .skip(1)
            .map(|(i, &c)| (i as f64) * c)
            .collect();
        Polynomial { coefficients }
    }
}

/// M 個の第一種 Chebyshev ノード `τ_j = cos(π (j + 1/2) / M)`（j = 0..M-1）。
/// 端点 ±1 は含まず、各値は開区間 (-1, 1)。
pub(crate) fn chebyshev_nodes(m: usize) -> Vec<f64> {
    (0..m)
        .map(|j| (PI * (j as f64 + 0.5) / m as f64).cos())
        .collect()
}

/// 値列 `(node_taus[j], node_values[j])` を `degree` 次の monomial Polynomial（変数 t）へ
/// Chebyshev 最小二乗フィットする。τ = (t − t_center) / half_width。
///
/// Chebyshev ノード上の離散直交性で係数を求め（正規方程式不要・好条件）、Chebyshev(τ) →
/// power(τ) → τ=(t−t_center)/half_width 代入で monomial(t) へ基底変換する（numerical-policy §A4）。
pub(crate) fn fit_chebyshev_monomial(
    node_taus: &[f64],
    node_values: &[f64],
    degree: usize,
    t_center: f64,
    half_width: f64,
) -> Polynomial {
    let m = node_taus.len();
    // 離散 Chebyshev 変換: c_0 = (1/M) Σ vⱼ T_0(τⱼ), c_n = (2/M) Σ vⱼ T_n(τⱼ) (n≥1)。
    let mut cheb = vec![0.0; degree + 1];
    for (&tau, &v) in node_taus.iter().zip(node_values.iter()) {
        cheb[0] += v; // v · T_0(τ) = v（T_0 = 1）
        if degree >= 1 {
            let mut t_prev = 1.0; // T_0(τ)
            let mut t_cur = tau; // T_1(τ)
            cheb[1] += v * t_cur;
            for c in cheb.iter_mut().skip(2) {
                let t_next = 2.0 * tau * t_cur - t_prev; // T_{n} = 2τ T_{n-1} − T_{n-2}
                *c += v * t_next;
                t_prev = t_cur;
                t_cur = t_next;
            }
        }
    }
    let inv_m = 1.0 / m as f64;
    for (n, c) in cheb.iter_mut().enumerate() {
        *c *= if n == 0 { inv_m } else { 2.0 * inv_m };
    }

    let power_tau = chebyshev_to_power(&cheb);
    let power_t = substitute_affine(&power_tau, t_center, half_width);
    Polynomial {
        coefficients: power_t,
    }
}

/// Chebyshev 係数（τ∈[-1,1]）を τ の単項式係数へ変換する（`T_n = 2τ T_{n-1} − T_{n-2}` で展開）。
fn chebyshev_to_power(cheb: &[f64]) -> Vec<f64> {
    let degree = cheb.len() - 1;
    let mut out = vec![0.0; degree + 1];
    out[0] += cheb[0]; // cheb[0] · T_0(τ)（T_0 = [1]）
    if degree >= 1 {
        let mut t_prev = vec![1.0]; // T_0(τ) の単項式係数
        let mut t_cur = vec![0.0, 1.0]; // T_1(τ)
        for (i, &a) in t_cur.iter().enumerate() {
            out[i] += cheb[1] * a;
        }
        for &ck in cheb.iter().skip(2) {
            // T_next = 2τ·T_cur − T_prev
            let mut t_next = vec![0.0; t_cur.len() + 1];
            for (i, &a) in t_cur.iter().enumerate() {
                t_next[i + 1] += 2.0 * a;
            }
            for (i, &a) in t_prev.iter().enumerate() {
                t_next[i] -= a;
            }
            for (i, &a) in t_next.iter().enumerate() {
                out[i] += ck * a;
            }
            t_prev = t_cur;
            t_cur = t_next;
        }
    }
    out
}

/// τ の単項式多項式に τ=(t − t_center)/half_width を代入し、t の単項式係数を返す。
/// `(t − t_center)^k = Σ_{i=0..k} C(k,i) tⁱ (−t_center)^{k−i}`。
fn substitute_affine(power_tau: &[f64], t_center: f64, half_width: f64) -> Vec<f64> {
    let degree = power_tau.len() - 1;
    let mut out = vec![0.0; degree + 1];
    for (k, &ak) in power_tau.iter().enumerate() {
        if ak == 0.0 {
            continue;
        }
        let coef = ak / half_width.powi(k as i32);
        // out の長さは degree+1 ≥ k+1。i は添字かつ二項展開の引数として使う。
        for (i, o) in out.iter_mut().take(k + 1).enumerate() {
            *o += coef * binomial(k, i) * (-t_center).powi((k - i) as i32);
        }
    }
    out
}

/// 二項係数 C(n, k)（小さい n のみ。乗算的に計算）。
fn binomial(n: usize, k: usize) -> f64 {
    let mut result = 1.0;
    for i in 0..k {
        result *= (n - i) as f64 / (i + 1) as f64;
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    // ----------------------------------------------------------------------
    // 浮動小数比較ヘルパ
    // ----------------------------------------------------------------------

    /// 絶対誤差比較。
    fn close_abs(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    /// 相対 or 絶対の緩い比較（大きい値での丸めを許容）。
    fn close_rel(a: f64, b: f64, tol: f64) -> bool {
        let scale = 1.0_f64.max(a.abs()).max(b.abs());
        (a - b).abs() < tol * scale
    }

    // ----------------------------------------------------------------------
    // 自前オラクル（実装側の関数は使わない）
    // ----------------------------------------------------------------------

    /// 素朴なべき乗和による多項式評価（Horner のクロスチェック用オラクル）。
    /// p(t) = Σ coeffs[i] · t^i。
    fn naive_poly_eval(coeffs: &[f64], t: f64) -> f64 {
        let mut acc = 0.0;
        for (i, &c) in coeffs.iter().enumerate() {
            acc += c * t.powi(i as i32);
        }
        acc
    }

    /// 中心差分による数値微分（解析微分のクロスチェック用オラクル）。
    fn central_diff(f: impl Fn(f64) -> f64, t: f64, h: f64) -> f64 {
        (f(t + h) - f(t - h)) / (2.0 * h)
    }

    /// 既知の monomial 多項式を素朴に評価するクロージャを返す。
    fn known_poly(coeffs: Vec<f64>) -> impl Fn(f64) -> f64 {
        move |t: f64| naive_poly_eval(&coeffs, t)
    }

    // ======================================================================
    // 契約1: Polynomial::eval は Horner で正しい多項式値
    // ======================================================================

    /// eval が素朴べき乗和に一致する（定数・1次・高次, 複数 t）。
    #[test]
    fn eval_matches_naive_power_sum() {
        let cases: Vec<Vec<f64>> = vec![
            vec![7.0],                             // 定数
            vec![0.0],                             // ゼロ定数
            vec![-2.0, 3.0],                       // 1次
            vec![1.5, 0.0, -0.7],                  // 2次（中間係数 0）
            vec![2.0, -3.0, 5.0, 4.0],             // 3次
            vec![0.1, -0.2, 0.3, -0.4, 0.5, -0.6], // 5次
        ];
        let ts = [-3.7, -1.0, -0.25, 0.0, 0.5, 1.0, 2.3, 10.0];
        for coeffs in &cases {
            let p = Polynomial {
                coefficients: coeffs.clone(),
            };
            for &t in &ts {
                let expected = naive_poly_eval(coeffs, t);
                // 相対許容（高次・大 t での丸めを吸収）。
                assert!(
                    close_rel(p.eval(t), expected, 1e-9),
                    "eval mismatch: coeffs={:?} t={} got={} want={}",
                    coeffs,
                    t,
                    p.eval(t),
                    expected
                );
            }
        }
    }

    /// 空係数 `coefficients=[]` は任意の t で 0 を返す（境界）。
    #[test]
    fn eval_empty_coefficients_is_zero() {
        let p = Polynomial {
            coefficients: vec![],
        };
        for &t in &[-5.0, 0.0, 1.0, 42.0] {
            assert!(
                close_abs(p.eval(t), 0.0, 1e-15),
                "empty poly must be 0 at t={t}"
            );
        }
    }

    // ======================================================================
    // 契約2: Polynomial::derivative は解析微分
    // ======================================================================

    /// 既知ケース: p = 2 -3t +5t² +4t³ の導関数は -3 +10t +12t²。
    /// 係数 [2,-3,5,4] -> [-3,10,12]（手計算済みの小さな既知ケース）。
    #[test]
    fn derivative_known_small_case() {
        let p = Polynomial {
            coefficients: vec![2.0, -3.0, 5.0, 4.0],
        };
        let d = p.derivative();
        let expected = [-3.0, 10.0, 12.0];
        assert_eq!(
            d.coefficients.len(),
            expected.len(),
            "derivative degree should drop by 1: got {:?}",
            d.coefficients
        );
        for (i, &e) in expected.iter().enumerate() {
            assert!(
                close_abs(d.coefficients[i], e, 1e-12),
                "derivative coeff[{i}] got={} want={e}",
                d.coefficients[i]
            );
        }
    }

    /// 一般則 i·a_i をテスト内で再計算してクロスチェック（マジック値を避ける）。
    #[test]
    fn derivative_coefficients_follow_i_times_a_i() {
        let coeffs = vec![1.1, -2.2, 3.3, -4.4, 5.5];
        let p = Polynomial {
            coefficients: coeffs.clone(),
        };
        let d = p.derivative();
        // p'(t) = Σ_{i>=1} i·a_i t^{i-1}
        let expected: Vec<f64> = (1..coeffs.len()).map(|i| (i as f64) * coeffs[i]).collect();
        assert_eq!(d.coefficients.len(), expected.len());
        for (i, &e) in expected.iter().enumerate() {
            assert!(
                close_abs(d.coefficients[i], e, 1e-12),
                "deriv coeff[{i}] got={} want={e}",
                d.coefficients[i]
            );
        }
    }

    /// 「derivative の eval」＝「元 eval の中心差分」で独立クロスチェック（複数 t, 数値微分 tol）。
    #[test]
    fn derivative_eval_matches_central_difference() {
        let p = Polynomial {
            coefficients: vec![0.5, -1.3, 2.7, -0.9, 0.4],
        };
        let d = p.derivative();
        let h = 1e-4;
        for &t in &[-2.0, -0.5, 0.0, 0.7, 1.5, 3.0] {
            let analytic = d.eval(t);
            let numeric = central_diff(|x| p.eval(x), t, h);
            assert!(
                close_rel(analytic, numeric, 1e-5),
                "derivative eval vs central diff: t={t} analytic={analytic} numeric={numeric}"
            );
        }
    }

    /// 定数多項式の導関数は 0（係数 [0] か空; eval が 0 を返せばよい）。
    #[test]
    fn derivative_of_constant_is_zero() {
        for c in [&vec![5.0][..], &vec![][..]] {
            let p = Polynomial {
                coefficients: c.to_vec(),
            };
            let d = p.derivative();
            for &t in &[-1.0, 0.0, 3.0] {
                assert!(
                    close_abs(d.eval(t), 0.0, 1e-12),
                    "derivative of constant {:?} must eval 0 at t={t}, got {}",
                    c,
                    d.eval(t)
                );
            }
        }
    }

    // ======================================================================
    // 契約3: chebyshev_nodes
    // ======================================================================

    /// 長さは M, 各値は (-1,1) に厳密に入る, ノードは相異なる。
    #[test]
    fn chebyshev_nodes_length_range_and_distinct() {
        for m in [1usize, 2, 3, 5, 8, 16] {
            let nodes = chebyshev_nodes(m);
            assert_eq!(nodes.len(), m, "length must equal M={m}");
            for &v in &nodes {
                assert!(v > -1.0 && v < 1.0, "node {v} must be strictly in (-1,1)");
            }
            // 相異なる（端点除外の cos なので厳密に distinct）。
            for i in 0..nodes.len() {
                for j in (i + 1)..nodes.len() {
                    assert!(
                        (nodes[i] - nodes[j]).abs() > 1e-12,
                        "nodes must be distinct: [{i}]={} [{j}]={}",
                        nodes[i],
                        nodes[j]
                    );
                }
            }
        }
    }

    /// M=1 のとき τ_0 = cos(π/2) = 0。
    #[test]
    fn chebyshev_nodes_single_is_zero() {
        let nodes = chebyshev_nodes(1);
        assert_eq!(nodes.len(), 1);
        assert!(
            close_abs(nodes[0], 0.0, 1e-12),
            "M=1 node must be 0, got {}",
            nodes[0]
        );
    }

    /// M=2 で {cos(π/4), cos(3π/4)} = {+√2/2, −√2/2}（順序非依存で集合一致）。
    #[test]
    fn chebyshev_nodes_two_known_values() {
        let nodes = chebyshev_nodes(2);
        assert_eq!(nodes.len(), 2);
        let want_pos = (std::f64::consts::PI / 4.0).cos(); // +√2/2
        let want_neg = (3.0 * std::f64::consts::PI / 4.0).cos(); // -√2/2
                                                                 // 順序を仮定せず、両期待値がそれぞれ一致するノードを持つことを確認。
        let has = |target: f64| nodes.iter().any(|&n| close_abs(n, target, 1e-12));
        assert!(has(want_pos), "missing node ≈ {want_pos}, got {nodes:?}");
        assert!(has(want_neg), "missing node ≈ {want_neg}, got {nodes:?}");
    }

    // ======================================================================
    // 契約4: fit_chebyshev_monomial の往復（最重要）
    // ======================================================================

    /// 既知 monomial 多項式 q を t_j = t_center + half_width·τ_j で評価し、
    /// フィット結果が区間内で q を復元する（eval 残差を主オラクル）。
    /// degree == deg(q), 複数の (t_center, half_width) を網羅。
    fn assert_fit_recovers(q_coeffs: &[f64], degree: usize, t_center: f64, half_width: f64) {
        let q = known_poly(q_coeffs.to_vec());
        // M >= degree+1。余裕を持って degree+3。
        let m = degree + 3;
        let taus = chebyshev_nodes(m);
        let values: Vec<f64> = taus
            .iter()
            .map(|&tau| q(t_center + half_width * tau))
            .collect();
        let fit = fit_chebyshev_monomial(&taus, &values, degree, t_center, half_width);

        // 出力は長さ degree+1 の monomial(t)。
        assert_eq!(
            fit.coefficients.len(),
            degree + 1,
            "fit output length must be degree+1={}",
            degree + 1
        );

        // 主オラクル: 区間 [t_center-half_width, t_center+half_width] の多数の t で eval 一致。
        let n = 50;
        let mut max_resid = 0.0_f64;
        for k in 0..=n {
            let frac = (k as f64) / (n as f64); // 0..1
            let t = (t_center - half_width) + 2.0 * half_width * frac;
            let resid = (fit.eval(t) - q(t)).abs();
            max_resid = max_resid.max(resid);
        }
        assert!(
            max_resid < 1e-6,
            "fit residual too large: q={:?} degree={degree} t_center={t_center} half_width={half_width} max_resid={max_resid}",
            q_coeffs
        );
    }

    /// 次数ちょうど（degree == deg(q)）での完全復元, t_center=0 と t_center≠0 の両方,
    /// half_width=1 と half_width≠1 の両方。
    #[test]
    fn fit_exact_degree_round_trip() {
        // (q_coeffs, degree)
        let polys: Vec<(Vec<f64>, usize)> = vec![
            (vec![1.5], 0),                        // 定数
            (vec![1.5, 2.0], 1),                   // 1次
            (vec![1.5, 2.0, -0.7], 2),             // 2次
            (vec![0.3, -1.1, 0.4, 0.9], 3),        // 3次
            (vec![-0.5, 0.2, -0.3, 0.4, -0.1], 4), // 4次
        ];
        // 複数の (t_center, half_width): クリーンと非クリーン, half_width≠1 を必ず含む。
        let frames = [(0.0, 1.0), (0.0, 2.5), (3.0, 2.5), (-4.2, 0.75)];
        for (q, degree) in &polys {
            for &(tc, hw) in &frames {
                assert_fit_recovers(q, *degree, tc, hw);
            }
        }
    }

    /// 過剰次数: degree が q の真の次数より高くても eval が q に一致（高次係数 ~0）。
    #[test]
    fn fit_overdetermined_degree_round_trip() {
        // q は 1次、degree=3 でフィット。
        let q = vec![2.0, -0.5];
        let degree = 3;
        let frames = [(0.0, 1.0), (3.0, 2.5)];
        for &(tc, hw) in &frames {
            assert_fit_recovers(&q, degree, tc, hw);
            // 高次係数（2次・3次）が ~0 であることも確認。
            let m = degree + 3;
            let taus = chebyshev_nodes(m);
            let qf = known_poly(q.clone());
            let values: Vec<f64> = taus.iter().map(|&tau| qf(tc + hw * tau)).collect();
            let fit = fit_chebyshev_monomial(&taus, &values, degree, tc, hw);
            assert!(
                close_abs(fit.coefficients[2], 0.0, 1e-6),
                "quadratic coeff should be ~0, got {}",
                fit.coefficients[2]
            );
            assert!(
                close_abs(fit.coefficients[3], 0.0, 1e-6),
                "cubic coeff should be ~0, got {}",
                fit.coefficients[3]
            );
        }
    }

    /// 係数レベルの復元も確認（基底変換の丸めを許容した緩い tol）。
    /// t_center=0 の単純フレームで monomial 係数が q の係数に一致。
    #[test]
    fn fit_recovers_monomial_coefficients() {
        let q = vec![1.5, 2.0, -0.7, 0.4];
        let degree = 3;
        let (tc, hw) = (0.0, 1.0);
        let m = degree + 3;
        let taus = chebyshev_nodes(m);
        let qf = known_poly(q.clone());
        let values: Vec<f64> = taus.iter().map(|&tau| qf(tc + hw * tau)).collect();
        let fit = fit_chebyshev_monomial(&taus, &values, degree, tc, hw);
        assert_eq!(fit.coefficients.len(), q.len());
        for (i, &want) in q.iter().enumerate() {
            assert!(
                close_abs(fit.coefficients[i], want, 1e-6),
                "monomial coeff[{i}] got={} want={want}",
                fit.coefficients[i]
            );
        }
    }

    // ======================================================================
    // 契約5（任意）: 非多項式関数の近似は区間内で残差が小さい
    // ======================================================================

    /// sin(t) サンプルを degree=4 でフィットすると区間内残差が小さい（緩い tol）。
    #[test]
    fn fit_approximates_smooth_function() {
        let (tc, hw) = (0.0, 1.0); // t ∈ [-1, 1]
        let degree = 4;
        let m = degree + 4;
        let taus = chebyshev_nodes(m);
        let values: Vec<f64> = taus.iter().map(|&tau| (tc + hw * tau).sin()).collect();
        let fit = fit_chebyshev_monomial(&taus, &values, degree, tc, hw);
        let n = 100;
        let mut max_resid = 0.0_f64;
        for k in 0..=n {
            let frac = (k as f64) / (n as f64);
            let t = (tc - hw) + 2.0 * hw * frac;
            max_resid = max_resid.max((fit.eval(t) - t.sin()).abs());
        }
        // degree=4 Chebyshev 近似で sin on [-1,1] の残差は十分小さい（緩い tol）。
        assert!(
            max_resid < 1e-3,
            "sin approximation residual too large: {max_resid}"
        );
    }
}
