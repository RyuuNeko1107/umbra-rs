//! 数値解法（`docs/numerical-policy.md` §A5）。
//!
//! - [`brent_root`] … 導関数不要・ブラケット必須の堅牢な求根（Newton 単独は禁止, conventions §11）。
//!   接触・合・地平線通過の標準解法。
//! - [`minimize_golden`] … 黄金分割による一次元最小化（粗ブラケット用）。
//!
//! 最大食時刻は「距離の最小化」ではなく `dm/dt = 0` の[`brent_root`]で解くのが正式手法
//! （numerical-policy §A5 / D2）。最小化はその粗いブラケット取りに使う。

use crate::error::SolverError;

/// Brent 法による求根。`f(a)` と `f(b)` が異符号（根をブラケット）である必要がある。
///
/// `tol` は区間幅 `|b − a|` の収束閾値（秒など、引数 `x` の単位）。`max_iter` は反復上限。
/// 区間幅が `tol` を下回るか `f` がちょうど 0 になった時点の推定値を返す。
pub fn brent_root<F>(
    mut f: F,
    mut a: f64,
    mut b: f64,
    tol: f64,
    max_iter: usize,
) -> Result<f64, SolverError>
where
    F: FnMut(f64) -> f64,
{
    let mut fa = f(a);
    let mut fb = f(b);
    if fa == 0.0 {
        return Ok(a);
    }
    if fb == 0.0 {
        return Ok(b);
    }
    if fa * fb > 0.0 {
        return Err(SolverError::RootNotBracketed);
    }
    // b を最良推定（|f| 最小）に保つ。
    if fa.abs() < fb.abs() {
        core::mem::swap(&mut a, &mut b);
        core::mem::swap(&mut fa, &mut fb);
    }
    let mut c = a;
    let mut fc = fa;
    let mut d = a; // !mflag のときのみ参照される。初期 mflag=true。
    let mut mflag = true;

    for _ in 0..max_iter {
        if fb == 0.0 || (b - a).abs() < tol {
            return Ok(b);
        }
        if !fa.is_finite() || !fb.is_finite() {
            return Err(SolverError::NumericalInstability);
        }

        // 逆二次補間、不可なら割線。
        let mut s = if fa != fc && fb != fc {
            a * fb * fc / ((fa - fb) * (fa - fc))
                + b * fa * fc / ((fb - fa) * (fb - fc))
                + c * fa * fb / ((fc - fa) * (fc - fb))
        } else {
            b - fb * (b - a) / (fb - fa)
        };

        // 補間が不適なら二分法へ切替。
        let lo = (3.0 * a + b) / 4.0;
        let s_between_lo_and_b = (s - lo) * (s - b) < 0.0;
        let reject = !s_between_lo_and_b
            || (mflag && (s - b).abs() >= (b - c).abs() / 2.0)
            || (!mflag && (s - b).abs() >= (c - d).abs() / 2.0)
            || (mflag && (b - c).abs() < tol)
            || (!mflag && (c - d).abs() < tol);
        if reject {
            s = (a + b) / 2.0;
            mflag = true;
        } else {
            mflag = false;
        }

        let fs = f(s);
        d = c;
        c = b;
        fc = fb;
        if fa * fs < 0.0 {
            b = s;
            fb = fs;
        } else {
            a = s;
            fa = fs;
        }
        if fa.abs() < fb.abs() {
            core::mem::swap(&mut a, &mut b);
            core::mem::swap(&mut fa, &mut fb);
        }
    }
    Err(SolverError::DidNotConverge)
}

/// 黄金分割による一次元最小化（区間 `[a, b]` で単峰を仮定）。最小点の推定値を返す。
///
/// 線形収束のため精密化には不向き（numerical-policy §A5: 最大食は `dm/dt = 0` 求根を正式手法とし、
/// 本関数は粗いブラケット取りに使う）。
pub fn minimize_golden<F>(mut f: F, mut a: f64, mut b: f64, tol: f64, max_iter: usize) -> f64
where
    F: FnMut(f64) -> f64,
{
    let gr = (5.0_f64.sqrt() - 1.0) / 2.0; // ≈0.618
    let mut c = b - gr * (b - a);
    let mut d = a + gr * (b - a);
    let mut fc = f(c);
    let mut fd = f(d);
    for _ in 0..max_iter {
        if (b - a).abs() < tol {
            break;
        }
        if fc < fd {
            b = d;
            d = c;
            fd = fc;
            c = b - gr * (b - a);
            fc = f(c);
        } else {
            a = c;
            c = d;
            fc = fd;
            d = a + gr * (b - a);
            fd = f(d);
        }
    }
    (a + b) / 2.0
}

#[cfg(test)]
mod tests {
    use super::*;
    use core::f64::consts::PI;

    #[test]
    fn brent_finds_cos_root() {
        // cos(x) = 0 の [0, 3] における根は π/2。
        let r = brent_root(|x| x.cos(), 0.0, 3.0, 1e-12, 100).unwrap();
        assert!((r - PI / 2.0).abs() < 1e-10, "r = {r}");
    }

    #[test]
    fn brent_finds_sqrt2() {
        let r = brent_root(|x| x * x - 2.0, 0.0, 2.0, 1e-12, 100).unwrap();
        assert!((r - 2.0_f64.sqrt()).abs() < 1e-10, "r = {r}");
    }

    #[test]
    fn brent_errors_when_not_bracketed() {
        let e = brent_root(|x| x * x + 1.0, 0.0, 2.0, 1e-12, 100).unwrap_err();
        assert_eq!(e, SolverError::RootNotBracketed);
    }

    #[test]
    fn brent_hits_exact_endpoint_root() {
        let r = brent_root(|x| x, -1.0, 0.0, 1e-12, 100).unwrap();
        assert_eq!(r, 0.0);
    }

    #[test]
    fn golden_minimizes_parabola() {
        let m = minimize_golden(|x| (x - 2.0) * (x - 2.0), 0.0, 5.0, 1e-10, 200);
        assert!((m - 2.0).abs() < 1e-6, "m = {m}");
    }

    #[test]
    fn golden_minimizes_cos() {
        // cos の [0, 2π] における最小は π。
        let m = minimize_golden(|x| x.cos(), 0.0, 2.0 * PI, 1e-10, 200);
        assert!((m - PI).abs() < 1e-5, "m = {m}");
    }
}
