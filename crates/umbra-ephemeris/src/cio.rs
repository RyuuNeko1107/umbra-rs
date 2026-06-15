//! CIO ベースのフレーム変換（ISSUE-035 のうち CIO 経路, `docs/algorithms/02-frames.md` §2.1）。
//!
//! NPB 行列（`frames::bias_precession_nutation_matrix_iau2006`）から CIP 座標 X,Y を取り出し、
//! CIO locator s（SOFA `iauS06`）を評価して celestial-to-intermediate 行列 `Q`（GCRS→CIRS,
//! SOFA `iauC2ixys`/`iauC2i06a`）を構成する。さらに ERA・極運動を合成して GCRS→ITRS
//! （SOFA `iauC2t06a`）を与える。CIO ベースで統一（分点 GST は使わない, 確定D4）。
//!
//! s 級数（iauS06）は小規模（多項式6 + 周期66項）の公開系列で、`docs/algorithms/02-frames.md`
//! §181 が packed 生成対象とする 1365 項章動とは別格のためインライン（出典: IERS Conventions
//! 2010 Table 5.2c/d, SOFA `s06.c` は参照のみで非移植。基本引数は IERS Conventions 2003）。

use crate::frames::{bias_precession_nutation_matrix_iau2006, era_rotation, polar_motion_matrix};
use crate::nutation::fundamental_arguments;
use umbra_core::constants::ARCSEC_TO_RAD;
use umbra_core::{Matrix3, Radians, TtInstant, Ut1Instant};

/// CIO locator s 級数の 1 項: 8 基本引数の整数乗数 + sin/cos 振幅（秒角、`e-6` 込み）。
struct SeriesTerm {
    /// 乗数 `[l, l', F, D, Ω, L_Ve, L_E, p_A]`（s06 の 8 引数、`fundamental_arguments` の部分集合）。
    nfa: [f64; 8],
    /// sin(ARG) の係数 \[″\]。
    s: f64,
    /// cos(ARG) の係数 \[″\]。
    c: f64,
}

// CIO locator s の多項式・周期係数（SOFA `iauS06` = IERS Conventions 2010 Table 5.2d）。
// 単位は秒角（`e-6` を含む）。SOFA `s06.c` は参照のみで非移植、係数は IERS 標準データ。
// `s = (w0 + (w1 + (w2 + (w3 + (w4 + w5·t)·t)·t)·t)·t)·DAS2R − X·Y/2`、w_k は SP[k] に
// 各次数の周期級数を加えたもの。

/// 多項式部 SP[0..5]（秒角）。
const SP: [f64; 6] = [
    94.00e-6,
    3808.65e-6,
    -122.68e-6,
    -72574.11e-6,
    27.98e-6,
    15.62e-6,
];

// t⁰ の周期項（33）。
const S0: &[SeriesTerm] = &[
    SeriesTerm {
        nfa: [0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0],
        s: -2640.73e-6,
        c: 0.39e-6,
    },
    SeriesTerm {
        nfa: [0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0],
        s: -63.53e-6,
        c: 0.02e-6,
    },
    SeriesTerm {
        nfa: [0.0, 0.0, 2.0, -2.0, 3.0, 0.0, 0.0, 0.0],
        s: -11.75e-6,
        c: -0.01e-6,
    },
    SeriesTerm {
        nfa: [0.0, 0.0, 2.0, -2.0, 1.0, 0.0, 0.0, 0.0],
        s: -11.21e-6,
        c: -0.01e-6,
    },
    SeriesTerm {
        nfa: [0.0, 0.0, 2.0, -2.0, 2.0, 0.0, 0.0, 0.0],
        s: 4.57e-6,
        c: 0.00e-6,
    },
    SeriesTerm {
        nfa: [0.0, 0.0, 2.0, 0.0, 3.0, 0.0, 0.0, 0.0],
        s: -2.02e-6,
        c: 0.00e-6,
    },
    SeriesTerm {
        nfa: [0.0, 0.0, 2.0, 0.0, 1.0, 0.0, 0.0, 0.0],
        s: -1.98e-6,
        c: 0.00e-6,
    },
    SeriesTerm {
        nfa: [0.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0],
        s: 1.72e-6,
        c: 0.00e-6,
    },
    SeriesTerm {
        nfa: [0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0],
        s: 1.41e-6,
        c: 0.01e-6,
    },
    SeriesTerm {
        nfa: [0.0, 1.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0],
        s: 1.26e-6,
        c: 0.01e-6,
    },
    SeriesTerm {
        nfa: [1.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0],
        s: 0.63e-6,
        c: 0.00e-6,
    },
    SeriesTerm {
        nfa: [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0],
        s: 0.63e-6,
        c: 0.00e-6,
    },
    SeriesTerm {
        nfa: [0.0, 1.0, 2.0, -2.0, 3.0, 0.0, 0.0, 0.0],
        s: -0.46e-6,
        c: 0.00e-6,
    },
    SeriesTerm {
        nfa: [0.0, 1.0, 2.0, -2.0, 1.0, 0.0, 0.0, 0.0],
        s: -0.45e-6,
        c: 0.00e-6,
    },
    SeriesTerm {
        nfa: [0.0, 0.0, 4.0, -4.0, 4.0, 0.0, 0.0, 0.0],
        s: -0.36e-6,
        c: 0.00e-6,
    },
    SeriesTerm {
        nfa: [0.0, 0.0, 1.0, -1.0, 1.0, -8.0, 12.0, 0.0],
        s: 0.24e-6,
        c: 0.12e-6,
    },
    SeriesTerm {
        nfa: [0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        s: -0.32e-6,
        c: 0.00e-6,
    },
    SeriesTerm {
        nfa: [0.0, 0.0, 2.0, 0.0, 2.0, 0.0, 0.0, 0.0],
        s: -0.28e-6,
        c: 0.00e-6,
    },
    SeriesTerm {
        nfa: [1.0, 0.0, 2.0, 0.0, 3.0, 0.0, 0.0, 0.0],
        s: -0.27e-6,
        c: 0.00e-6,
    },
    SeriesTerm {
        nfa: [1.0, 0.0, 2.0, 0.0, 1.0, 0.0, 0.0, 0.0],
        s: -0.26e-6,
        c: 0.00e-6,
    },
    SeriesTerm {
        nfa: [0.0, 0.0, 2.0, -2.0, 0.0, 0.0, 0.0, 0.0],
        s: 0.21e-6,
        c: 0.00e-6,
    },
    SeriesTerm {
        nfa: [0.0, 1.0, -2.0, 2.0, -3.0, 0.0, 0.0, 0.0],
        s: -0.19e-6,
        c: 0.00e-6,
    },
    SeriesTerm {
        nfa: [0.0, 1.0, -2.0, 2.0, -1.0, 0.0, 0.0, 0.0],
        s: -0.18e-6,
        c: 0.00e-6,
    },
    SeriesTerm {
        nfa: [0.0, 0.0, 0.0, 0.0, 0.0, 8.0, -13.0, -1.0],
        s: 0.10e-6,
        c: -0.05e-6,
    },
    SeriesTerm {
        nfa: [0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0],
        s: -0.15e-6,
        c: 0.00e-6,
    },
    SeriesTerm {
        nfa: [2.0, 0.0, -2.0, 0.0, -1.0, 0.0, 0.0, 0.0],
        s: 0.14e-6,
        c: 0.00e-6,
    },
    SeriesTerm {
        nfa: [0.0, 1.0, 2.0, -2.0, 2.0, 0.0, 0.0, 0.0],
        s: 0.14e-6,
        c: 0.00e-6,
    },
    SeriesTerm {
        nfa: [1.0, 0.0, 0.0, -2.0, 1.0, 0.0, 0.0, 0.0],
        s: -0.14e-6,
        c: 0.00e-6,
    },
    SeriesTerm {
        nfa: [1.0, 0.0, 0.0, -2.0, -1.0, 0.0, 0.0, 0.0],
        s: -0.14e-6,
        c: 0.00e-6,
    },
    SeriesTerm {
        nfa: [0.0, 0.0, 4.0, -2.0, 4.0, 0.0, 0.0, 0.0],
        s: -0.13e-6,
        c: 0.00e-6,
    },
    SeriesTerm {
        nfa: [0.0, 0.0, 2.0, -2.0, 4.0, 0.0, 0.0, 0.0],
        s: 0.11e-6,
        c: 0.00e-6,
    },
    SeriesTerm {
        nfa: [1.0, 0.0, -2.0, 0.0, -3.0, 0.0, 0.0, 0.0],
        s: -0.11e-6,
        c: 0.00e-6,
    },
    SeriesTerm {
        nfa: [1.0, 0.0, -2.0, 0.0, -1.0, 0.0, 0.0, 0.0],
        s: -0.11e-6,
        c: 0.00e-6,
    },
];
// t¹ の周期項（3）。
const S1: &[SeriesTerm] = &[
    SeriesTerm {
        nfa: [0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0],
        s: -0.07e-6,
        c: 3.57e-6,
    },
    SeriesTerm {
        nfa: [0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0],
        s: 1.73e-6,
        c: -0.03e-6,
    },
    SeriesTerm {
        nfa: [0.0, 0.0, 2.0, -2.0, 3.0, 0.0, 0.0, 0.0],
        s: 0.00e-6,
        c: 0.48e-6,
    },
];
// t² の周期項（25）。
const S2: &[SeriesTerm] = &[
    SeriesTerm {
        nfa: [0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0],
        s: 743.52e-6,
        c: -0.17e-6,
    },
    SeriesTerm {
        nfa: [0.0, 0.0, 2.0, -2.0, 2.0, 0.0, 0.0, 0.0],
        s: 56.91e-6,
        c: 0.06e-6,
    },
    SeriesTerm {
        nfa: [0.0, 0.0, 2.0, 0.0, 2.0, 0.0, 0.0, 0.0],
        s: 9.84e-6,
        c: -0.01e-6,
    },
    SeriesTerm {
        nfa: [0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0],
        s: -8.85e-6,
        c: 0.01e-6,
    },
    SeriesTerm {
        nfa: [0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        s: -6.38e-6,
        c: -0.05e-6,
    },
    SeriesTerm {
        nfa: [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        s: -3.07e-6,
        c: 0.00e-6,
    },
    SeriesTerm {
        nfa: [0.0, 1.0, 2.0, -2.0, 2.0, 0.0, 0.0, 0.0],
        s: 2.23e-6,
        c: 0.00e-6,
    },
    SeriesTerm {
        nfa: [0.0, 0.0, 2.0, 0.0, 1.0, 0.0, 0.0, 0.0],
        s: 1.67e-6,
        c: 0.00e-6,
    },
    SeriesTerm {
        nfa: [1.0, 0.0, 2.0, 0.0, 2.0, 0.0, 0.0, 0.0],
        s: 1.30e-6,
        c: 0.00e-6,
    },
    SeriesTerm {
        nfa: [0.0, 1.0, -2.0, 2.0, -2.0, 0.0, 0.0, 0.0],
        s: 0.93e-6,
        c: 0.00e-6,
    },
    SeriesTerm {
        nfa: [1.0, 0.0, 0.0, -2.0, 0.0, 0.0, 0.0, 0.0],
        s: 0.68e-6,
        c: 0.00e-6,
    },
    SeriesTerm {
        nfa: [0.0, 0.0, 2.0, -2.0, 1.0, 0.0, 0.0, 0.0],
        s: -0.55e-6,
        c: 0.00e-6,
    },
    SeriesTerm {
        nfa: [1.0, 0.0, -2.0, 0.0, -2.0, 0.0, 0.0, 0.0],
        s: 0.53e-6,
        c: 0.00e-6,
    },
    SeriesTerm {
        nfa: [0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0],
        s: -0.27e-6,
        c: 0.00e-6,
    },
    SeriesTerm {
        nfa: [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0],
        s: -0.27e-6,
        c: 0.00e-6,
    },
    SeriesTerm {
        nfa: [1.0, 0.0, -2.0, -2.0, -2.0, 0.0, 0.0, 0.0],
        s: -0.26e-6,
        c: 0.00e-6,
    },
    SeriesTerm {
        nfa: [1.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0],
        s: -0.25e-6,
        c: 0.00e-6,
    },
    SeriesTerm {
        nfa: [1.0, 0.0, 2.0, 0.0, 1.0, 0.0, 0.0, 0.0],
        s: 0.22e-6,
        c: 0.00e-6,
    },
    SeriesTerm {
        nfa: [2.0, 0.0, 0.0, -2.0, 0.0, 0.0, 0.0, 0.0],
        s: -0.21e-6,
        c: 0.00e-6,
    },
    SeriesTerm {
        nfa: [2.0, 0.0, -2.0, 0.0, -1.0, 0.0, 0.0, 0.0],
        s: 0.20e-6,
        c: 0.00e-6,
    },
    SeriesTerm {
        nfa: [0.0, 0.0, 2.0, 2.0, 2.0, 0.0, 0.0, 0.0],
        s: 0.17e-6,
        c: 0.00e-6,
    },
    SeriesTerm {
        nfa: [2.0, 0.0, 2.0, 0.0, 2.0, 0.0, 0.0, 0.0],
        s: 0.13e-6,
        c: 0.00e-6,
    },
    SeriesTerm {
        nfa: [2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        s: -0.13e-6,
        c: 0.00e-6,
    },
    SeriesTerm {
        nfa: [1.0, 0.0, 2.0, -2.0, 2.0, 0.0, 0.0, 0.0],
        s: -0.12e-6,
        c: 0.00e-6,
    },
    SeriesTerm {
        nfa: [0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        s: -0.11e-6,
        c: 0.00e-6,
    },
];
// t³ の周期項（4）。
const S3: &[SeriesTerm] = &[
    SeriesTerm {
        nfa: [0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0],
        s: 0.30e-6,
        c: -23.42e-6,
    },
    SeriesTerm {
        nfa: [0.0, 0.0, 2.0, -2.0, 2.0, 0.0, 0.0, 0.0],
        s: -0.03e-6,
        c: -1.46e-6,
    },
    SeriesTerm {
        nfa: [0.0, 0.0, 2.0, 0.0, 2.0, 0.0, 0.0, 0.0],
        s: -0.01e-6,
        c: -0.25e-6,
    },
    SeriesTerm {
        nfa: [0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0],
        s: 0.00e-6,
        c: 0.23e-6,
    },
];
// t⁴ の周期項（1）。
const S4: &[SeriesTerm] = &[SeriesTerm {
    nfa: [0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0],
    s: -0.26e-6,
    c: -0.01e-6,
}];

/// CIO locator s（SOFA `iauS06`）。`time_tt` と CIP 座標 X,Y から評価する \[rad\]。
fn cio_locator_s06(time_tt: TtInstant, x: f64, y: f64) -> f64 {
    let t = time_tt.jd2().julian_centuries_since_j2000();
    // s06 の 8 引数 = fundamental_arguments の [l,l',F,D,Ω,L_Ve,L_E,p_A]（idx 0..4,6,7,13）。
    let fa14 = fundamental_arguments(t);
    let fa = [
        fa14[0], fa14[1], fa14[2], fa14[3], fa14[4], fa14[6], fa14[7], fa14[13],
    ];
    let mut w = SP;
    // w[0..4] に t⁰..t⁴ の周期級数を加算（小振幅項から＝表末尾から、SOFA に倣う丸め抑制）。
    for (slot, terms) in w.iter_mut().zip([S0, S1, S2, S3, S4]) {
        for term in terms.iter().rev() {
            let arg: f64 = term.nfa.iter().zip(fa.iter()).map(|(n, f)| n * f).sum();
            *slot += term.s * arg.sin() + term.c * arg.cos();
        }
    }
    let poly = w[0] + t * (w[1] + t * (w[2] + t * (w[3] + t * (w[4] + t * w[5]))));
    poly * ARCSEC_TO_RAD - x * y / 2.0
}

/// CIP 座標 X,Y と CIO locator s から celestial-to-intermediate 行列を構成する。
/// SOFA `iauC2ixys`: `R3(−(E+s))·R2(d)·R3(E)`、E=atan2(Y,X)、d=atan√(r²/(1−r²))。
/// `docs/algorithms/02-frames.md` §2.1 (F1) の `a=1/(1+Z)` 行列形と等価。
fn celestial_to_intermediate(x: f64, y: f64, s: f64) -> Matrix3 {
    let r2 = x * x + y * y;
    // atan2(0,0)=0（IEEE/Rust 規約）。X=Y=0 でも e=0 になるため SOFA の r2>0 ガードは不要。
    let e = y.atan2(x);
    let d = (r2 / (1.0 - r2)).sqrt().atan();
    Matrix3::rotation_z(-(e + s))
        .mul_mat(&Matrix3::rotation_y(d))
        .mul_mat(&Matrix3::rotation_z(e))
}

/// GCRS→CIRS の celestial-to-intermediate 行列 `Q`（SOFA `iauC2i06a` 相当）。TT 引数。
/// NPB 行列（`iauPnm06a`）第3行から X=NPB\[2\]\[0\], Y=NPB\[2\]\[1\]、s は `iauS06`。
pub fn gcrs_to_cirs_matrix(time_tt: TtInstant) -> Matrix3 {
    let npb = bias_precession_nutation_matrix_iau2006(time_tt);
    let x = npb.rows[2][0];
    let y = npb.rows[2][1];
    let s = cio_locator_s06(time_tt, x, y);
    celestial_to_intermediate(x, y, s)
}

/// GCRS→ITRS の合成行列（SOFA `iauC2t06a` 相当）= 極運動 · R3(ERA) · C2I。
/// `time_tt` で歳差章動・s・s′、`time_ut1` で ERA、`xp,yp` で極運動。
pub fn gcrs_to_itrs_matrix(
    time_tt: TtInstant,
    time_ut1: Ut1Instant,
    xp: Radians,
    yp: Radians,
) -> Matrix3 {
    let c2i = gcrs_to_cirs_matrix(time_tt);
    let era = era_rotation(time_ut1);
    let pom = polar_motion_matrix(xp, yp, time_tt);
    pom.mul_mat(&era).mul_mat(&c2i)
}

#[cfg(test)]
mod tests {
    // SOFA/ERFA の検証値は f64 の表現可能桁数より多くの桁で配布されている。
    // provenance（一次ソースの逐語転記）を保つため桁を削らず、ここに限り
    // 過剰精度リント（余剰桁は最近接 f64 へ丸められ値は不変）を許可する。
    #![allow(clippy::excessive_precision)]

    use super::*;
    use umbra_core::JulianDate2;

    // ------------------------------------------------------------------
    // 一次オラクル: pyerfa（liberfa C ラッパ ＝ SOFA と数値同一）。
    //   取得: docker run --rm python:3.12-slim bash -c \
    //     "pip install -q pyerfa && python3 -c 'import erfa; ...'"
    //   erfa.__version__ = 2.0.1.5
    //   pyerfa は liberfa（SOFA 由来の C 実装）の薄い numpy ラッパで、本実装の
    //   純 Rust とは独立実装。戻り値は numpy 配列を float(...) 化して逐語転記。
    // 入力はすべて SOFA 流の 2要素 JD (date1=2400000.5, date2=MJD) ほか。
    // 各定数の provenance は使用箇所のコメントに erfa 関数名・入力・全要素を明記。
    // ------------------------------------------------------------------

    // SOFA の 2要素 JD は (date1, date2)。JulianDate2 へそのまま渡す。
    fn jd2(date1: f64, date2: f64) -> JulianDate2 {
        JulianDate2::new(date1, date2)
    }
    fn tt(date1: f64, date2: f64) -> TtInstant {
        TtInstant::from_jd2(jd2(date1, date2))
    }
    fn ut1(date1: f64, date2: f64) -> Ut1Instant {
        Ut1Instant::from_jd2(jd2(date1, date2))
    }

    /// 許容つきスカラ比較。
    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    /// 行列を期待値（行優先 \[行\]\[列\]）と要素ごとに比較。
    fn matrix_close(m: &Matrix3, expected: &[[f64; 3]; 3], tol: f64) -> bool {
        m.rows
            .iter()
            .zip(expected.iter())
            .all(|(mr, er)| mr.iter().zip(er.iter()).all(|(a, b)| (a - b).abs() < tol))
    }

    /// 直交性: M·Mᵀ ≈ I。
    fn assert_orthonormal(m: &Matrix3, tol: f64) {
        let prod = m.mul_mat(&m.transpose());
        for (i, row) in prod.rows.iter().enumerate() {
            for (j, &val) in row.iter().enumerate() {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(
                    (val - expected).abs() < tol,
                    "M·Mᵀ[{i}][{j}] = {val} (expected {expected})"
                );
            }
        }
    }

    /// 行列式（右手系の確認用、≈ +1 を期待）。
    fn det(m: &Matrix3) -> f64 {
        let r = &m.rows;
        r[0][0] * (r[1][1] * r[2][2] - r[1][2] * r[2][1])
            - r[0][1] * (r[1][0] * r[2][2] - r[1][2] * r[2][0])
            + r[0][2] * (r[1][0] * r[2][1] - r[1][1] * r[2][0])
    }

    // 許容（TOL）設計メモ（精度クリティカル）:
    //   gcrs_to_cirs/gcrs_to_itrs は CIP 座標 X,Y の章動由来項を含む。本実装の章動は
    //   IERS IAU 2000_R06 系列の直接評価、erfa の C2I（c2i06a）は nut00a の線形
    //   スケーリング近似に由来し、両者の章動表現差は X,Y で ~2e-11 rad。
    //   s06・ERA・極運動は同一式のため ~1e-13 で一致する。よって行列要素差は
    //   X,Y の representation 差（~2e-11）で律速される。一方、CIO 合成・s 評価・
    //   X,Y 抽出・ERA/極運動の合成順の実バグは行列要素を ≫1e-6 ずらす。
    //   要素ごと許容 TOL=1e-9 は representation 差（~2e-11）を吸収しつつ実バグを
    //   確実に検出する（NPB 突合の 1e-10 より一段緩めて C2I 合成段の丸めも吸収）。
    const TOL: f64 = 1e-9;

    // gcrs_to_itrs テストの極運動入力（erfa t_pom00 / c2t06a 検証と同一）。
    const XP: f64 = 2.55060238e-7;
    const YP: f64 = 1.860359247e-6;

    // ==================================================================
    // 1. gcrs_to_cirs vs erfa.c2i06a（GCRS→CIRS, celestial-to-intermediate）
    // ==================================================================

    /// SOFA iauC2i06a（pyerfa `erfa.c2i06a(d1, d2)`）多エポック突合。
    /// J2000 / SOFA-test / 1900 / 2100 の 4 エポックで C2I 全 9 要素が TOL 以内一致。
    /// 1900(t≈-1) / 2100(t≈+1) で歳差・章動・s の t¹ 項を両符号に励起する。
    ///
    /// オラクル: pyerfa 2.0.1.5 `erfa.c2i06a(date1, date2)`（liberfa C ラッパ、行優先 \[行\]\[列\]）。
    #[test]
    fn gcrs_to_cirs_matches_c2i06a_multi_epoch() {
        // erfa 2.0.1.5 erfa.c2i06a(2451545.0, 0.0) [TT]:
        const J2000: [[f64; 3]; 3] = [
            [
                0.9999999996369462,
                9.756652263881449e-09,
                2.694638043284609e-05,
            ],
            [
                -1.0511278236702282e-08,
                0.9999999996078677,
                2.8004720891691253e-05,
            ],
            [
                -2.694638014904722e-05,
                -2.800472116476493e-05,
                0.9999999992448141,
            ],
        ];
        // erfa 2.0.1.5 erfa.c2i06a(2400000.5, 53736.0) [TT]:
        const SOFA: [[f64; 3]; 3] = [
            [
                0.999999832303716,
                5.581121398368083e-10,
                -0.0005791308487740529,
            ],
            [
                -2.384253169895878e-08,
                0.9999999991917468,
                -4.020579392909148e-05,
            ],
            [
                0.0005791308482835292,
                4.020580099467486e-05,
                0.9999998314954629,
            ],
        ];
        // erfa 2.0.1.5 erfa.c2i06a(2400000.5, 15020.0) [TT] (≈ 1900, t≈-1 JC):
        const Y1900: [[f64; 3]; 3] = [
            [
                0.9999531110123181,
                -3.42094249674435e-07,
                0.00968378937552917,
            ],
            [
                -8.092428920476735e-07,
                0.9999999929323365,
                0.00011888932628313787,
            ],
            [
                -0.009683789347758761,
                -0.00011889158822070249,
                0.9999531039447093,
            ],
        ];
        // erfa 2.0.1.5 erfa.c2i06a(2400000.5, 88069.0) [TT] (≈ 2100, t≈+1 JC):
        const Y2100: [[f64; 3]; 3] = [
            [
                0.9999527538307167,
                3.3193606402241727e-07,
                -0.009720602155013537,
            ],
            [
                3.233043075834621e-07,
                0.9999999977281759,
                6.740581349921874e-05,
            ],
            [
                0.009720602155304459,
                -6.740577154529248e-05,
                0.9999527515588925,
            ],
        ];

        for (label, (d1, d2), expected) in [
            ("J2000", (2451545.0, 0.0), &J2000),
            ("SOFA", (2400000.5, 53736.0), &SOFA),
            ("1900", (2400000.5, 15020.0), &Y1900),
            ("2100", (2400000.5, 88069.0), &Y2100),
        ] {
            let m = gcrs_to_cirs_matrix(tt(d1, d2));
            assert!(
                matrix_close(&m, expected, TOL),
                "c2i06a mismatch at {label} ({d1}, {d2}); got {:?}, expected {expected:?}",
                m.rows
            );
        }
    }

    /// C2I は回転行列 → 直交かつ det≈+1（全 4 エポックで確認）。
    #[test]
    fn gcrs_to_cirs_is_orthonormal() {
        for &(d1, d2) in &[
            (2451545.0, 0.0),
            (2400000.5, 53736.0),
            (2400000.5, 15020.0),
            (2400000.5, 88069.0),
        ] {
            let m = gcrs_to_cirs_matrix(tt(d1, d2));
            assert_orthonormal(&m, 1e-12);
            assert!(
                close(det(&m), 1.0, 1e-12),
                "det = {} at ({d1}, {d2})",
                det(&m)
            );
        }
    }

    // ==================================================================
    // 2. gcrs_to_itrs vs erfa.c2t06a（GCRS→ITRS = 極運動·R3(ERA)·C2I）
    // ==================================================================

    /// SOFA iauC2t06a（pyerfa `erfa.c2t06a(tta,ttb, uta,utb, xp,yp)`）突合。
    /// 本テストでは UT1=TT（uta=tta, utb=ttb）、xp=2.55060238e-7, yp=1.860359247e-6 を使用。
    /// SOFA-test / 2100 の 2 ケースで全 9 要素が TOL 以内一致。
    /// 2100(t≈+1) で歳差・章動・s および ERA の永年項を励起する。
    ///
    /// オラクル: pyerfa 2.0.1.5 `erfa.c2t06a(d1,d2, d1,d2, xp,yp)`（行優先 \[行\]\[列\]）。
    /// 取得: docker run --rm python:3.12-slim 内で
    ///   xp,yp = 2.55060238e-7, 1.860359247e-6
    ///   erfa.c2t06a(d1,d2,d1,d2,xp,yp) for (d1,d2) in [(2400000.5,53736.0),(2400000.5,88069.0)]
    #[test]
    fn gcrs_to_itrs_matches_c2t06a() {
        // erfa 2.0.1.5 erfa.c2t06a(2400000.5,53736.0, 2400000.5,53736.0, XP,YP):
        const SOFA: [[f64; 3]; 3] = [
            [
                -0.18103321283058613,
                0.9834769806938598,
                6.555550962984977e-05,
            ],
            [
                -0.9834768134136221,
                -0.18103322036490946,
                0.0005749800844905839,
            ],
            [
                0.0005773474024748545,
                3.961816829646157e-05,
                0.9999998325501748,
            ],
        ];
        // erfa 2.0.1.5 erfa.c2t06a(2400000.5,88069.0, 2400000.5,88069.0, XP,YP) (≈ 2100):
        const Y2100: [[f64; 3]; 3] = [
            [-0.16429298597354, 0.9864101815991755, 0.0016638501674712722],
            [
                -0.9863637054911221,
                -0.16430139688544593,
                0.009575566370424785,
            ],
            [
                0.009718809069089283,
                -6.796302518408726e-05,
                0.9999527689502669,
            ],
        ];

        for (label, (d1, d2), expected) in [
            ("SOFA", (2400000.5, 53736.0), &SOFA),
            ("2100", (2400000.5, 88069.0), &Y2100),
        ] {
            // UT1=TT: time_ut1 に time_tt と同じ 2要素 JD を渡す。
            let m =
                gcrs_to_itrs_matrix(tt(d1, d2), ut1(d1, d2), Radians::new(XP), Radians::new(YP));
            assert!(
                matrix_close(&m, expected, TOL),
                "c2t06a mismatch at {label} ({d1}, {d2}); got {:?}, expected {expected:?}",
                m.rows
            );
        }
    }

    /// C2T も回転行列 → 直交かつ det≈+1（SOFA / 2100 の 2 ケース）。
    #[test]
    fn gcrs_to_itrs_is_orthonormal() {
        for &(d1, d2) in &[(2400000.5, 53736.0), (2400000.5, 88069.0)] {
            let m =
                gcrs_to_itrs_matrix(tt(d1, d2), ut1(d1, d2), Radians::new(XP), Radians::new(YP));
            assert_orthonormal(&m, 1e-12);
            assert!(
                close(det(&m), 1.0, 1e-12),
                "det = {} at ({d1}, {d2})",
                det(&m)
            );
        }
    }

    // ==================================================================
    // 3. cio_locator_s06 を erfa.s06（= xys06a の s）と直接突合。
    //    gcrs_to_cirs の 1e-9 行列許容では s（~µas）の係数符号 flip が埋もれるため、
    //    s を直接 µas 級で突合して s06 係数・乗数・多項式・−XY/2 を検証する。
    //    erfa の X,Y を渡して s06 系列を NPB 由来の X,Y 差から分離する。
    //    複数エポック（特に t=10 で t²/t³/t⁴ 群を増幅励起）で全係数を踏む。
    //    オラクル: pyerfa 2.0.1.5 `erfa.xys06a(d1,d2) -> (X, Y, s)`。
    // ==================================================================
    #[test]
    fn cio_locator_s06_matches_erfa() {
        // (d1, d2, X, Y, s)（erfa.xys06a, 全要素逐語転記）。
        const CASES: [(f64, f64, f64, f64, f64); 5] = [
            (
                2451545.0,
                0.0,
                -2.694638014904722e-05,
                -2.8004721164764934e-05,
                -1.0133965177563803e-08,
            ),
            (
                2400000.5,
                53736.0,
                0.0005791308482835291,
                4.0205800994674856e-05,
                -1.22003229416848e-08,
            ),
            (
                2400000.5,
                15020.0,
                -0.009683789347758761,
                -0.00011889158822070423,
                -2.3357979805368332e-07,
            ),
            (
                2400000.5,
                88069.0,
                0.009720602155304459,
                -6.740577154529248e-05,
                -4.315980180657891e-09,
            ),
            (
                2816795.0,
                0.0,
                0.09601386055497964,
                -0.010840769266888484,
                0.0001774778429708768,
            ),
        ];
        // s は同一係数・同一加算順のため erfa と ~1e-15 rad で一致する。許容 5e-14 rad
        // （= ~0.01 µas）は最小 cos 係数（0.01 µas）の符号 flip（~9.7e-14 rad）まで確実に検出する。
        const TOL_S: f64 = 5e-14;
        for (d1, d2, x, y, s_erfa) in CASES {
            let s = cio_locator_s06(tt(d1, d2), x, y);
            assert!(
                close(s, s_erfa, TOL_S),
                "s06 mismatch at ({d1}, {d2}): {s} vs erfa {s_erfa} (|diff| = {:e})",
                (s - s_erfa).abs()
            );
        }
    }
}
