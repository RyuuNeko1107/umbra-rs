//! 月 ELP2000-82B 評価（ISSUE-014）。
//!
//! 正本: `docs/algorithms/elp2000-82b-evaluation.md`（著者 Fortran `elp82b_1` から転記、
//! Chapront-Touzé & Chapront, IMCCE MCJCGF.9601）。
//!
//! S1（本コミット）= 基本引数 [`moon_arguments`]。S2 = 36 系列の総和＋組立＋J2000 回転。
//!
//! S1 の公開 IF はまだ評価器(S2)が consume しないため `dead_code` を一時許容（S2 で解除）。
#![allow(dead_code)]

/// 時刻 `t`（J2000 からのユリウス世紀 TDB, 無次元）における ELP2000-82B 基本引数（rad）。
/// 引数は**正規化しない**生の多項式和（t と共に増大する）。
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct MoonArguments {
    /// Delaunay 引数 D, l', l, F を完全多項式（t⁰..t⁴）で評価（主問題用）。
    pub delaunay_full: [f64; 4],
    /// Delaunay 引数 D, l', l, F を線形（t⁰, t¹）で評価（摂動用）。
    pub delaunay_lin: [f64; 4],
    /// zeta = W1 + 歳差, 線形（形状/潮汐/相対論/太陽離心率摂動用）。
    pub zeta_lin: f64,
    /// 8 惑星平均黄経 Me,Ve,Earth,Ma,Ju,Sa,Ur,Ne を線形評価（惑星摂動用）。
    pub planets_lin: [f64; 8],
    /// 月平均黄経 W1(t) 完全多項式（経度に加算）。
    pub w1: f64,
}

/// 秒角→rad 換算の逆数 `rad = 648000/π`（`elp82b_1` の `rad`、秒角/rad）。
const RAD: f64 = 648000.0 / core::f64::consts::PI;
/// 度→rad。
const DEG: f64 = core::f64::consts::PI / 180.0;
/// 歳差定数（秒角/世紀 → rad/世紀）。
const PRECES: f64 = 5029.0966 / RAD;

/// 度分秒（度・分・秒）を rad へ（定数項用）。
fn dms_to_rad(d: f64, m: f64, s: f64) -> f64 {
    (d + m / 60.0 + s / 3600.0) * DEG
}

/// 月平均黄経/近地点/昇交点 W1,W2,W3 の多項式係数（t⁰..t⁴, rad）。
fn w_coeffs() -> [[f64; 5]; 3] {
    [
        [
            dms_to_rad(218.0, 18.0, 59.95571),
            1732559343.73604 / RAD,
            -5.8883 / RAD,
            0.6604e-2 / RAD,
            -0.3169e-4 / RAD,
        ],
        [
            dms_to_rad(83.0, 21.0, 11.67475),
            14643420.2632 / RAD,
            -38.2776 / RAD,
            -0.45047e-1 / RAD,
            0.21301e-3 / RAD,
        ],
        [
            dms_to_rad(125.0, 2.0, 40.39816),
            -6967919.3622 / RAD,
            6.3622 / RAD,
            0.7625e-2 / RAD,
            -0.3586e-4 / RAD,
        ],
    ]
}

/// 地球(EMB)平均黄経 T の係数（t⁰..t⁴, rad）。
fn eart_coeffs() -> [f64; 5] {
    [
        dms_to_rad(100.0, 27.0, 59.22059),
        129597742.2758 / RAD,
        -0.0202 / RAD,
        0.9e-5 / RAD,
        0.15e-6 / RAD,
    ]
}

/// EMB 近日点 ϖ' の係数（t⁰..t⁴, rad）。
fn peri_coeffs() -> [f64; 5] {
    [
        dms_to_rad(102.0, 56.0, 14.42753),
        1161.2283 / RAD,
        0.5327 / RAD,
        -0.138e-3 / RAD,
        0.0,
    ]
}

/// 8 惑星平均黄経 Me,Ve,Earth,Ma,Ju,Sa,Ur,Ne の線形係数（t⁰, t¹, rad）。
fn planet_coeffs() -> [[f64; 2]; 8] {
    let eart = eart_coeffs();
    [
        [dms_to_rad(252.0, 15.0, 3.25986), 538101628.68898 / RAD],
        [dms_to_rad(181.0, 58.0, 47.28305), 210664136.43355 / RAD],
        [eart[0], eart[1]],
        [dms_to_rad(355.0, 25.0, 59.78866), 68905077.59284 / RAD],
        [dms_to_rad(34.0, 21.0, 5.34212), 10925660.42861 / RAD],
        [dms_to_rad(50.0, 4.0, 38.89694), 4399609.65932 / RAD],
        [dms_to_rad(314.0, 3.0, 18.01841), 1542481.19393 / RAD],
        [dms_to_rad(304.0, 20.0, 55.19575), 786550.32074 / RAD],
    ]
}

/// Delaunay 引数 D, l', l, F の多項式係数（各 t⁰..t⁴, rad）。
/// D=W1−T(+π), l'=T−ϖ', l=W1−W2, F=W1−W3（`elp82b_1` の del 構成）。
fn delaunay_coeffs() -> [[f64; 5]; 4] {
    let w = w_coeffs();
    let eart = eart_coeffs();
    let peri = peri_coeffs();
    let mut del = [[0.0_f64; 5]; 4];
    for k in 0..5 {
        del[0][k] = w[0][k] - eart[k]; // D
        del[1][k] = eart[k] - peri[k]; // l'
        del[2][k] = w[0][k] - w[1][k]; // l
        del[3][k] = w[0][k] - w[2][k]; // F
    }
    del[0][0] += core::f64::consts::PI; // D の定数項に +180°
    del
}

/// 多項式（昇べき係数）を t で評価。
fn poly(coeffs: &[f64], t: f64) -> f64 {
    coeffs.iter().rev().fold(0.0, |acc, &c| acc * t + c)
}

/// 時刻 `t`（世紀 TDB）における ELP2000-82B 基本引数を評価する。
pub(crate) fn moon_arguments(t: f64) -> MoonArguments {
    let del = delaunay_coeffs();
    let w = w_coeffs();
    let planets = planet_coeffs();

    let delaunay_full = [
        poly(&del[0], t),
        poly(&del[1], t),
        poly(&del[2], t),
        poly(&del[3], t),
    ];
    let delaunay_lin = [
        del[0][0] + del[0][1] * t,
        del[1][0] + del[1][1] * t,
        del[2][0] + del[2][1] * t,
        del[3][0] + del[3][1] * t,
    ];
    // zeta = W1 + 歳差（線形）。
    let zeta_lin = w[0][0] + (w[0][1] + PRECES) * t;
    let mut planets_lin = [0.0_f64; 8];
    for (out, p) in planets_lin.iter_mut().zip(planets.iter()) {
        *out = p[0] + p[1] * t;
    }
    let w1 = poly(&w[0], t);

    MoonArguments {
        delaunay_full,
        delaunay_lin,
        zeta_lin,
        planets_lin,
        w1,
    }
}

#[cfg(test)]
#[allow(clippy::excessive_precision)]
mod tests {
    use super::*;

    /// 絶対許容での近接表明（`clippy::float_cmp` 回避）。
    #[track_caller]
    fn assert_close(actual: f64, expected: f64, tol: f64, msg: &str) {
        let d = (actual - expected).abs();
        assert!(d <= tol, "{msg}: |{actual} − {expected}| = {d:e} > {tol:e}");
    }

    /// オラクル = 著者 Fortran elp82b_1 の独立 Python 移植（公開値 JD 2451555.5 に 5e-10 km
    /// で一致＝信頼可）。引数は生の多項式和（非正規化）。
    /// 許容 1e-7 rad: 最大 ~8400 rad で f64 丸めは ~1e-11、転記誤りは ≫1e-4 → 安全余裕大。
    const TOL: f64 = 1e-7;

    fn check(t: f64, df: [f64; 4], dl: [f64; 4], zl: f64, pl: [f64; 8], w1: f64) {
        let a = moon_arguments(t);
        for (i, (got, exp)) in a.delaunay_full.iter().zip(df.iter()).enumerate() {
            assert_close(*got, *exp, TOL, &format!("delaunay_full[{i}]"));
        }
        for (i, (got, exp)) in a.delaunay_lin.iter().zip(dl.iter()).enumerate() {
            assert_close(*got, *exp, TOL, &format!("delaunay_lin[{i}]"));
        }
        assert_close(a.zeta_lin, zl, TOL, "zeta_lin");
        for (i, (got, exp)) in a.planets_lin.iter().zip(pl.iter()).enumerate() {
            assert_close(*got, *exp, TOL, &format!("planets_lin[{i}]"));
        }
        assert_close(a.w1, w1, TOL, "w1");
    }

    /// J2000（t=0）: 全引数 = 定数項。delaunay_full と delaunay_lin は一致。
    #[test]
    fn arguments_at_j2000() {
        check(
            0.0,
            [
                5.198466741027443,
                -0.04312518020812495,
                2.3555558982657994,
                1.627905233371468,
            ],
            [
                5.198466741027443,
                -0.04312518020812495,
                2.3555558982657994,
                1.627905233371468,
            ],
            3.810344430588308,
            [
                4.4026088424029615,
                3.1761466969075944,
                1.753470343150658,
                6.203480913399945,
                0.5995464973886735,
                0.8740167565184808,
                5.481293871604991,
                5.311886286783467,
            ],
            3.810344430588308,
        );
    }

    /// 正の小 t（JD 2469000.5）。
    #[test]
    fn arguments_positive_small_t() {
        check(
            0.4779055441478439,
            [
                3719.1826843720487,
                300.22586198160286,
                3982.6834002283654,
                4032.0281250218522,
            ],
            [
                3719.1826908662197,
                300.22586259374236,
                3982.6833643369428,
                4032.0281385871735,
            ],
            4018.077899012051,
            [
                1251.157963497425,
                491.2747253468396,
                302.02514862670716,
                165.85320109933772,
                25.913771387767575,
                11.067698783672512,
                9.055147927152627,
                7.134285201422363,
            ],
            4018.0662403228034,
        );
    }

    /// 負の t（JD 2415020.5, t≈−1）: 符号と高次項を励起。
    #[test]
    fn arguments_negative_t() {
        check(
            -0.9999863107460644,
            [
                -7766.072324196238,
                -628.336482044857,
                -8326.221700713522,
                -8431.722864422874,
            ],
            [
                -7766.072295715539,
                -628.3364793636836,
                -8326.221857485394,
                -8431.72280503737,
            ],
            -8395.783783340637,
            [
                -2604.3519929219324,
                -1018.1384266982653,
                -626.5455135569244,
                -327.85318918664217,
                -52.36882490467052,
                -20.455600796733357,
                -1.9967636146801127,
                1.4986349243058288,
            ],
            -8395.759430504724,
        );
    }

    /// t=1.0（JD 2488070.0）: 全 t^k を等重みで励起。
    #[test]
    fn arguments_at_t_one() {
        check(
            1.0,
            [
                7776.575585135253,
                628.2588273084584,
                8331.047140130802,
                8435.09400396688,
            ],
            [
                7776.575613552785,
                628.2588299882799,
                8331.04698285382,
                8435.094063363911,
            ],
            8403.519457952856,
            [
                2613.1929229998136,
                1024.5047013180165,
                630.061055305306,
                340.2647240626296,
                53.568643006860725,
                22.203926300318486,
                12.959453728319344,
                9.125189850541922,
            ],
            8403.49504768908,
        );
    }

    /// t=0 では delaunay_full == delaunay_lin（定数項のみ）。full/lin 分離の起点を固定。
    #[test]
    fn full_equals_lin_at_j2000() {
        let a = moon_arguments(0.0);
        assert_eq!(
            a.delaunay_full, a.delaunay_lin,
            "at t=0 full and linear Delaunay must coincide"
        );
    }

    /// 遠方エポック t=10（JD 2816795.0）: t³ を ×1000、t⁴ を ×10000 に増幅し、
    /// 運用域では許容に埋もれる高次係数（W1/W3 の t⁴ ~1.5e-10 rad、ϖ' の t³ ~6.7e-10 rad）の
    /// 符号反転を delaunay_full / w1 で励起して殺す（章動 t=10 強化と同手法）。
    #[test]
    fn arguments_far_epoch_t10() {
        let a = moon_arguments(10.0);
        let df = [
            77718.96712035326,
            6282.976159171136,
            83289.2857667773,
            84336.28354258099,
        ];
        for (i, (got, exp)) in a.delaunay_full.iter().zip(df.iter()).enumerate() {
            assert_close(*got, *exp, TOL, &format!("t=10 delaunay_full[{i}]"));
        }
        assert_close(a.w1, 84000.65483792205, TOL, "t=10 w1");
        assert_close(a.zeta_lin, 84000.90147965326, TOL, "t=10 zeta_lin");
    }

    /// t≠0 では full と lin が相異（t²..t⁴ が full にのみ加わる）。"full==lin" 実装誤りを殺す。
    #[test]
    fn full_differs_from_lin_at_t1() {
        let a = moon_arguments(1.0);
        let diff_d = (a.delaunay_full[0] - a.delaunay_lin[0]).abs();
        assert!(
            diff_d > 1e-6,
            "t=1: delaunay_full[0] (D) must differ from delaunay_lin[0] (|diff| = {diff_d:e} <= 1e-6); \
             higher-order terms not applied?"
        );
    }
}
