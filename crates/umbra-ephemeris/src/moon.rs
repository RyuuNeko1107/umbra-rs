//! 月 ELP2000-82B 評価（ISSUE-014）。
//!
//! 正本: `docs/algorithms/elp2000-82b-evaluation.md`（著者 Fortran `elp82b_1` から転記、
//! Chapront-Touzé & Chapront, IMCCE MCJCGF.9601）。
//!
//! S1 = 基本引数 [`moon_arguments`]。S2 = 36 系列の総和＋組立＋J2000 回転
//! [`moon_geocentric_j2000`]。

use core::f64::consts::{FRAC_PI_2, TAU};
use std::sync::OnceLock;
use umbra_core::{TdbInstant, Vector3};

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

// ====================================================================
// S2: 36 系列の総和＋組立＋J2000 回転 → 地心直交（km）。
// ====================================================================

/// packed ELP2000-82B 係数（ISSUE-034/S0 生成物、flat little-endian f64）。
/// レイアウト: `[n_series, <系列 = [file, n_mult, n_coeff, n_terms, <項 = n_mult 乗数 + n_coeff 係数>]>...]`。
const PACKED: &[u8] = include_bytes!("../../../generated/elp2000-82b/elp2000-82b_moon.bin");

/// 月慣性質量比関連定数（`elp82b_1`）。
const AM: f64 = 0.074801329518;
const ALFA: f64 = 0.002571881335;
const DTASM: f64 = 2.0 * ALFA / (3.0 * AM);
/// 距離スケール基準（km）。`R = Σ · A0/ATH`。
const A0: f64 = 384747.9806448954;
const ATH: f64 = 384747.9806743165;
/// DE200/LE200 フィット補正（主問題の振幅補正に使用）。RAD/W1(t¹) を含むものは比で rad が相殺。
const W1_RATE: f64 = 1732559343.73604 / RAD;
const DELNU: f64 = (0.55604 / RAD) / W1_RATE;
const DELE: f64 = 0.01789 / RAD;
const DELG: f64 = -0.08066 / RAD;
const DELNP: f64 = (-0.06424 / RAD) / W1_RATE;
const DELEP: f64 = -0.12879 / RAD;
/// 平均黄道 of date → 平均黄道/J2000 慣性分点の回転（Laskar P,Q 多項式, rad）。
const PCOEF: [f64; 5] = [
    0.10180391e-4,
    0.47020439e-6,
    -0.5417367e-9,
    -0.2507948e-11,
    0.463486e-14,
];
const QCOEF: [f64; 5] = [
    -0.113469002e-3,
    0.12372674e-6,
    0.1265417e-8,
    -0.1371808e-11,
    -0.320334e-14,
];

/// 1 項: 整数乗数 + 係数（主問題=[A,B1..B6], 摂動=[φ(度), A]）。
struct PackedTerm {
    multipliers: Vec<i32>,
    coeffs: Vec<f64>,
}

/// 1 系列（ファイル 1..=36）。
struct PackedSeries {
    file: u8,
    terms: Vec<PackedTerm>,
}

/// 検証済み f64 を非負カウントへ（自前生成・verify-generated ゲート済みのため信頼）。
#[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
fn nint(value: f64) -> usize {
    value.round() as usize
}

/// 検証済み f64 を整数乗数 i32 へ。
#[allow(clippy::cast_possible_truncation)]
fn nint_i32(value: f64) -> i32 {
    value.round() as i32
}

/// packed を 1 度だけ復号する（xtask::elp::pack_model と byte-for-byte 対称）。
fn model() -> &'static [PackedSeries] {
    static MODEL: OnceLock<Vec<PackedSeries>> = OnceLock::new();
    MODEL.get_or_init(|| {
        let values: Vec<f64> = PACKED
            .chunks_exact(8)
            .map(|chunk| f64::from_le_bytes(chunk.try_into().expect("packed length multiple of 8")))
            .collect();
        let n_series = nint(values[0]);
        let mut idx = 1;
        let mut series = Vec::with_capacity(n_series);
        for _ in 0..n_series {
            let file = u8::try_from(nint(values[idx])).expect("file fits u8");
            let n_mult = nint(values[idx + 1]);
            let n_coeff = nint(values[idx + 2]);
            let n_terms = nint(values[idx + 3]);
            idx += 4;
            let mut terms = Vec::with_capacity(n_terms);
            for _ in 0..n_terms {
                let multipliers = values[idx..idx + n_mult]
                    .iter()
                    .map(|&v| nint_i32(v))
                    .collect();
                idx += n_mult;
                let coeffs = values[idx..idx + n_coeff].to_vec();
                idx += n_coeff;
                terms.push(PackedTerm {
                    multipliers,
                    coeffs,
                });
            }
            series.push(PackedSeries { file, terms });
        }
        series
    })
}

/// 摂動振幅の時間スケール（Poisson 項）: ×t（/t 系列）/ ×t²（/t² 系列）/ ×1。
fn amplitude_scale(file: u8, t: f64) -> f64 {
    match file {
        7..=9 | 13..=15 | 19..=21 | 25..=27 => t,
        34..=36 => t * t,
        _ => 1.0,
    }
}

/// 月の地心直交座標（mean dynamical ecliptic / inertial equinox of J2000, km）。
/// 補正前の幾何位置（光行時間・章動・光行差は ISSUE-015）。TDB 引数。
pub fn moon_geocentric_j2000(time_tdb: TdbInstant) -> Vector3 {
    let t = time_tdb.jd2().julian_centuries_since_j2000();
    let args = moon_arguments(t);
    let df = &args.delaunay_full;
    let dl = &args.delaunay_lin;
    let zeta = args.zeta_lin;
    let pl = &args.planets_lin;

    // r = [経度(秒角), 緯度(秒角), 距離(km)] の系列総和。
    let mut r = [0.0_f64; 3];
    for s in model() {
        let iv = usize::from((s.file - 1) % 3);
        if s.file <= 3 {
            // 主問題: 振幅は DE200/LE200 補正、引数は完全多項式 Delaunay。距離は cos（+π/2）。
            for term in &s.terms {
                let c = &term.coeffs;
                let tgv = c[1] + DTASM * c[5];
                let mut amp = c[0];
                if s.file == 3 {
                    amp -= 2.0 * amp * DELNU / 3.0;
                }
                let amp =
                    amp + tgv * (DELNP - AM * DELNU) + c[2] * DELG + c[3] * DELE + c[4] * DELEP;
                let mut y = 0.0;
                for (m, d) in term.multipliers.iter().zip(df.iter()) {
                    y += f64::from(*m) * d;
                }
                if iv == 2 {
                    y += FRAC_PI_2;
                }
                r[iv] += amp * (y % TAU).sin();
            }
        } else if (10..=21).contains(&s.file) {
            // 惑星摂動（表1/表2）。乗数 11、係数 [φ, A]。引数は線形。
            let scale = amplitude_scale(s.file, t);
            let table1 = s.file < 16;
            for term in &s.terms {
                let m = &term.multipliers;
                let mut y = term.coeffs[0] * DEG;
                if table1 {
                    // 表1: m[8]·D + m[9]·l + m[10]·F + Σ_{i=0..8} m[i]·planet[i]。
                    y += f64::from(m[8]) * dl[0]
                        + f64::from(m[9]) * dl[2]
                        + f64::from(m[10]) * dl[3];
                    for (mi, p) in m.iter().zip(pl.iter()) {
                        y += f64::from(*mi) * p;
                    }
                } else {
                    // 表2: Σ_{i=0..4} m[i+7]·del[i] + Σ_{i=0..7} m[i]·planet[i]。
                    for (mi, d) in m[7..11].iter().zip(dl.iter()) {
                        y += f64::from(*mi) * d;
                    }
                    for (mi, p) in m[0..7].iter().zip(pl.iter()) {
                        y += f64::from(*mi) * p;
                    }
                }
                r[iv] += term.coeffs[1] * scale * (y % TAU).sin();
            }
        } else {
            // 地球形状/潮汐/月形状/相対論/太陽離心率（4-9, 22-36）。乗数 5 = [iz, ilu0..3]。
            let scale = amplitude_scale(s.file, t);
            for term in &s.terms {
                let m = &term.multipliers;
                let mut y = term.coeffs[0] * DEG + f64::from(m[0]) * zeta;
                for (mi, d) in m[1..5].iter().zip(dl.iter()) {
                    y += f64::from(*mi) * d;
                }
                r[iv] += term.coeffs[1] * scale * (y % TAU).sin();
            }
        }
    }

    // 組立: 経度 V = Σ/rad + W1(t)、緯度 U = Σ/rad、距離 R = Σ·a0/ath。
    let v = r[0] / RAD + args.w1;
    let u = r[1] / RAD;
    let rr = r[2] * A0 / ATH;
    // 球面 → 平均黄道 of date 直交。
    let (sin_u, cos_u) = u.sin_cos();
    let (sin_v, cos_v) = v.sin_cos();
    let x1 = rr * cos_u;
    let x2 = x1 * sin_v;
    let x1 = x1 * cos_v;
    let x3 = rr * sin_u;
    // Laskar P,Q 回転で平均黄道/J2000 慣性分点へ。
    let pw = poly(&PCOEF, t) * t;
    let qw = poly(&QCOEF, t) * t;
    let ra = 2.0 * (1.0 - pw * pw - qw * qw).sqrt();
    let pwqw = 2.0 * pw * qw;
    let pw2 = 1.0 - 2.0 * pw * pw;
    let qw2 = 1.0 - 2.0 * qw * qw;
    let pw = pw * ra;
    let qw = qw * ra;
    Vector3::new(
        pw2 * x1 + pwqw * x2 + pw * x3,
        pwqw * x1 + qw2 * x2 - qw * x3,
        -pw * x1 + qw * x2 + (pw2 + qw2 - 1.0) * x3,
    )
}

#[cfg(test)]
#[allow(clippy::excessive_precision)]
mod tests {
    use super::*;
    use umbra_core::JulianDate2;

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

    // ==================================================================
    // S2: moon_geocentric_j2000 vs オラクル（著者 Fortran elp82b_1 の独立 Python 移植、
    // 公開値 JD 2451555.5 に 5e-10 km 一致＝信頼可）。出力 = mean ecliptic/J2000 慣性分点, km。
    // 1400–2100 を張り t²..t⁴ べき・×t/×t² Poisson 級数・Laskar 回転を励起。
    // 許容 2e-6 km(2 mm): Rust と移植は同一正本アルゴリズムで、差は f64 演算順のみ
    // （Rust は引数を delaunay_full で 1 度評価、移植はべき毎和。実測最大 6.65e-7 km＝0.66mm、
    //  1900 t≈−1.4 で最大。これは 37872 項 f64 総和の固有丸めで月モデル精度 0.40″≈750m の
    //  ~9 桁下）。2e-6 はこの床の ~3 倍で、固定 Docker イメージ上は決定的。これにより構造/
    //  論理バグに加え DE200/LE200 フィット補正（delnu/dele/delg/delnp/delep, ~m 級）と 1 次の
    //  Laskar 回転（pw2/qw2）の誤りも捕捉。残存の微小変異（DTASM 経由 tgv・回転 2 次 pw²・
    //  sin 周期性 %TAU→+TAU）は <床 or 等価で許容（docs/reviews/mutation-moon.md）。
    //  小振幅係数の転記精度は ISSUE-034 の round-trip+checksum が別途保証。
    // ==================================================================

    /// `(JD_TDB, [X_km, Y_km, Z_km])`。
    const ORACLE: [(f64, [f64; 3]); 7] = [
        (
            2305447.5,
            [-133215.08543039399, 359080.0506588529, 9791.334833276025],
        ),
        (
            2415020.5,
            [24466.25480643438, -367508.841233825, 7042.610859853953],
        ),
        (
            2444239.5,
            [43890.1382077410, 381188.7625491393, -31633.3768118953],
        ),
        (
            2451545.0,
            [-291608.23234374094, -274979.94408453995, 36271.17050269245],
        ),
        (
            2451555.5,
            [382979.76047304674, -68204.20174530093, -25987.71602589963],
        ),
        (
            2469000.5,
            [-361602.9853562691, 44996.9951025669, -30696.6531571589],
        ),
        (
            2488069.5,
            [-339520.7371925672, 151145.57646019128, 7061.378184640374],
        ),
    ];

    /// 許容 2e-6 km。実測 Rust↔移植差の最大 6.65e-7 km（1900, t≈−1.4）に ~3 倍の余裕。
    const TOL_KM: f64 = 2e-6;

    /// TDB の単一要素 JD から `TdbInstant` を構築（part2=0）。
    fn tdb(jd: f64) -> TdbInstant {
        TdbInstant::from_jd2(JulianDate2::new(jd, 0.0))
    }

    /// 全オラクル行で X,Y,Z を 1e-3 km 以内照合（系列群引数・距離位相・スケール・回転・単位を殺す）。
    #[test]
    fn moon_geocentric_j2000_matches_oracle() {
        for (jd, expected) in ORACLE {
            let v = moon_geocentric_j2000(tdb(jd));
            let actual = [v.x, v.y, v.z];
            for (i, (a, e)) in actual.iter().zip(expected.iter()).enumerate() {
                let comp = ["X", "Y", "Z"][i];
                assert_close(*a, *e, TOL_KM, &format!("JD{jd} {comp}_km"));
            }
        }
    }

    /// オーダーサニティ: ノルムが物理的な月距離域に入る（単位/スケール暴走を検出）。
    #[test]
    fn moon_distance_order_of_magnitude() {
        for (jd, _) in ORACLE {
            let norm = moon_geocentric_j2000(tdb(jd)).norm();
            assert!(
                (356000.0..407000.0).contains(&norm),
                "lunar distance at JD{jd} out of physical range [356000,407000) km: {norm}"
            );
        }
    }
}
