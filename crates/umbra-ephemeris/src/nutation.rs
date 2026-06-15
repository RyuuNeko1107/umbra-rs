//! IAU 2000_R06 章動の評価（ISSUE-035 スライス2）。
//!
//! ISSUE-040 が生成・コミットした packed 係数（`generated/nutation/nut_iau2000_r06.bin`）を
//! `include_bytes!` で取り込み、黄経章動 Δψ・黄道傾斜章動 Δε を評価する。これは eclipse の
//! フレーム連鎖 `iauPnm06a` が用いる 2006 整合章動（SOFA `iauNut06a` 相当）である。
//!
//! 基本引数（Delaunay 5 + 惑星黄経 8 + 一般歳差 1 = 14）の式は IERS Conventions 2003
//! （SOFA `iauFa*03`）の公開多項式。係数表（packed）の出典は
//! `data/coefficient-source/nutation/PROVENANCE.md`。`docs/algorithms/02-frames.md`。

use core::f64::consts::TAU;
use std::sync::OnceLock;
use umbra_core::constants::ARCSEC_TO_RAD;
use umbra_core::{Radians, TtInstant};

/// packed 章動係数（ISSUE-040 生成物、flat little-endian f64）。
/// レイアウト: `[n_psi0, n_psi1, n_eps0, n_eps1, <各項 = [m0..m13, sin_uas, cos_uas]>...]`。
const PACKED: &[u8] = include_bytes!("../../../generated/nutation/nut_iau2000_r06.bin");

/// 基本引数の個数（Delaunay 5 + 惑星黄経 8 + 一般歳差 1）。
const N_ARGS: usize = 14;
/// 1 回転あたりの秒角（Delaunay 引数の約分用）。
const TURNAS: f64 = 1_296_000.0;

/// 章動級数 1 項。乗数は内積で f64 として使うため f64 で保持（整数を無損失格納）。
struct Term {
    multipliers: [f64; N_ARGS],
    sin_amp_uas: f64,
    cos_amp_uas: f64,
}

/// 4 ブロックの章動モデル（packed から一度だけ復号）。
struct NutationModel {
    psi_constant: Vec<Term>,
    psi_rate: Vec<Term>,
    eps_constant: Vec<Term>,
    eps_rate: Vec<Term>,
}

/// packed の `i` 番目の f64（little-endian）。
fn packed_f64(index: usize) -> f64 {
    let start = index * 8;
    let octet: [u8; 8] = PACKED[start..start + 8]
        .try_into()
        .expect("packed nutation length is a multiple of 8 (ISSUE-040 verify-generated gate)");
    f64::from_le_bytes(octet)
}

/// packed ヘッダの項数（自前生成・verify-generated ゲート済みのため非負整数を保証）。
#[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
fn block_count(index: usize) -> usize {
    packed_f64(index).round() as usize
}

/// packed を 4 ブロックへ復号する（一度だけ実行）。
fn model() -> &'static NutationModel {
    static MODEL: OnceLock<NutationModel> = OnceLock::new();
    MODEL.get_or_init(|| {
        let counts = [
            block_count(0),
            block_count(1),
            block_count(2),
            block_count(3),
        ];
        let mut cursor = 4;
        let mut take_block = |count: usize| -> Vec<Term> {
            (0..count)
                .map(|_| {
                    let mut multipliers = [0.0; N_ARGS];
                    for slot in multipliers.iter_mut() {
                        *slot = packed_f64(cursor);
                        cursor += 1;
                    }
                    let sin_amp_uas = packed_f64(cursor);
                    cursor += 1;
                    let cos_amp_uas = packed_f64(cursor);
                    cursor += 1;
                    Term {
                        multipliers,
                        sin_amp_uas,
                        cos_amp_uas,
                    }
                })
                .collect()
        };
        NutationModel {
            psi_constant: take_block(counts[0]),
            psi_rate: take_block(counts[1]),
            eps_constant: take_block(counts[2]),
            eps_rate: take_block(counts[3]),
        }
    })
}

/// 秒角多項式（Horner、最低次から）をラジアンへ。Delaunay 引数は 1 回転で約分してから rad 化
/// （SOFA `iauFa*03` と同じく `fmod(.., TURNAS) * DAS2R`、桁落ち抑制）。
fn delaunay_arcsec(coeffs: &[f64; 5], t: f64) -> f64 {
    let mut acc = 0.0;
    for &c in coeffs.iter().rev() {
        acc = acc * t + c;
    }
    (acc % TURNAS) * ARCSEC_TO_RAD
}

/// 14 基本引数 `[l, l', F, D, Ω, L_Me, L_Ve, L_E, L_Ma, L_J, L_Sa, L_U, L_Ne, p_A]`（rad）。
/// Delaunay = IERS Conventions 2003（SOFA `iauFal03` 等、秒角多項式）。
/// 惑星黄経 = 平均黄経 `a + b·t`（rad、SOFA `iauFame03` 等）を 2π で約分。
/// 一般歳差 p_A = `(0.024381750 + 0.00000538691·t)·t`（rad、SOFA `iauFapa03`）。
pub(crate) fn fundamental_arguments(t: f64) -> [f64; N_ARGS] {
    // Delaunay 引数（秒角係数, IERS Conventions 2003）。
    let l = delaunay_arcsec(
        &[
            485868.249036,
            1717915923.2178,
            31.8792,
            0.051635,
            -0.00024470,
        ],
        t,
    );
    let lp = delaunay_arcsec(
        &[
            1287104.793048,
            129596581.0481,
            -0.5532,
            0.000136,
            -0.00001149,
        ],
        t,
    );
    let f = delaunay_arcsec(
        &[
            335779.526232,
            1739527262.8478,
            -12.7512,
            -0.001037,
            0.00000417,
        ],
        t,
    );
    let d = delaunay_arcsec(
        &[
            1072260.703692,
            1602961601.2090,
            -6.3706,
            0.006593,
            -0.00003169,
        ],
        t,
    );
    let om = delaunay_arcsec(
        &[450160.398036, -6962890.5431, 7.4722, 0.007702, -0.00005939],
        t,
    );
    // 惑星平均黄経（rad、2π で約分）。
    let plan = |a: f64, b: f64| (a + b * t) % TAU;
    let me = plan(4.402608842, 2608.7903141574);
    let ve = plan(3.176146697, 1021.3285546211);
    let ea = plan(1.753470314, 628.3075849991);
    let ma = plan(6.203480913, 334.0612426700);
    let ju = plan(0.599546497, 52.9690962641);
    let sa = plan(0.874016757, 21.3299104960);
    let ur = plan(5.481293872, 7.4781598567);
    let ne = plan(5.311886287, 3.8133035638);
    // 一般歳差（rad、約分しない小角）。
    let pa = (0.024381750 + 0.00000538691 * t) * t;
    [l, lp, f, d, om, me, ve, ea, ma, ju, sa, ur, ne, pa]
}

/// 1 ブロックの級数和（µas）。丸め誤差抑制のため振幅の小さい項（表末尾）から加算する
/// （SOFA `iauNut00a` の reverse 加算に倣う）。
fn sum_block(terms: &[Term], arg: &[f64; N_ARGS]) -> f64 {
    terms
        .iter()
        .rev()
        .map(|term| {
            // 基本引数は各 (-2π,2π) に約分済み（delaunay は mod TURNAS、惑星は mod 2π）。
            // 整数乗数との線形結合 ARG は高々 ~±250 rad に収まり、sin/cos は f64 で十分精確。
            // SOFA iauNut00a は未約分の生引数を使うため最終 fmod を要するが、本実装は引数約分済みで
            // 等価かつ高精度のため最終約分は不要（冗長な mod は入れない）。
            let phase: f64 = term
                .multipliers
                .iter()
                .zip(arg.iter())
                .map(|(m, a)| m * a)
                .sum();
            term.sin_amp_uas * phase.sin() + term.cos_amp_uas * phase.cos()
        })
        .sum()
}

/// IAU 2000_R06 章動 `(Δψ, Δε)`（SOFA `iauNut06a` 相当）。TT 引数。
/// 戻り値はラジアン。係数表は µas 単位（ISSUE-040）→ rad へ変換。
pub fn nutation_iau2006(time_tt: TtInstant) -> (Radians, Radians) {
    let t = time_tt.jd2().julian_centuries_since_j2000();
    let arg = fundamental_arguments(t);
    let model = model();
    let uas_to_rad = 1e-6 * ARCSEC_TO_RAD;
    let dpsi =
        (sum_block(&model.psi_constant, &arg) + t * sum_block(&model.psi_rate, &arg)) * uas_to_rad;
    let deps =
        (sum_block(&model.eps_constant, &arg) + t * sum_block(&model.eps_rate, &arg)) * uas_to_rad;
    (Radians::new(dpsi), Radians::new(deps))
}

#[cfg(test)]
mod tests {
    // SOFA/ERFA 由来の検証値は f64 の表現可能桁数を超える桁で配布される。provenance
    // （一次ソースの逐語転記）を保つため桁を削らず、ここに限り過剰精度リント
    // （余剰桁は最近接 f64 へ丸められ値は不変）を許可する。frames.rs と同様。
    #![allow(clippy::excessive_precision)]

    use super::*;
    use umbra_core::{JulianDate2, TtInstant};

    // ------------------------------------------------------------------
    // 一次オラクル: pyerfa（liberfa C ラッパ）`erfa.nut06a(d1, d2)`。
    //   独立計算経路 = 本実装の Rust（IERS Conventions 2010 IAU 2000_R06 章動系列表の
    //   直接評価）とは別実装（C SOFA/ERFA = nut00a への線形スケーリングで 2006 整合章動）。
    //
    //   取得手順（本テスト設計時に実行）:
    //     docker run --rm python:3.12-slim bash -c \
    //       "pip install -q pyerfa && python3 -c '...erfa.nut06a(d1,d2)...'"
    //   取得環境: erfa 2.0.1.5（pyerfa 同梱の liberfa）。
    //   pyerfa は numpy.float64 を返すため float(...) で Python float 化し repr を転記。
    //
    //   入力 (d1, d2) [TT, 2要素JD] と取得値 (dpsi, deps) [rad]:
    //     J2000      (2451545.0, 0.0)     dpsi=-6.754425598969512e-05  deps=-2.7970831192374137e-05
    //     SOFA-test  (2400000.5, 53736.0) dpsi=-9.630912025821214e-06  deps= 4.063238496887236e-05
    //     1900       (2400000.5, 15020.0) dpsi= 8.452092340677673e-05  deps=-1.1102991495414474e-05
    //     2100       (2400000.5, 88069.0) dpsi= 1.594261371114902e-05  deps= 4.152098077602096e-05
    //     t=10       (2816795.0, 0.0)     dpsi= 6.180294659490065e-05  deps=-3.467512718467922e-05
    // ------------------------------------------------------------------

    /// SOFA の 2要素 JD (date1, date2) を TT 瞬時として構築（frames.rs の `tt` 慣習踏襲）。
    fn tt(date1: f64, date2: f64) -> TtInstant {
        TtInstant::from_jd2(JulianDate2::new(date1, date2))
    }

    /// 許容つきスカラ比較（frames.rs の `close` 踏襲）。
    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    // 1 µas = 1e-6 秒角 = 1e-6 · (π/648000) rad ≈ 4.848e-12 rad。
    const MICROARCSEC_RAD: f64 = 1e-6 * (core::f64::consts::PI / 648_000.0);

    // 運用域（1900–2100, |t|≲1 世紀）の許容 [rad] = 5 µas。
    //
    // 本実装は IERS Conventions 2010 の IAU 2000_R06 章動系列表を**直接評価**する。一方
    // オラクル erfa.nut06a は nut00a への**線形スケーリング近似**で 2006 整合章動を作る。
    // 両者は表現が異なり、その差は章動振幅に比例（実測 ≈ 2e-7 × |Δψ|、すなわち erfa の
    // 線形スケーリングが落とす高次補正）。運用域の最大章動 ~17″≈8.5e-5 rad では差は最大
    // ~3.5 µas（実測: 1900 で Δψ 差 1.7e-11 rad、2100 で <0.7 µas, J2000/SOFA-test で <1 µas）。
    // 5 µas はこの representation 差を吸収しつつ、パース誤り・基本引数式/係数の誤り等の実バグ
    // （残差 ≫ arcsec ≈ 4.85e-6 rad、振幅でなく位相にも効く）は確実に検出する。これは
    // accuracy.md §2.1 のフレーム配分 0.05″=50000 µas の 4 桁下で、精度目標に影響しない。
    const TOL_RAD: f64 = 5.0 * MICROARCSEC_RAD;

    // 域外の極端エポック（t=10 JC ≈ 西暦3000）用の許容 [rad] = 100 µas。
    // representation 差は ∝振幅·（高次）で増え、t=10 では実測 ~47 µas（依然 0.05″ 予算の
    // 3 桁下）。運用域外のため厳密一致は要求せず、gross バグ（位相誤り→arcsec 級）の検出に留める。
    const TOL_FAR_RAD: f64 = 100.0 * MICROARCSEC_RAD;

    /// 単一エポックで erfa.nut06a 値（dpsi, deps）と `tol` 以内で一致することを表明する補助。
    fn assert_nut06a(label: &str, d1: f64, d2: f64, exp_dpsi: f64, exp_deps: f64, tol: f64) {
        let (dpsi, deps) = nutation_iau2006(tt(d1, d2));
        assert!(
            close(dpsi.0, exp_dpsi, tol),
            "{label}: Δψ = {} (expected {exp_dpsi}, |diff| = {:e} > tol {:e})",
            dpsi.0,
            (dpsi.0 - exp_dpsi).abs(),
            tol
        );
        assert!(
            close(deps.0, exp_deps, tol),
            "{label}: Δε = {} (expected {exp_deps}, |diff| = {:e} > tol {:e})",
            deps.0,
            (deps.0 - exp_deps).abs(),
            tol
        );
    }

    // ==================================================================
    // 1. multi-epoch erfa.nut06a 突合
    //    5 エポックで Δψ・Δε を独立計算経路の値と TOL_RAD 以内で照合する。
    //    1900 / 2100 / t=10 を含めることで基本引数の t 依存（Delaunay/惑星黄経/一般歳差の
    //    高次・永年項）を励起し、基本引数式や係数の誤りを検出する。
    // ==================================================================

    #[test]
    fn nut06a_matches_erfa_j2000() {
        // erfa 2.0.1.5  erfa.nut06a(2451545.0, 0.0)。
        assert_nut06a(
            "J2000",
            2451545.0,
            0.0,
            -6.754425598969512e-05,
            -2.7970831192374137e-05,
            TOL_RAD,
        );
    }

    #[test]
    fn nut06a_matches_erfa_sofa_test_epoch() {
        // erfa 2.0.1.5  erfa.nut06a(2400000.5, 53736.0)（T=0.06 JC, SOFA 標準テストエポック）。
        assert_nut06a(
            "SOFA-test",
            2400000.5,
            53736.0,
            -9.630912025821214e-06,
            4.063238496887236e-05,
            TOL_RAD,
        );
    }

    #[test]
    fn nut06a_matches_erfa_1900() {
        // erfa 2.0.1.5  erfa.nut06a(2400000.5, 15020.0)（≈1900, T≈-1 JC: 負の t を励起）。
        assert_nut06a(
            "1900",
            2400000.5,
            15020.0,
            8.452092340677673e-05,
            -1.1102991495414474e-05,
            TOL_RAD,
        );
    }

    #[test]
    fn nut06a_matches_erfa_2100() {
        // erfa 2.0.1.5  erfa.nut06a(2400000.5, 88069.0)（≈2100, T≈+1 JC）。
        assert_nut06a(
            "2100",
            2400000.5,
            88069.0,
            1.594261371114902e-05,
            4.152098077602096e-05,
            TOL_RAD,
        );
    }

    #[test]
    fn nut06a_matches_erfa_far_epoch_t10() {
        // erfa 2.0.1.5  erfa.nut06a(2816795.0, 0.0)（T=10 JC ≈ 西暦3000）。
        // 遠方エポックで基本引数の高次・永年項を強く励起し式の誤りを検出する。
        assert_nut06a(
            "t=10",
            2816795.0,
            0.0,
            6.180294659490065e-05,
            -3.467512718467922e-05,
            TOL_FAR_RAD,
        );
    }

    // ==================================================================
    // 2. オーダー・符号サニティ
    //    章動振幅は Δψ ~ 17″、Δε ~ 9″（~1e-4 rad）。J2000 では上記 erfa 値より
    //    Δψ<0, Δε<0。定数返却・スケール誤り・符号反転を粗く検出する。
    // ==================================================================

    #[test]
    fn nut06a_order_and_sign_at_j2000() {
        let (dpsi, deps) = nutation_iau2006(tt(2451545.0, 0.0));
        // J2000 で双方とも負（erfa 値 dpsi=-6.75e-5, deps=-2.80e-5 の符号）。
        assert!(dpsi.0 < 0.0, "Δψ(J2000) should be < 0, got {}", dpsi.0);
        assert!(deps.0 < 0.0, "Δε(J2000) should be < 0, got {}", deps.0);
        // 章動振幅オーダー ~1e-5〜1e-4 rad（最大振幅 ~17″≈8.2e-5 rad 級）。
        assert!(
            (1e-6..1e-3).contains(&dpsi.0.abs()),
            "|Δψ| out of nutation order: {}",
            dpsi.0.abs()
        );
        assert!(
            (1e-6..1e-3).contains(&deps.0.abs()),
            "|Δε| out of nutation order: {}",
            deps.0.abs()
        );
    }

    // ==================================================================
    // 3. 桁感（暴走検出）: 章動は最大でも ~20″ ≈ 1e-4 rad。全エポックで |Δψ|,|Δε| < 1e-3 rad。
    // ==================================================================

    #[test]
    fn nut06a_magnitude_bounded_across_epochs() {
        const BOUND: f64 = 1e-3; // ~206″。章動最大 ~20″ を 10 倍超で上回るのは暴走。
        for &(d1, d2) in &[
            (2451545.0, 0.0),
            (2400000.5, 53736.0),
            (2400000.5, 15020.0),
            (2400000.5, 88069.0),
            (2816795.0, 0.0),
        ] {
            let (dpsi, deps) = nutation_iau2006(tt(d1, d2));
            assert!(
                dpsi.0.abs() < BOUND,
                "|Δψ| = {} >= {BOUND} at ({d1},{d2})",
                dpsi.0.abs()
            );
            assert!(
                deps.0.abs() < BOUND,
                "|Δε| = {} >= {BOUND} at ({d1},{d2})",
                deps.0.abs()
            );
        }
    }

    // ==================================================================
    // 4. T 依存性（定数返却バグの検出）
    //    異なるエポックで Δψ が実際に変化すること。erfa 値より J2000 と 1900 では
    //    Δψ が符号も大きさも明確に異なる（-6.75e-5 vs +8.45e-5）。
    // ==================================================================

    #[test]
    fn nut06a_varies_with_epoch() {
        let (dpsi_j2000, deps_j2000) = nutation_iau2006(tt(2451545.0, 0.0));
        let (dpsi_1900, deps_1900) = nutation_iau2006(tt(2400000.5, 15020.0));
        // 定数を返す実装なら差は 0。章動の実周期変動で明確に異なるはず（≫ TOL_RAD）。
        assert!(
            (dpsi_j2000.0 - dpsi_1900.0).abs() > 1e-5,
            "Δψ did not change between J2000 and 1900: {} vs {}",
            dpsi_j2000.0,
            dpsi_1900.0
        );
        assert!(
            (deps_j2000.0 - deps_1900.0).abs() > 1e-6,
            "Δε did not change between J2000 and 1900: {} vs {}",
            deps_j2000.0,
            deps_1900.0
        );
    }

    // ==================================================================
    // 5. 基本引数 fundamental_arguments を erfa.fa*03（厳密・独立計算）と直接突合。
    //    nut06a 経由では representation 差・振幅で埋もれる高次係数や約分の誤りを、
    //    引数の生値で直接検出する（係数 sign flip → t=10 で arcsec 級ずれ）。
    //    順序は packed 乗数列順 [l,l',F,D,Ω,Me,Ve,E,Ma,Ju,Sa,Ur,Ne,p_A] と一致。
    //    オラクル: pyerfa 2.0.1.5  erfa.{fal03,falp03,faf03,fad03,faom03,fame03,fave03,
    //      fae03,fama03,faju03,fasa03,faur03,fane03,fapa03}(t)。生値（約分済み）を比較する
    //      ため reduction を入れない（未約分化する `%`→`+` 変異も巨大値として検出される）。
    // ==================================================================
    #[test]
    fn fundamental_arguments_match_erfa() {
        // (t, [14 erfa 値])。erfa.faXX03(t) を Python float 化して転記（erfa 2.0.1.5）。
        let cases: [(f64, [f64; N_ARGS]); 2] = [
            (
                10.0,
                [
                    5.664260211392868,
                    6.074036564661699,
                    3.3700944993124566,
                    2.2477998948536673,
                    -2.375541973345698,
                    4.520355006360376,
                    0.002383433993045969,
                    1.6440131254135508,
                    4.16132419346016,
                    2.502943334914818,
                    0.5448212728940547,
                    4.864668752844963,
                    5.745810081922478,
                    0.244356191,
                ],
            ),
            (
                1.0,
                [
                    5.8266042534984575,
                    6.223481898965822,
                    3.0593179383364326,
                    4.2753563474950695,
                    -0.15864395781697227,
                    5.671020519871966,
                    0.3454962478273771,
                    1.7425245951413046,
                    0.9727169953023775,
                    3.3031603036633115,
                    3.354371331461241,
                    0.3930831143408273,
                    2.8420045436204138,
                    0.02438713691,
                ],
            ),
        ];
        // 引数は厳密（representation 差なし）。実バグ（最小でも高次係数 sign flip ≈ 1e-6 rad）を
        // 確実に検出する水準。本実装と erfa は同一式・同一約分のため正常時は ~1e-13 rad 一致。
        const TOL: f64 = 1e-10;
        for (t, expected) in cases {
            let got = fundamental_arguments(t);
            for (k, (g, e)) in got.iter().zip(expected.iter()).enumerate() {
                assert!(
                    close(*g, *e, TOL),
                    "arg[{k}] at t={t}: {g} vs erfa {e} (|diff| = {:e} > TOL {:e})",
                    (g - e).abs(),
                    TOL
                );
            }
        }
    }
}
