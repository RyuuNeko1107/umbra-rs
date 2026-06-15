//! 地球姿勢フレーム変換の純式部分（ISSUE-035 のうち章動テーブル不要の部分）。
//!
//! `docs/algorithms/02-frames.md`。CIO ベースで統一する。本スライスは SOFA の
//! `iauObl06` / `iauEra00` / `iauPmat06` / `iauSp00` / `iauPom00` 相当の純式関数のみ。
//! 章動 (IAU2000A nut00a) を要する `gcrs_to_cirs_matrix` / `gcrs_to_itrs_matrix` は別スライス。
//!
//! 各係数は IAU2006 歳差（Capitaine, Wallace & Chapront 2003, A&A 412, 567）/ IAU2000
//! 地球回転角の **公開式**であり、出典を併記する（バルク係数表＝章動 nut00a は別途生成）。
//! SOFA C（`iauObl06`/`iauEra00`/`iauPfw06`/`iauFw2m`/`iauPom00`/`iauSp00`）は参照のみで移植しない。

use umbra_core::constants::{ARCSEC_TO_RAD, J2000_JD};
use umbra_core::{JulianDate2, Matrix3, Radians, TtInstant, Ut1Instant};

use core::f64::consts::TAU;

/// 2要素 JD から TT/TDB ユリウス世紀 T（J2000 基準）を桁落ちを避けて算出する。
/// エポック減算を整数部側で厳密に行う（`JulianDate2::julian_centuries_since_j2000`）。
fn julian_centuries(jd: JulianDate2) -> f64 {
    jd.julian_centuries_since_j2000()
}

/// 秒角単位の多項式 `c[0] + c[1]·t + … + c[n]·tⁿ`（Horner）をラジアンへ変換して返す。
fn arcsec_poly_to_rad(coeffs: &[f64], t: f64) -> f64 {
    let mut acc = 0.0;
    for &c in coeffs.iter().rev() {
        acc = acc * t + c;
    }
    acc * ARCSEC_TO_RAD
}

/// IAU2006 平均黄道傾斜 ε_A。SOFA `iauObl06` 相当。TT 引数。
///
/// `docs/algorithms/02-frames.md`。係数 \[″\]（Capitaine et al. 2003 / IAU2006）:
/// `84381.406 − 46.836769·t − 0.0001831·t² + 0.00200340·t³ − 5.76e-7·t⁴ − 4.34e-8·t⁵`
/// （t = J2000 からの TT ユリウス世紀）。
pub fn mean_obliquity_iau2006(time_tt: TtInstant) -> Radians {
    let t = julian_centuries(time_tt.jd2());
    Radians::new(arcsec_poly_to_rad(
        &[
            84_381.406,
            -46.836_769,
            -0.000_183_1,
            0.002_003_40,
            -0.000_000_576,
            -0.000_000_043_4,
        ],
        t,
    ))
}

/// Earth Rotation Angle（CIO ベース・UT1）。SOFA `iauEra00` 相当。\[0,2π) 正規化。
///
/// `ERA = 2π·(f + 0.7790572732640 + 0.00273781191135448·Tu)`（IERS Conventions 2010, eq. 5.15）。
/// ここで Tu = JD_UT1 − 2451545、f は UT1 の小数日部。巨大 JD と微小係数の積による桁落ちを
/// 避けるため、整数日由来の全回転を `f`（小数日）と分離してから 2π 正規化する（SOFA と同算法）。
fn earth_rotation_angle_rad(time_ut1: Ut1Instant) -> f64 {
    let jd = time_ut1.jd2();
    // 2要素 JD の小数日部の和（C fmod = Rust の `%` と同じく 0 方向丸め）。
    let f = (jd.part1 % 1.0) + (jd.part2 % 1.0);
    // J2000 からの経過日数 Tu（整数部側でエポック減算し桁落ちを避ける）。
    let tu = (jd.part1 - J2000_JD) + jd.part2;
    let theta = TAU * (f + 0.779_057_273_264_0 + 0.002_737_811_911_354_48 * tu);
    Radians::new(theta).normalized_two_pi().0
}

/// Earth Rotation Angle まわりの R3(ERA) 回転行列（CIRS→TIRS）。SOFA `iauEra00`＋R3 相当。
pub fn era_rotation(time_ut1: Ut1Instant) -> Matrix3 {
    Matrix3::rotation_z(earth_rotation_angle_rad(time_ut1))
}

/// frame bias + IAU2006 歳差の合成行列（GCRS→ "of date" 平均赤道系）。
/// SOFA `iauPmat06`（= `iauPfw06` の Fukushima-Williams 角 → `iauFw2m`）相当。TT 引数。
///
/// FW 角の係数 \[″\]（Capitaine et al. 2003 / IAU2006、`docs/algorithms/02-frames.md`）。
/// 合成は `R1(−ε_A)·R3(−ψ_b)·R1(φ_b)·R3(γ_b)`（SOFA `iauFw2m`）。
pub fn bias_precession_matrix_iau2006(time_tt: TtInstant) -> Matrix3 {
    let t = julian_centuries(time_tt.jd2());
    let gamb = arcsec_poly_to_rad(
        &[
            -0.052_928,
            10.556_378,
            0.493_204_4,
            -0.000_312_38,
            -0.000_002_788,
            0.000_000_026_0,
        ],
        t,
    );
    let phib = arcsec_poly_to_rad(
        &[
            84_381.412_819,
            -46.811_016,
            0.051_126_8,
            0.000_532_89,
            -0.000_000_440,
            -0.000_000_017_6,
        ],
        t,
    );
    let psib = arcsec_poly_to_rad(
        &[
            -0.041_775,
            5_038.481_484,
            1.558_417_5,
            -0.000_185_22,
            -0.000_026_452,
            -0.000_000_014_8,
        ],
        t,
    );
    let epsa = mean_obliquity_iau2006(time_tt).0;
    // R = R1(−ε_A)·R3(−ψ_b)·R1(φ_b)·R3(γ_b)（mul_mat は self·other）。
    Matrix3::rotation_x(-epsa)
        .mul_mat(&Matrix3::rotation_z(-psib))
        .mul_mat(&Matrix3::rotation_x(phib))
        .mul_mat(&Matrix3::rotation_z(gamb))
}

/// TIO locator s′。SOFA `iauSp00` 相当。`s' = −47e-6·t` \[″\]（t = TT ユリウス世紀）。
/// 永年ドリフトの線形近似（IERS Conventions 2010）。
fn tio_locator_rad(time_tt: TtInstant) -> f64 {
    let t = julian_centuries(time_tt.jd2());
    arcsec_poly_to_rad(&[0.0, -47.0e-6], t)
}

/// 極運動行列 TIRS→ITRS。SOFA `iauPom00(xp,yp,sp)` 相当。
/// s′（TIO locator）は内部で `iauSp00` 相当式から `time_tt` で算出する。
/// 合成は `R1(−yp)·R2(−xp)·R3(s')`（SOFA `iauPom00`）。
pub fn polar_motion_matrix(xp: Radians, yp: Radians, time_tt: TtInstant) -> Matrix3 {
    let sp = tio_locator_rad(time_tt);
    Matrix3::rotation_x(-yp.0)
        .mul_mat(&Matrix3::rotation_y(-xp.0))
        .mul_mat(&Matrix3::rotation_z(sp))
}

#[cfg(test)]
mod tests {
    // SOFA/ERFA の検証値は f64 の表現可能桁数より多くの桁で配布されている。
    // provenance（一次ソースの逐語転記）を保つため桁を削らず、ここに限り
    // 過剰精度リント（余剰桁は最近接 f64 へ丸められ値は不変）を許可する。
    #![allow(clippy::excessive_precision)]

    use super::*;
    use core::f64::consts::PI;
    use umbra_core::{JulianDate2, Matrix3, TtInstant, Ut1Instant, Vector3};

    // ------------------------------------------------------------------
    // 一次オラクル: ERFA t_erfa_c.c（SOFA 由来・数値同一の検証スイート）。
    //   source: https://github.com/liberfa/erfa  src/t_erfa_c.c (ERFA, master/v2.0.x)
    //   ERFA は SOFA の関数名のみを改名した派生で、IAU 公式 SOFA C の t_sofa_c.c と
    //   入力・期待出力・許容が同一（eraXxx ⇔ iauXxx）。本ライブラリは SOFA 相当式を実装する。
    // 各定数の provenance は使用箇所のコメントに SOFA 関数名・入力・期待値・出典を明記する。
    // ------------------------------------------------------------------

    // SOFA の 2要素 JD は (date1=2400000.5, date2=MJD)。JulianDate2 へそのまま渡す。
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

    /// 行列を期待値（行優先 [行][列]）と要素ごとに比較。
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

    // ==================================================================
    // 1. SOFA 値突合
    // ==================================================================

    /// SOFA iauObl06（ERFA t_obl06）:
    ///   input  : (date1=2400000.5, date2=54388.0)  [TT, 2要素JD]
    ///   expect : 0.4090749229387258204  [rad]
    ///   vvd tol: 1e-14
    ///   source : ERFA t_erfa_c.c / t_obl06; eq. iauObl06 in t_sofa_c.c
    ///
    /// なぜ TT か: iauObl06 は T = ((date1−2451545)+date2)/36525 [JC, TT] のみに依存する純多項式。
    /// 引数を UT1 と取り違えると本検証値からずれて検出できる。
    /// 許容: SOFA の vvd と同等の 1e-14（純多項式評価で十分達成可能）。
    #[test]
    fn obl06_matches_sofa() {
        let eps = mean_obliquity_iau2006(tt(2400000.5, 54388.0));
        assert!(
            close(eps.0, 0.4090749229387258204, 1e-14),
            "obl06 = {} (expected 0.4090749229387258204)",
            eps.0
        );
    }

    /// SOFA iauEra00（ERFA t_era00）:
    ///   input  : (dj1=2400000.5, dj2=54388.0)  [UT1, 2要素JD]
    ///   expect ERA : 0.4022837240028158102  [rad]
    ///   vvd tol: 1e-12
    ///   source : ERFA t_erfa_c.c / t_era00
    ///
    /// era_rotation は ERA まわりの R3 回転行列を返す（CIRS→TIRS）。
    /// よって期待行列は Matrix3::rotation_z(ERA)（R3 受動回転、umbra-core matrix.rs と同一規約）。
    /// 期待 ERA 値は SOFA の 2要素フラクショナルデイ算法による高精度値。
    ///
    /// 許容に関する設計メモ（精度クリティカル）:
    ///   ERA = 2π(0.7790572732640 + 1.00273781191135448·(JD_UT1−2451545)) を素朴に
    ///   単一 JD で評価すると、巨大 JD と微小係数の積で ~2e-12 rad の桁落ちが出て
    ///   期待値と 9〜10 有効数字目でずれる。SOFA は f=fmod(d1,1)+fmod(d2,1) の
    ///   フラクショナルデイ分割で評価しこれを回避している。本検証は SOFA の vvd と同じ
    ///   絶対許容 1e-12 を採用する（2要素算法なら ~1e-13、素朴式でもぎりぎり 1e-12 内）。
    #[test]
    fn era00_matches_sofa() {
        const ERA_EXPECT: f64 = 0.4022837240028158102;
        let m = era_rotation(ut1(2400000.5, 54388.0));
        let expected = Matrix3::rotation_z(ERA_EXPECT);
        assert!(
            matrix_close(&m, &expected.rows, 1e-12),
            "era_rotation != R3({ERA_EXPECT}); got {:?}",
            m.rows
        );
    }

    /// SOFA iauPmat06（ERFA t_pmat06）:
    ///   input  : (date1=2400000.5, date2=50123.9999)  [TT, 2要素JD]
    ///   expect rbp[i][j]（行優先）:
    ///     [ 0.9999995505176007047,  0.8695404617348208406e-3,  0.3779735201865589104e-3]
    ///     [-0.8695404723772031414e-3, 0.9999996219496027161,  -0.1361752497080270143e-6]
    ///     [-0.3779734957034089490e-3, -0.1924880847894457113e-6, 0.9999999285679971958]
    ///   vvd tol: 対角 1e-12 / 非対角 1e-14
    ///   source : ERFA t_erfa_c.c / t_pmat06
    ///
    /// なぜ TT か: bias-precession は T[JC,TT] の多項式（iauPfw06）に依存。UT1 取り違えで外れる。
    /// 許容: SOFA の最小許容（非対角 1e-14）を全要素に課す保守設定（純式実装で達成可能）。
    #[test]
    fn pmat06_matches_sofa() {
        const RBP: [[f64; 3]; 3] = [
            [
                0.9999995505176007047,
                0.8695404617348208406e-3,
                0.3779735201865589104e-3,
            ],
            [
                -0.8695404723772031414e-3,
                0.9999996219496027161,
                -0.1361752497080270143e-6,
            ],
            [
                -0.3779734957034089490e-3,
                -0.1924880847894457113e-6,
                0.9999999285679971958,
            ],
        ];
        let m = bias_precession_matrix_iau2006(tt(2400000.5, 50123.9999));
        assert!(
            matrix_close(&m, &RBP, 1e-14),
            "pmat06 mismatch; got {:?}",
            m.rows
        );
    }

    /// SOFA iauPom00 + iauSp00（ERFA t_pom00, t_sp00）:
    ///   polar_motion_matrix は内部で s′ を iauSp00 相当式から time_tt で算出する。
    ///   iauSp00:  s' = -47e-6 · T[JC,TT] · (π/648000)  [rad]
    ///   t_pom00 の検証は xp=2.55060238e-7, yp=1.860359247e-6,
    ///     sp=-0.1367174580728891460e-10 を直接渡したもの。
    ///   この sp は s' = -47e-6·T·DAS2R を逆算すると T=0.06 JC、すなわち
    ///     (date1=2400000.5, date2=53736.0) [TT] に一致する（独自検算で確認）。
    ///   したがって polar_motion_matrix(xp, yp, TT(2400000.5,53736.0)) は
    ///   t_pom00 の期待 rpom を再現しなければならない。
    ///   expect rpom[i][j]（行優先）:
    ///     [ 0.9999999999999674721,  -0.1367174580728846989e-10, 0.2550602379999972345e-6]
    ///     [ 0.1414624947957029801e-10, 0.9999999999982695317,  -0.1860359246998866389e-5]
    ///     [-0.2550602379741215021e-6,  0.1860359247002414021e-5, 0.9999999999982370039]
    ///   vvd tol: 対角 1e-12 / 非対角 1e-16
    ///   source : ERFA t_erfa_c.c / t_pom00, t_sp00 (sp00.c: -47e-6·T·DAS2R)
    ///
    /// なぜ s′ が TT か: iauSp00 は T[JC,TT] に比例（µas/世紀の永年ドリフト近似）。
    /// 許容: SOFA の最小許容（非対角 1e-16）を全要素に課す保守設定。
    #[test]
    fn pom00_with_internal_sp_matches_sofa() {
        const XP: f64 = 2.55060238e-7;
        const YP: f64 = 1.860359247e-6;
        const RPOM: [[f64; 3]; 3] = [
            [
                0.9999999999999674721,
                -0.1367174580728846989e-10,
                0.2550602379999972345e-6,
            ],
            [
                0.1414624947957029801e-10,
                0.9999999999982695317,
                -0.1860359246998866389e-5,
            ],
            [
                -0.2550602379741215021e-6,
                0.1860359247002414021e-5,
                0.9999999999982370039,
            ],
        ];
        // T=0.06 JC ⇒ date2 = 0.06*36525 + 2451545 - 2400000.5 = 53736.0。
        let m = polar_motion_matrix(Radians::new(XP), Radians::new(YP), tt(2400000.5, 53736.0));
        assert!(
            matrix_close(&m, &RPOM, 1e-16),
            "pom00 (internal sp) mismatch; got {:?}",
            m.rows
        );
    }

    /// iauSp00 値の単独突合（s′ が time_tt から正しく算出されることの直接確認）。
    /// polar_motion_matrix は s′ を公開しないため、s′ を R3(sp) 成分として持つ
    /// rpom の [0][1]/[1][0] 近傍で検証する代わりに、sp が小さい極限での rpom の
    /// 解析形（xp=yp=0）から s′ を取り出して突合する。
    ///   xp=yp=0 のとき rpom = R3(sp) なので rpom[0][1] = sin(sp) ≈ sp。
    ///   SOFA iauSp00(2400000.5, 52541.0) = -0.6216698469981019309e-11 [rad]
    ///     (ERFA t_sp00; sp00.c: s' = -47e-6·T·DAS2R)
    ///   この入力の date2=52541.0 を TT に与え、rpom[0][1] ≈ sin(sp) を突合。
    #[test]
    fn internal_sp00_matches_sofa_via_zero_polar_motion() {
        // iauSp00(2400000.5, 52541.0) expected:
        const SP_EXPECT: f64 = -0.6216698469981019309e-11;
        let m = polar_motion_matrix(Radians::new(0.0), Radians::new(0.0), tt(2400000.5, 52541.0));
        // xp=yp=0 ⇒ rpom = R3(sp): rows[0][1] = sin(sp) ≈ sp (|sp|~6e-12)。
        assert!(
            close(m.rows[0][1], SP_EXPECT.sin(), 1e-23),
            "internal sp (via rpom[0][1]) = {} (expected sin(sp)={})",
            m.rows[0][1],
            SP_EXPECT.sin()
        );
        // R3 形であること: rows[2][2] == 1, rows[2][0]=rows[2][1]=0。
        assert!(close(m.rows[2][2], 1.0, 1e-15));
        assert!(close(m.rows[2][0], 0.0, 1e-15));
        assert!(close(m.rows[2][1], 0.0, 1e-15));
    }

    // ==================================================================
    // 1b. ミューテーション堅牢化（遠方エポックで高次係数を励起・分割不変性）
    // ==================================================================
    //
    // SOFA 検証エポックは t≈0.04〜0.08 JC で、歳差・黄道傾斜の t⁴/t⁵ 係数の
    // 寄与が許容（~1e-14）未満になり、それら係数の符号取り違えを単一エポックでは
    // 検出できない。遠方エポック t=10 JC（西暦 ~3000）では t⁵ が 1e5 倍に増幅され
    // 高次係数を実効的に検証できる。これは多項式評価の単体検証であり、独立計算経路
    // （pyerfa = liberfa C ラッパ。本実装の Rust とは別実装）の値をオラクルに用いる。

    /// 遠方エポック t=10 JC（JD 2816795.0 = J2000+10 世紀）での黄道傾斜。
    /// 高次係数(t⁴,t⁵)を励起し符号取り違えを検出。
    /// オラクル: pyerfa 2.0.1.5  `erfa.obl06(2816795.0, 0.0) = 0.40683146498328676` \[rad\]。
    #[test]
    fn obl06_far_epoch_exercises_high_order_terms() {
        let eps = mean_obliquity_iau2006(tt(2816795.0, 0.0));
        assert!(
            close(eps.0, 0.40683146498328676, 1e-12),
            "obl06(t=10) = {} (expected 0.40683146498328676)",
            eps.0
        );
    }

    /// 遠方エポック t=10 JC での bias-precession 行列。FW角(gamb/phib/psib)の
    /// t⁴/t⁵ 係数を励起し符号取り違えを検出。
    /// オラクル: pyerfa 2.0.1.5  `erfa.pmat06(2816795.0, 0.0)`。
    #[test]
    fn pmat06_far_epoch_exercises_high_order_terms() {
        const RBP: [[f64; 3]; 3] = [
            [
                0.9702976997703714,
                -0.22205299729851868,
                -0.09599395924260168,
            ],
            [
                0.22205132252351573,
                0.9749747654752644,
                -0.010835905684833613,
            ],
            [
                0.09599783323535901,
                -0.010801531243165408,
                0.9953229339952532,
            ],
        ];
        let m = bias_precession_matrix_iau2006(tt(2816795.0, 0.0));
        assert!(
            matrix_close(&m, &RBP, 1e-12),
            "pmat06(t=10) mismatch; got {:?}",
            m.rows
        );
    }

    /// ERA は総 JD が同じなら part1/part2 の分割に依存しない（小数日部の合成が正しいことの検証）。
    /// 既存 SOFA 突合は整数 MJD（`part2 % 1.0 = 0`）のため `f = (part1%1)+(part2%1)` の
    /// part2 小数項を励起しない。本テストは part2 に小数を持たせ、同一総 JD の整数分割と
    /// 一致することを要求する（小数日合成の `+` を取り違えると分割で値が割れて検出される）。
    #[test]
    fn era_invariant_under_jd_split() {
        // 総 JD = 2454388.75 を 2 通りに分割（A: part2 に小数 0.75、B: 整数 part1）。
        let split_a = era_rotation(ut1(2451545.0, 2843.75));
        let split_b = era_rotation(ut1(2454388.75, 0.0));
        assert!(
            matrix_close(&split_a, &split_b.rows, 1e-13),
            "ERA differs across JD split: A={:?} B={:?}",
            split_a.rows,
            split_b.rows
        );
        // 念のため自明でない回転であること（恒等行列との比較で空テスト化を防ぐ）。
        assert!(
            !matrix_close(&split_a, &Matrix3::IDENTITY.rows, 1e-3),
            "ERA unexpectedly near identity"
        );
    }

    // ==================================================================
    // 2. 性質テスト
    // ==================================================================

    /// era_rotation は直交・det≈+1（回転行列）。
    #[test]
    fn era_rotation_is_orthonormal() {
        let m = era_rotation(ut1(2400000.5, 54388.0));
        assert_orthonormal(&m, 1e-12);
        assert!(close(det(&m), 1.0, 1e-12), "det = {}", det(&m));
    }

    /// era_rotation は R3(ERA): z 軸不変、左上 2×2 が回転、下端行/右端列が単位。
    #[test]
    fn era_rotation_is_r3_about_z() {
        let m = era_rotation(ut1(2400000.5, 54388.0));
        // z 軸ベクトルは不変。
        let z = m.mul_vec(Vector3::new(0.0, 0.0, 1.0));
        assert!(close(z.x, 0.0, 1e-14) && close(z.y, 0.0, 1e-14) && close(z.z, 1.0, 1e-14));
        // R3 形の固定要素。
        assert!(close(m.rows[2][2], 1.0, 1e-14));
        assert!(close(m.rows[0][2], 0.0, 1e-14));
        assert!(close(m.rows[1][2], 0.0, 1e-14));
        assert!(close(m.rows[2][0], 0.0, 1e-14));
        assert!(close(m.rows[2][1], 0.0, 1e-14));
        // R3 の符号規約（rows[0][1] = +sinθ, rows[1][0] = −sinθ）。
        assert!(close(m.rows[0][1], -m.rows[1][0], 1e-14));
        assert!(close(m.rows[0][0], m.rows[1][1], 1e-14));
    }

    /// ERA の [0,2π) 正規化: UT1 を大きくずらしても回転角が範囲内（行列が常に有効な回転）。
    #[test]
    fn era_rotation_normalized_over_large_ut1_range() {
        for &mjd in &[40000.0, 51544.5, 60000.0, 88069.0, 100000.0] {
            let m = era_rotation(ut1(2400000.5, mjd));
            // cos/sin が [-1,1] にあり直交であること = 正規化済み有効回転。
            assert_orthonormal(&m, 1e-11);
            let c = m.rows[0][0];
            let s = m.rows[0][1];
            assert!((-1.0..=1.0).contains(&c), "cos out of range: {c}");
            assert!((-1.0..=1.0).contains(&s), "sin out of range: {s}");
            // 角度を atan2 で復元し [0,2π) に入ること。
            let theta = s.atan2(c).rem_euclid(2.0 * PI);
            assert!((0.0..2.0 * PI).contains(&theta), "theta = {theta}");
        }
    }

    /// bias_precession_matrix_iau2006 は直交・det≈+1。
    #[test]
    fn bias_precession_is_orthonormal() {
        let m = bias_precession_matrix_iau2006(tt(2400000.5, 50123.9999));
        assert_orthonormal(&m, 1e-12);
        assert!(close(det(&m), 1.0, 1e-12), "det = {}", det(&m));
    }

    /// polar_motion_matrix は直交・det≈+1。
    #[test]
    fn polar_motion_is_orthonormal() {
        let m = polar_motion_matrix(
            Radians::new(2.55060238e-7),
            Radians::new(1.860359247e-6),
            tt(2400000.5, 53736.0),
        );
        assert_orthonormal(&m, 1e-12);
        assert!(close(det(&m), 1.0, 1e-12), "det = {}", det(&m));
    }

    /// mean_obliquity_iau2006 は J2000(TT) で ≈ 23°26′21.4″ ≈ 0.4090928 rad（オーダーサニティ）。
    /// 厳密値は obl06_matches_sofa で別途担保。ここでは桁・近傍のみ確認。
    #[test]
    fn obliquity_at_j2000_is_about_23_4_deg() {
        // J2000.0 TT = (2451545.0, 0.0)。
        let eps = mean_obliquity_iau2006(tt(2451545.0, 0.0));
        // iauObl06 の定数項 84381.406″ → 0.40909260060... rad。
        assert!(
            close(eps.0, 0.4090926006005829, 1e-9),
            "eps(J2000) = {} rad (expected ~0.4090926)",
            eps.0
        );
        // 度に直して 23.4393° 近傍。
        let deg = eps.0 * 180.0 / PI;
        assert!(close(deg, 23.439_28, 1e-3), "eps = {deg} deg");
    }

    // ==================================================================
    // 3. 時刻系取り違え回帰
    // ==================================================================

    /// 時刻系取り違えガード（ERA=UT1 / obliquity=TT の identity を交差確認）。
    ///
    /// SOFA の era00 と obl06 は同じ 2要素 JD (2400000.5, 54388.0) に対し
    /// **全く異なる量**を返す（ERA=0.40228372..., obl=0.40907492...）。
    /// 実装が UT1/TT を取り違えても、各関数は受け取った Instant 型で固定されるため
    /// コンパイル時に弾かれるが、内部で別の時刻量（例: ERA に T[JC] を使う等）を
    /// 取り違えた場合は SOFA 突合で外れる。本テストは両者が一致しないこと自体を
    /// 明示し、取り違え実装が双方の SOFA 値を同時に満たすことは不可能であることを示す。
    #[test]
    fn era_and_obliquity_differ_at_same_jd() {
        let m_era = era_rotation(ut1(2400000.5, 54388.0));
        // ERA を R3 角として復元。
        let era_angle = m_era.rows[0][1]
            .atan2(m_era.rows[0][0])
            .rem_euclid(2.0 * PI);
        let obl = mean_obliquity_iau2006(tt(2400000.5, 54388.0)).0;
        // 0.4022837 vs 0.4090749 — 約 6.8e-3 rad 異なる。取り違え実装は両立不能。
        assert!(
            (era_angle - obl).abs() > 1e-3,
            "ERA ({era_angle}) and obliquity ({obl}) unexpectedly equal — possible time-scale swap"
        );
    }
}
