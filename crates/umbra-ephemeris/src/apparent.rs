//! 見かけ地心位置（ISSUE-015）。
//!
//! S1 = 幾何地心位置のフレーム調和: 黄道座標の暦（太陽 VSOP87D = 黄道 of date、
//! 月 ELP2000-82B = 黄道 J2000）を共通の **GCRS**（ICRS 軸・地心）へ載せる（`*_geocentric_gcrs`）。
//! S2（本コミット）= 光行時間補正（`*_light_time_corrected_gcrs`, SOFA `iauAtciq` の light-time
//! ステップ相当）: 天体 = 放射時刻 t−τ・観測者 = 観測時刻 t を一貫させた幾何地心ベクトルを返す。
//! 光行差（S3, `iauAb`）→ 歳差章動（S4, GCRS→CIRS）は後続スライス（SOFA `iauAtciq` 順, ISSUE-015）。
//!
//! 出力フレーム = GCRS（`docs/issues/ISSUE-015` 確定 / iauAtciq 標準）。入力は `TtInstant`
//! （位置計算標準, conventions §6）。暦評価は TT≈TDB 近似（差 ≲2ms, metadata 帰属外）。
//!
//! 注（ELP の J2000 慣性分点）: ELP2000-82B の出力は「平均力学黄道・**慣性**分点 of J2000」。
//! 本スライスはこれを IAU2006 の J2000 平均黄道・分点と同一視し、J2000 での `ecliptic_to_gcrs_matrix`
//! （= ERFA `ecm06(J2000)ᵀ`、frame bias + J2000 黄道傾斜を含む）で GCRS へ回す。慣性 vs 回転分点の
//! 微小オフセット（~0.1″）は既知の近似で、M10 の JPL DE 差分で確定する（R06 章動と同じ実用判断）。

use crate::frames::ecliptic_to_gcrs_matrix;
use crate::moon::moon_geocentric_j2000;
use crate::sun::sun_geocentric_ecliptic_of_date;
use umbra_core::constants::{ASTRONOMICAL_UNIT_KM, J2000_JD};
use umbra_core::{JulianDate2, TdbInstant, TtInstant, Vector3};

/// 太陽の幾何地心位置（GCRS, km）。補正前（光行時間・光行差は後続）。TT 入力。
/// VSOP87D（黄道 of date, AU）を **観測日**の黄道→GCRS 行列で回転し km 化する。
pub fn sun_geocentric_gcrs(time_tt: TtInstant) -> Vector3 {
    let tdb = TdbInstant::from_jd2(time_tt.jd2());
    let ecl_km = sun_geocentric_ecliptic_of_date(tdb).scale(ASTRONOMICAL_UNIT_KM);
    ecliptic_to_gcrs_matrix(time_tt).mul_vec(ecl_km)
}

/// 月の幾何地心位置（GCRS, km）。補正前。TT 入力。
/// ELP2000-82B（黄道 J2000, km）を **J2000** の黄道→GCRS 行列で回転する（暦が J2000 黄道系のため
/// 行列は観測日でなく J2000 固定）。
pub fn moon_geocentric_gcrs(time_tt: TtInstant) -> Vector3 {
    let tdb = TdbInstant::from_jd2(time_tt.jd2());
    let ecl_j2000 = moon_geocentric_j2000(tdb);
    let m = ecliptic_to_gcrs_matrix(TtInstant::from_jd2(JulianDate2::new(J2000_JD, 0.0)));
    m.mul_vec(ecl_j2000)
}

/// 光行時間補正後の幾何地心位置（GCRS, km）＋収束した光行時間 τ（秒）。
#[derive(Debug, Clone, Copy)]
pub struct LightTimeCorrected {
    /// SOFA「astrometric」幾何ベクトル: 天体 = 放射時刻 t−τ、観測者 = 観測時刻 t（GCRS, km）。
    pub position_gcrs: Vector3,
    /// 収束した光行時間 τ（秒）。
    pub light_time_seconds: f64,
}

/// 光行時間補正の本体（SOFA `iauAtciq` の light-time ステップ相当）。
///
/// 観測時刻 `t` における幾何地心ベクトル `B_geo(t') = 天体(t') − 地球(t')`（GCRS, km）を
/// `body_gcrs` が供給する。出力は**天体 = 放射時刻 t−τ・観測者 = 観測時刻 t** の幾何ベクトル
///
/// ```text
///   s = B_geo(t−τ) + ( E(t−τ) − E(t) )      （一次近似 E(t−τ)−E(t) ≈ −v_E·τ）
/// ```
///
/// 第2項は地球が光行時間 τ の間に動いた変位で、角度にして約 `v_E/c ≈ 20.5″`（距離によらず一定）。
/// これを含めることで `s` が SOFA「astrometric」幾何ベクトル（地球運動分を除いた、後段の恒星光行差
/// `iauAb` が乗る前の量）になり、S3 で純粋な光行差を当てても二重計上にならない（ISSUE-015 D3）。
/// 素朴に `B_geo(t−τ)` だけを返すと地球運動分が混ざり、後段光行差と数十″の誤差を生む。
///
/// `v_E` は VSOP87D 解析微分（黄道 of date, km/s）を観測日行列で GCRS へ回したもの。曲率項
/// `½·a·τ² ≈ 0.7 km`（太陽, ≈0.001″）は予算外として一次近似で省略する。
fn light_time_correct(
    time_tt: TtInstant,
    body_gcrs: impl Fn(TtInstant) -> Vector3,
) -> LightTimeCorrected {
    let c = umbra_core::constants::SPEED_OF_LIGHT_KM_S;
    let tdb = TdbInstant::from_jd2(time_tt.jd2());
    let v_e_gcrs = ecliptic_to_gcrs_matrix(time_tt).mul_vec(
        crate::sun::earth_heliocentric_velocity_ecliptic_of_date(tdb),
    );
    let p0 = body_gcrs(time_tt);
    let mut tau = p0.norm() / c;
    let mut position = p0;
    for _ in 0..5 {
        let emit = TtInstant::from_jd2(time_tt.jd2().add_days(-tau / 86400.0));
        let s = body_gcrs(emit) + v_e_gcrs.scale(-tau);
        let next = s.norm() / c;
        position = s;
        let converged = (next - tau).abs() < 1e-6;
        tau = next;
        if converged {
            break;
        }
    }
    LightTimeCorrected {
        position_gcrs: position,
        light_time_seconds: tau,
    }
}

/// 太陽の光行時間補正後の幾何地心位置（GCRS, km）と τ。
pub fn sun_light_time_corrected_gcrs(time_tt: TtInstant) -> LightTimeCorrected {
    light_time_correct(time_tt, sun_geocentric_gcrs)
}

/// 月の光行時間補正後の幾何地心位置（GCRS, km）と τ。
pub fn moon_light_time_corrected_gcrs(time_tt: TtInstant) -> LightTimeCorrected {
    light_time_correct(time_tt, moon_geocentric_gcrs)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// 観測時刻 TtInstant を 1要素 JD（小数日 0）から構築。
    fn tt(jd: f64) -> TtInstant {
        TtInstant::from_jd2(JulianDate2::new(jd, 0.0))
    }

    /// TT≈TDB として同 JD の TdbInstant を構築（暦呼び出し用）。
    fn tdb(jd: f64) -> TdbInstant {
        TdbInstant::from_jd2(JulianDate2::new(jd, 0.0))
    }

    /// 許容つきスカラ比較（clippy::float_cmp 回避）。
    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    /// 2ベクトルの要素ごと一致（合成同一性検証用）。
    fn vec_close(a: Vector3, b: Vector3, tol: f64) -> bool {
        close(a.x, b.x, tol) && close(a.y, b.y, tol) && close(a.z, b.z, tol)
    }

    /// J2000 における黄道→GCRS 行列（月変換の固定行列）。
    fn matrix_at_j2000() -> umbra_core::Matrix3 {
        ecliptic_to_gcrs_matrix(tt(J2000_JD))
    }

    // 太陽: sun_gcrs = ecliptic_to_gcrs_matrix(time_tt) · (黄道ベクトル×AU)。回転はノルム不変。

    /// (b) 合成同一性: 戻り値 = 観測日の行列 · (黄道×AU)。行列日・スケール・順序の取り違えを殺す。
    #[test]
    fn sun_gcrs_equals_observation_date_matrix_times_ecliptic() {
        for &jd in &[J2000_JD, 2469807.0] {
            let ecl_km = sun_geocentric_ecliptic_of_date(tdb(jd)).scale(ASTRONOMICAL_UNIT_KM);
            let expected = ecliptic_to_gcrs_matrix(tt(jd)).mul_vec(ecl_km);
            let got = sun_geocentric_gcrs(tt(jd));
            assert!(
                vec_close(got, expected, 1e-3),
                "sun_gcrs(jd={jd}) = {got:?}, expected {expected:?}"
            );
        }
    }

    /// (a) ノルム保存: |sun_gcrs| == |黄道×AU|。
    #[test]
    fn sun_gcrs_preserves_norm() {
        for &jd in &[J2000_JD, 2469807.0] {
            let ecl_norm = sun_geocentric_ecliptic_of_date(tdb(jd))
                .scale(ASTRONOMICAL_UNIT_KM)
                .norm();
            let gcrs_norm = sun_geocentric_gcrs(tt(jd)).norm();
            assert!(
                close(gcrs_norm, ecl_norm, ecl_norm * 1e-6),
                "sun norm(jd={jd}): gcrs={gcrs_norm}, ecliptic={ecl_norm}"
            );
        }
    }

    /// (c) オーダーサニティ: 太陽地心距離 1.4e8〜1.6e8 km（≈1 AU）。
    #[test]
    fn sun_gcrs_distance_order_of_magnitude() {
        for &jd in &[J2000_JD, 2469807.0] {
            let r = sun_geocentric_gcrs(tt(jd)).norm();
            assert!(
                (1.4e8..1.6e8).contains(&r),
                "sun distance(jd={jd}) = {r} km out of [1.4e8, 1.6e8]"
            );
        }
    }

    // 月: moon_gcrs = ecliptic_to_gcrs_matrix(J2000) · (月 J2000 黄道)。行列は常に J2000。

    /// (b) 合成同一性: 戻り値 = J2000 の行列 · (月 J2000 黄道)。行列を観測日で取る誤りを殺す。
    #[test]
    fn moon_gcrs_equals_j2000_matrix_times_ecliptic() {
        for &jd in &[J2000_JD, 2469807.0] {
            let expected = matrix_at_j2000().mul_vec(moon_geocentric_j2000(tdb(jd)));
            let got = moon_geocentric_gcrs(tt(jd));
            assert!(
                vec_close(got, expected, 1e-6),
                "moon_gcrs(jd={jd}) = {got:?}, expected {expected:?}"
            );
        }
    }

    /// (a) ノルム保存: |moon_gcrs| == |月 J2000 黄道|。
    #[test]
    fn moon_gcrs_preserves_norm() {
        for &jd in &[J2000_JD, 2469807.0] {
            let ecl_norm = moon_geocentric_j2000(tdb(jd)).norm();
            let gcrs_norm = moon_geocentric_gcrs(tt(jd)).norm();
            assert!(
                close(gcrs_norm, ecl_norm, ecl_norm * 1e-6),
                "moon norm(jd={jd}): gcrs={gcrs_norm}, ecliptic={ecl_norm}"
            );
        }
    }

    /// (c) オーダーサニティ: 月地心距離 356000〜407000 km。
    #[test]
    fn moon_gcrs_distance_order_of_magnitude() {
        for &jd in &[J2000_JD, 2469807.0] {
            let r = moon_geocentric_gcrs(tt(jd)).norm();
            assert!(
                (356_000.0..407_000.0).contains(&r),
                "moon distance(jd={jd}) = {r} km out of [356000, 407000]"
            );
        }
    }

    // ============================================================
    // S2: 光行時間補正 (sun/moon_light_time_corrected_gcrs)
    // ============================================================

    use umbra_core::constants::SPEED_OF_LIGHT_KM_S;

    /// 2ベクトル間の角（rad）。微小角でも acos のクランプで安定。
    fn angle_between(a: Vector3, b: Vector3) -> f64 {
        let c = a.dot(b) / (a.norm() * b.norm());
        c.clamp(-1.0, 1.0).acos()
    }

    /// 秒角 → ラジアン。
    fn arcsec_to_rad(s: f64) -> f64 {
        s * std::f64::consts::PI / (180.0 * 3600.0)
    }

    /// 出力 τ を所与に t−τ を作る（仕様: add_days(-τ/86400)）。
    fn retarded(time_tt: TtInstant, tau_s: f64) -> TtInstant {
        TtInstant::from_jd2(time_tt.jd2().add_days(-tau_s / 86400.0))
    }

    /// 第2項 −v_E·τ を GCRS で組む（一次近似; 観測日行列で回転）。
    fn earth_displacement_term_gcrs(time_tt: TtInstant, tau_s: f64) -> Vector3 {
        let tdb_now = TdbInstant::from_jd2(time_tt.jd2());
        let v_ecl = crate::sun::earth_heliocentric_velocity_ecliptic_of_date(tdb_now);
        ecliptic_to_gcrs_matrix(time_tt)
            .mul_vec(v_ecl)
            .scale(-tau_s)
    }

    // ---- A. 合成同一性 ----

    #[test]
    fn sun_light_time_corrected_equals_definition_formula() {
        for &jd in &[J2000_JD, 2469807.0] {
            let out = sun_light_time_corrected_gcrs(tt(jd));
            let tau = out.light_time_seconds;
            let b_retarded = sun_geocentric_gcrs(retarded(tt(jd), tau));
            let term = earth_displacement_term_gcrs(tt(jd), tau);
            let expected = b_retarded + term;
            assert!(
                vec_close(out.position_gcrs, expected, 50.0),
                "sun s(jd={jd}) = {:?}, expected {:?}",
                out.position_gcrs,
                expected
            );
        }
    }

    #[test]
    fn moon_light_time_corrected_equals_definition_formula() {
        for &jd in &[J2000_JD, 2469807.0] {
            let out = moon_light_time_corrected_gcrs(tt(jd));
            let tau = out.light_time_seconds;
            let b_retarded = moon_geocentric_gcrs(retarded(tt(jd), tau));
            let term = earth_displacement_term_gcrs(tt(jd), tau);
            let expected = b_retarded + term;
            assert!(
                vec_close(out.position_gcrs, expected, 1.0),
                "moon s(jd={jd}) = {:?}, expected {:?}",
                out.position_gcrs,
                expected
            );
        }
    }

    // ---- B. 第2項 −v_E·τ の有無・符号（★最重要回帰）----

    #[test]
    fn sun_correction_includes_earth_displacement_term() {
        for &jd in &[J2000_JD, 2469807.0] {
            let out = sun_light_time_corrected_gcrs(tt(jd));
            let naive = sun_geocentric_gcrs(retarded(tt(jd), out.light_time_seconds));
            let theta = angle_between(out.position_gcrs, naive);
            assert!(
                theta > arcsec_to_rad(15.0) && theta < arcsec_to_rad(26.0),
                "sun displacement angle(jd={jd}) = {} arcsec, want ~20.5",
                theta / arcsec_to_rad(1.0)
            );
        }
    }

    // 月は太陽と違い v_E と月視線の角度が任意（v_E ⊥ 視線が成り立たない）ため、第2項の視線直交
    // 成分＝角度ずれは幾何依存で [0, 20.5″] を取りうる（J2000 では 11.4″）。よって月では角度でなく
    // **変位の大きさ |out−naive| = |v_E·τ|** を不変量にする（幾何非依存・第2項欠落で 0km）。
    #[test]
    fn moon_correction_includes_earth_displacement_term() {
        for &jd in &[J2000_JD, 2469807.0] {
            let out = moon_light_time_corrected_gcrs(tt(jd));
            let naive = moon_geocentric_gcrs(retarded(tt(jd), out.light_time_seconds));
            let disp = (out.position_gcrs - naive).norm();
            // |v_E·τ|: v_E≈29.3–30.3 km/s × τ_moon≈1.19–1.36 s ≈ 35–41 km。
            assert!(
                (30.0..45.0).contains(&disp),
                "moon displacement(jd={jd}) = {disp} km, want ~|v_E·tau|≈38"
            );
        }
    }

    #[test]
    fn sun_earth_displacement_term_has_correct_sign() {
        for &jd in &[J2000_JD, 2469807.0] {
            let out = sun_light_time_corrected_gcrs(tt(jd));
            let naive = sun_geocentric_gcrs(retarded(tt(jd), out.light_time_seconds));
            let got_term = out.position_gcrs - naive;
            let expected_term = earth_displacement_term_gcrs(tt(jd), out.light_time_seconds);
            let cos = got_term.dot(expected_term) / (got_term.norm() * expected_term.norm());
            assert!(
                cos > 0.9,
                "sun displacement direction(jd={jd}) cos = {cos}, want >0.9"
            );
        }
    }

    // ---- C. 光行時間 τ ----

    #[test]
    fn sun_light_time_about_499_seconds() {
        for &jd in &[J2000_JD, 2469807.0] {
            let lt = sun_light_time_corrected_gcrs(tt(jd)).light_time_seconds;
            // 近日点 r≈1.471e8 km (τ≈490.7s) 〜遠日点 r≈1.521e8 km (τ≈507.4s) を覆う。
            assert!((488.0..510.0).contains(&lt), "sun tau(jd={jd}) = {lt} s");
        }
    }

    #[test]
    fn moon_light_time_about_one_second() {
        for &jd in &[J2000_JD, 2469807.0] {
            let lt = moon_light_time_corrected_gcrs(tt(jd)).light_time_seconds;
            assert!((1.1..1.4).contains(&lt), "moon tau(jd={jd}) = {lt} s");
        }
    }

    #[test]
    fn sun_and_moon_light_time_differ_by_orders_of_magnitude() {
        for &jd in &[J2000_JD, 2469807.0] {
            let ts = sun_light_time_corrected_gcrs(tt(jd)).light_time_seconds;
            let tm = moon_light_time_corrected_gcrs(tt(jd)).light_time_seconds;
            let ratio = ts / tm;
            assert!(
                (300.0..450.0).contains(&ratio),
                "tau ratio(jd={jd}) = {ratio}"
            );
        }
    }

    #[test]
    fn light_time_consistent_with_output_norm() {
        for &jd in &[J2000_JD, 2469807.0] {
            let sun = sun_light_time_corrected_gcrs(tt(jd));
            assert!(
                close(
                    sun.light_time_seconds,
                    sun.position_gcrs.norm() / SPEED_OF_LIGHT_KM_S,
                    1e-4
                ),
                "sun lt/norm mismatch(jd={jd})"
            );
            let moon = moon_light_time_corrected_gcrs(tt(jd));
            assert!(
                close(
                    moon.light_time_seconds,
                    moon.position_gcrs.norm() / SPEED_OF_LIGHT_KM_S,
                    1e-4
                ),
                "moon lt/norm mismatch(jd={jd})"
            );
        }
    }

    // ---- D. 反復収束（不動点残差）----

    #[test]
    fn sun_iteration_converges_within_tolerance() {
        for &jd in &[J2000_JD, 2469807.0] {
            let out = sun_light_time_corrected_gcrs(tt(jd));
            let tau = out.light_time_seconds;
            let b = sun_geocentric_gcrs(retarded(tt(jd), tau));
            let s = b + earth_displacement_term_gcrs(tt(jd), tau);
            let rhs = s.norm() / SPEED_OF_LIGHT_KM_S;
            assert!(
                (tau - rhs).abs() < 1e-6,
                "sun fixed-point residual(jd={jd}) = {}",
                (tau - rhs).abs()
            );
        }
    }

    #[test]
    fn moon_iteration_converges_within_tolerance() {
        for &jd in &[J2000_JD, 2469807.0] {
            let out = moon_light_time_corrected_gcrs(tt(jd));
            let tau = out.light_time_seconds;
            let b = moon_geocentric_gcrs(retarded(tt(jd), tau));
            let s = b + earth_displacement_term_gcrs(tt(jd), tau);
            let rhs = s.norm() / SPEED_OF_LIGHT_KM_S;
            assert!(
                (tau - rhs).abs() < 1e-6,
                "moon fixed-point residual(jd={jd}) = {}",
                (tau - rhs).abs()
            );
        }
    }

    #[test]
    fn light_time_results_are_finite() {
        for &jd in &[J2000_JD, 2469807.0] {
            for out in [
                sun_light_time_corrected_gcrs(tt(jd)),
                moon_light_time_corrected_gcrs(tt(jd)),
            ] {
                assert!(out.light_time_seconds.is_finite());
                assert!(
                    out.position_gcrs.x.is_finite()
                        && out.position_gcrs.y.is_finite()
                        && out.position_gcrs.z.is_finite()
                );
            }
        }
    }

    // ---- E. 補正後距離オーダー ----

    #[test]
    fn sun_corrected_distance_order_of_magnitude() {
        for &jd in &[J2000_JD, 2469807.0] {
            let r = sun_light_time_corrected_gcrs(tt(jd)).position_gcrs.norm();
            assert!(
                (1.4e8..1.6e8).contains(&r),
                "sun corrected distance(jd={jd}) = {r}"
            );
        }
    }

    #[test]
    fn moon_corrected_distance_order_of_magnitude() {
        for &jd in &[J2000_JD, 2469807.0] {
            let r = moon_light_time_corrected_gcrs(tt(jd)).position_gcrs.norm();
            assert!(
                (356_000.0..407_000.0).contains(&r),
                "moon corrected distance(jd={jd}) = {r}"
            );
        }
    }
}
