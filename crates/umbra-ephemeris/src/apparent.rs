//! 見かけ地心位置（ISSUE-015）。
//!
//! S1（本コミット）= 幾何地心位置のフレーム調和: 黄道座標の暦（太陽 VSOP87D = 黄道 of date、
//! 月 ELP2000-82B = 黄道 J2000）を共通の **GCRS**（ICRS 軸・地心）へ載せる。光行時間 → 光行差 →
//! 歳差章動（GCRS→CIRS）の補正は後続スライス（SOFA `iauAtciq` 順, ISSUE-015 S2–S4）。
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
}
