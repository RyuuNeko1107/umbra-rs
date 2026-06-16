//! 見かけ地心位置（ISSUE-015）。
//!
//! S1 = 幾何地心位置のフレーム調和: 黄道座標の暦（太陽 VSOP87D = 黄道 of date、
//! 月 ELP2000-82B = 黄道 J2000）を共通の **GCRS**（ICRS 軸・地心）へ載せる（`*_geocentric_gcrs`）。
//! S2 = 光行時間補正（`*_light_time_corrected_gcrs`, SOFA `iauAtciq` の light-time ステップ）:
//! 天体 = 放射時刻 t−τ・観測者 = 観測時刻 t を一貫させた幾何地心ベクトル。
//! S3 = 恒星光行差（`*_aberrated_gcrs`, SOFA `iauAb` 逐語, GCRS）。
//! S4（本コミット）= 歳差章動 GCRS→CIRS（`*_apparent_cirs`, `cio::gcrs_to_cirs_matrix` 連鎖）。
//! `*_apparent_cirs` が**フル見かけ位置**（光行時間＋光行差＋歳差章動, 偏向 iauLd 既定 OFF）。
//! 補正順序は SOFA `iauAtciq` / D3（光行差まで GCRS、章動を最後にまとめて回す）。
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
use umbra_core::{JulianDate2, TdbInstant, TtInstant, UnitVector3, Vector3};

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

/// 恒星光行差（年周光行差）を astrometric 単位方向へ適用する（SOFA `iauAb` 逐語）。
///
/// `pnat` = 補正前単位方向、`v` = 観測者速度を c で無次元化（GCRS）、`s_au` = 太陽-観測者距離 \[AU\]、
/// `bm1` = √(1−|v|²)。`w2 = SRS/s_au` は `iauAb` 内蔵の微小項で、角度依存の太陽光偏向 `iauLd`
/// （既定 OFF）とは別物。戻り値は光行差後の単位方向。
fn apply_iau_ab(pnat: UnitVector3, v: Vector3, s_au: f64, bm1: f64) -> UnitVector3 {
    let p = pnat.get();
    let pdv = p.dot(v);
    let w1 = 1.0 + pdv / (1.0 + bm1);
    let w2 = umbra_core::constants::SRS / s_au;
    let aberrated = Vector3 {
        x: p.x * bm1 + w1 * v.x + w2 * (v.x - pdv * p.x),
        y: p.y * bm1 + w1 * v.y + w2 * (v.y - pdv * p.y),
        z: p.z * bm1 + w1 * v.z + w2 * (v.z - pdv * p.z),
    };
    aberrated
        .normalized()
        .expect("aberrated vector is non-zero (|p·bm1 + ...| ≈ 1)")
}

/// astrometric ベクトル（S2 出力, GCRS km）に恒星光行差を適用した見かけ地心位置（GCRS, km）。
/// 観測者速度 = 地球日心速度（黄道 of date → 観測日行列で GCRS）/c。太陽-観測者距離は地球日心
/// 距離 R \[AU\]（`iauAb` の `s`）。光行差は方向のみ変えるため距離 |s2| を保つ。
fn aberrated_gcrs(time_tt: TtInstant, astrometric: Vector3) -> Vector3 {
    let dist = astrometric.norm();
    let pnat = astrometric
        .normalized()
        .expect("astrometric vector is non-zero");
    let tdb = TdbInstant::from_jd2(time_tt.jd2());
    let v = ecliptic_to_gcrs_matrix(time_tt)
        .mul_vec(crate::sun::earth_heliocentric_velocity_ecliptic_of_date(
            tdb,
        ))
        .scale(1.0 / umbra_core::constants::SPEED_OF_LIGHT_KM_S);
    let s_au = crate::sun::earth_heliocentric_lbr(tdb).2;
    let bm1 = (1.0 - v.dot(v)).sqrt();
    apply_iau_ab(pnat, v, s_au, bm1).get().scale(dist)
}

/// 太陽の見かけ地心位置（光行時間＋恒星光行差, GCRS, km）。歳差章動（S4）前・偏向は既定 OFF。
pub fn sun_aberrated_gcrs(time_tt: TtInstant) -> Vector3 {
    aberrated_gcrs(
        time_tt,
        sun_light_time_corrected_gcrs(time_tt).position_gcrs,
    )
}

/// 月の見かけ地心位置（光行時間＋恒星光行差, GCRS, km）。歳差章動（S4）前・偏向は既定 OFF。
pub fn moon_aberrated_gcrs(time_tt: TtInstant) -> Vector3 {
    aberrated_gcrs(
        time_tt,
        moon_light_time_corrected_gcrs(time_tt).position_gcrs,
    )
}

/// 太陽の見かけ地心位置（CIRS, km）。光行時間＋恒星光行差＋歳差章動を適用したフル見かけ位置
/// （SOFA `iauAtciq` 相当, 偏向 iauLd は既定 OFF）。S3 の見かけ GCRS 位置に GCRS→CIRS 回転
/// （`cio::gcrs_to_cirs_matrix` = frame bias + IAU2006 歳差 + IAU2000A 章動, CIO ベース）を当てる。
/// 補正順序 D3: 光行差まで GCRS で適用済み、章動は最後にまとめて回す。回転ゆえ距離保存。
pub fn sun_apparent_cirs(time_tt: TtInstant) -> Vector3 {
    crate::cio::gcrs_to_cirs_matrix(time_tt).mul_vec(sun_aberrated_gcrs(time_tt))
}

/// 月の見かけ地心位置（CIRS, km）。光行時間＋恒星光行差＋歳差章動を適用したフル見かけ位置
/// （SOFA `iauAtciq` 相当, 偏向 iauLd は既定 OFF）。詳細は [`sun_apparent_cirs`] と同型。
pub fn moon_apparent_cirs(time_tt: TtInstant) -> Vector3 {
    crate::cio::gcrs_to_cirs_matrix(time_tt).mul_vec(moon_aberrated_gcrs(time_tt))
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

    // ============================================================
    // S3: 恒星光行差 (sun/moon_aberrated_gcrs) — SOFA iauAb 逐語
    // ============================================================
    //
    // 一次オラクル: pyerfa 2.0.1.5 `erfa.ab(pnat, v, s, bm1)`（liberfa = SOFA 由来 C,
    //   独立実装）。本実装は ab.c 逐語のため単位方向で tol 1e-12 厳密一致するはず。
    //   provenance（入力 pnat,v,s,bm1 と erfa version, 全要素）は EXPECTED 定数の直上に転記する。
    //   erfa.ab に渡す s は本実装が w2=SRS/s に使う値（earth_heliocentric_lbr(tdb).2）と同一。

    // pyerfa 2.0.1.5 `erfa.ab(pnat, v, s, bm1)` 出力（liberfa = SOFA ab.c, 独立実装）。
    // 入力（J2000, 本実装が生成した f64 を逐語転記）:
    //   v   = [-9.93866806674108765e-5, -1.67390893828670883e-5, -7.25658131851178186e-6]（共通）
    //   s   = 9.83327681910549090e-1（= earth_heliocentric_lbr(tdb).2, AU）
    //   bm1 = 9.99999994894716249e-1
    //   pnat(SUN)  = [ 1.80138755190583838e-1, -9.02474816600209229e-1, -3.91266193633957426e-1]
    //   pnat(MOON) = [-7.24548996157033720e-1, -6.62772556607322594e-1, -1.89106558257581575e-1]
    /// erfa.ab の単位方向期待値（太陽, J2000）。pyerfa 出力を逐語転記（桁保持のため allow）。
    #[allow(clippy::excessive_precision)]
    const SUN_AB_DIR_J2000: [f64; 3] = [
        1.800_393_599_401_254_32e-1,
        -9.024_915_127_553_109_21e-1,
        -3.912_734_316_012_022_04e-1,
    ];
    /// erfa.ab の単位方向期待値（月, J2000）。pyerfa 出力を逐語転記（桁保持のため allow）。
    #[allow(clippy::excessive_precision)]
    const MOON_AB_DIR_J2000: [f64; 3] = [
        -7.245_871_723_897_979_73e-1,
        -6.627_333_073_004_310_07e-1,
        -1.890_978_397_623_566_68e-1,
    ];

    /// 単位方向を要素ごとに期待値と比較。
    fn unit_close(v: Vector3, expected: [f64; 3], tol: f64) -> bool {
        let u = v.normalized().expect("non-zero").get();
        close(u.x, expected[0], tol) && close(u.y, expected[1], tol) && close(u.z, expected[2], tol)
    }

    // ---- (a) erfa.ab 厳密一致（主オラクル, tol 1e-12）----

    #[test]
    fn sun_aberrated_matches_erfa_ab() {
        let out = sun_aberrated_gcrs(tt(J2000_JD));
        assert!(
            unit_close(out, SUN_AB_DIR_J2000, 1e-12),
            "sun ab dir = {:?}, expected {SUN_AB_DIR_J2000:?}",
            out.normalized().map(|u| u.get())
        );
    }

    #[test]
    fn moon_aberrated_matches_erfa_ab() {
        let out = moon_aberrated_gcrs(tt(J2000_JD));
        assert!(
            unit_close(out, MOON_AB_DIR_J2000, 1e-12),
            "moon ab dir = {:?}, expected {MOON_AB_DIR_J2000:?}",
            out.normalized().map(|u| u.get())
        );
    }

    // `apply_iau_ab` を**増幅速度**（|v|≈0.06, s=0.5）で erfa.ab と直接突合。実エポックでは
    // w2=SRS/s 微項（相対論補正 ~0.004″）の単位ベクトル寄与が ~1e-12 で tol に埋もれ、w2 項
    // 内部（SRS/s 除算・(v−pdv·pnat)・符号）のミューテーションが生存する。純関数を非物理だが
    // 有効な大入力で直接検証して w2 項を ~1e-9 に励起し捕捉する（s=0.5 で SRS/s と SRS*s/SRS%s を区別）。
    /// erfa.ab(pnat,v,0.5,bm1) 期待値（増幅入力）。pyerfa 2.0.1.5 出力。
    /// pnat=unit([0.3,-0.8,0.5]), v=[0.05,-0.02,0.03], s=0.5, bm1=√(1−|v|²)。
    #[allow(clippy::excessive_precision)]
    const AB_AMPLIFIED_PPR: [f64; 3] = [
        3.379_296_290_574_336_23e-1,
        -7.903_261_551_162_078_51e-1,
        5.110_656_849_606_880_49e-1,
    ];

    #[test]
    fn apply_iau_ab_matches_erfa_ab_amplified() {
        let pnat = (Vector3 {
            x: 0.3,
            y: -0.8,
            z: 0.5,
        })
        .normalized()
        .expect("non-zero");
        let v = Vector3 {
            x: 0.05,
            y: -0.02,
            z: 0.03,
        };
        let s_au = 0.5;
        let bm1 = (1.0 - v.dot(v)).sqrt();
        let got = apply_iau_ab(pnat, v, s_au, bm1).get();
        assert!(
            unit_close(got, AB_AMPLIFIED_PPR, 1e-12),
            "apply_iau_ab amplified = {got:?}, expected {AB_AMPLIFIED_PPR:?}"
        );
    }

    // ---- (b) 距離不変: |out| == |s_S2| ----

    #[test]
    fn sun_aberration_preserves_distance() {
        for &jd in &[J2000_JD, 2469807.0] {
            let s2 = sun_light_time_corrected_gcrs(tt(jd)).position_gcrs.norm();
            let out = sun_aberrated_gcrs(tt(jd)).norm();
            assert!(
                close(out, s2, s2 * 1e-9),
                "sun |out|={out} |s2|={s2} (jd={jd})"
            );
        }
    }

    #[test]
    fn moon_aberration_preserves_distance() {
        for &jd in &[J2000_JD, 2469807.0] {
            let s2 = moon_light_time_corrected_gcrs(tt(jd)).position_gcrs.norm();
            let out = moon_aberrated_gcrs(tt(jd)).norm();
            assert!(
                close(out, s2, s2 * 1e-9),
                "moon |out|={out} |s2|={s2} (jd={jd})"
            );
        }
    }

    // ---- (c) 光行差角: 太陽 ≈ 20.5″（apex 満角, v⊥pnat）----

    #[test]
    fn sun_aberration_angle_about_20_arcsec() {
        for &jd in &[J2000_JD, 2469807.0] {
            let s2 = sun_light_time_corrected_gcrs(tt(jd)).position_gcrs;
            let out = sun_aberrated_gcrs(tt(jd));
            let theta = angle_between(out, s2);
            assert!(
                theta > arcsec_to_rad(15.0) && theta < arcsec_to_rad(26.0),
                "sun aberration angle(jd={jd}) = {} arcsec, want ~20.5",
                theta / arcsec_to_rad(1.0)
            );
        }
    }

    // 月は pnat と v_E の角が任意 → aberration 角 ∈ [0,20.5]″（幾何依存）。下限は pin せず
    // 物理上限のみサニティ。月の向き正しさは erfa.ab 厳密一致が担保。
    #[test]
    fn moon_aberration_angle_within_physical_bound() {
        for &jd in &[J2000_JD, 2469807.0] {
            let s2 = moon_light_time_corrected_gcrs(tt(jd)).position_gcrs;
            let out = moon_aberrated_gcrs(tt(jd));
            let theta = angle_between(out, s2);
            assert!(
                theta <= arcsec_to_rad(21.0),
                "moon aberration angle(jd={jd}) = {} arcsec exceeds 20.5 bound",
                theta / arcsec_to_rad(1.0)
            );
        }
    }

    // ---- (d)/(e) シフト向き = apex(+v) 横成分（太陽。bm1/w1/w2 符号は erfa.ab が捕捉）----

    #[test]
    fn aberration_shift_is_along_apex_direction() {
        for &jd in &[J2000_JD, 2469807.0] {
            let s2 = sun_light_time_corrected_gcrs(tt(jd)).position_gcrs;
            let out = sun_aberrated_gcrs(tt(jd));
            let tdb_now = TdbInstant::from_jd2(tt(jd).jd2());
            let v = ecliptic_to_gcrs_matrix(tt(jd))
                .mul_vec(crate::sun::earth_heliocentric_velocity_ecliptic_of_date(
                    tdb_now,
                ))
                .scale(1.0 / SPEED_OF_LIGHT_KM_S);
            let pnat = s2.scale(1.0 / s2.norm());
            let v_perp = v - pnat.scale(v.dot(pnat));
            let shift = out.scale(s2.norm() / out.norm()) - s2;
            let cos = shift.dot(v_perp) / (shift.norm() * v_perp.norm());
            assert!(
                cos > 0.9,
                "sun shift vs apex cos(jd={jd}) = {cos}, want >0.9"
            );
        }
    }

    // ---- (f) finite ----

    #[test]
    fn aberrated_results_are_finite() {
        for &jd in &[J2000_JD, 2469807.0] {
            for v in [sun_aberrated_gcrs(tt(jd)), moon_aberrated_gcrs(tt(jd))] {
                assert!(v.x.is_finite() && v.y.is_finite() && v.z.is_finite());
            }
        }
    }

    // ---- オーダーサニティ（距離保持の独立確認）----

    #[test]
    fn sun_aberrated_distance_order_of_magnitude() {
        for &jd in &[J2000_JD, 2469807.0] {
            let r = sun_aberrated_gcrs(tt(jd)).norm();
            assert!(
                (1.4e8..1.6e8).contains(&r),
                "sun ab distance(jd={jd}) = {r}"
            );
        }
    }

    #[test]
    fn moon_aberrated_distance_order_of_magnitude() {
        for &jd in &[J2000_JD, 2469807.0] {
            let r = moon_aberrated_gcrs(tt(jd)).norm();
            assert!(
                (356_000.0..407_000.0).contains(&r),
                "moon ab distance(jd={jd}) = {r}"
            );
        }
    }

    // ============================================================
    // S4: 歳差章動 GCRS→CIRS (sun/moon_apparent_cirs) — SOFA iauC2i06a 合成
    // ============================================================
    //
    // 契約: apparent_cirs(t) = gcrs_to_cirs_matrix(t) · aberrated_gcrs(t)（回転のみ→距離保存）。
    // end-to-end オラクル: CIRS 見かけ単位方向 = c2i06a(2451545.0,0.0) × S3 aberrated 単位方向。
    //   c2i06a は cio.rs で erfa 検証済み、aberrated 単位方向は S3 で erfa.ab 検証済み。
    //   slice 境界: JPL DE 突合(M10)・erfa.atci13(別 API)・偏向 iauLd(既定 OFF) は対象外。

    use crate::cio::gcrs_to_cirs_matrix;

    // pyerfa 2.0.1.5: C2I = erfa.c2i06a(2451545.0, 0.0)[行優先]; v = C2I @ ab; v/|v|。
    //   ab = SUN_AB_DIR_J2000 / MOON_AB_DIR_J2000（S3 検証済み aberrated 単位方向）。
    /// CIRS 見かけ単位方向の期待値（太陽, J2000）。
    #[allow(clippy::excessive_precision)]
    const SUN_APPARENT_CIRS_DIR_J2000: [f64; 3] = [
        1.800_288_076_667_243_83e-1,
        -9.025_024_717_971_028_56e-1,
        -3.912_530_086_915_844_70e-1,
    ];
    /// CIRS 見かけ単位方向の期待値（月, J2000）。
    #[allow(clippy::excessive_precision)]
    const MOON_APPARENT_CIRS_DIR_J2000: [f64; 3] = [
        -7.245_922_740_951_216_13e-1,
        -6.627_385_950_564_382_84e-1,
        -1.890_597_549_566_767_62e-1,
    ];

    // end-to-end の tol は cio.rs の c2i06a 行列許容 1e-9（実装章動 R06 直接評価 vs erfa nut06a
    // スケーリング近似の representation 差 ~2e-11 を吸収）に整合。実装 gcrs_to_cirs_matrix と
    // erfa c2i06a の差が単位方向に乗るため 1e-12 でなく 1e-9（≈2e-4″, 0.10″ 予算の 3 桁下）。
    const CIRS_DIR_TOL: f64 = 1e-9;

    // ---- (a) end-to-end erfa オラクル一致（主オラクル, J2000）----

    #[test]
    fn sun_apparent_cirs_matches_erfa_c2i_ab() {
        let out = sun_apparent_cirs(tt(J2000_JD));
        assert!(
            unit_close(out, SUN_APPARENT_CIRS_DIR_J2000, CIRS_DIR_TOL),
            "sun apparent CIRS dir = {:?}, expected {SUN_APPARENT_CIRS_DIR_J2000:?}",
            out.normalized().map(|u| u.get())
        );
    }

    #[test]
    fn moon_apparent_cirs_matches_erfa_c2i_ab() {
        let out = moon_apparent_cirs(tt(J2000_JD));
        assert!(
            unit_close(out, MOON_APPARENT_CIRS_DIR_J2000, CIRS_DIR_TOL),
            "moon apparent CIRS dir = {:?}, expected {MOON_APPARENT_CIRS_DIR_J2000:?}",
            out.normalized().map(|u| u.get())
        );
    }

    // ---- (b) 合成同一性: out == gcrs_to_cirs_matrix(t) · aberrated_gcrs(t)（2エポック）----

    #[test]
    fn sun_apparent_cirs_equals_matrix_times_aberrated() {
        for &jd in &[J2000_JD, 2469807.0] {
            let expected = gcrs_to_cirs_matrix(tt(jd)).mul_vec(sun_aberrated_gcrs(tt(jd)));
            let got = sun_apparent_cirs(tt(jd));
            assert!(
                vec_close(got, expected, 1e-3),
                "sun apparent CIRS(jd={jd}) = {got:?}, expected {expected:?}"
            );
        }
    }

    #[test]
    fn moon_apparent_cirs_equals_matrix_times_aberrated() {
        for &jd in &[J2000_JD, 2469807.0] {
            let expected = gcrs_to_cirs_matrix(tt(jd)).mul_vec(moon_aberrated_gcrs(tt(jd)));
            let got = moon_apparent_cirs(tt(jd));
            assert!(
                vec_close(got, expected, 1e-6),
                "moon apparent CIRS(jd={jd}) = {got:?}, expected {expected:?}"
            );
        }
    }

    // ---- (c) 距離保存: |out| == |aberrated_gcrs|（回転のみ, 2エポック）----

    #[test]
    fn sun_apparent_cirs_preserves_distance() {
        for &jd in &[J2000_JD, 2469807.0] {
            let s = sun_aberrated_gcrs(tt(jd)).norm();
            let out = sun_apparent_cirs(tt(jd)).norm();
            assert!(
                close(out, s, s * 1e-9),
                "sun |apparent|={out} |aberrated|={s} (jd={jd})"
            );
        }
    }

    #[test]
    fn moon_apparent_cirs_preserves_distance() {
        for &jd in &[J2000_JD, 2469807.0] {
            let s = moon_aberrated_gcrs(tt(jd)).norm();
            let out = moon_apparent_cirs(tt(jd)).norm();
            assert!(
                close(out, s, s * 1e-9),
                "moon |apparent|={out} |aberrated|={s} (jd={jd})"
            );
        }
    }

    // ---- (d) finite ----

    #[test]
    fn apparent_cirs_results_are_finite() {
        for &jd in &[J2000_JD, 2469807.0] {
            for v in [sun_apparent_cirs(tt(jd)), moon_apparent_cirs(tt(jd))] {
                assert!(v.x.is_finite() && v.y.is_finite() && v.z.is_finite());
            }
        }
    }

    // ---- (e) オーダーサニティ ----

    #[test]
    fn sun_apparent_cirs_distance_order_of_magnitude() {
        for &jd in &[J2000_JD, 2469807.0] {
            let r = sun_apparent_cirs(tt(jd)).norm();
            assert!(
                (1.4e8..1.6e8).contains(&r),
                "sun apparent CIRS distance(jd={jd}) = {r}"
            );
        }
    }

    #[test]
    fn moon_apparent_cirs_distance_order_of_magnitude() {
        for &jd in &[J2000_JD, 2469807.0] {
            let r = moon_apparent_cirs(tt(jd)).norm();
            assert!(
                (356_000.0..407_000.0).contains(&r),
                "moon apparent CIRS distance(jd={jd}) = {r}"
            );
        }
    }

    // ---- (f) 回転を実際に適用している（恒等＝aberrated そのまま返す不完全実装を red に）----
    // J2000 では C2I の回転角が frame bias 支配で ≈8″ と小さく恒等(0″)と紛れる。歳差が累積する
    // 2469807.0(~2050) では実回転が数百″（太陽 433″/月 995″, CIP 軸との位置で幾何依存）になるため、
    // 恒等(0″)・bias 級(≈8″)と実回転を安全に弁別する閾値 100″ で「回転が適用されていること」を要求する。

    #[test]
    fn sun_apparent_cirs_rotation_is_applied_not_identity() {
        let jd = 2469807.0;
        let theta = angle_between(sun_apparent_cirs(tt(jd)), sun_aberrated_gcrs(tt(jd)));
        assert!(
            theta > arcsec_to_rad(100.0),
            "sun apparent CIRS(jd={jd}) rotated only {} arcsec from GCRS (no-rotation impl?)",
            theta / arcsec_to_rad(1.0)
        );
    }

    #[test]
    fn moon_apparent_cirs_rotation_is_applied_not_identity() {
        let jd = 2469807.0;
        let theta = angle_between(moon_apparent_cirs(tt(jd)), moon_aberrated_gcrs(tt(jd)));
        assert!(
            theta > arcsec_to_rad(100.0),
            "moon apparent CIRS(jd={jd}) rotated only {} arcsec from GCRS (no-rotation impl?)",
            theta / arcsec_to_rad(1.0)
        );
    }
}
