//! 地心合の精解（ISSUE-017）。
//!
//! ISSUE-016 の各新月候補窓内で、太陽と月の**地心合**（黄経差 Δλ または赤経差 Δα のゼロ点）を
//! Brent 法で精密に解く層。`solve_conjunction` が候補窓を受け取り、合の種別 [`ConjunctionKind`]
//! （黄経合 / 赤経合）に応じた角度差の符号変化を窓内で粗走査でブラケットし、`solve_zero_in_window`
//! 経由で Brent 求根する。返り値 [`Conjunction`] は合時刻 `time_tt`・種別・合時刻の月-太陽角距離
//! `separation`（早期棄却 ISSUE-018 の入力）を持つ。
//!
//! 角度差 Δ(t) は見かけ位置（ISSUE-015）から構成する: 黄経合は GCRS 見かけ位置を平均黄道 of date へ
//! 回した λ、赤経合は CIRS 見かけ位置の α。`Δ = (angle_moon − angle_sun)` を `[-π,π)` 連続化し（±1日窓で
//! 折返しなし・月が太陽を追い越すため単調増加）、窓内を粗走査で符号変化ブラケット→ Brent 求根する
//! （Newton 単独禁止, conventions §11）。`separation` は合時刻の月-太陽角距離（acos クランプ）。
//!
//! 注: 仕様の `eph: &impl Ephemeris` 引数は省略した。現アーキの `apparent::*`（見かけ位置）が VSOP/ELP
//! 直結で Ephemeris ジェネリックでないため（ISSUE-037 と同じ繰延）。暦の差し替えは ISSUE-043 で統合する。

// 粗走査の分割数・収束反復の整数→f64 変換は小さな添字のみ（精度クリティカルな天文量ではない）。
#![allow(clippy::cast_precision_loss)]
// solve_conjunction / solve_zero_in_window は ISSUE-018（候補棄却）/023 が消費する。結線され次第
// この許容は外す（candidates.rs と同手順）。
#![allow(dead_code)]

use umbra_core::{brent_root, JulianDate2, Radians, SolverError, TtInstant, Vector3};
use umbra_ephemeris::apparent::{
    moon_aberrated_gcrs, moon_apparent_cirs, sun_aberrated_gcrs, sun_apparent_cirs,
};
use umbra_ephemeris::frames::ecliptic_to_gcrs_matrix;

use crate::candidates::NewMoonCandidate;
use crate::error::EclipseError;

/// 粗走査の分割数。窓（~±1 日）で Δ 角は単調 1 通過ゆえ十分（偽陰性防止に余裕を持たせる）。
const COARSE_SCAN_INTERVALS: usize = 16;

/// 合の種別（どの座標で合を定義するか）。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum ConjunctionKind {
    /// `Δλ = λ_moon − λ_sun = 0`（黄経合, 既定。見かけ位置を平均黄道 of date へ射影）。
    EclipticLongitude,
    /// `Δα = α_moon − α_sun = 0`（赤経合, CIRS 見かけ位置の赤経）。
    RightAscension,
}

/// 地心合の精解結果。
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct Conjunction {
    /// 合の TT 時刻。
    pub time_tt: TtInstant,
    /// 合の種別。
    pub kind: ConjunctionKind,
    /// 合時刻の月-太陽角距離（早期棄却 ISSUE-018 の入力）。
    pub separation: Radians,
}

/// 求根設定（ISSUE-008 Brent へ渡す）。
#[derive(Debug, Clone, Copy)]
pub(crate) struct RootConfig {
    /// Brent の独立変数（日）の収束許容。**実用下限 ≈1e-9 日**: 求根は単一 f64 の TT-JD（≈2.46e6,
    /// ULP≈4.6e-10 日）で動くため、これより小さい許容は区間幅条件が成立せず収束は `f==0` か
    /// 反復上限（`DidNotConverge`）に依存する。±1〜2s 目標には 1e-7〜1e-9 日で十分。
    pub x_tolerance_days: f64,
    /// 反復上限。
    pub max_iterations: usize,
}

/// 単一 f64 の TT-JD から TtInstant（求根は単一 f64 JD 空間で動く。2 要素分割は collapse され、
/// ±1 日窓・JD≈2.46e6 での再構成誤差 ≈ULP≈4.6e-10 日 ≈0.04ms は ±1〜2s 目標の内側で許容）。
fn tt(jd: f64) -> TtInstant {
    TtInstant::from_jd2(JulianDate2::from_jd(jd))
}

/// 赤経 α = atan2(y, x)（CIRS 見かけ位置, rad）。
fn right_ascension(cirs: Vector3) -> f64 {
    cirs.y.atan2(cirs.x)
}

/// 黄経 λ（平均黄道 of date, rad）。GCRS 見かけ位置を `ecliptic_to_gcrs_matrix(t).transpose()`
/// （= GCRS→黄道 of date）で回し λ = atan2(y, x)。共通の章動・分点オフセットは Δλ で相殺する。
fn ecliptic_longitude(aberrated_gcrs: Vector3, t: TtInstant) -> f64 {
    let ecl = ecliptic_to_gcrs_matrix(t)
        .transpose()
        .mul_vec(aberrated_gcrs);
    ecl.y.atan2(ecl.x)
}

/// 種別に応じた連続化角度差 `Δ(t) = (angle_moon − angle_sun)`（[-π, π)）。合時刻でゼロになる。
fn angle_difference(kind: ConjunctionKind, t: TtInstant) -> f64 {
    let (a_moon, a_sun) = match kind {
        ConjunctionKind::EclipticLongitude => (
            ecliptic_longitude(moon_aberrated_gcrs(t), t),
            ecliptic_longitude(sun_aberrated_gcrs(t), t),
        ),
        ConjunctionKind::RightAscension => (
            right_ascension(moon_apparent_cirs(t)),
            right_ascension(sun_apparent_cirs(t)),
        ),
    };
    Radians::new(a_moon - a_sun).normalized_signed().0
}

/// 太陽-月の見かけ地心方向（CIRS）のなす角（離角, rad, acos クランプ）。
fn elongation(t: TtInstant) -> f64 {
    let s = sun_apparent_cirs(t);
    let m = moon_apparent_cirs(t);
    (s.dot(m) / (s.norm() * m.norm())).clamp(-1.0, 1.0).acos()
}

/// 窓 `[t0_jd, t1_jd]`（TT-JD）内で連続関数 `f` のゼロ点を粗走査でブラケット→ Brent 求根し、根の
/// TT-JD を返す。符号変化なし → `Err(Solver(RootNotBracketed))`。Newton 単独は使わない。
pub(crate) fn solve_zero_in_window<F: FnMut(f64) -> f64>(
    mut f: F,
    t0_jd: f64,
    t1_jd: f64,
    config: RootConfig,
) -> Result<f64, EclipseError> {
    let n = COARSE_SCAN_INTERVALS;
    let mut prev_jd = t0_jd;
    let mut prev_f = f(prev_jd);
    if prev_f == 0.0 {
        return Ok(prev_jd);
    }
    for i in 1..=n {
        let frac = i as f64 / n as f64;
        let cur_jd = t0_jd + (t1_jd - t0_jd) * frac;
        let cur_f = f(cur_jd);
        if cur_f == 0.0 {
            return Ok(cur_jd);
        }
        if prev_f * cur_f < 0.0 {
            // 符号変化したサブ区間 [prev_jd, cur_jd] を Brent で精解。
            let root = brent_root(
                &mut f,
                prev_jd,
                cur_jd,
                config.x_tolerance_days,
                config.max_iterations,
            )?;
            return Ok(root);
        }
        prev_jd = cur_jd;
        prev_f = cur_f;
    }
    Err(EclipseError::Solver(SolverError::RootNotBracketed))
}

/// 候補窓内で太陽-月の地心合を精解する（ISSUE-016 の窓 → 合の TT 時刻）。
pub(crate) fn solve_conjunction(
    candidate: &NewMoonCandidate,
    kind: ConjunctionKind,
    config: RootConfig,
) -> Result<Conjunction, EclipseError> {
    let t0 = candidate.search_window.start.jd2().jd();
    let t1 = candidate.search_window.end.jd2().jd();
    let root_jd = solve_zero_in_window(|jd| angle_difference(kind, tt(jd)), t0, t1, config)?;
    let time_tt = tt(root_jd);
    let separation = Radians(elongation(time_tt));
    Ok(Conjunction {
        time_tt,
        kind,
        separation,
    })
}

#[cfg(test)]
mod tests {
    // 実装本体（同モジュール直下）が定義する型・関数。impl 担当の `use` 構成に依存しないよう、
    // テスト側で必要なシンボルを `super::*` で取り込む（candidates.rs と同手順）。
    use super::*;

    use crate::candidates::{new_moon_candidates, NewMoonCandidate};
    use crate::error::EclipseError;

    use umbra_core::{
        JulianDate2, Radians, SolverError, TimeRange, TtInstant, UtcInstant, Vector3,
    };
    use umbra_ephemeris::apparent::{
        moon_aberrated_gcrs, moon_apparent_cirs, sun_aberrated_gcrs, sun_apparent_cirs,
    };
    use umbra_ephemeris::frames::ecliptic_to_gcrs_matrix;

    // ============================================================
    // 共通ヘルパ・定数
    // ============================================================

    /// 許容つきスカラ比較（clippy::float_cmp 回避）。
    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    /// 単一 JD（UTC スケール）から UtcInstant。
    fn utc_from_jd(jd: f64) -> UtcInstant {
        UtcInstant::from_jd2(JulianDate2::from_jd(jd))
    }

    /// UTC の `[start_jd, end_jd]` 範囲。
    fn utc_range(start_jd: f64, end_jd: f64) -> TimeRange<UtcInstant> {
        TimeRange {
            start: utc_from_jd(start_jd),
            end: utc_from_jd(end_jd),
        }
    }

    /// 単一 JD（TT スケール）から TtInstant。
    fn tt_from_jd(jd: f64) -> TtInstant {
        TtInstant::from_jd2(JulianDate2::from_jd(jd))
    }

    /// TtInstant を JD（TT）で取り出す。
    fn jd_of(t: TtInstant) -> f64 {
        t.jd2().jd()
    }

    /// 緩めの収束許容を持つ標準 RootConfig（実 ephemeris 用）。
    /// x_tolerance_days = 1e-7 日 ≈ 8.6 ms。角速度 Δλ≈12.2°/日 → 角度残差は ~0.04″ 級。
    fn config_tight() -> RootConfig {
        RootConfig {
            x_tolerance_days: 1e-7,
            max_iterations: 100,
        }
    }

    /// 機構テスト（合成関数）用の RootConfig。
    fn config_for(x_tol: f64) -> RootConfig {
        RootConfig {
            x_tolerance_days: x_tol,
            max_iterations: 200,
        }
    }

    // ---- 角度オラクル（テスト側で独立に計算。実装の内部関数に依存しない）----

    /// CIRS 見かけ位置の赤経 α = atan2(y, x)（rad）。
    fn right_ascension_of(cirs: Vector3) -> f64 {
        cirs.y.atan2(cirs.x)
    }

    /// 平均黄道 of date の黄経 λ（rad）。GCRS 見かけ位置を `ecliptic_to_gcrs_matrix(t).transpose()`
    /// （= GCRS→黄道 of date）で黄道系へ回し λ = atan2(y, x)。
    fn ecliptic_longitude_of(aberrated_gcrs: Vector3, t: TtInstant) -> f64 {
        let ecl = ecliptic_to_gcrs_matrix(t)
            .transpose()
            .mul_vec(aberrated_gcrs);
        ecl.y.atan2(ecl.x)
    }

    /// 黄経合の角度差 Δλ = λ_moon − λ_sun（連続化, [-π, π)）。±1日窓で |Δλ|<13° → 折返しなし。
    fn delta_lambda(t: TtInstant) -> f64 {
        let lam_moon = ecliptic_longitude_of(moon_aberrated_gcrs(t), t);
        let lam_sun = ecliptic_longitude_of(sun_aberrated_gcrs(t), t);
        Radians::new(lam_moon - lam_sun).normalized_signed().0
    }

    /// 赤経合の角度差 Δα = α_moon − α_sun（連続化, [-π, π)）。
    fn delta_alpha(t: TtInstant) -> f64 {
        let a_moon = right_ascension_of(moon_apparent_cirs(t));
        let a_sun = right_ascension_of(sun_apparent_cirs(t));
        Radians::new(a_moon - a_sun).normalized_signed().0
    }

    /// 種別に応じた角度差（テスト用オラクル, rad）。合時刻でゼロになるはずの量。
    fn angle_difference(kind: ConjunctionKind, t: TtInstant) -> f64 {
        match kind {
            ConjunctionKind::EclipticLongitude => delta_lambda(t),
            ConjunctionKind::RightAscension => delta_alpha(t),
        }
    }

    /// 太陽・月の見かけ地心方向（CIRS）のなす角（離角, rad）。実装と同じ acos クランプ。
    fn elongation_rad(t: TtInstant) -> f64 {
        let s = sun_apparent_cirs(t);
        let m = moon_apparent_cirs(t);
        (s.dot(m) / (s.norm() * m.norm())).clamp(-1.0, 1.0).acos()
    }

    /// 区間 `[lo, hi]`（TT-JD）を細かく走査し、離角が最小となる TT-JD を返す（独立オラクル7用）。
    fn elongation_min_jd(lo: f64, hi: f64, step: f64) -> f64 {
        let mut best_jd = lo;
        let mut best = f64::INFINITY;
        let mut jd = lo;
        while jd <= hi {
            let e = elongation_rad(tt_from_jd(jd));
            if e < best {
                best = e;
                best_jd = jd;
            }
            jd += step;
        }
        best_jd
    }

    // ---- 2017-08-21 皆既日食の候補を取り出す ----

    /// 代表 JD（UTC）。2017-08-21 ≈ 2457987.0、皆既帯は 18:30 UTC 付近 ≈ 2457987.27。
    const JD_2017_08_21: f64 = 2_457_987.0;
    /// 2017-08 の朔候補を確実に含む UTC-JD 範囲（朔は 2017-08-21 18:30 UTC 付近）。
    const JD_2017_08_RANGE_LO: f64 = 2_457_980.0;
    const JD_2017_08_RANGE_HI: f64 = 2_457_995.0;

    /// 2017-08-21 の朔（皆既日食）に対応する候補を返す。approx_tt が 2017-08-21 に最も近い候補。
    fn candidate_2017_08_21() -> NewMoonCandidate {
        let candidates = new_moon_candidates(utc_range(JD_2017_08_RANGE_LO, JD_2017_08_RANGE_HI))
            .expect("post-1972 range");
        // approx_tt(TT-JD) が 2017-08-21(UTC-JD≈2457987.0) に最も近い候補を選ぶ。
        // ΔT≈69 s（≈8e-4 日）は朔の選別に影響しない。
        candidates
            .into_iter()
            .min_by(|a, b| {
                (jd_of(a.approx_tt) - JD_2017_08_21)
                    .abs()
                    .total_cmp(&(jd_of(b.approx_tt) - JD_2017_08_21).abs())
            })
            .expect("non-empty candidates near 2017-08")
    }

    /// 候補窓の TT-JD 区間 `[start, end]`。
    fn window_jd(c: &NewMoonCandidate) -> (f64, f64) {
        (
            c.search_window.start.jd2().jd(),
            c.search_window.end.jd2().jd(),
        )
    }

    // ============================================================
    // solve_zero_in_window（機構, 合成関数で検証）
    // ============================================================

    /// 契約1: 既知のゼロ点を解く。f(jd)=jd−target（線形, 根=target）を target 内包の窓で解くと、
    /// 返り値（根の TT-JD）が target に x_tolerance_days まで一致する。求根機構の基本健全性。
    #[test]
    fn solve_zero_in_window_finds_known_linear_root() {
        let target = 2_458_000.3;
        let tol = 1e-9;
        let root = solve_zero_in_window(
            |jd| jd - target,
            target - 1.0,
            target + 1.0,
            config_for(tol),
        )
        .expect("linear root inside window must be found");
        assert!(
            close(root, target, tol.max(1e-9) * 10.0),
            "root {root} != target {target} within tolerance"
        );
    }

    /// 契約2a: 区間で f が同符号（定数 1.0, 根なし）→ Err(Solver(RootNotBracketed))。
    /// 符号変化が無い窓を求根しようとすると粗走査がブラケットに失敗する。
    #[test]
    fn solve_zero_in_window_constant_positive_is_not_bracketed() {
        let result = solve_zero_in_window(|_jd| 1.0, 2_458_000.0, 2_458_001.0, config_for(1e-9));
        assert!(
            matches!(
                result,
                Err(EclipseError::Solver(SolverError::RootNotBracketed))
            ),
            "constant-sign function must yield RootNotBracketed, got {result:?}"
        );
    }

    /// 契約2b: 区間内に根が無い線形（f(jd)=jd−target, target が窓の外）→ 同符号 → RootNotBracketed。
    #[test]
    fn solve_zero_in_window_root_outside_window_is_not_bracketed() {
        // 根 target は窓 [lo, hi] の外（右側）。窓内で f は常に負。
        let lo = 2_458_000.0;
        let hi = 2_458_001.0;
        let target = hi + 5.0;
        let result = solve_zero_in_window(|jd| jd - target, lo, hi, config_for(1e-9));
        assert!(
            matches!(
                result,
                Err(EclipseError::Solver(SolverError::RootNotBracketed))
            ),
            "root outside window must yield RootNotBracketed, got {result:?}"
        );
    }

    /// 契約3: 単調でない単根（f(jd)=(jd−target)³）。導関数が根で 0 でも、粗走査が符号変化を捉え
    /// Brent で解ける。返り値は target に緩く一致（立方は根近傍で平坦のため tol は緩め）。
    #[test]
    fn solve_zero_in_window_cubic_single_root() {
        let target = 2_458_000.6;
        let root = solve_zero_in_window(
            |jd| (jd - target).powi(3),
            target - 1.0,
            target + 1.0,
            config_for(1e-9),
        )
        .expect("cubic single root must be found");
        // 立方の根近傍は平坦なので x の収束は緩い。1e-3 日（≈1.4分）以内で十分。
        assert!(
            close(root, target, 1e-3),
            "cubic root {root} not near target {target}"
        );
    }

    /// 契約4: 小さい x_tolerance_days（1e-7 日 ≈ 8.6ms）でも線形根が tol 内に収束する。
    /// 求根の収束許容が x_tolerance_days に追従していること。
    #[test]
    fn solve_zero_in_window_converges_to_small_tolerance() {
        let target = 2_458_123.456;
        let x_tol = 1e-7;
        let root = solve_zero_in_window(
            |jd| jd - target,
            target - 1.0,
            target + 1.0,
            config_for(x_tol),
        )
        .expect("root must be found");
        // 線形なら根は機械精度級に求まる。x_tol の数倍を上限にする。
        assert!(
            close(root, target, x_tol * 10.0),
            "root {root} not within {x_tol} of target {target}"
        );
    }

    /// 機構の健全性: 根が窓 [t0, t1] の内側に返ること（線形・窓中央でない根でも端を返さない）。
    #[test]
    fn solve_zero_in_window_root_lies_inside_window() {
        let t0 = 2_458_000.0;
        let t1 = 2_458_002.0;
        let target = t0 + 0.7; // 窓中央でない
        let root =
            solve_zero_in_window(|jd| jd - target, t0, t1, config_for(1e-9)).expect("root inside");
        assert!(
            t0 <= root && root <= t1,
            "root {root} not inside window [{t0}, {t1}]"
        );
    }

    /// 機構: 粗走査が窓の**後半**まで覆うこと。根を窓 [t0, t1] の後半（t0+1.7, 幅2日）に置くと、走査が
    /// 窓幅を正しく刻んで初めて符号変化を捉えられる。走査の刻み・終端が前半までしか届かないと根を
    /// 取りこぼし RootNotBracketed になる（偽陰性）。前半のみの既存テストと合わせ全窓を縛る。
    #[test]
    fn solve_zero_in_window_finds_root_in_second_half() {
        let t0 = 2_458_000.0;
        let t1 = 2_458_002.0;
        let target = t0 + 1.7; // 窓の後半
        let root = solve_zero_in_window(|jd| jd - target, t0, t1, config_for(1e-9))
            .expect("root in second half of window must be found");
        assert!(
            close(root, target, 1e-7),
            "second-half root {root} != target {target}"
        );
    }

    // ============================================================
    // solve_conjunction（合, 実 ephemeris で検証）
    // ============================================================

    /// 契約5（黄経合）: 返り値 time_tt で黄経差 |Δλ| がほぼ 0。テスト側の独立オラクル
    /// `delta_lambda` で評価する（実装内部関数は使わない）。tol は緩めに 5e-6 rad（≈1″）。
    #[test]
    fn solve_conjunction_ecliptic_longitude_is_zero_at_result() {
        let cand = candidate_2017_08_21();
        let conj = solve_conjunction(&cand, ConjunctionKind::EclipticLongitude, config_tight())
            .expect("conjunction must solve for 2017-08-21 candidate");
        let d = delta_lambda(conj.time_tt).abs();
        assert!(
            d < 5e-6,
            "|Δλ| at conjunction = {d} rad (~{}″) not ≈ 0",
            d * 180.0 * 3600.0 / std::f64::consts::PI
        );
    }

    /// 契約5（赤経合）: 返り値 time_tt で赤経差 |Δα| がほぼ 0（独立オラクル `delta_alpha`）。
    #[test]
    fn solve_conjunction_right_ascension_is_zero_at_result() {
        let cand = candidate_2017_08_21();
        let conj = solve_conjunction(&cand, ConjunctionKind::RightAscension, config_tight())
            .expect("RA conjunction must solve for 2017-08-21 candidate");
        let d = delta_alpha(conj.time_tt).abs();
        assert!(
            d < 5e-6,
            "|Δα| at conjunction = {d} rad (~{}″) not ≈ 0",
            d * 180.0 * 3600.0 / std::f64::consts::PI
        );
    }

    /// 契約6: time_tt が候補の検索窓 `[start, end]`（TT-JD）内にある（両 kind）。
    #[test]
    fn solve_conjunction_result_is_within_candidate_window() {
        let cand = candidate_2017_08_21();
        let (lo, hi) = window_jd(&cand);
        for kind in [
            ConjunctionKind::EclipticLongitude,
            ConjunctionKind::RightAscension,
        ] {
            let conj =
                solve_conjunction(&cand, kind, config_tight()).expect("conjunction must solve");
            let t = jd_of(conj.time_tt);
            assert!(
                lo <= t && t <= hi,
                "{kind:?}: time_tt {t} not within window [{lo}, {hi}]"
            );
        }
    }

    /// 契約7（独立オラクル, 離角極小）: 合時刻は太陽-月見かけ離角の極小時刻の近傍にある（= 同じ朔・
    /// 正しいイベントを解いている独立確認）。**合（Δλ/Δα=0）と離角極小（角距離最接近）は厳密一致しない**:
    /// 月が黄緯/赤緯方向に動くため最接近は合から ±十数分ずれる（2017 は月が交点付近で緯度運動が速く、
    /// 赤経合で実測 ~12.7 分）。よって閾値は物理的に妥当な 20 分とする（合のゼロ点厳密性は契約5の
    /// |Δ|≈0 が別途担保。本テストは「正しい朔を捕まえているか」の粗い独立チェック）。
    #[test]
    fn solve_conjunction_time_is_near_elongation_minimum() {
        let cand = candidate_2017_08_21();
        let (lo, hi) = window_jd(&cand);
        let min_jd = elongation_min_jd(lo, hi, 5.0 / 1440.0); // 5 分刻み（閾値 20 分に対し十分）
        for kind in [
            ConjunctionKind::EclipticLongitude,
            ConjunctionKind::RightAscension,
        ] {
            let conj =
                solve_conjunction(&cand, kind, config_tight()).expect("conjunction must solve");
            let diff_days = (jd_of(conj.time_tt) - min_jd).abs();
            assert!(
                diff_days < 20.0 / 1440.0,
                "{kind:?}: conjunction time differs from elongation minimum by {} min (>20)",
                diff_days * 1440.0
            );
        }
    }

    /// 契約8: 2017-08-21 の皆既日食候補では separation が小さい（< ~1.5° = 0.026 rad, 皆既ゆえ近接）。
    /// かつ separation はテスト側で独立計算した合時刻離角（acos クランプ）に一致する。
    /// さらに合時刻が 2017-08-21 18:30 UTC 付近にあること（TT-JD ≈ 2457987.27 + ΔT≈8e-4 日）。
    #[test]
    fn solve_conjunction_2017_total_eclipse_has_small_separation() {
        let cand = candidate_2017_08_21();
        let conj = solve_conjunction(&cand, ConjunctionKind::EclipticLongitude, config_tight())
            .expect("conjunction must solve");
        // 皆既日食 → 合時刻の月-太陽角距離は小さい。
        assert!(
            conj.separation.0 < 0.026,
            "separation {} rad (~{}°) not < 1.5° for total eclipse",
            conj.separation.0,
            conj.separation.0 * 180.0 / std::f64::consts::PI
        );
        // separation は合時刻の独立離角（acos クランプ）と一致するはず。
        let oracle = elongation_rad(conj.time_tt);
        assert!(
            close(conj.separation.0, oracle, 1e-6),
            "separation {} != independent elongation {oracle} at conjunction",
            conj.separation.0
        );
        // 合時刻は 2017-08-21 18:30 UTC 付近（UTC-JD≈2457987.27, TT は +ΔT≈8e-4 日）。
        // 窓半幅 ~1 日に対し ±0.05 日（≈1.2h）で固定すれば年代・朔の取り違えを殺せる。
        let approx_tt_jd = 2_457_987.27 + 69.0 / 86_400.0;
        assert!(
            close(jd_of(conj.time_tt), approx_tt_jd, 0.05),
            "conjunction time {} not near 2017-08-21 18:30 UTC",
            jd_of(conj.time_tt)
        );
    }

    /// 契約9: 両 kind がともに解け、time_tt が互いに ~30 分以内（黄経合と赤経合の定義差）。
    #[test]
    fn solve_conjunction_both_kinds_agree_within_30_minutes() {
        let cand = candidate_2017_08_21();
        let ecl = solve_conjunction(&cand, ConjunctionKind::EclipticLongitude, config_tight())
            .expect("ecliptic conjunction must solve");
        let ra = solve_conjunction(&cand, ConjunctionKind::RightAscension, config_tight())
            .expect("RA conjunction must solve");
        // kind フィールドが取り違えられていないこと。
        assert_eq!(ecl.kind, ConjunctionKind::EclipticLongitude);
        assert_eq!(ra.kind, ConjunctionKind::RightAscension);
        let diff_days = (jd_of(ecl.time_tt) - jd_of(ra.time_tt)).abs();
        assert!(
            diff_days < 30.0 / 1440.0,
            "ecliptic vs RA conjunction differ by {} min (>30)",
            diff_days * 1440.0
        );
    }

    /// 契約10（収束許容）: x_tolerance_days を小さく設定しても合が解ける（より厳しい設定でも収束する）。
    /// 緩い設定と厳しい設定で得た合時刻が x_tolerance に整合する範囲で一致する。
    #[test]
    fn solve_conjunction_converges_under_tight_tolerance() {
        let cand = candidate_2017_08_21();
        let loose = RootConfig {
            x_tolerance_days: 1e-4,
            max_iterations: 100,
        };
        let tight = RootConfig {
            x_tolerance_days: 1e-9,
            max_iterations: 200,
        };
        let c_loose = solve_conjunction(&cand, ConjunctionKind::EclipticLongitude, loose)
            .expect("loose conjunction must solve");
        let c_tight = solve_conjunction(&cand, ConjunctionKind::EclipticLongitude, tight)
            .expect("tight conjunction must solve");
        // 両者の差は緩い側の許容（1e-4 日）程度に収まるはず。
        let diff_days = (jd_of(c_loose.time_tt) - jd_of(c_tight.time_tt)).abs();
        assert!(
            diff_days < 1e-3,
            "loose vs tight conjunction differ by {diff_days} days (tolerance not honored?)"
        );
        // 厳しい設定では角度差がより 0 に近い（独立オラクル）。
        assert!(
            delta_lambda(c_tight.time_tt).abs() < 5e-6,
            "tight |Δλ| not ≈ 0"
        );
    }

    /// 健全性: 返り値の time_tt と separation が有限（NaN/Inf を返さない）。
    #[test]
    fn solve_conjunction_results_are_finite() {
        let cand = candidate_2017_08_21();
        for kind in [
            ConjunctionKind::EclipticLongitude,
            ConjunctionKind::RightAscension,
        ] {
            let conj =
                solve_conjunction(&cand, kind, config_tight()).expect("conjunction must solve");
            assert!(
                jd_of(conj.time_tt).is_finite(),
                "{kind:?}: time_tt non-finite"
            );
            assert!(
                conj.separation.0.is_finite(),
                "{kind:?}: separation non-finite"
            );
            // separation は角距離なので非負。
            assert!(
                conj.separation.0 >= 0.0,
                "{kind:?}: separation {} negative",
                conj.separation.0
            );
        }
    }

    /// 機構と合の整合（メタ）: 合時刻 time_tt で、種別に対応する角度差の符号が窓端で反転している
    /// （f(t0)<0<f(t1) もしくは逆）。実装がブラケットを取れる前提条件が実 ephemeris で満たされること。
    /// これは「窓内に必ず合がある」という ISSUE-016 契約の独立確認でもある。
    #[test]
    fn angle_difference_changes_sign_across_window_for_both_kinds() {
        let cand = candidate_2017_08_21();
        let lo = cand.search_window.start;
        let hi = cand.search_window.end;
        for kind in [
            ConjunctionKind::EclipticLongitude,
            ConjunctionKind::RightAscension,
        ] {
            let f_lo = angle_difference(kind, lo);
            let f_hi = angle_difference(kind, hi);
            assert!(
                f_lo * f_hi < 0.0,
                "{kind:?}: angle difference does not change sign across window (f_lo={f_lo}, f_hi={f_hi})"
            );
        }
    }
}
