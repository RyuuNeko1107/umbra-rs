//! 全球接触 P1/U1/U4/P4 時刻求解（`docs/algorithms/08-global.md` §C, ISSUE-023 S6b-i,
//! Explanatory Supplement to the Astronomical Almanac §11 / Meeus *Astronomical Algorithms* Ch.54）。
//!
//! 影錐の縁が地球縁に最初/最後に外接する全球接触の**時刻**を求める（地球を 1 天体として扱い、
//! 観測者は登場しない＝局地接触 `local_contacts.rs` の全球版）:
//! - 半影縁の外接（P1=部分食開始 / P4=部分食終了）: `x² + y² = (ρ_g + l1)²`。
//! - 本影/反本影縁の外接（U1=中心食開始 / U4=中心食終了）: `x² + y² = (ρ_g + |l2|)²`。
//!
//! `ρ_g` は地球縁の基本面射影半径（球近似 1.0。扁平込み有効半径は要確認, algorithms 08 §B/§C）。
//! 求根対象は連続関数 `g_P(t) = (x²+y²) − (ρ_g+l1)²` / `g_U(t) = (x²+y²) − (ρ_g+|l2|)²`。探索窓
//! （`source.fit_interval()`）を粗走査で符号変化区間に分割 → 各区間を Brent（ISSUE-008, 無条件
//! Newton 禁止・conventions §11）。**U1/U4 は中心食でなければ None**（本影が地球縁に届かず符号変化なし）。
//!
//! 本層は接触**時刻のみ**（TT/UTC）。地表点（`GlobalContact.position`）は S6b-ii で付与する。

// solve_global_contacts（pub(crate)）は ISSUE-043 S6c（classify_global 結線）が消費するまで未使用。
#![allow(dead_code)]
// 粗走査の分割数（窓秒数 / 刻み）の整数⇔f64 変換は小さな添字のみ（local_contacts.rs と同様）。
#![allow(
    clippy::cast_precision_loss,
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss
)]

use umbra_core::time::tt_to_utc;
use umbra_core::{brent_root, JulianDate2, TtInstant};

use crate::conjunction::RootConfig;
use crate::error::EclipseError;
use crate::local_contacts::ContactInstant;
use crate::source::BesselianSource;

/// 地球縁の基本面射影半径 ρ_g（Re）。球近似 1.0（扁平込み有効半径は要確認, algorithms 08 §B/§C）。
const EARTH_LIMB_RADIUS: f64 = 1.0;
/// 粗走査の刻み（SI 秒）。全球接触対（特に U1/U4＝中心食帯の出入り）を取りこぼさない細かさ
/// （偽陰性不可・architecture §3）。local_contacts.rs と同値。
const CONTACT_SCAN_STEP_SECONDS: f64 = 30.0;
/// 1 日 = 86400 SI 秒。
const SECONDS_PER_DAY: f64 = 86_400.0;

/// 全球接触の時刻集合（P1/U1/U4/P4, TT+UTC）。U1/U4 は中心食でなければ `None`。
/// 地表点（`GlobalContact.position`）は S6b-ii で付与する（本層は時刻のみ）。
#[derive(Clone, Copy, Debug, PartialEq, Default)]
pub(crate) struct GlobalContactTimes {
    /// P1: 部分食開始（半影縁が地球縁に最初に外接）。
    pub p1: Option<ContactInstant>,
    /// U1: 中心食開始（本影縁が地球縁に最初に外接, 中心食のみ）。
    pub u1: Option<ContactInstant>,
    /// U4: 中心食終了（本影縁が地球縁に最後に外接, 中心食のみ）。
    pub u4: Option<ContactInstant>,
    /// P4: 部分食終了（半影縁が地球縁に最後に外接）。
    pub p4: Option<ContactInstant>,
}

/// 全球接触 P1/U1/U4/P4 の時刻を `source.fit_interval()` 内で求解する。
///
/// - `source`: 瞬時ベッセル要素の供給源（ISSUE-037, 既定 `DirectBesselianSource`）。
/// - `config`: Brent 求根設定（root_tolerance は目標の 1/10 以下）。
///
/// U1/U4 は中心食でなければ `None`。日食でない（半影も地球縁に届かない）窓では P1/P4 も `None`。
/// 本影走査は P1・P4 の**両方**が見つかった場合のみ行う（窓が部分食区間を片側で切ると P1 または
/// P4 のみ Some・U1/U4 は走査せず None になりうる。通常は窓が全接触を内部に括る前提）。
pub(crate) fn solve_global_contacts<B: BesselianSource>(
    source: &B,
    config: RootConfig,
) -> Result<GlobalContactTimes, EclipseError> {
    let iv = source.fit_interval();
    let t0 = iv.start.jd2().jd();
    let t1 = iv.end.jd2().jd();

    // 半影外接 P1/P4: g_P = (x²+y²) − (ρ_g+l1)²。+→−（外側→内側）= P1, −→+（内側→外側）= P4。
    let penumbra = scan_sign_change_roots(|jd| g_penumbra(source, jd), t0, t1, config)?;
    let p1_jd = penumbra.iter().find(|r| !r.ascending).map(|r| r.time_jd);
    let p4_jd = penumbra
        .iter()
        .rev()
        .find(|r| r.ascending)
        .map(|r| r.time_jd);

    // 本影外接 U1/U4: g_U = (x²+y²) − (ρ_g+|l2|)²。中心食でなければ符号変化なし → None。
    // 部分食区間 [P1,P4] に絞って走査する（無駄な走査と取りこぼしを避ける, local_contacts と同方針）。
    let (mut u1_jd, mut u4_jd) = (None, None);
    if let (Some(a), Some(b)) = (p1_jd, p4_jd) {
        let umbra = scan_sign_change_roots(|jd| g_umbra(source, jd), a, b, config)?;
        u1_jd = umbra.iter().find(|r| !r.ascending).map(|r| r.time_jd);
        u4_jd = umbra.iter().rev().find(|r| r.ascending).map(|r| r.time_jd);
    }

    Ok(GlobalContactTimes {
        p1: contact_at(p1_jd)?,
        u1: contact_at(u1_jd)?,
        u4: contact_at(u4_jd)?,
        p4: contact_at(p4_jd)?,
    })
}

/// 粗走査で見つけた符号変化の根（TT-JD）と、その点で関数が「−→+（昇）」か「+→−（降）」か。
#[derive(Clone, Copy, Debug)]
struct SignChangeRoot {
    time_jd: f64,
    /// `true`: 関数が負→正（昇）。`false`: 正→負（降）。
    ascending: bool,
}

/// 単一 TT-JD から TtInstant（求根は単一 f64 JD 空間。local_contacts.rs と同方針）。
fn tt(jd: f64) -> TtInstant {
    TtInstant::from_jd2(JulianDate2::from_jd(jd))
}

/// 半影外接 `g_P(t) = (x²+y²) − (ρ_g+l1)²`（Expl. Suppl. §11 / Meeus Ch.54）。
/// gamma=√(x²+y²) が半影縁 ρ_g+l1 と一致する時刻が P1/P4（地球を 1 天体として扱う全球版）。
fn g_penumbra<B: BesselianSource>(source: &B, jd: f64) -> Result<f64, EclipseError> {
    let e = source.at(tt(jd))?;
    let r = EARTH_LIMB_RADIUS + e.l1;
    Ok(e.x * e.x + e.y * e.y - r * r)
}

/// 本影/反本影外接 `g_U(t) = (x²+y²) − (ρ_g+|l2|)²`（l2 は符号付き → |l2|）。
fn g_umbra<B: BesselianSource>(source: &B, jd: f64) -> Result<f64, EclipseError> {
    let e = source.at(tt(jd))?;
    let r = EARTH_LIMB_RADIUS + e.l2.abs();
    Ok(e.x * e.x + e.y * e.y - r * r)
}

/// 粗走査の分割数（窓幅 / 刻み, 最低 2）。**走査解像度のみ**を決め、接触検出の正否には影響しない
/// （n が十分大きければ符号変化を捉える）。`span`/`n` の算術変異は等価（細かい n でも根は同じ）/timeout
/// （巨大 n）になるため `mutation.yml` で `--exclude-re 'in scan_point_count'` 除外する
/// （docs/reviews/mutation-local-contacts.md と同方針）。
fn scan_point_count(t0_jd: f64, t1_jd: f64, step_seconds: f64) -> usize {
    let span_seconds = (t1_jd - t0_jd) * SECONDS_PER_DAY;
    (span_seconds / step_seconds).ceil().max(2.0) as usize
}

/// `[t0_jd, t1_jd]`（TT-JD）を一定刻みで粗走査し、`f` の符号変化区間を Brent で精解して全根を返す。
/// 接触が無ければ空 Vec（その窓に該当接触なし＝食なし/部分食は正常, エラーにしない）。
/// 偽陰性回避のため刻みは [`CONTACT_SCAN_STEP_SECONDS`]。無条件 Newton 禁止（conventions §11）。
fn scan_sign_change_roots<F>(
    mut f: F,
    t0_jd: f64,
    t1_jd: f64,
    config: RootConfig,
) -> Result<Vec<SignChangeRoot>, EclipseError>
where
    F: FnMut(f64) -> Result<f64, EclipseError>,
{
    let n = scan_point_count(t0_jd, t1_jd, CONTACT_SCAN_STEP_SECONDS);

    let mut roots = Vec::new();
    let mut prev_jd = t0_jd;
    let mut prev_f = f(prev_jd)?;
    for i in 1..=n {
        let frac = i as f64 / n as f64;
        let cur_jd = t0_jd + (t1_jd - t0_jd) * frac;
        let cur_f = f(cur_jd)?;
        if prev_f * cur_f < 0.0 {
            let ascending = cur_f > prev_f;
            // Brent でサブ区間 [prev_jd, cur_jd] を精解。f のエラーは捕捉して伝播する。
            let mut eval_err: Option<EclipseError> = None;
            let root = brent_root(
                &mut |jd| match f(jd) {
                    Ok(v) => v,
                    Err(e) => {
                        eval_err = Some(e);
                        0.0
                    }
                },
                prev_jd,
                cur_jd,
                config.x_tolerance_days,
                config.max_iterations,
            );
            if let Some(e) = eval_err {
                return Err(e);
            }
            roots.push(SignChangeRoot {
                time_jd: root?,
                ascending,
            });
        }
        prev_jd = cur_jd;
        prev_f = cur_f;
    }
    Ok(roots)
}

/// TT-JD（あれば）から `ContactInstant`（TT + UTC）を作る。UTC は ΔT 経由（ISSUE-007）。
fn contact_at(jd: Option<f64>) -> Result<Option<ContactInstant>, EclipseError> {
    match jd {
        Some(jd) => {
            let time_tt = tt(jd);
            let time_utc = tt_to_utc(time_tt)?;
            Ok(Some(ContactInstant { time_tt, time_utc }))
        }
        None => Ok(None),
    }
}

#[cfg(test)]
mod tests {
    //! ISSUE-023 S6b-i 受け入れテスト（strict）。
    //!
    //! ## オラクル戦略（追認回避）
    //! 主オラクルは **独立内部再計算**: 全球接触の定義関数
    //!   `g_P(t) = (x²+y²) − (ρ_g + l1)²`（P1/P4・半影縁の外接）
    //!   `g_U(t) = (x²+y²) − (ρ_g + |l2|)²`（U1/U4・本影縁の外接, |l2| 符号付き）
    //! を、テスト側で `source.at(t)`（ISSUE-037, 検証済プリミティブ）から **別経路** で
    //! 組み直し、`solve_global_contacts` が返した各接触時刻でその g がほぼゼロになることを縛る。
    //! ここで使う供給源は solver の内部関数ではなく公開プリミティブなので、実装本体をコピーした
    //! 「追認」にはならない（solver 本体は g の合成・粗走査・Brent・接触割当を担う）。
    //! `ρ_g = 1.0`（地球縁の球近似射影半径, 扁平込みは要確認で deferred）。
    //!
    //! 幾何の出典: Explanatory Supplement to the Astronomical Almanac §11 /
    //! Meeus *Astronomical Algorithms* 2nd ed. Ch.54（global.rs の greatest/classify と整合）。
    //!
    //! 補助オラクル: (2) 接触条件ゼロ点（tight）、(3) 区間内外の符号（独立, P1<midpoint<P4 で
    //! g_P<0・窓外で g_P>0）、(4) TT/UTC 整合 time_utc == tt_to_utc(time_tt)、
    //! (5) NASA ballpark（時刻の桁・時帯のみ・絶対基準にしない）、(6) 非中心 ⇒ U1/U4=None・
    //! P1/P4=Some（時変合成供給源, gamma_min を本影/半影限界の間に置く）、
    //! (7) 食なし窓 ⇒ 全 None。

    use super::*;

    use umbra_core::constants::{EARTH_EQUATORIAL_RADIUS_M, SOLAR_RADIUS_KM};
    use umbra_core::{EspenakMeeusDeltaT, JulianDate2, Radians, TimeInterval, TtInstant};

    use crate::besselian::InstantaneousBesselianElements;
    use crate::source::{BesselianSource, DirectBesselianSource};

    /// 1 日 = 86400 SI 秒。秒↔日換算に使う。
    const SECONDS_PER_DAY: f64 = 86_400.0;
    /// 太陽物理半径[km]（global.rs / local_contacts.rs テストと同一）。
    const R_SUN: f64 = SOLAR_RADIUS_KM;
    /// 月半径[km]（k·Re, IAU 慣習 k=0.2725076・同上）。
    const R_MOON: f64 = 0.2725076 * (EARTH_EQUATORIAL_RADIUS_M / 1000.0);
    /// 地球縁の射影半径（球近似 1.0・spec §C / global_contacts ドキュメント）。
    const RHO_G: f64 = 1.0;

    /// 2017-08-21 中心皆既を内部に括る探索窓（TT-JD 2457986〜2457988, global.rs と同形）。
    /// 最大食 TT-JD≈2457986.768 が内部にあり、P1/U1/U4/P4 全接触を覆う。
    fn window_2017() -> TimeInterval<TtInstant> {
        TimeInterval {
            start: TtInstant::from_jd2(JulianDate2::new(2_457_986.0, 0.0)),
            end: TtInstant::from_jd2(JulianDate2::new(2_457_988.0, 0.0)),
        }
    }

    /// 単一 TT-JD（f64）から TtInstant。
    fn tt_jd(jd: f64) -> TtInstant {
        TtInstant::from_jd2(JulianDate2::from_jd(jd))
    }

    /// TtInstant の TT-JD。
    fn jd_of(t: TtInstant) -> f64 {
        t.jd2().jd()
    }

    /// 2017-08-21 中心皆既を括る `DirectBesselianSource`（fit_interval=探索窓）。
    fn make_source_2017<'d>(
        dt: &'d EspenakMeeusDeltaT,
    ) -> DirectBesselianSource<'d, EspenakMeeusDeltaT> {
        DirectBesselianSource::new(R_SUN, R_MOON, dt, window_2017())
    }

    /// 標準の Brent 設定。x_tolerance_days = 1e-9 日（≈8.6e-5 s）は接触 ±2s 目標の
    /// 1/10（≤0.2s）を十分に下回る。max_iterations は余裕を見て 200（local_contacts と同一）。
    fn config_tight() -> RootConfig {
        RootConfig {
            x_tolerance_days: 1e-9,
            max_iterations: 200,
        }
    }

    // ============================================================
    // 独立オラクル: 接触条件関数 g_P / g_U（テスト側で別経路再構成）
    // ============================================================

    /// 半影外接 `g_P(t) = (x²+y²) − (ρ_g + l1)²`（P1/P4）を独立に評価する。
    /// `source.at` は検証済プリミティブ。solver 内部関数には依存しない。
    fn g_penumbra<B: BesselianSource>(source: &B, t: TtInstant) -> f64 {
        let e = source
            .at(t)
            .expect("source.at should succeed near 2017 eclipse");
        let r = RHO_G + e.l1;
        (e.x * e.x + e.y * e.y) - r * r
    }

    /// 本影外接 `g_U(t) = (x²+y²) − (ρ_g + |l2|)²`（U1/U4, l2 符号付き → |l2|）。
    fn g_umbra<B: BesselianSource>(source: &B, t: TtInstant) -> f64 {
        let e = source
            .at(t)
            .expect("source.at should succeed near 2017 eclipse");
        let r = RHO_G + e.l2.abs();
        (e.x * e.x + e.y * e.y) - r * r
    }

    /// 接触時刻 t で g がほぼゼロであることの許容。
    /// g = (x²+y²)−R²、接触付近 dg/dt ≈ 2·gamma·d(gamma)/dt ~ O(1) Re²/day。
    /// root_tolerance ≤0.2s = 2.3e-6 day に対し残差 |g| ≲ |dg/dt|·tol ~ 数×1e-6 Re²。
    /// 期待残差の一桁上の 1e-5 Re² を安全側の上限とする（local_contacts.rs と同値）。
    const G_ROOT_TOL: f64 = 1e-5;

    // ============================================================
    // 非中心 / 食なし 用の時変合成供給源（global.rs の AxisMissSource と同形）
    // ============================================================

    /// gamma=√(x²+y²)=|x| が中心で内部極小を持つ時変合成供給源。
    /// x = x_min + K·(jd−center)², y=0, l1/l2 固定。極小値 x_min を選んで
    /// 半影は届く（g_P<0）が本影は届かない（g_U>0）等の状況を作る。
    struct SyntheticGammaSource {
        center_jd: f64,
        x_min: f64,
        k: f64,
        l1: f64,
        l2: f64,
        window: TimeInterval<TtInstant>,
    }

    impl BesselianSource for SyntheticGammaSource {
        fn at(&self, t: TtInstant) -> Result<InstantaneousBesselianElements, EclipseError> {
            let dj = t.jd2().jd() - self.center_jd;
            Ok(InstantaneousBesselianElements {
                x: self.x_min + self.k * dj * dj,
                y: 0.0,
                declination: Radians(0.0),
                mu: Radians(0.0),
                l1: self.l1,
                l2: self.l2,
                tan_f1: 0.0047,
                tan_f2: 0.0046,
                time_tt: t,
            })
        }

        fn fit_interval(&self) -> TimeInterval<TtInstant> {
            self.window
        }
    }

    /// center ± half_day を窓に持つ合成供給源を作る。極小（=x_min の gamma）が内部に来る。
    fn synthetic_source(
        x_min: f64,
        k: f64,
        l1: f64,
        l2: f64,
        half_day: f64,
    ) -> SyntheticGammaSource {
        let center_jd = 2_457_986.768;
        SyntheticGammaSource {
            center_jd,
            x_min,
            k,
            l1,
            l2,
            window: TimeInterval {
                start: tt_jd(center_jd - half_day),
                end: tt_jd(center_jd + half_day),
            },
        }
    }

    // ============================================================
    // テスト本体
    // ============================================================

    /// 主オラクル(1)（全接触存在＋順序）: 2017-08-21 中心皆既で P1/U1/U4/P4 が全て Some、
    /// かつ time_tt が P1 < U1 < U4 < P4（半影は本影より広く先に触れ後に離れる）。
    #[test]
    fn central_2017_all_four_present_and_ordered() {
        let dt = EspenakMeeusDeltaT;
        let src = make_source_2017(&dt);
        let times = solve_global_contacts(&src, config_tight())
            .expect("2017 central eclipse should yield global contacts");

        let p1 = times.p1.expect("P1 must be Some (partial begins)");
        let u1 = times.u1.expect("U1 must be Some (central begins)");
        let u4 = times.u4.expect("U4 must be Some (central ends)");
        let p4 = times.p4.expect("P4 must be Some (partial ends)");

        let (p1j, u1j, u4j, p4j) = (
            jd_of(p1.time_tt),
            jd_of(u1.time_tt),
            jd_of(u4.time_tt),
            jd_of(p4.time_tt),
        );
        assert!(p1j < u1j, "expected P1 < U1, got {p1j} !< {u1j}");
        assert!(u1j < u4j, "expected U1 < U4, got {u1j} !< {u4j}");
        assert!(u4j < p4j, "expected U4 < P4, got {u4j} !< {p4j}");
    }

    /// 主オラクル(2)（接触条件ゼロ点・tight）: 返った各時刻で独立再構成した g がほぼゼロ。
    /// P1/P4 で g_P ≈ 0、U1/U4 で g_U ≈ 0（|g| < 1e-6 Re²）。接触の **定義** を直接縛る。
    #[test]
    fn contacts_are_zeros_of_independent_contact_functions() {
        let dt = EspenakMeeusDeltaT;
        let src = make_source_2017(&dt);
        let times = solve_global_contacts(&src, config_tight())
            .expect("2017 central eclipse should yield global contacts");

        let p1 = times.p1.expect("P1");
        let p4 = times.p4.expect("P4");
        for (label, c) in [("P1", p1), ("P4", p4)] {
            let g = g_penumbra(&src, c.time_tt);
            assert!(
                g.abs() < G_ROOT_TOL,
                "{label}: g_P({}) = {g} not ≈ 0 (tol {G_ROOT_TOL})",
                jd_of(c.time_tt)
            );
        }

        let u1 = times.u1.expect("U1");
        let u4 = times.u4.expect("U4");
        for (label, c) in [("U1", u1), ("U4", u4)] {
            let g = g_umbra(&src, c.time_tt);
            assert!(
                g.abs() < G_ROOT_TOL,
                "{label}: g_U({}) = {g} not ≈ 0 (tol {G_ROOT_TOL})",
                jd_of(c.time_tt)
            );
        }
    }

    /// 主オラクル(3)（区間内外の符号・独立）: P1/P4 の中点で g_P < 0（半影接触区間の内側）、
    /// [P1,P4] の外側（P1−300s / P4+300s, 窓内へクランプ）で g_P > 0。U1/U4 についても同様に
    /// 中点で g_U < 0・外側で g_U > 0。接触が区間境界であることを符号で独立に縛る。
    #[test]
    fn contact_intervals_have_correct_signs() {
        let dt = EspenakMeeusDeltaT;
        let src = make_source_2017(&dt);
        let w = window_2017();
        let (lo, hi) = (jd_of(w.start), jd_of(w.end));
        let times = solve_global_contacts(&src, config_tight())
            .expect("2017 central eclipse should yield global contacts");

        let pad = 300.0 / SECONDS_PER_DAY; // 300 秒（粗走査刻みより十分大）

        // 半影 P1/P4。
        let p1 = jd_of(times.p1.expect("P1").time_tt);
        let p4 = jd_of(times.p4.expect("P4").time_tt);
        let mid_p = 0.5 * (p1 + p4);
        assert!(
            g_penumbra(&src, tt_jd(mid_p)) < 0.0,
            "g_P at P1/P4 midpoint must be < 0 (inside penumbra-touch interval)"
        );
        let before_p1 = (p1 - pad).max(lo);
        let after_p4 = (p4 + pad).min(hi);
        assert!(
            g_penumbra(&src, tt_jd(before_p1)) > 0.0,
            "g_P just before P1 must be > 0 (outside penumbra interval)"
        );
        assert!(
            g_penumbra(&src, tt_jd(after_p4)) > 0.0,
            "g_P just after P4 must be > 0 (outside penumbra interval)"
        );

        // 本影 U1/U4（中心食）。
        let u1 = jd_of(times.u1.expect("U1").time_tt);
        let u4 = jd_of(times.u4.expect("U4").time_tt);
        let mid_u = 0.5 * (u1 + u4);
        assert!(
            g_umbra(&src, tt_jd(mid_u)) < 0.0,
            "g_U at U1/U4 midpoint must be < 0 (inside umbra-touch interval)"
        );
        let before_u1 = (u1 - pad).max(lo);
        let after_u4 = (u4 + pad).min(hi);
        assert!(
            g_umbra(&src, tt_jd(before_u1)) > 0.0,
            "g_U just before U1 must be > 0 (outside umbra interval)"
        );
        assert!(
            g_umbra(&src, tt_jd(after_u4)) > 0.0,
            "g_U just after U4 must be > 0 (outside umbra interval)"
        );
    }

    /// 補助オラクル(4)（TT/UTC 整合）: 全接触で time_utc == tt_to_utc(time_tt)（過去日食は変換可）。
    #[test]
    fn contact_utc_matches_tt_to_utc() {
        let dt = EspenakMeeusDeltaT;
        let src = make_source_2017(&dt);
        let times = solve_global_contacts(&src, config_tight())
            .expect("2017 central eclipse should yield global contacts");

        for (label, c) in [
            ("p1", times.p1),
            ("u1", times.u1),
            ("u4", times.u4),
            ("p4", times.p4),
        ] {
            if let Some(c) = c {
                let want = umbra_core::time::tt_to_utc(c.time_tt)
                    .expect("2017 is post-1972, tt_to_utc must succeed");
                assert_eq!(
                    c.time_utc, want,
                    "{label}: time_utc must equal tt_to_utc(time_tt)"
                );
            }
        }
    }

    /// 補助オラクル(5)（NASA ballpark・緩め・第二義）: 2017-08-21 全球接触の時刻が
    /// 既知の時帯（UTC）に収まる。NASA 公表 P1≈15:46 / U1≈16:48 / U4≈20:02 / P4≈21:04 UTC を
    /// **WIDE な時レンジ**で括る（絶対基準にしない・k/ΔT 慣習差で秒値は再現しない）。
    #[test]
    fn global_contacts_2017_utc_hours_in_nasa_ballpark() {
        let dt = EspenakMeeusDeltaT;
        let src = make_source_2017(&dt);
        let times = solve_global_contacts(&src, config_tight())
            .expect("2017 central eclipse should yield global contacts");

        let hour_of = |c: ContactInstant| {
            let (_y, _mo, _d, hh, mm, _ss) = c.time_utc.to_gregorian();
            f64::from(hh) + f64::from(mm) / 60.0
        };

        let p1 = hour_of(times.p1.expect("P1"));
        let u1 = hour_of(times.u1.expect("U1"));
        let u4 = hour_of(times.u4.expect("U4"));
        let p4 = hour_of(times.p4.expect("P4"));

        assert!(
            (15.0..=17.0).contains(&p1),
            "P1 UTC hour {p1} not in ballpark [15,17] (NASA≈15:46)"
        );
        assert!(
            (16.0..=18.0).contains(&u1),
            "U1 UTC hour {u1} not in ballpark [16,18] (NASA≈16:48)"
        );
        assert!(
            (19.0..=21.0).contains(&u4),
            "U4 UTC hour {u4} not in ballpark [19,21] (NASA≈20:02)"
        );
        assert!(
            (20.0..=22.0).contains(&p4),
            "P4 UTC hour {p4} not in ballpark [20,22] (NASA≈21:04)"
        );
    }

    /// 補助オラクル(6)（非中心 ⇒ U1/U4=None・P1/P4=Some）: 時変合成供給源で gamma の内部極小を
    /// 本影限界と半影限界の **間** に置く。gamma_min = x_min = 1.2、l1=0.54, l2=−0.01:
    ///   1+|l2| = 1.01 < gamma_min(1.2) < 1+l1 = 1.54。
    /// 中心: g_P = 1.2² − 1.54² = 1.44 − 2.3716 < 0（半影は届く）、
    ///       g_U = 1.2² − 1.01² = 1.44 − 1.0201 > 0（本影は届かない）。
    /// 端: x_edge = x_min + k·half_day² = 1.2 + 50·0.15² = 1.2 + 1.125 = 2.325 なので
    ///       g_P = 2.325² − 1.54² = 5.405625 − 2.3716 > 0（半影が縁から離れる ⇒ g_P が窓内でゼロ交差
    ///       ⇒ P1/P4 が存在）、g_U = 2.325² − 1.01² > 0（本影は終始届かない ⇒ U1/U4=None）。
    /// 窓幅は g_P が縁で正になるよう half_day=0.15 を選ぶ（50·half_day² > 0.34 ⇒ half_day > 0.0826）。
    #[test]
    fn non_central_has_no_umbral_contacts_but_has_penumbral() {
        // gamma_min=1.2 は 1+|l2|=1.01 と 1+l1=1.54 の間。窓は center ± 0.15 day（端で半影が縁を離れる）。
        let half_day = 0.15;
        let src = synthetic_source(1.2, 50.0, 0.54, -0.01, half_day);

        // 前提の独立確認(中心): g_P<0（半影到達）・g_U>0（本影未到達）。
        // g_P(center)=1.44−2.3716<0, g_U(center)=1.44−1.0201>0。
        let center_jd = 2_457_986.768;
        let center = tt_jd(center_jd);
        assert!(
            g_penumbra(&src, center) < 0.0,
            "precondition: g_P at gamma_min(1.2) must be < 0 (penumbra reaches limb)"
        );
        assert!(
            g_umbra(&src, center) > 0.0,
            "precondition: g_U at gamma_min(1.2) must be > 0 (umbra never reaches limb)"
        );

        // 前提の独立確認(両端): g_P>0（半影が縁を離れる ⇒ 窓内に P1/P4 のゼロ交差が存在する）。
        // x_edge=2.325 ⇒ g_P(edge)=5.405625−2.3716>0。本影は終始届かないので U1/U4 は None。
        for edge in [tt_jd(center_jd - half_day), tt_jd(center_jd + half_day)] {
            assert!(
                g_penumbra(&src, edge) > 0.0,
                "precondition: g_P at window edge must be > 0 (penumbra off the limb \
                 so P1/P4 exist as zero crossings)"
            );
        }

        let times = solve_global_contacts(&src, config_tight())
            .expect("non-central penumbral eclipse must yield Ok (P1/P4 present)");
        assert!(
            times.u1.is_none() && times.u4.is_none(),
            "non-central: U1/U4 must be None (umbra misses limb), got u1={:?} u4={:?}",
            times.u1,
            times.u4
        );
        assert!(
            times.p1.is_some() && times.p4.is_some(),
            "non-central: P1/P4 must be Some (penumbra touches limb)"
        );
    }

    /// 補助オラクル(7)（食なし窓 ⇒ 全 None）: gamma が全域で 1+l1 を超える合成供給源
    /// （x_min=2.0 ⇒ gamma_min=2.0 > 1+0.54=1.54）。g_P が一度も負にならず符号変化が無いので
    /// P1/U1/U4/P4 全て None（「窓内に該当接触なし＝食なし」は正常で Ok）。
    #[test]
    fn no_eclipse_window_yields_all_none() {
        let src = synthetic_source(2.0, 50.0, 0.54, -0.01, 0.05);

        // 前提: center で g_P>0（半影すら届かない＝食なし）。
        let center = tt_jd(2_457_986.768);
        assert!(
            g_penumbra(&src, center) > 0.0,
            "precondition: g_P at gamma_min(2.0) must be > 0 (no eclipse at all)"
        );

        let times = solve_global_contacts(&src, config_tight())
            .expect("no-eclipse window must be Ok (all None), not Err");
        assert_eq!(
            times,
            GlobalContactTimes::default(),
            "no-eclipse window must yield all-None contacts, got {times:?}"
        );
    }

    /// 主オラクル(A)（接触種別の遷移方向）: 各接触時刻 t の前後 ±ε（8s 相当 day, root_tolerance
    /// 1e-9 day より大・走査刻み 30s より小）で独立オラクル g の符号方向を縛る。中点/窓外の符号テスト
    /// では捉えられない P1↔P4 / U1↔U4 取り違えを弾く（local_contacts の遷移方向テストと同形）:
    ///   P1: g_P(p1−ε)>0 && g_P(p1+ε)<0（外側→内側＝部分食開始）
    ///   P4: g_P(p4−ε)<0 && g_P(p4+ε)>0（内側→外側＝部分食終了）
    ///   U1: g_U(u1−ε)>0 && g_U(u1+ε)<0（本影が地球縁に触れる＝中心食開始）
    ///   U4: g_U(u4−ε)<0 && g_U(u4+ε)>0（本影が地球縁を離れる＝中心食終了）
    #[test]
    fn contacts_have_correct_sign_transitions() {
        let dt = EspenakMeeusDeltaT;
        let src = make_source_2017(&dt);
        let times = solve_global_contacts(&src, config_tight())
            .expect("2017 central eclipse should yield global contacts");

        let eps = 8.0 / SECONDS_PER_DAY;
        let before = |t: TtInstant| tt_jd(jd_of(t) - eps);
        let after = |t: TtInstant| tt_jd(jd_of(t) + eps);

        let p1 = times.p1.expect("P1").time_tt;
        let p4 = times.p4.expect("P4").time_tt;
        let u1 = times.u1.expect("U1").time_tt;
        let u4 = times.u4.expect("U4").time_tt;

        // P1: 部分食開始（外側 g_P>0 → 内側 g_P<0）。
        let p1_b = g_penumbra(&src, before(p1));
        let p1_a = g_penumbra(&src, after(p1));
        assert!(
            p1_b > 0.0 && p1_a < 0.0,
            "P1 must be penumbra +→− (g_P before={p1_b}, after={p1_a})"
        );

        // P4: 部分食終了（内側 g_P<0 → 外側 g_P>0）。
        let p4_b = g_penumbra(&src, before(p4));
        let p4_a = g_penumbra(&src, after(p4));
        assert!(
            p4_b < 0.0 && p4_a > 0.0,
            "P4 must be penumbra −→+ (g_P before={p4_b}, after={p4_a})"
        );

        // U1: 中心食開始（本影 g_U>0 → g_U<0）。
        let u1_b = g_umbra(&src, before(u1));
        let u1_a = g_umbra(&src, after(u1));
        assert!(
            u1_b > 0.0 && u1_a < 0.0,
            "U1 must be umbra +→− (g_U before={u1_b}, after={u1_a})"
        );

        // U4: 中心食終了（本影 g_U<0 → g_U>0）。
        let u4_b = g_umbra(&src, before(u4));
        let u4_a = g_umbra(&src, after(u4));
        assert!(
            u4_b < 0.0 && u4_a > 0.0,
            "U4 must be umbra −→+ (g_U before={u4_b}, after={u4_a})"
        );
    }

    /// 主オラクル(B)（窓幅不変性・偽陰性ガード）: 全接触を内部に括る既定窓（2 day 幅）と、より狭い窓
    /// （greatest ±6h）の 2 通りで解き、P1/U1/U4/P4（time_tt）が 2 秒以内で一致することを縛る。
    /// `solve_global_contacts` は `source.fit_interval()` を窓に使うため、狭い窓は別の
    /// `DirectBesselianSource`（fit_interval=狭め）として作る。2017 の部分食は greatest（TT-JD
    /// ≈2457987.2685）周り ~±2.6h なので ±6h（=±0.25 day）は全接触を内部に含む安全な狭め窓。窓幅で
    /// 粗走査分割数（窓秒/30s）が変わるため、刻み感度の偽陰性を実効的に縛る（local_contacts と同形）。
    #[test]
    fn global_contacts_window_width_invariant() {
        let dt = EspenakMeeusDeltaT;

        // 既定（広い）窓: 2457986〜2457988（P1..P4 を内部に括る）。
        let wide_src = make_source_2017(&dt);
        let wide = solve_global_contacts(&wide_src, config_tight())
            .expect("wide window should yield global contacts");

        // 狭い窓: greatest TT-JD≈2457987.2685（≈18:25 UTC + ΔT）の ±6h（±0.25 day）。
        // 部分食 ~±2.6h を内部に括る安全な狭め窓。
        let greatest = 2_457_986.5 + 0.768_532_222_222_222_2;
        let narrow_iv = TimeInterval {
            start: tt_jd(greatest - 0.25),
            end: tt_jd(greatest + 0.25),
        };
        let narrow_src = DirectBesselianSource::new(R_SUN, R_MOON, &dt, narrow_iv);
        let narrow = solve_global_contacts(&narrow_src, config_tight())
            .expect("narrow window should still yield global contacts");

        let pick = |t: &GlobalContactTimes| {
            [
                jd_of(t.p1.expect("P1").time_tt),
                jd_of(t.u1.expect("U1").time_tt),
                jd_of(t.u4.expect("U4").time_tt),
                jd_of(t.p4.expect("P4").time_tt),
            ]
        };
        let (a, b) = (pick(&wide), pick(&narrow));
        for (k, label) in ["P1", "U1", "U4", "P4"].iter().enumerate() {
            assert!(
                (a[k] - b[k]).abs() * SECONDS_PER_DAY < 2.0,
                "{label}: wide vs narrow window differ by {} s",
                (a[k] - b[k]).abs() * SECONDS_PER_DAY
            );
        }
    }
}
