//! 局地最大食 solver（`docs/issues/ISSUE-026`、`docs/conventions.md` §8,
//! Explanatory Supplement §11 / Meeus *Astronomical Algorithms* Ch.54）。
//!
//! 観測地点で太陽-月の基本面上中心間距離 `m(t) = √((ξ−x)² + (η−y)²)` が最小となる最大食の時刻と
//! その最小値 `m_min` を求める。**D2: 最小化対象は m²=u²+v²（u=ξ−x, v=η−y）に統一**し
//! （√ の中心線尖点の微分特異を回避、ISSUE-025 の g=m²−L² と整合）、正式手法は
//! **d(m²)/dt = 2(u·u'+v·v') = 0 の Brent 求根**（ISSUE-008, 無条件 Newton 禁止・conventions §11）。
//! 粗ブラケット（最小付近の括り出し）は粗走査で行う（黄金分割 ISSUE-009 は併用に降格）。
//!
//! 瞬時要素は時間微分（x',y'）を持たないため、d(m²)/dt は m²(t) の中心差分で評価する
//! （numerical-policy §A2/§A5）。
//!
//! 注: 食分・食面積（ISSUE-027）と高度・方位・可視性（ISSUE-028）は本 issue の非目的。本層は
//! 最大食時刻と m_min の算出が責務で、magnitude/obscuration/contact は ISSUE-027/028/043 で充足する。

// solve_local_maximum（pub(crate)）は ISSUE-043（EclipseEngine 結線）が消費するまで未使用。
#![allow(dead_code)]
// 粗走査の分割数（窓秒数 / 刻み）の整数⇔f64 変換は小さな添字のみ（天文量ではない, conjunction.rs 同様）。
#![allow(
    clippy::cast_precision_loss,
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss
)]

use umbra_core::ellipsoid::GeocentricObserver;
use umbra_core::time::tt_to_utc;
use umbra_core::{
    brent_root, JulianDate2, Radians, SolverError, TimeInterval, TtInstant, UtcInstant,
};

use crate::conjunction::RootConfig;
use crate::error::EclipseError;
use crate::projection::project_observer_to_fundamental;
use crate::source::BesselianSource;

/// 最大食の粗ブラケット走査の刻み（SI 秒）。m²(t) の極小付近を 3 点で括るための粗走査。
/// 部分食継続（数時間）の中の単一谷を確実に括れる細かさ。
const MAXIMUM_SCAN_STEP_SECONDS: f64 = 60.0;
/// d(m²)/dt を中心差分で評価する刻み（SI 秒）。瞬時要素が時間微分を持たないための数値微分
/// （numerical-policy §A2）。極小近傍で符号が安定する程度に小さく、量子化に埋もれない程度に大きく。
const DERIVATIVE_STEP_SECONDS: f64 = 2.0;
/// 1 日 = 86400 SI 秒。
const SECONDS_PER_DAY: f64 = 86_400.0;

/// 局地最大食の結果（時刻 TT+UTC と最小中心間距離）。
///
/// 食分・食面積（ISSUE-027）と高度・方位・可視性（ISSUE-028）は本 issue 非目的のため持たない
/// （ISSUE-027/028/043 で拡張）。`min_separation` は食分算出（ISSUE-027）の入力。
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LocalMaximum {
    /// 最大食の TT 時刻（幾何相対の一級値, conventions §6）。
    pub time_tt: TtInstant,
    /// 最大食の UTC 時刻（TT→UTC は ΔT 経由 ISSUE-007）。
    pub time_utc: UtcInstant,
    /// 最大食時点の基本面上中心間距離 m_min（Re）。
    pub min_separation: f64,
}

/// 観測地点の局地最大食を探索窓内で求解する。
///
/// - `source`: 瞬時ベッセル要素の供給源（ISSUE-037）。
/// - `observer`/`east_longitude`: 観測者（ρsinφ′/ρcosφ′ ＋ 東経正, ISSUE-024）。
/// - `search`: 探索窓（全球最大 ±マージンの TT 区間）。
/// - `config`: Brent 求根設定。
pub(crate) fn solve_local_maximum<B: BesselianSource>(
    source: &B,
    observer: &GeocentricObserver,
    east_longitude: Radians,
    search: TimeInterval<TtInstant>,
    config: RootConfig,
) -> Result<LocalMaximum, EclipseError> {
    let t0 = search.start.jd2().jd();
    let t1 = search.end.jd2().jd();

    let m2 = |jd: f64| m2_at(source, observer, east_longitude, jd);
    // d(m²)/dt の符号関数＝中心差分の**非正規化分子** `m²(jd+h) − m²(jd−h)`（瞬時要素が x',y' を
    // 持たないための数値微分, numerical-policy §A2）。Brent は符号・ゼロ点しか使わず、正の定数 1/(2h)
    // で割っても根は不変なので正規化は省く（`/(2h)` の等価変異を排除し、差分の符号を担う `-` を
    // load-bearing にする）。極小で 0、最小より前で負（減少）、後で正（増加）。
    let h = DERIVATIVE_STEP_SECONDS / SECONDS_PER_DAY;
    let dm2_sign = |jd: f64| -> Result<f64, EclipseError> { Ok(m2(jd + h)? - m2(jd - h)?) };

    // 1. 粗走査で m²(t) の最小サンプルを探し、内部 3 点ブラケット [t_{i-1}, t_{i+1}] を作る。
    let n = scan_point_count(t0, t1, MAXIMUM_SCAN_STEP_SECONDS);
    let jd_at = |i: usize| t0 + (t1 - t0) * (i as f64 / n as f64);

    let mut min_i = 0usize;
    let mut min_val = m2(t0)?;
    for i in 1..=n {
        let v = m2(jd_at(i))?;
        if v < min_val {
            min_val = v;
            min_i = i;
        }
    }
    // 最小が窓端＝窓内に内部極小なし（単調 or 定数=平底）。ブラケットできない（conventions §11）。
    if min_i == 0 || min_i == n {
        return Err(EclipseError::Solver(SolverError::RootNotBracketed));
    }

    // 2. （D2 正式手法）d(m²)/dt = 0 を Brent 求根（ISSUE-008, 無条件 Newton 禁止）。
    //    粗ブラケット端で導関数が異符号（減少→増加）でなければブラケット不成立。
    let bracket_lo = jd_at(min_i - 1);
    let bracket_hi = jd_at(min_i + 1);
    let mut eval_err: Option<EclipseError> = None;
    let root = brent_root(
        &mut |jd| match dm2_sign(jd) {
            Ok(v) => v,
            Err(e) => {
                eval_err = Some(e);
                0.0
            }
        },
        bracket_lo,
        bracket_hi,
        config.x_tolerance_days,
        config.max_iterations,
    );
    if let Some(e) = eval_err {
        return Err(e);
    }
    let time_jd = root?;

    // 3. 最大食時点の m_min と TT/UTC。
    let time_tt = tt(time_jd);
    let time_utc = tt_to_utc(time_tt)?;
    let min_separation = m2(time_jd)?.max(0.0).sqrt();

    Ok(LocalMaximum {
        time_tt,
        time_utc,
        min_separation,
    })
}

/// 粗ブラケットの分割数（窓幅 / 刻み, 最低 2）。**走査解像度のみ**を決め、最大食の検出正否には
/// 影響しない（n が十分大きければ単一谷を括れる）。`span`/`n` の算術には振る舞い契約が無く、ここの
/// 算術変異は等価（細かい n でも同じ最大食）/timeout（巨大 n）になるため `mutation.yml` で
/// `--exclude-re 'in scan_point_count'` 除外する（docs/reviews/mutation-local-maximum.md）。
fn scan_point_count(t0_jd: f64, t1_jd: f64, step_seconds: f64) -> usize {
    let span_seconds = (t1_jd - t0_jd) * SECONDS_PER_DAY;
    (span_seconds / step_seconds).ceil().max(2.0) as usize
}

/// 単一 TT-JD から TtInstant（求根は単一 f64 JD 空間。conjunction.rs / local_contacts.rs と同方針）。
fn tt(jd: f64) -> TtInstant {
    TtInstant::from_jd2(JulianDate2::from_jd(jd))
}

/// 観測者基本面座標 (ξ,η) と影軸交点 (x,y) から `m²(t) = (ξ−x)² + (η−y)²`（Re²）を評価する。
/// (ξ,η)=`project_observer_to_fundamental`（ISSUE-024）、(x,y)=`source.at`（ISSUE-037）。
fn m2_at<B: BesselianSource>(
    source: &B,
    observer: &GeocentricObserver,
    east_longitude: Radians,
    jd: f64,
) -> Result<f64, EclipseError> {
    let e = source.at(tt(jd))?;
    let p = project_observer_to_fundamental(observer, east_longitude, &e);
    let du = p.xi - e.x;
    let dv = p.eta - e.y;
    Ok(du * du + dv * dv)
}

#[cfg(test)]
mod tests {
    //! ISSUE-026 受け入れテスト（strict・局地最大食 solver）。
    //!
    //! ## オラクル戦略（追認回避）
    //! 主オラクルは **独立内部再計算による「最小性」検証**。solver の内部手法
    //! （d(m²)/dt=0 の Brent 求根・粗ブラケット・中心差分）には一切依存せず、外部観測可能な
    //! 振る舞い「返った time_tt で m²(t) が局地最小」を、検証済プリミティブ
    //! `project_observer_to_fundamental`（ISSUE-024）＋ `source.at(t)`（ISSUE-037）から
    //! **別経路で組んだ m²(t)** で縛る。射影・供給源は solver 内部関数ではなく公開プリミティブ
    //! なので、実装本体をコピーした「追認」にはならない（solver 本体は m² の合成・粗ブラケット・
    //! Brent 求根・最小判定を担う）。
    //!
    //! 幾何の出典: Explanatory Supplement to the Astronomical Almanac §11 /
    //! Meeus *Astronomical Algorithms* 2nd ed. Ch.54。
    //!   m² = (ξ−x)² + (η−y)²（観測者基本面座標 (ξ,η) と影軸交点 (x,y) の距離²）。
    //!
    //! テスト観点:
    //! (0) 前提メタ: 窓内に厳密な単一局地最小（減少→増加）が存在することを独立走査で確認。
    //! (1) 最小性（最強）: m²(t_max) ≤ m²(t_max±δ)、δ∈{1s,10s,60s} 両側で増加。
    //! (2) min_separation² 整合: 返った min_separation² が独立再計算 m²(t_max) に一致。
    //! (3) d(m²)/dt≈0: 中心差分の導関数を二階微分で割り等価時刻誤差 <0.5s に換算して縛る。
    //! (4) 中心線上: 中心食地点で m_min が小さい（≈0, <0.05 Re）。
    //! (5) 接触の内側: solve_local_contacts の C1/C4 に対し c1 < max < c4。
    //! (6) 窓幅不変: ±3h / ±2h で time_tt が 2 秒内一致。
    //! (7) TT/UTC 整合: time_utc == tt_to_utc(time_tt)。
    //! (8) 異常系: 最小を含まない単調窓（最大食より十分後で m が単調増加）→ Solver(RootNotBracketed)。
    //!     窓内 m² が単調・端で最小＝ d(m²)/dt をブラケットできない（機構的にブラケット失敗）。
    //!     テスト側で「窓内 m² が単調・端で最小」を独立確認してから縛る。DegenerateGeometry は
    //!     影錐の幾何退化専用で兄弟 solver と不整合のため使わない。
    //! (8b) 平底/退化: 定数 m²（時不変供給源）でも panic せず Solver(RootNotBracketed) を返す。
    //!     spec の MockEphemeris 平底受入を、Mock 静的・m が物理的に単一極小ゆえ「定数 m² 退化」
    //!     として再解釈（最小が窓端→ブラケット不成立と同経路）。
    //! (9) 実日食 ballpark（第二義・緩め）: 2017-08-21 中心食地点の最大食が当日 18時台 UTC・
    //!     部分食の中央付近。flaky な厳密秒値は避ける。
    //!
    //! NASA/USNO の地点別最大食時刻・食分は絶対基準にせず、ballpark 整合のみに用いる。

    use super::*;

    use umbra_core::constants::{EARTH_EQUATORIAL_RADIUS_M, SOLAR_RADIUS_KM};
    use umbra_core::ellipsoid::{observer_geocentric, Ellipsoid, GeocentricObserver};
    use umbra_core::{
        EspenakMeeusDeltaT, JulianDate2, Radians, SolverError, TimeInterval, TtInstant,
    };

    use crate::besselian::InstantaneousBesselianElements;
    use crate::local_contacts::solve_local_contacts;
    use crate::projection::project_observer_to_fundamental;
    use crate::source::{BesselianSource, DirectBesselianSource};

    const WGS84: Ellipsoid = Ellipsoid::WGS84;
    const R_SUN: f64 = SOLAR_RADIUS_KM;
    // k·Re（月/地球赤道半径比、IAU 慣習 k=0.2725076）。source.rs / local_contacts.rs と同一。
    const R_MOON: f64 = 0.2725076 * (EARTH_EQUATORIAL_RADIUS_M / 1000.0);

    /// 1 日 = 86400 SI 秒。
    const SECONDS_PER_DAY: f64 = 86_400.0;

    /// TT を 2 要素 JD から構築（source.rs / local_contacts.rs テストと同形）。
    fn tt(jd1: f64, jd2: f64) -> TtInstant {
        TtInstant::from_jd2(JulianDate2::new(jd1, jd2))
    }

    /// 単一 JD（TT）から TtInstant。
    fn tt_jd(jd: f64) -> TtInstant {
        TtInstant::from_jd2(JulianDate2::from_jd(jd))
    }

    /// TtInstant の TT-JD。
    fn jd_of(t: TtInstant) -> f64 {
        t.jd2().jd()
    }

    /// 2017-08-21 皆既日食 最大食付近の TT エポック（besselian.rs / local_contacts.rs と同一）。
    /// 最大食 ≈ 2017-08-21 18:25 UTC（TT-JD ≈ 2457986.5 + 0.76853）。
    fn tt_2017_max() -> TtInstant {
        tt(2_457_986.5, 7.685_322_222_222_222e-1)
    }

    /// この供給源の妥当区間（2017 を含む広め）。`fit_interval` 広告用で at() は外でも評価可。
    fn source_interval() -> TimeInterval<TtInstant> {
        TimeInterval {
            start: tt(2_457_985.0, 0.0),
            end: tt(2_457_989.0, 0.0),
        }
    }

    fn make_source(dt: &EspenakMeeusDeltaT) -> DirectBesselianSource<'_, EspenakMeeusDeltaT> {
        DirectBesselianSource::new(R_SUN, R_MOON, dt, source_interval())
    }

    /// 標準の Brent 設定。x_tolerance_days = 1e-9 日（≈8.6e-5 s）は最大食時刻 ±1〜2s 目標の
    /// 1/10（≤0.2s）を十分に下回る。max_iterations は余裕を見て 200。local_contacts.rs と同一。
    fn config_tight() -> RootConfig {
        RootConfig {
            x_tolerance_days: 1e-9,
            max_iterations: 200,
        }
    }

    /// 最大食 ±half_hours の探索窓（部分食全体を覆う）。local_contacts.rs と同形。
    fn window_around(center: TtInstant, half_hours: f64) -> TimeInterval<TtInstant> {
        let c = jd_of(center);
        let half = half_hours / 24.0;
        TimeInterval {
            start: tt_jd(c - half),
            end: tt_jd(c + half),
        }
    }

    // ============================================================
    // 観測地点（2017-08-21・local_contacts.rs と同一・プローブ確認済）
    // ============================================================
    //
    // 経度は東経正（conventions §3）。北米西経は負で与える（西経吸収は ISSUE-024 で検証済）。

    /// 中心食地点（皆既帯中心付近, 皆既最長 ~2.6 分）。lat=37.5°N, lon=西経89.2°。
    /// 中心線上では最大食で m_min がほぼ 0 になる。主オラクルの基準地点。
    fn central_observer() -> (GeocentricObserver, Radians) {
        let lat = 37.5_f64.to_radians();
        let east_lon = Radians::new((-89.2_f64).to_radians()); // 西経 89.2° → 東経負
        (observer_geocentric(&WGS84, lat, 200.0), east_lon)
    }

    /// 部分食のみの地点（高緯度・本影外）。lat=60°N, lon=西経100°。
    /// 部分食地点でも最大食（幾何最小）は必ず存在する（非 Option）ことの縛りに使う。
    fn partial_observer() -> (GeocentricObserver, Radians) {
        let lat = 60.0_f64.to_radians();
        let east_lon = Radians::new((-100.0_f64).to_radians());
        (observer_geocentric(&WGS84, lat, 200.0), east_lon)
    }

    // ============================================================
    // 独立オラクル: m²(t) をテスト側で別経路再構成
    // ============================================================

    /// 観測者の基本面座標 (ξ,η) と影軸交点 (x,y) から `m²(t) = (ξ−x)² + (η−y)²` を独立に評価する。
    /// `project_observer_to_fundamental`（ISSUE-024 検証済）と `source.at`（ISSUE-037 検証済）から
    /// 組むのみ。solver 内部関数には依存しない（追認回避）。
    /// 出典: Explanatory Supplement §11 / Meeus Ch.54。
    fn m2_at<B: BesselianSource>(
        source: &B,
        observer: &GeocentricObserver,
        east_longitude: Radians,
        t: TtInstant,
    ) -> f64 {
        let e = source
            .at(t)
            .expect("source.at should succeed near 2017 eclipse");
        let p = project_observer_to_fundamental(observer, east_longitude, &e);
        let du = p.xi - e.x;
        let dv = p.eta - e.y;
        du * du + dv * dv
    }

    /// 任意秒オフセットだけ動かした時刻で m² を評価するショートカット。
    fn m2_shift<B: BesselianSource>(
        source: &B,
        observer: &GeocentricObserver,
        east_longitude: Radians,
        t: TtInstant,
        delta_s: f64,
    ) -> f64 {
        let shifted = tt_jd(jd_of(t) + delta_s / SECONDS_PER_DAY);
        m2_at(source, observer, east_longitude, shifted)
    }

    // ============================================================
    // 退化供給源: 時刻によらず一定の瞬時要素を返す（平底/定数 m² 再現用）
    // ============================================================

    /// 時刻に依存せず常に同一の瞬時要素を返す供給源。観測者・経度を固定すれば m²(t) も時刻不変＝定数。
    /// spec の MockEphemeris「皆既帯平底」を、Mock が静的で時変 m²(t) を作れない／m が物理的に単一極小を
    /// 持つ（平坦なのは食分であって m ではない）ことを踏まえ「定数 m² 退化」として再解釈するための装置。
    /// 内部に極小がない（窓全域で d(m²)/dt がゼロ→ブラケット不成立）状況を作る。
    struct ConstantSource {
        elems: InstantaneousBesselianElements,
    }

    impl BesselianSource for ConstantSource {
        fn at(&self, _t: TtInstant) -> Result<InstantaneousBesselianElements, EclipseError> {
            // 時刻ラベルだけ要求時刻にしても m² には影響しない（射影は declination/mu のみ参照）。
            Ok(self.elems)
        }

        fn fit_interval(&self) -> TimeInterval<TtInstant> {
            // 任意の妥当区間（2017 を含む広め・ほかの供給源と同形）。
            source_interval()
        }
    }

    /// 退化供給源の要素を最小構成する（有限ダミー値）。`projection.rs`/`local_contacts.rs` の
    /// `elems_with` 相当。射影で実際に使われるのは declination/mu と影軸交点 x,y のみだが、
    /// 全フィールドに有限値を入れて健全に構成する。
    fn constant_elems() -> InstantaneousBesselianElements {
        InstantaneousBesselianElements {
            x: 0.10,
            y: 0.20,
            declination: Radians::new(0.207),
            mu: Radians::new(1.0),
            l1: 0.54,
            l2: -0.009,
            tan_f1: 0.004_65,
            tan_f2: 0.004_63,
            time_tt: tt(2_451_545.0, 0.0),
        }
    }

    // ============================================================
    // テスト本体
    // ============================================================

    /// 前提の独立確認（メタ・観点0）: 中心食地点の探索窓内で、独立再構成 m²(t) が
    /// 「減少 → 増加」の単一谷を持つ（厳密な局地最小が一意に存在する）。これは後続の
    /// 最小性・接触内側テストの妥当性を支える独立チェック（実装の取りこぼし追認でないこと）。
    #[test]
    fn central_site_window_contains_a_single_interior_minimum() {
        let dt = EspenakMeeusDeltaT;
        let src = make_source(&dt);
        let (obs, lon) = central_observer();
        let w = window_around(tt_2017_max(), 3.0);

        // 窓を細かく走査し m²(t) の符号反転（差分の正→負・負→正）回数を数える。
        let n = 4000;
        let (lo, hi) = (jd_of(w.start), jd_of(w.end));
        let mut sign_changes = 0;
        let mut prev = m2_at(&src, &obs, lon, tt_jd(lo));
        let mut prev_diff_sign = 0.0_f64;
        let mut argmin_jd = lo;
        let mut min_val = prev;
        for k in 1..=n {
            let jd = lo + (hi - lo) * (f64::from(k) / f64::from(n));
            let cur = m2_at(&src, &obs, lon, tt_jd(jd));
            let diff = cur - prev;
            if diff != 0.0 {
                let s = diff.signum();
                if prev_diff_sign != 0.0 && s != prev_diff_sign {
                    sign_changes += 1;
                }
                prev_diff_sign = s;
            }
            if cur < min_val {
                min_val = cur;
                argmin_jd = jd;
            }
            prev = cur;
        }
        // 単一谷なら傾きの符号反転は 1 回（負→正）のみ。
        assert_eq!(
            sign_changes, 1,
            "central site: m²(t) should have exactly one interior minimum (slope sign change), got {sign_changes}"
        );
        // 最小は窓の内部（両端から十分離れている）。
        assert!(
            argmin_jd > lo + (hi - lo) * 0.05 && argmin_jd < hi - (hi - lo) * 0.05,
            "interior minimum should be away from window edges (argmin_jd={argmin_jd})"
        );
    }

    /// 主オラクル1（最小性・最強, 観点1）: 中心食地点で返った time_tt について、独立再構成 m²(t) が
    /// 両側 δ∈{1s,10s,60s} で増加する（局地最小）。solver の内部手法に依存しない外部縛り。
    #[test]
    fn central_maximum_is_a_local_minimum_of_independent_m2() {
        let dt = EspenakMeeusDeltaT;
        let src = make_source(&dt);
        let (obs, lon) = central_observer();
        let max = solve_local_maximum(
            &src,
            &obs,
            lon,
            window_around(tt_2017_max(), 3.0),
            config_tight(),
        )
        .expect("central site should yield a local maximum");

        let m2_center = m2_at(&src, &obs, lon, max.time_tt);
        // 同時刻判定の数値揺らぎを避けるため、谷底でも明確に増加するスケールを選ぶ。
        // 1s/10s/60s いずれも局地最小なら m²(t±δ) ≥ m²(t)（厳密谷では >）。
        for &delta in &[1.0_f64, 10.0, 60.0] {
            let plus = m2_shift(&src, &obs, lon, max.time_tt, delta);
            let minus = m2_shift(&src, &obs, lon, max.time_tt, -delta);
            assert!(
                plus >= m2_center,
                "m²(t_max+{delta}s)={plus} must be ≥ m²(t_max)={m2_center} (local min)"
            );
            assert!(
                minus >= m2_center,
                "m²(t_max−{delta}s)={minus} must be ≥ m²(t_max)={m2_center} (local min)"
            );
        }
        // 谷であること（平底でない・少なくとも一方の 60s で狭義増加）を確認し、最大値（最小食）取り違えを弾く。
        let p60 = m2_shift(&src, &obs, lon, max.time_tt, 60.0);
        let m60 = m2_shift(&src, &obs, lon, max.time_tt, -60.0);
        assert!(
            p60 > m2_center && m60 > m2_center,
            "60s away on both sides must strictly increase (rejects flat/maximum mix-up): \
             +60s={p60}, −60s={m60}, center={m2_center}"
        );
    }

    /// 主オラクル2（min_separation 整合, 観点2）: 返った min_separation² が独立再計算 m²(t_max) に一致。
    /// min_separation が「その時刻の中心間距離」であることを直接縛る（単位 Re）。
    #[test]
    fn min_separation_squared_matches_independent_m2() {
        let dt = EspenakMeeusDeltaT;
        let src = make_source(&dt);
        let (obs, lon) = central_observer();
        let max = solve_local_maximum(
            &src,
            &obs,
            lon,
            window_around(tt_2017_max(), 3.0),
            config_tight(),
        )
        .expect("central site should yield a local maximum");

        let m2_oracle = m2_at(&src, &obs, lon, max.time_tt);
        let m2_returned = max.min_separation * max.min_separation;
        // min_separation は非負（距離）。
        assert!(
            max.min_separation >= 0.0,
            "min_separation must be non-negative, got {}",
            max.min_separation
        );
        // 1e-9 Re²。最大食付近で |d(m²)/dt| は小さい（谷底）ので、x_tolerance による
        // 時刻差は m² に二次でしか効かず、整合は非常に厳しく取れる。
        assert!(
            (m2_returned - m2_oracle).abs() < 1e-9,
            "min_separation²={m2_returned} must match independent m²(t_max)={m2_oracle} (tol 1e-9 Re²)"
        );
    }

    /// 主オラクル3（d(m²)/dt≈0 の時刻誤差換算, 観点3）: 返った time_tt での中心差分 d(m²)/dt を
    /// 二階微分 d²(m²)/dt² で割った「最小からの等価時刻ずれ」が <0.5s（目標 ±1〜2s の 1/10 級）。
    /// Newton 1 ステップ Δt = −(d/dt)/(d²/dt²) を最小からのずれの推定として用いる
    /// （ISSUE-025 の残差→時刻誤差換算と同方針。ここでは導関数ゼロ点なので二階微分で割る）。
    #[test]
    fn central_derivative_zero_implies_subsecond_time_error() {
        let dt = EspenakMeeusDeltaT;
        let src = make_source(&dt);
        let (obs, lon) = central_observer();
        let max = solve_local_maximum(
            &src,
            &obs,
            lon,
            window_around(tt_2017_max(), 3.0),
            config_tight(),
        )
        .expect("central site should yield a local maximum");

        // 中心差分を 2 つのステップ幅 h∈{1s,4s} で評価し、どちらでも等価時刻誤差 <0.5s を確認する。
        // 単一スケールだと谷底の打切り誤差で偶然通る恐れがあるため、刻みを変えても結論が不変なことを縛る
        // （R2）。単位 Re²/s, Re²/s²。
        let f0 = m2_at(&src, &obs, lon, max.time_tt);
        for &h in &[1.0_f64, 4.0] {
            let fp = m2_shift(&src, &obs, lon, max.time_tt, h);
            let fm = m2_shift(&src, &obs, lon, max.time_tt, -h);
            let first = (fp - fm) / (2.0 * h); // d(m²)/dt
            let second = (fp - 2.0 * f0 + fm) / (h * h); // d²(m²)/dt²
            assert!(
                second > 0.0,
                "d²(m²)/dt² must be positive at a minimum (h={h}s, got {second}); not a minimum"
            );
            let time_error_s = (first / second).abs(); // Newton 1 step → 最小からのずれ[s]
            assert!(
                time_error_s < 0.5,
                "equivalent time error from minimum {time_error_s} s exceeds 0.5 s \
                 (h={h}s, d/dt={first}, d²/dt²={second})"
            );
        }
    }

    /// 観点4（中心線上）: 中心食地点では最大食の m_min が小さい（≈0）。
    /// 影軸が観測者をほぼ通過する中心線なので最小中心間距離は本影半径以下。
    /// 格子・地点近似で <0.05 Re（≈320km）に緩める（projection.rs の外部ピンと同基準）。
    #[test]
    fn central_site_min_separation_is_near_zero() {
        let dt = EspenakMeeusDeltaT;
        let src = make_source(&dt);
        let (obs, lon) = central_observer();
        let max = solve_local_maximum(
            &src,
            &obs,
            lon,
            window_around(tt_2017_max(), 3.0),
            config_tight(),
        )
        .expect("central site should yield a local maximum");
        assert!(
            max.min_separation < 0.05,
            "central (centerline) site min_separation should be ≈0, got {} Re",
            max.min_separation
        );
    }

    /// 観点5（接触の内側）: 中心食地点で最大食は C1 と C4 の間に入る（c1 < max < c4）。
    /// 接触は ISSUE-025（独立に検証済の別 solver）から得る。最大食が部分食区間の内側にある
    /// という幾何制約を別経路で縛る。
    #[test]
    fn central_maximum_is_between_c1_and_c4() {
        let dt = EspenakMeeusDeltaT;
        let src = make_source(&dt);
        let (obs, lon) = central_observer();
        let w = window_around(tt_2017_max(), 3.0);

        let max = solve_local_maximum(&src, &obs, lon, w, config_tight())
            .expect("central site should yield a local maximum");
        let contacts = solve_local_contacts(&src, &obs, lon, w, config_tight())
            .expect("central site should yield a contact set");

        let c1 = jd_of(contacts.c1.expect("central site must have C1").time_tt);
        let c4 = jd_of(contacts.c4.expect("central site must have C4").time_tt);
        let tmax = jd_of(max.time_tt);
        assert!(
            c1 < tmax && tmax < c4,
            "maximum must lie inside partial phase: c1={c1} < max={tmax} < c4={c4}"
        );
    }

    /// 観点6（窓幅不変）: 探索窓を ±3h / ±2h / 非対称（−1h〜+3h）で解き、time_tt が 2 秒以内で一致。
    /// 窓幅・窓位相が変わると実装内部の粗ブラケット分割（分割位相）が変わるため、刻み感度・位相依存の
    /// 偽陰性を実効的に縛る（R5）。いずれの窓も最小を内部に含む（半幅/前後幅 > 部分食半継続 ~80 分）。
    #[test]
    fn maximum_time_is_window_width_invariant() {
        let dt = EspenakMeeusDeltaT;
        let src = make_source(&dt);
        let (obs, lon) = central_observer();

        let solve_window = |w: TimeInterval<TtInstant>| {
            solve_local_maximum(&src, &obs, lon, w, config_tight())
                .expect("central site should yield a local maximum")
        };

        let wide = solve_window(window_around(tt_2017_max(), 3.0)); // ±3h
        let narrow = solve_window(window_around(tt_2017_max(), 2.0)); // ±2h（半幅 120 分 > 部分食半継続）
                                                                      // 非対称窓（最大食 −1h 〜 +3h）。window_around は対称なので TimeInterval を直接構築する。
        let center = jd_of(tt_2017_max());
        let asymmetric = solve_window(TimeInterval {
            start: tt_jd(center - 1.0 / 24.0),
            end: tt_jd(center + 3.0 / 24.0),
        });

        let diff_narrow = (jd_of(wide.time_tt) - jd_of(narrow.time_tt)).abs() * SECONDS_PER_DAY;
        assert!(
            diff_narrow < 2.0,
            "max time differs between ±3h and ±2h windows by {diff_narrow} s (tol 2 s); \
             coarse-bracket step sensitivity",
        );
        let diff_asym = (jd_of(wide.time_tt) - jd_of(asymmetric.time_tt)).abs() * SECONDS_PER_DAY;
        assert!(
            diff_asym < 2.0,
            "max time differs between ±3h and asymmetric (−1h..+3h) windows by {diff_asym} s \
             (tol 2 s); coarse-bracket split-phase dependence",
        );
    }

    /// 観点7（TT/UTC 整合）: time_utc == tt_to_utc(time_tt)（2017 は post-1972 で変換可能）。
    #[test]
    fn maximum_utc_matches_tt_to_utc() {
        let dt = EspenakMeeusDeltaT;
        let src = make_source(&dt);
        let (obs, lon) = central_observer();
        let max = solve_local_maximum(
            &src,
            &obs,
            lon,
            window_around(tt_2017_max(), 3.0),
            config_tight(),
        )
        .expect("central site should yield a local maximum");
        let want = umbra_core::time::tt_to_utc(max.time_tt)
            .expect("2017 is post-1972, tt_to_utc must succeed");
        assert_eq!(max.time_utc, want, "time_utc must equal tt_to_utc(time_tt)");
    }

    /// 観点1 系（部分食地点でも最大食は存在）: 高緯度の部分食のみ地点でも最大食（幾何最小）が
    /// 返り（非 Option）、独立再構成 m²(t) の局地最小であること（両側 60s で増加）。
    /// 部分食地点では m_min は本影半径より大きい（≈0 でない）が、最小性は同じく成立する。
    #[test]
    fn partial_site_still_has_a_local_minimum() {
        let dt = EspenakMeeusDeltaT;
        let src = make_source(&dt);
        let (obs, lon) = partial_observer();
        let max = solve_local_maximum(
            &src,
            &obs,
            lon,
            window_around(tt_2017_max(), 3.0),
            config_tight(),
        )
        .expect("partial site must still have a (geometric) maximum eclipse (non-Option)");

        let center = m2_at(&src, &obs, lon, max.time_tt);
        let p = m2_shift(&src, &obs, lon, max.time_tt, 60.0);
        let m = m2_shift(&src, &obs, lon, max.time_tt, -60.0);
        assert!(
            p > center && m > center,
            "partial site maximum must be a strict local min of m²: +60s={p}, −60s={m}, center={center}"
        );
        // 部分食地点なので min_separation は中心線より大きい（本影を外す）。0 ではない健全性。
        assert!(
            max.min_separation > 0.05,
            "partial (high-lat) site min_separation should be well above 0, got {} Re",
            max.min_separation
        );
    }

    /// 観点1 系・前提独立確認（R4）: `partial_observer`(60N,−100°) が本当に「部分食地点」であること
    /// を、独立検証済の別 solver `solve_local_contacts`（ISSUE-025）で確認する。本影を外す地点では
    /// 内接 C2/C3 が存在せず（None）、外接 C1/C4 のみ存在する（Some）。これにより
    /// 「部分食地点でも最大食（幾何最小）は非 Option で存在する」という分類の前提が独立に立つ。
    #[test]
    fn partial_site_has_no_central_contacts_only_partial() {
        let dt = EspenakMeeusDeltaT;
        let src = make_source(&dt);
        let (obs, lon) = partial_observer();
        let w = window_around(tt_2017_max(), 3.0);

        let contacts = solve_local_contacts(&src, &obs, lon, w, config_tight())
            .expect("partial site should yield a contact set");

        assert!(
            contacts.c2.is_none() && contacts.c3.is_none(),
            "partial site must have no central (inscribed) contacts: c2={:?}, c3={:?}",
            contacts.c2,
            contacts.c3
        );
        assert!(
            contacts.c1.is_some() && contacts.c4.is_some(),
            "partial site must still have outer contacts C1/C4: c1={:?}, c4={:?}",
            contacts.c1,
            contacts.c4
        );
    }

    /// 観点8（異常系・単調窓）: 最小を含まない単調な窓では Solver(RootNotBracketed) を返す。
    /// 最大食より十分後（部分食終了 C4 の後）の窓を取ると m は単調増加（月影が遠ざかる）。
    /// このとき最小は窓端で内部極小がない＝ d(m²)/dt=0 をブラケットできない（機構的ブラケット失敗）。
    /// DegenerateGeometry は影錐の幾何退化専用で兄弟 solver と不整合のため使わない（契約 M2）。
    /// テスト側オラクルで「窓内 m² が単調増加・最小は左端」を独立確認してから縛る
    /// （窓前提が成り立つことの保証）。
    #[test]
    fn monotone_window_without_interior_minimum_is_root_not_bracketed() {
        let dt = EspenakMeeusDeltaT;
        let src = make_source(&dt);
        let (obs, lon) = central_observer();

        // 最大食の +2h を起点に +1h 後までの窓（部分食はとうに終了し m は単調増加）。
        let start = tt_jd(jd_of(tt_2017_max()) + 2.0 / 24.0);
        let end = tt_jd(jd_of(tt_2017_max()) + 3.0 / 24.0);
        let w = TimeInterval { start, end };

        // 前提の独立確認: 窓全域で m² が単調増加（隣接サンプルが常に増える）→ 内部に極小なし。
        let n = 400;
        let (lo, hi) = (jd_of(w.start), jd_of(w.end));
        let mut prev = m2_at(&src, &obs, lon, tt_jd(lo));
        let mut monotone_increasing = true;
        for k in 1..=n {
            let jd = lo + (hi - lo) * (f64::from(k) / f64::from(n));
            let cur = m2_at(&src, &obs, lon, tt_jd(jd));
            if cur <= prev {
                monotone_increasing = false;
                break;
            }
            prev = cur;
        }
        assert!(
            monotone_increasing,
            "precondition: chosen window must have monotonically increasing m² (no interior minimum)"
        );

        let result = solve_local_maximum(&src, &obs, lon, w, config_tight());
        assert_eq!(
            result,
            Err(EclipseError::Solver(SolverError::RootNotBracketed)),
            "monotone window (minimum at edge, none interior) must return Solver(RootNotBracketed), got {result:?}"
        );
    }

    /// 観点8b（平底/定数 m² 退化・graceful）: 時刻不変供給源では m²(t) が定数になり、窓内に
    /// 内部極小がない（粗走査の最小が窓端／d(m²)/dt をブラケットできない）。このとき実装は
    /// **panic せず** `Err(Solver(RootNotBracketed))` を返す（単調窓と同経路の堅牢性）。
    ///
    /// 契約 M1 の再解釈: spec は「皆既帯で m が |L2| を下回る平底（m がほぼ一定）」の MockEphemeris
    /// 受入を要求するが、(i) 皆既中も m は単一極小を持ち平坦なのは食分(obscuration)であって m ではない、
    /// (ii) MockEphemeris は静的で時変 m²(t) を作れない。よって平底＝「m² が(ほぼ)定数になる退化」と
    /// 再解釈し、定数 m² でも graceful に RootNotBracketed を返す堅牢性で縛る（最小が一意でない平底でも
    /// 破綻しないという spec の核心）。
    #[test]
    fn constant_m2_flat_bottom_is_graceful_root_not_bracketed() {
        let src = ConstantSource {
            elems: constant_elems(),
        };
        let (obs, lon) = central_observer();

        // 前提の独立確認: この供給源では m²(t) が時刻によらず一定（定数）→ 内部極小なし。
        // 窓内の複数サンプルで m² が同一値であることを確認してから縛る。
        let w = window_around(tt_2017_max(), 3.0);
        let (lo, hi) = (jd_of(w.start), jd_of(w.end));
        let base = m2_at(&src, &obs, lon, tt_jd(lo));
        for k in 1..=8 {
            let jd = lo + (hi - lo) * (f64::from(k) / 8.0);
            let v = m2_at(&src, &obs, lon, tt_jd(jd));
            assert_eq!(
                v, base,
                "precondition: constant source must yield time-invariant m² (flat bottom), got {v} != {base}"
            );
        }

        // 定数 m²（平底/退化）でも実装は panic せず graceful にブラケット失敗を返す。
        let result = solve_local_maximum(&src, &obs, lon, w, config_tight());
        assert_eq!(
            result,
            Err(EclipseError::Solver(SolverError::RootNotBracketed)),
            "constant m² (flat bottom / degenerate) must gracefully return Solver(RootNotBracketed), got {result:?}"
        );
    }

    /// 観点9（実日食 ballpark, 第二義・緩め）: 2017-08-21 中心食地点(37.5N)の最大食が
    ///  - 当日 2017-08-21 にある（TT-JD 2457986.5〜2457987.6）
    ///  - UTC 18 時台（北米中緯度の最大食は 18:xx UTC）
    ///  - 部分食の中央付近（(max−c1) と (c4−max) がどちらも部分食半継続オーダー）
    ///
    /// flaky なハードコード秒値は使わず、桁・継続スケールのみで縛る。
    /// 出典: 一般に知られる 2017-08-21 皆既日食の最大食時刻スケール（絶対基準にしない）。
    #[test]
    fn central_site_2017_maximum_in_ballpark() {
        let dt = EspenakMeeusDeltaT;
        let src = make_source(&dt);
        let (obs, lon) = central_observer();
        let w = window_around(tt_2017_max(), 3.0);
        let max = solve_local_maximum(&src, &obs, lon, w, config_tight())
            .expect("central site should yield a local maximum");
        let contacts =
            solve_local_contacts(&src, &obs, lon, w, config_tight()).expect("contact set");

        // 当日（TT-JD）。
        let tmax = jd_of(max.time_tt);
        assert!(
            (2_457_986.5..2_457_987.6).contains(&tmax),
            "max TT-JD {tmax} not on 2017-08-21"
        );

        // UTC 18 時台（グレゴリオ暦へ変換して hour+minute を取り出す）。
        let (_y, _mo, _d, hh, mm, _ss) = max.time_utc.to_gregorian();
        let hour_utc = f64::from(hh) + f64::from(mm) / 60.0;
        assert!(
            (17.5..=19.0).contains(&hour_utc),
            "max UTC hour {hour_utc} not ~18:xx (north-american mid-latitude greatest eclipse)"
        );

        // 部分食の中央付近: max が C1 と C4 の間で、両側マージンが部分食半継続の妥当割合。
        let c1 = jd_of(contacts.c1.expect("C1").time_tt);
        let c4 = jd_of(contacts.c4.expect("C4").time_tt);
        let before = (tmax - c1) / (c4 - c1);
        assert!(
            (0.3..=0.7).contains(&before),
            "maximum should sit near the middle of the partial phase, fraction-from-c1={before}"
        );
    }
}
