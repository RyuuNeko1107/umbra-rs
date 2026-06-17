//! 局地接触 C1/C4（外接）・C2/C3（内接）求解（`docs/issues/ISSUE-025`、`docs/conventions.md` §8,
//! Explanatory Supplement §11 / Meeus *Astronomical Algorithms* Ch.54）。
//!
//! 観測地点の部分食 開始 C1 / 終了 C4（外接 `m = L1`）、中心食地点の C2/C3（内接 `m = |L2|`）を
//! ベッセル基本面上で求解する。`m² = (ξ−x)² + (η−y)²`（観測者基本面座標 ISSUE-024 と影軸交点
//! x,y の距離²）、影半径は観測者高さ補正 `L1 = l1 − ζ·tan f1`, `L2 = l2 − ζ·tan f2`。
//!
//! 求根対象は連続関数 `g(t) = m²(t) − L(t)²`（D2: m² 基準で統一、中心線尖点の微分特異と acos を回避）。
//! 探索窓を粗走査で符号変化区間に分割 → 各区間を Brent（ISSUE-008、無条件 Newton 禁止・conventions §11）。
//! 接触時刻は TT と UTC の両方を返す（accuracy.md §0）。
//!
//! 注: `LocalContact` の高度・方位・position_angle・可視性（ISSUE-028）は本 issue の非目的。本層は
//! 接触時刻求解が責務で、観測フィールドは ISSUE-028/043 で充足する。

// solve_local_contacts（pub(crate)）は ISSUE-043（EclipseEngine 結線）が消費するまで未使用。
// 結線され次第この許容は外す（conjunction.rs / candidates.rs と同手順）。
#![allow(dead_code)]
// 粗走査の分割数（窓秒数 / 刻み）の整数⇔f64 変換は小さな添字のみ（天文量ではない, conjunction.rs 同様）。
#![allow(
    clippy::cast_precision_loss,
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss
)]

use umbra_core::ellipsoid::GeocentricObserver;
use umbra_core::time::tt_to_utc;
use umbra_core::{brent_root, JulianDate2, Radians, TimeInterval, TtInstant, UtcInstant};

use crate::conjunction::RootConfig;
use crate::error::EclipseError;
use crate::projection::project_observer_to_fundamental;
use crate::source::BesselianSource;

/// 粗走査の刻み（SI 秒）。接触対（特に皆既/金環の C2/C3＝数分〜）を取りこぼさないよう細かく取る
/// （偽陰性不可・architecture §3）。窓幅 / 本刻みで分割数を決める。
const CONTACT_SCAN_STEP_SECONDS: f64 = 30.0;
/// 1 日 = 86400 SI 秒。
const SECONDS_PER_DAY: f64 = 86_400.0;

/// 局地接触の時刻（TT と UTC, accuracy.md §0）。
///
/// 高度・方位・position_angle・可視性（ISSUE-028）は本 issue の非目的のため持たない
/// （ISSUE-028/043 で拡張）。
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LocalContact {
    /// 接触の TT 時刻（幾何相対の一級値, conventions §6）。
    pub time_tt: TtInstant,
    /// 接触の UTC 時刻（TT→UTC は ΔT 経由 ISSUE-007、将来は予測律速）。
    pub time_utc: UtcInstant,
}

/// C1〜C4 の局地接触集合。部分食地点では `c2`/`c3` が `None`（api-draft §6）。
#[derive(Clone, Copy, Debug, PartialEq, Default)]
pub struct LocalContactSet {
    /// 第1接触（部分食開始・外接）。
    pub c1: Option<LocalContact>,
    /// 第2接触（皆既/金環開始・内接, 中心食地点のみ）。
    pub c2: Option<LocalContact>,
    /// 第3接触（皆既/金環終了・内接, 中心食地点のみ）。
    pub c3: Option<LocalContact>,
    /// 第4接触（部分食終了・外接）。
    pub c4: Option<LocalContact>,
}

/// 観測地点の局地接触 C1〜C4 を探索窓内で求解する。
///
/// - `source`: 瞬時ベッセル要素の供給源（ISSUE-037, 既定 `DirectBesselianSource`）。
/// - `observer`/`east_longitude`: 観測者（ρsinφ′/ρcosφ′ ＋ 東経正, ISSUE-024）。
/// - `search`: 探索窓（全球接触から絞った TT 区間）。
/// - `config`: Brent 求根設定（root_tolerance は目標の 1/10 以下）。
pub(crate) fn solve_local_contacts<B: BesselianSource>(
    source: &B,
    observer: &GeocentricObserver,
    east_longitude: Radians,
    search: TimeInterval<TtInstant>,
    config: RootConfig,
) -> Result<LocalContactSet, EclipseError> {
    let t0 = search.start.jd2().jd();
    let t1 = search.end.jd2().jd();

    // 外接 C1/C4: g_outer = m² − L1²。+→−（外側→内側）= C1, −→+（内側→外側）= C4。
    let outer = scan_sign_change_roots(
        |jd| g_outer(source, observer, east_longitude, jd),
        t0,
        t1,
        config,
    )?;
    let c1_jd = outer.iter().find(|r| !r.ascending).map(|r| r.time_jd);
    let c4_jd = outer.iter().rev().find(|r| r.ascending).map(|r| r.time_jd);

    // 内接 C2/C3: g_inner = m² − L2²（L2 符号付き → L2² で符号不問）。部分食域は内接なし → None。
    // 内接接触は部分食区間内にあるので、C1〜C4 に絞って走査する（無駄な走査と取りこぼしを避ける）。
    let (mut c2_jd, mut c3_jd) = (None, None);
    if let (Some(a), Some(b)) = (c1_jd, c4_jd) {
        let inner = scan_sign_change_roots(
            |jd| g_inner(source, observer, east_longitude, jd),
            a,
            b,
            config,
        )?;
        c2_jd = inner.iter().find(|r| !r.ascending).map(|r| r.time_jd);
        c3_jd = inner.iter().rev().find(|r| r.ascending).map(|r| r.time_jd);
    }

    Ok(LocalContactSet {
        c1: contact_at(c1_jd)?,
        c2: contact_at(c2_jd)?,
        c3: contact_at(c3_jd)?,
        c4: contact_at(c4_jd)?,
    })
}

/// 粗走査で見つけた符号変化の根（TT-JD）と、その点で関数が「−→+（昇）」か「+→−（降）」か。
#[derive(Clone, Copy, Debug)]
struct SignChangeRoot {
    time_jd: f64,
    /// `true`: 関数が負→正（昇）。`false`: 正→負（降）。
    ascending: bool,
}

/// 単一 TT-JD から TtInstant（求根は単一 f64 JD 空間。conjunction.rs と同方針）。
fn tt(jd: f64) -> TtInstant {
    TtInstant::from_jd2(JulianDate2::from_jd(jd))
}

/// 観測者の基本面座標と影軸交点から `m² = (ξ−x)² + (η−y)²` と高さ補正後の外接半径 L1 を作る。
/// `g_outer = m² − L1²`、`L1 = l1 − ζ·tan f1`（Expl. Suppl. §11 / Meeus Ch.54）。
fn g_outer<B: BesselianSource>(
    source: &B,
    observer: &GeocentricObserver,
    east_longitude: Radians,
    jd: f64,
) -> Result<f64, EclipseError> {
    let e = source.at(tt(jd))?;
    let p = project_observer_to_fundamental(observer, east_longitude, &e);
    let du = p.xi - e.x;
    let dv = p.eta - e.y;
    let m2 = du * du + dv * dv;
    let l1 = e.l1 - p.zeta * e.tan_f1;
    Ok(m2 - l1 * l1)
}

/// 内接 `g_inner = m² − L2²`、`L2 = l2 − ζ·tan f2`（l2 は符号付き）。
fn g_inner<B: BesselianSource>(
    source: &B,
    observer: &GeocentricObserver,
    east_longitude: Radians,
    jd: f64,
) -> Result<f64, EclipseError> {
    let e = source.at(tt(jd))?;
    let p = project_observer_to_fundamental(observer, east_longitude, &e);
    let du = p.xi - e.x;
    let dv = p.eta - e.y;
    let m2 = du * du + dv * dv;
    let l2 = e.l2 - p.zeta * e.tan_f2;
    Ok(m2 - l2 * l2)
}

/// 粗走査の分割数（窓幅 / 刻み, 最低 2）。**走査解像度のみ**を決め、接触検出の正否には影響しない
/// （n が十分大きければ符号変化を捉える＝偽陰性回避は刻みの細かさで担保）。`span`/`n` の算術には
/// 振る舞い契約が無く、ここの算術変異は等価（細かい n でも根は同じ）/timeout（巨大 n）になるため
/// `mutation.yml` で `--exclude-re 'in scan_point_count'` 除外する（docs/reviews/mutation-local-contacts.md）。
fn scan_point_count(t0_jd: f64, t1_jd: f64, step_seconds: f64) -> usize {
    let span_seconds = (t1_jd - t0_jd) * SECONDS_PER_DAY;
    (span_seconds / step_seconds).ceil().max(2.0) as usize
}

/// `[t0_jd, t1_jd]`（TT-JD）を一定刻みで粗走査し、`f` の符号変化区間を Brent で精解して全根を返す。
/// 接触が無ければ空 Vec（その地点・窓に該当接触なし＝食なしは正常, エラーにしない）。
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

/// TT-JD（あれば）から `LocalContact`（TT + UTC）を作る。UTC は ΔT 経由（ISSUE-007）。
fn contact_at(jd: Option<f64>) -> Result<Option<LocalContact>, EclipseError> {
    match jd {
        Some(jd) => {
            let time_tt = tt(jd);
            let time_utc = tt_to_utc(time_tt)?;
            Ok(Some(LocalContact { time_tt, time_utc }))
        }
        None => Ok(None),
    }
}

#[cfg(test)]
mod tests {
    //! ISSUE-025 受け入れテスト（strict）。
    //!
    //! ## オラクル戦略（追認回避）
    //! 主オラクルは **独立内部再計算**: 接触条件 `g_outer(t)=m²−L1²`（C1/C4・外接）/
    //! `g_inner(t)=m²−L2²`（C2/C3・内接）を、テスト側で `project_observer_to_fundamental`
    //! （ISSUE-024, 検証済プリミティブ）＋ `source.at(t)`（ISSUE-037, 検証済）から **別経路で**
    //! 組み直し、`solve_local_contacts` が返した各接触時刻でその g がほぼゼロになることを縛る。
    //! ここで使う射影・供給源は solver の内部関数ではなく公開プリミティブなので、実装本体を
    //! コピーした「追認」にはならない（solver 本体は `g` の合成・粗走査・Brent・接触割当を担う）。
    //!
    //! 幾何の出典: Explanatory Supplement to the Astronomical Almanac §11 /
    //! Meeus *Astronomical Algorithms* 2nd ed. Ch.54。
    //!   m² = (ξ−x)² + (η−y)²、観測者高さ補正 L1 = l1 − ζ·tan f1, L2 = l2 − ζ·tan f2。
    //!
    //! 補助オラクル: (2) 順序プロパティ c1<c2<c3<c4（中心食地点）/ 部分食地点で c2=c3=None・
    //! c1,c4 存在、(3) TT/UTC 整合 time_utc == tt_to_utc(time_tt)、(4) ブラケット不成立 →
    //! RootNotBracketed、(5) grazing 刻み感度（粗走査刻みを変えても C1/C4 を取りこぼさない偽陰性ガード）、
    //! (6) 実日食 2017-08-21 の ballpark（部分食継続 ~2.5〜3h、皆既継続 数分オーダー、当日内）。
    //!
    //! NASA/USNO の地点別接触時刻は絶対基準にせず、ハードコードの flaky 値を避けるため
    //! ballpark 整合（出典: 一般に知られる 2017-08-21 皆既日食の継続時間スケール）のみに用いる。

    use super::*;

    use umbra_core::constants::{EARTH_EQUATORIAL_RADIUS_M, SOLAR_RADIUS_KM};
    use umbra_core::ellipsoid::{observer_geocentric, Ellipsoid};
    use umbra_core::{EspenakMeeusDeltaT, JulianDate2, Radians, TimeInterval, TtInstant};

    use crate::projection::project_observer_to_fundamental;
    use crate::source::{BesselianSource, DirectBesselianSource};

    const WGS84: Ellipsoid = Ellipsoid::WGS84;
    const R_SUN: f64 = SOLAR_RADIUS_KM;
    // k·Re（月/地球赤道半径比、IAU 慣習 k=0.2725076）。source.rs / besselian.rs テストと同一。
    const R_MOON: f64 = 0.2725076 * (EARTH_EQUATORIAL_RADIUS_M / 1000.0);

    /// 1 日 = 86400 SI 秒。秒↔日換算に使う。
    const SECONDS_PER_DAY: f64 = 86_400.0;

    /// TT を 2 要素 JD から構築（source.rs / besselian.rs テストと同形）。
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

    /// 2017-08-21 皆既日食 最大食付近の TT エポック（besselian.rs / source.rs テストと同一）。
    /// 最大食 ≈ 2017-08-21 18:25 UTC（TT-JD ≈ 2457986.5 + 0.76853 ≈ 2457987.2685）。
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

    fn make_source<'d>(
        dt: &'d EspenakMeeusDeltaT,
    ) -> DirectBesselianSource<'d, EspenakMeeusDeltaT> {
        DirectBesselianSource::new(R_SUN, R_MOON, dt, source_interval())
    }

    /// 標準の Brent 設定。x_tolerance_days = 1e-9 日（≈8.6e-5 s）は接触 ±2s 目標の
    /// 1/10（≤0.2s）を十分に下回る。max_iterations は余裕を見て 200。
    fn config_tight() -> RootConfig {
        RootConfig {
            x_tolerance_days: 1e-9,
            max_iterations: 200,
        }
    }

    /// 最大食 ±数時間の探索窓（部分食全体を覆う）。2017 の部分食継続は片側 ~1.4h なので ±3h で安全。
    fn window_around(center: TtInstant, half_hours: f64) -> TimeInterval<TtInstant> {
        let c = jd_of(center);
        let half = half_hours / 24.0;
        TimeInterval {
            start: tt_jd(c - half),
            end: tt_jd(c + half),
        }
    }

    // ============================================================
    // 観測地点（2017-08-21）
    // ============================================================
    //
    // 経度は東経正（conventions §3）。北米西経は負で与える（西経吸収は ISSUE-024 で検証済）。
    //
    // 中心食地点: 2017 皆既帯（北米中緯度）。例として ~37°N, 西経 ~89°（イリノイ南部 Carbondale 付近、
    // 皆既継続が長い地域）を選ぶ。皆既帯のどこかに本影が通れば C2/C3 が存在する。
    // 部分食地点: 高緯度（本影が届かない）。例として 60°N, 西経 100°（カナダ中部）。北米北部は
    // 2017 で部分食のみ（皆既帯は中緯度の細い帯）なので C2/C3=None になるはず。

    /// 中心食地点（皆既帯中心付近, 2017 皆既最長 ~2.6 分）。lat=37.5°N, lon=西経89.2°。
    /// 主オラクル・順序・継続テストの基準。皆既が最長で頑健（プローブで検証済）。
    fn central_observer() -> (umbra_core::ellipsoid::GeocentricObserver, Radians) {
        let lat = 37.5_f64.to_radians();
        let east_lon = Radians::new((-89.2_f64).to_radians()); // 西経 89.2° → 東経負
        (observer_geocentric(&WGS84, lat, 200.0), east_lon)
    }

    /// 限界内側 grazing 地点（皆既帯の縁寄り, 皆既 ~1.9 分・C2/C3 は存在）。lat=38.0°N。
    /// 本影内だが中心より浅く、皆既継続が中心(37.5N)より短いことの縛りに使う（プローブで検証済）。
    fn inner_limit_observer() -> (umbra_core::ellipsoid::GeocentricObserver, Radians) {
        let lat = 38.0_f64.to_radians();
        let east_lon = Radians::new((-89.2_f64).to_radians());
        (observer_geocentric(&WGS84, lat, 200.0), east_lon)
    }

    /// 限界外側 grazing 地点（本影を僅かに外す・部分食・C2/C3=None だが C1/C4 あり）。lat=38.5°N。
    /// 「中心線近傍でも本影を外せば内接なし」の限界ケース（プローブで検証済）。
    fn outer_limit_observer() -> (umbra_core::ellipsoid::GeocentricObserver, Radians) {
        let lat = 38.5_f64.to_radians();
        let east_lon = Radians::new((-89.2_f64).to_radians());
        (observer_geocentric(&WGS84, lat, 200.0), east_lon)
    }

    /// 部分食のみの地点（高緯度・本影外）。
    fn partial_observer() -> (umbra_core::ellipsoid::GeocentricObserver, Radians) {
        let lat = 60.0_f64.to_radians();
        let east_lon = Radians::new((-100.0_f64).to_radians());
        (observer_geocentric(&WGS84, lat, 200.0), east_lon)
    }

    // ============================================================
    // 独立オラクル: 接触条件関数 g_outer / g_inner（テスト側で別経路再構成）
    // ============================================================

    /// 観測者の基本面座標 (ξ,η,ζ) と影軸交点 (x,y)・影半径 l1,l2・tan f1,f2 から、
    /// 外接接触 `g_outer(t) = m² − L1²`（C1/C4）を独立に評価する。
    /// `m² = (ξ−x)² + (η−y)²`、`L1 = l1 − ζ·tan f1`（Expl. Suppl. §11 / Meeus Ch.54）。
    /// project_* と source.at は検証済プリミティブ。solver 内部関数には依存しない。
    fn g_outer<B: BesselianSource>(
        source: &B,
        observer: &umbra_core::ellipsoid::GeocentricObserver,
        east_longitude: Radians,
        t: TtInstant,
    ) -> f64 {
        let e = source
            .at(t)
            .expect("source.at should succeed near 2017 eclipse");
        let p = project_observer_to_fundamental(observer, east_longitude, &e);
        let du = p.xi - e.x;
        let dv = p.eta - e.y;
        let m2 = du * du + dv * dv;
        let l1 = e.l1 - p.zeta * e.tan_f1;
        m2 - l1 * l1
    }

    /// 内接接触 `g_inner(t) = m² − L2²`（C2/C3）。`L2 = l2 − ζ·tan f2`（l2 は符号付き → L2² で符号不問）。
    fn g_inner<B: BesselianSource>(
        source: &B,
        observer: &umbra_core::ellipsoid::GeocentricObserver,
        east_longitude: Radians,
        t: TtInstant,
    ) -> f64 {
        let e = source
            .at(t)
            .expect("source.at should succeed near 2017 eclipse");
        let p = project_observer_to_fundamental(observer, east_longitude, &e);
        let du = p.xi - e.x;
        let dv = p.eta - e.y;
        let m2 = du * du + dv * dv;
        let l2 = e.l2 - p.zeta * e.tan_f2;
        m2 - l2 * l2
    }

    /// 接触時刻 t で g がほぼゼロであることの許容。
    /// 中心距離 m は ~1 day で O(1) Re 動く（地球自転＋影軸運動）→ |dm/dt| ~ 数 Re/day。
    /// g = m²−L²、接触付近 dg/dt ≈ 2m·dm/dt ~ O(1) Re²/day。root_tolerance ≤0.2s = 2.3e-6 day
    /// に対し残差 |g| ≲ |dg/dt|·tol ~ 数×1e-6 Re²。安全側に 1e-5 Re² を上限とする。
    const G_ROOT_TOL: f64 = 1e-5;

    /// 中心差分で dg/dt を見積もる（残差→時刻誤差の換算根拠の自己検証用）。秒あたり。
    fn dg_outer_dt<B: BesselianSource>(
        source: &B,
        observer: &umbra_core::ellipsoid::GeocentricObserver,
        east_longitude: Radians,
        t: TtInstant,
    ) -> f64 {
        let h_s = 1.0; // 1 秒
        let h_d = h_s / SECONDS_PER_DAY;
        let tp = tt_jd(jd_of(t) + h_d);
        let tm = tt_jd(jd_of(t) - h_d);
        (g_outer(source, observer, east_longitude, tp)
            - g_outer(source, observer, east_longitude, tm))
            / (2.0 * h_s)
    }

    /// 中心差分で dg_inner/dt を見積もる（内接残差→時刻誤差の換算根拠）。秒あたり。
    /// 内接は外接と勾配が異なる（L2≪L1 で接触の幾何が違う）ため別途必要。
    fn dg_inner_dt<B: BesselianSource>(
        source: &B,
        observer: &umbra_core::ellipsoid::GeocentricObserver,
        east_longitude: Radians,
        t: TtInstant,
    ) -> f64 {
        let h_s = 1.0; // 1 秒
        let h_d = h_s / SECONDS_PER_DAY;
        let tp = tt_jd(jd_of(t) + h_d);
        let tm = tt_jd(jd_of(t) - h_d);
        (g_inner(source, observer, east_longitude, tp)
            - g_inner(source, observer, east_longitude, tm))
            / (2.0 * h_s)
    }

    /// 遷移方向テスト用の ε（day 換算）。8 秒 ≈ 9.26e-5 day。
    /// root_tolerance 1e-9 day（≈8.6e-5 s）より十分大きく、粗走査刻み 30s より小さい。
    /// これにより接触点を挟んだ符号判定が「数値ゼロ近傍の揺らぎ」でなく確かな符号反転を見る。
    const EPS_SECONDS: f64 = 8.0;
    fn eps_day() -> f64 {
        EPS_SECONDS / SECONDS_PER_DAY
    }

    // ============================================================
    // テスト本体
    // ============================================================

    /// 前提の独立確認（メタ）: 中心食地点では探索窓内で g_outer が「正→負→正」（C1, C4 を挟む）に
    /// 加え g_inner も符号反転する（C2/C3 が存在）。これは「選んだ地点・窓が実際に皆既を含む」ことの
    /// 独立チェックで、後続の None/順序テストの妥当性を支える。実 ephemeris を使う。
    #[test]
    fn central_site_window_brackets_both_outer_and_inner_contacts() {
        let dt = EspenakMeeusDeltaT;
        let src = make_source(&dt);
        let (obs, lon) = central_observer();
        let w = window_around(tt_2017_max(), 3.0);

        // 窓を細かく走査し g_outer / g_inner の符号反転回数を数える。
        let n = 2000;
        let (lo, hi) = (jd_of(w.start), jd_of(w.end));
        let mut outer_changes = 0;
        let mut inner_changes = 0;
        let mut prev_o = g_outer(&src, &obs, lon, tt_jd(lo));
        let mut prev_i = g_inner(&src, &obs, lon, tt_jd(lo));
        for k in 1..=n {
            let jd = lo + (hi - lo) * (f64::from(k) / f64::from(n));
            let o = g_outer(&src, &obs, lon, tt_jd(jd));
            let i = g_inner(&src, &obs, lon, tt_jd(jd));
            if prev_o * o < 0.0 {
                outer_changes += 1;
            }
            if prev_i * i < 0.0 {
                inner_changes += 1;
            }
            prev_o = o;
            prev_i = i;
        }
        assert_eq!(
            outer_changes, 2,
            "central site: g_outer should change sign twice (C1, C4), got {outer_changes}"
        );
        assert_eq!(
            inner_changes, 2,
            "central site: g_inner should change sign twice (C2, C3), got {inner_changes}"
        );
    }

    /// 主オラクル1（外接ゼロ点）: 中心食地点で返った C1/C4 時刻で、独立再構成 g_outer ≈ 0。
    /// 「正しい m²=L1² を解いているか」を直接縛る。
    #[test]
    fn central_c1_c4_are_zeros_of_independent_g_outer() {
        let dt = EspenakMeeusDeltaT;
        let src = make_source(&dt);
        let (obs, lon) = central_observer();
        let set = solve_local_contacts(
            &src,
            &obs,
            lon,
            window_around(tt_2017_max(), 3.0),
            config_tight(),
        )
        .expect("central site should yield a contact set");

        let c1 = set.c1.expect("central site must have C1");
        let c4 = set.c4.expect("central site must have C4");
        for (label, c) in [("C1", c1), ("C4", c4)] {
            let g = g_outer(&src, &obs, lon, c.time_tt);
            assert!(
                g.abs() < G_ROOT_TOL,
                "{label}: g_outer({}) = {g} not ≈ 0 (tol {G_ROOT_TOL})",
                jd_of(c.time_tt)
            );
        }
    }

    /// 主オラクル1（内接ゼロ点）: 中心食地点で返った C2/C3 時刻で、独立再構成 g_inner ≈ 0。
    #[test]
    fn central_c2_c3_are_zeros_of_independent_g_inner() {
        let dt = EspenakMeeusDeltaT;
        let src = make_source(&dt);
        let (obs, lon) = central_observer();
        let set = solve_local_contacts(
            &src,
            &obs,
            lon,
            window_around(tt_2017_max(), 3.0),
            config_tight(),
        )
        .expect("central site should yield a contact set");

        let c2 = set.c2.expect("central site must have C2");
        let c3 = set.c3.expect("central site must have C3");
        for (label, c) in [("C2", c2), ("C3", c3)] {
            let g = g_inner(&src, &obs, lon, c.time_tt);
            assert!(
                g.abs() < G_ROOT_TOL,
                "{label}: g_inner({}) = {g} not ≈ 0 (tol {G_ROOT_TOL})",
                jd_of(c.time_tt)
            );
        }
    }

    /// 主オラクル1 を時刻誤差に換算した強い縛り: C1 のゼロ点残差を局所勾配 dg/dt で割った
    /// 「等価時刻誤差」が root_tolerance スケール（≤0.2s に余裕を見て 0.5s）以内。
    /// g の絶対許容が緩すぎないことを、勾配経由で時刻 ±0.5s に翻訳して保証する。
    #[test]
    fn central_c1_residual_implies_subsecond_time_error() {
        let dt = EspenakMeeusDeltaT;
        let src = make_source(&dt);
        let (obs, lon) = central_observer();
        let set = solve_local_contacts(
            &src,
            &obs,
            lon,
            window_around(tt_2017_max(), 3.0),
            config_tight(),
        )
        .expect("central site should yield a contact set");
        let c1 = set.c1.expect("central site must have C1");

        let g = g_outer(&src, &obs, lon, c1.time_tt);
        let slope = dg_outer_dt(&src, &obs, lon, c1.time_tt); // Re²/s
        assert!(
            slope.abs() > 1e-9,
            "dg/dt at C1 is ~0 ({slope}); cannot translate residual to time"
        );
        let time_error_s = (g / slope).abs();
        assert!(
            time_error_s < 0.5,
            "C1 equivalent time error {time_error_s} s exceeds 0.5 s (g={g}, dg/dt={slope})"
        );
    }

    /// 補助オラクル2（順序プロパティ, L8）: 中心食地点で c1 < c2 < c3 < c4（最大食は c2,c3 の間）。
    #[test]
    fn central_contacts_are_time_ordered() {
        let dt = EspenakMeeusDeltaT;
        let src = make_source(&dt);
        let (obs, lon) = central_observer();
        let set = solve_local_contacts(
            &src,
            &obs,
            lon,
            window_around(tt_2017_max(), 3.0),
            config_tight(),
        )
        .expect("central site should yield a contact set");

        let c1 = jd_of(set.c1.expect("C1").time_tt);
        let c2 = jd_of(set.c2.expect("C2").time_tt);
        let c3 = jd_of(set.c3.expect("C3").time_tt);
        let c4 = jd_of(set.c4.expect("C4").time_tt);
        assert!(c1 < c2, "expected c1 < c2, got {c1} !< {c2}");
        assert!(c2 < c3, "expected c2 < c3, got {c2} !< {c3}");
        assert!(c3 < c4, "expected c3 < c4, got {c3} !< {c4}");
    }

    /// 補助オラクル2（None 分岐）: 部分食のみの地点では c2 == None && c3 == None、c1/c4 は存在。
    #[test]
    fn partial_site_has_no_inner_contacts_but_has_outer() {
        let dt = EspenakMeeusDeltaT;
        let src = make_source(&dt);
        let (obs, lon) = partial_observer();
        let set = solve_local_contacts(
            &src,
            &obs,
            lon,
            window_around(tt_2017_max(), 3.0),
            config_tight(),
        )
        .expect("partial site should yield a contact set");

        assert!(
            set.c2.is_none(),
            "partial site must have c2 == None, got {:?}",
            set.c2
        );
        assert!(
            set.c3.is_none(),
            "partial site must have c3 == None, got {:?}",
            set.c3
        );
        assert!(
            set.c1.is_some(),
            "partial site must still have C1 (partial eclipse begins)"
        );
        assert!(
            set.c4.is_some(),
            "partial site must still have C4 (partial eclipse ends)"
        );
    }

    /// 補助オラクル2 + 主オラクル1 の合わせ: 部分食地点の C1/C4 も g_outer ≈ 0、かつ c1 < c4。
    #[test]
    fn partial_site_outer_contacts_are_zeros_and_ordered() {
        let dt = EspenakMeeusDeltaT;
        let src = make_source(&dt);
        let (obs, lon) = partial_observer();
        let set = solve_local_contacts(
            &src,
            &obs,
            lon,
            window_around(tt_2017_max(), 3.0),
            config_tight(),
        )
        .expect("partial site should yield a contact set");

        let c1 = set.c1.expect("C1");
        let c4 = set.c4.expect("C4");
        for (label, c) in [("C1", c1), ("C4", c4)] {
            let g = g_outer(&src, &obs, lon, c.time_tt);
            assert!(g.abs() < G_ROOT_TOL, "{label}: g_outer = {g} not ≈ 0");
        }
        assert!(
            jd_of(c1.time_tt) < jd_of(c4.time_tt),
            "partial site: expected c1 < c4"
        );
    }

    /// 補助オラクル3（TT/UTC 整合）: 全接触で time_utc == tt_to_utc(time_tt)（過去日食は変換可能）。
    #[test]
    fn contact_utc_matches_tt_to_utc() {
        let dt = EspenakMeeusDeltaT;
        let src = make_source(&dt);
        let (obs, lon) = central_observer();
        let set = solve_local_contacts(
            &src,
            &obs,
            lon,
            window_around(tt_2017_max(), 3.0),
            config_tight(),
        )
        .expect("central site should yield a contact set");

        for (label, c) in [
            ("c1", set.c1),
            ("c2", set.c2),
            ("c3", set.c3),
            ("c4", set.c4),
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

    /// 補助オラクル4（契約: 食なしは正常）: 接触を一切含まない探索窓（日食の数日前で月影が遠い）では
    /// **`Ok` で全接触 None** を返す（「窓内に該当接触なし＝その地点は食なし」は正常でありエラーにしない）。
    /// `RootNotBracketed` は返さない。テスト側オラクルで両端 g_outer>0 を前提確認してから縛る。
    #[test]
    fn window_without_contacts_returns_all_none() {
        let dt = EspenakMeeusDeltaT;
        let src = make_source(&dt);
        let (obs, lon) = central_observer();
        // 最大食の 2 日前を中心に ±30 分。この窓には接触が無い（部分食の数時間内に入らない）。
        let center = tt_jd(jd_of(tt_2017_max()) - 2.0);
        let w = window_around(center, 0.5);

        // テスト側オラクルで「この窓は全域 外接 g_outer>0（部分食域の外）」を独立に確認（前提保証）。
        // 両端だけでなく内部も走査し、符号反転が無い＝接触が無いことを担保する。
        let n = 200;
        let (lo, hi) = (jd_of(w.start), jd_of(w.end));
        for k in 0..=n {
            let jd = lo + (hi - lo) * (f64::from(k) / f64::from(n));
            let g = g_outer(&src, &obs, lon, tt_jd(jd));
            assert!(
                g > 0.0,
                "precondition: window must be fully outside partial phase (g_outer>0 everywhere), got {g} at jd {jd}"
            );
        }

        let set = solve_local_contacts(&src, &obs, lon, w, config_tight())
            .expect("window without contacts is a no-eclipse case and must be Ok, not Err");
        assert_eq!(
            set,
            LocalContactSet::default(),
            "window without contacts must yield all-None set, got {set:?}"
        );
    }

    /// 補助オラクル5（grazing 刻み感度・偽陰性ガード）: 部分食の C1/C4 は外接（掠め接触に相当する
    /// 浅い符号変化）。**粗走査の刻みに過敏でないこと**を、探索窓を広めに取った設定で C1/C4 を
    /// 取りこぼさずに検出できることで縛る。実装が粗走査刻みを十分細かく取らないと、外接の浅い
    /// 符号変化を飛ばして RootNotBracketed や None になる（偽陰性）。
    ///
    /// 中心食地点の C1〜C4（外接 2 + 内接 2）が、最大食から離れた広い窓（±3h）でも全て検出される
    /// ことを要求する。窓が広い＝粗走査の 1 区間が長い＝接触を飛ばしやすい設定なので偽陰性ガードになる。
    #[test]
    fn coarse_scan_does_not_miss_contacts_in_wide_window() {
        let dt = EspenakMeeusDeltaT;
        let src = make_source(&dt);
        let (obs, lon) = central_observer();
        // 広い窓（±3h, 6 時間幅）。粗走査刻みが粗いと外接の浅い符号変化を飛ばす危険が最大。
        let set = solve_local_contacts(
            &src,
            &obs,
            lon,
            window_around(tt_2017_max(), 3.0),
            config_tight(),
        )
        .expect("wide window should still resolve all contacts");
        assert!(
            set.c1.is_some(),
            "C1 missed in wide window (coarse-scan false negative)"
        );
        assert!(set.c2.is_some(), "C2 missed in wide window");
        assert!(set.c3.is_some(), "C3 missed in wide window");
        assert!(set.c4.is_some(), "C4 missed in wide window");
    }

    /// 補助オラクル6（実日食 ballpark, 第二義・緩め）: 2017-08-21 中心食地点(37.5N)で
    ///  - 部分食継続 (C4 − C1) が ~2〜3.5 時間（北米中緯度の典型）
    ///  - 皆既継続 (C3 − C2) が **1〜5 分**（37.5N は ~2.6 分。下限>1分で C2≈C3 の取り違えを弾く）
    ///  - 全接触が 2017-08-21（当日, UTC-JD 2457986.5〜2457987.5）にある
    ///
    /// flaky なハードコード秒値は使わず、桁・継続スケールのみで縛る。
    #[test]
    fn central_site_2017_durations_are_in_ballpark() {
        let dt = EspenakMeeusDeltaT;
        let src = make_source(&dt);
        let (obs, lon) = central_observer();
        let set = solve_local_contacts(
            &src,
            &obs,
            lon,
            window_around(tt_2017_max(), 3.0),
            config_tight(),
        )
        .expect("central site should yield a contact set");

        let c1 = jd_of(set.c1.expect("C1").time_tt);
        let c2 = jd_of(set.c2.expect("C2").time_tt);
        let c3 = jd_of(set.c3.expect("C3").time_tt);
        let c4 = jd_of(set.c4.expect("C4").time_tt);

        let partial_hours = (c4 - c1) * 24.0;
        assert!(
            (2.0..=3.5).contains(&partial_hours),
            "partial duration {partial_hours} h not in 2.0–3.5 h ballpark"
        );
        let total_minutes = (c3 - c2) * 24.0 * 60.0;
        assert!(
            (1.0..=5.0).contains(&total_minutes),
            "totality duration {total_minutes} min not in 1–5 min ballpark (>1 min rejects C2≈C3 mix-up)"
        );

        // 2017-08-21 当日（UTC-JD 2457986.5 ≤ 当日 < 2457987.5、TT は +ΔT≈8e-4 日内で同日扱い）。
        for (label, jd) in [("c1", c1), ("c2", c2), ("c3", c3), ("c4", c4)] {
            assert!(
                (2_457_986.5..2_457_987.6).contains(&jd),
                "{label} TT-JD {jd} not on 2017-08-21"
            );
        }
    }

    /// [要修正-1] 接触種別の遷移方向を独立に縛る（順序だけでは検出できない C1↔C4 / C2↔C3 取り違えを弾く）。
    /// 中心食地点(37.5N)で返った各接触時刻 t の前後 ±ε（8s 相当 day, root_tolerance 1e-9 day より大・
    /// 走査刻み 30s より小）で、独立オラクル g の符号方向を検証する:
    ///   C1: g_outer(c1−ε)>0 && g_outer(c1+ε)<0（外側→内側＝部分食開始）
    ///   C4: g_outer(c4−ε)<0 && g_outer(c4+ε)>0（内側→外側＝部分食終了）
    ///   C2: g_inner(c2−ε)>0 && g_inner(c2+ε)<0
    ///   C3: g_inner(c3−ε)<0 && g_inner(c3+ε)>0
    #[test]
    fn contact_kinds_have_correct_sign_transitions() {
        let dt = EspenakMeeusDeltaT;
        let src = make_source(&dt);
        let (obs, lon) = central_observer();
        let set = solve_local_contacts(
            &src,
            &obs,
            lon,
            window_around(tt_2017_max(), 3.0),
            config_tight(),
        )
        .expect("central site should yield a contact set");

        let eps = eps_day();
        let before = |t: TtInstant| tt_jd(jd_of(t) - eps);
        let after = |t: TtInstant| tt_jd(jd_of(t) + eps);

        let c1 = set.c1.expect("C1").time_tt;
        let c4 = set.c4.expect("C4").time_tt;
        let c2 = set.c2.expect("C2").time_tt;
        let c3 = set.c3.expect("C3").time_tt;

        // C1: 部分食開始（外側 g_outer>0 → 内側 g_outer<0）。
        let c1_b = g_outer(&src, &obs, lon, before(c1));
        let c1_a = g_outer(&src, &obs, lon, after(c1));
        assert!(
            c1_b > 0.0 && c1_a < 0.0,
            "C1 must be outer +→− (g_outer before={c1_b}, after={c1_a})"
        );

        // C4: 部分食終了（内側 g_outer<0 → 外側 g_outer>0）。
        let c4_b = g_outer(&src, &obs, lon, before(c4));
        let c4_a = g_outer(&src, &obs, lon, after(c4));
        assert!(
            c4_b < 0.0 && c4_a > 0.0,
            "C4 must be outer −→+ (g_outer before={c4_b}, after={c4_a})"
        );

        // C2: 皆既開始（内接 g_inner>0 → g_inner<0）。
        let c2_b = g_inner(&src, &obs, lon, before(c2));
        let c2_a = g_inner(&src, &obs, lon, after(c2));
        assert!(
            c2_b > 0.0 && c2_a < 0.0,
            "C2 must be inner +→− (g_inner before={c2_b}, after={c2_a})"
        );

        // C3: 皆既終了（内接 g_inner<0 → g_inner>0）。
        let c3_b = g_inner(&src, &obs, lon, before(c3));
        let c3_a = g_inner(&src, &obs, lon, after(c3));
        assert!(
            c3_b < 0.0 && c3_a > 0.0,
            "C3 must be inner −→+ (g_inner before={c3_b}, after={c3_a})"
        );
    }

    /// [要修正-2] 偽陰性ガード（複数窓幅で同一接触）: 探索窓を ±3h / ±2h / ±1.5h の 3 通りで解き、
    /// 得られる C1〜C4（time_tt）が互いに小さな許容（2 秒）内で一致することを検証する。
    /// 窓幅が変わると実装内部の粗走査分割数（窓秒/30s）が変わるため、刻み感度の偽陰性を実効的に縛る。
    /// 3 窓とも半幅 > 90 分とし、中心食地点の部分食半継続（~80 分）より大きくして C1/C4 を必ず含める
    /// （半幅が部分食半継続より小さいと窓外の接触が None になり、窓幅不変性ではなく窓外契約に落ちてしまう）。
    #[test]
    fn contacts_are_window_width_invariant() {
        let dt = EspenakMeeusDeltaT;
        let src = make_source(&dt);
        let (obs, lon) = central_observer();

        let solve = |half_hours: f64| {
            solve_local_contacts(
                &src,
                &obs,
                lon,
                window_around(tt_2017_max(), half_hours),
                config_tight(),
            )
            .expect("central site should yield a contact set")
        };

        // 37.5N の部分食半継続は実測 >90 分（部分食全体 >3h）。3 窓とも半幅 ≥120 分で C1/C4 を確実に含める。
        let wide = solve(3.0); // ±3h（半幅 180 分）
        let mid = solve(2.5); // ±2.5h（半幅 150 分）
        let narrow = solve(2.0); // ±2h（半幅 120 分 > 部分食半継続）

        // 2 秒 = 2/86400 day。窓幅由来の刻み感度を弾く（root_tolerance 1e-9 day より十分緩い実用許容）。
        let tol = 2.0 / SECONDS_PER_DAY;
        let pick = |s: &LocalContactSet| {
            [
                jd_of(s.c1.expect("C1").time_tt),
                jd_of(s.c2.expect("C2").time_tt),
                jd_of(s.c3.expect("C3").time_tt),
                jd_of(s.c4.expect("C4").time_tt),
            ]
        };
        let (a, b, c) = (pick(&wide), pick(&mid), pick(&narrow));
        for (k, label) in ["c1", "c2", "c3", "c4"].iter().enumerate() {
            assert!(
                (a[k] - b[k]).abs() < tol,
                "{label}: ±3h vs ±2.5h differ by {} s",
                (a[k] - b[k]).abs() * SECONDS_PER_DAY
            );
            assert!(
                (a[k] - c[k]).abs() < tol,
                "{label}: ±3h vs ±2h differ by {} s",
                (a[k] - c[k]).abs() * SECONDS_PER_DAY
            );
        }
    }

    /// [要修正-3a] 限界内側 grazing 地点(38.0N): C2/C3 が存在し、皆既継続 (c3−c2) が
    /// 中心(37.5N)より短いこと、かつ c1<c2<c3<c4。本影の縁寄りで皆既が浅いケースの縛り。
    #[test]
    fn inner_limit_site_has_shorter_totality_than_center() {
        let dt = EspenakMeeusDeltaT;
        let src = make_source(&dt);

        let solve = |obs, lon| {
            solve_local_contacts(
                &src,
                &obs,
                lon,
                window_around(tt_2017_max(), 3.0),
                config_tight(),
            )
            .expect("site should yield a contact set")
        };

        let (cobs, clon) = central_observer();
        let center = solve(cobs, clon);
        let (iobs, ilon) = inner_limit_observer();
        let inner = solve(iobs, ilon);

        // 内側 grazing でも C1〜C4 すべて存在し時刻順。
        let i1 = jd_of(inner.c1.expect("inner-limit C1").time_tt);
        let i2 = jd_of(inner.c2.expect("inner-limit C2").time_tt);
        let i3 = jd_of(inner.c3.expect("inner-limit C3").time_tt);
        let i4 = jd_of(inner.c4.expect("inner-limit C4").time_tt);
        assert!(i1 < i2 && i2 < i3 && i3 < i4, "inner-limit not ordered");

        let center_tot = jd_of(center.c3.expect("center C3").time_tt)
            - jd_of(center.c2.expect("center C2").time_tt);
        let inner_tot = i3 - i2;
        assert!(
            inner_tot > 0.0,
            "inner-limit totality must be positive, got {} min",
            inner_tot * 24.0 * 60.0
        );
        assert!(
            inner_tot < center_tot,
            "inner-limit(38.0N) totality {} min must be shorter than center(37.5N) {} min",
            inner_tot * 24.0 * 60.0,
            center_tot * 24.0 * 60.0
        );
    }

    /// [要修正-3b] 限界外側 grazing 地点(38.5N): 本影を僅かに外す → **c2==None && c3==None だが
    /// c1/c4 は Some**（部分食・本影を僅かに外す）。「中心線近傍でも本影を外せば内接なし」の限界ケース。
    #[test]
    fn outer_limit_site_has_no_inner_but_has_outer() {
        let dt = EspenakMeeusDeltaT;
        let src = make_source(&dt);
        let (obs, lon) = outer_limit_observer();
        let set = solve_local_contacts(
            &src,
            &obs,
            lon,
            window_around(tt_2017_max(), 3.0),
            config_tight(),
        )
        .expect("outer-limit site should yield a contact set");

        assert!(
            set.c2.is_none() && set.c3.is_none(),
            "outer-limit(38.5N) just outside umbra must have c2=c3=None, got c2={:?} c3={:?}",
            set.c2,
            set.c3
        );
        let c1 = set.c1.expect("outer-limit C1 (partial begins)");
        let c4 = set.c4.expect("outer-limit C4 (partial ends)");
        assert!(
            jd_of(c1.time_tt) < jd_of(c4.time_tt),
            "outer-limit: expected c1 < c4"
        );
    }

    /// [要修正-4] 部分食地点の前提を独立確認するメタテスト: 部分食地点(60N,−100°)で窓全域を走査し
    /// **g_outer が 2 回符号反転（C1/C4 存在）かつ g_inner>0 が全域**（内接接触なし＝真に部分食）を
    /// 独立に確認する。これで `partial_site_has_no_inner_contacts_but_has_outer` の None が
    /// 「実装の取りこぼし」を追認していないことを担保する。
    #[test]
    fn partial_site_window_has_outer_brackets_and_no_inner() {
        let dt = EspenakMeeusDeltaT;
        let src = make_source(&dt);
        let (obs, lon) = partial_observer();
        let w = window_around(tt_2017_max(), 3.0);

        let n = 2000;
        let (lo, hi) = (jd_of(w.start), jd_of(w.end));
        let mut outer_changes = 0;
        let mut min_inner = f64::INFINITY;
        let mut prev_o = g_outer(&src, &obs, lon, tt_jd(lo));
        for k in 1..=n {
            let jd = lo + (hi - lo) * (f64::from(k) / f64::from(n));
            let o = g_outer(&src, &obs, lon, tt_jd(jd));
            let i = g_inner(&src, &obs, lon, tt_jd(jd));
            if prev_o * o < 0.0 {
                outer_changes += 1;
            }
            if i < min_inner {
                min_inner = i;
            }
            prev_o = o;
        }
        assert_eq!(
            outer_changes, 2,
            "partial site: g_outer should change sign twice (C1, C4), got {outer_changes}"
        );
        assert!(
            min_inner > 0.0,
            "partial site: g_inner must stay > 0 over whole window (no inner contact), min was {min_inner}"
        );
    }

    /// [推奨] C2 の残差→等価時刻誤差: 中心食地点で C2 の g_inner 残差を局所勾配 dg_inner/dt で
    /// 割った等価時刻誤差が 0.5s 以内であること。内接は外接と勾配が異なるため別途必要
    /// （central_c1_residual_implies_subsecond_time_error の内接版）。
    #[test]
    fn central_c2_residual_implies_subsecond_time_error() {
        let dt = EspenakMeeusDeltaT;
        let src = make_source(&dt);
        let (obs, lon) = central_observer();
        let set = solve_local_contacts(
            &src,
            &obs,
            lon,
            window_around(tt_2017_max(), 3.0),
            config_tight(),
        )
        .expect("central site should yield a contact set");
        let c2 = set.c2.expect("central site must have C2");

        let g = g_inner(&src, &obs, lon, c2.time_tt);
        let slope = dg_inner_dt(&src, &obs, lon, c2.time_tt); // Re²/s
        assert!(
            slope.abs() > 1e-9,
            "dg_inner/dt at C2 is ~0 ({slope}); cannot translate residual to time"
        );
        let time_error_s = (g / slope).abs();
        assert!(
            time_error_s < 0.5,
            "C2 equivalent time error {time_error_s} s exceeds 0.5 s (g={g}, dg_inner/dt={slope})"
        );
    }
}
