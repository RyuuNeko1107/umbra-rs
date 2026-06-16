//! ベッセル要素の供給源抽象（`docs/issues/ISSUE-037`, `docs/api-draft.md §3.3`,
//! `docs/architecture.md §6.1`）。
//!
//! [`BesselianSource`] は任意 TT 時刻で [`InstantaneousBesselianElements`] を供給する抽象で、
//! 直接評価（[`DirectBesselianSource`]・fit 誤差ゼロ）と多項式近似（ISSUE-022, `BesselianPolynomial`）を
//! `&dyn BesselianSource` で差し替え可能にする（局地 solver が供給源にジェネリック, architecture §6.1）。
//! 供給源を替えても単位・座標系は不変（x,y,l1,l2 = Re 無次元, d,μ = rad, FundamentalPlane, TT 基準）。
//!
//! [`DirectBesselianSource`] は各 `at()` 呼出で [`besselian_elements_at`] を再評価する
//! （= ISSUE-015 見かけ位置 → ISSUE-019 影円錐 → ISSUE-020 基底 → ISSUE-021 瞬時要素のパイプラインを
//! 毎回回す）。多項式 fit 誤差を持たない（精度最優先, architecture §6.1）。Standard 局地計算の
//! 既定供給源。
//!
//! 注: ISSUE-037 の `InstantaneousEvaluator<E: Ephemeris>`（`EngineConfig`/`TimeScales`/
//! `AstrometryOptions` 統合・暦ジェネリック）は、それら公開型（ISSUE-042/043/015 公開API）と
//! `apparent` 層の Ephemeris ジェネリック化が整うまで繰延。現アーキ（`apparent::*_apparent_cirs` が
//! VSOP/ELP 直結）に合わせ、本スライスは半径・ΔT モデルを保持する直接供給源を提供する。

use umbra_core::deltat::DeltaTModel;
use umbra_core::{TimeInterval, TtInstant};

use crate::besselian::{besselian_elements_at, InstantaneousBesselianElements};
use crate::error::EclipseError;

/// 任意 TT 時刻で瞬時ベッセル要素を供給する抽象。
///
/// 直接評価（[`DirectBesselianSource`], fit 誤差ゼロ, ISSUE-037）と多項式近似
/// （`BesselianPolynomial`, ISSUE-022）を同一契約で差し替え可能にする（`&dyn BesselianSource`,
/// architecture §6.1）。
pub trait BesselianSource {
    /// 時刻 `time`（TT）における瞬時ベッセル要素を返す。
    fn at(&self, time: TtInstant) -> Result<InstantaneousBesselianElements, EclipseError>;
    /// この供給源が妥当な TT 区間（多項式は fit 区間、直接評価は構築時に渡された推奨範囲）。
    fn fit_interval(&self) -> TimeInterval<TtInstant>;
}

/// 各 `at()` 呼出で暦を直接再評価する供給源（**fit 誤差ゼロ**, ISSUE-037）。
///
/// [`besselian_elements_at`]（ISSUE-021）を保持した半径・ΔT モデルで毎回評価する。
/// `at()` は区間外でも評価可能（暦が有効な範囲なら値を返す）。`fit_interval` は推奨範囲を広告するのみ。
#[derive(Debug, Clone, Copy)]
pub struct DirectBesselianSource<'d, M: DeltaTModel> {
    /// 太陽物理半径[km]。
    r_sun_km: f64,
    /// 月半径[km]（= k·Re）。
    r_moon_km: f64,
    /// ΔT モデル（μ の UT1 変換に使用）。
    delta_t: &'d M,
    /// この供給源の妥当（推奨）TT 区間。
    interval: TimeInterval<TtInstant>,
}

impl<'d, M: DeltaTModel> DirectBesselianSource<'d, M> {
    /// `r_sun_km` = 太陽物理半径[km], `r_moon_km` = 月半径[km]（= k·Re）,
    /// `delta_t` = ΔT モデル, `interval` = この供給源の妥当（推奨）TT 区間。
    pub fn new(
        r_sun_km: f64,
        r_moon_km: f64,
        delta_t: &'d M,
        interval: TimeInterval<TtInstant>,
    ) -> Self {
        Self {
            r_sun_km,
            r_moon_km,
            delta_t,
            interval,
        }
    }
}

impl<M: DeltaTModel> BesselianSource for DirectBesselianSource<'_, M> {
    /// 各呼出で [`besselian_elements_at`] を再評価（fit 誤差ゼロ）。式は一切加えない。
    fn at(&self, time: TtInstant) -> Result<InstantaneousBesselianElements, EclipseError> {
        besselian_elements_at(time, self.r_sun_km, self.r_moon_km, self.delta_t)
    }

    /// 構築時に渡された推奨区間をそのまま返す。
    fn fit_interval(&self) -> TimeInterval<TtInstant> {
        self.interval
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use umbra_core::constants::{EARTH_EQUATORIAL_RADIUS_M, SOLAR_RADIUS_KM};
    use umbra_core::{EspenakMeeusDeltaT, JulianDate2, TimeInterval, TtInstant};

    use crate::besselian::besselian_elements_at;

    const R_SUN: f64 = SOLAR_RADIUS_KM;
    const R_MOON: f64 = 0.2725076 * (EARTH_EQUATORIAL_RADIUS_M / 1000.0);

    /// TT 時刻を 2 要素 JD から構築するヘルパ。
    fn tt(jd1: f64, jd2: f64) -> TtInstant {
        TtInstant::from_jd2(JulianDate2::new(jd1, jd2))
    }

    /// 2017-08-21 最大食付近の TT（besselian.rs テストと同一エポック）。
    fn tt_2017_max() -> TtInstant {
        tt(2_457_986.5, 7.685_322_222_222_222e-1)
    }

    /// J2000.0（TT）。
    fn tt_j2000() -> TtInstant {
        tt(2_451_545.0, 0.0)
    }

    /// 値同一性検証で使う「複数の異なる時刻」。実太陽・月位置で評価可能なエポックを散らす。
    fn sample_times() -> [TtInstant; 4] {
        [
            tt_2017_max(),
            tt_j2000(),
            tt(2_458_000.5, 0.25),
            tt(2_459_731.75, 0.0),
        ]
    }

    /// `InstantaneousBesselianElements` の全フィールドが厳密（`==`）一致することを表明する。
    /// 同一計算経路の同一性検証なので tol を使わずビット/厳密一致を要求する（契約1）。
    fn assert_exact_eq(
        got: &InstantaneousBesselianElements,
        want: &InstantaneousBesselianElements,
        label: &str,
    ) {
        assert_eq!(got.x, want.x, "{label}: x");
        assert_eq!(got.y, want.y, "{label}: y");
        assert_eq!(
            got.declination.0, want.declination.0,
            "{label}: declination"
        );
        assert_eq!(got.mu.0, want.mu.0, "{label}: mu");
        assert_eq!(got.l1, want.l1, "{label}: l1");
        assert_eq!(got.l2, want.l2, "{label}: l2");
        assert_eq!(got.tan_f1, want.tan_f1, "{label}: tan_f1");
        assert_eq!(got.tan_f2, want.tan_f2, "{label}: tan_f2");
        assert_eq!(got.time_tt, want.time_tt, "{label}: time_tt");
        // 派生値も一致（全フィールド一致なら自明だが gamma 経路を踏む）。
        assert_eq!(got.gamma(), want.gamma(), "{label}: gamma");
    }

    /// 区間内の任意時刻を含む、テスト用の妥当区間。
    fn interval() -> TimeInterval<TtInstant> {
        TimeInterval {
            start: tt(2_451_545.0, 0.0),
            end: tt(2_460_000.0, 0.0),
        }
    }

    /// 契約1（最重要）: `DirectBesselianSource::at(t)` は同じ引数で直接呼んだ
    /// `besselian_elements_at(t, r_sun, r_moon, &dt)` と全フィールド厳密一致する。
    /// オーケストレーション層が式を一切変えないことを複数時刻で担保する。
    #[test]
    fn at_matches_besselian_elements_at_for_multiple_times() {
        let dt = EspenakMeeusDeltaT;
        let src = DirectBesselianSource::new(R_SUN, R_MOON, &dt, interval());
        for t in sample_times() {
            let got = src
                .at(t)
                .expect("source at() should succeed at real positions");
            let want = besselian_elements_at(t, R_SUN, R_MOON, &dt)
                .expect("oracle besselian_elements_at should succeed");
            assert_exact_eq(&got, &want, "at()==besselian_elements_at");
        }
    }

    /// 契約1 系: 渡した半径引数（r_sun, r_moon）がそのまま暦評価に伝わることを、
    /// 別の半径値でも at()==besselian_elements_at が成り立つことで担保する。
    #[test]
    fn at_uses_provided_radii() {
        let dt = EspenakMeeusDeltaT;
        // 既定値と異なる半径でも、供給源は同じ引数で直接呼んだ結果と一致しなければならない。
        let r_sun = R_SUN * 1.01;
        let r_moon = R_MOON * 0.99;
        let src = DirectBesselianSource::new(r_sun, r_moon, &dt, interval());
        let t = tt_2017_max();
        let got = src.at(t).expect("at() should succeed");
        let want = besselian_elements_at(t, r_sun, r_moon, &dt).expect("oracle should succeed");
        assert_exact_eq(&got, &want, "at() honors custom radii");
    }

    /// 契約2: `new(.., interval)` に渡した `interval` が `fit_interval()` でそのまま
    /// （start/end とも）返る。
    #[test]
    fn fit_interval_returns_constructed_interval() {
        let dt = EspenakMeeusDeltaT;
        let iv = interval();
        let src = DirectBesselianSource::new(R_SUN, R_MOON, &dt, iv);
        let got = src.fit_interval();
        assert_eq!(got, iv, "fit_interval should echo the constructed interval");
        assert_eq!(got.start, iv.start, "fit_interval start");
        assert_eq!(got.end, iv.end, "fit_interval end");
    }

    /// 契約3: `fit_interval` の外の時刻でも at() はエラーにせず、
    /// besselian_elements_at と同じ値を返す（区間は推奨範囲を広告するだけ）。
    #[test]
    fn at_evaluates_outside_fit_interval() {
        let dt = EspenakMeeusDeltaT;
        // 区間を狭く取り、検証時刻が外側に来るようにする。
        let narrow = TimeInterval {
            start: tt(2_457_986.0, 0.0),
            end: tt(2_457_987.0, 0.0),
        };
        let src = DirectBesselianSource::new(R_SUN, R_MOON, &dt, narrow);

        // 区間外（過去側）の J2000 で評価する。
        let outside = tt_j2000();
        assert!(
            outside.jd2().jd() < narrow.start.jd2().jd(),
            "test time must be outside (before) the interval"
        );
        let got = src
            .at(outside)
            .expect("at() must not error outside the advertised interval");
        let want =
            besselian_elements_at(outside, R_SUN, R_MOON, &dt).expect("oracle should succeed");
        assert_exact_eq(
            &got,
            &want,
            "at() outside interval == besselian_elements_at",
        );
    }

    /// 契約4: `&dyn BesselianSource` 経由で at()/fit_interval() が呼べる（object-safe）。
    /// 多項式版と `&dyn BesselianSource` で差し替え可能であることの担保。
    #[test]
    fn usable_through_dyn_trait_object() {
        let dt = EspenakMeeusDeltaT;
        let iv = interval();
        let concrete = DirectBesselianSource::new(R_SUN, R_MOON, &dt, iv);
        let s: &dyn BesselianSource = &concrete;

        // fit_interval() を trait オブジェクト越しに。
        assert_eq!(s.fit_interval(), iv, "fit_interval via &dyn");

        // at() を trait オブジェクト越しに。直接呼びと厳密一致。
        let t = tt_2017_max();
        let got = s.at(t).expect("at() via &dyn should succeed");
        let want = besselian_elements_at(t, R_SUN, R_MOON, &dt).expect("oracle should succeed");
        assert_exact_eq(&got, &want, "at() via &dyn == besselian_elements_at");
    }
}
