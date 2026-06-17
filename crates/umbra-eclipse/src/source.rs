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
//! [`InstantaneousEvaluator`]`<E: Ephemeris>`（ISSUE-043 S3）は任意の暦バックエンドを駆動する
//! 暦ジェネリックな供給源で、`apparent` 層の Ephemeris ジェネリック化（ISSUE-043 S2
//! `apparent_cirs<E>` + `AstrometryOptions`）の上に構築される。`DirectBesselianSource` は
//! VSOP/ELP 直結（具象）の既定供給源、`InstantaneousEvaluator` は暦差し替え・Mock 幾何検証用。

use umbra_core::deltat::{tt_to_ut1, DeltaTModel};
use umbra_core::{TimeInterval, TtInstant};
use umbra_ephemeris::{apparent_cirs, AstrometryOptions, Body, Ephemeris};

use crate::besselian::{
    besselian_elements_at, instantaneous_from_cirs, InstantaneousBesselianElements,
};
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

/// 任意 [`Ephemeris`] バックエンドから各 `at()` で瞬時ベッセル要素を直接評価する供給源
/// （**fit 誤差ゼロ**, ISSUE-043 S3）。
///
/// [`DirectBesselianSource`] が VSOP/ELP 直結の `besselian_elements_at` を呼ぶのに対し、本型は
/// [`umbra_ephemeris::apparent_cirs`]`<E>`（ISSUE-043 S2・暦ジェネリック）で太陽・月の見かけ CIRS
/// 位置を `options` に従って得て、共有組立 [`instantaneous_from_cirs`] へ渡す。`AnalyticalEphemeris`
/// ＋ [`AstrometryOptions::standard`] では [`besselian_elements_at`] と一致し（apparent_cirs の回帰
/// ブリッジ済）、[`MockEphemeris`](umbra_ephemeris::MockEphemeris) ＋ [`AstrometryOptions::geometric`]
/// では幾何配置の検証に使える（EclipseEngine の Mock CI 経路, 受け入れ §77）。
///
/// `ephemeris`・`delta_t` は借用（EclipseEngine が保持する E・D を各呼出で借りる想定）。
#[derive(Debug, Clone, Copy)]
pub struct InstantaneousEvaluator<'e, 'd, E: Ephemeris, M: DeltaTModel> {
    /// 暦バックエンド（借用）。
    ephemeris: &'e E,
    /// ΔT モデル（借用, μ の UT1 変換に使用）。
    delta_t: &'d M,
    /// 太陽物理半径[km]。
    r_sun_km: f64,
    /// 月半径[km]（= k·Re）。
    r_moon_km: f64,
    /// 見かけ補正フラグ（標準は全 ON、Mock 幾何検証は全 OFF）。
    options: AstrometryOptions,
    /// この供給源の妥当（推奨）TT 区間。
    interval: TimeInterval<TtInstant>,
}

impl<'e, 'd, E: Ephemeris, M: DeltaTModel> InstantaneousEvaluator<'e, 'd, E, M> {
    /// `ephemeris` = 暦バックエンド（借用）, `delta_t` = ΔT モデル（借用）,
    /// `r_sun_km` = 太陽物理半径, `r_moon_km` = 月半径（= k·Re）, `options` = 見かけ補正フラグ,
    /// `interval` = この供給源の妥当（推奨）TT 区間。
    pub fn new(
        ephemeris: &'e E,
        delta_t: &'d M,
        r_sun_km: f64,
        r_moon_km: f64,
        options: AstrometryOptions,
        interval: TimeInterval<TtInstant>,
    ) -> Self {
        Self {
            ephemeris,
            delta_t,
            r_sun_km,
            r_moon_km,
            options,
            interval,
        }
    }
}

impl<E: Ephemeris, M: DeltaTModel> BesselianSource for InstantaneousEvaluator<'_, '_, E, M> {
    /// 各呼出で `apparent_cirs<E>` から太陽・月の見かけ CIRS 位置を得て瞬時要素を組み立てる
    /// （fit 誤差ゼロ）。暦の `EphemerisError` は `?` で [`EclipseError::Ephemeris`] へ透過。
    fn at(&self, time: TtInstant) -> Result<InstantaneousBesselianElements, EclipseError> {
        let sun = apparent_cirs(self.ephemeris, Body::Sun, time, self.options)?;
        let moon = apparent_cirs(self.ephemeris, Body::Moon, time, self.options)?;
        let ut1 = tt_to_ut1(time, self.delta_t);
        instantaneous_from_cirs(sun, moon, self.r_sun_km, self.r_moon_km, time, ut1)
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

// ============================================================
// ISSUE-043 S3: InstantaneousEvaluator<E: Ephemeris>（暦ジェネリックな直接供給源）
// ============================================================
//
// `InstantaneousEvaluator` は `besselian_elements_at`（VSOP/ELP 直結）を `apparent_cirs<E>` で
// 暦ジェネリック化したもの。回帰の核心は「AnalyticalEphemeris + standard では
// besselian_elements_at と**ビット級厳密一致**する」こと（S2 回帰ブリッジ済みゆえ）。
// Mock 統合では幾何配置がエンジン供給源として正しく流れることを設計オラクルで縛る。
#[cfg(test)]
mod evaluator_tests {
    use super::*;

    use umbra_core::constants::{EARTH_EQUATORIAL_RADIUS_M, SOLAR_RADIUS_KM};
    use umbra_core::{EspenakMeeusDeltaT, JulianDate2, TimeInterval, TtInstant};
    use umbra_ephemeris::{AnalyticalEphemeris, AstrometryOptions, EphemerisError, MockEphemeris};

    use crate::besselian::{besselian_elements_at, InstantaneousBesselianElements};

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

    /// 回帰用に複数の異なる実エポックを散らす（DirectBesselianSource テストと同一集合）。
    fn sample_times() -> [TtInstant; 4] {
        [
            tt_2017_max(),
            tt_j2000(),
            tt(2_458_000.5, 0.25),
            tt(2_459_731.75, 0.0),
        ]
    }

    /// 区間内の任意時刻を含む、テスト用の妥当区間。
    fn interval() -> TimeInterval<TtInstant> {
        TimeInterval {
            start: tt(2_451_545.0, 0.0),
            end: tt(2_460_000.0, 0.0),
        }
    }

    /// `InstantaneousBesselianElements` の全フィールド厳密（`==`）一致を表明する。
    /// 同一計算経路の同一性検証なので tol を使わずビット/厳密一致を要求する（回帰の核心）。
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
        // 派生値も一致（gamma 経路を踏む）。
        assert_eq!(got.gamma(), want.gamma(), "{label}: gamma");
    }

    /// 回帰（最重要）: 評価器(AnalyticalEphemeris + standard).at(t) は
    /// `besselian_elements_at(t, R_SUN, R_MOON, &EspenakMeeusDeltaT)` と**全フィールド厳密一致**する。
    /// apparent_cirs(&Analytical, standard) == 具象 *_apparent_cirs（S2 回帰ブリッジ済）ゆえ、
    /// 評価器の組立 == besselian_elements_at がビット級で一致するはず。
    /// 殺す変異: apparent 取得元の取り違え・options 既定値改変・組立順や引数差し替え・μ の UT1 変換漏れ。
    #[test]
    fn evaluator_analytical_standard_matches_besselian_elements_at() {
        let eph = AnalyticalEphemeris::new();
        let dt = EspenakMeeusDeltaT;
        let eval = InstantaneousEvaluator::new(
            &eph,
            &dt,
            R_SUN,
            R_MOON,
            AstrometryOptions::standard(),
            interval(),
        );
        for t in sample_times() {
            let got = eval
                .at(t)
                .expect("evaluator at() should succeed at real positions");
            let want = besselian_elements_at(t, R_SUN, R_MOON, &dt)
                .expect("oracle besselian_elements_at should succeed");
            assert_exact_eq(
                &got,
                &want,
                "evaluator(Analytical,standard)==besselian_elements_at",
            );
        }
    }

    /// 回帰系: 渡した半径引数（r_sun, r_moon）が暦評価へそのまま伝わる。既定と異なる半径でも
    /// 評価器(Analytical+standard) == besselian_elements_at（同一半径引数）が成り立つことで縛る。
    /// 殺す変異: r_sun/r_moon のハードコード化・取り違え・無視。
    #[test]
    fn evaluator_honors_provided_radii() {
        let eph = AnalyticalEphemeris::new();
        let dt = EspenakMeeusDeltaT;
        let r_sun = R_SUN * 1.01;
        let r_moon = R_MOON * 0.99;
        let eval = InstantaneousEvaluator::new(
            &eph,
            &dt,
            r_sun,
            r_moon,
            AstrometryOptions::standard(),
            interval(),
        );
        let t = tt_2017_max();
        let got = eval.at(t).expect("at() should succeed");
        let want = besselian_elements_at(t, r_sun, r_moon, &dt).expect("oracle should succeed");
        assert_exact_eq(&got, &want, "evaluator honors custom radii");
    }

    /// 実日食ゲート: 評価器(Analytical+standard).at(2017 最大食) の gamma が
    /// NASA 公表 gamma≈0.4367 と 4 桁一致（[0.43,0.44]）。apparent 経路の実日食検証。
    /// 殺す変異: apparent 補正経路の破壊・基本面射影/Re 正規化の誤り・time_tt ずれ。
    #[test]
    fn evaluator_2017_total_eclipse_gamma_matches_nasa() {
        let eph = AnalyticalEphemeris::new();
        let dt = EspenakMeeusDeltaT;
        let eval = InstantaneousEvaluator::new(
            &eph,
            &dt,
            R_SUN,
            R_MOON,
            AstrometryOptions::standard(),
            interval(),
        );
        let e = eval.at(tt_2017_max()).expect("at() should succeed");
        // NASA gamma=0.4367 を [0.43,0.44] で締める（モデル差 ΔT/k/平均月縁の余裕）。
        assert!(
            (0.43..0.44).contains(&e.gamma()),
            "gamma = {} (NASA 0.4367)",
            e.gamma()
        );
        // 念のため time_tt ラベルが入力 TT を保持。
        assert_eq!(e.time_tt, tt_2017_max(), "time_tt label preserved");
    }

    /// Mock 統合（central_total, geometric）: 影軸が地心近傍を貫く中心皆既配置ゆえ gamma≈0、
    /// 本影は皆既ゆえ l2<0。Mock 幾何 besselian がエンジン供給源として流れることを縛る。
    /// 設計オラクル（besselian.rs central_total テストと整合）: gamma<1e-6, l2<0, l1>0,
    /// |l2| が実日食域 [0.005,0.05) Re。geometric ゆえ velocity 不要で Ok。
    /// 殺す変異: Mock 配置の取り違え・geometric で誤って velocity を要求・l2 符号反転。
    #[test]
    fn evaluator_mock_central_total_geometric_design_oracle() {
        let eph = MockEphemeris::central_total();
        let dt = EspenakMeeusDeltaT;
        let eval = InstantaneousEvaluator::new(
            &eph,
            &dt,
            R_SUN,
            R_MOON,
            AstrometryOptions::geometric(),
            interval(),
        );
        let e = eval
            .at(tt_j2000())
            .expect("geometric Mock evaluation should succeed (no velocity needed)");
        assert!(e.gamma() < 1e-6, "central_total gamma = {}", e.gamma());
        assert!(e.l2 < 0.0, "central_total l2 = {} (皆既は負)", e.l2);
        assert!(e.l1 > 0.0, "l1 = {} (半影は正)", e.l1);
        assert!(
            (0.005..0.05).contains(&e.l2.abs()),
            "|l2| = {} (実日食域)",
            e.l2.abs()
        );
    }

    /// Mock 統合（clear_annular, geometric）: 中心配置 gamma≈0 だが遠地点で月が小さく金環 ⇒ l2>0。
    /// l2 符号オラクルを central_total（l2<0）と対で縛る（符号反転変異を殺す）。
    #[test]
    fn evaluator_mock_clear_annular_geometric_positive_l2() {
        let eph = MockEphemeris::clear_annular();
        let dt = EspenakMeeusDeltaT;
        let eval = InstantaneousEvaluator::new(
            &eph,
            &dt,
            R_SUN,
            R_MOON,
            AstrometryOptions::geometric(),
            interval(),
        );
        let e = eval.at(tt_j2000()).expect("geometric Mock should succeed");
        assert!(e.gamma() < 1e-6, "clear_annular gamma = {}", e.gamma());
        assert!(e.l2 > 0.0, "clear_annular l2 = {} (金環は正)", e.l2);
        assert!((0.005..0.05).contains(&e.l2.abs()), "|l2| = {}", e.l2.abs());
    }

    /// Mock 統合（clear_partial, geometric）: 影軸が地球縁を外す部分食配置 ⇒ gamma が大きい
    /// （[1.0,1.55)）。central_total（gamma≈0）と対で gamma の大小を縛る。
    /// 殺す変異: 影軸射影の縮退・gamma の x/y 片成分落とし。
    #[test]
    fn evaluator_mock_clear_partial_geometric_large_gamma() {
        let eph = MockEphemeris::clear_partial();
        let dt = EspenakMeeusDeltaT;
        let eval = InstantaneousEvaluator::new(
            &eph,
            &dt,
            R_SUN,
            R_MOON,
            AstrometryOptions::geometric(),
            interval(),
        );
        let e = eval.at(tt_j2000()).expect("geometric Mock should succeed");
        assert!(
            (1.0..1.55).contains(&e.gamma()),
            "clear_partial gamma = {}",
            e.gamma()
        );
    }

    /// Mock 統合（shadow_misses_earth, geometric）: 影軸を地球から大きく外す ⇒ gamma>1.55。
    /// 部分食（[1.0,1.55)）より更に大きいことで配置の単調性を縛る。
    #[test]
    fn evaluator_mock_shadow_miss_geometric_very_large_gamma() {
        let eph = MockEphemeris::shadow_misses_earth();
        let dt = EspenakMeeusDeltaT;
        let eval = InstantaneousEvaluator::new(
            &eph,
            &dt,
            R_SUN,
            R_MOON,
            AstrometryOptions::geometric(),
            interval(),
        );
        let e = eval.at(tt_j2000()).expect("geometric Mock should succeed");
        assert!(
            e.gamma() > 1.55,
            "shadow_misses_earth gamma = {}",
            e.gamma()
        );
    }

    /// fit_interval 透過: 構築時 interval を start/end とも そのまま返す。
    /// 殺す変異: 区間の取り違え・start/end 入れ替え・別区間の生成。
    #[test]
    fn evaluator_fit_interval_returns_constructed_interval() {
        let eph = MockEphemeris::central_total();
        let dt = EspenakMeeusDeltaT;
        let iv = interval();
        let eval = InstantaneousEvaluator::new(
            &eph,
            &dt,
            R_SUN,
            R_MOON,
            AstrometryOptions::geometric(),
            iv,
        );
        let got = eval.fit_interval();
        assert_eq!(got, iv, "fit_interval should echo the constructed interval");
        assert_eq!(got.start, iv.start, "fit_interval start");
        assert_eq!(got.end, iv.end, "fit_interval end");
    }

    /// エラー透過: Mock は velocity None。standard（light_time/aberration ON）ゆえ apparent_cirs が
    /// `EphemerisError::DataUnavailable` を返し、評価器の at() は `?` で
    /// `EclipseError::Ephemeris(EphemerisError::DataUnavailable)` に透過ラップして返す。
    /// 殺す変異: エラーの握り潰し・unwrap 化・別 variant への誤分類。
    #[test]
    fn evaluator_mock_standard_propagates_data_unavailable() {
        let eph = MockEphemeris::central_total();
        let dt = EspenakMeeusDeltaT;
        let eval = InstantaneousEvaluator::new(
            &eph,
            &dt,
            R_SUN,
            R_MOON,
            AstrometryOptions::standard(),
            interval(),
        );
        let r = eval.at(tt_j2000());
        assert!(
            matches!(
                r,
                Err(EclipseError::Ephemeris(EphemerisError::DataUnavailable))
            ),
            "expected Err(Ephemeris(DataUnavailable)), got {r:?}"
        );
    }

    /// BesselianSource trait 経由: `&dyn BesselianSource` 越しに at()/fit_interval() を呼べる
    /// （DirectBesselianSource と差し替え可能・object-safe）。at() は具象呼びと厳密一致。
    /// 殺す変異: trait 実装の欠落・object-safety 違反・dyn 越しの値ずれ。
    #[test]
    fn evaluator_usable_through_dyn_trait_object() {
        let eph = AnalyticalEphemeris::new();
        let dt = EspenakMeeusDeltaT;
        let iv = interval();
        let concrete = InstantaneousEvaluator::new(
            &eph,
            &dt,
            R_SUN,
            R_MOON,
            AstrometryOptions::standard(),
            iv,
        );
        let s: &dyn BesselianSource = &concrete;

        // fit_interval() を trait オブジェクト越しに。
        assert_eq!(s.fit_interval(), iv, "fit_interval via &dyn");

        // at() を trait オブジェクト越しに。besselian_elements_at と厳密一致。
        let t = tt_2017_max();
        let got = s.at(t).expect("at() via &dyn should succeed");
        let want = besselian_elements_at(t, R_SUN, R_MOON, &dt).expect("oracle should succeed");
        assert_exact_eq(
            &got,
            &want,
            "evaluator at() via &dyn == besselian_elements_at",
        );
    }

    /// options 影響: 同一 Analytical 暦で standard と geometric の at() は結果が変わる
    /// （見かけ補正の有無が x/y/mu に効く）。geometric は velocity 不要ゆえ Analytical でも Ok。
    /// 殺す変異: options 引数の無視・常に standard/geometric を使うハードコード。
    #[test]
    fn evaluator_options_change_result_on_analytical() {
        let eph = AnalyticalEphemeris::new();
        let dt = EspenakMeeusDeltaT;
        let t = tt_2017_max();
        let std_eval = InstantaneousEvaluator::new(
            &eph,
            &dt,
            R_SUN,
            R_MOON,
            AstrometryOptions::standard(),
            interval(),
        );
        let geo_eval = InstantaneousEvaluator::new(
            &eph,
            &dt,
            R_SUN,
            R_MOON,
            AstrometryOptions::geometric(),
            interval(),
        );
        let std = std_eval.at(t).expect("standard at() should succeed");
        let geo = geo_eval
            .at(t)
            .expect("geometric at() should succeed (no velocity needed)");
        // 見かけ補正（光行時間＋光行差）の有無で少なくとも 1 フィールドは有意に変わる。
        let differs = (std.x - geo.x).abs() > 1e-9
            || (std.y - geo.y).abs() > 1e-9
            || (std.mu.0 - geo.mu.0).abs() > 1e-9;
        assert!(
            differs,
            "standard vs geometric should differ: std=({},{},{}) geo=({},{},{})",
            std.x, std.y, std.mu.0, geo.x, geo.y, geo.mu.0
        );
    }
}
