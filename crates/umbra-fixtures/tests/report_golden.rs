//! ISSUE-030 S30d 受け入れテスト（strict / `report_against_golden` オーケストレーション）。
//!
//! 本ファイルは `umbra-fixtures` の **公開 API のみ**を対象とした統合テスト（tests/ 配下）。
//! 対象は本スライスの **エンドツーエンドのゴールデン照合オーケストレーション**:
//! - `trait GoldenComputer`（食計算・局地条件計算を注入する純インターフェース）。
//! - `report_against_golden(computer, golden, profile) -> Result<GoldenReport, EclipseError>`
//!   （S30b/S30c の純比較・集計を end-to-end に束ねる純オーケストレーション）。
//! - `struct GoldenReport`（全球＋地点別レポート＋発見/取りこぼし/比較地点の件数）。
//!
//! ## テスト戦略（mutation-resistant）
//! オーケストレーションを **遅いエンジンを走らせずに** 縛るため、計算能力を `GoldenComputer`
//! trait で注入する。FAST テストは結果を完全制御する **モック** computer を使い、件数・配線・
//! エラー伝播・空入力を独立に固定する。エンジン実走は **SLOW テスト 1 件のみ**（実エンジン）。
//!
//! モックが返す `SolarEclipse` / `LocalCircumstances` は computed == golden になるよう構築し
//! （誤差が予測可能）、または既知の固定誤差を仕込んで集計統計を縛る。
//!
//! ## 期待される RED（実装前）
//! `GoldenComputer` / `GoldenReport` / `report_against_golden` はまだ存在しないため、
//! 本ファイルは **未解決インポート（E0432）でコンパイル不能 = RED** になる。これが想定どおりの赤。

#![allow(clippy::excessive_precision)]

use umbra_core::{Degrees, JulianDate2, Kilometers, Observer, TimeInterval, TtInstant, UtcInstant};
use umbra_eclipse::{
    AccuracyProfile, BesselFitError, BesselianPolynomial, CalculationMetadata, EclipseError,
    EclipseMagnitude, GlobalCircumstances, GlobalContact, GreatestEclipse, LocalCircumstances,
    LocalContact, LocalContactSet, Obscuration, Polynomial, SolarEclipse, SolarEclipseKind,
    UtcRange, Visibility,
};
use umbra_geo::GeoPoint;

use umbra_fixtures::{
    golden_eclipses, report_against_golden, GoldenComputer, GoldenContact, GoldenEclipse,
    GoldenLocation, GoldenReport, LocationClass, OracleSource, ToleranceProfile,
};

/// 計算統計の f64 一致に用いる厳密許容。
const EPS: f64 = 1e-9;

// ============================================================
// 構築ヘルパ（report_global.rs / report_local.rs / results.rs のパターンをミラー）
// ============================================================

/// UTC 瞬時を整数引数で組む。
fn utc(year: i32, month: u8, day: u8, hour: u8, minute: u8, second: f64) -> UtcInstant {
    UtcInstant::from_gregorian(year, month, day, hour, minute, second).expect("有効な UTC 日時")
}

/// TT 瞬時を 2 要素 JD で組む。
fn tt(jd1: f64, jd2: f64) -> TtInstant {
    TtInstant::from_jd2(JulianDate2::new(jd1, jd2))
}

/// 地表点（lat, lon）を度から組む。
fn geo(lat: f64, lon: f64) -> GeoPoint {
    GeoPoint::from_degrees(lat, lon).expect("有効な地表点")
}

/// 固定の最小 BesselianPolynomial（フィラー）。
fn minimal_bessel() -> BesselianPolynomial {
    let c = |v: f64| Polynomial {
        coefficients: vec![v],
    };
    BesselianPolynomial {
        epoch_tt: tt(2_451_545.0, 0.0),
        x: c(0.20),
        y: c(-0.30),
        d: c(0.2070),
        mu: c(1.2),
        l1: c(0.5400),
        l2: c(-0.0090),
        tan_f1: 0.004_65,
        tan_f2: 0.004_63,
        fit_interval: TimeInterval {
            start: tt(2_451_544.9, 0.0),
            end: tt(2_451_545.1, 0.0),
        },
        fit_error: BesselFitError {
            max_x: 1.0e-7,
            max_y: 2.0e-7,
            max_l1: 3.0e-7,
            max_l2: 4.0e-7,
        },
    }
}

/// 固定のメタデータ（フィラー）。
fn metadata() -> CalculationMetadata {
    CalculationMetadata {
        library_version: "0.1.0".to_string(),
        ephemeris_model: "ELP/MPP02+VSOP87D".to_string(),
        ephemeris_version: "2024a".to_string(),
        delta_t_model: "EspenakMeeus".to_string(),
        delta_t_uncertainty_seconds: 0.5,
        earth_model: "WGS84".to_string(),
        lunar_radius_model: "IauMean".to_string(),
        accuracy_profile: AccuracyProfile::Standard,
        generated_at: utc(2026, 6, 18, 0, 0, 0.0),
    }
}

/// 固定の全球接触点（フィラー）。
fn global_contact() -> GlobalContact {
    GlobalContact {
        time_utc: utc(2024, 4, 8, 16, 0, 0.0),
        time_tt: tt(2_460_409.0, 0.01),
        position: geo(30.0, -100.0),
    }
}

/// 合成 `computed: SolarEclipse`。compare_global が読む 4 値（greatest.time_tt /
/// greatest.time_utc / greatest.magnitude / global.gamma）を引数で指定し、残りは固定フィラー。
fn computed_eclipse(
    greatest_tt: TtInstant,
    greatest_utc: UtcInstant,
    gamma: f64,
    magnitude: f64,
) -> SolarEclipse {
    let greatest = GreatestEclipse {
        time_utc: greatest_utc,
        time_tt: greatest_tt,
        position: geo(25.0, -104.0),
        magnitude: EclipseMagnitude(magnitude),
        obscuration: Obscuration(1.0),
        path_width: Some(Kilometers(197.0)),
        central_duration: Some(268.0),
        sun_altitude: Degrees(70.3),
    };
    let global = GlobalCircumstances {
        kind: SolarEclipseKind::Total,
        partial_begin: Some(global_contact()),
        central_begin: Some(global_contact()),
        greatest,
        central_end: Some(global_contact()),
        partial_end: Some(global_contact()),
        gamma,
    };
    SolarEclipse {
        event_key: "computed".to_string(),
        kind: SolarEclipseKind::Total,
        global,
        bessel: minimal_bessel(),
        metadata: metadata(),
    }
}

/// 計算側の局地接触（時刻 2 値以外はフィラー）。
fn local_contact(time_utc: UtcInstant, time_tt: TtInstant) -> LocalContact {
    LocalContact {
        time_utc,
        time_tt,
        sun_altitude: Degrees(40.0),
        sun_azimuth: Degrees(200.0),
        position_angle: Degrees(300.0),
        visible: true,
    }
}

/// 合成 `computed: LocalCircumstances`。compare_local が読む値を引数で指定し、metadata はフィラー。
fn computed_local(
    contacts: LocalContactSet,
    magnitude: f64,
    obscuration: f64,
    max_alt: f64,
    visibility: Visibility,
) -> LocalCircumstances {
    LocalCircumstances {
        contacts,
        magnitude: EclipseMagnitude(magnitude),
        obscuration: Obscuration(obscuration),
        maximum_altitude: Degrees(max_alt),
        visibility,
        metadata: metadata(),
    }
}

/// ゴールデン接触（UTC・任意 TT・高度）。
fn golden_contact(
    time_utc: UtcInstant,
    time_tt: Option<TtInstant>,
    altitude_deg: f64,
) -> GoldenContact {
    GoldenContact {
        time_utc,
        time_tt,
        altitude_deg,
    }
}

/// 合成 `golden: GoldenLocation`（compare_local が読む値を引数で指定、残りはフィラー）。
#[allow(clippy::too_many_arguments)]
fn golden_location(
    c1: Option<GoldenContact>,
    c2: Option<GoldenContact>,
    maximum: GoldenContact,
    c3: Option<GoldenContact>,
    c4: Option<GoldenContact>,
    magnitude: f64,
    obscuration: f64,
    max_altitude_deg: f64,
    visibility_expected: Visibility,
) -> GoldenLocation {
    GoldenLocation {
        name: "x".to_string(),
        latitude_deg: 35.0,
        east_longitude_deg: 139.0,
        elevation_m: 0.0,
        location_class: LocationClass::Centerline,
        c1,
        c2,
        maximum,
        c3,
        c4,
        magnitude,
        obscuration,
        max_altitude_deg,
        max_azimuth_deg: 200.0,
        visibility_expected,
    }
}

/// 合成 `golden: GoldenEclipse`。`event_key` を引数化（モックの分岐キー）、全球 4 値と locations を指定。
fn golden_eclipse(
    event_key: &str,
    greatest_utc: UtcInstant,
    greatest_tt: Option<TtInstant>,
    gamma: f64,
    magnitude: f64,
    locations: Vec<GoldenLocation>,
) -> GoldenEclipse {
    GoldenEclipse {
        event_key: event_key.to_string(),
        kind_expected: SolarEclipseKind::Total,
        greatest_time_utc: greatest_utc,
        greatest_time_tt: greatest_tt,
        gamma,
        magnitude,
        delta_t_seconds: None,
        locations,
        source: OracleSource {
            name: "".into(),
            url: "".into(),
            retrieved: "".into(),
            delta_t_convention: "".into(),
            k_convention: "".into(),
            license_note: "".into(),
        },
    }
}

// ============================================================
// FAST テスト用モック GoldenComputer
// ============================================================

/// モック computer。`event_key` の接頭辞で分岐:
/// - `"found"`: `eclipse_on` が `Ok(Some(eclipse))`。`eclipse` は `eclipse_for` が返す合成値。
/// - `"missing"`: `eclipse_on` が `Ok(None)`（取りこぼし）。
/// - `"err"`: `eclipse_on` が `Err(NotImplemented)`。
///
/// `local_at` は通常 `Ok(local)`（`local_for` の合成値）を返すが、`local_err == true` のときは
/// `Err(DegenerateGeometry)` を返す（local_at エラー伝播テスト用）。
///
/// `eclipse_offset_seconds` / `local_*` で computed-vs-golden の既知誤差を制御する。
struct MockComputer {
    /// `eclipse_for` が返す greatest 時刻に足す既知秒差（computed − golden を作る）。
    eclipse_offset_seconds: f64,
    /// `local_at` が `Err` を返すか（local_at エラー伝播テスト用）。
    local_err: bool,
}

impl MockComputer {
    /// 既知誤差なし（computed == golden になるよう golden から合成）・local 正常。
    fn aligned() -> Self {
        MockComputer {
            eclipse_offset_seconds: 0.0,
            local_err: false,
        }
    }

    /// golden の全球値から「offset 秒ずらした」computed SolarEclipse を作る。
    /// greatest TT は golden の greatest_time_tt（あれば）に offset 秒を足す。
    /// gamma / magnitude は golden 一致（誤差 0）にして、誤差を greatest_seconds に単離する。
    fn eclipse_for(&self, golden: &GoldenEclipse) -> SolarEclipse {
        let base_tt = golden
            .greatest_time_tt
            .unwrap_or_else(|| tt(2_451_545.0, 0.0));
        let off_jd = self.eclipse_offset_seconds / 86_400.0;
        // offset は part2 に足して 2 要素を保つ（単一 f64 の jd()+off は JD≈2.45e6 で桁落ちし、
        // days_since で 1.0s が 0.99999..s になるため。julian §2要素表現）。
        let base = base_tt.jd2();
        let shifted_tt = TtInstant::from_jd2(JulianDate2::new(base.part1, base.part2 + off_jd));
        computed_eclipse(
            shifted_tt,
            golden.greatest_time_utc,
            golden.gamma,
            golden.magnitude,
        )
    }

    /// golden 地点の値に一致する computed LocalCircumstances（誤差 0）を作る。
    /// 各接触は golden の maximum 時刻に揃え（両 Some の接触のみ contact 誤差 0）、
    /// magnitude / obscuration / max_altitude / visibility も golden 一致にする。
    fn local_for(&self, loc: &GoldenLocation) -> LocalCircumstances {
        // golden の maximum TT（なければフィラー）に揃えた computed maximum。
        let max_tt = loc.maximum.time_tt.unwrap_or_else(|| tt(2_451_545.0, 0.0));
        let mk = |gc: &Option<GoldenContact>| -> Option<LocalContact> {
            gc.as_ref().map(|g| {
                let t_tt = g.time_tt.unwrap_or_else(|| tt(2_451_545.0, 0.0));
                local_contact(g.time_utc, t_tt)
            })
        };
        let contacts = LocalContactSet {
            c1: mk(&loc.c1),
            c2: mk(&loc.c2),
            maximum: local_contact(loc.maximum.time_utc, max_tt),
            c3: mk(&loc.c3),
            c4: mk(&loc.c4),
        };
        computed_local(
            contacts,
            loc.magnitude,
            loc.obscuration,
            loc.max_altitude_deg,
            loc.visibility_expected,
        )
    }
}

impl GoldenComputer for MockComputer {
    fn eclipse_on(&self, golden: &GoldenEclipse) -> Result<Option<SolarEclipse>, EclipseError> {
        if golden.event_key.starts_with("err") {
            Err(EclipseError::NotImplemented)
        } else if golden.event_key.starts_with("missing") {
            Ok(None)
        } else {
            // "found"（およびその他）→ Some。
            Ok(Some(self.eclipse_for(golden)))
        }
    }

    fn local_at(
        &self,
        _eclipse: &SolarEclipse,
        location: &GoldenLocation,
    ) -> Result<LocalCircumstances, EclipseError> {
        if self.local_err {
            Err(EclipseError::DegenerateGeometry)
        } else {
            Ok(self.local_for(location))
        }
    }
}

/// 既定の golden 地点（全 Some 接触・全部 golden 一致用に TT 付与）。
fn loc_with_n_contacts() -> GoldenLocation {
    let g_utc = utc(2024, 4, 8, 18, 0, 0.0);
    let base = tt(2_451_545.0, 0.0);
    let gc = |alt: f64| golden_contact(g_utc, Some(base), alt);
    golden_location(
        Some(gc(10.0)),
        Some(gc(20.0)),
        gc(50.0),
        Some(gc(30.0)),
        Some(gc(40.0)),
        1.0,
        1.0,
        50.0,
        Visibility::FullyVisible,
    )
}

// ============================================================
// FAST テスト（MockComputer 経由・mutation-resistant）
// ============================================================

/// 受け入れ「`eclipses_found`/`eclipses_missing` は computer の Some/None 数で数える」。
/// golden 3 件: 2 件 Some・1 件 None → found==2, missing==1。
/// 殺す変異: found/missing の数え間違い・入れ替え、missing を増分しない。
#[test]
fn report_against_golden_counts_found_and_missing() {
    let g_utc = utc(2024, 4, 8, 18, 0, 0.0);
    let base = tt(2_451_545.0, 0.0);
    let golden = vec![
        golden_eclipse("found-a", g_utc, Some(base), 0.4, 1.03, vec![]),
        golden_eclipse("missing-b", g_utc, Some(base), 0.4, 1.03, vec![]),
        golden_eclipse("found-c", g_utc, Some(base), 0.4, 1.03, vec![]),
    ];
    let report: GoldenReport = report_against_golden(
        &MockComputer::aligned(),
        &golden,
        &ToleranceProfile::standard(),
    )
    .expect("モックはエラーを返さない");
    assert_eq!(report.eclipses_found, 2, "Some を返した golden 数 = 2");
    assert_eq!(report.eclipses_missing, 1, "None を返した golden 数 = 1");
}

/// 受け入れ「`locations_compared` は **found** golden の地点総数のみ（missing の地点は数えない）」。
/// found-a に 2 地点・found-c に 3 地点・missing-b にも 2 地点（数えてはならない）→ 5。
/// 殺す変異: missing の地点を数える、ループ誤り、オフ・バイ・ワン。
#[test]
fn report_against_golden_counts_locations_across_found() {
    let g_utc = utc(2024, 4, 8, 18, 0, 0.0);
    let base = tt(2_451_545.0, 0.0);
    let loc = loc_with_n_contacts();
    let golden = vec![
        golden_eclipse(
            "found-a",
            g_utc,
            Some(base),
            0.4,
            1.03,
            vec![loc.clone(), loc.clone()], // 2 地点
        ),
        golden_eclipse(
            "missing-b",
            g_utc,
            Some(base),
            0.4,
            1.03,
            vec![loc.clone(), loc.clone()], // 2 地点（数えてはならない）
        ),
        golden_eclipse(
            "found-c",
            g_utc,
            Some(base),
            0.4,
            1.03,
            vec![loc.clone(), loc.clone(), loc.clone()], // 3 地点
        ),
    ];
    let report = report_against_golden(
        &MockComputer::aligned(),
        &golden,
        &ToleranceProfile::standard(),
    )
    .expect("モックはエラーを返さない");
    assert_eq!(report.eclipses_found, 2, "found は 2 件");
    assert_eq!(
        report.locations_compared, 5,
        "found 地点総数 = 2 + 3 = 5（missing の 2 地点は数えない）"
    );
}

/// 受け入れ「compare 結果が aggregate に正しく配線され、全球/地点別が混ざらない」。
/// モックに既知の全球誤差 +1.0 s を仕込み（gamma/magnitude は誤差 0）、found 2 件。
/// → `global.greatest.n == eclipses_found(==2)` かつ `global.greatest.max_abs ≈ 1.0`。
/// 各 found に 2 地点（全 golden 一致＝誤差 0）→ `local.maximum.n == locations_compared(==4)`。
/// 殺す変異: compare 結果を aggregate に渡さない、global/local の取り違え、誤った vec 配線。
#[test]
fn report_against_golden_aggregates_global_and_local() {
    let g_utc = utc(2024, 4, 8, 18, 0, 0.0);
    let base = tt(2_451_545.0, 0.0);
    let loc = loc_with_n_contacts();
    let golden = vec![
        golden_eclipse(
            "found-a",
            g_utc,
            Some(base),
            0.4,
            1.03,
            vec![loc.clone(), loc.clone()],
        ),
        golden_eclipse(
            "found-b",
            g_utc,
            Some(base),
            0.4,
            1.03,
            vec![loc.clone(), loc.clone()],
        ),
    ];
    // 既知全球誤差 +1.0 s（gamma/magnitude は一致＝誤差 0）。local は一致（誤差 0）。
    let computer = MockComputer {
        eclipse_offset_seconds: 1.0,
        local_err: false,
    };
    let report = report_against_golden(&computer, &golden, &ToleranceProfile::standard())
        .expect("モックはエラーを返さない");

    assert_eq!(report.eclipses_found, 2, "found は 2 件");
    assert_eq!(report.locations_compared, 4, "found 地点総数 = 2 + 2 = 4");

    // 全球: 各 found の greatest_seconds 誤差 ≈ +1.0 が n=2 件集計される。
    assert_eq!(
        report.global.greatest.n, 2,
        "global.greatest.n は found 数（2）= 全球誤差を集計に渡している"
    );
    assert!(
        (report.global.greatest.max_abs - 1.0).abs() < EPS,
        "global.greatest.max_abs ≈ 1.0（仕込んだ既知誤差）, got {}",
        report.global.greatest.max_abs
    );

    // 地点別: maximum 接触は 1 地点 1 件 → n = locations_compared(4)。誤差 0。
    assert_eq!(
        report.local.maximum.n, 4,
        "local.maximum.n は比較地点数（4）= 地点別誤差を集計に渡している"
    );
}

/// 受け入れ「`eclipse_on` の `Err` は `?` で関数全体に伝播する」。
/// golden 3 件のうち 1 件が "err"（`eclipse_on` が Err）→ `report_against_golden` は Err。
/// 殺す変異: エラーを握り潰す、unwrap、err を None 扱いにする。
#[test]
fn report_against_golden_propagates_eclipse_on_error() {
    let g_utc = utc(2024, 4, 8, 18, 0, 0.0);
    let base = tt(2_451_545.0, 0.0);
    let golden = vec![
        golden_eclipse("found-a", g_utc, Some(base), 0.4, 1.03, vec![]),
        golden_eclipse("err-b", g_utc, Some(base), 0.4, 1.03, vec![]),
        golden_eclipse("found-c", g_utc, Some(base), 0.4, 1.03, vec![]),
    ];
    let r = report_against_golden(
        &MockComputer::aligned(),
        &golden,
        &ToleranceProfile::standard(),
    );
    assert!(
        r.is_err(),
        "eclipse_on の Err が伝播して全体が Err になる, got {r:?}"
    );
}

/// 受け入れ「`local_at` の `Err` も `?` で関数全体に伝播する」。
/// found golden に 1 地点。computer は `eclipse_on` で Ok(Some) を返すが `local_at` で Err。
/// → `report_against_golden` は Err。
/// 殺す変異: local_at のエラーを無視する、unwrap、地点エラーを握り潰す。
#[test]
fn report_against_golden_propagates_local_at_error() {
    let g_utc = utc(2024, 4, 8, 18, 0, 0.0);
    let base = tt(2_451_545.0, 0.0);
    let golden = vec![golden_eclipse(
        "found-a",
        g_utc,
        Some(base),
        0.4,
        1.03,
        vec![loc_with_n_contacts()],
    )];
    let computer = MockComputer {
        eclipse_offset_seconds: 0.0,
        local_err: true, // local_at が Err を返す
    };
    let r = report_against_golden(&computer, &golden, &ToleranceProfile::standard());
    assert!(
        r.is_err(),
        "local_at の Err が伝播して全体が Err になる, got {r:?}"
    );
}

/// 受け入れ「空 golden は空レポート（全件数 0・空集計の vacuous pass）」。
/// 空スライス → found==0, missing==0, locations_compared==0、global/local は vacuous pass==true。
/// 殺す変異: 空でパニック、空で件数が 0 にならない、空で pass==false。
#[test]
fn report_against_golden_empty_is_empty_report() {
    let golden: Vec<GoldenEclipse> = vec![];
    let report = report_against_golden(
        &MockComputer::aligned(),
        &golden,
        &ToleranceProfile::standard(),
    )
    .expect("空入力はエラーにしない");
    assert_eq!(report.eclipses_found, 0, "空 → found==0");
    assert_eq!(report.eclipses_missing, 0, "空 → missing==0");
    assert_eq!(report.locations_compared, 0, "空 → locations_compared==0");
    assert!(report.global.pass, "空 → global は vacuous pass==true");
    assert!(report.local.pass, "空 → local は vacuous pass==true");
    assert_eq!(report.global.greatest.n, 0, "空 → global.greatest.n==0");
    assert_eq!(report.local.maximum.n, 0, "空 → local.maximum.n==0");
}

/// 受け入れ「missing(None) golden は global/local 集計に **何も寄与しない**（地点があっても）」。
/// found-a（地点 0）と missing-b（地点 3・数えてはならない）→ global.greatest.n==1（found のみ）、
/// local.maximum.n==0（missing の地点は集計しない）、locations_compared==0。
/// 殺す変異: missing の地点を処理する、missing を集計 vec に積む。
#[test]
fn report_against_golden_missing_does_not_add_errors() {
    let g_utc = utc(2024, 4, 8, 18, 0, 0.0);
    let base = tt(2_451_545.0, 0.0);
    let loc = loc_with_n_contacts();
    let golden = vec![
        golden_eclipse("found-a", g_utc, Some(base), 0.4, 1.03, vec![]),
        golden_eclipse(
            "missing-b",
            g_utc,
            Some(base),
            0.4,
            1.03,
            vec![loc.clone(), loc.clone(), loc.clone()], // 3 地点（無視されるべき）
        ),
    ];
    let report = report_against_golden(
        &MockComputer::aligned(),
        &golden,
        &ToleranceProfile::standard(),
    )
    .expect("モックはエラーを返さない");
    assert_eq!(report.eclipses_found, 1, "found は 1 件");
    assert_eq!(report.eclipses_missing, 1, "missing は 1 件");
    assert_eq!(
        report.locations_compared, 0,
        "found の地点は 0、missing の 3 地点は数えない → 0"
    );
    assert_eq!(
        report.global.greatest.n, 1,
        "global.greatest.n は found のみ（1）"
    );
    assert_eq!(
        report.local.maximum.n, 0,
        "local.maximum.n は found 地点のみ（0）= missing の地点を処理しない"
    );
}

// ============================================================
// SLOW 統合テスト（実エンジン・1 件のみ）
// ============================================================

/// 実エンジンを包む `GoldenComputer`（SLOW テスト専用）。
struct EngineComputer {
    engine: umbra_eclipse::StandardEngine,
}

impl GoldenComputer for EngineComputer {
    fn eclipse_on(&self, golden: &GoldenEclipse) -> Result<Option<SolarEclipse>, EclipseError> {
        // golden の最大食 UTC を中心に ±0.5 日窓を探索し、先頭の食を採る。
        let jd = golden.greatest_time_utc.jd2().jd();
        let start = UtcInstant::from_jd2(JulianDate2::from_jd(jd - 0.5));
        let end = UtcInstant::from_jd2(JulianDate2::from_jd(jd + 0.5));
        Ok(self
            .engine
            .search(UtcRange { start, end })?
            .into_iter()
            .next())
    }

    fn local_at(
        &self,
        eclipse: &SolarEclipse,
        loc: &GoldenLocation,
    ) -> Result<LocalCircumstances, EclipseError> {
        let obs =
            Observer::from_degrees(loc.latitude_deg, loc.east_longitude_deg, loc.elevation_m)?;
        self.engine.local_circumstances(eclipse, obs)
    }
}

/// 【SLOW・1 件】実エンジンで 1 ゴールデンを end-to-end 照合し、**構造**（件数・統計の populated）
/// を縛る。`pass` は実エンジン-vs-オラクルの本来の検証結果（false でも正当）ゆえ **assert しない**。
/// 実 `search`（数分）＋数回の `local_circumstances` を 1 ゴールデンに限定して実行時間を抑える。
/// 殺す変異: found/missing/locations の数え違い、global/local 集計の未配線。
// SLOW
#[test]
fn report_against_golden_real_engine_one_golden() {
    let computer = EngineComputer {
        engine: umbra_eclipse::standard_engine(umbra_ephemeris::bundled_time_data()),
    };
    let goldens = golden_eclipses();
    // 2017-08-21-total を優先・無ければ先頭。
    let target = goldens
        .iter()
        .find(|g| g.event_key == "2017-08-21-total")
        .or_else(|| goldens.first())
        .expect("ゴールデンが 1 件以上ある");

    let report = report_against_golden(
        &computer,
        std::slice::from_ref(target),
        &ToleranceProfile::standard(),
    )
    .expect("実エンジン照合は Ok を返す");

    assert_eq!(report.eclipses_found, 1, "1 ゴールデン → found==1");
    assert_eq!(report.eclipses_missing, 0, "取りこぼしなし → missing==0");
    assert_eq!(
        report.locations_compared,
        target.locations.len(),
        "比較地点数 = ゴールデンの地点数"
    );
    // 統計が populated（pass の真偽は問わない）。
    assert_eq!(report.global.greatest.n, 1, "global.greatest.n==1");
    assert_eq!(
        report.local.maximum.n,
        target.locations.len(),
        "local.maximum.n = 地点数"
    );
}
