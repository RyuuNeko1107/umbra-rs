//! ISSUE-030 S30f 受け入れテスト（strict / `cargo xtask validate` ゴールデン照合コマンド）。
//!
//! 本ファイルは `xtask` crate の **公開 API のみ**を対象とした統合テスト（tests/ 配下）。
//! 対象は本スライスが導入する `xtask::validate` モジュール:
//! - `parse_format` / `parse_accuracy`（純引数解釈・既定値・未知/欠落エラー）。
//! - `validate_report`（注入 `GoldenComputer` でゴールデン照合し text/json 文字列を返す純粋寄り関数）。
//! - `EngineGoldenComputer`（実エンジンを包む `GoldenComputer`。SLOW テスト 1 件のみで縛る）。
//!
//! ## テスト戦略（mutation-resistant）
//! `validate_report` の配線（照合→集計→レンダ）を **遅いエンジンを走らせずに** 縛るため、計算能力を
//! `umbra_fixtures::GoldenComputer` trait で注入する。FAST テストは結果を完全制御する **モック**
//! computer を使い、レンダ種別（text/json）・件数の通し・エラー伝播を独立に固定する。実エンジン経路は
//! **SLOW テスト 1 件のみ**（`EngineGoldenComputer` × 1 ゴールデン）。
//!
//! 合成 `SolarEclipse` / `LocalCircumstances` / `GoldenEclipse` 等の構築ヘルパは
//! `crates/umbra-fixtures/tests/report_golden.rs` のパターンをミラーする（computed == golden 一致）。
//!
//! ## 期待される RED（実装前）
//! `xtask::validate` モジュールも `parse_format` / `parse_accuracy` / `validate_report` /
//! `EngineGoldenComputer` も存在しないため、本ファイルは **未解決インポート/モジュール（E0432/E0433）
//! でコンパイル不能 = RED** になる。これが想定どおりの赤。

#![allow(clippy::excessive_precision)]

use serde_json::Value;

use umbra_core::{Degrees, JulianDate2, Kilometers, TimeInterval, TtInstant, UtcInstant};
use umbra_eclipse::{
    AccuracyProfile, BesselFitError, BesselianPolynomial, CalculationMetadata, EclipseError,
    EclipseMagnitude, GlobalCircumstances, GlobalContact, GreatestEclipse, LocalCircumstances,
    LocalContact, LocalContactSet, Obscuration, Polynomial, SolarEclipse, SolarEclipseKind,
    Visibility,
};
use umbra_geo::GeoPoint;

use umbra_fixtures::{
    golden_eclipses, GoldenComputer, GoldenContact, GoldenEclipse, GoldenLocation, LocationClass,
    OracleSource, ToleranceProfile,
};

use xtask::validate::{
    parse_accuracy, parse_format, validate_report, AccuracyArg, EngineGoldenComputer, ReportFormat,
};

// ============================================================
// 引数列ヘルパ
// ============================================================

/// &str スライス → Vec<String>（parse_format/parse_accuracy の引数列構築ヘルパ）。
fn args(parts: &[&str]) -> Vec<String> {
    parts.iter().map(|s| s.to_string()).collect()
}

// ============================================================
// 構築ヘルパ（report_golden.rs のパターンをミラー）
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

/// 合成 `computed: SolarEclipse`（compare_global が読む 4 値を引数、残りはフィラー）。
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

/// 合成 `computed: LocalCircumstances`（compare_local が読む値を引数、metadata はフィラー）。
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

/// 合成 `golden: GoldenLocation`（compare_local が読む値を引数、残りはフィラー）。
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

/// 合成 `golden: GoldenEclipse`。`event_key` を引数化（モックの分岐キー）。
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
// FAST テスト用モック GoldenComputer
// ============================================================

/// モック computer。`event_key` の接頭辞で分岐:
/// - `"found"`: `eclipse_on` → `Ok(Some(eclipse))`（golden から合成・誤差 0）。
/// - `"missing"`: `eclipse_on` → `Ok(None)`（取りこぼし）。
/// - `"err"`: `eclipse_on` → `Err(NotImplemented)`（エラー伝播テスト用）。
///
/// `local_at` は常に `Ok(local)`（golden 一致の合成値）を返す。
struct MockComputer;

impl MockComputer {
    /// golden の全球値に一致する computed SolarEclipse（誤差 0）を作る。
    fn eclipse_for(golden: &GoldenEclipse) -> SolarEclipse {
        let base_tt = golden
            .greatest_time_tt
            .unwrap_or_else(|| tt(2_451_545.0, 0.0));
        computed_eclipse(
            base_tt,
            golden.greatest_time_utc,
            golden.gamma,
            golden.magnitude,
        )
    }

    /// golden 地点の値に一致する computed LocalCircumstances（誤差 0）を作る。
    fn local_for(loc: &GoldenLocation) -> LocalCircumstances {
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
            Ok(Some(Self::eclipse_for(golden)))
        }
    }

    fn local_at(
        &self,
        _eclipse: &SolarEclipse,
        location: &GoldenLocation,
    ) -> Result<LocalCircumstances, EclipseError> {
        Ok(Self::local_for(location))
    }
}

// ============================================================
// parse_format（FAST・純引数解釈）
// ============================================================

/// 受け入れ「`--format` 既定 Text・text→Text・json→Json」。
/// 殺す変異: 既定値の取り違え、text/json のマッピング入れ替え、フラグ名比較の改変。
#[test]
fn parse_format_defaults_and_values() {
    assert_eq!(
        parse_format(&[]).expect("既定は Ok(Text)"),
        ReportFormat::Text,
        "引数なし → 既定 Text"
    );
    assert_eq!(
        parse_format(&args(&["--format", "text"])).expect("text は Ok"),
        ReportFormat::Text,
        "text → Text"
    );
    assert_eq!(
        parse_format(&args(&["--format", "json"])).expect("json は Ok"),
        ReportFormat::Json,
        "json → Json"
    );
}

/// 受け入れ「未知値は Err・値欠落は Err(MissingArgument)」。
/// 殺す変異: 未知値を既定にフォールバックする、欠落を黙って既定にする、エラー種別の取り違え。
#[test]
fn parse_format_unknown_and_missing_are_errors() {
    let unknown = parse_format(&args(&["--format", "xml"]));
    assert!(
        matches!(
            unknown,
            Err(xtask::error::XtaskError::InvalidArgument { ref flag, ref value })
                if flag == "--format" && value == "xml"
        ),
        "未知値は InvalidArgument{{flag:--format, value:xml}}, got {unknown:?}"
    );

    let missing = parse_format(&args(&["--format"])).expect_err("値欠落は Err");
    assert!(
        matches!(missing, xtask::error::XtaskError::MissingArgument(_)),
        "値欠落は MissingArgument, got {missing:?}"
    );
}

// ============================================================
// parse_accuracy（FAST・純引数解釈）
// ============================================================

/// 受け入れ「`--accuracy` 既定 Standard・standard→Standard・reference→Reference」。
/// 殺す変異: 既定値の取り違え、standard/reference のマッピング入れ替え。
#[test]
fn parse_accuracy_defaults_and_values() {
    assert_eq!(
        parse_accuracy(&[]).expect("既定は Ok(Standard)"),
        AccuracyArg::Standard,
        "引数なし → 既定 Standard"
    );
    assert_eq!(
        parse_accuracy(&args(&["--accuracy", "standard"])).expect("standard は Ok"),
        AccuracyArg::Standard,
        "standard → Standard"
    );
    assert_eq!(
        parse_accuracy(&args(&["--accuracy", "reference"])).expect("reference は Ok"),
        AccuracyArg::Reference,
        "reference → Reference"
    );
}

/// 受け入れ「未知値は Err・値欠落は Err(MissingArgument)」。
/// 殺す変異: 未知値を既定にフォールバックする、欠落を黙って既定にする、エラー種別の取り違え。
#[test]
fn parse_accuracy_unknown_and_missing_are_errors() {
    let unknown = parse_accuracy(&args(&["--accuracy", "bogus"]));
    assert!(
        matches!(
            unknown,
            Err(xtask::error::XtaskError::InvalidArgument { ref flag, ref value })
                if flag == "--accuracy" && value == "bogus"
        ),
        "未知値は InvalidArgument{{flag:--accuracy, value:bogus}}, got {unknown:?}"
    );

    let missing = parse_accuracy(&args(&["--accuracy"])).expect_err("値欠落は Err");
    assert!(
        matches!(missing, xtask::error::XtaskError::MissingArgument(_)),
        "値欠落は MissingArgument, got {missing:?}"
    );
}

// ============================================================
// validate_report（FAST・MockComputer 経由）
// ============================================================

/// 受け入れ「Text 形式は人間可読レポート（render_text 相当）を返す」。
/// `render_text` は "GLOBAL"/"LOCAL"/"found" を出すので、小文字化後に "global"/"local"/"found" を含む。
/// found 2 件の件数行（"2 found"）も含む。
/// 殺す変異: Json 経路に分岐する、render_text を呼ばない、空文字列を返す、件数を通さない。
#[test]
fn validate_report_text_contains_report() {
    let g_utc = utc(2024, 4, 8, 18, 0, 0.0);
    let base = tt(2_451_545.0, 0.0);
    let golden = vec![
        golden_eclipse("found-a", g_utc, Some(base), 0.4, 1.03, vec![]),
        golden_eclipse("found-b", g_utc, Some(base), 0.4, 1.03, vec![]),
    ];
    let s = validate_report(
        &MockComputer,
        &golden,
        &ToleranceProfile::standard(),
        ReportFormat::Text,
    )
    .expect("Text レンダは Ok");
    let lower = s.to_lowercase();
    assert!(lower.contains("global"), "Text に global マーカ: {s}");
    assert!(lower.contains("local"), "Text に local マーカ: {s}");
    assert!(lower.contains("found"), "Text に found マーカ: {s}");
    assert!(
        s.contains("2 found"),
        "Text に found 件数 2 が通っている: {s}"
    );
}

/// 受け入れ「Json 形式は 5 トップレベルキーを持つ JSON オブジェクトで、件数が正しく通る」。
/// golden を found 2・missing 1・found 側合計 3 地点で構築 →
/// eclipses_found==2, eclipses_missing==1, locations_compared==3。
/// 殺す変異: Text 経路に分岐する、render_json を呼ばない、件数を取り違える/入れ替える、
///   report_against_golden に誤データを渡す。
#[test]
fn validate_report_json_is_object_with_keys_and_counts() {
    let g_utc = utc(2024, 4, 8, 18, 0, 0.0);
    let base = tt(2_451_545.0, 0.0);
    let loc = loc_with_n_contacts();
    let golden = vec![
        // found-a: 1 地点
        golden_eclipse("found-a", g_utc, Some(base), 0.4, 1.03, vec![loc.clone()]),
        // missing-b: 2 地点（数えてはならない）
        golden_eclipse(
            "missing-b",
            g_utc,
            Some(base),
            0.4,
            1.03,
            vec![loc.clone(), loc.clone()],
        ),
        // found-c: 2 地点
        golden_eclipse(
            "found-c",
            g_utc,
            Some(base),
            0.4,
            1.03,
            vec![loc.clone(), loc.clone()],
        ),
    ];
    let s = validate_report(
        &MockComputer,
        &golden,
        &ToleranceProfile::standard(),
        ReportFormat::Json,
    )
    .expect("Json レンダは Ok");

    let v: Value = serde_json::from_str(&s).expect("出力は妥当な JSON");
    assert!(v.is_object(), "トップレベルは JSON オブジェクト: {s}");
    let obj = v.as_object().expect("オブジェクト");
    for key in [
        "global",
        "local",
        "eclipses_found",
        "eclipses_missing",
        "locations_compared",
    ] {
        assert!(obj.contains_key(key), "トップレベルキー {key} が存在: {s}");
    }
    assert_eq!(
        v["eclipses_found"], 2,
        "found==2（Some を返した golden 数）"
    );
    assert_eq!(
        v["eclipses_missing"], 1,
        "missing==1（None を返した golden 数）"
    );
    assert_eq!(
        v["locations_compared"], 3,
        "locations_compared==3（found 側 1+2、missing の 2 地点は数えない）"
    );
}

/// 受け入れ「computer の `eclipse_on` Err は `?` で validate_report 全体に伝播する」。
/// golden に "err" を 1 件含める → validate_report は Err（EclipseError → XtaskError::Eclipse）。
/// 殺す変異: エラーを握り潰す、unwrap、err を None 扱いにする、Ok を返す。
#[test]
fn validate_report_propagates_computer_error() {
    let g_utc = utc(2024, 4, 8, 18, 0, 0.0);
    let base = tt(2_451_545.0, 0.0);
    let golden = vec![
        golden_eclipse("found-a", g_utc, Some(base), 0.4, 1.03, vec![]),
        golden_eclipse("err-b", g_utc, Some(base), 0.4, 1.03, vec![]),
    ];
    let r = validate_report(
        &MockComputer,
        &golden,
        &ToleranceProfile::standard(),
        ReportFormat::Json,
    );
    assert!(
        matches!(r, Err(xtask::error::XtaskError::Eclipse(_))),
        "eclipse_on の EclipseError が XtaskError::Eclipse に写像されて伝播する, got {r:?}"
    );
}

// ============================================================
// SLOW 統合テスト（実エンジン・1 ゴールデンのみ）
// ============================================================

/// 【SLOW・1 件】`EngineGoldenComputer`（実エンジン）で 1 ゴールデンを end-to-end 照合し、Json の
/// 件数（found==1・missing==0・locations_compared==地点数）を縛る。`pass` は実エンジン-vs-オラクルの
/// 本来の検証結果（false でも正当）ゆえ **assert しない**。
/// 実 `search`（数分）＋数回の `local_circumstances` を 1 ゴールデンに限定して実行時間を抑える。
/// 殺す変異: found/missing/locations の数え違い、validate_report が件数を通さない、エンジン未配線。
// SLOW
#[test]
fn validate_report_real_engine_one_golden_json() {
    let computer = EngineGoldenComputer::new(AccuracyArg::Standard);
    let goldens = golden_eclipses();
    // 2017-08-21-total を優先・無ければ先頭。
    let target = goldens
        .iter()
        .find(|g| g.event_key == "2017-08-21-total")
        .or_else(|| goldens.first())
        .expect("ゴールデンが 1 件以上ある");

    let s = validate_report(
        &computer,
        std::slice::from_ref(target),
        &ToleranceProfile::standard(),
        ReportFormat::Json,
    )
    .expect("実エンジン照合の Json レンダは Ok");

    let v: Value = serde_json::from_str(&s).expect("出力は妥当な JSON");
    assert_eq!(v["eclipses_found"], 1, "1 ゴールデン → found==1");
    assert_eq!(v["eclipses_missing"], 0, "取りこぼしなし → missing==0");
    assert_eq!(
        v["locations_compared"],
        target.locations.len(),
        "比較地点数 = ゴールデンの地点数"
    );
}
