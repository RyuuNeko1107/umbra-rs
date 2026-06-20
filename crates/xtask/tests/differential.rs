//! ISSUE-030 DE 差分・誤差層分解ハーネス受け入れテスト（strict / `cargo xtask differential`）。
//!
//! 本ファイルは `xtask` crate の **公開 API のみ**を対象とした統合テスト（tests/ 配下）。
//! 対象は本スライスが導入する `xtask::differential` モジュール:
//! - `differential_report`（注入 2 `GoldenComputer`〔analytical 役・DE 役〕でゴールデンを層分解し
//!   text/json 文字列を返す純粋寄り関数。`report_differential` →
//!   `render_differential_text`/`render_differential_json` の結線）。
//! - `JplGoldenComputer::from_spk_path`（実 DE440s SPK を包む `GoldenComputer`。異常系 FAST・
//!   正常系 SLOW 1 件のみ）。
//!
//! ## テスト戦略（mutation-resistant・負荷配分）
//! `differential_report` の配線（層分解→レンダ種別分岐→エラー伝播）を **遅い実エンジン/DE を走らせず
//! に** 縛るため、計算能力を `umbra_fixtures::GoldenComputer` trait で 2 系統注入する。FAST テストは
//! 結果を完全制御する **モック** computer（analytical 役・DE 役）を使い、固定の既知誤差で層分解値・
//! レンダ種別（text/json）・被覆カウント・エラー伝播を独立に固定する。実 DE 経路は
//! **SLOW テスト 1 件のみ**（`JplGoldenComputer` × `EngineGoldenComputer` × 1 ゴールデン）で、
//! de440s.bsp 不在時は `de_diff.rs` と同様 eprintln してスキップする。
//!
//! 合成 `SolarEclipse` / `LocalCircumstances` / `GoldenEclipse` 等の構築ヘルパは
//! `crates/xtask/tests/validate.rs` / `crates/umbra-fixtures/tests/report_golden.rs` のパターンを
//! ミラーする。
//!
//! ## 期待される RED（実装前）
//! `xtask::differential` モジュールも `differential_report` / `JplGoldenComputer` も存在しないため、
//! 本ファイルは **未解決インポート/モジュール（E0432/E0433）でコンパイル不能 = RED** になる。
//! これが想定どおりの赤。

#![allow(clippy::excessive_precision)]

use std::path::Path;

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
    OracleSource,
};

use xtask::differential::{differential_report, JplGoldenComputer};
use xtask::validate::{AccuracyArg, EngineGoldenComputer, ReportFormat};

// ============================================================
// 構築ヘルパ（validate.rs / report_golden.rs のパターンをミラー）
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

/// 合成 `computed: SolarEclipse`（greatest の時刻 2 値・gamma・magnitude を引数、残りフィラー）。
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

/// 計算側の局地接触（最大食 1 点のみ・C1〜C4 は None で済ませる単純構成）。
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

/// 合成 `computed: LocalCircumstances`（最大食時刻＋スカラ metric を引数、C1〜C4 None・metadata フィラー）。
fn computed_local(
    max_utc: UtcInstant,
    max_tt: TtInstant,
    magnitude: f64,
    obscuration: f64,
    max_alt: f64,
    visibility: Visibility,
) -> LocalCircumstances {
    let contacts = LocalContactSet {
        c1: None,
        c2: None,
        maximum: local_contact(max_utc, max_tt),
        c3: None,
        c4: None,
    };
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

/// 合成 `golden: GoldenLocation`（最大食のみ・C1〜C4 None の単純構成。スカラ値を引数化）。
fn golden_location(
    maximum: GoldenContact,
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
        c1: None,
        c2: None,
        maximum,
        c3: None,
        c4: None,
        magnitude,
        obscuration,
        max_altitude_deg,
        max_azimuth_deg: 200.0,
        visibility_expected,
    }
}

/// 合成 `golden: GoldenEclipse`。`event_key`（モックの分岐キー）と全球値・地点列を引数化。TT 付与。
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

/// 1 地点（最大食のみ・全 metric は magnitude=1.0/obsc=1.0/alt=50.0・FullyVisible）。
fn loc_simple() -> GoldenLocation {
    let g_utc = utc(2024, 4, 8, 18, 0, 0.0);
    let base = tt(2_451_545.0, 0.0);
    golden_location(
        golden_contact(g_utc, Some(base), 50.0),
        1.0,
        1.0,
        50.0,
        Visibility::FullyVisible,
    )
}

// ============================================================
// FAST テスト用モック GoldenComputer（固定の既知誤差を注入）
// ============================================================

/// 既知の固定オフセットを golden 値に加えて返すモック computer。
///
/// 役割（analytical / DE）を 1 つの型でパラメタライズ: `greatest_offset_s`（全球最大食時刻に
/// 加える秒・TT 側へ）/ `magnitude_offset`（全球＆地点食分に加える）/ `max_offset_s`（地点最大食
/// 接触時刻に加える秒・TT 側へ）。これにより層分解値（ephemeris=a−d / geometry=d−g / total=a−g）を
/// テスト側で算術的に予言できる。
///
/// `event_key` の接頭辞で分岐:
/// - `"err"`: `eclipse_on` → `Err(NotImplemented)`（エラー伝播テスト用）。
/// - `"missing"`: `eclipse_on` → `Ok(None)`（取りこぼし）。
/// - その他: `eclipse_on` → `Ok(Some(..))`（golden + offset の合成値）。
struct OffsetComputer {
    greatest_offset_s: f64,
    magnitude_offset: f64,
    max_offset_s: f64,
}

/// 1 日の秒数（時刻オフセットを JD 日数へ換算する）。
const SECONDS_PER_DAY: f64 = 86_400.0;

impl OffsetComputer {
    /// 秒オフセットを TT 瞬時へ加える（`add_days` で 2 要素 JD の桁落ちを避ける）。
    fn shift_tt(base: TtInstant, offset_s: f64) -> TtInstant {
        TtInstant::from_jd2(base.jd2().add_days(offset_s / SECONDS_PER_DAY))
    }
}

impl GoldenComputer for OffsetComputer {
    fn eclipse_on(&self, golden: &GoldenEclipse) -> Result<Option<SolarEclipse>, EclipseError> {
        if golden.event_key.starts_with("err") {
            return Err(EclipseError::NotImplemented);
        }
        if golden.event_key.starts_with("missing") {
            return Ok(None);
        }
        let base_tt = golden
            .greatest_time_tt
            .unwrap_or_else(|| tt(2_451_545.0, 0.0));
        Ok(Some(computed_eclipse(
            Self::shift_tt(base_tt, self.greatest_offset_s),
            golden.greatest_time_utc,
            golden.gamma,
            golden.magnitude + self.magnitude_offset,
        )))
    }

    fn local_at(
        &self,
        _eclipse: &SolarEclipse,
        location: &GoldenLocation,
    ) -> Result<LocalCircumstances, EclipseError> {
        let base_tt = location
            .maximum
            .time_tt
            .unwrap_or_else(|| tt(2_451_545.0, 0.0));
        Ok(computed_local(
            location.maximum.time_utc,
            Self::shift_tt(base_tt, self.max_offset_s),
            location.magnitude + self.magnitude_offset,
            location.obscuration,
            location.max_altitude_deg,
            location.visibility_expected,
        ))
    }
}

// ============================================================
// differential_report（FAST・注入モック）
// ============================================================

/// 受け入れ「Text 形式は層分解レポート（render_differential_text 相当）を返す」。
/// `render_differential_text` は 3 層ラベル "ephemeris"/"geometry"/"total" と被覆カウント行
/// （"compared"）を出す。脆い完全一致は避け含有のサニティで縛る。found 2 件 → "2 compared"。
/// 殺す変異: Json 経路に分岐する、render_differential_text を呼ばない、空文字列を返す、
///   層分解せず被覆カウントを通さない。
#[test]
fn differential_report_text_contains_layers_and_counts() {
    let g_utc = utc(2024, 4, 8, 18, 0, 0.0);
    let base = tt(2_451_545.0, 0.0);
    let golden = vec![
        golden_eclipse("a", g_utc, Some(base), 0.4, 1.03, vec![loc_simple()]),
        golden_eclipse("b", g_utc, Some(base), 0.4, 1.03, vec![loc_simple()]),
    ];
    let analytical = OffsetComputer {
        greatest_offset_s: 0.0,
        magnitude_offset: 0.0,
        max_offset_s: 0.0,
    };
    let de = OffsetComputer {
        greatest_offset_s: 0.0,
        magnitude_offset: 0.0,
        max_offset_s: 0.0,
    };
    let s = differential_report(&analytical, &de, &golden, ReportFormat::Text)
        .expect("Text レンダは Ok");
    let lower = s.to_lowercase();
    assert!(
        lower.contains("ephemeris"),
        "Text に ephemeris 層ラベル: {s}"
    );
    assert!(lower.contains("geometry"), "Text に geometry 層ラベル: {s}");
    assert!(lower.contains("total"), "Text に total 層ラベル: {s}");
    assert!(
        s.contains("2 compared"),
        "Text に層分解できた食件数 2 が通っている: {s}"
    );
}

/// 受け入れ「Json 形式は妥当な JSON・末尾改行・主要キーを持ち、被覆カウントが正しく通る」。
/// golden を found 2・missing(de 側 None) 1・found 側合計 3 地点で構築 →
/// eclipses_compared==2, eclipses_missing==1, locations_compared==3。
/// 殺す変異: Text 経路に分岐する、render_differential_json を呼ばない（末尾改行欠落）、
///   被覆カウントを取り違える/入れ替える、report_differential に誤データを渡す。
#[test]
fn differential_report_json_is_object_with_keys_and_counts() {
    let g_utc = utc(2024, 4, 8, 18, 0, 0.0);
    let base = tt(2_451_545.0, 0.0);
    let loc = loc_simple();
    let golden = vec![
        // a: 両エンジン Some・1 地点
        golden_eclipse("a", g_utc, Some(base), 0.4, 1.03, vec![loc.clone()]),
        // missing-b: DE 役が None を返す → eclipses_missing。地点 2 は数えない。
        golden_eclipse(
            "missing-b",
            g_utc,
            Some(base),
            0.4,
            1.03,
            vec![loc.clone(), loc.clone()],
        ),
        // c: 両エンジン Some・2 地点
        golden_eclipse(
            "c",
            g_utc,
            Some(base),
            0.4,
            1.03,
            vec![loc.clone(), loc.clone()],
        ),
    ];
    // analytical は全件 Some。DE は "missing-b" で None を返すよう、同じ接頭辞分岐を使う。
    let analytical = OffsetComputer {
        greatest_offset_s: 0.0,
        magnitude_offset: 0.0,
        max_offset_s: 0.0,
    };
    let de = OffsetComputer {
        greatest_offset_s: 0.0,
        magnitude_offset: 0.0,
        max_offset_s: 0.0,
    };
    let s = differential_report(&analytical, &de, &golden, ReportFormat::Json)
        .expect("Json レンダは Ok");

    assert!(s.ends_with('\n'), "Json は末尾改行で終わる: {s:?}");
    let v: Value = serde_json::from_str(&s).expect("出力は妥当な JSON");
    assert!(v.is_object(), "トップレベルは JSON オブジェクト: {s}");
    let obj = v.as_object().expect("オブジェクト");
    for key in [
        "global_greatest_seconds",
        "global_magnitude",
        "local_maximum_seconds",
        "eclipses_compared",
        "eclipses_missing",
        "locations_compared",
    ] {
        assert!(obj.contains_key(key), "トップレベルキー {key} が存在: {s}");
    }
    assert_eq!(
        v["eclipses_compared"], 2,
        "compared==2（両エンジン Some の golden 数）"
    );
    assert_eq!(
        v["eclipses_missing"], 1,
        "missing==1（DE 役が None を返した golden 数）"
    );
    assert_eq!(
        v["locations_compared"], 3,
        "locations_compared==3（found 側 1+2、missing の 2 地点は数えない）"
    );
}

/// 受け入れ「注入した既知オフセットが 3 層（ephemeris/geometry/total）に算術どおり現れる」。
/// analytical は golden 最大食時刻に +3.0 s、DE は +1.0 s（共に TT 側）を加える（1 地点・接触は最大食のみ）。
/// 期待（符号付き = computed − reference を絶対値統計化, n=1 なので max==mean==p95）:
/// - global_greatest_seconds.ephemeris(a−d) = +3.0−(+1.0) = 2.0
/// - global_greatest_seconds.geometry(d−g)  = +1.0
/// - global_greatest_seconds.total(a−g)     = +3.0
///
/// 殺す変異: ephemeris/geometry/total の取り違え、a/d/g の引数順入れ替え、層分解せず golden 比較に
///   倒す、ErrorStats の max_abs を別フィールドに割り当てる。
#[test]
fn differential_report_json_layer_values_match_injected_offsets() {
    let g_utc = utc(2024, 4, 8, 18, 0, 0.0);
    let base = tt(2_451_545.0, 0.0);
    let golden = vec![golden_eclipse(
        "a",
        g_utc,
        Some(base),
        0.4,
        1.03,
        vec![loc_simple()],
    )];
    let analytical = OffsetComputer {
        greatest_offset_s: 3.0,
        magnitude_offset: 0.0,
        max_offset_s: 0.0,
    };
    let de = OffsetComputer {
        greatest_offset_s: 1.0,
        magnitude_offset: 0.0,
        max_offset_s: 0.0,
    };
    let s = differential_report(&analytical, &de, &golden, ReportFormat::Json)
        .expect("Json レンダは Ok");
    let v: Value = serde_json::from_str(&s).expect("妥当な JSON");

    let layer = &v["global_greatest_seconds"];
    let eps = 1e-6; // shift_tt の JD 換算往復で ~1e-9 級の丸めが入りうるため緩めの EPS。
    let eph = layer["ephemeris"]["max_abs"].as_f64().expect("eph max_abs");
    let geo = layer["geometry"]["max_abs"].as_f64().expect("geo max_abs");
    let tot = layer["total"]["max_abs"].as_f64().expect("tot max_abs");
    assert!(
        (eph - 2.0).abs() < eps,
        "ephemeris(a−d)=3.0−1.0=2.0 だが {eph}（層の取り違え/引数順）: {s}"
    );
    assert!((geo - 1.0).abs() < eps, "geometry(d−g)=1.0 だが {geo}: {s}");
    assert!((tot - 3.0).abs() < eps, "total(a−g)=3.0 だが {tot}: {s}");
}

/// 受け入れ「analytical/DE どちらの `eclipse_on` Err も differential_report 全体へ伝播する」。
/// golden に "err" を 1 件含める → どちらの computer 経由でも Err（EclipseError → XtaskError::Eclipse）。
/// 殺す変異: エラーを握り潰す、unwrap、err を None 扱いにする、片方の computer のエラーだけ見る、Ok を返す。
#[test]
fn differential_report_propagates_computer_error() {
    let g_utc = utc(2024, 4, 8, 18, 0, 0.0);
    let base = tt(2_451_545.0, 0.0);
    let golden = vec![
        golden_eclipse("a", g_utc, Some(base), 0.4, 1.03, vec![]),
        golden_eclipse("err-b", g_utc, Some(base), 0.4, 1.03, vec![]),
    ];
    let analytical = OffsetComputer {
        greatest_offset_s: 0.0,
        magnitude_offset: 0.0,
        max_offset_s: 0.0,
    };
    let de = OffsetComputer {
        greatest_offset_s: 0.0,
        magnitude_offset: 0.0,
        max_offset_s: 0.0,
    };
    let r = differential_report(&analytical, &de, &golden, ReportFormat::Json);
    assert!(
        matches!(r, Err(xtask::error::XtaskError::Eclipse(_))),
        "eclipse_on の EclipseError が XtaskError::Eclipse に写像されて伝播する, got {r:?}"
    );
}

// ============================================================
// JplGoldenComputer::from_spk_path（FAST 異常系）
// ============================================================

/// 受け入れ「存在しない SPK パスは Err」。変種（MalformedSpk/DataUnavailable のどちらにマップされても）
/// は問わず、**Ok にならない**ことだけを縛る（descriptive エラーは実装裁量）。
/// 殺す変異: 不在パスでも Ok を返す、エラーを握り潰してダミーエンジンを返す。
#[test]
fn jpl_golden_computer_from_nonexistent_spk_is_err() {
    let path = Path::new("/no/such/de440s/__definitely_missing__.bsp");
    let r = JplGoldenComputer::from_spk_path(path, AccuracyArg::Standard);
    assert!(
        r.is_err(),
        "存在しない SPK パスは Err であるべき（Ok にしてはならない）"
    );
}

// ============================================================
// SLOW 統合テスト（実 DE440s × 解析エンジン・1 ゴールデンのみ）
// ============================================================

/// 実 DE440s（リポジトリ root の data/spk）。CARGO_MANIFEST_DIR は crates/xtask。
const DE440S_PATH: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../data/spk/de440s.bsp");

/// 【SLOW・1 件】実 `JplGoldenComputer`（DE440s）× `EngineGoldenComputer`（解析暦）で 1 ゴールデンを
/// end-to-end に層分解し、Json が妥当・主要キー・compared==1 を満たすことを縛る。層の数値（実エンジン
/// -vs-DE の本来の差分）は assert しない（精度検証は別テストの領分）。
/// 重い実 `search`（数分）＋少数 `local_circumstances` を 1 ゴールデンに限定して実行時間を抑える。
/// de440s.bsp 不在時は de_diff.rs と同様 eprintln してスキップ（return）。
/// 殺す変異: JplGoldenComputer が DE エンジンを配線していない、differential_report が DE 側を呼ばない、
///   compared の数え違い、Json レンダ未配線。
// SLOW
#[test]
fn differential_report_real_de_one_golden_json() {
    let spk = Path::new(DE440S_PATH);
    let de = match JplGoldenComputer::from_spk_path(spk, AccuracyArg::Standard) {
        Ok(j) => j,
        Err(e) => {
            eprintln!(
                "skip differential_report_real_de_one_golden_json: \
                 {DE440S_PATH} を読めない（{e:?}）。実 DE440s は CI 非同梱（ISSUE-036）。"
            );
            return;
        }
    };
    let analytical = EngineGoldenComputer::new(AccuracyArg::Standard);

    let goldens = golden_eclipses();
    // 2017-08-21-total を優先・無ければ先頭。地点は重いので 1 件に切り詰める。
    let mut target = goldens
        .iter()
        .find(|g| g.event_key == "2017-08-21-total")
        .or_else(|| goldens.first())
        .expect("ゴールデンが 1 件以上ある")
        .clone();
    target.locations.truncate(1);

    let s = differential_report(
        &analytical,
        &de,
        std::slice::from_ref(&target),
        ReportFormat::Json,
    )
    .expect("実 DE × 解析暦の層分解 Json レンダは Ok");

    assert!(!s.is_empty(), "レポート文字列は非空");
    let v: Value = serde_json::from_str(&s).expect("出力は妥当な JSON");
    for key in [
        "global_greatest_seconds",
        "eclipses_compared",
        "locations_compared",
    ] {
        assert!(v.get(key).is_some(), "トップレベルキー {key} が存在: {s}");
    }
    assert_eq!(
        v["eclipses_compared"], 1,
        "両エンジンが食を返せば 1 ゴールデン → compared==1"
    );
}
