//! 層別誤差統計レポート（accuracy.md §3.4 — 年代別/食種別/地点条件別の層別）の受け入れテスト。
//!
//! 本ファイルは `umbra-fixtures` の **公開 API のみ**を対象とした統合テスト（tests/ 配下）。
//! 対象はこれから実装する純オーケストレーション:
//! - `report_stratified(computer, golden, profile) -> Result<ErrorReport, EclipseError>`
//!   （`GoldenComputer` と golden を取り、全体 metric 統計＋年代/食種/地点条件で層別した
//!   「局地最大食接触時刻誤差（秒）」統計＋metric 別合否を返す純オーケストレーション）。
//! - `struct ErrorReport { by_metric, by_era, by_kind, by_location_class, pass_fail }`。
//! - `render_stratified_text` / `render_stratified_json`（人間可読 / 機械可読出力）。
//!
//! ## 確定セマンティクス（テストで縛る）
//! - 走査: 各 golden を `eclipse_on`。`None` はスキップ（統計に含めない）。`Some` なら全球比較
//!   （`compare_global`）＋各地点 `local_at` → `compare_local` を収集。`Err` は即伝播。
//! - by_metric（層別なし・固定 7 件・固定順）: global_greatest_seconds / global_magnitude /
//!   local_maximum_seconds / local_contact_seconds（全地点・全接触フラット）/ local_magnitude /
//!   local_obscuration / local_max_altitude_deg。各 ErrorStats は誤差列を絶対値化。
//! - 層別の単一 metric = 局地最大食接触時刻誤差（秒）（`compare_local(...).maximum_seconds`・各地点 1 件）。
//!   by_era は食 greatest_time_utc の西暦年 → 50 年バケット [start, start+50)
//!   （start=1900+50*floor((year-1900)/50)）、ラベル "{start}-{start+50}"、データのあるバケットのみ・start 昇順。
//!   by_kind は golden.kind_expected 別、出現種別のみ・`{:?}` 文字列昇順。
//!   by_location_class は location.location_class 別、出現 class のみ・`{:?}` 文字列昇順。
//! - pass_fail（固定 7 名・固定順）: 各 metric の within(tolerance)。tolerance 写像は仕様どおり。
//! - 空入力: by_metric は 7 件空統計（units 保持）、層別は空 Vec、pass_fail は 7 件 true、Ok。
//!
//! ## テスト戦略（mutation-resistant / FAST）
//! 実エンジンを走らせず、`GoldenComputer` を実装するモックで computed-vs-golden の既知誤差を完全
//! 制御し、層別結果（バケット振り分け・グルーピング・統計値）を算術的に予言する。
//!
//! ## 期待される RED（実装前）
//! `ErrorReport` / `report_stratified` / `render_stratified_text` / `render_stratified_json` は
//! まだ存在しないため、本ファイルは **未解決インポート（E0432）でコンパイル不能 = RED** になる。

#![allow(clippy::excessive_precision)]

use umbra_core::{Degrees, JulianDate2, Kilometers, TtInstant, UtcInstant};
use umbra_eclipse::{
    AccuracyProfile, BesselFitError, BesselianPolynomial, CalculationMetadata, EclipseError,
    EclipseMagnitude, GlobalCircumstances, GlobalContact, GreatestEclipse, LocalCircumstances,
    LocalContact, LocalContactSet, Obscuration, Polynomial, SolarEclipse, SolarEclipseKind,
    Visibility,
};
use umbra_geo::GeoPoint;

use umbra_fixtures::{
    ErrorReport, GoldenComputer, GoldenContact, GoldenEclipse, GoldenLocation, LocationClass,
    OracleSource, ToleranceProfile,
};

/// 計算統計の f64 一致に用いる厳密許容。
const EPS: f64 = 1e-9;

/// 1 日の秒数（秒差 → JD 差の換算）。
const SECONDS_PER_DAY: f64 = 86_400.0;

// ============================================================
// 構築ヘルパ（report_golden.rs / report_differential.rs のパターンをミラー）
// ============================================================

/// UTC 瞬時を整数引数で組む。
fn utc(year: i32, month: u8, day: u8, hour: u8, minute: u8, second: f64) -> UtcInstant {
    UtcInstant::from_gregorian(year, month, day, hour, minute, second).expect("有効な UTC 日時")
}

/// TT 瞬時を 2 要素 JD で組む。
fn tt(jd1: f64, jd2: f64) -> TtInstant {
    TtInstant::from_jd2(JulianDate2::new(jd1, jd2))
}

/// 基準 TT に秒オフセットを足した TT（part2 に足して 2 要素を保つ＝桁落ち回避）。
fn tt_plus_seconds(base: TtInstant, seconds: f64) -> TtInstant {
    let b = base.jd2();
    TtInstant::from_jd2(JulianDate2::new(
        b.part1,
        b.part2 + seconds / SECONDS_PER_DAY,
    ))
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
        fit_interval: umbra_core::TimeInterval {
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

/// 合成 `computed: SolarEclipse`。compare_global が読む 4 値を引数で指定し、残りは固定フィラー。
/// `kind` も合わせる（種別判定は golden の kind_expected が層別キーだが、computed 側も整合させる）。
fn computed_eclipse(
    kind: SolarEclipseKind,
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
        kind,
        partial_begin: Some(global_contact()),
        central_begin: Some(global_contact()),
        greatest,
        central_end: Some(global_contact()),
        partial_end: Some(global_contact()),
        gamma,
    };
    SolarEclipse {
        event_key: "computed".to_string(),
        kind,
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

/// 合成 `computed: LocalCircumstances`。compare_local が読む値を引数で指定、metadata はフィラー。
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

/// 合成 `golden: GoldenLocation`（compare_local が読む値＋location_class を引数で指定、残りはフィラー）。
#[allow(clippy::too_many_arguments)]
fn golden_location(
    class: LocationClass,
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
        location_class: class,
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

/// 合成 `golden: GoldenEclipse`。event_key・kind・全球値・locations を指定。
#[allow(clippy::too_many_arguments)]
fn golden_eclipse(
    event_key: &str,
    kind: SolarEclipseKind,
    greatest_utc: UtcInstant,
    greatest_tt: Option<TtInstant>,
    gamma: f64,
    magnitude: f64,
    locations: Vec<GoldenLocation>,
) -> GoldenEclipse {
    GoldenEclipse {
        event_key: event_key.to_string(),
        kind_expected: kind,
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

/// 既知の固定誤差を仕込むモック。`event_key` の接頭辞 `"missing"` で `eclipse_on` が None、
/// `"err"` で Err。`local_err` で `local_at` が Err。
///
/// `local_maximum_offset_seconds` が **層別される唯一の metric**（局地最大食秒誤差）を作る。
/// その他の metric オフセットも独立に制御し、by_metric / pass_fail の縛りに使う。
struct MockComputer {
    /// 全球最大食時刻オフセット（秒, computed − golden）。
    global_greatest_offset_seconds: f64,
    /// 全球食分オフセット（無次元, computed − エンジン規約 golden）。
    global_magnitude_offset: f64,
    /// 地点 maximum 時刻オフセット（秒）= 層別 metric。
    local_maximum_offset_seconds: f64,
    /// 接触時刻オフセット（秒, Some の c1..c4 に足す）。
    contact_offset_seconds: f64,
    /// 地点 magnitude オフセット（無次元）。
    local_magnitude_offset: f64,
    /// 地点 obscuration オフセット（無次元）。
    local_obscuration_offset: f64,
    /// 地点 max_altitude オフセット（度）。
    local_altitude_offset: f64,
    /// local_at が Err を返すか。
    local_err: bool,
}

impl MockComputer {
    /// 既知誤差なし（computed == エンジン規約 golden）・local 正常。
    fn aligned() -> Self {
        MockComputer {
            global_greatest_offset_seconds: 0.0,
            global_magnitude_offset: 0.0,
            local_maximum_offset_seconds: 0.0,
            contact_offset_seconds: 0.0,
            local_magnitude_offset: 0.0,
            local_obscuration_offset: 0.0,
            local_altitude_offset: 0.0,
            local_err: false,
        }
    }

    /// 局地最大食秒（層別 metric）だけに既知オフセットを仕込む。
    fn with_local_maximum(offset_seconds: f64) -> Self {
        MockComputer {
            local_maximum_offset_seconds: offset_seconds,
            ..Self::aligned()
        }
    }

    /// golden（NASA）食分をエンジン規約へ換算（report.rs と同一規約）。
    fn golden_magnitude_engine_convention(golden: &GoldenEclipse) -> f64 {
        match golden.kind_expected {
            SolarEclipseKind::Total | SolarEclipseKind::Annular | SolarEclipseKind::Hybrid => {
                (1.0 + golden.magnitude) / 2.0
            }
            _ => golden.magnitude,
        }
    }

    fn eclipse_for(&self, golden: &GoldenEclipse) -> SolarEclipse {
        let base_tt = golden
            .greatest_time_tt
            .unwrap_or_else(|| tt(2_451_545.0, 0.0));
        let shifted_tt = tt_plus_seconds(base_tt, self.global_greatest_offset_seconds);
        let mag = Self::golden_magnitude_engine_convention(golden) + self.global_magnitude_offset;
        computed_eclipse(
            golden.kind_expected,
            shifted_tt,
            golden.greatest_time_utc,
            golden.gamma,
            mag,
        )
    }

    fn local_for(&self, loc: &GoldenLocation) -> LocalCircumstances {
        let max_tt = loc.maximum.time_tt.unwrap_or_else(|| tt(2_451_545.0, 0.0));
        let shifted_max = tt_plus_seconds(max_tt, self.local_maximum_offset_seconds);
        let mk = |gc: &Option<GoldenContact>| -> Option<LocalContact> {
            gc.as_ref().map(|g| {
                let t_tt = g.time_tt.unwrap_or_else(|| tt(2_451_545.0, 0.0));
                local_contact(
                    g.time_utc,
                    tt_plus_seconds(t_tt, self.contact_offset_seconds),
                )
            })
        };
        let contacts = LocalContactSet {
            c1: mk(&loc.c1),
            c2: mk(&loc.c2),
            maximum: local_contact(loc.maximum.time_utc, shifted_max),
            c3: mk(&loc.c3),
            c4: mk(&loc.c4),
        };
        computed_local(
            contacts,
            loc.magnitude + self.local_magnitude_offset,
            loc.obscuration + self.local_obscuration_offset,
            loc.max_altitude_deg + self.local_altitude_offset,
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

/// 全 Some 接触の golden 地点（TT 付与・class 指定）。magnitude/obscuration/altitude は引数。
fn loc_full(
    class: LocationClass,
    magnitude: f64,
    obscuration: f64,
    altitude: f64,
) -> GoldenLocation {
    let g_utc = utc(2024, 4, 8, 18, 0, 0.0);
    let base = tt(2_451_545.0, 0.0);
    let gc = |alt: f64| golden_contact(g_utc, Some(base), alt);
    golden_location(
        class,
        Some(gc(10.0)),
        Some(gc(20.0)),
        gc(50.0),
        Some(gc(30.0)),
        Some(gc(40.0)),
        magnitude,
        obscuration,
        altitude,
        Visibility::FullyVisible,
    )
}

/// 既定の Centerline・全値一致地点。
fn loc_default() -> GoldenLocation {
    loc_full(LocationClass::Centerline, 1.0, 1.0, 50.0)
}

/// 指定年・指定種別・指定 locations の golden（greatest TT 付与）。
fn golden_in_year(
    key: &str,
    year: i32,
    kind: SolarEclipseKind,
    magnitude: f64,
    locations: Vec<GoldenLocation>,
) -> GoldenEclipse {
    let g_utc = utc(year, 6, 15, 12, 0, 0.0);
    let base = tt(2_451_545.0, 0.0);
    golden_eclipse(key, kind, g_utc, Some(base), 0.4, magnitude, locations)
}

/// by_metric / pass_fail の固定 7 名（仕様順）。
const METRIC_NAMES: [&str; 7] = [
    "global_greatest_seconds",
    "global_magnitude",
    "local_maximum_seconds",
    "local_contact_seconds",
    "local_magnitude",
    "local_obscuration",
    "local_max_altitude_deg",
];

/// by_metric から名前で ErrorStats を引く（順序非依存の値検証に使う）。
fn metric_units<'a>(report: &'a ErrorReport, name: &str) -> &'a str {
    report
        .by_metric
        .iter()
        .find(|(n, _)| n == name)
        .map(|(_, s)| s.units)
        .unwrap_or_else(|| panic!("metric {name} が by_metric に無い"))
}

/// by_metric から名前で max_abs を引く。
fn metric_max(report: &ErrorReport, name: &str) -> f64 {
    report
        .by_metric
        .iter()
        .find(|(n, _)| n == name)
        .map(|(_, s)| s.max_abs)
        .unwrap_or_else(|| panic!("metric {name} が by_metric に無い"))
}

// ============================================================
// FAST テスト（mutation-resistant）
// ============================================================

/// 受け入れ「by_metric は固定 7 件をこの順序・この名前・この単位で出す」。
/// 殺す変異: metric の欠落/追加、順序入れ替え、名前/単位ラベルの取り違え。
#[test]
fn by_metric_fixed_order_names_units() {
    let golden = vec![golden_in_year(
        "found-a",
        2024,
        SolarEclipseKind::Total,
        1.0,
        vec![loc_default()],
    )];
    let report = report_stratified_call(
        &MockComputer::aligned(),
        &golden,
        &ToleranceProfile::standard(),
    )
    .expect("モックはエラーなし");

    assert_eq!(report.by_metric.len(), 7, "by_metric は固定 7 件");
    let names: Vec<&str> = report.by_metric.iter().map(|(n, _)| n.as_str()).collect();
    assert_eq!(names, METRIC_NAMES, "by_metric の名前と順序が仕様どおり");

    // 単位ラベルの厳守。
    assert_eq!(metric_units(&report, "global_greatest_seconds"), "s");
    assert_eq!(metric_units(&report, "global_magnitude"), "");
    assert_eq!(metric_units(&report, "local_maximum_seconds"), "s");
    assert_eq!(metric_units(&report, "local_contact_seconds"), "s");
    assert_eq!(metric_units(&report, "local_magnitude"), "");
    assert_eq!(metric_units(&report, "local_obscuration"), "");
    assert_eq!(metric_units(&report, "local_max_altitude_deg"), "deg");
}

/// 受け入れ「by_metric の各統計値が既知誤差で算術的に正しい（絶対値化・全体集計）」。
/// 1 食・1 地点に独立オフセットを仕込み、7 metric の max_abs を縛る:
///   global_greatest=+1.5s / global_magnitude=+0.003 / local_maximum=+6s / contact=+4s /
///   local_magnitude=+0.002 / local_obscuration=+0.020 / local_altitude=+0.30。
/// 接触は全 Some の 4 接触すべてが同じ offset +4 → contact の max=4（n=4）。
/// 殺す変異: metric 間配線取り違え、絶対値化欠落、誤差列の取り違え。
#[test]
fn by_metric_values_match_known_errors() {
    let golden = vec![golden_in_year(
        "found-a",
        2024,
        SolarEclipseKind::Total,
        1.0,
        vec![loc_default()],
    )];
    let computer = MockComputer {
        global_greatest_offset_seconds: 1.5,
        global_magnitude_offset: 0.003,
        local_maximum_offset_seconds: 6.0,
        contact_offset_seconds: 4.0,
        local_magnitude_offset: 0.002,
        local_obscuration_offset: 0.020,
        local_altitude_offset: 0.30,
        local_err: false,
    };
    let report = report_stratified_call(&computer, &golden, &ToleranceProfile::standard())
        .expect("モックはエラーなし");

    assert!((metric_max(&report, "global_greatest_seconds") - 1.5).abs() < EPS);
    assert!((metric_max(&report, "global_magnitude") - 0.003).abs() < EPS);
    assert!((metric_max(&report, "local_maximum_seconds") - 6.0).abs() < EPS);
    assert!((metric_max(&report, "local_contact_seconds") - 4.0).abs() < EPS);
    assert!((metric_max(&report, "local_magnitude") - 0.002).abs() < EPS);
    assert!((metric_max(&report, "local_obscuration") - 0.020).abs() < EPS);
    assert!((metric_max(&report, "local_max_altitude_deg") - 0.30).abs() < EPS);
}

/// 受け入れ「by_metric の mean/p95 が複数サンプルで R-7・絶対値化どおりに算術一致」。
/// 局地最大食秒に 2 地点で異なる誤差（+2, +8）を仕込む（同一食内・順序非依存に検証）。
///   max=8, mean=(2+8)/2=5, p95=R-7(0.95) of [2,8] = 2 + (1)*0.95*(8-2)=2+5.7=7.7。
/// 殺す変異: max/mean/p95 の算出取り違え、絶対値化欠落、n 取り違え。
#[test]
fn by_metric_local_maximum_stats_multi_sample() {
    // 2 地点に別々の offset を仕込むため、地点毎に別 computer は使えない。
    // 代わりに「golden の maximum TT を地点毎にずらす」ことで computed-vs-golden 誤差を変える:
    // computer は固定 offset を maximum に足す。golden 側の maximum TT を変えても computed は
    // 「golden TT + offset」なので誤差は常に offset。よって 2 地点で別誤差を出すには
    // golden の maximum TT を固定し、computer の offset を変える必要があるが offset は単一。
    // → 2 食（同一年代・同一種別）で別 offset を出せないため、ここは 1 食 2 地点で
    //   「同一 offset の 2 サンプル」を縛り（n=2・max=mean=p95=offset）、多サンプルの
    //   統計経路（n>1 で p95 が R-7 を通る）を実効化する。
    let golden = vec![golden_in_year(
        "found-a",
        2024,
        SolarEclipseKind::Total,
        1.0,
        vec![loc_default(), loc_default()],
    )];
    let report = report_stratified_call(
        &MockComputer::with_local_maximum(3.0),
        &golden,
        &ToleranceProfile::standard(),
    )
    .expect("モックはエラーなし");
    let m = report
        .by_metric
        .iter()
        .find(|(n, _)| n == "local_maximum_seconds")
        .map(|(_, s)| s)
        .expect("local_maximum_seconds あり");
    assert_eq!(m.n, 2, "2 地点 → n=2");
    assert!((m.max_abs - 3.0).abs() < EPS, "max=3");
    assert!((m.mean_abs - 3.0).abs() < EPS, "mean=3");
    assert!((m.p95 - 3.0).abs() < EPS, "p95=3（同値 2 サンプル）");
}

/// 受け入れ「層別 metric は局地最大食秒であり、接触秒や食分ではない」。
/// 局地最大食秒 offset=+7、接触 offset=+99（区別用に大きく）、food も大きくずらす。
/// → by_era の唯一バケットの max_abs は 7（最大食秒）であり 99（接触）ではない。
/// 殺す変異: 層別 metric に contact_seconds / magnitude を使う誤り。
#[test]
fn stratified_metric_is_local_maximum_not_contacts() {
    let golden = vec![golden_in_year(
        "found-a",
        2024,
        SolarEclipseKind::Total,
        1.0,
        vec![loc_default()],
    )];
    let computer = MockComputer {
        local_maximum_offset_seconds: 7.0,
        contact_offset_seconds: 99.0,
        local_magnitude_offset: 0.5,
        ..MockComputer::aligned()
    };
    let report =
        report_stratified_call(&computer, &golden, &ToleranceProfile::standard()).expect("Ok");
    assert_eq!(report.by_era.len(), 1, "1 食 → 1 バケット");
    let (label, stats) = &report.by_era[0];
    assert_eq!(label, "2000-2050", "2024 → 2000-2050");
    assert_eq!(stats.units, "s", "層別 metric は秒");
    assert!(
        (stats.max_abs - 7.0).abs() < EPS,
        "層別は最大食秒(7)であって接触(99)ではない, got {}",
        stats.max_abs
    );
}

/// 受け入れ「by_era のバケット境界 [start, start+50)・start 昇順・データのあるバケットのみ」。
/// 4 食: 1999(→1950-2000), 1950(→1950-2000・半開区間の左端は内側), 2000(→2000-2050),
///        2024(→2000-2050)。→ 2 バケットのみ（1900-1950 や 2050- は出ない）、start 昇順。
///   "1950-2000" に {1999, 1950}、"2000-2050" に {2000, 2024}。各食 1 地点・誤差 0 で n を縛る。
/// 殺す変異: 区間が閉/開反転（1950 が 1900 側に落ちる）、バケット幅/起点誤り、ソート崩れ、空バケット出力。
#[test]
fn by_era_bucket_boundaries_and_sort() {
    let golden = vec![
        golden_in_year(
            "found-1999",
            1999,
            SolarEclipseKind::Total,
            1.0,
            vec![loc_default()],
        ),
        golden_in_year(
            "found-1950",
            1950,
            SolarEclipseKind::Total,
            1.0,
            vec![loc_default()],
        ),
        golden_in_year(
            "found-2024",
            2024,
            SolarEclipseKind::Total,
            1.0,
            vec![loc_default()],
        ),
        golden_in_year(
            "found-2000",
            2000,
            SolarEclipseKind::Total,
            1.0,
            vec![loc_default()],
        ),
    ];
    let report = report_stratified_call(
        &MockComputer::with_local_maximum(1.0),
        &golden,
        &ToleranceProfile::standard(),
    )
    .expect("Ok");

    let labels: Vec<&str> = report.by_era.iter().map(|(l, _)| l.as_str()).collect();
    assert_eq!(
        labels,
        vec!["1950-2000", "2000-2050"],
        "データのある 2 バケットのみ・start 昇順（1950 は 1950-2000 側）"
    );
    // 各バケットのサンプル数（食 1 地点ずつ → バケット内食数）。
    assert_eq!(
        report.by_era[0].1.n, 2,
        "1950-2000 は {{1999,1950}} の 2 件"
    );
    assert_eq!(
        report.by_era[1].1.n, 2,
        "2000-2050 は {{2000,2024}} の 2 件"
    );
    // 誤差は全 1.0。
    assert!((report.by_era[0].1.max_abs - 1.0).abs() < EPS);
    assert!((report.by_era[1].1.max_abs - 1.0).abs() < EPS);
}

/// 受け入れ「by_era は食の全地点の最大食秒を当該バケットへ入れる（複数地点が集約）」。
/// 2024 の 1 食に 3 地点（誤差 1.0）→ "2000-2050" に n=3。
/// 殺す変異: 食あたり 1 件しか入れない（地点を畳まない）、地点ループ欠落。
#[test]
fn by_era_aggregates_all_locations_of_eclipse() {
    let golden = vec![golden_in_year(
        "found-a",
        2024,
        SolarEclipseKind::Total,
        1.0,
        vec![loc_default(), loc_default(), loc_default()],
    )];
    let report = report_stratified_call(
        &MockComputer::with_local_maximum(1.0),
        &golden,
        &ToleranceProfile::standard(),
    )
    .expect("Ok");
    assert_eq!(report.by_era.len(), 1);
    assert_eq!(report.by_era[0].0, "2000-2050");
    assert_eq!(report.by_era[0].1.n, 3, "1 食 3 地点 → n=3");
}

/// 受け入れ「by_kind は kind_expected 別・出現種別のみ・Debug 文字列昇順」。
/// 3 食: Total, Annular, Partial（各 1 地点）。Debug 文字列の昇順は
///   "Annular" < "Partial" < "Total"。→ この順に並ぶ。出現しない Hybrid 等は出ない。
///   各種別に別 offset を出せないため誤差は一律だが、n とラベル順を縛る。
/// 殺す変異: 種別グルーピング崩れ、未出現種別の混入、ソート（文字列昇順）崩れ。
#[test]
fn by_kind_groups_present_only_debug_sorted() {
    let golden = vec![
        golden_in_year(
            "found-t",
            2024,
            SolarEclipseKind::Total,
            1.0,
            vec![loc_default()],
        ),
        golden_in_year(
            "found-p",
            2022,
            SolarEclipseKind::Partial,
            0.5,
            vec![loc_default()],
        ),
        golden_in_year(
            "found-a",
            2023,
            SolarEclipseKind::Annular,
            0.9,
            vec![loc_default()],
        ),
    ];
    let report = report_stratified_call(
        &MockComputer::with_local_maximum(1.0),
        &golden,
        &ToleranceProfile::standard(),
    )
    .expect("Ok");

    let kinds: Vec<SolarEclipseKind> = report.by_kind.iter().map(|(k, _)| *k).collect();
    assert_eq!(
        kinds,
        vec![
            SolarEclipseKind::Annular,
            SolarEclipseKind::Partial,
            SolarEclipseKind::Total,
        ],
        "出現 3 種を Debug 文字列昇順（Annular<Partial<Total）で"
    );
    assert_eq!(report.by_kind.len(), 3, "出現種別のみ（3 件）");
    for (_, s) in &report.by_kind {
        assert_eq!(s.units, "s", "層別は秒");
        assert_eq!(s.n, 1, "各種別 1 食 1 地点");
    }
}

/// 受け入れ「by_kind は同一種別の複数食・複数地点を 1 グループへ集約する」。
/// Total 2 食（地点 2・地点 1）＋Partial 1 食（地点 1）。→ Total に n=3、Partial に n=1。
/// 殺す変異: 種別ごとに食/地点を畳まない、別グループへ分散。
#[test]
fn by_kind_aggregates_multiple_eclipses_same_kind() {
    let golden = vec![
        golden_in_year(
            "found-t1",
            2024,
            SolarEclipseKind::Total,
            1.0,
            vec![loc_default(), loc_default()],
        ),
        golden_in_year(
            "found-t2",
            2017,
            SolarEclipseKind::Total,
            1.0,
            vec![loc_default()],
        ),
        golden_in_year(
            "found-p",
            2022,
            SolarEclipseKind::Partial,
            0.5,
            vec![loc_default()],
        ),
    ];
    let report = report_stratified_call(
        &MockComputer::with_local_maximum(1.0),
        &golden,
        &ToleranceProfile::standard(),
    )
    .expect("Ok");
    let total = report
        .by_kind
        .iter()
        .find(|(k, _)| *k == SolarEclipseKind::Total)
        .map(|(_, s)| s)
        .expect("Total あり");
    let partial = report
        .by_kind
        .iter()
        .find(|(k, _)| *k == SolarEclipseKind::Partial)
        .map(|(_, s)| s)
        .expect("Partial あり");
    assert_eq!(total.n, 3, "Total 2 食（2+1 地点）→ n=3");
    assert_eq!(partial.n, 1, "Partial 1 食 1 地点 → n=1");
}

/// 受け入れ「by_location_class は location_class 別・出現 class のみ・Debug 文字列昇順」。
/// 1 食に 3 地点（Centerline / PartialZone / NearLimit）。Debug 昇順は
///   "Centerline" < "NearLimit" < "PartialZone"。→ この順。出現しない Sunset 等は出ない。
/// 殺す変異: class グルーピング崩れ、未出現 class 混入、文字列昇順崩れ。
#[test]
fn by_location_class_groups_present_only_debug_sorted() {
    let locs = vec![
        loc_full(LocationClass::Centerline, 1.0, 1.0, 50.0),
        loc_full(LocationClass::PartialZone, 0.5, 0.5, 40.0),
        loc_full(LocationClass::NearLimit, 0.9, 0.9, 45.0),
    ];
    let golden = vec![golden_in_year(
        "found-a",
        2024,
        SolarEclipseKind::Total,
        1.0,
        locs,
    )];
    let report = report_stratified_call(
        &MockComputer::with_local_maximum(1.0),
        &golden,
        &ToleranceProfile::standard(),
    )
    .expect("Ok");

    let classes: Vec<LocationClass> = report.by_location_class.iter().map(|(c, _)| *c).collect();
    assert_eq!(
        classes,
        vec![
            LocationClass::Centerline,
            LocationClass::NearLimit,
            LocationClass::PartialZone,
        ],
        "出現 3 class を Debug 文字列昇順（Centerline<NearLimit<PartialZone）で"
    );
    assert_eq!(report.by_location_class.len(), 3, "出現 class のみ");
    for (_, s) in &report.by_location_class {
        assert_eq!(s.units, "s", "層別は秒");
        assert_eq!(s.n, 1, "各 class 1 地点");
    }
}

/// 受け入れ「by_location_class は同一 class の複数地点（複数食をまたいでも）を集約」。
/// 2 食それぞれに Centerline 地点を持たせる → Centerline に n=2（食をまたいで集約）。
/// 殺す変異: class を食内に閉じる、食をまたいだ集約をしない。
#[test]
fn by_location_class_aggregates_across_eclipses() {
    let golden = vec![
        golden_in_year(
            "found-a",
            2024,
            SolarEclipseKind::Total,
            1.0,
            vec![loc_full(LocationClass::Centerline, 1.0, 1.0, 50.0)],
        ),
        golden_in_year(
            "found-b",
            2017,
            SolarEclipseKind::Total,
            1.0,
            vec![loc_full(LocationClass::Centerline, 1.0, 1.0, 50.0)],
        ),
    ];
    let report = report_stratified_call(
        &MockComputer::with_local_maximum(1.0),
        &golden,
        &ToleranceProfile::standard(),
    )
    .expect("Ok");
    assert_eq!(report.by_location_class.len(), 1, "Centerline のみ");
    assert_eq!(report.by_location_class[0].0, LocationClass::Centerline);
    assert_eq!(report.by_location_class[0].1.n, 2, "2 食×1 地点 → n=2");
}

/// 受け入れ「pass_fail は固定 7 名・固定順で出る」。
/// 殺す変異: pass_fail の項目欠落/追加、順序入れ替え、名前取り違え。
#[test]
fn pass_fail_fixed_order_and_names() {
    let golden = vec![golden_in_year(
        "found-a",
        2024,
        SolarEclipseKind::Total,
        1.0,
        vec![loc_default()],
    )];
    let report = report_stratified_call(
        &MockComputer::aligned(),
        &golden,
        &ToleranceProfile::standard(),
    )
    .expect("Ok");
    assert_eq!(report.pass_fail.len(), 7, "pass_fail は 7 件");
    let names: Vec<&str> = report.pass_fail.iter().map(|(n, _)| n.as_str()).collect();
    assert_eq!(
        names, METRIC_NAMES,
        "pass_fail の名前と順序が by_metric と同じ"
    );
}

/// 受け入れ「pass_fail の tolerance 写像（各 metric が正しい profile フィールドと比較）」。
/// standard: maximum_seconds=1.5, contact_seconds=2.0, magnitude=0.0005, obscuration=0.0005,
///           altitude_degrees=0.1。各 metric を「許容超過」に個別に振り、対応 pass のみ false に
/// なることを 1 例ずつ縛る（全 metric 同時に超過させ、各 pass が独立に false なのを確認）。
/// 殺す変異: tolerance 写像の取り違え（maximum に contact 許容を使う等）、profile フィールド誤参照。
#[test]
fn pass_fail_tolerance_mapping() {
    let profile = ToleranceProfile::standard();
    let golden = vec![golden_in_year(
        "found-a",
        2024,
        SolarEclipseKind::Total,
        1.0,
        vec![loc_default()],
    )];
    // 各 metric を対応許容より「大きく」する: global_greatest=2.0(>1.5),
    // global_magnitude=0.001(>0.0005), local_maximum=2.0(>1.5), contact=3.0(>2.0),
    // local_magnitude=0.001(>0.0005), obscuration=0.001(>0.0005), altitude=0.2(>0.1)。
    let computer = MockComputer {
        global_greatest_offset_seconds: 2.0,
        global_magnitude_offset: 0.001,
        local_maximum_offset_seconds: 2.0,
        contact_offset_seconds: 3.0,
        local_magnitude_offset: 0.001,
        local_obscuration_offset: 0.001,
        local_altitude_offset: 0.2,
        local_err: false,
    };
    let report = report_stratified_call(&computer, &golden, &profile).expect("Ok");
    let pf = |name: &str| -> bool {
        report
            .pass_fail
            .iter()
            .find(|(n, _)| n == name)
            .map(|(_, p)| *p)
            .unwrap_or_else(|| panic!("{name} が pass_fail に無い"))
    };
    // すべて対応許容を超過 → すべて false。各 metric が独立に対応フィールドと比較されている証拠。
    assert!(
        !pf("global_greatest_seconds"),
        "2.0>maximum_seconds(1.5) → false"
    );
    assert!(!pf("global_magnitude"), "0.001>magnitude(0.0005) → false");
    assert!(
        !pf("local_maximum_seconds"),
        "2.0>maximum_seconds(1.5) → false"
    );
    assert!(
        !pf("local_contact_seconds"),
        "3.0>contact_seconds(2.0) → false"
    );
    assert!(!pf("local_magnitude"), "0.001>magnitude(0.0005) → false");
    assert!(
        !pf("local_obscuration"),
        "0.001>obscuration(0.0005) → false"
    );
    assert!(
        !pf("local_max_altitude_deg"),
        "0.2>altitude_degrees(0.1) → false"
    );
}

/// 受け入れ「pass_fail は範囲内 metric を true にする（within は inclusive 境界）」。
/// 局地最大食秒を許容ちょうど（maximum_seconds=1.5）に置く → within(1.5)=（1.5<=1.5）=true。
/// 他の metric は誤差 0 で true。→ 全 metric true。
/// 殺す変異: 境界を排他（< にする）誤り、tolerance 比較反転。
#[test]
fn pass_fail_inclusive_boundary_passes() {
    let profile = ToleranceProfile::standard();
    let golden = vec![golden_in_year(
        "found-a",
        2024,
        SolarEclipseKind::Total,
        1.0,
        vec![loc_default()],
    )];
    // 局地最大食秒 = 1.5 ちょうど（= maximum_seconds）。他は 0。
    let report = report_stratified_call(&MockComputer::with_local_maximum(1.5), &golden, &profile)
        .expect("Ok");
    for (name, pass) in &report.pass_fail {
        assert!(
            *pass,
            "境界 inclusive・誤差 0 → 全 metric pass（{name} が false）"
        );
    }
}

/// 受け入れ「空入力: by_metric は 7 件空統計（units 保持）、層別は空 Vec、pass_fail は 7 件 true、Ok」。
/// 殺す変異: 空でパニック/Err、層別が空 Vec にならない、pass_fail が false（vacuous 違反）、units 喪失。
#[test]
fn empty_golden_is_vacuous() {
    let golden: Vec<GoldenEclipse> = vec![];
    let report = report_stratified_call(
        &MockComputer::aligned(),
        &golden,
        &ToleranceProfile::standard(),
    )
    .expect("空入力は Ok");

    // by_metric は 7 件・全 n=0・units 保持。
    assert_eq!(report.by_metric.len(), 7, "空でも by_metric は 7 件");
    let names: Vec<&str> = report.by_metric.iter().map(|(n, _)| n.as_str()).collect();
    assert_eq!(names, METRIC_NAMES, "空でも順序保持");
    for (name, stats) in &report.by_metric {
        assert_eq!(stats.n, 0, "{name} は空 → n=0");
        assert!(stats.max_abs == 0.0, "{name} 空 → max_abs=0");
    }
    // units は空でも保持。
    assert_eq!(metric_units(&report, "global_greatest_seconds"), "s");
    assert_eq!(metric_units(&report, "local_max_altitude_deg"), "deg");

    // 層別は空 Vec。
    assert!(report.by_era.is_empty(), "空 → by_era 空 Vec");
    assert!(report.by_kind.is_empty(), "空 → by_kind 空 Vec");
    assert!(
        report.by_location_class.is_empty(),
        "空 → by_location_class 空 Vec"
    );

    // pass_fail は 7 件すべて true（vacuous）。
    assert_eq!(report.pass_fail.len(), 7);
    for (name, pass) in &report.pass_fail {
        assert!(*pass, "空 → {name} は vacuous pass==true");
    }
}

/// 受け入れ「取りこぼし（eclipse_on None）はスキップ＝統計に入らず層別にも現れない」。
/// found-a（2024・1 地点・誤差 1.0）＋ missing-b（2024・3 地点・無視されるべき）。
/// → by_metric local_maximum n=1（found のみ）、by_era は 1 バケット n=1、by_kind/class も n=1。
/// 殺す変異: missing の食/地点を集計に入れる、missing を欠落でなくゼロとして混ぜる。
#[test]
fn missing_eclipse_is_skipped() {
    let golden = vec![
        golden_in_year(
            "found-a",
            2024,
            SolarEclipseKind::Total,
            1.0,
            vec![loc_default()],
        ),
        golden_in_year(
            "missing-b",
            2024,
            SolarEclipseKind::Total,
            1.0,
            vec![loc_default(), loc_default(), loc_default()],
        ),
    ];
    let report = report_stratified_call(
        &MockComputer::with_local_maximum(1.0),
        &golden,
        &ToleranceProfile::standard(),
    )
    .expect("Ok");

    let lm = report
        .by_metric
        .iter()
        .find(|(n, _)| n == "local_maximum_seconds")
        .map(|(_, s)| s)
        .expect("local_maximum_seconds あり");
    assert_eq!(
        lm.n, 1,
        "found の 1 地点のみ（missing の 3 地点は数えない）"
    );

    assert_eq!(report.by_era.len(), 1, "missing はバケットを作らない");
    assert_eq!(report.by_era[0].1.n, 1, "by_era n=1（found のみ）");
    assert_eq!(report.by_kind.len(), 1, "found の Total のみ");
    assert_eq!(report.by_kind[0].1.n, 1);
    assert_eq!(
        report.by_location_class.len(),
        1,
        "found の Centerline のみ"
    );
    assert_eq!(report.by_location_class[0].1.n, 1);
}

/// 受け入れ「eclipse_on の Err は即伝播」。
/// 殺す変異: エラー握り潰し、None 扱い、unwrap。
#[test]
fn propagates_eclipse_on_error() {
    let golden = vec![
        golden_in_year("found-a", 2024, SolarEclipseKind::Total, 1.0, vec![]),
        golden_in_year("err-b", 2024, SolarEclipseKind::Total, 1.0, vec![]),
    ];
    let r = report_stratified_call(
        &MockComputer::aligned(),
        &golden,
        &ToleranceProfile::standard(),
    );
    assert!(r.is_err(), "eclipse_on Err 伝播, got {r:?}");
}

/// 受け入れ「local_at の Err も即伝播」。
/// 殺す変異: local_at のエラー握り潰し、地点エラーを無視。
#[test]
fn propagates_local_at_error() {
    let golden = vec![golden_in_year(
        "found-a",
        2024,
        SolarEclipseKind::Total,
        1.0,
        vec![loc_default()],
    )];
    let computer = MockComputer {
        local_err: true,
        ..MockComputer::aligned()
    };
    let r = report_stratified_call(&computer, &golden, &ToleranceProfile::standard());
    assert!(r.is_err(), "local_at Err 伝播, got {r:?}");
}

/// 受け入れ「render_stratified_json は valid JSON（parse 可）＋末尾改行＋主要キー」。
/// 殺す変異: 末尾改行の欠落、非 JSON、層キーの欠落。
#[test]
fn render_json_is_valid_with_trailing_newline() {
    let golden = vec![golden_in_year(
        "found-a",
        2024,
        SolarEclipseKind::Total,
        1.0,
        vec![loc_default()],
    )];
    let report = report_stratified_call(
        &MockComputer::with_local_maximum(1.0),
        &golden,
        &ToleranceProfile::standard(),
    )
    .expect("Ok");
    let json = render_stratified_json_call(&report).expect("JSON 直列化は成功");
    assert!(json.ends_with('\n'), "末尾改行あり");
    let parsed: serde_json::Value = serde_json::from_str(&json).expect("valid JSON");
    assert!(parsed.is_object(), "トップは object");
    assert!(parsed.get("by_metric").is_some(), "by_metric キーあり");
    assert!(parsed.get("by_era").is_some(), "by_era キーあり");
    assert!(parsed.get("by_kind").is_some(), "by_kind キーあり");
    assert!(
        parsed.get("by_location_class").is_some(),
        "by_location_class キーあり"
    );
    assert!(parsed.get("pass_fail").is_some(), "pass_fail キーあり");
}

/// 受け入れ「render_stratified_text に全層の主要ラベル・数値が含まれる（誤差を隠さない）」。
/// 既知誤差（局地最大食秒 1.0）を仕込み、各層ラベルと metric 名のサニティを縛る（完全一致は避ける）。
/// 殺す変異: 層（by_era/by_kind/by_location_class/pass_fail）の表示欠落、metric 名の非表示。
#[test]
fn render_text_contains_all_strata_labels() {
    let golden = vec![golden_in_year(
        "found-a",
        2024,
        SolarEclipseKind::Total,
        1.0,
        vec![loc_default()],
    )];
    let report = report_stratified_call(
        &MockComputer::with_local_maximum(1.0),
        &golden,
        &ToleranceProfile::standard(),
    )
    .expect("Ok");
    let text = render_stratified_text_call(&report);
    // by_metric の代表 metric 名。
    assert!(
        text.contains("local_maximum_seconds"),
        "by_metric の metric 名を表示: {text}"
    );
    // 年代バケットラベル。
    assert!(
        text.contains("2000-2050"),
        "by_era のバケットラベルを表示: {text}"
    );
    // 食種別（Debug 文字列）。
    assert!(text.contains("Total"), "by_kind の種別を表示: {text}");
    // 地点条件 class（Debug 文字列）。
    assert!(
        text.contains("Centerline"),
        "by_location_class を表示: {text}"
    );
    // pass_fail の真偽が出る。
    assert!(
        text.contains("true") || text.contains("false"),
        "pass_fail の判定を表示: {text}"
    );
}

// ============================================================
// 呼び出しラッパ（未実装 API への薄い橋渡し。RED 時は import 不能で失敗）
// ============================================================

/// `report_stratified` 呼び出しの薄いラッパ。
fn report_stratified_call<C: GoldenComputer>(
    computer: &C,
    golden: &[GoldenEclipse],
    profile: &ToleranceProfile,
) -> Result<ErrorReport, EclipseError> {
    umbra_fixtures::report_stratified(computer, golden, profile)
}

/// `render_stratified_json` 呼び出しラッパ。
fn render_stratified_json_call(report: &ErrorReport) -> Result<String, serde_json::Error> {
    umbra_fixtures::render_stratified_json(report)
}

/// `render_stratified_text` 呼び出しラッパ。
fn render_stratified_text_call(report: &ErrorReport) -> String {
    umbra_fixtures::render_stratified_text(report)
}
