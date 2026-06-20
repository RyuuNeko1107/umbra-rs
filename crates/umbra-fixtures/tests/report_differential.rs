//! DE 差分・誤差層分解レポート（accuracy.md §4 — 同一パイプライン 2 エンジン差分法）の受け入れテスト。
//!
//! 本ファイルは `umbra-fixtures` の **公開 API のみ**を対象とした統合テスト（tests/ 配下）。
//! 対象はこれから実装する純オーケストレーション:
//! - `report_differential(analytical, de, golden) -> Result<DifferentialReport, EclipseError>`
//!   （2 つの `GoldenComputer`（analytical 役・DE 役）と golden を取り、各 metric を
//!   ephemeris / geometry / total の 3 層に分解する純オーケストレーション）。
//! - `struct LayeredError { ephemeris, geometry, total }`（1 metric の層分解 = 各層 `ErrorStats`）。
//! - `struct DifferentialReport`（全球＋地点別 metric 毎の `LayeredError` ＋被覆カウント）。
//! - `render_differential_text` / `render_differential_json`（人間可読 / 機械可読出力）。
//!
//! ## 層分解の式（符号付き = computed − reference、集計時に ErrorStats が絶対値化）
//! - ephemeris 層 = analytical − DE（同一パイプライン → 暦差のみ）
//! - geometry 層  = DE − golden（DE 入力 → 幾何/数値＋慣習差）
//! - total 層     = analytical − golden（実測総誤差）
//! - 恒等性: 各サンプルで（符号付き値が）total == ephemeris + geometry が厳密成立
//!   （同じ 3 実数 a, d, g から (a−d)+(d−g)=a−g）。
//!
//! ## テスト戦略（mutation-resistant / FAST）
//! 実エンジンを走らせず、`GoldenComputer` を実装するモックを 2 つ（analytical 役・DE 役）用意し、
//! 返す `SolarEclipse` / `LocalCircumstances` を既知の固定値で構築する。analytical・DE・golden の
//! 3 者の値を完全制御し、各層の max/mean/p95 が手計算で確定するように仕込む。
//!
//! ## 期待される RED（実装前）
//! `LayeredError` / `DifferentialReport` / `report_differential` /
//! `render_differential_text` / `render_differential_json` はまだ存在しないため、本ファイルは
//! **未解決インポート（E0432）でコンパイル不能 = RED** になる。これが想定どおりの赤。

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
    DifferentialReport, GoldenComputer, GoldenContact, GoldenEclipse, GoldenLocation, LayeredError,
    LocationClass, OracleSource,
};

/// 計算統計の f64 一致に用いる厳密許容。
const EPS: f64 = 1e-9;

/// 1 日の秒数（秒差 → JD 差の換算）。
const SECONDS_PER_DAY: f64 = 86_400.0;

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

/// 基準 TT に秒オフセットを足した TT（part2 に足して 2 要素を保つ＝桁落ち回避）。
fn tt_plus_seconds(base: TtInstant, seconds: f64) -> TtInstant {
    let b = base.jd2();
    TtInstant::from_jd2(JulianDate2::new(
        b.part1,
        b.part2 + seconds / SECONDS_PER_DAY,
    ))
}

/// 基準 UTC に秒オフセットを足した UTC（part2 に足して 2 要素を保つ）。
fn utc_plus_seconds(base: UtcInstant, seconds: f64) -> UtcInstant {
    let b = base.jd2();
    UtcInstant::from_jd2(JulianDate2::new(
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

/// 合成 `golden: GoldenEclipse`。`event_key` を引数化、種別・全球値・locations を指定。
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
// FAST テスト用モック GoldenComputer（analytical 役 / DE 役）
// ============================================================

/// 1 接触の存在制御フラグ（analytical / DE で個別に Some/None を切る。golden 側の存在は
/// 渡す `GoldenLocation` 自体の Some/None で決まるためここには持たない）。
#[derive(Clone, Copy)]
struct ContactPresence {
    analytical: bool,
    de: bool,
}

/// 既知の固定誤差を仕込むモック computer。
///
/// `role` が `Analytical` か `De` かで「golden に足すオフセット」を切り替える。これにより
/// 1 つの golden から analytical 値・DE 値を独立に合成し、3 者（a, d, g）を完全制御する。
///
/// - 全球時刻: golden の greatest TT に `greatest_offset_seconds[role]` を足す。
/// - 全球食分: golden の **エンジン規約食分**（中心食なら (1+ratio)/2）に
///   `magnitude_offset[role]` を足す（→ analytical/DE はともにエンジン出力なので換算後を基準にする）。
/// - 地点 maximum 時刻 / magnitude / obscuration / max_altitude: 同様にオフセットを足す。
/// - 接触 c1..c4: `contact_presence` で role 毎に Some/None を切り、Some の接触時刻には
///   maximum と同じ時刻オフセットを足す。
///
/// `eclipse_missing` / `local_err` でオーケストレーションの分岐（取りこぼし・エラー伝播）を制御。
#[derive(Clone, Copy)]
enum Role {
    Analytical,
    De,
}

struct MockComputer {
    role: Role,
    /// greatest TT 時刻オフセット（秒, computed − golden）。
    greatest_offset_seconds: f64,
    /// greatest UTC 時刻オフセット（秒, computed − golden）。golden が TT 無しのとき UTC 経路を縛る。
    greatest_utc_offset_seconds: f64,
    /// 全球食分オフセット（無次元, computed − エンジン規約 golden）。
    global_magnitude_offset: f64,
    /// 地点 maximum 時刻オフセット（秒）。
    local_maximum_offset_seconds: f64,
    /// 地点 magnitude オフセット（無次元）。
    local_magnitude_offset: f64,
    /// 地点 obscuration オフセット（無次元）。
    local_obscuration_offset: f64,
    /// 地点 max_altitude オフセット（度）。
    local_altitude_offset: f64,
    /// 接触時刻オフセット（秒, Some の c1..c4 に足す）。
    contact_offset_seconds: f64,
    /// 各接触の Some/None 制御（c1, c2, c3, c4）。None なら全 role で素直に従う。
    contact_presence: Option<[ContactPresence; 4]>,
    /// eclipse_on が None を返すか。
    eclipse_missing: bool,
    /// eclipse_on が Err を返すか。
    eclipse_err: bool,
    /// local_at が Err を返すか。
    local_err: bool,
}

impl MockComputer {
    /// 既知誤差なし（computed == エンジン規約 golden）の analytical 役 / DE 役を作る。
    fn aligned(role: Role) -> Self {
        MockComputer {
            role,
            greatest_offset_seconds: 0.0,
            greatest_utc_offset_seconds: 0.0,
            global_magnitude_offset: 0.0,
            local_maximum_offset_seconds: 0.0,
            local_magnitude_offset: 0.0,
            local_obscuration_offset: 0.0,
            local_altitude_offset: 0.0,
            contact_offset_seconds: 0.0,
            contact_presence: None,
            eclipse_missing: false,
            eclipse_err: false,
            local_err: false,
        }
    }

    /// golden（NASA）食分をエンジン規約へ換算（report.rs と同一規約）。
    /// 中心食（Total/Annular/Hybrid）は (1+ratio)/2、それ以外は素通し。
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
        let shifted_tt = tt_plus_seconds(base_tt, self.greatest_offset_seconds);
        let shifted_utc =
            utc_plus_seconds(golden.greatest_time_utc, self.greatest_utc_offset_seconds);
        // 全球食分は「エンジン規約 golden」に offset を足す（analytical/DE はともにエンジン出力）。
        let mag = Self::golden_magnitude_engine_convention(golden) + self.global_magnitude_offset;
        computed_eclipse(shifted_tt, shifted_utc, golden.gamma, mag)
    }

    /// role に対応する接触の存在（presence 指定がなければ golden に従う）。
    fn present(&self, idx: usize, golden_has: bool) -> bool {
        match self.contact_presence {
            None => golden_has,
            Some(p) => match self.role {
                Role::Analytical => p[idx].analytical,
                Role::De => p[idx].de,
            },
        }
    }

    fn local_for(&self, loc: &GoldenLocation) -> LocalCircumstances {
        let max_tt = loc.maximum.time_tt.unwrap_or_else(|| tt(2_451_545.0, 0.0));
        let shifted_max = tt_plus_seconds(max_tt, self.local_maximum_offset_seconds);
        let mk = |idx: usize, gc: &Option<GoldenContact>| -> Option<LocalContact> {
            if !self.present(idx, gc.is_some()) {
                return None;
            }
            // Some の接触は golden 時刻（あれば）に contact offset を足す。golden が None でも
            // presence で Some を強制した場合はフィラー基準。
            let (base_utc, base_tt) = match gc.as_ref() {
                Some(g) => (
                    g.time_utc,
                    g.time_tt.unwrap_or_else(|| tt(2_451_545.0, 0.0)),
                ),
                None => (loc.maximum.time_utc, tt(2_451_545.0, 0.0)),
            };
            Some(local_contact(
                base_utc,
                tt_plus_seconds(base_tt, self.contact_offset_seconds),
            ))
        };
        let contacts = LocalContactSet {
            c1: mk(0, &loc.c1),
            c2: mk(1, &loc.c2),
            maximum: local_contact(loc.maximum.time_utc, shifted_max),
            c3: mk(2, &loc.c3),
            c4: mk(3, &loc.c4),
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
        if self.eclipse_err {
            Err(EclipseError::NotImplemented)
        } else if self.eclipse_missing {
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

/// 全 Some 接触の golden 地点（TT 付与）。magnitude/obscuration/altitude は引数。
fn loc_full(magnitude: f64, obscuration: f64, altitude: f64) -> GoldenLocation {
    let g_utc = utc(2024, 4, 8, 18, 0, 0.0);
    let base = tt(2_451_545.0, 0.0);
    let gc = |alt: f64| golden_contact(g_utc, Some(base), alt);
    golden_location(
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

/// 単純な全 Some・全値 golden 一致用の既定地点。
fn loc_default() -> GoldenLocation {
    loc_full(1.0, 1.0, 50.0)
}

/// 1 件の Total golden（greatest TT 付与・指定 locations）。
fn total_golden(key: &str, magnitude: f64, locations: Vec<GoldenLocation>) -> GoldenEclipse {
    let g_utc = utc(2024, 4, 8, 18, 0, 0.0);
    let base = tt(2_451_545.0, 0.0);
    golden_eclipse(
        key,
        SolarEclipseKind::Total,
        g_utc,
        Some(base),
        0.4,
        magnitude,
        locations,
    )
}

// ============================================================
// FAST テスト（mutation-resistant）
// ============================================================

/// 受け入れ「層分解の式（ephemeris=a−d, geometry=d−g, total=a−g）が greatest 時刻で正しい」。
/// 仕込み（秒, computed − golden）: analytical=+5.0, DE=+2.0 → golden 一致。
///   ephemeris = a−d = +3.0 → |3.0|
///   geometry  = d−g = +2.0 → |2.0|
///   total     = a−g = +5.0 → |5.0|
/// 単一サンプルなので max=mean=p95=各層の絶対値。
/// 殺す変異: 層の式の取り違え（a−g を ephemeris に入れる等）、引き算の方向反転。
#[test]
fn layered_split_global_greatest_seconds() {
    let golden = vec![total_golden("e", 1.0, vec![])];
    let analytical = MockComputer {
        greatest_offset_seconds: 5.0,
        ..MockComputer::aligned(Role::Analytical)
    };
    let de = MockComputer {
        greatest_offset_seconds: 2.0,
        ..MockComputer::aligned(Role::De)
    };
    let report: DifferentialReport =
        report_differential_call(&analytical, &de, &golden).expect("モックはエラーなし");

    let l: &LayeredError = &report.global_greatest_seconds;
    assert_eq!(l.ephemeris.units, "s", "greatest は秒");
    assert!(
        (l.ephemeris.max_abs - 3.0).abs() < EPS,
        "ephemeris=a−d=5−2=3, got {}",
        l.ephemeris.max_abs
    );
    assert!(
        (l.geometry.max_abs - 2.0).abs() < EPS,
        "geometry=d−g=2−0=2, got {}",
        l.geometry.max_abs
    );
    assert!(
        (l.total.max_abs - 5.0).abs() < EPS,
        "total=a−g=5−0=5, got {}",
        l.total.max_abs
    );
    // n は found 数（1）。
    assert_eq!(l.ephemeris.n, 1, "ephemeris.n=1");
    assert_eq!(l.geometry.n, 1, "geometry.n=1");
    assert_eq!(l.total.n, 1, "total.n=1");
}

/// 受け入れ「恒等性: 符号付き total == ephemeris + geometry が各サンプルで成立」。
/// 2 サンプルで符号も含めて崩れない設定にする（負の値・絶対値化前の恒等性を統計経由で検証）。
///   サンプル1: a=+10, d=+4 → eph=+6, geo=+4, tot=+10
///   サンプル2: a=−10, d=−4 → eph=−6, geo=−4, tot=−10（同符号で並ぶよう offset を逆向きに）
/// ここでは各 computer のオフセットは固定なので、2 つの golden で同じ符号差が出る。
/// 絶対値化後でも total(=10) == ephemeris(=6)+geometry(=4) が max で成立する。
/// 殺す変異: 層を独立に計算して恒等性を破る式（例: total を別経路で算出して a−g にならない）。
#[test]
fn layered_identity_total_equals_eph_plus_geo() {
    let golden = vec![
        total_golden("e1", 1.0, vec![]),
        total_golden("e2", 1.0, vec![]),
    ];
    let analytical = MockComputer {
        greatest_offset_seconds: 10.0,
        ..MockComputer::aligned(Role::Analytical)
    };
    let de = MockComputer {
        greatest_offset_seconds: 4.0,
        ..MockComputer::aligned(Role::De)
    };
    let report = report_differential_call(&analytical, &de, &golden).expect("モックはエラーなし");
    let l = &report.global_greatest_seconds;
    // 各層の max は全サンプル同一値（6 / 4 / 10）。
    assert!((l.ephemeris.max_abs - 6.0).abs() < EPS, "eph=6");
    assert!((l.geometry.max_abs - 4.0).abs() < EPS, "geo=4");
    assert!((l.total.max_abs - 10.0).abs() < EPS, "tot=10");
    // 恒等性（絶対値・同符号なので max でも成立）。
    assert!(
        (l.total.max_abs - (l.ephemeris.max_abs + l.geometry.max_abs)).abs() < EPS,
        "total == ephemeris + geometry"
    );
    assert_eq!(l.total.n, 2, "2 サンプル集計");
}

/// 受け入れ「全球食分の golden 換算が中心食で効く（geometry/total に出て ephemeris には出ない）」。
/// Total golden の NASA 食分 ratio=0.94 → エンジン規約 = (1+0.94)/2 = 0.97。
/// analytical/DE はともにエンジン規約 0.97 に offset を足す:
///   analytical offset=+0.01 → a=0.98、DE offset=0.00 → d=0.97、エンジン規約 golden g=0.97。
///   ephemeris = a−d = +0.01
///   geometry  = d−g = 0.00  （DE が換算済み golden と一致）
///   total     = a−g = +0.01
/// 換算が無いと g=0.94 になり geometry=0.03・total=0.04 になってしまう（区別が付く）。
/// 殺す変異: 換算の欠落（geometry/total に 0.03 級の差が出る）、ephemeris にも換算を掛ける誤り。
#[test]
fn layered_global_magnitude_golden_conversion_central() {
    let golden = vec![total_golden("e", 0.94, vec![])]; // NASA ratio 0.94 → エンジン規約 0.97
    let analytical = MockComputer {
        global_magnitude_offset: 0.01,
        ..MockComputer::aligned(Role::Analytical)
    };
    let de = MockComputer::aligned(Role::De); // offset 0 → d = エンジン規約 0.97
    let report = report_differential_call(&analytical, &de, &golden).expect("モックはエラーなし");
    let m = &report.global_magnitude;
    assert_eq!(m.ephemeris.units, "", "magnitude は無次元");
    assert!(
        (m.ephemeris.max_abs - 0.01).abs() < EPS,
        "ephemeris=a−d=0.98−0.97=0.01, got {}",
        m.ephemeris.max_abs
    );
    assert!(
        m.geometry.max_abs < EPS,
        "geometry=d−g=0.97−0.97≈0（換算が効いている）, got {}",
        m.geometry.max_abs
    );
    assert!(
        (m.total.max_abs - 0.01).abs() < EPS,
        "total=a−g=0.98−0.97=0.01, got {}",
        m.total.max_abs
    );
}

/// 受け入れ「部分食では golden 換算しない（素通し）— geometry/total に換算が出ない」。
/// Partial golden の食分 0.50（素通し）。analytical offset=+0.02、DE offset=0。
///   geometry = d−g = 0.50−0.50 = 0、total = a−g = 0.52−0.50 = 0.02、ephemeris = a−d = 0.02。
/// 換算を誤って掛けると g=(1+0.5)/2=0.75 になり geometry/total が壊れる（区別が付く）。
/// 殺す変異: 部分食にも (1+ratio)/2 を掛ける誤り。
#[test]
fn layered_global_magnitude_partial_passthrough() {
    let g_utc = utc(2022, 10, 25, 11, 0, 0.0);
    let base = tt(2_451_545.0, 0.0);
    let golden = vec![golden_eclipse(
        "p",
        SolarEclipseKind::Partial,
        g_utc,
        Some(base),
        1.07,
        0.50,
        vec![],
    )];
    let analytical = MockComputer {
        global_magnitude_offset: 0.02,
        ..MockComputer::aligned(Role::Analytical)
    };
    let de = MockComputer::aligned(Role::De);
    let report = report_differential_call(&analytical, &de, &golden).expect("モックはエラーなし");
    let m = &report.global_magnitude;
    assert!(
        m.geometry.max_abs < EPS,
        "部分食は素通し → geometry≈0, got {}",
        m.geometry.max_abs
    );
    assert!(
        (m.total.max_abs - 0.02).abs() < EPS,
        "total=0.52−0.50=0.02, got {}",
        m.total.max_abs
    );
}

/// 受け入れ「地点別 metric（maximum 秒 / magnitude / obscuration / 高度）の層分解」。
/// 1 地点に既知オフセットを仕込み、各 metric で 3 層を縛る。
///   maximum 秒: a=+6, d=+2 → eph=4, geo=2, tot=6
///   magnitude:  a=+0.003, d=+0.001 → eph=0.002, geo=0.001, tot=0.003
///   obscuration:a=+0.030, d=+0.010 → eph=0.020, geo=0.010, tot=0.030
///   altitude:   a=+0.40, d=+0.10 → eph=0.30, geo=0.10, tot=0.40
/// 殺す変異: metric 間の配線取り違え、層の式取り違え、単位ラベル取り違え。
#[test]
fn layered_local_metrics_split() {
    let golden = vec![total_golden("e", 1.0, vec![loc_default()])];
    let analytical = MockComputer {
        local_maximum_offset_seconds: 6.0,
        local_magnitude_offset: 0.003,
        local_obscuration_offset: 0.030,
        local_altitude_offset: 0.40,
        ..MockComputer::aligned(Role::Analytical)
    };
    let de = MockComputer {
        local_maximum_offset_seconds: 2.0,
        local_magnitude_offset: 0.001,
        local_obscuration_offset: 0.010,
        local_altitude_offset: 0.10,
        ..MockComputer::aligned(Role::De)
    };
    let report = report_differential_call(&analytical, &de, &golden).expect("モックはエラーなし");

    let mx = &report.local_maximum_seconds;
    assert_eq!(mx.ephemeris.units, "s");
    assert!((mx.ephemeris.max_abs - 4.0).abs() < EPS, "max秒 eph=4");
    assert!((mx.geometry.max_abs - 2.0).abs() < EPS, "max秒 geo=2");
    assert!((mx.total.max_abs - 6.0).abs() < EPS, "max秒 tot=6");

    let mg = &report.local_magnitude;
    assert_eq!(mg.ephemeris.units, "");
    assert!((mg.ephemeris.max_abs - 0.002).abs() < EPS, "mag eph=0.002");
    assert!((mg.geometry.max_abs - 0.001).abs() < EPS, "mag geo=0.001");
    assert!((mg.total.max_abs - 0.003).abs() < EPS, "mag tot=0.003");

    let ob = &report.local_obscuration;
    assert_eq!(ob.ephemeris.units, "");
    assert!((ob.ephemeris.max_abs - 0.020).abs() < EPS, "obsc eph=0.020");
    assert!((ob.geometry.max_abs - 0.010).abs() < EPS, "obsc geo=0.010");
    assert!((ob.total.max_abs - 0.030).abs() < EPS, "obsc tot=0.030");

    let al = &report.local_max_altitude_deg;
    assert_eq!(al.ephemeris.units, "deg");
    assert!((al.ephemeris.max_abs - 0.30).abs() < EPS, "alt eph=0.30");
    assert!((al.geometry.max_abs - 0.10).abs() < EPS, "alt geo=0.10");
    assert!((al.total.max_abs - 0.40).abs() < EPS, "alt tot=0.40");
}

/// 受け入れ「接触秒は analytical・DE・golden の 3 者すべて Some の接触のみ層分解に寄与」。
/// 1 地点・golden は c1..c4 全 Some。presence を:
///   c1: a=Some, d=Some, g=Some → 寄与する（3 者そろう）
///   c2: a=None, d=Some, g=Some → 全層から落ちる
///   c3: a=Some, d=None, g=Some → 全層から落ちる
///   c4: a=Some, d=Some, g=Some → 寄与する
/// → 寄与は c1, c4 の 2 接触のみ（全層 n=2）。接触 offset: a=+3, d=+1 → eph=2, geo=1, tot=3。
/// 殺す変異: 2 者 Some で寄与させる、None 接触を 0 として混ぜる、片方の層だけ落とす。
#[test]
fn layered_local_contacts_require_three_some() {
    // golden は全 Some。
    let loc = loc_default();
    let golden = vec![total_golden("e", 1.0, vec![loc])];
    let presence = [
        ContactPresence {
            analytical: true,
            de: true,
        }, // c1: a,d Some（golden も Some）→ 3 者
        ContactPresence {
            analytical: false,
            de: true,
        }, // c2: a None → 落ちる
        ContactPresence {
            analytical: true,
            de: false,
        }, // c3: d None → 落ちる
        ContactPresence {
            analytical: true,
            de: true,
        }, // c4: a,d Some（golden も Some）→ 3 者
    ];
    let analytical = MockComputer {
        contact_offset_seconds: 3.0,
        contact_presence: Some(presence),
        ..MockComputer::aligned(Role::Analytical)
    };
    let de = MockComputer {
        contact_offset_seconds: 1.0,
        contact_presence: Some(presence),
        ..MockComputer::aligned(Role::De)
    };
    let report = report_differential_call(&analytical, &de, &golden).expect("モックはエラーなし");
    let c = &report.local_contact_seconds;
    assert_eq!(c.ephemeris.units, "s", "接触は秒");
    // 3 者そろう c1, c4 の 2 接触のみ → 全層 n=2。
    assert_eq!(c.ephemeris.n, 2, "ephemeris.n=2（c1,c4 のみ）");
    assert_eq!(c.geometry.n, 2, "geometry.n=2（c1,c4 のみ）");
    assert_eq!(c.total.n, 2, "total.n=2（c1,c4 のみ）");
    assert!(
        (c.ephemeris.max_abs - 2.0).abs() < EPS,
        "接触 eph=a−d=3−1=2"
    );
    assert!((c.geometry.max_abs - 1.0).abs() < EPS, "接触 geo=d−g=1−0=1");
    assert!((c.total.max_abs - 3.0).abs() < EPS, "接触 tot=a−g=3−0=3");
}

/// 受け入れ「golden が接触を持たない（None）なら 3 者そろわず接触は全層から落ちる」。
/// golden の c2/c3 を None にし、a/d は Some。→ c2/c3 は寄与しない。c1/c4 のみ（n=2）。
/// 殺す変異: golden None を無視して a/d だけで寄与させる。
#[test]
fn layered_local_contacts_golden_none_excluded() {
    let g_utc = utc(2024, 4, 8, 18, 0, 0.0);
    let base = tt(2_451_545.0, 0.0);
    let gc = |alt: f64| golden_contact(g_utc, Some(base), alt);
    // c2, c3 を golden で None に。
    let loc = golden_location(
        Some(gc(10.0)),
        None,
        gc(50.0),
        None,
        Some(gc(40.0)),
        1.0,
        1.0,
        50.0,
        Visibility::FullyVisible,
    );
    let golden = vec![total_golden("e", 1.0, vec![loc])];
    // presence 指定なし → a/d は golden の Some/None に従う（c2/c3 は両者 None）。
    let analytical = MockComputer {
        contact_offset_seconds: 4.0,
        ..MockComputer::aligned(Role::Analytical)
    };
    let de = MockComputer {
        contact_offset_seconds: 1.0,
        ..MockComputer::aligned(Role::De)
    };
    let report = report_differential_call(&analytical, &de, &golden).expect("モックはエラーなし");
    let c = &report.local_contact_seconds;
    assert_eq!(c.total.n, 2, "golden Some の c1,c4 のみ → n=2");
    assert!((c.total.max_abs - 4.0).abs() < EPS, "接触 tot=4−0=4");
}

/// 受け入れ「eclipses_missing: analytical None・DE None・両方 None のいずれでもスキップ＋計上」。
/// 4 件: (a Some, d Some)=比較、(a None, d Some)=missing、(a Some, d None)=missing、
///        (a None, d None)=missing。→ eclipses_compared=1, eclipses_missing=3。
/// missing の食は地点を評価しない。
/// 殺す変異: どれかの None を比較対象にしてしまう、missing を数え落とす。
#[test]
fn missing_when_either_engine_returns_none() {
    // role 毎に eclipse_missing を出し分けるため、event_key で分岐するモックを使う。
    // ここでは presence ベースではなく、各 golden を別 computer で処理できないので、
    // event_key で missing を判定する専用モックを使う。
    struct KeyedMock {
        role: Role,
    }
    impl GoldenComputer for KeyedMock {
        fn eclipse_on(&self, golden: &GoldenEclipse) -> Result<Option<SolarEclipse>, EclipseError> {
            let miss = match self.role {
                Role::Analytical => golden.event_key.contains("aNone"),
                Role::De => golden.event_key.contains("dNone"),
            };
            if miss {
                Ok(None)
            } else {
                let base = tt(2_451_545.0, 0.0);
                Ok(Some(computed_eclipse(
                    base,
                    golden.greatest_time_utc,
                    golden.gamma,
                    (1.0 + golden.magnitude) / 2.0,
                )))
            }
        }
        fn local_at(
            &self,
            _e: &SolarEclipse,
            loc: &GoldenLocation,
        ) -> Result<LocalCircumstances, EclipseError> {
            let max_tt = loc.maximum.time_tt.unwrap_or_else(|| tt(2_451_545.0, 0.0));
            Ok(computed_local(
                LocalContactSet {
                    c1: None,
                    c2: None,
                    maximum: local_contact(loc.maximum.time_utc, max_tt),
                    c3: None,
                    c4: None,
                },
                loc.magnitude,
                loc.obscuration,
                loc.max_altitude_deg,
                loc.visibility_expected,
            ))
        }
    }
    let golden = vec![
        total_golden("both-some", 1.0, vec![loc_default()]),
        total_golden("aNone-1", 1.0, vec![loc_default()]),
        total_golden("dNone-1", 1.0, vec![loc_default()]),
        total_golden("aNone-dNone-1", 1.0, vec![loc_default()]),
    ];
    let report = report_differential_call(
        &KeyedMock {
            role: Role::Analytical,
        },
        &KeyedMock { role: Role::De },
        &golden,
    )
    .expect("モックはエラーなし");
    assert_eq!(report.eclipses_compared, 1, "両 Some は 1 件のみ");
    assert_eq!(report.eclipses_missing, 3, "どちらか None は 3 件");
    assert_eq!(
        report.locations_compared, 1,
        "比較対象食(1件×1地点)のみ → 1（missing 食の地点は数えない）"
    );
}

/// 受け入れ「eclipses_compared / locations_compared の計数（複数地点）」。
/// 両 Some の 2 食（地点 2・地点 3）→ compared=2, locations=5。
/// 殺す変異: compared/locations の数え違い、missing 混入。
#[test]
fn counts_compared_and_locations() {
    let golden = vec![
        total_golden("e1", 1.0, vec![loc_default(), loc_default()]),
        total_golden("e2", 1.0, vec![loc_default(), loc_default(), loc_default()]),
    ];
    let report = report_differential_call(
        &MockComputer::aligned(Role::Analytical),
        &MockComputer::aligned(Role::De),
        &golden,
    )
    .expect("モックはエラーなし");
    assert_eq!(report.eclipses_compared, 2, "両 Some 2 食");
    assert_eq!(report.locations_compared, 5, "2+3=5 地点");
    // 地点 metric の n も地点数に一致（maximum は 1 地点 1 件）。
    assert_eq!(report.local_maximum_seconds.total.n, 5, "max秒 total.n=5");
}

/// 受け入れ「空 golden は全 LayeredError が空統計（n=0・全 0.0・units 保持）＋カウント 0＋Ok」。
/// 殺す変異: 空でパニック、units を落とす、カウントが 0 にならない、Err にする。
#[test]
fn empty_golden_is_vacuous() {
    let golden: Vec<GoldenEclipse> = vec![];
    let report = report_differential_call(
        &MockComputer::aligned(Role::Analytical),
        &MockComputer::aligned(Role::De),
        &golden,
    )
    .expect("空入力は Ok");
    assert_eq!(report.eclipses_compared, 0);
    assert_eq!(report.eclipses_missing, 0);
    assert_eq!(report.locations_compared, 0);
    // 全層が n=0・0.0、units は保持。
    let check = |l: &LayeredError, u: &str| {
        for s in [&l.ephemeris, &l.geometry, &l.total] {
            assert_eq!(s.n, 0, "空 → n=0");
            assert!(s.max_abs == 0.0, "空 → max_abs=0");
            assert!(s.mean_abs == 0.0, "空 → mean_abs=0");
            assert!(s.p95 == 0.0, "空 → p95=0");
            assert_eq!(s.units, u, "空でも units 保持");
        }
    };
    check(&report.global_greatest_seconds, "s");
    check(&report.global_magnitude, "");
    check(&report.local_maximum_seconds, "s");
    check(&report.local_contact_seconds, "s");
    check(&report.local_magnitude, "");
    check(&report.local_obscuration, "");
    check(&report.local_max_altitude_deg, "deg");
}

/// 受け入れ「analytical の eclipse_on の Err は即伝播」。
/// 殺す変異: エラー握り潰し、None 扱い。
#[test]
fn propagates_analytical_eclipse_error() {
    let golden = vec![total_golden("e", 1.0, vec![])];
    let analytical = MockComputer {
        eclipse_err: true,
        ..MockComputer::aligned(Role::Analytical)
    };
    let de = MockComputer::aligned(Role::De);
    let r = report_differential_call(&analytical, &de, &golden);
    assert!(r.is_err(), "analytical eclipse_on Err 伝播, got {r:?}");
}

/// 受け入れ「DE の eclipse_on の Err も即伝播」。
/// 殺す変異: DE 側エラーだけ握り潰す。
#[test]
fn propagates_de_eclipse_error() {
    let golden = vec![total_golden("e", 1.0, vec![])];
    let analytical = MockComputer::aligned(Role::Analytical);
    let de = MockComputer {
        eclipse_err: true,
        ..MockComputer::aligned(Role::De)
    };
    let r = report_differential_call(&analytical, &de, &golden);
    assert!(r.is_err(), "DE eclipse_on Err 伝播, got {r:?}");
}

/// 受け入れ「local_at の Err も即伝播（analytical / DE どちらでも）」。
/// 殺す変異: local_at のエラー握り潰し。
#[test]
fn propagates_local_at_error() {
    let golden = vec![total_golden("e", 1.0, vec![loc_default()])];
    // analytical 側 local_err。
    let a_err = MockComputer {
        local_err: true,
        ..MockComputer::aligned(Role::Analytical)
    };
    let de = MockComputer::aligned(Role::De);
    assert!(
        report_differential_call(&a_err, &de, &golden).is_err(),
        "analytical local_at Err 伝播"
    );
    // DE 側 local_err。
    let analytical = MockComputer::aligned(Role::Analytical);
    let d_err = MockComputer {
        local_err: true,
        ..MockComputer::aligned(Role::De)
    };
    assert!(
        report_differential_call(&analytical, &d_err, &golden).is_err(),
        "DE local_at Err 伝播"
    );
}

/// 受け入れ「render_differential_json は valid JSON（parse 可）＋末尾改行」。
/// 殺す変異: 末尾改行の欠落、非 JSON 文字列。
#[test]
fn render_json_is_valid_with_trailing_newline() {
    let golden = vec![total_golden("e", 1.0, vec![loc_default()])];
    let analytical = MockComputer {
        greatest_offset_seconds: 1.5,
        ..MockComputer::aligned(Role::Analytical)
    };
    let de = MockComputer {
        greatest_offset_seconds: 0.5,
        ..MockComputer::aligned(Role::De)
    };
    let report = report_differential_call(&analytical, &de, &golden).expect("Ok");
    let json = render_differential_json_call(&report).expect("JSON 直列化は成功");
    assert!(json.ends_with('\n'), "末尾改行あり");
    let parsed: serde_json::Value = serde_json::from_str(&json).expect("valid JSON");
    assert!(parsed.is_object(), "トップは object");
    // 主要キーが存在（脆い完全一致は避け、構造の存在のみ）。
    assert!(
        parsed.get("global_greatest_seconds").is_some(),
        "global_greatest_seconds キーあり"
    );
    assert!(
        parsed
            .get("global_greatest_seconds")
            .and_then(|v| v.get("ephemeris"))
            .is_some(),
        "層分解 ephemeris キーあり"
    );
    assert!(
        parsed.get("eclipses_compared").is_some(),
        "compared キーあり"
    );
}

/// 受け入れ「render_differential_text に主要数値・ラベルが含まれる（被覆カウント＋各層 n/max）」。
/// 既知誤差（greatest eph=1.0）を仕込み、ラベルと数値のサニティのみ縛る（完全一致は避ける）。
/// 殺す変異: 層の表示欠落（誤差を隠さない）、被覆カウントの非表示。
#[test]
fn render_text_contains_labels_and_numbers() {
    let golden = vec![total_golden("e", 1.0, vec![loc_default()])];
    let analytical = MockComputer {
        greatest_offset_seconds: 2.0,
        ..MockComputer::aligned(Role::Analytical)
    };
    let de = MockComputer {
        greatest_offset_seconds: 1.0,
        ..MockComputer::aligned(Role::De)
    };
    let report = report_differential_call(&analytical, &de, &golden).expect("Ok");
    let text = render_differential_text_call(&report);
    // 3 層のラベルが出る（誤差を隠さない）。
    assert!(
        text.contains("ephemeris"),
        "ephemeris 層ラベルを表示: {text}"
    );
    assert!(text.contains("geometry"), "geometry 層ラベルを表示: {text}");
    assert!(text.contains("total"), "total 層ラベルを表示: {text}");
    // 被覆カウント（compared=1）が出る。
    assert!(text.contains('1'), "被覆カウント等の数値が出る: {text}");
    // metric ラベルの少なくとも 1 つ。
    assert!(
        text.contains("greatest") || text.contains("maximum"),
        "metric ラベルが出る: {text}"
    );
}

/// 受け入れ「時刻誤差は TT 優先: golden が TT を持てば TT 差、持たなければ UTC 差（geometry 層）」。
/// 同一 UTC・異なる TT の golden 2 種で geometry(=DE−golden) の値が切り替わることを縛る。
/// - golden が TT を持つ場合: DE の greatest TT は golden TT に offset を足したものなので
///   TT 差 = offset（DE offset を 0 にして DE TT == golden TT → geometry≈0）。
/// - golden が TT を持たない場合: UTC 差で評価。DE の greatest UTC は golden UTC 一致なので
///   UTC 差≈0。これも geometry≈0 になるが、TT を持つ場合に DE の TT をずらすと TT 差が出る一方
///   UTC は一致 → TT 経路が使われていることを区別できる。
///
/// 殺す変異: TT/UTC 経路の取り違え（TT があるのに UTC 差を使う）。
#[test]
fn time_error_prefers_tt_when_golden_has_tt() {
    // golden が TT を持つ: DE の TT を golden TT から +3s ずらし、UTC は一致させる。
    // → geometry(d−g) は TT 差 = 3（UTC 差なら 0 になるはずなので区別できる）。
    let g_utc = utc(2024, 4, 8, 18, 0, 0.0);
    let g_tt = tt(2_451_545.0, 0.0);
    let golden_with_tt = vec![golden_eclipse(
        "tt",
        SolarEclipseKind::Total,
        g_utc,
        Some(g_tt),
        0.4,
        1.0,
        vec![],
    )];
    // DE: greatest TT を +3s（UTC は computed_eclipse が golden UTC をそのまま使う＝一致）。
    let de = MockComputer {
        greatest_offset_seconds: 3.0,
        ..MockComputer::aligned(Role::De)
    };
    // analytical: DE と同じ +3s → ephemeris(a−d)=0 に単離（geometry のみ見る）。
    let analytical = MockComputer {
        greatest_offset_seconds: 3.0,
        ..MockComputer::aligned(Role::Analytical)
    };
    let report = report_differential_call(&analytical, &de, &golden_with_tt).expect("Ok");
    let l = &report.global_greatest_seconds;
    assert!(
        (l.geometry.max_abs - 3.0).abs() < EPS,
        "TT を持つ → geometry は TT 差 3（UTC 差なら 0）, got {}",
        l.geometry.max_abs
    );
    assert!(
        l.ephemeris.max_abs < EPS,
        "a,d 同 offset → ephemeris≈0, got {}",
        l.ephemeris.max_abs
    );

    // golden が TT を持たない: UTC 差で評価。DE/analytical の UTC は golden 一致なので geometry≈0。
    let golden_no_tt = vec![golden_eclipse(
        "no-tt",
        SolarEclipseKind::Total,
        g_utc,
        None,
        0.4,
        1.0,
        vec![],
    )];
    let report2 = report_differential_call(&analytical, &de, &golden_no_tt).expect("Ok");
    assert!(
        report2.global_greatest_seconds.geometry.max_abs < EPS,
        "TT 無し → UTC 差で評価、UTC 一致なので geometry≈0, got {}",
        report2.global_greatest_seconds.geometry.max_abs
    );
}

/// 受け入れ「golden が TT を持たないとき、ephemeris 層も UTC 経路で a_utc−d_utc を測る（非零を縛る）」。
/// golden は TT 無し。analytical/DE の greatest **UTC** を独立にずらす（TT はずらさない）:
///   analytical UTC offset=+4s, DE UTC offset=+1s（golden UTC 基準）。
///   ephemeris = a_utc−d_utc = 4−1 = 3、geometry = d_utc−g_utc = 1、total = a_utc−g_utc = 4。
/// 全層が UTC 経路で算出され、恒等性 total(4)=eph(3)+geo(1) も保たれる。
/// 殺す変異: `computed_pair_time_error_seconds` の else（UTC）分岐の削除/0 化/TT 流用、引き算方向反転。
#[test]
fn time_error_uses_utc_when_golden_has_no_tt() {
    let g_utc = utc(2024, 4, 8, 18, 0, 0.0);
    let golden = vec![golden_eclipse(
        "no-tt-utc",
        SolarEclipseKind::Total,
        g_utc,
        None, // TT 無し → UTC 経路
        0.4,
        1.0,
        vec![],
    )];
    let analytical = MockComputer {
        greatest_utc_offset_seconds: 4.0,
        ..MockComputer::aligned(Role::Analytical)
    };
    let de = MockComputer {
        greatest_utc_offset_seconds: 1.0,
        ..MockComputer::aligned(Role::De)
    };
    let report = report_differential_call(&analytical, &de, &golden).expect("Ok");
    let l = &report.global_greatest_seconds;
    assert!(
        (l.ephemeris.max_abs - 3.0).abs() < EPS,
        "ephemeris(UTC)=a_utc−d_utc=4−1=3, got {}",
        l.ephemeris.max_abs
    );
    assert!(
        (l.geometry.max_abs - 1.0).abs() < EPS,
        "geometry(UTC)=d_utc−g_utc=1−0=1, got {}",
        l.geometry.max_abs
    );
    assert!(
        (l.total.max_abs - 4.0).abs() < EPS,
        "total(UTC)=a_utc−g_utc=4−0=4, got {}",
        l.total.max_abs
    );
}

// ============================================================
// 呼び出しラッパ（未実装 API への薄い橋渡し。RED 時は import 不能で失敗）
// ============================================================

/// `report_differential` 呼び出しの薄いラッパ。
fn report_differential_call<A: GoldenComputer, D: GoldenComputer>(
    analytical: &A,
    de: &D,
    golden: &[GoldenEclipse],
) -> Result<DifferentialReport, EclipseError> {
    umbra_fixtures::report_differential(analytical, de, golden)
}

/// `render_differential_json` 呼び出しラッパ。
fn render_differential_json_call(report: &DifferentialReport) -> Result<String, serde_json::Error> {
    umbra_fixtures::render_differential_json(report)
}

/// `render_differential_text` 呼び出しラッパ。
fn render_differential_text_call(report: &DifferentialReport) -> String {
    umbra_fixtures::render_differential_text(report)
}
