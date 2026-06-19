//! ISSUE-030 S30b 受け入れテスト（strict / 純全球比較: compare_global + aggregate_global）。
//!
//! 本ファイルは `umbra-fixtures` の **公開 API のみ**を対象とした統合テスト（tests/ 配下）。
//! 対象は本スライスの **純全球比較**のみ:
//! - `compare_global(computed, golden) -> GlobalErrors`（符号付き = computed − golden）。
//! - `aggregate_global(errors, profile) -> GlobalReport`（metric 別 `ErrorStats` ＋ 合否）。
//!
//! ## スコープ外（後続スライス・本ファイルでは検証しない）
//! エンジン駆動の `report_against_golden`・地点別比較・JSON/CLI 出力は後続スライス。
//! 本スライスは **エンジンを一切走らせず**、合成（手構築）値のみで純比較を縛る。
//!
//! ## オラクル戦略（mutation-resistant）
//! 期待値はすべてテスト側で手計算した literal。各フィールド・各分岐・符号・境界を独立に縛る。
//! `compare_global` の符号は **computed − golden** に固定（被減数/減数の入れ替えを撃破）。
//! `ErrorStats::from_errors` は内部で abs を取るため、集計統計は |誤差| 上の値（手計算 R-7 p95）。
//!
//! ## 期待される RED（実装前）
//! `GlobalErrors` / `compare_global` / `GlobalReport` / `aggregate_global` はまだ存在しないため、
//! 本ファイルは **未解決インポート（E0432/E0425）でコンパイル不能 = RED** になる。これが想定どおりの赤。

#![allow(clippy::excessive_precision)]

use umbra_core::{Degrees, JulianDate2, Kilometers, TimeInterval, TtInstant, UtcInstant};
use umbra_eclipse::{
    AccuracyProfile, BesselFitError, BesselianPolynomial, CalculationMetadata, EclipseMagnitude,
    GlobalCircumstances, GlobalContact, GreatestEclipse, Obscuration, Polynomial, SolarEclipse,
    SolarEclipseKind,
};
use umbra_geo::GeoPoint;

use umbra_fixtures::{
    aggregate_global, compare_global, ErrorStats, GlobalErrors, GlobalReport, GoldenEclipse,
    OracleSource, ToleranceProfile,
};

/// 計算統計の f64 一致に用いる厳密許容。
const EPS: f64 = 1e-12;
/// TT 秒差（×86400 換算）の許容。`Δs/86400` の往復で浮動小数の丸めが出るため緩め。
const SEC_EPS: f64 = 1e-6;

// ============================================================
// 構築ヘルパ（results.rs の test ヘルパをミラー）
// ============================================================

/// UTC 瞬時を整数引数で組む。
fn utc(year: i32, month: u8, day: u8, hour: u8, minute: u8, second: f64) -> UtcInstant {
    UtcInstant::from_gregorian(year, month, day, hour, minute, second).expect("有効な UTC 日時")
}

/// TT 瞬時を 2 要素 JD で組む。`jd2().jd()` は part1+part2。
fn tt(jd1: f64, jd2: f64) -> TtInstant {
    TtInstant::from_jd2(JulianDate2::new(jd1, jd2))
}

/// 地表点（lat, lon）を度から組む。
fn geo(lat: f64, lon: f64) -> GeoPoint {
    GeoPoint::from_degrees(lat, lon).expect("有効な地表点")
}

/// 固定の最小 BesselianPolynomial（フィラー。compare_global は読まない）。
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

/// 固定の全球接触点（フィラー。compare_global は読まない）。
fn contact() -> GlobalContact {
    GlobalContact {
        time_utc: utc(2024, 4, 8, 16, 0, 0.0),
        time_tt: tt(2_460_409.0, 0.01),
        position: geo(30.0, -100.0),
    }
}

/// 合成 `computed: SolarEclipse`。compare_global が読む 4 値（greatest.time_tt /
/// greatest.time_utc / greatest.magnitude / global.gamma）だけを引数で指定し、残りは固定フィラー。
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
        partial_begin: Some(contact()),
        central_begin: Some(contact()),
        greatest,
        central_end: Some(contact()),
        partial_end: Some(contact()),
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

/// 合成 `golden: GoldenEclipse`。比較に効く 4 値（greatest_time_utc / greatest_time_tt /
/// gamma / magnitude）だけを指定し、locations は空・最小 OracleSource。
fn golden_eclipse(
    greatest_utc: UtcInstant,
    greatest_tt: Option<TtInstant>,
    gamma: f64,
    magnitude: f64,
) -> GoldenEclipse {
    GoldenEclipse {
        event_key: "x".into(),
        kind_expected: SolarEclipseKind::Total,
        greatest_time_utc: greatest_utc,
        greatest_time_tt: greatest_tt,
        gamma,
        magnitude,
        delta_t_seconds: None,
        locations: vec![],
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
// compare_global — greatest_seconds（TT 優先・符号付き ×86400）
// ============================================================

/// 受け入れ「greatest_seconds: golden に TT があれば TT を優先・符号は computed − golden」。
/// golden TT = base、computed TT = base + 1.5 s → greatest_seconds ≈ +1.5（computed が後）。
/// 第二ケース: computed TT = base − 2.0 s → ≈ −2.0。computed の UTC は **わざと遠い別値**にして、
/// TT 分岐が UTC を読んでいないことを証明する。
/// 殺す変異: 符号反転（golden − computed）, TT があるのに UTC を使う, ×86400 欠落, 被減数/減数入替。
#[test]
fn compare_global_greatest_seconds_tt_preferred_signed() {
    let base = tt(2_451_545.0, 0.0);
    let g_utc = utc(2024, 4, 8, 18, 0, 0.0);
    let golden = golden_eclipse(g_utc, Some(base), 0.40, 1.0300);

    // computed TT は golden TT より +1.5 s 後。UTC は遠い別日（TT 優先の証明）。
    let computed_after = computed_eclipse(
        tt(2_451_545.0, 1.5 / 86400.0),
        utc(2000, 1, 1, 0, 0, 0.0),
        0.40,
        1.0300,
    );
    let e_after = compare_global(&computed_after, &golden);
    assert!(
        (e_after.greatest_seconds - 1.5).abs() < SEC_EPS,
        "computed が +1.5 s 後 → greatest_seconds ≈ +1.5（computed − golden, TT 優先）, got {}",
        e_after.greatest_seconds
    );

    // computed TT は golden TT より −2.0 s 前 → 符号は負。
    let computed_before = computed_eclipse(
        tt(2_451_545.0, -2.0 / 86400.0),
        utc(2000, 1, 1, 0, 0, 0.0),
        0.40,
        1.0300,
    );
    let e_before = compare_global(&computed_before, &golden);
    assert!(
        (e_before.greatest_seconds + 2.0).abs() < SEC_EPS,
        "computed が −2.0 s 前 → greatest_seconds ≈ −2.0（符号: computed − golden）, got {}",
        e_before.greatest_seconds
    );
}

/// 受け入れ「greatest_seconds: golden の TT が None なら UTC へフォールバック」。
/// golden TT=None、computed UTC は golden UTC より +3.0 s 後 → greatest_seconds ≈ +3.0。
/// computed TT は **無関係な別値**にして、この分岐で TT が読まれていないことを証明する。
/// 殺す変異: 常に TT を読む（computed TT を拾ってしまう）, 誤フォールバック, 符号反転。
#[test]
fn compare_global_greatest_seconds_falls_back_to_utc_when_tt_none() {
    let g_utc = utc(2024, 4, 8, 18, 0, 0.0);
    let golden = golden_eclipse(g_utc, None, 0.40, 1.0300);

    // computed UTC は golden UTC より +3.0 s 後。computed TT は無関係（読まれないはず）。
    let computed = computed_eclipse(
        tt(2_400_000.0, 0.0),
        utc(2024, 4, 8, 18, 0, 3.0),
        0.40,
        1.0300,
    );
    let e = compare_global(&computed, &golden);
    assert!(
        (e.greatest_seconds - 3.0).abs() < SEC_EPS,
        "TT=None → UTC フォールバックで +3.0 s（computed − golden）, got {}",
        e.greatest_seconds
    );
}

// ============================================================
// compare_global — gamma / magnitude（符号付き）
// ============================================================

/// 受け入れ「gamma: computed − golden の符号付き差」。
/// computed 0.43, golden 0.40 → +0.03。負ケース computed 0.10, golden 0.50 → −0.40。
/// TT は両ケースで一致させ gamma だけを動かす（gamma フィールド読み取りを単離）。
/// 殺す変異: 符号反転, 別フィールド（magnitude/time）読み取り。
#[test]
fn compare_global_gamma_signed() {
    let base = tt(2_451_545.0, 0.0);
    let g_utc = utc(2024, 4, 8, 18, 0, 0.0);

    let golden_pos = golden_eclipse(g_utc, Some(base), 0.40, 1.0300);
    let computed_pos = computed_eclipse(base, g_utc, 0.43, 1.0300);
    let e_pos = compare_global(&computed_pos, &golden_pos);
    assert!(
        (e_pos.gamma - 0.03).abs() < EPS,
        "gamma = 0.43 − 0.40 = +0.03, got {}",
        e_pos.gamma
    );

    let golden_neg = golden_eclipse(g_utc, Some(base), 0.50, 1.0300);
    let computed_neg = computed_eclipse(base, g_utc, 0.10, 1.0300);
    let e_neg = compare_global(&computed_neg, &golden_neg);
    assert!(
        (e_neg.gamma + 0.40).abs() < EPS,
        "gamma = 0.10 − 0.50 = −0.40, got {}",
        e_neg.gamma
    );
}

/// 受け入れ「magnitude: computed.magnitude.0 − golden.magnitude の符号付き差」。
/// computed 1.0306, golden 1.0300 → +0.0006。負ケース computed 1.0290, golden 1.0300 → −0.0010。
/// 殺す変異: 符号反転, gamma など別フィールド読み取り, EclipseMagnitude(.0) の取り違え。
#[test]
fn compare_global_magnitude_signed() {
    let base = tt(2_451_545.0, 0.0);
    let g_utc = utc(2024, 4, 8, 18, 0, 0.0);

    let golden = golden_eclipse(g_utc, Some(base), 0.40, 1.0300);
    let computed_pos = computed_eclipse(base, g_utc, 0.40, 1.0306);
    let e_pos = compare_global(&computed_pos, &golden);
    assert!(
        (e_pos.magnitude - 0.0006).abs() < EPS,
        "magnitude = 1.0306 − 1.0300 = +0.0006, got {}",
        e_pos.magnitude
    );

    let computed_neg = computed_eclipse(base, g_utc, 0.40, 1.0290);
    let e_neg = compare_global(&computed_neg, &golden);
    assert!(
        (e_neg.magnitude + 0.0010).abs() < EPS,
        "magnitude = 1.0290 − 1.0300 = −0.0010, got {}",
        e_neg.magnitude
    );
}

/// 受け入れ「フィールド独立: 3 誤差が互いに異なる非ゼロ値で各フィールドに正しく配線される」。
/// time(+1.5 s) / gamma(+0.03) / magnitude(+0.0006) を **同時に異値**で与え、3 フィールドを独立検証。
/// 殺す変異: フィールド間の取り違え（time→gamma, gamma→magnitude 等のミスルーティング）。
#[test]
fn compare_global_fields_are_independent() {
    let base = tt(2_451_545.0, 0.0);
    let g_utc = utc(2024, 4, 8, 18, 0, 0.0);
    let golden = golden_eclipse(g_utc, Some(base), 0.40, 1.0300);

    let computed = computed_eclipse(
        tt(2_451_545.0, 1.5 / 86400.0),
        utc(2000, 1, 1, 0, 0, 0.0),
        0.43,
        1.0306,
    );
    let e = compare_global(&computed, &golden);
    assert!(
        (e.greatest_seconds - 1.5).abs() < SEC_EPS,
        "greatest_seconds は +1.5（gamma/magnitude と混線しない）, got {}",
        e.greatest_seconds
    );
    assert!(
        (e.gamma - 0.03).abs() < EPS,
        "gamma は +0.03（time/magnitude と混線しない）, got {}",
        e.gamma
    );
    assert!(
        (e.magnitude - 0.0006).abs() < EPS,
        "magnitude は +0.0006（time/gamma と混線しない）, got {}",
        e.magnitude
    );
}

// ============================================================
// aggregate_global — metric 別 ErrorStats ＋ 合否
// ============================================================

/// 受け入れ「metric 別統計の構築（abs 集計・単位・R-7 p95）」。
/// greatest_seconds=[+1.0,−2.0,+0.5], gamma=[+0.4,−0.1,+0.2], magnitude=[+0.0003,−0.0002,+0.0001]。
/// greatest: |e|=[1,2,0.5] → n=3, max_abs=2.0, mean_abs=(1+2+0.5)/3, 昇順[0.5,1,2] の R-7 p95:
///   h=(3-1)*0.95=1.9, lo=1 → 1.0 + 0.9*(2.0-1.0)=1.9。units="s"。
/// gamma: |e|=[0.4,0.1,0.2] → max_abs=0.4, units="Re"。
/// magnitude: |e|=[0.0003,0.0002,0.0001] → max_abs=0.0003, units=""。
/// 殺す変異: metric→stats の配線違い, 単位取り違え, 符号付き（abs 忘れ）投入で max が変わる。
#[test]
fn aggregate_global_builds_per_metric_stats() {
    let errors = [
        GlobalErrors {
            greatest_seconds: 1.0,
            gamma: 0.4,
            magnitude: 0.0003,
        },
        GlobalErrors {
            greatest_seconds: -2.0,
            gamma: -0.1,
            magnitude: -0.0002,
        },
        GlobalErrors {
            greatest_seconds: 0.5,
            gamma: 0.2,
            magnitude: 0.0001,
        },
    ];
    let report: GlobalReport = aggregate_global(&errors, &ToleranceProfile::standard());

    // greatest（units "s"）。
    assert_eq!(report.greatest.n, 3, "greatest.n must be 3");
    assert!(
        (report.greatest.max_abs - 2.0).abs() < EPS,
        "greatest.max_abs must be 2.0 (abs of -2.0), got {}",
        report.greatest.max_abs
    );
    let greatest_mean = (1.0 + 2.0 + 0.5) / 3.0;
    assert!(
        (report.greatest.mean_abs - greatest_mean).abs() < 1e-9,
        "greatest.mean_abs must be (1+2+0.5)/3, got {}",
        report.greatest.mean_abs
    );
    assert!(
        (report.greatest.p95 - 1.9).abs() < EPS,
        "greatest.p95 (R-7 over [0.5,1,2]) must be 1.9, got {}",
        report.greatest.p95
    );
    assert_eq!(report.greatest.units, "s", "greatest.units must be 's'");

    // gamma（units "Re", 非 gated だが統計は出る）。
    assert_eq!(report.gamma.n, 3, "gamma.n must be 3");
    assert!(
        (report.gamma.max_abs - 0.4).abs() < EPS,
        "gamma.max_abs must be 0.4, got {}",
        report.gamma.max_abs
    );
    assert_eq!(report.gamma.units, "Re", "gamma.units must be 'Re'");

    // magnitude（units ""）。
    assert_eq!(report.magnitude.n, 3, "magnitude.n must be 3");
    assert!(
        (report.magnitude.max_abs - 0.0003).abs() < EPS,
        "magnitude.max_abs must be 0.0003, got {}",
        report.magnitude.max_abs
    );
    assert_eq!(report.magnitude.units, "", "magnitude.units must be ''");
}

/// 受け入れ「pass=true でも統計を必ず出す（accuracy.md §3.4）」。
/// 全誤差が standard（greatest ±1.5 s, magnitude ±0.0005）内 → pass==true。
/// それでも各 ErrorStats は n>0・非ゼロ（誤差を隠さない）。
/// 殺す変異: pass が統計をゼロ化する / 空にする。
#[test]
fn aggregate_global_pass_true_within_tolerance() {
    let errors = [
        GlobalErrors {
            greatest_seconds: 1.0,
            gamma: 0.5,
            magnitude: 0.0003,
        },
        GlobalErrors {
            greatest_seconds: -0.8,
            gamma: -0.2,
            magnitude: -0.0002,
        },
    ];
    let report = aggregate_global(&errors, &ToleranceProfile::standard());
    assert!(report.pass, "全誤差が standard 内 → pass==true");
    // pass でも統計は populated。
    assert_eq!(report.greatest.n, 2, "greatest.n は 2（統計は出す）");
    assert!(
        report.greatest.max_abs > 0.0,
        "greatest.max_abs は非ゼロ（誤差を隠さない）, got {}",
        report.greatest.max_abs
    );
    assert_eq!(report.magnitude.n, 2, "magnitude.n は 2（統計は出す）");
    assert!(
        report.magnitude.max_abs > 0.0,
        "magnitude.max_abs は非ゼロ, got {}",
        report.magnitude.max_abs
    );
}

/// 受け入れ「greatest が許容超過なら pass==false」。
/// greatest_seconds=2.0（> standard.maximum_seconds=1.5）, magnitude は許容内。
/// 殺す変異: greatest を gating しない, maximum_seconds 以外の許容フィールドを誤参照。
#[test]
fn aggregate_global_pass_false_when_greatest_exceeds() {
    let errors = [GlobalErrors {
        greatest_seconds: 2.0, // > 1.5
        gamma: 0.0,
        magnitude: 0.0001, // < 0.0005（OK）
    }];
    let report = aggregate_global(&errors, &ToleranceProfile::standard());
    assert!(
        !report.pass,
        "greatest 2.0 s > 1.5 s → pass==false（greatest を gating する）"
    );
}

/// 受け入れ「magnitude が許容超過なら pass==false」。
/// magnitude=0.001（> standard.magnitude=0.0005）, greatest は許容内。
/// 殺す変異: magnitude を gating しない, magnitude 許容フィールドの誤参照。
#[test]
fn aggregate_global_pass_false_when_magnitude_exceeds() {
    let errors = [GlobalErrors {
        greatest_seconds: 0.5, // < 1.5（OK）
        gamma: 0.0,
        magnitude: 0.001, // > 0.0005
    }];
    let report = aggregate_global(&errors, &ToleranceProfile::standard());
    assert!(
        !report.pass,
        "magnitude 0.001 > 0.0005 → pass==false（magnitude を gating する）"
    );
}

/// 受け入れ「gamma は gating しない（統計のみ）」。
/// gamma=10.0（巨大）だが greatest/magnitude は許容内 → pass==true。
/// それでも gamma.max_abs は 10.0 を反映（報告はする）。
/// 殺す変異: pass 判定に gamma を AND する変異。
#[test]
fn aggregate_global_gamma_not_gated() {
    let errors = [GlobalErrors {
        greatest_seconds: 0.5, // OK
        gamma: 10.0,           // 巨大だが非 gated
        magnitude: 0.0001,     // OK
    }];
    let report = aggregate_global(&errors, &ToleranceProfile::standard());
    assert!(
        report.pass,
        "gamma が巨大でも greatest/magnitude が OK なら pass==true（gamma は非 gated）"
    );
    assert!(
        (report.gamma.max_abs - 10.0).abs() < EPS,
        "gamma.max_abs は 10.0 を反映（統計は報告する）, got {}",
        report.gamma.max_abs
    );
}

/// 受け入れ「空列は vacuous pass・ゼロ統計」。
/// aggregate_global(&[], standard) → 3 ErrorStats すべて n=0・全 0.0、pass==true（0.0 <= 正の許容）。
/// 殺す変異: 空でパニック, 空で pass==false。
#[test]
fn aggregate_global_empty_is_vacuous_pass_zero_stats() {
    let report = aggregate_global(&[], &ToleranceProfile::standard());
    assert!(report.pass, "空列は vacuous pass==true");

    let zero = |s: &ErrorStats, label: &str| {
        assert_eq!(s.n, 0, "{label}.n must be 0 for empty");
        assert!(
            s.max_abs == 0.0,
            "{label}.max_abs must be 0.0, got {}",
            s.max_abs
        );
        assert!(
            s.mean_abs == 0.0 && !s.mean_abs.is_nan(),
            "{label}.mean_abs must be 0.0 (not NaN), got {}",
            s.mean_abs
        );
        assert!(s.p95 == 0.0, "{label}.p95 must be 0.0, got {}", s.p95);
    };
    zero(&report.greatest, "greatest");
    zero(&report.gamma, "gamma");
    zero(&report.magnitude, "magnitude");
}
