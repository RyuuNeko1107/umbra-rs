//! 1900-2100 全食スイープの**自己カタログ集計＋完備性突合**（純関数）の受け入れテスト。
//!
//! 本ファイルは `umbra-fixtures` の **公開 API のみ**を対象とした統合テスト（tests/ 配下）。
//! 対象はこれから実装する純集計（実装は別担当・本ファイルはテストのみ）:
//! - `summarize_sweep(eclipses: &[SolarEclipse], expected: CatalogCounts) -> SweepSummary`
//!   （検出食群の自己集計＝種別件数・|gamma|/magnitude の範囲統計・NASA 4 区分への完備性突合）。
//! - `struct SweepSummary { total, by_kind, gamma_abs, magnitude, completeness }`。
//! - `struct RangeStats { n, min, max, mean }`（素の値の min/max/mean・絶対値化しない）。
//! - `struct CatalogCounts { total, annular, hybrid, partial }`（NASA 4 区分）。
//! - `struct CompletenessReport { detected, expected, all_match }`。
//! - `render_sweep_text` / `render_sweep_json`（人間可読 / 機械可読出力）。
//!
//! ## 確定セマンティクス（テストで縛る）
//! 1. total = `eclipses.len()`。
//! 2. by_kind = 各 `eclipse.kind`（raw エンジン種別）の件数・**出現種別のみ**・`{:?}` 文字列昇順。
//!    NonCentral も raw のまま別カウント（畳まない）。
//! 3. gamma_abs = 各 `eclipse.global.gamma` の**絶対値** |gamma| の min/max/mean。空は n=0・全 0.0。
//! 4. magnitude = 各 `eclipse.global.greatest.magnitude.0` の min/max/mean（**素の値**）。空は 0。
//! 5. completeness.detected: NASA 4 区分へ畳む（total←Total+NonCentralTotal,
//!    annular←Annular+NonCentralAnnular, hybrid←Hybrid, partial←Partial）。
//! 6. completeness.expected = 引数 `expected`（そのまま）。
//! 7. completeness.all_match = `detected == expected`（4 区分すべて一致で true）。
//! 8. 空 eclipses: total=0・by_kind 空 Vec・range n=0/0.0・detected 全 0・all_match は expected の全 0 性。
//! 9. render_sweep_json: serde pretty＋末尾改行。render_sweep_text: 人間可読・全項目を漏れなく。
//!
//! ## テスト戦略（mutation-resistant / FAST）
//! 実エンジンを走らせず、種別・gamma・magnitude を既知値で制御した合成 `SolarEclipse` 群を作り
//! （report_stratified.rs のヘルパをミラー）、集計を算術的に予言する。
//!
//! ## 期待される RED（実装前）
//! `SweepSummary` / `summarize_sweep` / `render_sweep_*` / `RangeStats` / `CatalogCounts` /
//! `CompletenessReport` はまだ存在しないため、本ファイルは **未解決インポート（E0432）で
//! コンパイル不能 = RED** になる。これが想定どおりの赤。

#![allow(clippy::excessive_precision)]

use umbra_core::{Degrees, JulianDate2, Kilometers, TtInstant, UtcInstant};
use umbra_eclipse::{
    AccuracyProfile, BesselFitError, BesselianPolynomial, CalculationMetadata, EclipseMagnitude,
    GlobalCircumstances, GlobalContact, GreatestEclipse, Obscuration, Polynomial, SolarEclipse,
    SolarEclipseKind,
};
use umbra_geo::GeoPoint;

use umbra_fixtures::{
    summarize_sweep, CatalogCounts, CompletenessReport, RangeStats, SweepSummary,
};

/// 統計の f64 一致に用いる厳密許容。
const EPS: f64 = 1e-9;

// ============================================================
// 構築ヘルパ（report_stratified.rs のパターンをミラー）
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

/// 合成 `computed: SolarEclipse`。集計が読む 3 値（kind / global.gamma /
/// greatest.magnitude）を引数で指定し、残りは固定フィラー。時刻はスイープ集計に無関係なので固定。
fn eclipse_with(kind: SolarEclipseKind, gamma: f64, magnitude: f64) -> SolarEclipse {
    let greatest = GreatestEclipse {
        time_utc: utc(2024, 4, 8, 18, 0, 0.0),
        time_tt: tt(2_460_409.0, 0.25),
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

/// `CatalogCounts` を組む短縮。
fn counts(total: usize, annular: usize, hybrid: usize, partial: usize) -> CatalogCounts {
    CatalogCounts {
        total,
        annular,
        hybrid,
        partial,
    }
}

/// by_kind から指定種別の件数を引く（出現しなければ None）。
fn kind_count(summary: &SweepSummary, kind: SolarEclipseKind) -> Option<usize> {
    summary
        .by_kind
        .iter()
        .find(|(k, _)| *k == kind)
        .map(|(_, c)| *c)
}

// ============================================================
// total
// ============================================================

/// 受け入れ「total = eclipses.len()」（複数件で len と一致）。
/// 殺す変異: total を len 以外（定数・by_kind 合計の誤算）にする、+1/-1 のオフバイワン。
#[test]
fn total_equals_len() {
    let eclipses = vec![
        eclipse_with(SolarEclipseKind::Total, 0.1, 1.05),
        eclipse_with(SolarEclipseKind::Annular, 0.2, 0.95),
        eclipse_with(SolarEclipseKind::Partial, 0.3, 0.50),
    ];
    let summary = summarize_sweep(&eclipses, counts(0, 0, 0, 0));
    assert_eq!(summary.total, 3, "total は eclipses.len()");
}

// ============================================================
// by_kind（raw 種別・出現のみ・Debug 文字列昇順・NonCentral 別カウント）
// ============================================================

/// 受け入れ「by_kind は各 raw 種別の件数・出現種別のみ・Debug 文字列昇順」。
/// Total×2, Annular×1, Partial×3 を投入 → 出現 3 種が "Annular"<"Partial"<"Total" の順で
/// 件数 (1, 3, 2)。出現しない Hybrid 等は現れない。
/// 殺す変異: 件数の取り違え、出現しない種別の混入、ソート（文字列昇順）崩れ、種別グルーピング崩れ。
#[test]
fn by_kind_counts_present_only_debug_sorted() {
    let eclipses = vec![
        eclipse_with(SolarEclipseKind::Total, 0.1, 1.05),
        eclipse_with(SolarEclipseKind::Partial, 0.2, 0.40),
        eclipse_with(SolarEclipseKind::Annular, 0.3, 0.95),
        eclipse_with(SolarEclipseKind::Partial, 0.4, 0.42),
        eclipse_with(SolarEclipseKind::Total, 0.5, 1.02),
        eclipse_with(SolarEclipseKind::Partial, 0.6, 0.38),
    ];
    let summary = summarize_sweep(&eclipses, counts(0, 0, 0, 0));

    let pairs: Vec<(SolarEclipseKind, usize)> = summary.by_kind.clone();
    assert_eq!(
        pairs,
        vec![
            (SolarEclipseKind::Annular, 1),
            (SolarEclipseKind::Partial, 3),
            (SolarEclipseKind::Total, 2),
        ],
        "出現 3 種を Debug 文字列昇順（Annular<Partial<Total）で・件数つき"
    );
    assert_eq!(summary.by_kind.len(), 3, "出現種別のみ（3 件）");
}

/// 受け入れ「by_kind は NonCentral を raw のまま別カウント（畳まない）かつ全種の Debug 昇順」。
/// 6 種すべて 1 件ずつ投入 → Annular<Hybrid<NonCentralAnnular<NonCentralTotal<Partial<Total。
/// NonCentralTotal/NonCentralAnnular が Total/Annular に**畳まれず**独立に現れることを縛る。
/// 殺す変異: by_kind で NonCentral を中心種別へ畳む（完備性 detected と取り違え）、Debug 昇順崩れ。
#[test]
fn by_kind_noncentral_kept_raw_and_full_sort() {
    let eclipses = vec![
        eclipse_with(SolarEclipseKind::Partial, 0.1, 0.40),
        eclipse_with(SolarEclipseKind::Total, 0.2, 1.05),
        eclipse_with(SolarEclipseKind::NonCentralTotal, 0.3, 1.00),
        eclipse_with(SolarEclipseKind::Hybrid, 0.4, 1.01),
        eclipse_with(SolarEclipseKind::NonCentralAnnular, 0.5, 0.99),
        eclipse_with(SolarEclipseKind::Annular, 0.6, 0.96),
    ];
    let summary = summarize_sweep(&eclipses, counts(0, 0, 0, 0));

    let kinds: Vec<SolarEclipseKind> = summary.by_kind.iter().map(|(k, _)| *k).collect();
    assert_eq!(
        kinds,
        vec![
            SolarEclipseKind::Annular,
            SolarEclipseKind::Hybrid,
            SolarEclipseKind::NonCentralAnnular,
            SolarEclipseKind::NonCentralTotal,
            SolarEclipseKind::Partial,
            SolarEclipseKind::Total,
        ],
        "6 種 raw を Debug 文字列昇順（NonCentral は畳まず別キー）"
    );
    // NonCentral が独立 1 件ずつ（中心種別へ畳まれていない）。
    assert_eq!(
        kind_count(&summary, SolarEclipseKind::NonCentralTotal),
        Some(1),
        "NonCentralTotal は raw で 1 件（Total に畳まない）"
    );
    assert_eq!(
        kind_count(&summary, SolarEclipseKind::NonCentralAnnular),
        Some(1),
        "NonCentralAnnular は raw で 1 件（Annular に畳まない）"
    );
    assert_eq!(
        kind_count(&summary, SolarEclipseKind::Total),
        Some(1),
        "Total は raw で 1 件（NonCentralTotal を取り込まない）"
    );
}

/// 受け入れ「by_kind の件数は同一種別の複数件を 1 グループへ集約する」。
/// Total×4・Annular×1 → Total に 4・Annular に 1。
/// 殺す変異: 件数を 1 に丸める（出現フラグ化）、別グループへ分散、件数の加算を or に取り違え。
#[test]
fn by_kind_aggregates_counts_same_kind() {
    let eclipses = vec![
        eclipse_with(SolarEclipseKind::Total, 0.1, 1.05),
        eclipse_with(SolarEclipseKind::Total, 0.2, 1.04),
        eclipse_with(SolarEclipseKind::Annular, 0.3, 0.95),
        eclipse_with(SolarEclipseKind::Total, 0.4, 1.03),
        eclipse_with(SolarEclipseKind::Total, 0.5, 1.02),
    ];
    let summary = summarize_sweep(&eclipses, counts(0, 0, 0, 0));
    assert_eq!(
        kind_count(&summary, SolarEclipseKind::Total),
        Some(4),
        "Total 4 件を集約"
    );
    assert_eq!(
        kind_count(&summary, SolarEclipseKind::Annular),
        Some(1),
        "Annular 1 件"
    );
    assert_eq!(summary.by_kind.len(), 2, "出現 2 種のみ");
}

// ============================================================
// gamma_abs（|gamma| の min/max/mean）
// ============================================================

/// 受け入れ「gamma_abs は |gamma| の min/max/mean（負 gamma も絶対値化される）」。
/// gamma = {-0.9, 0.3, -0.1, 0.7} → |gamma| = {0.9, 0.3, 0.1, 0.7}。
///   n=4・min=0.1・max=0.9・mean=(0.9+0.3+0.1+0.7)/4=0.5。
/// 殺す変異: 絶対値化の欠落（負を含めると min が負・max がずれる）、min/max 取り違え、mean の n 誤り。
#[test]
fn gamma_abs_is_absolute_min_max_mean() {
    let eclipses = vec![
        eclipse_with(SolarEclipseKind::Total, -0.9, 1.0),
        eclipse_with(SolarEclipseKind::Total, 0.3, 1.0),
        eclipse_with(SolarEclipseKind::Total, -0.1, 1.0),
        eclipse_with(SolarEclipseKind::Total, 0.7, 1.0),
    ];
    let summary = summarize_sweep(&eclipses, counts(0, 0, 0, 0));
    let g = &summary.gamma_abs;
    assert_eq!(g.n, 4, "n=4");
    assert!((g.min - 0.1).abs() < EPS, "min(|gamma|)=0.1, got {}", g.min);
    assert!((g.max - 0.9).abs() < EPS, "max(|gamma|)=0.9, got {}", g.max);
    assert!(
        (g.mean - 0.5).abs() < EPS,
        "mean(|gamma|)=0.5（絶対値化後）, got {}",
        g.mean
    );
}

// ============================================================
// magnitude（greatest.magnitude.0 の素の値 min/max/mean）
// ============================================================

/// 受け入れ「magnitude は greatest.magnitude.0 の素の値 min/max/mean（絶対値化しない）」。
/// magnitude = {0.40, 1.05, 0.95} → n=3・min=0.40・max=1.05・mean=(0.40+1.05+0.95)/3=0.8。
/// 殺す変異: 絶対値化の誤挿入、min/max 取り違え、gamma 列との取り違え（配線ミス）、mean の n 誤り。
#[test]
fn magnitude_is_raw_min_max_mean() {
    let eclipses = vec![
        eclipse_with(SolarEclipseKind::Partial, 0.1, 0.40),
        eclipse_with(SolarEclipseKind::Total, 0.2, 1.05),
        eclipse_with(SolarEclipseKind::Annular, 0.3, 0.95),
    ];
    let summary = summarize_sweep(&eclipses, counts(0, 0, 0, 0));
    let m = &summary.magnitude;
    assert_eq!(m.n, 3, "n=3");
    assert!((m.min - 0.40).abs() < EPS, "min=0.40, got {}", m.min);
    assert!((m.max - 1.05).abs() < EPS, "max=1.05, got {}", m.max);
    assert!((m.mean - 0.8).abs() < EPS, "mean=0.8, got {}", m.mean);
}

/// 受け入れ「gamma_abs と magnitude は別列（配線が独立）」。
/// gamma=-2.0（|·|=2.0）・magnitude=0.30 の 1 件 → gamma_abs.max=2.0, magnitude.max=0.30。
/// 殺す変異: gamma_abs と magnitude のソースを同一列に取り違える（両者が同値になる変異を殺す）。
#[test]
fn gamma_and_magnitude_are_distinct_sources() {
    let eclipses = vec![eclipse_with(SolarEclipseKind::Total, -2.0, 0.30)];
    let summary = summarize_sweep(&eclipses, counts(0, 0, 0, 0));
    assert!(
        (summary.gamma_abs.max - 2.0).abs() < EPS,
        "gamma_abs は |gamma|=2.0"
    );
    assert!(
        (summary.magnitude.max - 0.30).abs() < EPS,
        "magnitude は 0.30（gamma と独立）"
    );
}

// ============================================================
// completeness（NASA 4 区分への畳み込み・expected・all_match）
// ============================================================

/// 受け入れ「detected は NASA 4 区分へ畳む（NonCentralTotal→total, NonCentralAnnular→annular,
/// Hybrid→hybrid, Partial→partial）。各畳み込みを 1 件ずつ独立に縛る」。
/// 投入: Total×1・NonCentralTotal×1（→total=2）/ Annular×1・NonCentralAnnular×1（→annular=2）/
///       Hybrid×1（→hybrid=1）/ Partial×1（→partial=1）。
/// 殺す変異: NonCentralTotal を total に足さない、NonCentralAnnular を annular に足さない、
///   Hybrid/Partial の配線取り違え、区分間の取り違え（NonCentralTotal→annular 等）。
#[test]
fn completeness_detected_folds_into_four_categories() {
    let eclipses = vec![
        eclipse_with(SolarEclipseKind::Total, 0.1, 1.05),
        eclipse_with(SolarEclipseKind::NonCentralTotal, 0.2, 1.00),
        eclipse_with(SolarEclipseKind::Annular, 0.3, 0.95),
        eclipse_with(SolarEclipseKind::NonCentralAnnular, 0.4, 0.99),
        eclipse_with(SolarEclipseKind::Hybrid, 0.5, 1.01),
        eclipse_with(SolarEclipseKind::Partial, 0.6, 0.40),
    ];
    // expected は detected と完全一致させて all_match=true も同時に縛る。
    let expected = counts(2, 2, 1, 1);
    let summary = summarize_sweep(&eclipses, expected);
    let d = summary.completeness.detected;
    assert_eq!(d.total, 2, "total ← Total(1)+NonCentralTotal(1)");
    assert_eq!(d.annular, 2, "annular ← Annular(1)+NonCentralAnnular(1)");
    assert_eq!(d.hybrid, 1, "hybrid ← Hybrid(1)");
    assert_eq!(d.partial, 1, "partial ← Partial(1)");
    assert!(
        summary.completeness.all_match,
        "detected==expected → all_match=true"
    );
}

/// 受け入れ「completeness.expected は引数をそのまま保持する」。
/// detected と異なる expected を渡し、expected フィールドが引数の値そのものであることを縛る。
/// 殺す変異: expected を detected で上書き、expected の区分入れ替え、引数の取り違え。
#[test]
fn completeness_expected_is_passthrough() {
    let eclipses = vec![eclipse_with(SolarEclipseKind::Total, 0.1, 1.05)];
    let expected = counts(7, 5, 3, 9);
    let summary = summarize_sweep(&eclipses, expected);
    assert_eq!(
        summary.completeness.expected,
        counts(7, 5, 3, 9),
        "expected は引数そのまま"
    );
}

/// 受け入れ「all_match は detected==expected で true（全 4 区分一致）」。
/// detected=(1,1,1,1) になる 4 件を投入し expected も (1,1,1,1) → all_match=true。
/// 殺す変異: all_match を常に false/true 固定、比較の反転。
#[test]
fn all_match_true_when_all_four_equal() {
    let eclipses = vec![
        eclipse_with(SolarEclipseKind::Total, 0.1, 1.05),
        eclipse_with(SolarEclipseKind::Annular, 0.2, 0.95),
        eclipse_with(SolarEclipseKind::Hybrid, 0.3, 1.01),
        eclipse_with(SolarEclipseKind::Partial, 0.4, 0.40),
    ];
    let summary = summarize_sweep(&eclipses, counts(1, 1, 1, 1));
    assert!(
        summary.completeness.all_match,
        "4 区分一致 → all_match=true"
    );
}

/// 受け入れ「all_match は total 区分が 1 でもずれれば false」。
/// detected total=1（Total 1 件）だが expected total=2 → all_match=false。他 3 区分は一致。
/// 殺す変異: all_match の比較で total 区分を無視する。
#[test]
fn all_match_false_when_total_differs() {
    let eclipses = vec![
        eclipse_with(SolarEclipseKind::Total, 0.1, 1.05),
        eclipse_with(SolarEclipseKind::Annular, 0.2, 0.95),
        eclipse_with(SolarEclipseKind::Hybrid, 0.3, 1.01),
        eclipse_with(SolarEclipseKind::Partial, 0.4, 0.40),
    ];
    // detected=(1,1,1,1)・expected total だけ +1。
    let summary = summarize_sweep(&eclipses, counts(2, 1, 1, 1));
    assert!(
        !summary.completeness.all_match,
        "total 区分のずれ → all_match=false"
    );
}

/// 受け入れ「all_match は annular 区分が 1 でもずれれば false」。
/// 殺す変異: all_match の比較で annular 区分を無視する。
#[test]
fn all_match_false_when_annular_differs() {
    let eclipses = vec![
        eclipse_with(SolarEclipseKind::Total, 0.1, 1.05),
        eclipse_with(SolarEclipseKind::Annular, 0.2, 0.95),
        eclipse_with(SolarEclipseKind::Hybrid, 0.3, 1.01),
        eclipse_with(SolarEclipseKind::Partial, 0.4, 0.40),
    ];
    let summary = summarize_sweep(&eclipses, counts(1, 2, 1, 1));
    assert!(
        !summary.completeness.all_match,
        "annular 区分のずれ → all_match=false"
    );
}

/// 受け入れ「all_match は hybrid 区分が 1 でもずれれば false」。
/// 殺す変異: all_match の比較で hybrid 区分を無視する。
#[test]
fn all_match_false_when_hybrid_differs() {
    let eclipses = vec![
        eclipse_with(SolarEclipseKind::Total, 0.1, 1.05),
        eclipse_with(SolarEclipseKind::Annular, 0.2, 0.95),
        eclipse_with(SolarEclipseKind::Hybrid, 0.3, 1.01),
        eclipse_with(SolarEclipseKind::Partial, 0.4, 0.40),
    ];
    let summary = summarize_sweep(&eclipses, counts(1, 1, 2, 1));
    assert!(
        !summary.completeness.all_match,
        "hybrid 区分のずれ → all_match=false"
    );
}

/// 受け入れ「all_match は partial 区分が 1 でもずれれば false」。
/// 殺す変異: all_match の比較で partial 区分を無視する。
#[test]
fn all_match_false_when_partial_differs() {
    let eclipses = vec![
        eclipse_with(SolarEclipseKind::Total, 0.1, 1.05),
        eclipse_with(SolarEclipseKind::Annular, 0.2, 0.95),
        eclipse_with(SolarEclipseKind::Hybrid, 0.3, 1.01),
        eclipse_with(SolarEclipseKind::Partial, 0.4, 0.40),
    ];
    let summary = summarize_sweep(&eclipses, counts(1, 1, 1, 2));
    assert!(
        !summary.completeness.all_match,
        "partial 区分のずれ → all_match=false"
    );
}

// ============================================================
// 空入力（vacuous）
// ============================================================

/// 受け入れ「空 eclipses: total=0・by_kind 空 Vec・range n=0/全 0.0・detected 全 0。
/// expected が全 0 なら all_match=true」。
/// 殺す変異: 空でパニック、by_kind が非空、range が NaN/非 0、detected が非 0、空時 all_match を false 固定。
#[test]
fn empty_eclipses_all_zero_expected_zero_is_match() {
    let eclipses: Vec<SolarEclipse> = vec![];
    let summary = summarize_sweep(&eclipses, counts(0, 0, 0, 0));

    assert_eq!(summary.total, 0, "total=0");
    assert!(summary.by_kind.is_empty(), "by_kind 空 Vec");

    for (label, r) in [
        ("gamma_abs", &summary.gamma_abs),
        ("magnitude", &summary.magnitude),
    ] {
        assert_eq!(r.n, 0, "{label} は n=0");
        assert!(r.min == 0.0, "{label} 空 → min=0.0");
        assert!(r.max == 0.0, "{label} 空 → max=0.0");
        assert!(r.mean == 0.0, "{label} 空 → mean=0.0");
    }

    let d = summary.completeness.detected;
    assert_eq!(
        (d.total, d.annular, d.hybrid, d.partial),
        (0, 0, 0, 0),
        "detected 全 0"
    );
    assert!(
        summary.completeness.all_match,
        "空・expected 全 0 → all_match=true（vacuous）"
    );
}

/// 受け入れ「空 eclipses で expected が非 0 なら all_match=false（取りこぼし検知）」。
/// detected 全 0 vs expected total=1 → 不一致。
/// 殺す変異: 空入力で all_match を常に true 固定（expected を無視する）。
#[test]
fn empty_eclipses_nonzero_expected_is_mismatch() {
    let eclipses: Vec<SolarEclipse> = vec![];
    let summary = summarize_sweep(&eclipses, counts(1, 0, 0, 0));
    assert_eq!(summary.completeness.detected.total, 0, "detected total=0");
    assert!(
        !summary.completeness.all_match,
        "空・expected 非 0 → all_match=false"
    );
}

// ============================================================
// render_sweep_json / render_sweep_text
// ============================================================

/// 受け入れ「render_sweep_json は valid JSON（parse 可）＋末尾改行＋主要キー」。
/// 殺す変異: 末尾改行の欠落、非 JSON、主要キー（total/by_kind/gamma_abs/magnitude/completeness）の欠落。
#[test]
fn render_json_is_valid_with_trailing_newline() {
    let eclipses = vec![
        eclipse_with(SolarEclipseKind::Total, -0.3, 1.05),
        eclipse_with(SolarEclipseKind::Partial, 0.6, 0.40),
    ];
    let summary = summarize_sweep(&eclipses, counts(1, 0, 0, 1));
    let json = render_sweep_json_call(&summary).expect("JSON 直列化は成功");
    assert!(json.ends_with('\n'), "末尾改行あり");
    let parsed: serde_json::Value = serde_json::from_str(&json).expect("valid JSON");
    assert!(parsed.is_object(), "トップは object");
    assert!(parsed.get("total").is_some(), "total キーあり");
    assert!(parsed.get("by_kind").is_some(), "by_kind キーあり");
    assert!(parsed.get("gamma_abs").is_some(), "gamma_abs キーあり");
    assert!(parsed.get("magnitude").is_some(), "magnitude キーあり");
    assert!(
        parsed.get("completeness").is_some(),
        "completeness キーあり"
    );
}

/// 受け入れ「render_sweep_text に全項目の主要ラベル・数値が含まれる（完全一致は避ける）」。
/// 既知の集計（total=2・Total 1/Partial 1・gamma/magnitude・detected vs expected・all_match）を
/// 仕込み、各主要ラベルのサニティを縛る。
/// 殺す変異: total/by_kind/gamma/magnitude/completeness（detected/expected/all_match）の表示欠落。
#[test]
fn render_text_contains_all_sections() {
    let eclipses = vec![
        eclipse_with(SolarEclipseKind::Total, -0.3, 1.05),
        eclipse_with(SolarEclipseKind::Partial, 0.6, 0.40),
    ];
    let summary = summarize_sweep(&eclipses, counts(1, 0, 0, 1));
    let text = render_sweep_text_call(&summary);

    // total の表示。
    assert!(text.contains("total"), "total ラベルを表示: {text}");
    // by_kind の種別（Debug 文字列）。
    assert!(
        text.contains("Total"),
        "by_kind の種別 Total を表示: {text}"
    );
    assert!(
        text.contains("Partial"),
        "by_kind の種別 Partial を表示: {text}"
    );
    // gamma / magnitude のラベル。
    assert!(
        text.contains("gamma"),
        "gamma（gamma_abs）ラベルを表示: {text}"
    );
    assert!(text.contains("magnitude"), "magnitude ラベルを表示: {text}");
    // completeness の detected / expected / all_match。
    assert!(
        text.contains("detected"),
        "completeness の detected を表示: {text}"
    );
    assert!(
        text.contains("expected"),
        "completeness の expected を表示: {text}"
    );
    assert!(
        text.contains("all_match"),
        "completeness の all_match を表示: {text}"
    );
    // all_match の真偽値（detected≠expected を仕込んだので false が出るはず）。
    assert!(
        text.contains("true") || text.contains("false"),
        "all_match の判定（真偽）を表示: {text}"
    );
}

// ============================================================
// 呼び出しラッパ（未実装 API への薄い橋渡し。RED 時は import 不能で失敗）
// ============================================================

/// `render_sweep_json` 呼び出しラッパ。
fn render_sweep_json_call(summary: &SweepSummary) -> Result<String, serde_json::Error> {
    umbra_fixtures::render_sweep_json(summary)
}

/// `render_sweep_text` 呼び出しラッパ。
fn render_sweep_text_call(summary: &SweepSummary) -> String {
    umbra_fixtures::render_sweep_text(summary)
}

// `RangeStats` / `CompletenessReport` を import だけで終わらせずに 1 度は型として触れて
// 「未使用 import」警告を防ぎつつ、フィールドの存在を縛る（コンパイル時の型チェック）。
#[test]
fn types_are_constructible_and_fielded() {
    let r = RangeStats {
        n: 1,
        min: -1.0,
        max: 2.0,
        mean: 0.5,
    };
    assert_eq!(r.n, 1);
    assert!((r.min - -1.0).abs() < EPS);
    assert!((r.max - 2.0).abs() < EPS);
    assert!((r.mean - 0.5).abs() < EPS);

    let report = CompletenessReport {
        detected: counts(1, 2, 3, 4),
        expected: counts(1, 2, 3, 4),
        all_match: true,
    };
    assert_eq!(report.detected, report.expected);
    assert!(report.all_match);
}
