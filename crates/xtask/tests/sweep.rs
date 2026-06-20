//! 全食スイープ実走ランナー受け入れテスト（strict / `cargo xtask sweep`）。
//!
//! 本ファイルは `xtask` crate の **公開 API のみ**を対象とした統合テスト（tests/ 配下）。
//! 対象は本スライスが導入する `xtask::sweep` モジュール:
//! - `parse_range`（`--from`/`--to` 西暦年 → 各年 1/1 0h UTC の [start, end) 対・既定 1900/2100・
//!   非整数/逆転は `XtaskError::InvalidArgument`）。
//! - `parse_expected_counts`（`--expected-total/annular/hybrid/partial` → `CatalogCounts`・既定全 0・
//!   非整数は `XtaskError::InvalidArgument`）。
//! - `sweep_report`（検出食群を `summarize_sweep` → `render_sweep_text`/`render_sweep_json`・
//!   エンジン非依存の純粋寄り結線）。
//! - `run_sweep`（実エンジン実走の薄い印字部分。SLOW 経路は `sweep_report`＋実 `search` で縛る）。
//!
//! ## テスト戦略（mutation-resistant・負荷配分）
//! FAST はすべて合成データ・引数解釈で縛る（実エンジン非実走）。`sweep_report` は種別・gamma・
//! magnitude を既知値で制御した合成 `SolarEclipse` 群（`crates/umbra-fixtures/tests/sweep.rs` の
//! `eclipse_with` ヘルパをミラー）で集計結果を予言する。`parse_range`/`parse_expected_counts` は
//! `validate.rs` の `parse_format` 流儀（`--flag value`）をミラーした引数列で縛る。
//! SLOW は実 `search` を **狭い窓 1 つ**（2024-03-15〜2024-05-01・2024-04-08 の皆既 1 件）に限定。
//!
//! ## 期待される RED（実装前）
//! `xtask::sweep` モジュールも `parse_range` / `parse_expected_counts` / `sweep_report` /
//! `run_sweep` も存在しないため、本ファイルは **未解決インポート/モジュール（E0432/E0433）で
//! コンパイル不能 = RED** になる。これが想定どおりの赤。

#![allow(clippy::excessive_precision)]

use serde_json::Value;

use umbra_core::{Degrees, JulianDate2, Kilometers, TtInstant, UtcInstant};
use umbra_eclipse::{
    standard_engine, AccuracyProfile, BesselFitError, BesselianPolynomial, CalculationMetadata,
    EclipseEngine, EclipseMagnitude, GlobalCircumstances, GlobalContact, GreatestEclipse,
    Obscuration, Polynomial, SolarEclipse, SolarEclipseKind, UtcRange,
};
use umbra_ephemeris::bundled_time_data;
use umbra_geo::GeoPoint;

use umbra_fixtures::CatalogCounts;

use xtask::sweep::{parse_expected_counts, parse_range, sweep_report};
use xtask::validate::ReportFormat;

/// 1 日の秒数（JD 差 → 秒の換算・年境界の時刻一致確認に使う）。
const SECONDS_PER_DAY: f64 = 86_400.0;

/// 年境界 jd2 一致の許容（秒換算）。gregorian→jd2 は決定的なので極小で縛る。
const TIME_EPS_S: f64 = 1e-6;

/// 範囲統計の f64 一致に用いる厳密許容。
const EPS: f64 = 1e-9;

// ============================================================
// 引数列ヘルパ
// ============================================================

/// &str スライス → Vec<String>（parse_range/parse_expected_counts の引数列構築ヘルパ）。
fn args(parts: &[&str]) -> Vec<String> {
    parts.iter().map(|s| s.to_string()).collect()
}

/// UTC 瞬時を整数引数で組む（年境界期待値・合成食の時刻に使う）。
fn utc(year: i32, month: u8, day: u8, hour: u8, minute: u8, second: f64) -> UtcInstant {
    UtcInstant::from_gregorian(year, month, day, hour, minute, second).expect("有効な UTC 日時")
}

/// 2 つの UTC 瞬時の差（秒・桁落ち回避の days_since 経由）。年境界一致の検証に使う。
fn utc_diff_seconds(a: UtcInstant, b: UtcInstant) -> f64 {
    a.jd2().days_since(b.jd2()) * SECONDS_PER_DAY
}

// ============================================================
// 合成 SolarEclipse 構築ヘルパ（umbra-fixtures/tests/sweep.rs の eclipse_with をミラー）
// ============================================================

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

/// 合成 `SolarEclipse`。集計が読む 3 値（kind / global.gamma / greatest.magnitude）を引数で
/// 指定し、残りは固定フィラー。時刻はスイープ集計に無関係なので固定。
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

// ============================================================
// parse_range（FAST・純引数解釈）
// ============================================================

/// 受け入れ「フラグ無し → 既定 (1900-01-01, 2100-01-01) の各年 1/1 0h UTC」。
/// 殺す変異: 既定値の取り違え（from/to の年数・1900/2100 以外）、1/1 0h 以外の境界、start/end 入れ替え。
#[test]
fn parse_range_defaults_to_1900_2100() {
    let (start, end) = parse_range(&[]).expect("既定は Ok");
    assert!(
        utc_diff_seconds(start, utc(1900, 1, 1, 0, 0, 0.0)).abs() < TIME_EPS_S,
        "既定 start は 1900-01-01 0h UTC"
    );
    assert!(
        utc_diff_seconds(end, utc(2100, 1, 1, 0, 0, 0.0)).abs() < TIME_EPS_S,
        "既定 end は 2100-01-01 0h UTC"
    );
}

/// 受け入れ「`--from 2000 --to 2050` → (2000-01-01, 2050-01-01) の各年 1/1 0h UTC」。
/// 殺す変異: from/to を読まず既定にフォールバック、from↔to 取り違え、年→境界変換のオフセット混入。
#[test]
fn parse_range_reads_from_and_to() {
    let (start, end) =
        parse_range(&args(&["--from", "2000", "--to", "2050"])).expect("整数年は Ok");
    assert!(
        utc_diff_seconds(start, utc(2000, 1, 1, 0, 0, 0.0)).abs() < TIME_EPS_S,
        "start は 2000-01-01 0h UTC"
    );
    assert!(
        utc_diff_seconds(end, utc(2050, 1, 1, 0, 0, 0.0)).abs() < TIME_EPS_S,
        "end は 2050-01-01 0h UTC"
    );
}

/// 受け入れ「`--from` だけ指定なら `--to` は既定 2100 を保つ（部分指定の独立性）」。
/// 殺す変異: 片方指定で他方を巻き添えにする（両方 from にする等）、既定 to を無視する。
#[test]
fn parse_range_partial_keeps_other_default() {
    let (start, end) = parse_range(&args(&["--from", "1950"])).expect("from のみは Ok");
    assert!(
        utc_diff_seconds(start, utc(1950, 1, 1, 0, 0, 0.0)).abs() < TIME_EPS_S,
        "start は 1950-01-01 0h UTC"
    );
    assert!(
        utc_diff_seconds(end, utc(2100, 1, 1, 0, 0, 0.0)).abs() < TIME_EPS_S,
        "end は既定 2100-01-01 0h UTC（from だけ指定でも保つ）"
    );
}

/// 受け入れ「非整数の `--from` は InvalidArgument{flag:--from, value:abc}」。
/// 殺す変異: 非整数を 0 等へ握り潰す、エラー種別の取り違え、flag/value フィールドの取り違え。
#[test]
fn parse_range_non_integer_from_is_error() {
    let r = parse_range(&args(&["--from", "abc"]));
    assert!(
        matches!(
            r,
            Err(xtask::error::XtaskError::InvalidArgument { ref flag, ref value })
                if flag == "--from" && value == "abc"
        ),
        "非整数 --from は InvalidArgument{{flag:--from, value:abc}}, got {r:?}"
    );
}

/// 受け入れ「非整数の `--to` も InvalidArgument{flag:--to, ..}」（from とは別フラグ名）。
/// 殺す変異: --to の検証欠落（既定へフォールバック）、flag 名を --from に取り違える。
#[test]
fn parse_range_non_integer_to_is_error() {
    let r = parse_range(&args(&["--to", "zzz"]));
    assert!(
        matches!(
            r,
            Err(xtask::error::XtaskError::InvalidArgument { ref flag, ref value })
                if flag == "--to" && value == "zzz"
        ),
        "非整数 --to は InvalidArgument{{flag:--to, value:zzz}}, got {r:?}"
    );
}

/// 受け入れ「from>to（範囲逆転）は InvalidArgument（Ok にしてはならない）」。
/// エラーの flag/value 文字列は実装裁量（逆転は from/to のどちらに帰属させても可）ゆえ問わない。
/// 殺す変異: 逆転を許して空/負の範囲を返す、検証を欠落させる。
#[test]
fn parse_range_from_greater_than_to_is_error() {
    let r = parse_range(&args(&["--from", "2100", "--to", "2000"]));
    assert!(
        matches!(r, Err(xtask::error::XtaskError::InvalidArgument { .. })),
        "from>to は InvalidArgument（逆転を Ok にしてはならない）, got {r:?}"
    );
}

/// 受け入れ「from==to は Ok（空区間 [y, y) は逆転ではない・境界 inclusive 下端）」。
/// 殺す変異: 比較を `>=` にして等値を誤って弾く（境界の取り違え）。
#[test]
fn parse_range_from_equals_to_is_ok() {
    let (start, end) = parse_range(&args(&["--from", "2000", "--to", "2000"]))
        .expect("from==to は空区間として Ok");
    assert!(
        utc_diff_seconds(start, end).abs() < TIME_EPS_S,
        "from==to は start==end（空区間）"
    );
}

// ============================================================
// parse_expected_counts（FAST・純引数解釈）
// ============================================================

/// 受け入れ「フラグ無し → 既定 (0,0,0,0)」。
/// 殺す変異: 既定を非 0 にする、どれかのフィールドを別既定にする。
#[test]
fn parse_expected_counts_defaults_to_zero() {
    let c = parse_expected_counts(&[]).expect("既定は Ok");
    assert_eq!(c, counts(0, 0, 0, 0), "フラグ無しは全 0");
}

/// 受け入れ「4 フラグ全指定 → CatalogCounts{5,3,1,7}（各フラグが正しい区分へ通る）」。
/// 殺す変異: total/annular/hybrid/partial の配線取り違え、フラグ名比較の改変、値の読み落とし。
#[test]
fn parse_expected_counts_reads_all_four() {
    let c = parse_expected_counts(&args(&[
        "--expected-total",
        "5",
        "--expected-annular",
        "3",
        "--expected-hybrid",
        "1",
        "--expected-partial",
        "7",
    ]))
    .expect("4 フラグ全指定は Ok");
    assert_eq!(
        c,
        counts(5, 3, 1, 7),
        "total=5 annular=3 hybrid=1 partial=7 が各区分へ正しく通る"
    );
}

/// 受け入れ「一部だけ指定なら他区分は 0 を保つ（部分指定の独立性）」。
/// `--expected-annular 4` のみ → annular=4・他 3 区分は 0。
/// 殺す変異: 未指定区分を非 0 既定にする、指定値を全区分へ波及させる。
#[test]
fn parse_expected_counts_partial_keeps_others_zero() {
    let c = parse_expected_counts(&args(&["--expected-annular", "4"])).expect("一部指定は Ok");
    assert_eq!(
        c,
        counts(0, 4, 0, 0),
        "annular=4・他区分は 0（部分指定の独立性）"
    );
}

/// 受け入れ「非整数の期待件数は InvalidArgument」。`--expected-total foo` → Err。
/// flag/value の文字列は実装裁量だが、Ok にしてはならない。
/// 殺す変異: 非整数を 0 へ握り潰す、検証欠落、エラー種別の取り違え。
#[test]
fn parse_expected_counts_non_integer_is_error() {
    let r = parse_expected_counts(&args(&["--expected-total", "foo"]));
    assert!(
        matches!(r, Err(xtask::error::XtaskError::InvalidArgument { .. })),
        "非整数の期待件数は InvalidArgument（Ok にしてはならない）, got {r:?}"
    );
}

// ============================================================
// sweep_report（FAST・合成 eclipses で縛る）
// ============================================================

/// 受け入れ「Text 形式は summarize_sweep→render_sweep_text の人間可読サマリを返す」。
/// render_sweep_text は total/by_kind 種別/gamma/magnitude/completeness ラベルを出す（含有サニティ）。
/// 合成: Total×1・Partial×1（total=2・by_kind に Total/Partial）。脆い完全一致は避ける。
/// 殺す変異: Json 経路へ分岐する、render_sweep_text を呼ばない、空文字列を返す、集計を通さない。
#[test]
fn sweep_report_text_contains_all_sections() {
    let eclipses = vec![
        eclipse_with(SolarEclipseKind::Total, -0.3, 1.05),
        eclipse_with(SolarEclipseKind::Partial, 0.6, 0.40),
    ];
    let s =
        sweep_report(&eclipses, counts(1, 0, 0, 1), ReportFormat::Text).expect("Text レンダは Ok");
    assert!(s.contains("total"), "total ラベル: {s}");
    assert!(s.contains("Total"), "by_kind の種別 Total: {s}");
    assert!(s.contains("Partial"), "by_kind の種別 Partial: {s}");
    assert!(s.contains("gamma"), "gamma_abs ラベル: {s}");
    assert!(s.contains("magnitude"), "magnitude ラベル: {s}");
    assert!(s.contains("completeness"), "completeness ラベル: {s}");
}

/// 受け入れ「Json 形式は妥当な JSON・末尾改行・主要キーを持ち、集計値が正しく通る」。
/// 合成: Total×2・Annular×1・Partial×1（total=4）。gamma を制御して gamma_abs.max を予言。
/// 殺す変異: Text 経路へ分岐する、render_sweep_json を呼ばない（末尾改行欠落）、total を取り違える、
///   summarize_sweep に誤データを渡す。
#[test]
fn sweep_report_json_is_object_with_keys_and_counts() {
    let eclipses = vec![
        eclipse_with(SolarEclipseKind::Total, 0.1, 1.05),
        eclipse_with(SolarEclipseKind::Total, -0.9, 1.02),
        eclipse_with(SolarEclipseKind::Annular, 0.2, 0.95),
        eclipse_with(SolarEclipseKind::Partial, 0.3, 0.40),
    ];
    let s =
        sweep_report(&eclipses, counts(2, 1, 0, 1), ReportFormat::Json).expect("Json レンダは Ok");

    assert!(s.ends_with('\n'), "Json は末尾改行で終わる: {s:?}");
    let v: Value = serde_json::from_str(&s).expect("出力は妥当な JSON");
    assert!(v.is_object(), "トップレベルは JSON オブジェクト: {s}");
    let obj = v.as_object().expect("オブジェクト");
    for key in ["total", "by_kind", "gamma_abs", "magnitude", "completeness"] {
        assert!(obj.contains_key(key), "トップレベルキー {key} が存在: {s}");
    }
    assert_eq!(v["total"], 4, "total==4（eclipses.len()）");
    let gmax = v["gamma_abs"]["max"]
        .as_f64()
        .expect("gamma_abs.max は数値");
    assert!(
        (gmax - 0.9).abs() < EPS,
        "gamma_abs.max==|−0.9|=0.9（合成 gamma を通している）, got {gmax}"
    );
}

/// 受け入れ「expected を渡すと completeness.expected と all_match が反映される（一致で true）」。
/// 合成 detected=(1,1,1,1) に expected=(1,1,1,1) → all_match=true。
/// 殺す変異: expected を detected で上書き、all_match を false 固定、summarize へ expected を渡さない。
#[test]
fn sweep_report_json_completeness_all_match_true_when_expected_matches() {
    let eclipses = vec![
        eclipse_with(SolarEclipseKind::Total, 0.1, 1.05),
        eclipse_with(SolarEclipseKind::Annular, 0.2, 0.95),
        eclipse_with(SolarEclipseKind::Hybrid, 0.3, 1.01),
        eclipse_with(SolarEclipseKind::Partial, 0.4, 0.40),
    ];
    let s =
        sweep_report(&eclipses, counts(1, 1, 1, 1), ReportFormat::Json).expect("Json レンダは Ok");
    let v: Value = serde_json::from_str(&s).expect("妥当な JSON");
    let comp = &v["completeness"];
    assert_eq!(comp["expected"]["total"], 1, "expected.total は引数 1");
    assert_eq!(comp["expected"]["annular"], 1, "expected.annular は引数 1");
    assert_eq!(comp["expected"]["hybrid"], 1, "expected.hybrid は引数 1");
    assert_eq!(comp["expected"]["partial"], 1, "expected.partial は引数 1");
    assert_eq!(
        comp["all_match"], true,
        "detected==expected=(1,1,1,1) → all_match=true"
    );
}

/// 受け入れ「expected が detected とずれれば all_match=false（完備性突合が反映される）」。
/// 合成 detected total=1（Total 1 件）だが expected total=2 → all_match=false。
/// 殺す変異: all_match を true 固定、expected を無視して常に一致扱いにする。
#[test]
fn sweep_report_json_completeness_all_match_false_when_expected_differs() {
    let eclipses = vec![eclipse_with(SolarEclipseKind::Total, 0.1, 1.05)];
    let s =
        sweep_report(&eclipses, counts(2, 0, 0, 0), ReportFormat::Json).expect("Json レンダは Ok");
    let v: Value = serde_json::from_str(&s).expect("妥当な JSON");
    assert_eq!(
        v["completeness"]["all_match"], false,
        "detected total=1 ≠ expected total=2 → all_match=false"
    );
}

// ============================================================
// SLOW 統合テスト（実エンジン・狭い窓 1 つ・2024-04-08 皆既 1 件）
// ============================================================

/// 【SLOW・1 件】解析暦エンジン（`standard_engine(bundled_time_data())`）で **狭い窓**
/// （2024-03-15〜2024-05-01・2024-04-08 の皆既を含む）を実 `search` → `sweep_report(Json)` し、
/// 妥当 JSON・total>=1・by_kind に Total を含むことを縛る。`run_sweep` 自体は stdout 印字で戻り値が
/// `()` のため直接呼ばず、`sweep_report`＋実 `search` の経路で結線を担保する（薄い印字部分は手動確認）。
/// 実 `search` は重いので **窓 1 つ**に限定し実行時間を抑える。
/// 殺す変異: sweep_report が search 結果を集計しない、Json レンダ未配線、total を通さない、
///   by_kind に検出種別を載せない。
// SLOW
#[test]
fn sweep_report_real_engine_narrow_window_json() {
    let engine: EclipseEngine<_, _, _> = standard_engine(bundled_time_data());
    let range = UtcRange {
        start: utc(2024, 3, 15, 0, 0, 0.0),
        end: utc(2024, 5, 1, 0, 0, 0.0),
    };
    let eclipses = engine.search(range).expect("狭い窓の search は Ok");

    // 2024-04-08 の皆既を含む窓ゆえ少なくとも 1 件は検出される（期待件数は未指定＝全 0）。
    let s = sweep_report(&eclipses, counts(0, 0, 0, 0), ReportFormat::Json)
        .expect("実 search 結果の Json レンダは Ok");

    let v: Value = serde_json::from_str(&s).expect("出力は妥当な JSON");
    let total = v["total"].as_u64().expect("total は整数");
    assert!(
        total >= 1,
        "2024-04-08 皆既を含む窓 → total>=1, got {total}"
    );

    // by_kind は配列で、Total 種別（皆既）を 1 種以上含む。
    let by_kind = v["by_kind"].as_array().expect("by_kind は配列");
    let has_total = by_kind.iter().any(|pair| {
        pair.as_array()
            .and_then(|kv| kv.first())
            .and_then(|k| k.as_str())
            == Some("Total")
    });
    assert!(
        has_total,
        "2024-04-08 は皆既 → by_kind に Total を含む: {s}"
    );
}
