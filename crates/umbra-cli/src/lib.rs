//! `umbra` CLI ライブラリ（ISSUE-031 `umbra search`）。
//!
//! 薄い CLI ラッパ: 引数解釈（clap）・日付パース・`EclipseEngine::search` 呼び出し・整形出力。
//! 計算は umbra-eclipse が担保。本クレートは境界（引数・パース・出力・エラー/終了コード）が責務。
//!
//! S31a 範囲: `umbra search`（text 出力）。`--format json`（serde 横断配線）は後続スライス。

use clap::{Args, Parser, Subcommand, ValueEnum};
use umbra_core::{jd2_to_gregorian, EspenakMeeusDeltaT, UtcInstant};
use umbra_eclipse::{
    standard_engine, EclipseEngine, EclipseError, EngineConfig, SolarEclipse, SolarEclipseKind,
    UtcRange,
};
use umbra_ephemeris::{bundled_time_data, AnalyticalEphemeris};

/// `umbra` CLI ルート。
#[derive(Debug, Parser)]
#[command(name = "umbra", about = "umbra-rs solar eclipse CLI", version)]
pub struct Cli {
    /// 実行するサブコマンド。
    #[command(subcommand)]
    pub command: Command,
}

/// サブコマンド（S31a は `search` のみ。local/path 等は後続 issue）。
#[derive(Debug, Subcommand)]
pub enum Command {
    /// 期間内の太陽食を列挙する（`EclipseEngine::search`）。
    Search(SearchArgs),
}

/// `umbra search` の引数。
#[derive(Debug, Args)]
pub struct SearchArgs {
    /// 開始日（`YYYY-MM-DD`, UTC・境界含む）。
    #[arg(long)]
    pub from: String,
    /// 終了日（`YYYY-MM-DD`, UTC・境界含む）。
    #[arg(long)]
    pub to: String,
    /// 精度プロファイル（既定 standard）。
    #[arg(long, value_enum, default_value_t = AccuracyArg::Standard)]
    pub accuracy: AccuracyArg,
    /// 種別フィルタ（既定 all）。
    #[arg(long, value_enum, default_value_t = KindFilter::All)]
    pub kind: KindFilter,
}

/// 精度プロファイル引数（公開 2 層, api-draft §3.1）。
#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
pub enum AccuracyArg {
    /// 標準（既定）。
    Standard,
    /// 参照（高精度・低速）。
    Reference,
}

/// 種別フィルタ引数。
#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
pub enum KindFilter {
    /// すべて。
    All,
    /// 皆既のみ（非中心皆既含む）。
    Total,
    /// 金環のみ（非中心金環含む）。
    Annular,
    /// 部分のみ。
    Partial,
    /// ハイブリッドのみ。
    Hybrid,
}

/// CLI 境界のエラー（終了コード非 0・メッセージ用）。
#[derive(Debug, thiserror::Error)]
pub enum CliError {
    /// 日付パース失敗（`YYYY-MM-DD` 以外・暦範囲外など）。
    #[error("invalid date '{0}' (expected YYYY-MM-DD)")]
    InvalidDate(String),
    /// `--from` が `--to` より後。
    #[error("--from ({from}) must be on or before --to ({to})")]
    RangeOrder {
        /// 開始日（入力文字列）。
        from: String,
        /// 終了日（入力文字列）。
        to: String,
    },
    /// エンジン側エラー（探索・時刻系・暦・範囲外など, 透過）。
    #[error(transparent)]
    Eclipse(#[from] EclipseError),
}

/// `YYYY-MM-DD`（UTC 0:00:00）を [`UtcInstant`] にパースする。
///
/// `-` 区切りちょうど 3 要素・各要素が 10 進整数で、暦として妥当（月/日範囲）なら成功。
/// それ以外は [`CliError::InvalidDate`]（入力文字列をそのまま保持）。負年・RFC3339 は S31a 非対応。
pub fn parse_date(text: &str) -> Result<UtcInstant, CliError> {
    let invalid = || CliError::InvalidDate(text.to_string());
    let parts: Vec<&str> = text.split('-').collect();
    if parts.len() != 3 {
        return Err(invalid());
    }
    let year: i32 = parts[0].parse().map_err(|_| invalid())?;
    let month: u8 = parts[1].parse().map_err(|_| invalid())?;
    let day: u8 = parts[2].parse().map_err(|_| invalid())?;
    // 月/日範囲などの暦妥当性は from_gregorian（calendar）が検証する（month∉1..=12 は OutOfRange）。
    UtcInstant::from_gregorian(year, month, day, 0, 0, 0.0).map_err(|_| invalid())
}

/// 種別フィルタが日食種別に合致するか（`All` は常に true）。
///
/// 皆既/金環フィルタは**非中心**（[`SolarEclipseKind::NonCentralTotal`]/`NonCentralAnnular`）を含む
/// （確定仕様, ISSUE-031）。`Partial`/`Hybrid` はそれぞれ単一種別のみ。
pub fn kind_matches(filter: KindFilter, kind: SolarEclipseKind) -> bool {
    use SolarEclipseKind::{Annular, Hybrid, NonCentralAnnular, NonCentralTotal, Partial, Total};
    match filter {
        KindFilter::All => true,
        KindFilter::Total => matches!(kind, Total | NonCentralTotal),
        KindFilter::Annular => matches!(kind, Annular | NonCentralAnnular),
        KindFilter::Partial => matches!(kind, Partial),
        KindFilter::Hybrid => matches!(kind, Hybrid),
    }
}

/// 日食リストを人間可読 text に整形する。各日食の event_key・種別・最大食 **UTC と TT の両方**
/// （accuracy.md §0）・gamma・食分・食面積＋計算メタデータ（暦/ΔT モデル名・ΔT 不確実性帯）。
/// 空リストは「該当なし」を表す 1 行（具体的日食情報を含まない・非 panic）。
pub fn format_search_text(eclipses: &[SolarEclipse]) -> String {
    if eclipses.is_empty() {
        return "No solar eclipses found in the given range.\n".to_string();
    }
    let mut out = String::new();
    for e in eclipses {
        let g = &e.global.greatest;
        let (y, mo, d, h, mi, s) = g.time_utc.to_gregorian();
        // TT も暦形式で併記（accuracy.md §0: UTC+TT 両方）。TT-JD をグレゴリオ暦に直す。
        let (ty, tmo, td, th, tmi, ts) = jd2_to_gregorian(g.time_tt.jd2());
        let m = &e.metadata;
        out.push_str(&format!(
            "{key}  {kind:?}\n",
            key = e.event_key,
            kind = e.kind
        ));
        out.push_str(&format!(
            "  greatest: {y:04}-{mo:02}-{d:02} {h:02}:{mi:02}:{s:04.1} UTC  /  \
             {ty:04}-{tmo:02}-{td:02} {th:02}:{tmi:02}:{ts:04.1} TT\n",
        ));
        out.push_str(&format!(
            "  gamma: {gamma:.4}  magnitude: {mag:.4}  obscuration: {obsc:.4}\n",
            gamma = e.global.gamma,
            mag = g.magnitude.0,
            obsc = g.obscuration.0,
        ));
        out.push_str(&format!(
            "  ephemeris: {em} {ev}  ΔT: {dt} (±{unc:.2}s)  accuracy: {acc:?}\n",
            em = m.ephemeris_model,
            ev = m.ephemeris_version,
            dt = m.delta_t_model,
            unc = m.delta_t_uncertainty_seconds,
            acc = m.accuracy_profile,
        ));
    }
    out
}

/// `umbra search` を実行し、整形済み text 出力を返す（日付パース→エンジン構築→search→
/// kind フィルタ→整形）。出力は呼び出し側（main）が印字する（テスト容易性のため String 返し）。
///
/// 不正日付・`from > to` は **search を呼ぶ前に** fast-fail（[`CliError::InvalidDate`]/
/// [`CliError::RangeOrder`]）。エンジンは同梱データ（[`bundled_time_data`]・実行時ネットワークなし）。
pub fn run_search(args: &SearchArgs) -> Result<String, CliError> {
    let from = parse_date(&args.from)?;
    let to = parse_date(&args.to)?;
    if from.jd2().jd() > to.jd2().jd() {
        return Err(CliError::RangeOrder {
            from: args.from.clone(),
            to: args.to.clone(),
        });
    }

    let time = bundled_time_data();
    let range = UtcRange {
        start: from,
        end: to,
    };
    // 精度プロファイルでエンジンを構築（Standard は standard_engine ショートカット、Reference は
    // 同梱 EOP を複製して reference config で直接構築。いずれも StandardEngine 型）。
    let eclipses = match args.accuracy {
        AccuracyArg::Standard => standard_engine(time).search(range)?,
        AccuracyArg::Reference => {
            let earth_orientation = time.eop().clone();
            EclipseEngine::new(
                AnalyticalEphemeris::new(),
                EspenakMeeusDeltaT,
                earth_orientation,
                time,
                EngineConfig::reference(),
            )
            .search(range)?
        }
    };

    let filtered: Vec<SolarEclipse> = eclipses
        .into_iter()
        .filter(|eclipse| kind_matches(args.kind, eclipse.kind))
        .collect();
    Ok(format_search_text(&filtered))
}

#[cfg(test)]
mod tests {
    //! ISSUE-031 S31a 受け入れテスト（standard・`umbra search` text 出力）。
    //!
    //! ## オラクル戦略（実装方針に立ち入らず、確定仕様の公開 IF だけを縛る）
    //! - **parse_date**: 期待値を検証済みプリミティブ `UtcInstant::from_gregorian` で独立に組む
    //!   （実装の文字列処理を写経しない）。妥当日付は round-trip 一致、明確な不正は `InvalidDate`。
    //! - **kind_matches**: 仕様の真理値表（皆既/金環フィルタは非中心を含む）を手で記述。
    //! - **format_search_text**: 既知 fixture（皆既/部分の 2 件）の値が出力に出る存在確認
    //!   （厳密整形でなく内容＝部分文字列）。空リストは非 panic で「該当なし」を表す。
    //! - **run_search**: fast-fail 2 件（from>to／不正 from は search を呼ぶ前に Err）＝高速。
    //!   正常系 1 件のみ SLOW（2017-08 実日食 search を実走）。
    //!
    //! ## red 設計（本体未実装）
    //! `parse_date`/`kind_matches`/`format_search_text`/`run_search` は本体 `unimplemented!`。
    //! テストは戻り値・出力内容を要求するため `unimplemented!` の panic で red。fast-fail テストも
    //! 入口の `unimplemented!` で panic する（red 段階では search 非実走＝速い）。

    #![allow(clippy::excessive_precision)]

    use super::*;

    use umbra_core::{Degrees, JulianDate2, Kilometers, TimeInterval, TtInstant, UtcInstant};
    use umbra_eclipse::{
        AccuracyProfile, BesselFitError, BesselianPolynomial, CalculationMetadata,
        EclipseMagnitude, GlobalCircumstances, GlobalContact, GreatestEclipse, Obscuration,
        Polynomial, SolarEclipse, SolarEclipseKind,
    };
    use umbra_geo::GeoPoint;

    // ------------------------------------------------------------------
    // 構築ヘルパ（engine.rs / results.rs の minimal_* パターンを再掲）
    // ------------------------------------------------------------------

    /// UTC 瞬時を整数引数で組む。
    fn utc(year: i32, month: u8, day: u8, hour: u8, minute: u8, second: f64) -> UtcInstant {
        UtcInstant::from_gregorian(year, month, day, hour, minute, second).expect("有効な UTC 日時")
    }

    /// TT 瞬時を 2 要素 JD で組む（UTC と区別できる別スケール値）。
    fn tt(jd1: f64, jd2: f64) -> TtInstant {
        TtInstant::from_jd2(JulianDate2::new(jd1, jd2))
    }

    /// 地表点（lat, lon）を度から組む。
    fn geo(lat: f64, lon: f64) -> GeoPoint {
        GeoPoint::from_degrees(lat, lon).expect("有効な地表点")
    }

    /// 全球接触点（時刻 TT/UTC ＋ 地表点）。
    fn contact(h: u8, lat: f64) -> GlobalContact {
        GlobalContact {
            time_utc: utc(2024, 4, 8, h, 0, 0.0),
            time_tt: tt(2_460_409.0, f64::from(h) * 0.01),
            position: geo(lat, -100.0),
        }
    }

    /// 最小 BesselianPolynomial（results.rs の minimal_bessel パターン）。
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

    /// 既知の計算メタデータ（レシピ全フィールド非空・format の存在確認に使う識別文字列）。
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

    /// 皆既（中心食）の fixture SolarEclipse（format の主力・既知値で内容を縛る）。
    /// event_key・最大食 UTC/TT・gamma・magnitude・obscuration を **互いに異なる既知値**で持つ。
    fn total_eclipse() -> SolarEclipse {
        let greatest = GreatestEclipse {
            time_utc: utc(2024, 4, 8, 18, 17, 0.0),
            time_tt: tt(2_460_409.0, 0.123),
            position: geo(25.0, -104.0),
            magnitude: EclipseMagnitude(1.0566),
            obscuration: Obscuration(1.0),
            path_width: Some(Kilometers(197.0)),
            central_duration: Some(268.0),
            sun_altitude: Degrees(70.3),
        };
        let global = GlobalCircumstances {
            kind: SolarEclipseKind::Total,
            partial_begin: Some(contact(16, 10.0)),
            central_begin: Some(contact(17, 20.0)),
            greatest,
            central_end: Some(contact(19, 40.0)),
            partial_end: Some(contact(20, 50.0)),
            gamma: 0.3431,
        };
        SolarEclipse {
            event_key: "2024-04-08#1252".to_string(),
            kind: SolarEclipseKind::Total,
            global,
            bessel: minimal_bessel(),
            metadata: metadata(),
        }
    }

    /// 部分食の fixture SolarEclipse（中心食用フィールドが None・magnitude<1）。
    fn partial_eclipse() -> SolarEclipse {
        let greatest = GreatestEclipse {
            time_utc: utc(2025, 3, 29, 10, 47, 0.0),
            time_tt: tt(2_460_763.0, 0.456),
            position: geo(61.0, -77.0),
            magnitude: EclipseMagnitude(0.938),
            obscuration: Obscuration(0.91),
            path_width: None,
            central_duration: None,
            sun_altitude: Degrees(12.0),
        };
        let global = GlobalCircumstances {
            kind: SolarEclipseKind::Partial,
            partial_begin: Some(contact(10, 60.0)),
            central_begin: None,
            greatest,
            central_end: None,
            partial_end: Some(contact(12, 62.0)),
            gamma: 1.21,
        };
        SolarEclipse {
            event_key: "2025-03-29#1264".to_string(),
            kind: SolarEclipseKind::Partial,
            global,
            bessel: minimal_bessel(),
            metadata: metadata(),
        }
    }

    // ==================================================================
    // 1. parse_date（FAST・独立オラクル = from_gregorian）
    // ==================================================================

    /// 妥当な `YYYY-MM-DD` 3 件が UTC 0:00:00 の `UtcInstant`（`from_gregorian` 独立オラクル）と
    /// round-trip 一致する。日・月・年がすべて正しく配線されることを 3 件で縛る。
    /// 殺す変異: 年↔月↔日の取り違え、時刻を 0:00 以外にする、月/日のオフ・バイ・ワン。
    #[test]
    fn parse_date_valid_round_trips_via_from_gregorian() {
        assert_eq!(
            parse_date("2024-01-15").expect("妥当日付はパース成功"),
            utc(2024, 1, 15, 0, 0, 0.0),
            "2024-01-15 は UTC 0:00 の from_gregorian と一致"
        );
        assert_eq!(
            parse_date("2017-08-21").expect("妥当日付はパース成功"),
            utc(2017, 8, 21, 0, 0, 0.0),
            "2017-08-21（皆既日）も一致"
        );
        assert_eq!(
            parse_date("2000-12-31").expect("妥当日付はパース成功"),
            utc(2000, 12, 31, 0, 0, 0.0),
            "2000-12-31（年末・桁境界）も一致"
        );
    }

    /// パース結果の時刻成分が **00:00:00** であることを `to_gregorian` で直接確認する
    /// （round-trip では時刻ずれが UtcInstant 等値の中に埋もれうるため、時/分/秒を陽に縛る）。
    /// 殺す変異: 正午（12:00, JD 慣習）など 0:00 以外の時刻に固定する。
    #[test]
    fn parse_date_sets_midnight_utc() {
        let (y, mo, d, h, mi, s) = parse_date("2024-06-19")
            .expect("妥当日付はパース成功")
            .to_gregorian();
        assert_eq!((y, mo, d), (2024, 6, 19), "暦日が一致");
        assert_eq!(h, 0, "時=0");
        assert_eq!(mi, 0, "分=0");
        assert!(s.abs() < 1e-6, "秒=0（got {s}）");
    }

    /// 月範囲外 `"2024-13-01"` は `Err(CliError::InvalidDate(入力文字列))`（入力をそのまま保持）。
    /// 殺す変異: 範囲外月を黙って受理する、InvalidDate に入力文字列でなく別文字列を載せる。
    #[test]
    fn parse_date_month_out_of_range_is_invalid_date() {
        let r = parse_date("2024-13-01");
        match r {
            Err(CliError::InvalidDate(s)) => {
                assert_eq!(s, "2024-13-01", "InvalidDate は入力文字列を保持");
            }
            other => panic!("expected Err(InvalidDate(\"2024-13-01\")), got {other:?}"),
        }
    }

    /// 非日付文字列 `"abc"` は `Err(CliError::InvalidDate("abc"))`。
    /// 殺す変異: 数値以外を 0 等にフォールバックして Ok を返す。
    #[test]
    fn parse_date_garbage_is_invalid_date() {
        let r = parse_date("abc");
        assert!(
            matches!(r, Err(CliError::InvalidDate(ref s)) if s == "abc"),
            "expected Err(InvalidDate(\"abc\")), got {r:?}"
        );
    }

    /// 空文字列 `""` は `Err(CliError::InvalidDate(""))`。
    /// 殺す変異: 空入力を既定日付（エポック等）にフォールバックして Ok を返す。
    #[test]
    fn parse_date_empty_is_invalid_date() {
        let r = parse_date("");
        assert!(
            matches!(r, Err(CliError::InvalidDate(ref s)) if s.is_empty()),
            "expected Err(InvalidDate(\"\")), got {r:?}"
        );
    }

    // ==================================================================
    // 2. kind_matches（FAST・真理値表 5×6, 非中心含む）
    // ==================================================================

    /// 全 6 種別を返すヘルパ（真理値表の走査に使う）。
    fn all_kinds() -> [SolarEclipseKind; 6] {
        [
            SolarEclipseKind::Partial,
            SolarEclipseKind::Annular,
            SolarEclipseKind::Total,
            SolarEclipseKind::Hybrid,
            SolarEclipseKind::NonCentralTotal,
            SolarEclipseKind::NonCentralAnnular,
        ]
    }

    /// `All` は任意の種別で true（6 種すべて）。
    /// 殺す変異: All を特定種別だけ true にする、常に false にする。
    #[test]
    fn kind_matches_all_accepts_every_kind() {
        for k in all_kinds() {
            assert!(
                kind_matches(KindFilter::All, k),
                "All は {k:?} を受理するべき"
            );
        }
    }

    /// `Total` は `Total` と `NonCentralTotal` のみ true、他 4 種は false（非中心皆既を含む確定仕様）。
    /// 殺す変異: NonCentralTotal を除外する、Annular/部分まで通す、常に true/false。
    #[test]
    fn kind_matches_total_includes_non_central_total_only() {
        for k in all_kinds() {
            let expected = matches!(
                k,
                SolarEclipseKind::Total | SolarEclipseKind::NonCentralTotal
            );
            assert_eq!(
                kind_matches(KindFilter::Total, k),
                expected,
                "Total フィルタ × {k:?} の真理値"
            );
        }
    }

    /// `Annular` は `Annular` と `NonCentralAnnular` のみ true、他は false（非中心金環を含む確定仕様）。
    /// 殺す変異: NonCentralAnnular を除外する、Total を通す、常に true/false。
    #[test]
    fn kind_matches_annular_includes_non_central_annular_only() {
        for k in all_kinds() {
            let expected = matches!(
                k,
                SolarEclipseKind::Annular | SolarEclipseKind::NonCentralAnnular
            );
            assert_eq!(
                kind_matches(KindFilter::Annular, k),
                expected,
                "Annular フィルタ × {k:?} の真理値"
            );
        }
    }

    /// `Partial` は `Partial` のみ true、他 5 種（非中心含む）は false。
    /// 殺す変異: 非中心や Hybrid まで Partial に含める、常に true/false。
    #[test]
    fn kind_matches_partial_matches_partial_only() {
        for k in all_kinds() {
            let expected = matches!(k, SolarEclipseKind::Partial);
            assert_eq!(
                kind_matches(KindFilter::Partial, k),
                expected,
                "Partial フィルタ × {k:?} の真理値"
            );
        }
    }

    /// `Hybrid` は `Hybrid` のみ true、他 5 種は false。
    /// 殺す変異: Total/Annular を Hybrid に含める、常に true/false。
    #[test]
    fn kind_matches_hybrid_matches_hybrid_only() {
        for k in all_kinds() {
            let expected = matches!(k, SolarEclipseKind::Hybrid);
            assert_eq!(
                kind_matches(KindFilter::Hybrid, k),
                expected,
                "Hybrid フィルタ × {k:?} の真理値"
            );
        }
    }

    // ==================================================================
    // 3. format_search_text（FAST・既知 fixture の内容存在で縛る）
    // ==================================================================

    /// 皆既 fixture を整形すると、event_key・種別名・最大食 UTC と TT の両方・gamma・食分・
    /// 食面積・計算メタデータ（ephemeris/ΔT モデル名・ΔT 不確実性帯）が出力に含まれる。
    /// **厳密レイアウトは縛らず内容（部分文字列）の存在で縛る**（accuracy.md §0: UTC+TT 併記必須）。
    /// 殺す変異: event_key/種別/gamma/食分/食面積/メタデータの欠落、TT を出さず UTC のみ出す。
    #[test]
    fn format_search_text_contains_all_key_fields_for_total() {
        let out = format_search_text(&[total_eclipse()]);

        // event_key（安定キー）。
        assert!(
            out.contains("2024-04-08#1252"),
            "event_key が出力に含まれる: {out}"
        );
        // 種別名（皆既を示す表記。Debug 名 "Total" を最小契約とする）。
        assert!(out.contains("Total"), "種別 Total が出力に含まれる: {out}");
        // 最大食 UTC の日付（暦日）。
        assert!(
            out.contains("2024-04-08"),
            "最大食 UTC 日付が出力に含まれる: {out}"
        );
        // 最大食 TT も併記される（UTC+TT 必須, accuracy.md §0）。TT は "TT" ラベルで縛る。
        assert!(
            out.contains("TT"),
            "最大食 TT（TT ラベル併記）が出力に含まれる: {out}"
        );
        // gamma 値（0.3431 の有効数字。先頭桁列を部分一致で）。
        assert!(
            out.contains("0.3431"),
            "gamma=0.3431 が出力に含まれる: {out}"
        );
        // 食分 magnitude（1.0566）。
        assert!(
            out.contains("1.0566"),
            "食分 1.0566 が出力に含まれる: {out}"
        );
        // 食面積 obscuration（1.0 を含む浮動小数表記。"1.0" 部分一致）。
        assert!(out.contains("1.0"), "食面積（1.0）が出力に含まれる: {out}");
        // 計算メタデータ: ephemeris モデル名・ΔT モデル名。
        assert!(
            out.contains("ELP/MPP02+VSOP87D"),
            "ephemeris モデル名が出力に含まれる: {out}"
        );
        assert!(
            out.contains("EspenakMeeus"),
            "ΔT モデル名が出力に含まれる: {out}"
        );
        // ΔT 不確実性帯（0.5 秒）。
        assert!(
            out.contains("0.5"),
            "ΔT 不確実性帯（0.5）が出力に含まれる: {out}"
        );
    }

    /// 複数 fixture（皆既＋部分）を整形すると、両方の event_key・両方の種別が出力に含まれる
    /// （リスト全要素が出力される＝先頭だけ・末尾だけ出す変異を撃破）。
    /// 殺す変異: 入力リストの一部のみ整形する、2 件目を落とす。
    #[test]
    fn format_search_text_includes_every_eclipse_in_list() {
        let out = format_search_text(&[total_eclipse(), partial_eclipse()]);
        assert!(
            out.contains("2024-04-08#1252"),
            "1 件目（皆既）の event_key: {out}"
        );
        assert!(
            out.contains("2025-03-29#1264"),
            "2 件目（部分）の event_key: {out}"
        );
        assert!(out.contains("Total"), "1 件目の種別 Total: {out}");
        assert!(out.contains("Partial"), "2 件目の種別 Partial: {out}");
    }

    /// 空リストは panic せず「該当なし」を表す出力を返す（空文字 or "no eclipses" 等）。
    /// fixture の event_key/種別名のような具体的日食情報を含まないことだけ確認する。
    /// 殺す変異: 空リストで panic する、空でないのに架空の日食を 1 件捏造して出す。
    #[test]
    fn format_search_text_empty_list_is_non_panicking() {
        let out = format_search_text(&[]);
        // 非 panic（ここに到達すれば panic していない）。具体的日食情報を含まない。
        assert!(
            !out.contains("2024-04-08#1252"),
            "空リスト出力に具体的日食 event_key は含まれない: {out}"
        );
    }

    // ==================================================================
    // 4. run_search（fast-fail 2 件 = FAST / 正常系 1 件 = SLOW）
    // ==================================================================

    /// `from > to`（妥当日付だが逆順）は search を呼ぶ前に `Err(CliError::RangeOrder{..})`
    /// （fast-fail）。RangeOrder には入力された from/to 文字列が載る。
    /// このテストは search を実走しないため高速（red 段階では入口の unimplemented! で panic）。
    /// 殺す変異: 逆順を検出せず search に渡す、RangeOrder でなく別 variant を返す、from/to を入替。
    #[test]
    fn run_search_from_after_to_is_range_order_error() {
        let args = SearchArgs {
            from: "2020-12-31".to_string(),
            to: "2020-01-01".to_string(),
            accuracy: AccuracyArg::Standard,
            kind: KindFilter::All,
        };
        let r = run_search(&args);
        match r {
            Err(CliError::RangeOrder { from, to }) => {
                assert_eq!(from, "2020-12-31", "RangeOrder.from に入力 from");
                assert_eq!(to, "2020-01-01", "RangeOrder.to に入力 to");
            }
            other => panic!("expected Err(RangeOrder {{..}}), got {other:?}"),
        }
    }

    /// `from` が不正日付なら search を呼ぶ前に `Err(CliError::InvalidDate)`（fast-fail）。
    /// to は妥当でも from の不正で即エラー（search 非実走＝高速）。
    /// 殺す変異: 不正 from を黙って受理して search に進む、InvalidDate でなく別 variant を返す。
    #[test]
    fn run_search_invalid_from_is_invalid_date_error() {
        let args = SearchArgs {
            from: "not-a-date".to_string(),
            to: "2020-12-31".to_string(),
            accuracy: AccuracyArg::Standard,
            kind: KindFilter::All,
        };
        let r = run_search(&args);
        assert!(
            matches!(r, Err(CliError::InvalidDate(ref s)) if s == "not-a-date"),
            "expected Err(InvalidDate(\"not-a-date\")), got {r:?}"
        );
    }

    /// 【SLOW・1 件】正常系: 2017-08 の範囲・Standard・All で `Ok(output)`。output に 2017 皆既の
    /// event_key 日付 `"2017-08-21"` と皆既を示す種別表記 `"Total"` が含まれる（物理事実オラクル）。
    /// 内部で `EclipseEngine::search`（実日食解, ≈120s）を実走するため SLOW。red 段階では入口の
    /// `unimplemented!` で panic（search 非実走）するので速い。
    /// 殺す変異: 範囲を search に渡さない、kind フィルタが All で日食を落とす、整形を呼ばない。
    #[test]
    fn run_search_2017_total_eclipse_appears_in_output() {
        let args = SearchArgs {
            from: "2017-08-01".to_string(),
            to: "2017-09-01".to_string(),
            accuracy: AccuracyArg::Standard,
            kind: KindFilter::All,
        };
        let out = run_search(&args).expect("2017-08 の探索は成功する");
        assert!(
            out.contains("2017-08-21"),
            "2017 皆既の event_key 日付が出力に含まれる: {out}"
        );
        assert!(
            out.contains("Total"),
            "皆既を示す種別表記 Total が出力に含まれる: {out}"
        );
    }
}
