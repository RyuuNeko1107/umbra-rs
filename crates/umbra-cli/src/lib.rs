//! `umbra` CLI ライブラリ（ISSUE-031 `umbra search` / ISSUE-032 `umbra local` / ISSUE-048 `umbra path`）。
//!
//! 薄い CLI ラッパ: 引数解釈（clap）・日付パース・`EclipseEngine` 呼び出し・整形出力。
//! 計算は umbra-eclipse が担保。本クレートは境界（引数・パース・出力・エラー/終了コード）が責務。
//!
//! - `umbra search`（`--format <text|json>`）: S31a text、S31b json（serde 横断配線・
//!   `SolarEclipse` 推移閉包に Serialize を通し `serde_json` で整形）。
//! - `umbra local`（S32a・text）: 指定日・指定地点の局地条件（`EclipseEngine::local_circumstances`）。
//!   西経入力吸収（`Observer::from_degrees`）・UTC オフセット表示・可視性 6 値。`--format json` は S32b。
//! - `umbra path`（`--format <text|geojson>`・ISSUE-048）: 指定日の日食の経路（`EclipseEngine::path`）。
//!   text は中心線・南北限界線・部分食域の要約、geojson は `EclipsePath::to_geojson` の生出力。

use clap::{Args, Parser, Subcommand, ValueEnum};
use umbra_core::{
    jd2_to_gregorian_deciseconds, DomainError, EspenakMeeusDeltaT, JulianDate2, Observer,
    UtcInstant,
};
use umbra_eclipse::{
    standard_engine, EclipseEngine, EclipseError, EclipsePath, EngineConfig, LocalCircumstances,
    LocalContact, PathOptions, RefractionModel, SolarEclipse, SolarEclipseKind, UtcRange,
    VisibleSolarEclipse,
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

/// サブコマンド（`search` = ISSUE-031 ／ `local` = ISSUE-032 ／ `path` = ISSUE-048。
/// `bessel`/`inspect`/`validate`/`bench` は後続 issue）。
#[derive(Debug, Subcommand)]
pub enum Command {
    /// 期間内の太陽食を列挙する（`EclipseEngine::search`）。
    Search(SearchArgs),
    /// 指定日・指定地点の局地条件を表示する（`EclipseEngine::local_circumstances`）。
    Local(LocalArgs),
    /// 指定日の日食の経路（中心線・限界線・部分食域）を表示する（`EclipseEngine::path`）。
    Path(PathArgs),
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
    /// 出力形式（既定 text）。
    #[arg(long, value_enum, default_value_t = FormatArg::Text)]
    pub format: FormatArg,
}

/// `umbra path` の引数（ISSUE-048）。
#[derive(Debug, Args)]
pub struct PathArgs {
    /// 対象日（`YYYY-MM-DD`, UTC）。当日に起こる日食の経路を求める。
    #[arg(long)]
    pub date: String,
    /// 精度プロファイル（既定 standard・`search`/`local` と同一）。
    #[arg(long, value_enum, default_value_t = AccuracyArg::Standard)]
    pub accuracy: AccuracyArg,
    /// 出力形式（既定 text）。
    #[arg(long, value_enum, default_value_t = PathFormatArg::Text)]
    pub format: PathFormatArg,
    /// サンプル間隔 \[s\]（既定は `PathOptions::default()`）。正の有限値のみ。
    #[arg(long, default_value_t = PathOptions::default().sample_interval_seconds)]
    pub interval: f64,
    /// 指定時は限界線・部分食域を計算しない（`include_limits=false`）。
    #[arg(long)]
    pub no_limits: bool,
}

/// `umbra path` の出力形式（ISSUE-048）。`search`/`local` の [`FormatArg`] とは別型
/// （経路の機械可読形は汎用 JSON ではなく **GeoJSON** ＝ `EclipsePath::to_geojson` だから）。
#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
pub enum PathFormatArg {
    /// 人間可読テキスト（既定）。
    Text,
    /// GeoJSON（`EclipsePath::to_geojson` の生出力）。
    Geojson,
}

/// 出力形式引数（S31b）。
#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
pub enum FormatArg {
    /// 人間可読テキスト（既定）。
    Text,
    /// JSON（`SolarEclipse` の serde・配列。機械可読・列挙は `{type:..}` タグ付き）。
    Json,
}

/// `umbra local` の引数（S32a。`--format` は S32b で追加・本スライスは text のみ）。
#[derive(Debug, Args)]
pub struct LocalArgs {
    /// 対象日（`YYYY-MM-DD`, UTC）。当日に起こる日食の局地条件を求める。
    #[arg(long)]
    pub date: String,
    /// 測地緯度（度, [-90, 90]）。負＝南緯（`allow_hyphen_values` で受理）。
    #[arg(long, allow_hyphen_values = true)]
    pub lat: f64,
    /// 経度（度・東経正。**負＝西経**も受理し東経へ正規化吸収, conventions §3）。
    #[arg(long, allow_hyphen_values = true)]
    pub lon: f64,
    /// 楕円体高（m, 既定 0）。
    #[arg(long, default_value_t = 0.0, allow_hyphen_values = true)]
    pub elevation: f64,
    /// ローカル時刻表示用の UTC オフセット（例 `+09:00` / `-0500` / `Z`）。内部計算は UTC/TT 不変。
    #[arg(long)]
    pub timezone: Option<String>,
    /// 精度プロファイル（既定 standard）。
    #[arg(long, value_enum, default_value_t = AccuracyArg::Standard)]
    pub accuracy: AccuracyArg,
    /// 大気差モデル（既定 standard・conventions §7 / EngineConfig 既定と一致）。
    #[arg(long, value_enum, default_value_t = RefractionArg::Standard)]
    pub refraction: RefractionArg,
    /// 出力形式（既定 text。S32b で json 追加）。
    #[arg(long, value_enum, default_value_t = FormatArg::Text)]
    pub format: FormatArg,
}

/// 大気差モデル引数（S32a）。`umbra_eclipse::RefractionModel` に対応。
#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
pub enum RefractionArg {
    /// 大気差なし（幾何学的高度のみ）。
    None,
    /// 標準大気差（既定・Saemundsson）。
    Standard,
}

impl From<RefractionArg> for RefractionModel {
    fn from(arg: RefractionArg) -> Self {
        match arg {
            RefractionArg::None => RefractionModel::None,
            RefractionArg::Standard => RefractionModel::Standard,
        }
    }
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
    /// JSON 整形失敗（`--format json`・serde_json 由来, 透過）。
    #[error(transparent)]
    Json(#[from] serde_json::Error),
    /// `--timezone` パース失敗（UTC オフセット形式以外）。
    #[error("invalid timezone '{0}' (expected UTC offset like +09:00, -0500, or Z)")]
    InvalidTimezone(String),
    /// 入力の定義域違反（緯度/経度範囲外など, 透過）。
    #[error(transparent)]
    Domain(#[from] DomainError),
    /// `--interval` が正の有限値でない（`umbra path`・ISSUE-048 §4）。
    #[error("invalid --interval {0} (expected a positive, finite number of seconds)")]
    InvalidInterval(f64),
    /// `--interval` が小さすぎて経路のサンプル数が上限を超える（`umbra path`・ISSUE-048 §4）。
    #[error(
        "--interval {interval} is too small: it would need about {estimated_samples:.0} path samples (limit {MAX_PATH_SAMPLES})"
    )]
    IntervalTooSmall {
        /// 指定された間隔 \[s\]。
        interval: f64,
        /// 推定サンプル数（P1〜P4 の秒数 ÷ 間隔）。
        estimated_samples: f64,
    },
}

/// `umbra path` が許す経路サンプル数の上限。
///
/// **正本は [`umbra_eclipse::MAX_PATH_SAMPLES`]**（ISSUE-049 でエンジン側に入力検証を置いた）。
/// CLI 側は数値を重複定義せず再輸出するだけにして、上限が二重管理で食い違わないようにする。
/// CLI の検査は「`path` を呼ぶ前に、より親切なメッセージで早期に弾く」ためのもので、
/// **安全性の正本はエンジン側**（CLI を経由しない呼び出し元も守られる）。
pub use umbra_eclipse::MAX_PATH_SAMPLES;

/// `--interval` が経路サンプル数の上限に収まるかを検査する（ISSUE-048 §4）。
///
/// `span_seconds` は P1〜P4 の秒数。`None`（P1/P4 が欠ける・起こり得ない想定）のときは
/// **検査しない**（捏造した span で判定しない）。`interval` は呼び出し前に正の有限値であることが
/// 保証されている前提。
pub(crate) fn check_path_sample_count(
    span_seconds: Option<f64>,
    interval: f64,
) -> Result<(), CliError> {
    let Some(span) = span_seconds else {
        return Ok(());
    };
    let estimated_samples = span / interval;
    if estimated_samples > MAX_PATH_SAMPLES {
        return Err(CliError::IntervalTooSmall {
            interval,
            estimated_samples,
        });
    }
    Ok(())
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
        // 表示は 0.1 秒丸め（暦往復 ±eps による "15:59:60.0" 不正表記を回避・共通ヘルパ）。
        let (y, mo, d, h, mi, s) = jd2_to_gregorian_deciseconds(g.time_utc.jd2());
        // TT も暦形式で併記（accuracy.md §0: UTC+TT 両方）。TT-JD をグレゴリオ暦に直す。
        let (ty, tmo, td, th, tmi, ts) = jd2_to_gregorian_deciseconds(g.time_tt.jd2());
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

/// 日食リストを JSON（`SolarEclipse` の serde 配列）に整形する（S31b）。
///
/// 出力は pretty-print された JSON **配列**（1 日食 1 要素・入力順）＋末尾改行。空リストは
/// `[]`（該当なしは空配列・エラーにしない, api-draft §3.2）。各日食は自身の `metadata`
/// （暦/ΔT モデル・不確実性帯）を含む。時刻は `{iso, jd}`（自己記述＋可逆）、列挙は
/// `{type:..}`（A7 タグ付き）、数値の単位はフィールド名（`path_width_km` 等, A7）で出力する。
/// 改変・丸めで精度を捏造しない（コア値を素通し, accuracy.md §0）。
pub fn format_search_json(eclipses: &[SolarEclipse]) -> Result<String, CliError> {
    let mut out = serde_json::to_string_pretty(eclipses)?;
    out.push('\n');
    Ok(out)
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
    match args.format {
        FormatArg::Text => Ok(format_search_text(&filtered)),
        FormatArg::Json => format_search_json(&filtered),
    }
}

/// `umbra path` を実行し、整形済み出力を返す（ISSUE-048）。日付パース→`--interval` 検証→
/// エンジン構築→当日の `search`→`path`→整形。出力は呼び出し側（main）が印字する。
///
/// 不正日付・正の有限でない `--interval` は **エンジンを呼ぶ前に** fast-fail
/// （[`CliError::InvalidDate`] / [`CliError::InvalidInterval`]）。当日に日食が無いのは**成功**で、
/// text は 1 行の告知・geojson は**空 `FeatureCollection`**（`null` は GeoJSON として妥当でない）。
pub fn run_path(args: &PathArgs) -> Result<String, CliError> {
    let date = parse_date(&args.date)?;
    // 正の有限値のみ受理する。`is_finite()` が NaN・±∞ を弾き、`> 0.0` が 0・負を弾く
    // （NaN は `is_finite()` 側で落ちるので、順序比較を否定する必要はない）。
    if !args.interval.is_finite() || args.interval <= 0.0 {
        return Err(CliError::InvalidInterval(args.interval));
    }

    let time = bundled_time_data();
    let config = match args.accuracy {
        AccuracyArg::Standard => EngineConfig::standard(),
        AccuracyArg::Reference => EngineConfig::reference(),
    };
    let earth_orientation = time.eop().clone();
    let engine = EclipseEngine::new(
        AnalyticalEphemeris::new(),
        EspenakMeeusDeltaT,
        earth_orientation,
        time,
        config,
    );

    // 指定日（UTC 暦日 [date, date+1 日)）に起こる日食を探索する（`run_local` と同型）。
    let day_end = UtcInstant::from_jd2(JulianDate2::from_jd(date.jd2().jd() + 1.0));
    let range = UtcRange {
        start: date,
        end: day_end,
    };
    let Some(eclipse) = engine.search(range)?.into_iter().next() else {
        // 該当日食なし（エラーにしない・架空の経路を作らない）。
        return Ok(match args.format {
            PathFormatArg::Text => format!(
                "No solar eclipse on {}.
",
                args.date
            ),
            PathFormatArg::Geojson => {
                let mut out = serde_json::to_string_pretty(&serde_json::json!({
                    "type": "FeatureCollection",
                    "features": [],
                }))?;
                out.push('\n');
                out
            }
        });
    };

    // P1〜P4 の秒数でサンプル数を見積もり、上限超えなら `path` を呼ぶ前に弾く（ISSUE-048 §4）。
    let span_seconds = match (&eclipse.global.partial_begin, &eclipse.global.partial_end) {
        (Some(p1), Some(p4)) => Some((p4.time_tt.jd2().jd() - p1.time_tt.jd2().jd()) * 86_400.0),
        _ => None,
    };
    check_path_sample_count(span_seconds, args.interval)?;

    let options = PathOptions {
        sample_interval_seconds: args.interval,
        include_limits: !args.no_limits,
        ..PathOptions::default()
    };
    let path = engine.path(&eclipse, options)?;

    match args.format {
        PathFormatArg::Text => Ok(format_path_text(&eclipse, &path)),
        PathFormatArg::Geojson => {
            let mut out = path.to_geojson()?;
            out.push('\n');
            Ok(out)
        }
    }
}

/// 経路を人間可読テキストへ整形する（ISSUE-048 §5・ラベルは仕様の表で固定）。
///
/// **`None` の要素は行を省略せず「なし」と明示**する（conventions §11「未提供を隠さない」）。
/// 時刻は 0.1 秒丸めの共通ヘルパ（暦往復 ±eps による不正表記の回避）。
pub fn format_path_text(eclipse: &SolarEclipse, path: &EclipsePath) -> String {
    let g = &eclipse.global.greatest;
    let (y, mo, d, h, mi, sec) = jd2_to_gregorian_deciseconds(g.time_utc.jd2());
    let mut out = format!(
        "{key}
",
        key = eclipse.event_key
    );
    out.push_str(&format!(
        "  種別: {kind:?}
",
        kind = eclipse.kind
    ));
    out.push_str(&format!(
        "  最大食: {y:04}-{mo:02}-{d:02} {h:02}:{mi:02}:{sec:04.1} UTC           lat {lat:.4}° lon {lon:.4}°
",
        lat = path.greatest_point.lat.degrees().0,
        lon = path.greatest_point.lon.degrees().0,
    ));

    out.push_str(&match &path.center_line {
        Some(line) if !line.points.is_empty() => {
            let first = line.points[0];
            let last = line.points[line.points.len() - 1];
            format!(
                "  中心線: {n} 点  始 lat {flat:.4}° lon {flon:.4}°  終 lat {llat:.4}° lon {llon:.4}°
",
                n = line.points.len(),
                flat = first.lat.degrees().0,
                flon = first.lon.degrees().0,
                llat = last.lat.degrees().0,
                llon = last.lon.degrees().0,
            )
        }
        // 点列が空なら「0 点」ではなく「なし」（空の成功出力を作らない）。
        _ => "  中心線: なし
".to_string(),
    });
    out.push_str(&format_limit_line(
        "北限",
        path.northern_limit.as_ref().map(|l| l.points.len()),
    ));
    out.push_str(&format_limit_line(
        "南限",
        path.southern_limit.as_ref().map(|l| l.points.len()),
    ));
    out.push_str(&match &path.partial_limit {
        Some(poly) if !poly.rings.is_empty() => format!(
            "  部分食域: {rings} 環  外環 {outer} 頂点
",
            rings = poly.rings.len(),
            outer = poly.rings[0].len(),
        ),
        _ => "  部分食域: なし
"
        .to_string(),
    });
    out.push_str(&format!(
        "  samples: {n}
",
        n = path.samples.len()
    ));
    out
}

/// 限界線 1 本を整形する（ISSUE-048 §5）。`None`・空点列はいずれも「なし」
/// （行を省略しない＝未提供を隠さない・conventions §11）。
fn format_limit_line(label: &str, points: Option<usize>) -> String {
    match points {
        Some(n) if n > 0 => format!(
            "  {label}: {n} 点
"
        ),
        _ => format!(
            "  {label}: なし
"
        ),
    }
}

/// UTC オフセット文字列を **符号付き分**へパースする（`umbra local --timezone`・S32a）。
///
/// 受理: `"Z"`（=0）, `"+HH:MM"` / `"-HH:MM"` / `"+HHMM"` / `"-HHMM"`（HH∈[00,23], MM∈[00,59]）。
/// それ以外は [`CliError::InvalidTimezone`]（入力文字列を保持）。表示専用で内部 UTC/TT は不変。
pub fn parse_utc_offset(text: &str) -> Result<i32, CliError> {
    let invalid = || CliError::InvalidTimezone(text.to_string());
    if text == "Z" {
        return Ok(0);
    }
    let (sign, rest) = if let Some(r) = text.strip_prefix('+') {
        (1, r)
    } else if let Some(r) = text.strip_prefix('-') {
        (-1, r)
    } else {
        return Err(invalid());
    };
    // "HH:MM"（コロンあり）または "HHMM"（コロンなし 4 桁）。
    let (hh, mm) = if let Some((h, m)) = rest.split_once(':') {
        (h, m)
    } else if rest.len() == 4 {
        rest.split_at(2)
    } else {
        return Err(invalid());
    };
    if hh.len() != 2 || mm.len() != 2 {
        return Err(invalid());
    }
    let hours: i32 = hh.parse().map_err(|_| invalid())?;
    let minutes: i32 = mm.parse().map_err(|_| invalid())?;
    if !(0..=23).contains(&hours) || !(0..=59).contains(&minutes) {
        return Err(invalid());
    }
    Ok(sign * (hours * 60 + minutes))
}

/// 局地接触 1 点を整形する（`time_utc` UTC ＋ `time_tt` TT ＋任意でローカル時刻・高度方位）。
/// `None`（部分食地点の C2/C3 など）は em ダッシュ `—` で描き、架空時刻を捏造しない。
fn format_contact(
    label: &str,
    contact: Option<&LocalContact>,
    timezone: Option<(&str, i32)>,
) -> String {
    let Some(c) = contact else {
        return format!("  {label}: —\n");
    };
    // 表示は 0.1 秒丸め（暦往復 ±eps による "15:59:60.0" 不正表記を回避・共通ヘルパ）。
    let (y, mo, d, h, mi, s) = jd2_to_gregorian_deciseconds(c.time_utc.jd2());
    let (ty, tmo, td, th, tmi, ts) = jd2_to_gregorian_deciseconds(c.time_tt.jd2());
    let mut line = format!(
        "  {label}: {y:04}-{mo:02}-{d:02} {h:02}:{mi:02}:{s:04.1} UTC / \
         {ty:04}-{tmo:02}-{td:02} {th:02}:{tmi:02}:{ts:04.1} TT",
    );
    if let Some((tz_label, offset_min)) = timezone {
        // ローカル時刻は UTC を表示のためだけにオフセット（内部 UTC/TT は不変, conventions §6）。
        let local_jd = c.time_utc.jd2().jd() + f64::from(offset_min) / 1440.0;
        let (ly, lmo, ld, lh, lmi, ls) =
            jd2_to_gregorian_deciseconds(JulianDate2::from_jd(local_jd));
        line.push_str(&format!(
            " / {ly:04}-{lmo:02}-{ld:02} {lh:02}:{lmi:02}:{ls:04.1} {tz_label}"
        ));
    }
    line.push_str(&format!(
        "  (alt {:.1}°, az {:.1}°)\n",
        c.sun_altitude.0, c.sun_azimuth.0
    ));
    line
}

/// 局地条件を人間可読 text に整形する（S32a）。接触 C1/C2/最大/C3/C4 を時系列で（各 **UTC+TT**・
/// `timezone` 指定時はローカル時刻も併記, accuracy.md §0）、食分・食面積・最大高度・**可視性 6 値**・
/// 計算メタデータ（暦/ΔT モデル名・ΔT 不確実性帯）を出す。部分食地点の C2/C3（None）は `—`。
pub fn format_local_text(circ: &LocalCircumstances, timezone: Option<(&str, i32)>) -> String {
    let cs = &circ.contacts;
    let mut out = String::new();
    out.push_str(&format_contact("C1 ", cs.c1.as_ref(), timezone));
    out.push_str(&format_contact("C2 ", cs.c2.as_ref(), timezone));
    out.push_str(&format_contact("max", Some(&cs.maximum), timezone));
    out.push_str(&format_contact("C3 ", cs.c3.as_ref(), timezone));
    out.push_str(&format_contact("C4 ", cs.c4.as_ref(), timezone));
    out.push_str(&format!(
        "  magnitude: {mag:.4}  obscuration: {obsc:.4}  max altitude: {alt:.1}°\n",
        mag = circ.magnitude.0,
        obsc = circ.obscuration.0,
        alt = circ.maximum_altitude.0,
    ));
    out.push_str(&format!("  visibility: {:?}\n", circ.visibility));
    let m = &circ.metadata;
    out.push_str(&format!(
        "  ephemeris: {em} {ev}  ΔT: {dt} (±{unc:.2}s)  accuracy: {acc:?}\n",
        em = m.ephemeris_model,
        ev = m.ephemeris_version,
        dt = m.delta_t_model,
        unc = m.delta_t_uncertainty_seconds,
        acc = m.accuracy_profile,
    ));
    out
}

/// `umbra local` を実行し、整形済み text 出力を返す（S32a）。
///
/// 不正日付・緯度経度範囲外・不正 timezone は **エンジン実走前に** fast-fail
/// （[`CliError::InvalidDate`]/[`CliError::Domain`]/[`CliError::InvalidTimezone`]）。
/// `--date` の UTC 暦日 `[date, date+1日)` を `search` し、見つかった日食に
/// `local_circumstances(observer)` を適用する。該当日食なしは「食なし」を返す（エラーにしない）。
/// 西経入力は [`Observer::from_degrees`] が東経へ正規化吸収する（conventions §3）。
pub fn run_local(args: &LocalArgs) -> Result<String, CliError> {
    let date = parse_date(&args.date)?;
    let observer = Observer::from_degrees(args.lat, args.lon, args.elevation)?;
    let timezone: Option<(&str, i32)> = match &args.timezone {
        Some(tz) => Some((tz.as_str(), parse_utc_offset(tz)?)),
        None => None,
    };

    // 大気差を反映したエンジン設定（精度プロファイル＋ refraction 上書き）。
    let mut config = match args.accuracy {
        AccuracyArg::Standard => EngineConfig::standard(),
        AccuracyArg::Reference => EngineConfig::reference(),
    };
    config.refraction = args.refraction.into();

    let time = bundled_time_data();
    let earth_orientation = time.eop().clone();
    let engine = EclipseEngine::new(
        AnalyticalEphemeris::new(),
        EspenakMeeusDeltaT,
        earth_orientation,
        time,
        config,
    );

    // 指定日（UTC 暦日 [date, date+1 日)）に起こる日食を探索し、あれば局地条件を求める。
    let day_end = UtcInstant::from_jd2(JulianDate2::from_jd(date.jd2().jd() + 1.0));
    let range = UtcRange {
        start: date,
        end: day_end,
    };
    let found: Option<VisibleSolarEclipse> = match engine.search(range)?.into_iter().next() {
        Some(eclipse) => {
            let local = engine.local_circumstances(&eclipse, observer)?;
            Some(VisibleSolarEclipse { eclipse, local })
        }
        None => None,
    };

    match args.format {
        FormatArg::Text => Ok(match &found {
            Some(vse) => {
                let mut out = format!(
                    "{key}  {kind:?}\n",
                    key = vse.eclipse.event_key,
                    kind = vse.eclipse.kind
                );
                out.push_str(&format_local_text(&vse.local, timezone));
                out
            }
            // 該当日食なし（エラーにしない・架空 event_key を出さない）。
            None => format!("No solar eclipse on {} at this location.\n", args.date),
        }),
        FormatArg::Json => format_local_json(found.as_ref()),
    }
}

/// 局地条件を JSON に整形する（S32b）。該当日食あり→`VisibleSolarEclipse`（`{eclipse, local}`）の
/// pretty JSON、なし→JSON `null`。いずれも末尾改行付き。時刻は `{iso, jd}`、列挙は `{type:..}`、
/// 数値の単位はフィールド名（`maximum_altitude_deg` 等, A7）。コア値を素通し（accuracy.md §0）。
pub fn format_local_json(found: Option<&VisibleSolarEclipse>) -> Result<String, CliError> {
    let mut out = serde_json::to_string_pretty(&found)?;
    out.push('\n');
    Ok(out)
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
            format: FormatArg::Text,
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
            format: FormatArg::Text,
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
            format: FormatArg::Text,
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

    // ==================================================================
    // === S31b: --format json ===
    // ==================================================================
    // S31b 受け入れテスト（standard・`umbra search --format json`）。
    //
    // ## オラクル戦略
    // 出力文字列を `serde_json::from_str::<serde_json::Value>` でパースし、**パースされた
    // Value** に対して構造を assert する（生 JSON の文字列一致は構造検証に使わない。
    // 部分文字列は二次シグナルとしてのみ許容）。値は既知 fixture（total_eclipse/partial_eclipse）
    // 由来で、観測可能な JSON 契約（凍結仕様）をそのまま縛る。
    //
    // ## red 設計（本体未実装）
    // `FormatArg`・`format_search_json`・`SearchArgs.format`・serde 配線は本スライスで導入予定で
    // 現状未実装。テストはコンパイル時点で未解決シンボル（red）。

    use serde_json::Value;

    /// `format_search_json` 出力をパースして JSON 配列 Value を返すヘルパ（パース成功＝有効 JSON）。
    fn parse_json_array(eclipses: &[SolarEclipse]) -> Vec<Value> {
        let s = format_search_json(eclipses).expect("JSON 整形は成功する");
        let v: Value = serde_json::from_str(&s).expect("出力は有効な JSON（パース成功必須）");
        v.as_array().expect("トップレベルは JSON 配列").clone()
    }

    /// 皆既 fixture 1 件の JSON 契約を全面的に縛る（event_key・kind タグ・gamma・greatest の
    /// time_utc(Z)/time_tt(非 Z)/jd・position・magnitude/obscuration（透過数値）・
    /// path_width_km/central_duration_seconds/sun_altitude_deg（改名値）・metadata・bessel）。
    /// 殺す変異: kind をタグ無し（bare string / untagged）にする、gamma/greatest 各値の欠落・誤配線、
    ///   magnitude/obscuration をオブジェクト（newtype 透過しない）にする、フィールド改名漏れ、
    ///   TT iso に Z を付ける、jd チャネルを落とす、metadata/bessel フィールド欠落。
    #[test]
    fn format_search_json_total_fixture_contract() {
        let arr = parse_json_array(&[total_eclipse()]);
        assert_eq!(arr.len(), 1, "皆既 1 件 → 配列長 1");
        let obj = &arr[0];

        // event_key（安定キー・文字列）。
        assert_eq!(
            obj["event_key"],
            Value::from("2024-04-08#1252"),
            "event_key（皆既）"
        );

        // kind は内部タグ付き enum オブジェクト {type:"Total"}（bare string でない）。
        assert_eq!(
            obj["kind"]["type"],
            Value::from("Total"),
            "kind は {{type:\"Total\"}} のタグ付きオブジェクト"
        );
        assert!(
            !obj["kind"].is_string(),
            "kind は文字列でなくオブジェクト（bare string 回帰を撃破）: {}",
            obj["kind"]
        );

        // global.
        let global = &obj["global"];
        assert_eq!(
            global["kind"]["type"],
            Value::from("Total"),
            "global.kind タグ（皆既）"
        );
        assert_eq!(
            global["gamma"].as_f64().expect("gamma は数値"),
            0.3431,
            "global.gamma（皆既）"
        );

        // 接触点 Option: 皆既は全て Some（オブジェクト）。
        assert!(
            global["partial_begin"].is_object(),
            "partial_begin は Some（オブジェクト）: {}",
            global["partial_begin"]
        );
        assert!(
            global["central_begin"].is_object(),
            "central_begin は Some（皆既）: {}",
            global["central_begin"]
        );
        assert!(
            global["central_end"].is_object(),
            "central_end は Some（皆既）: {}",
            global["central_end"]
        );
        assert!(
            global["partial_end"].is_object(),
            "partial_end は Some（皆既）: {}",
            global["partial_end"]
        );

        // P1 接触: contact(16, 10.0) → lat 10.0, time_utc iso は 2024-04-08T16:00:00 で始まる。
        assert_eq!(
            global["partial_begin"]["position"]["lat_deg"]
                .as_f64()
                .expect("P1 lat_deg は数値"),
            10.0,
            "partial_begin.position.lat_deg == 10.0"
        );
        let p1_iso = global["partial_begin"]["time_utc"]["iso"]
            .as_str()
            .expect("P1 time_utc.iso は文字列");
        assert!(
            p1_iso.starts_with("2024-04-08T16:00:00"),
            "partial_begin.time_utc.iso は 2024-04-08T16:00:00 で始まる: {p1_iso}"
        );

        // greatest.
        let greatest = &global["greatest"];

        // greatest.time_utc: {iso, jd:{part1,part2}}。iso は UTC → 末尾 Z。
        let utc_iso = greatest["time_utc"]["iso"]
            .as_str()
            .expect("greatest.time_utc.iso は文字列");
        assert_eq!(
            utc_iso, "2024-04-08T18:17:00.0Z",
            "greatest.time_utc.iso（UTC・末尾 Z）"
        );
        assert!(
            greatest["time_utc"]["jd"]["part1"].is_number(),
            "time_utc.jd.part1 は数値（lossless チャネル）: {}",
            greatest["time_utc"]["jd"]["part1"]
        );
        assert!(
            greatest["time_utc"]["jd"]["part2"].is_number(),
            "time_utc.jd.part2 は数値: {}",
            greatest["time_utc"]["jd"]["part2"]
        );

        // greatest.time_tt: {iso, jd}。TT iso は UTC でないため末尾 Z を持たない。
        let tt_iso = greatest["time_tt"]["iso"]
            .as_str()
            .expect("greatest.time_tt.iso は文字列");
        assert!(
            !tt_iso.ends_with('Z'),
            "TT iso は末尾 Z を持たない（TT は UTC でない）: {tt_iso}"
        );
        assert!(
            greatest["time_tt"]["jd"]["part1"].is_number(),
            "time_tt.jd.part1 は数値: {}",
            greatest["time_tt"]["jd"]["part1"]
        );
        assert!(
            greatest["time_tt"]["jd"]["part2"].is_number(),
            "time_tt.jd.part2 は数値: {}",
            greatest["time_tt"]["jd"]["part2"]
        );

        // greatest.position == {lat_deg:25.0, lon_deg:-104.0}。
        assert_eq!(
            greatest["position"]["lat_deg"]
                .as_f64()
                .expect("position.lat_deg は数値"),
            25.0,
            "greatest.position.lat_deg"
        );
        assert_eq!(
            greatest["position"]["lon_deg"]
                .as_f64()
                .expect("position.lon_deg は数値"),
            -104.0,
            "greatest.position.lon_deg"
        );

        // magnitude/obscuration は透過 newtype = bare number（オブジェクトでない）。
        assert_eq!(
            greatest["magnitude"].as_f64().expect("magnitude は数値"),
            1.0566,
            "greatest.magnitude（透過数値）"
        );
        assert!(
            greatest["magnitude"].is_number(),
            "magnitude は bare number（newtype 透過）: {}",
            greatest["magnitude"]
        );
        assert_eq!(
            greatest["obscuration"]
                .as_f64()
                .expect("obscuration は数値"),
            1.0,
            "greatest.obscuration（透過数値）"
        );
        assert!(
            greatest["obscuration"].is_number(),
            "obscuration は bare number（newtype 透過）: {}",
            greatest["obscuration"]
        );

        // 改名フィールド（単位サフィックス付き）。
        assert_eq!(
            greatest["path_width_km"]
                .as_f64()
                .expect("path_width_km は数値（皆既）"),
            197.0,
            "greatest.path_width_km（皆既）"
        );
        assert_eq!(
            greatest["central_duration_seconds"]
                .as_f64()
                .expect("central_duration_seconds は数値（皆既）"),
            268.0,
            "greatest.central_duration_seconds（皆既）"
        );
        assert_eq!(
            greatest["sun_altitude_deg"]
                .as_f64()
                .expect("sun_altitude_deg は数値"),
            70.3,
            "greatest.sun_altitude_deg（皆既）"
        );

        // metadata.
        let metadata = &obj["metadata"];
        assert_eq!(
            metadata["ephemeris_model"],
            Value::from("ELP/MPP02+VSOP87D"),
            "metadata.ephemeris_model"
        );
        assert_eq!(
            metadata["delta_t_model"],
            Value::from("EspenakMeeus"),
            "metadata.delta_t_model"
        );
        assert_eq!(
            metadata["delta_t_uncertainty_seconds"]
                .as_f64()
                .expect("delta_t_uncertainty_seconds は数値"),
            0.5,
            "metadata.delta_t_uncertainty_seconds"
        );
        assert_eq!(
            metadata["accuracy_profile"]["type"],
            Value::from("Standard"),
            "metadata.accuracy_profile は {{type:\"Standard\"}} タグ付き"
        );
        assert!(
            metadata["generated_at"]["iso"].is_string(),
            "metadata.generated_at.iso は文字列（UtcInstant オブジェクト形）: {}",
            metadata["generated_at"]
        );

        // bessel.
        let bessel = &obj["bessel"];
        assert_eq!(
            bessel["tan_f1"].as_f64().expect("tan_f1 は数値"),
            0.00465,
            "bessel.tan_f1"
        );
        assert_eq!(
            bessel["x"]["coefficients"],
            Value::from(vec![Value::from(0.20)]),
            "bessel.x.coefficients == [0.20]"
        );
        assert_eq!(
            bessel["fit_error"]["max_x"]
                .as_f64()
                .expect("fit_error.max_x は数値"),
            1.0e-7,
            "bessel.fit_error.max_x"
        );
        assert!(
            bessel["fit_interval"]["start"]["jd"]["part1"].is_number(),
            "bessel.fit_interval.start.jd.part1 は数値（TtInstant オブジェクト形）: {}",
            bessel["fit_interval"]["start"]
        );
    }

    /// 部分食 fixture の特有契約（中心食フィールドが null・magnitude<1・kind タグ Partial）。
    /// 殺す変異: None を null でなく省略/0 で出す、path_width_km/central_duration_seconds を皆既値で埋める、
    ///   central_begin/central_end の null を欠落させる、kind タグを Total/未タグにする。
    #[test]
    fn format_search_json_partial_fixture_contract() {
        let arr = parse_json_array(&[partial_eclipse()]);
        assert_eq!(arr.len(), 1, "部分 1 件 → 配列長 1");
        let obj = &arr[0];

        assert_eq!(
            obj["event_key"],
            Value::from("2025-03-29#1264"),
            "event_key（部分）"
        );
        assert_eq!(
            obj["kind"]["type"],
            Value::from("Partial"),
            "kind は {{type:\"Partial\"}}"
        );

        let greatest = &obj["global"]["greatest"];
        assert_eq!(
            greatest["magnitude"].as_f64().expect("magnitude は数値"),
            0.938,
            "magnitude（部分）"
        );
        // 中心食フィールドは None → null。
        assert_eq!(
            greatest["path_width_km"],
            Value::Null,
            "path_width_km は null（部分は中心食でない）"
        );
        assert_eq!(
            greatest["central_duration_seconds"],
            Value::Null,
            "central_duration_seconds は null（部分）"
        );

        let global = &obj["global"];
        assert_eq!(
            global["central_begin"],
            Value::Null,
            "global.central_begin は null（部分）"
        );
        assert_eq!(
            global["central_end"],
            Value::Null,
            "global.central_end は null（部分）"
        );
    }

    /// 複数 fixture（皆既→部分）で配列長 2・順序保存（全要素が順序通り出力される）。
    /// 殺す変異: 1 件目だけ/末尾だけ出す、順序を入れ替える、2 件目を落とす。
    #[test]
    fn format_search_json_preserves_order_and_every_element() {
        let arr = parse_json_array(&[total_eclipse(), partial_eclipse()]);
        assert_eq!(arr.len(), 2, "2 件入力 → 配列長 2");
        assert_eq!(
            arr[0]["event_key"],
            Value::from("2024-04-08#1252"),
            "要素 0 は皆既 event_key（順序保存）"
        );
        assert_eq!(
            arr[1]["event_key"],
            Value::from("2025-03-29#1264"),
            "要素 1 は部分 event_key（順序保存）"
        );
    }

    /// 空入力は長さ 0 の JSON 配列（空出力・架空日食を捏造しない）。パース成功必須。
    /// 殺す変異: 空入力で架空の日食を 1 件出す、null/オブジェクトを出す、パース不能な出力。
    #[test]
    fn format_search_json_empty_input_is_empty_array() {
        let arr = parse_json_array(&[]);
        assert_eq!(arr.len(), 0, "空入力 → 長さ 0 の配列（捏造しない）");
    }

    /// kind タグの安定性: kind は key "type"・値はバリアント名のオブジェクト（Total/Partial を網羅）。
    /// bare string / untagged 表現への回帰を陽に撃破する（time_utc 等の他フィールドと独立に縛る）。
    /// 殺す変異: kind を bare string にする、untagged にする、タグ key を "type" 以外にする、
    ///   バリアント名を別表記（小文字/別名）にする。
    #[test]
    fn format_search_json_kind_is_internally_tagged_object() {
        let total = parse_json_array(&[total_eclipse()]);
        assert!(
            total[0]["kind"].is_object(),
            "kind はオブジェクト（bare string でない）: {}",
            total[0]["kind"]
        );
        assert_eq!(
            total[0]["kind"]["type"],
            Value::from("Total"),
            "kind.type == バリアント名 \"Total\""
        );

        let partial = parse_json_array(&[partial_eclipse()]);
        assert!(
            partial[0]["kind"].is_object(),
            "kind はオブジェクト（部分）: {}",
            partial[0]["kind"]
        );
        assert_eq!(
            partial[0]["kind"]["type"],
            Value::from("Partial"),
            "kind.type == バリアント名 \"Partial\""
        );
    }

    /// run_search ディスパッチ: format=Json でも入力検証（不正日付）を fast-fail でバイパスしない。
    /// 殺す変異: JSON 経路で不正日付検証を飛ばす、別 variant を返す。
    #[test]
    fn run_search_json_invalid_from_still_fast_fails() {
        let args = SearchArgs {
            from: "not-a-date".to_string(),
            to: "2020-12-31".to_string(),
            accuracy: AccuracyArg::Standard,
            kind: KindFilter::All,
            format: FormatArg::Json,
        };
        let r = run_search(&args);
        assert!(
            matches!(r, Err(CliError::InvalidDate(ref s)) if s == "not-a-date"),
            "JSON 経路でも不正 from は InvalidDate で fast-fail: {r:?}"
        );
    }

    /// run_search ディスパッチ: format=Json でも from>to を fast-fail（RangeOrder）でバイパスしない。
    /// 殺す変異: JSON 経路で順序検証を飛ばす、別 variant を返す。
    #[test]
    fn run_search_json_from_after_to_still_range_order() {
        let args = SearchArgs {
            from: "2020-12-31".to_string(),
            to: "2020-01-01".to_string(),
            accuracy: AccuracyArg::Standard,
            kind: KindFilter::All,
            format: FormatArg::Json,
        };
        let r = run_search(&args);
        match r {
            Err(CliError::RangeOrder { from, to }) => {
                assert_eq!(from, "2020-12-31", "RangeOrder.from（JSON 経路）");
                assert_eq!(to, "2020-01-01", "RangeOrder.to（JSON 経路）");
            }
            other => panic!("expected Err(RangeOrder {{..}}) on JSON path, got {other:?}"),
        }
    }

    /// 【SLOW・1 件】正常系: 2017-08・Standard・All・format=Json で `Ok(s)`。s は有効な JSON 配列で
    /// パースでき、少なくとも 1 要素の `["kind"]["type"]` が `"Total"`（2017 皆既）。
    /// 内部で実エンジン（≈分）を実走するため SLOW。red 段階では未解決シンボルでコンパイル不能。
    /// 殺す変異: JSON 経路でエンジンを実走しない、整形を text に流す、kind タグを落とす。
    // SLOW
    #[test]
    fn run_search_json_2017_total_eclipse_is_valid_json_array() {
        let args = SearchArgs {
            from: "2017-08-01".to_string(),
            to: "2017-09-01".to_string(),
            accuracy: AccuracyArg::Standard,
            kind: KindFilter::All,
            format: FormatArg::Json,
        };
        let s = run_search(&args).expect("2017-08 の JSON 探索は成功する");
        let v: Value = serde_json::from_str(&s).expect("出力は有効な JSON（パース成功必須）");
        let arr = v.as_array().expect("トップレベルは JSON 配列");
        assert!(!arr.is_empty(), "2017-08 には少なくとも 1 件の日食");
        assert!(
            arr.iter().any(|e| e["kind"]["type"] == "Total"),
            "少なくとも 1 要素の kind.type が \"Total\"（2017 皆既）: {s}"
        );
    }

    // ==================================================================
    // === S32a: umbra local (text) ===
    // ==================================================================
    // ISSUE-032 S32a 受け入れテスト（standard・`umbra local` text 出力・--format なし）。
    //
    // ## オラクル戦略（実装方針に立ち入らず、確定仕様の公開 IF だけを縛る）
    // - **parse_utc_offset**: 期待分値を独立算術（例 9*60）で組む（実装の文字列処理を写経しない）。
    //   明確に妥当な offset は Ok(分)、明確なゴミ（"foo"/""）は InvalidTimezone（入力保持）。
    // - **format_local_text**: 既知 fixture（中心地点 FullyVisible／部分地点）を構造体リテラルで
    //   組み、出力に既知値が部分文字列として出る存在確認（厳密レイアウトは縛らない）。
    //   C2/C3 が None の地点は em ダッシュ "—" で描かれ、架空の接触時刻を捏造しないことを縛る。
    //   timezone Some 時はローカル時刻ラベルが追加されつつ UTC/TT が温存されることを縛る。
    // - **run_local（fast-fail）**: 不正 date／緯度範囲外（Domain）／不正 timezone は **エンジン実走前**
    //   に Err（高速）。red 段階では未解決シンボルでコンパイル不能。
    // - **西経吸収**: コアの正規化契約（−100° ≡ 260°）を Observer 等値で固定（CLI が依存する契約）。
    // - **run_local（SLOW・1〜2 件）**: 実エンジンで 2024-04-08 の皆既路上地点を解く（≈分）。
    //
    // ## red 設計（本体未実装）
    // `Command::Local`/`LocalArgs`/`RefractionArg`/`run_local`/`parse_utc_offset`/`format_local_text`、
    // および `CliError::{InvalidTimezone, Domain}` は本スライスで導入予定で現状未実装。
    // テストはコンパイル時点で未解決シンボル（red）。

    use umbra_core::{DomainError, Observer};
    use umbra_eclipse::{LocalCircumstances, LocalContact, LocalContactSet, Visibility};

    // ------------------------------------------------------------------
    // S32a 構築ヘルパ（局地接触・局地条件 fixture）
    // ------------------------------------------------------------------

    /// 局地接触点を既知値で組む（time_utc/time_tt は互いに区別できる別値）。
    fn local_contact(
        utc_h: u8,
        utc_mi: u8,
        tt_frac: f64,
        alt: f64,
        az: f64,
        pa: f64,
        visible: bool,
    ) -> LocalContact {
        LocalContact {
            time_utc: utc(2024, 4, 8, utc_h, utc_mi, 0.0),
            time_tt: tt(2_460_409.0, tt_frac),
            sun_altitude: Degrees(alt),
            sun_azimuth: Degrees(az),
            position_angle: Degrees(pa),
            visible,
        }
    }

    /// 中心地点（FullyVisible）の局地条件 fixture。c1..c4 すべて Some・各々別時刻。
    fn fully_visible_circ() -> LocalCircumstances {
        LocalCircumstances {
            contacts: LocalContactSet {
                c1: Some(local_contact(17, 18, 0.111, 60.0, 120.0, 250.0, true)),
                c2: Some(local_contact(18, 32, 0.222, 68.0, 150.0, 260.0, true)),
                maximum: local_contact(18, 34, 0.234, 70.5, 152.0, 265.0, true),
                c3: Some(local_contact(18, 36, 0.246, 69.0, 154.0, 80.0, true)),
                c4: Some(local_contact(20, 1, 0.333, 45.0, 230.0, 90.0, true)),
            },
            magnitude: EclipseMagnitude(1.0123),
            obscuration: Obscuration(1.0),
            maximum_altitude: Degrees(70.5),
            visibility: Visibility::FullyVisible,
            metadata: metadata(),
        }
    }

    /// 部分地点（PartialVisible）の局地条件 fixture。c2/c3 は None（中心接触なし）。
    fn partial_visible_circ() -> LocalCircumstances {
        LocalCircumstances {
            contacts: LocalContactSet {
                c1: Some(local_contact(17, 45, 0.150, 30.0, 100.0, 240.0, true)),
                c2: None,
                maximum: local_contact(19, 0, 0.270, 22.0, 200.0, 255.0, true),
                c3: None,
                c4: Some(local_contact(20, 15, 0.390, 8.0, 250.0, 95.0, true)),
            },
            magnitude: EclipseMagnitude(0.62),
            obscuration: Obscuration(0.51),
            maximum_altitude: Degrees(22.0),
            visibility: Visibility::PartialVisible,
            metadata: metadata(),
        }
    }

    // ------------------------------------------------------------------
    // 1. parse_utc_offset（FAST・独立算術オラクル）
    // ------------------------------------------------------------------

    /// 妥当な UTC offset 各表記が **符号付き分**へ正しくパースされる（期待値は独立算術で構成）。
    /// "+HH:MM"/"-HH:MM"/"+HHMM"/"-HHMM"/"Z" の全表記と符号・コロン有無を 1 件ずつ縛る。
    /// 殺す変異: 符号を無視する、コロン有無で分岐を誤る、時↔分の取り違え、Z を 0 にしない。
    #[test]
    fn parse_utc_offset_valid_forms_to_signed_minutes() {
        assert_eq!(
            parse_utc_offset("+09:00").expect("妥当 offset"),
            9 * 60,
            "+09:00 = +540 分"
        );
        assert_eq!(
            parse_utc_offset("-05:30").expect("妥当 offset"),
            -(5 * 60 + 30),
            "-05:30 = -330 分"
        );
        assert_eq!(parse_utc_offset("Z").expect("Z は 0"), 0, "Z = 0 分");
        assert_eq!(
            parse_utc_offset("+00:00").expect("妥当 offset"),
            0,
            "+00:00 = 0 分"
        );
        assert_eq!(
            parse_utc_offset("+0000").expect("妥当 offset"),
            0,
            "+0000 = 0 分"
        );
        assert_eq!(
            parse_utc_offset("+0900").expect("妥当 offset"),
            9 * 60,
            "+0900 = +540 分（コロン無し）"
        );
        assert_eq!(
            parse_utc_offset("-0500").expect("妥当 offset"),
            -(5 * 60),
            "-0500 = -300 分（コロン無し）"
        );
    }

    /// 明確なゴミ文字列は `Err(CliError::InvalidTimezone(入力))`（入力文字列をそのまま保持）。
    /// 殺す変異: 不正入力を 0 等にフォールバックして Ok を返す、InvalidTimezone に別文字列を載せる。
    #[test]
    fn parse_utc_offset_garbage_is_invalid_timezone() {
        match parse_utc_offset("foo") {
            Err(CliError::InvalidTimezone(s)) => {
                assert_eq!(s, "foo", "InvalidTimezone は入力 \"foo\" を保持");
            }
            other => panic!("expected Err(InvalidTimezone(\"foo\")), got {other:?}"),
        }
        match parse_utc_offset("") {
            Err(CliError::InvalidTimezone(s)) => {
                assert_eq!(s, "", "InvalidTimezone は空入力を保持");
            }
            other => panic!("expected Err(InvalidTimezone(\"\")), got {other:?}"),
        }
    }

    // ------------------------------------------------------------------
    // 2. format_local_text（FAST・既知 fixture の内容存在で縛る）
    // ------------------------------------------------------------------

    /// 中心地点（FullyVisible・timezone なし）を整形すると、各接触の UTC 暦日と TT ラベルが
    /// 併記され、食分・食面積・最大高度・可視性名・計算メタデータ（暦/ΔT モデル名・不確実性帯）が
    /// 出力に含まれる（accuracy.md §0: UTC+TT 併記必須）。
    /// 殺す変異: TT を接触行から落とし UTC のみ出す、食分/食面積/最大高度/可視性名/メタデータの欠落。
    #[test]
    fn format_local_text_central_site_contains_key_fields() {
        let out = format_local_text(&fully_visible_circ(), None);

        // 接触の UTC 暦日（2024-04-08）。
        assert!(
            out.contains("2024-04-08"),
            "接触の UTC 暦日が出力に含まれる: {out}"
        );
        // TT 併記（TT ラベル）。
        assert!(out.contains("TT"), "接触の TT 併記が出力に含まれる: {out}");
        // 食分（1.0123）。
        assert!(
            out.contains("1.0123"),
            "食分 1.0123 が出力に含まれる: {out}"
        );
        // 食面積（1.0）。
        assert!(out.contains("1.0"), "食面積 1.0 が出力に含まれる: {out}");
        // 最大高度（70.5）。
        assert!(
            out.contains("70.5"),
            "最大高度 70.5 が出力に含まれる: {out}"
        );
        // 可視性名（FullyVisible）。
        assert!(
            out.contains("FullyVisible"),
            "可視性名 FullyVisible が出力に含まれる: {out}"
        );
        // 計算メタデータ: 暦モデル名・ΔT モデル名・ΔT 不確実性帯。
        assert!(
            out.contains("ELP/MPP02+VSOP87D"),
            "ephemeris モデル名が出力に含まれる: {out}"
        );
        assert!(
            out.contains("EspenakMeeus"),
            "ΔT モデル名が出力に含まれる: {out}"
        );
        assert!(
            out.contains("0.5"),
            "ΔT 不確実性帯（0.5）が出力に含まれる: {out}"
        );
    }

    /// 部分地点（c2/c3 が None）を整形すると、欠落接触は em ダッシュ "—"（U+2014）で描かれ、
    /// 架空の中心接触時刻を捏造しない。非 Option の maximum は常に出る。可視性名も出る（非 panic）。
    /// 殺す変異: c2/c3 の None を "—" で表示しない、None なのに架空時刻を出す、maximum を落とす、
    ///   部分地点で panic する。
    #[test]
    fn format_local_text_partial_site_shows_dash_for_none_contacts() {
        let circ = partial_visible_circ();
        // fixture 前提: maximum ≈ 19:00 UTC（架空 c2/c3 時刻と区別する材料）。to_gregorian は丸め
        // 境界で 18:59:59.9995 を返しうる（S31b で判明した暦往復の ±eps）ため、暦成分の厳密比較で
        // なく JD レベルで ±1 分以内を確認する。
        let max_jd = circ.contacts.maximum.time_utc.jd2().jd();
        let expected_jd = utc(2024, 4, 8, 19, 0, 0.0).jd2().jd();
        assert!(
            (max_jd - expected_jd).abs() < 1.0 / 1440.0,
            "fixture 前提: maximum ≈ 19:00 UTC（|Δ| < 1 分）"
        );

        let out = format_local_text(&circ, None);

        // 欠落接触（c2/c3）は em ダッシュで描かれる。
        assert!(
            out.contains('—'),
            "None 接触は em ダッシュ '—'(U+2014) で描かれる: {out}"
        );
        // 非 Option の maximum は常に存在（最大食 UTC 暦日が出る）。
        assert!(
            out.contains("2024-04-08"),
            "maximum（非 Option）の UTC 暦日が出力に含まれる: {out}"
        );
        // 可視性名（PartialVisible）。
        assert!(
            out.contains("PartialVisible"),
            "部分地点の可視性名 PartialVisible が出力に含まれる: {out}"
        );
    }

    /// timezone Some 時はローカル時刻ラベル（"+09:00"）が追加されつつ、内部の UTC/TT は温存される。
    /// 殺す変異: timezone ラベルを描かない、ローカル表示で UTC を消す、TT を落とす。
    #[test]
    fn format_local_text_with_timezone_adds_label_and_keeps_utc_tt() {
        let out = format_local_text(&fully_visible_circ(), Some(("+09:00", 540)));

        // ローカル時刻ラベルが追加される。
        assert!(
            out.contains("+09:00"),
            "timezone ラベル +09:00 が出力に含まれる: {out}"
        );
        // UTC は温存（接触の UTC 暦日が依然として出る）。
        assert!(
            out.contains("2024-04-08"),
            "ローカル表示でも UTC 暦日は温存される: {out}"
        );
        // TT も温存。
        assert!(
            out.contains("TT"),
            "ローカル表示でも TT 併記は温存される: {out}"
        );
    }

    // ------------------------------------------------------------------
    // 3. run_local（fast-fail = FAST / 正常系 = SLOW）
    // ------------------------------------------------------------------

    /// 不正 `--date` はエンジン実走前に `Err(CliError::InvalidDate(入力))`（fast-fail）。
    /// 殺す変異: 不正 date を黙って受理してエンジンに進む、InvalidDate でなく別 variant を返す。
    #[test]
    fn run_local_invalid_date_is_invalid_date_error() {
        let args = LocalArgs {
            date: "not-a-date".to_string(),
            lat: 35.0,
            lon: 139.0,
            elevation: 0.0,
            timezone: None,
            accuracy: AccuracyArg::Standard,
            refraction: RefractionArg::Standard,
            format: FormatArg::Text,
        };
        let r = run_local(&args);
        assert!(
            matches!(r, Err(CliError::InvalidDate(ref s)) if s == "not-a-date"),
            "expected Err(InvalidDate(\"not-a-date\")), got {r:?}"
        );
    }

    /// 緯度範囲外（91.0）は Observer 構築で弾かれ `Err(CliError::Domain(_))`（エンジン実走前）。
    /// 殺す変異: 緯度範囲チェックを迂回してエンジンに進む、Domain でなく別 variant を返す。
    #[test]
    fn run_local_latitude_out_of_range_is_domain_error() {
        let args = LocalArgs {
            date: "2024-04-08".to_string(),
            lat: 91.0,
            lon: 0.0,
            elevation: 0.0,
            timezone: None,
            accuracy: AccuracyArg::Standard,
            refraction: RefractionArg::Standard,
            format: FormatArg::Text,
        };
        let r = run_local(&args);
        assert!(
            matches!(r, Err(CliError::Domain(_))),
            "expected Err(Domain(_)) for lat=91.0, got {r:?}"
        );
    }

    /// 不正 `--timezone`（"foo"）は date/lat/lon が妥当でもエンジン実走前に
    /// `Err(CliError::InvalidTimezone("foo"))`。
    /// 殺す変異: timezone 検証をエンジン後に回す/省く、InvalidTimezone でなく別 variant を返す。
    #[test]
    fn run_local_invalid_timezone_is_invalid_timezone_error() {
        let args = LocalArgs {
            date: "2024-04-08".to_string(),
            lat: 35.0,
            lon: 139.0,
            elevation: 0.0,
            timezone: Some("foo".to_string()),
            accuracy: AccuracyArg::Standard,
            refraction: RefractionArg::Standard,
            format: FormatArg::Text,
        };
        let r = run_local(&args);
        assert!(
            matches!(r, Err(CliError::InvalidTimezone(ref s)) if s == "foo"),
            "expected Err(InvalidTimezone(\"foo\")), got {r:?}"
        );
    }

    // ------------------------------------------------------------------
    // 4. 西経吸収（FAST・コア正規化契約・型レベル）
    // ------------------------------------------------------------------

    /// CLI が依存するコア契約: 西経 −100° は東経 260° と正規化後に等価な Observer になる。
    /// 殺す変異: 経度正規化を行わず西経入力が別地点になる（CLI の lon 吸収契約破壊）。
    #[test]
    fn observer_west_longitude_absorbs_to_equivalent_east() {
        let west = Observer::from_degrees(35.0, -100.0, 0.0).expect("妥当 Observer");
        let east = Observer::from_degrees(35.0, 260.0, 0.0).expect("妥当 Observer");
        assert_eq!(
            west, east,
            "西経 −100° ≡ 東経 260°（正規化後の Observer 等値）"
        );
        // DomainError 型の存在も固定（lat 範囲外は Domain で surface する契約の土台）。
        let _ = DomainError::OutOfRange {
            what: "geodetic latitude",
        };
    }

    /// clap が **負の緯度（南緯）・負の経度（西経）** を `--lat`/`--lon` の値として受理する
    /// （`-34.0` をフラグと誤認しない）。`allow_hyphen_values` の脱落＝南半球/西経の入力拒否を撃破。
    /// 殺す変異: `--lat` の `allow_hyphen_values` を外す（負緯度が clap パースエラーになる）。
    #[test]
    fn local_args_accept_negative_lat_lon_via_clap() {
        let cli = Cli::try_parse_from([
            "umbra",
            "local",
            "--date",
            "2024-04-08",
            "--lat",
            "-34.0",
            "--lon",
            "-58.0",
        ])
        .expect("負の緯度・経度は clap が値として受理する");
        match cli.command {
            Command::Local(args) => {
                assert_eq!(args.lat, -34.0, "南緯 -34.0 が値として渡る");
                assert_eq!(args.lon, -58.0, "西経 -58.0 が値として渡る");
            }
            other => panic!("expected Command::Local, got {other:?}"),
        }
    }

    // ------------------------------------------------------------------
    // 5. run_local（SLOW・実エンジン）
    // ------------------------------------------------------------------

    /// 【SLOW・1 件】正常系: 2024-04-08・米中部 Texas の皆既路上地点・timezone -05:00 で `Ok(s)`。
    /// s は日食日付 "2024-04-08"・可視性名のいずれか・timezone ラベル "-05:00" を含む（構造的存在
    /// のみ・計算値は過度に縛らない＝コアの精度は別所でテスト済み）。実エンジン実走で SLOW。
    /// 殺す変異: 観測地点をエンジンに渡さない、整形を呼ばない、timezone ラベルを出さない。
    // SLOW
    #[test]
    fn run_local_2024_total_path_site_produces_output() {
        let args = LocalArgs {
            date: "2024-04-08".to_string(),
            lat: 30.0,
            lon: -98.0,
            elevation: 0.0,
            timezone: Some("-05:00".to_string()),
            accuracy: AccuracyArg::Standard,
            refraction: RefractionArg::Standard,
            format: FormatArg::Text,
        };
        let out = run_local(&args).expect("2024-04-08 Texas の局地計算は成功する");
        assert!(
            out.contains("2024-04-08"),
            "日食日付 2024-04-08 が出力に含まれる: {out}"
        );
        // 6 値のいずれかの可視性名が出る（どれかは地点依存・少なくとも 1 つ）。
        let visibility_names = [
            "NotVisible",
            "BelowHorizon",
            "SunriseEclipse",
            "SunsetEclipse",
            "PartialVisible",
            "FullyVisible",
        ];
        assert!(
            visibility_names.iter().any(|n| out.contains(n)),
            "可視性名（6 値のいずれか）が出力に含まれる: {out}"
        );
        assert!(
            out.contains("-05:00"),
            "timezone ラベル -05:00 が出力に含まれる: {out}"
        );
    }

    /// 【SLOW・1 件】日食なし日付（2024-06-15）は `Ok(s)`・s は「日食なし」を示し、架空の event_key
    /// （`#` 付き安定キー）を捏造しない。実エンジン実走で SLOW。
    /// 殺す変異: 日食なし日に架空イベントを出す、no-eclipse メッセージを出さず panic/Err。
    // SLOW
    #[test]
    fn run_local_no_eclipse_date_reports_no_eclipse() {
        let args = LocalArgs {
            date: "2024-06-15".to_string(),
            lat: 35.0,
            lon: 139.0,
            elevation: 0.0,
            timezone: None,
            accuracy: AccuracyArg::Standard,
            refraction: RefractionArg::Standard,
            format: FormatArg::Text,
        };
        let out = run_local(&args).expect("日食なし日付でも Ok（メッセージ）を返す");
        // 架空の安定キー（"YYYY-MM-DD#NNNN" の '#'）を捏造しない。
        assert!(
            !out.contains('#'),
            "日食なし日に架空 event_key（'#'）を捏造しない: {out}"
        );
    }

    // ==================================================================
    // === S32b: umbra local --format json ===
    // ==================================================================
    // ISSUE-032 S32b 受け入れテスト（standard・`umbra local --format json`）。
    //
    // ## オラクル戦略
    // 出力文字列を `serde_json::from_str::<serde_json::Value>` でパースし、**パースされた
    // Value** に対して構造を assert する（生 JSON のレイアウト一致は使わない。部分文字列は
    // 二次シグナルとしてのみ許容）。値は既知 fixture（fully_visible_circ / partial_visible_circ /
    // total_eclipse）由来で、凍結された観測可能 JSON 契約をそのまま縛る。S31a→S31b と同形で
    // local 結果型グラフ（VisibleSolarEclipse → LocalCircumstances → LocalContactSet/
    // LocalContact など）に Serialize を横断配線し serde_json で整形する。
    //
    // ## red 設計（本体未実装）
    // `LocalArgs.format`・`format_local_json`・local 結果型の Serialize 配線は本スライスで導入予定で
    // 現状未実装。テストはコンパイル時点で未解決シンボル／未実装トレイト境界（red）。

    use umbra_eclipse::VisibleSolarEclipse;

    /// `format_local_json(Some(..))` 出力をパースして JSON Value を返すヘルパ（パース成功＝有効 JSON）。
    fn parse_local_json_value(vse: &VisibleSolarEclipse) -> Value {
        let s = format_local_json(Some(vse)).expect("JSON 整形は成功する");
        serde_json::from_str(&s).expect("出力は有効な JSON（パース成功必須）")
    }

    /// FullyVisible fixture の JSON 契約を全面的に縛る。`obj["eclipse"]`（SolarEclipse 形・軽い確認）と
    /// `obj["local"]`（magnitude/obscuration の透過数値・maximum_altitude_deg 改名・visibility の
    /// 内部タグ・metadata・contacts.maximum の instant 形＋`_deg` 改名キー＋未サフィックスキー不在＋
    /// visible bool・c1..c4 全 Some）を検証する。
    /// 殺す変異: eclipse/local の二分割を崩す、visibility を bare string にする、
    ///   maximum_altitude を `_deg` 改名しない、maximum 接触の alt/az/PA を `_deg` 改名しない、
    ///   maximum を null/省略にする、c1..c4 を Some で出さない、time_utc/time_tt の iso/jd を落とす、
    ///   UTC iso の末尾 Z を落とす・TT iso に Z を付ける、magnitude/obscuration をオブジェクト化する。
    #[test]
    fn format_local_json_some_fully_visible_contract() {
        let vse = VisibleSolarEclipse {
            eclipse: total_eclipse(),
            local: fully_visible_circ(),
        };
        let obj = parse_local_json_value(&vse);
        assert!(
            obj.is_object(),
            "Some(..) は JSON オブジェクトにシリアライズ: {obj}"
        );

        // --- eclipse（S31b で確立済みの SolarEclipse 形・軽い確認のみ）---
        let eclipse = &obj["eclipse"];
        assert!(
            eclipse["event_key"].is_string(),
            "eclipse.event_key は文字列: {eclipse}"
        );
        assert_eq!(
            eclipse["event_key"],
            Value::from("2024-04-08#1252"),
            "eclipse.event_key（total_eclipse）"
        );
        assert!(
            !eclipse["kind"]["type"].is_null(),
            "eclipse.kind.type が存在（タグ付き enum）: {eclipse}"
        );

        // --- local ---
        let local = &obj["local"];

        // magnitude / obscuration は透過 newtype = bare number。
        assert!(
            local["magnitude"].is_number(),
            "local.magnitude は bare number（newtype 透過）: {}",
            local["magnitude"]
        );
        assert_eq!(
            local["magnitude"].as_f64().expect("magnitude は数値"),
            1.0123,
            "local.magnitude（FullyVisible fixture）"
        );
        assert!(
            local["obscuration"].is_number(),
            "local.obscuration は bare number（newtype 透過）: {}",
            local["obscuration"]
        );
        assert_eq!(
            local["obscuration"].as_f64().expect("obscuration は数値"),
            1.0,
            "local.obscuration（FullyVisible fixture）"
        );

        // maximum_altitude_deg（A7 単位サフィックス改名）。未改名 maximum_altitude は不在。
        assert_eq!(
            local["maximum_altitude_deg"]
                .as_f64()
                .expect("maximum_altitude_deg は数値"),
            70.5,
            "local.maximum_altitude_deg（改名・単位付き）"
        );
        assert!(
            local["maximum_altitude"].is_null(),
            "未改名キー maximum_altitude は不在（_deg 改名漏れを撃破）: {}",
            local["maximum_altitude"]
        );

        // visibility は内部タグ付き enum オブジェクト {type:"FullyVisible"}（bare string でない）。
        assert!(
            local["visibility"].is_object(),
            "local.visibility はオブジェクト（bare string でない）: {}",
            local["visibility"]
        );
        assert_eq!(
            local["visibility"]["type"],
            Value::from("FullyVisible"),
            "local.visibility は {{type:\"FullyVisible\"}} タグ付き"
        );

        // metadata（metadata() 由来・S31b と同形）。
        let metadata = &local["metadata"];
        assert_eq!(
            metadata["ephemeris_model"],
            Value::from("ELP/MPP02+VSOP87D"),
            "local.metadata.ephemeris_model"
        );
        assert_eq!(
            metadata["accuracy_profile"]["type"],
            Value::from("Standard"),
            "local.metadata.accuracy_profile は {{type:\"Standard\"}} タグ付き"
        );

        // --- contacts ---
        let contacts = &local["contacts"];

        // maximum は非 Option（常にオブジェクト）。
        let maximum = &contacts["maximum"];
        assert!(
            maximum.is_object(),
            "contacts.maximum は常にオブジェクト（非 Option）: {maximum}"
        );

        // maximum.time_utc / time_tt は {iso, jd:{part1,part2}} の instant 形。
        let utc_iso = maximum["time_utc"]["iso"]
            .as_str()
            .expect("maximum.time_utc.iso は文字列");
        assert!(
            utc_iso.ends_with('Z'),
            "maximum.time_utc.iso は UTC（末尾 Z）: {utc_iso}"
        );
        assert!(
            maximum["time_utc"]["jd"]["part1"].is_number(),
            "maximum.time_utc.jd.part1 は数値（lossless チャネル）: {}",
            maximum["time_utc"]["jd"]["part1"]
        );
        assert!(
            maximum["time_utc"]["jd"]["part2"].is_number(),
            "maximum.time_utc.jd.part2 は数値: {}",
            maximum["time_utc"]["jd"]["part2"]
        );
        let tt_iso = maximum["time_tt"]["iso"]
            .as_str()
            .expect("maximum.time_tt.iso は文字列");
        assert!(
            !tt_iso.ends_with('Z'),
            "maximum.time_tt.iso は末尾 Z を持たない（TT は UTC でない）: {tt_iso}"
        );
        assert!(
            maximum["time_tt"]["jd"]["part1"].is_number(),
            "maximum.time_tt.jd.part1 は数値: {}",
            maximum["time_tt"]["jd"]["part1"]
        );
        assert!(
            maximum["time_tt"]["jd"]["part2"].is_number(),
            "maximum.time_tt.jd.part2 は数値: {}",
            maximum["time_tt"]["jd"]["part2"]
        );

        // maximum の角度は `_deg` 改名（A7）・visible は bool。
        assert!(
            maximum["sun_altitude_deg"].is_number(),
            "maximum.sun_altitude_deg は数値（改名）: {}",
            maximum["sun_altitude_deg"]
        );
        assert!(
            maximum["sun_azimuth_deg"].is_number(),
            "maximum.sun_azimuth_deg は数値（改名）: {}",
            maximum["sun_azimuth_deg"]
        );
        assert!(
            maximum["position_angle_deg"].is_number(),
            "maximum.position_angle_deg は数値（改名）: {}",
            maximum["position_angle_deg"]
        );
        assert!(
            maximum["visible"].is_boolean(),
            "maximum.visible は bool: {}",
            maximum["visible"]
        );
        // 未改名キー（_deg なし）は maximum 接触上に不在。
        assert!(
            maximum["sun_altitude"].is_null(),
            "未改名キー sun_altitude は不在（_deg 改名漏れを撃破）: {}",
            maximum["sun_altitude"]
        );
        assert!(
            maximum["sun_azimuth"].is_null(),
            "未改名キー sun_azimuth は不在: {}",
            maximum["sun_azimuth"]
        );
        assert!(
            maximum["position_angle"].is_null(),
            "未改名キー position_angle は不在: {}",
            maximum["position_angle"]
        );

        // FullyVisible fixture: c1..c4 すべて Some（オブジェクト）。
        assert!(
            contacts["c1"].is_object(),
            "contacts.c1 は Some（オブジェクト）: {}",
            contacts["c1"]
        );
        assert!(
            contacts["c2"].is_object(),
            "contacts.c2 は Some（オブジェクト）: {}",
            contacts["c2"]
        );
        assert!(
            contacts["c3"].is_object(),
            "contacts.c3 は Some（オブジェクト）: {}",
            contacts["c3"]
        );
        assert!(
            contacts["c4"].is_object(),
            "contacts.c4 は Some（オブジェクト）: {}",
            contacts["c4"]
        );
    }

    /// PartialVisible fixture: 内側接触 c2/c3（None）が JSON null・c1/maximum/c4 はオブジェクト・
    /// visibility タグは PartialVisible。
    /// 殺す変異: None の内側接触を null でなく省略/0/架空時刻で出す、c1/maximum/c4 を落とす、
    ///   visibility タグを FullyVisible/未タグにする。
    #[test]
    fn format_local_json_partial_visible_has_null_inner_contacts() {
        let vse = VisibleSolarEclipse {
            eclipse: total_eclipse(),
            local: partial_visible_circ(),
        };
        let obj = parse_local_json_value(&vse);
        let local = &obj["local"];
        let contacts = &local["contacts"];

        // c2/c3 は None → JSON null。
        assert_eq!(
            contacts["c2"],
            Value::Null,
            "contacts.c2 は null（部分地点は中心接触なし）"
        );
        assert_eq!(
            contacts["c3"],
            Value::Null,
            "contacts.c3 は null（部分地点）"
        );

        // c1 / maximum / c4 はオブジェクト。
        assert!(
            contacts["c1"].is_object(),
            "contacts.c1 は Some（部分地点でも C1 あり）: {}",
            contacts["c1"]
        );
        assert!(
            contacts["maximum"].is_object(),
            "contacts.maximum は常にオブジェクト: {}",
            contacts["maximum"]
        );
        assert!(
            contacts["c4"].is_object(),
            "contacts.c4 は Some（部分地点でも C4 あり）: {}",
            contacts["c4"]
        );

        // visibility タグ。
        assert_eq!(
            local["visibility"]["type"],
            Value::from("PartialVisible"),
            "local.visibility は {{type:\"PartialVisible\"}}（部分地点）"
        );
    }

    /// 日食なし（None）は JSON リテラル `null` にシリアライズされ、出力は有効な JSON。
    /// 殺す変異: None を null でなく空オブジェクト/空文字/架空イベントで出す、パース不能な出力。
    #[test]
    fn format_local_json_none_is_json_null() {
        let s = format_local_json(None).expect("None でも JSON 整形は成功する");
        let v: Value = serde_json::from_str(&s).expect("出力は有効な JSON（パース成功必須）");
        assert!(
            v.is_null(),
            "format_local_json(None) は JSON null（日食なしを null で表す）: {v}"
        );
    }

    /// visibility タグ安定性: `local["visibility"]` は key "type"・値はバリアント名のオブジェクト
    /// （FullyVisible / PartialVisible を網羅）。bare string / untagged への回帰を陽に撃破する。
    /// 殺す変異: visibility を bare string にする、untagged にする、タグ key を "type" 以外にする、
    ///   バリアント名を別表記にする。
    #[test]
    fn format_local_json_visibility_is_internally_tagged() {
        let full = VisibleSolarEclipse {
            eclipse: total_eclipse(),
            local: fully_visible_circ(),
        };
        let full_obj = parse_local_json_value(&full);
        let full_vis = &full_obj["local"]["visibility"];
        assert!(
            full_vis.is_object(),
            "visibility はオブジェクト（bare string でない）: {full_vis}"
        );
        assert_eq!(
            full_vis["type"],
            Value::from("FullyVisible"),
            "visibility.type == バリアント名 \"FullyVisible\""
        );

        let partial = VisibleSolarEclipse {
            eclipse: total_eclipse(),
            local: partial_visible_circ(),
        };
        let partial_obj = parse_local_json_value(&partial);
        let partial_vis = &partial_obj["local"]["visibility"];
        assert!(
            partial_vis.is_object(),
            "visibility はオブジェクト（部分地点）: {partial_vis}"
        );
        assert_eq!(
            partial_vis["type"],
            Value::from("PartialVisible"),
            "visibility.type == バリアント名 \"PartialVisible\""
        );
    }

    /// run_local ディスパッチ: format=Json でも不正日付検証を fast-fail でバイパスしない
    /// （エンジン非実走＝高速）。
    /// 殺す変異: JSON 経路で不正日付検証を飛ばす、別 variant を返す。
    #[test]
    fn run_local_json_invalid_date_still_fast_fails() {
        let args = LocalArgs {
            date: "not-a-date".to_string(),
            lat: 35.0,
            lon: 139.0,
            elevation: 0.0,
            timezone: None,
            accuracy: AccuracyArg::Standard,
            refraction: RefractionArg::Standard,
            format: FormatArg::Json,
        };
        let r = run_local(&args);
        assert!(
            matches!(r, Err(CliError::InvalidDate(ref s)) if s == "not-a-date"),
            "JSON 経路でも不正 date は InvalidDate で fast-fail: {r:?}"
        );
    }

    /// run_local ディスパッチ: format=Json でも不正 timezone 検証を fast-fail でバイパスしない
    /// （date/lat/lon は妥当・エンジン非実走＝高速）。text 経路の対称テスト。
    /// 殺す変異: JSON 経路で timezone 検証をエンジン後に回す/省く、別 variant を返す。
    #[test]
    fn run_local_json_invalid_timezone_still_fast_fails() {
        let args = LocalArgs {
            date: "2024-04-08".to_string(),
            lat: 35.0,
            lon: 139.0,
            elevation: 0.0,
            timezone: Some("foo".to_string()),
            accuracy: AccuracyArg::Standard,
            refraction: RefractionArg::Standard,
            format: FormatArg::Json,
        };
        let r = run_local(&args);
        assert!(
            matches!(r, Err(CliError::InvalidTimezone(ref s)) if s == "foo"),
            "JSON 経路でも不正 timezone は InvalidTimezone で fast-fail: {r:?}"
        );
    }

    /// 【SLOW・1 件】正常系: 2024-04-08・米中部 Texas 路上地点・timezone -05:00・format=Json で
    /// `Ok(s)`。s は有効な JSON オブジェクトでパースでき、`["eclipse"]["event_key"]` が文字列・
    /// `["local"]["visibility"]["type"]` が存在する（計算値は過度に縛らない）。実エンジン実走で SLOW。
    /// 殺す変異: JSON 経路でエンジンを実走しない、整形を text に流す、eclipse/local 二分割を崩す、
    ///   visibility タグを落とす。
    // SLOW
    #[test]
    fn run_local_json_2024_total_path_site_is_valid_json() {
        let args = LocalArgs {
            date: "2024-04-08".to_string(),
            lat: 30.0,
            lon: -98.0,
            elevation: 0.0,
            timezone: Some("-05:00".to_string()),
            accuracy: AccuracyArg::Standard,
            refraction: RefractionArg::Standard,
            format: FormatArg::Json,
        };
        let s = run_local(&args).expect("2024-04-08 Texas の JSON 局地計算は成功する");
        let v: Value = serde_json::from_str(&s).expect("出力は有効な JSON（パース成功必須）");
        assert!(v.is_object(), "トップレベルは JSON オブジェクト: {s}");
        assert!(
            v["eclipse"]["event_key"].is_string(),
            "eclipse.event_key は文字列: {s}"
        );
        assert!(
            !v["local"]["visibility"]["type"].is_null(),
            "local.visibility.type が存在: {s}"
        );
    }

    // ==================================================================
    // === text 0.1s 丸め（:60.0 回帰防止） ===
    // ==================================================================
    // text 整形は秒を `{:04.1}` で印字するため、暦往復のドリフトで jd2_to_gregorian が
    // 16:00:00 に対し 15:59:59.9995 を返すと "15:59:60.0"（不正な :60 秒）が出る。
    // text 整形器は共有の 0.1 秒丸めヘルパ経由で秒を桁上げし、:60 を絶対に出さないことを縛る。
    //
    // ## オラクル戦略
    // 既知ドリフト時刻 utc(2024,4,8,16,0,0.0)（jd2_to_gregorian で 15:59:59.9995）を最大食/最大
    // 接触に据えた fixture を構造体リテラルで組み、**観測可能な text 出力**に "16:00:00" を含み
    // ":60"（および "15:59:60"）を含まないことを縛る。JSON 経路は S31b/S32b 済みなので再検証しない。
    //
    // ## red 設計
    // 整形器が依然 raw jd2_to_gregorian を使うため、出力に "15:59:60.0" が現れ実行時に FAIL する
    // （コンパイルは通る＝既存公開 IF のまま）。

    /// 既知ドリフト時刻を最大食に据えた皆既 fixture（total_eclipse の greatest 時刻のみ差し替え）。
    fn drift_greatest_eclipse() -> SolarEclipse {
        let greatest = GreatestEclipse {
            // jd2_to_gregorian で 15:59:59.9995 を返す既知ドリフト分境界。
            time_utc: utc(2024, 4, 8, 16, 0, 0.0),
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

    /// 既知ドリフト時刻を最大接触に据えた局地条件 fixture（maximum を 16:00 に差し替え）。
    fn drift_maximum_circ() -> LocalCircumstances {
        LocalCircumstances {
            contacts: LocalContactSet {
                c1: Some(local_contact(15, 30, 0.111, 60.0, 120.0, 250.0, true)),
                c2: None,
                // jd2_to_gregorian で 15:59:59.9995 を返す既知ドリフト分境界。
                maximum: local_contact(16, 0, 0.234, 70.5, 152.0, 265.0, true),
                c3: None,
                c4: Some(local_contact(16, 30, 0.333, 45.0, 230.0, 90.0, true)),
            },
            magnitude: EclipseMagnitude(1.0123),
            obscuration: Obscuration(1.0),
            maximum_altitude: Degrees(70.5),
            visibility: Visibility::FullyVisible,
            metadata: metadata(),
        }
    }

    /// format_search_text の最大食 UTC が 16:00:00 のとき、出力に "16:00:00" を含み ":60"
    /// （まして "15:59:60"）を含まない。raw jd2_to_gregorian の 15:59:59.9995 を `{:04.1}` で
    /// 印字すると "15:59:60.0" が出る回帰を撃破する。
    /// 殺す変異: text 経路で raw jd2_to_gregorian を使う（"15:59:60.0" 出力）、丸めヘルパを経由しない。
    #[test]
    fn format_search_text_rounds_drift_minute_no_sixty_seconds() {
        let out = format_search_text(&[drift_greatest_eclipse()]);
        assert!(
            out.contains("16:00:00"),
            "最大食 UTC は 16:00:00 に丸められる: {out}"
        );
        assert!(
            !out.contains(":60"),
            ":60 秒（不正な秒）を出力に含まない: {out}"
        );
        assert!(
            !out.contains("15:59:60"),
            "15:59:60 を出力に含まない（ドリフト丸め回帰）: {out}"
        );
    }

    /// format_local_text の最大接触 UTC が 16:00:00 のとき、出力に "16:00:00" を含み ":60" を
    /// 含まない。局地 text 経路でも raw jd2_to_gregorian による "15:59:60.0" 回帰を撃破する。
    /// 殺す変異: 局地 text 経路で raw jd2_to_gregorian を使う、丸めヘルパを経由しない。
    #[test]
    fn format_local_text_rounds_drift_minute_no_sixty_seconds() {
        let out = format_local_text(&drift_maximum_circ(), None);
        assert!(
            out.contains("16:00:00"),
            "最大接触 UTC は 16:00:00 に丸められる: {out}"
        );
        assert!(
            !out.contains(":60"),
            ":60 秒（不正な秒）を出力に含まない: {out}"
        );
        assert!(
            !out.contains("15:59:60"),
            "15:59:60 を出力に含まない（ドリフト丸め回帰）: {out}"
        );
    }

    // ==================================================================
    // === ISSUE-048: umbra path（text / geojson） ===
    // ==================================================================
    // ISSUE-048 受け入れテスト（`umbra path`・`run_path`）。
    //
    // ## オラクル戦略（外部表のハードコード禁止, conventions §11）
    // - 数値オラクルは **同じ公開 API から導出**する: 参照 `EclipsePath` を
    //   `standard_engine(...).search(1 日範囲) → path(options)` で独立に組み立て、run_path の
    //   出力（text の件数・geojson の role 集合）と突き合わせる。NASA 等の外部表は一切使わない。
    // - geojson は文字列一致でなく `serde_json::Value` へパースして構造を assert する
    //   （S31b/S32b と同じ方針）。
    // - 「日食の無い日」は 2024-06-15 を使う。選定理由: 2024 年の日食は 4 月と 10 月のみで 6 月中旬
    //   には起こらない、かつ **エンジン自身**（`engine_path("2024-06-15", ..) == None`）で検証する
    //   （テスト内で判定を実行するので外部表に依存しない）。既存 S32a
    //   `run_local_no_eclipse_date_reports_no_eclipse` と同一日付。
    // - text のラベル語彙は確定仕様 §5/§7（「中心線: なし」）を最小契約とし、
    //   `中心線`/`北限`/`南限`/`部分食域`/`samples` を含む **行が存在すること**を縛る
    //   （行ごと落とす実装を撃破するため、行の存在と内容を分けて assert する）。
    //
    // ## red 設計（本体未実装）
    // `PathArgs`/`PathFormatArg`/`run_path`/`CliError::InvalidInterval` は未導入のため、
    // 本節はコンパイル時点で未解決シンボル＝red。

    use umbra_eclipse::{EclipsePath, PathOptions};

    /// `PathArgs` を既定値（standard・limits あり・既定サンプル間隔）で組むヘルパ。
    /// 間隔の既定値は外部から持ち込まず `PathOptions::default()` から取る（確定仕様 §公開 IF）。
    fn path_args(date: &str, format: PathFormatArg) -> PathArgs {
        PathArgs {
            date: date.to_string(),
            accuracy: AccuracyArg::Standard,
            format,
            interval: PathOptions::default().sample_interval_seconds,
            no_limits: false,
        }
    }

    /// 参照 `EclipsePath` を公開 API だけで独立に組む（run_path の実装を写経しない独立オラクル）。
    /// 当日に日食が無ければ `None`。
    fn engine_path(date: &str, options: PathOptions) -> Option<EclipsePath> {
        let start = parse_date(date).expect("妥当日付");
        let end = UtcInstant::from_jd2(JulianDate2::from_jd(start.jd2().jd() + 1.0));
        let engine = standard_engine(bundled_time_data());
        let eclipse = engine
            .search(UtcRange { start, end })
            .expect("1 日範囲の探索は成功する")
            .into_iter()
            .next()?;
        Some(engine.path(&eclipse, options).expect("経路計算は成功する"))
    }

    /// `key` を含む **行** を返す（無ければ panic）。行ごと欠落する実装を明示的に落とすための補助。
    fn line_containing<'a>(out: &'a str, key: &str) -> &'a str {
        out.lines()
            .find(|l| l.contains(key))
            .unwrap_or_else(|| panic!("出力に '{key}' を含む行が無い: {out}"))
    }

    /// geojson 文字列をパースして `Value` を返す（パース成功＝妥当な JSON）。
    fn parse_geojson(s: &str) -> Value {
        serde_json::from_str(s).expect("geojson 出力は妥当な JSON（パース成功必須）")
    }

    // ------------------------------------------------------------------
    // 1. fast-fail（エンジン非実走＝FAST）
    // ------------------------------------------------------------------

    /// 【ISSUE-048 §確定仕様 7】不正日付は `search` を呼ぶ前に `Err(CliError::InvalidDate(入力))`。
    /// 殺す変異: 不正日付を黙って既定日に落として search へ進む、InvalidDate に入力文字列を載せない。
    #[test]
    fn run_path_invalid_date_is_invalid_date_error() {
        let args = path_args("not-a-date", PathFormatArg::Text);
        let r = run_path(&args);
        assert!(
            matches!(r, Err(CliError::InvalidDate(ref s)) if s == "not-a-date"),
            "expected Err(InvalidDate(\"not-a-date\")), got {r:?}"
        );
    }

    /// 【ISSUE-048 §確定仕様 4】`--interval 0` は `path` を呼ぶ前に
    /// `Err(CliError::InvalidInterval(0.0))`（**違反入力そのもの**を載せる）。
    /// 殺す変異: 0 を既定値 60 s に黙って差し替える、InvalidInterval に 0 でない値を載せる、
    ///   検査を engine 実走の後に置く。
    #[test]
    fn run_path_zero_interval_is_invalid_interval() {
        let mut args = path_args("2024-04-08", PathFormatArg::Text);
        args.interval = 0.0;
        let r = run_path(&args);
        match r {
            Err(CliError::InvalidInterval(v)) => {
                assert_eq!(v, 0.0, "InvalidInterval は違反入力 0.0 を保持");
            }
            other => panic!("expected Err(InvalidInterval(0.0)), got {other:?}"),
        }
    }

    /// 【ISSUE-048 §確定仕様 4】負の `--interval` も `InvalidInterval`（違反入力を保持）。
    /// **`--format geojson` でも同じ**（fast-fail は出力形式に依存しない）。
    /// 殺す変異: 負値を絶対値に補正する、geojson 経路だけ検査を飛ばす、別 variant を返す。
    #[test]
    fn run_path_negative_interval_is_invalid_interval_even_for_geojson() {
        let mut args = path_args("2024-04-08", PathFormatArg::Geojson);
        args.interval = -60.0;
        let r = run_path(&args);
        match r {
            Err(CliError::InvalidInterval(v)) => {
                assert_eq!(v, -60.0, "InvalidInterval は違反入力 -60.0 を保持");
            }
            other => panic!("expected Err(InvalidInterval(-60.0)), got {other:?}"),
        }
    }

    /// 【ISSUE-048 §確定仕様 4】`NaN` の `--interval` は非有限として `InvalidInterval`
    /// （載る値も NaN）。`v > 0.0` だけの素朴な検査は NaN を弾くが、`!(v <= 0.0)` 型の実装は
    /// 通してしまうため、NaN を独立に縛る。
    /// 殺す変異: `!(interval <= 0.0)` で判定して NaN を受理する、NaN を既定値に差し替える。
    #[test]
    fn run_path_nan_interval_is_invalid_interval() {
        let mut args = path_args("2024-04-08", PathFormatArg::Text);
        args.interval = f64::NAN;
        let r = run_path(&args);
        match r {
            Err(CliError::InvalidInterval(v)) => {
                assert!(
                    v.is_nan(),
                    "InvalidInterval は違反入力 NaN を保持（got {v}）"
                );
            }
            other => panic!("expected Err(InvalidInterval(NaN)), got {other:?}"),
        }
    }

    /// 【ISSUE-048 §確定仕様 4】`+∞` / `-∞` の `--interval` は非有限として `InvalidInterval`。
    /// 正の無限大は「非正」検査だけでは通ってしまうため、有限性検査の存在を独立に縛る。
    /// 殺す変異: `interval <= 0.0` のみ検査して +∞ を受理する（エンジンが無限ループ/発散する）。
    #[test]
    fn run_path_infinite_interval_is_invalid_interval() {
        for bad in [f64::INFINITY, f64::NEG_INFINITY] {
            let mut args = path_args("2024-04-08", PathFormatArg::Text);
            args.interval = bad;
            let r = run_path(&args);
            match r {
                Err(CliError::InvalidInterval(v)) => {
                    assert_eq!(v, bad, "InvalidInterval は違反入力 {bad} を保持");
                }
                other => panic!("expected Err(InvalidInterval({bad})), got {other:?}"),
            }
        }
    }

    // ------------------------------------------------------------------
    // 2. 日食の無い日（成功・SLOW＝1 日範囲の実 search）
    // ------------------------------------------------------------------

    /// 【SLOW】【ISSUE-048 §確定仕様 2】日食の無い日（2024-06-15）の text は **成功**し、
    /// `No solar eclipse on 2024-06-15.` の告知 1 行。架空の安定キー（'#'）を捏造しない。
    /// 日付の選定はテスト内でエンジン自身に確認させる（外部表を使わない）。
    /// 殺す変異: 日食なしで Err を返す（撤回された NoEclipse を含む）、空文字を返す、
    ///   架空イベントを 1 件でっち上げる、日付を告知文に載せない。
    // SLOW
    #[test]
    fn run_path_no_eclipse_date_text_announces_success() {
        // エンジン自身による検証: この日に日食は無い。
        assert!(
            engine_path("2024-06-15", PathOptions::default()).is_none(),
            "2024-06-15 はエンジン判定で日食なし（オラクルの自己検証）"
        );

        let out = run_path(&path_args("2024-06-15", PathFormatArg::Text))
            .expect("日食なし日も Ok（エラーにしない）");
        assert!(
            out.contains("No solar eclipse on 2024-06-15."),
            "告知 1 行（日付入り）が出力に含まれる: {out}"
        );
        assert!(
            !out.contains('#'),
            "架空 event_key（'#'）を捏造しない: {out}"
        );
    }

    /// 【SLOW】【ISSUE-048 §確定仕様 2】日食の無い日の geojson は **空 `FeatureCollection`**
    /// （`null` でも空文字でもない）。パース可能・`type == "FeatureCollection"`・`features == []`。
    /// 殺す変異: `null` を返す、空文字/空 JSON `{}` を返す、Err にする、
    ///   feature を 0 件でなく最大食点 1 件だけ捏造する。
    // SLOW
    #[test]
    fn run_path_no_eclipse_date_geojson_is_empty_feature_collection() {
        let out = run_path(&path_args("2024-06-15", PathFormatArg::Geojson))
            .expect("日食なし日も Ok（エラーにしない）");
        let v = parse_geojson(&out);
        assert_eq!(
            v["type"],
            Value::from("FeatureCollection"),
            "type は FeatureCollection: {out}"
        );
        let features = v["features"]
            .as_array()
            .unwrap_or_else(|| panic!("features は配列: {out}"));
        assert!(
            features.is_empty(),
            "日食なし日の features は空（0 件）: {out}"
        );
    }

    // ------------------------------------------------------------------
    // 3. geojson 出力（role 集合 == EclipsePath の Some 要素）
    // ------------------------------------------------------------------

    /// 【SLOW】【ISSUE-048 §確定仕様 6・§受け入れテスト戦略】`--format geojson` の出力は妥当な JSON
    /// かつ `FeatureCollection` で、`properties.role` の集合が **同一 options で独立に計算した
    /// `EclipsePath` の `Some` 要素と厳密一致**する（`greatest` は常に存在、`center_line` /
    /// `northern_limit` / `southern_limit` / `partial_limit` は Some のときだけ）。
    /// オラクルは外部表でなく同じ公開 API（`engine.path`）から導出する。
    /// 殺す変異: `to_geojson` を呼ばず自前で再整形する、None の要素にも feature を出す（座標捏造）、
    ///   Some の要素の feature を落とす、role 名を綴り違いにする、FeatureCollection でなく生配列を出す。
    // SLOW
    #[test]
    fn run_path_geojson_roles_match_some_elements_of_path() {
        let args = path_args("2024-04-08", PathFormatArg::Geojson);
        let options = PathOptions {
            sample_interval_seconds: args.interval,
            include_limits: true,
            ..PathOptions::default()
        };
        let reference = engine_path("2024-04-08", options).expect("2024-04-08 には日食がある");

        // 期待 role 集合（独立オラクル）。
        let mut expected: Vec<&str> = vec!["greatest"];
        if reference.center_line.is_some() {
            expected.push("center_line");
        }
        if reference.northern_limit.is_some() {
            expected.push("northern_limit");
        }
        if reference.southern_limit.is_some() {
            expected.push("southern_limit");
        }
        if reference.partial_limit.is_some() {
            expected.push("partial_limit");
        }
        expected.sort_unstable();

        let out = run_path(&args).expect("2024-04-08 の経路計算は成功する");
        let v = parse_geojson(&out);
        assert_eq!(
            v["type"],
            Value::from("FeatureCollection"),
            "type は FeatureCollection: {out}"
        );
        let features = v["features"]
            .as_array()
            .unwrap_or_else(|| panic!("features は配列: {out}"));
        let mut actual: Vec<String> = features
            .iter()
            .map(|f| {
                f["properties"]["role"]
                    .as_str()
                    .unwrap_or_else(|| panic!("各 feature は properties.role（文字列）を持つ: {f}"))
                    .to_string()
            })
            .collect();
        actual.sort();
        assert_eq!(
            actual,
            expected
                .iter()
                .map(|s| (*s).to_string())
                .collect::<Vec<_>>(),
            "role 集合は EclipsePath の Some 要素と一致する: {out}"
        );
    }

    // ------------------------------------------------------------------
    // 4. text 出力（None を「なし」と明示・--no-limits の lockstep）
    // ------------------------------------------------------------------

    /// 【SLOW】【ISSUE-048 §確定仕様 4/5・M9.7 lockstep】`--no-limits`（`include_limits=false`）では
    /// 北限・南限・部分食域が **行ごと残ったまま「なし」** になり、`samples` 件数が 0 になる。
    /// 中心線は残る（`no_limits` が中心線まで消さないことを同時に縛る）。
    /// 独立オラクル: 同じ options の `engine.path` が実際に limits=None・samples 空であることを
    /// テスト内で確認してから text を縛る。
    /// 殺す変異: `no_limits` を `PathOptions` に配線しない（限界線が出てしまう）、
    ///   「なし」の行を丸ごと省略する、samples 件数を limits 無効時も非 0 で出す、
    ///   `no_limits` で中心線まで落とす。
    // SLOW
    #[test]
    fn run_path_no_limits_reports_none_and_empty_samples() {
        let mut args = path_args("2024-04-08", PathFormatArg::Text);
        args.no_limits = true;
        let options = PathOptions {
            sample_interval_seconds: args.interval,
            include_limits: false,
            ..PathOptions::default()
        };
        let reference = engine_path("2024-04-08", options).expect("2024-04-08 には日食がある");
        assert!(
            reference.northern_limit.is_none() && reference.southern_limit.is_none(),
            "オラクル自己検証: include_limits=false で限界線は None"
        );
        assert!(
            reference.partial_limit.is_none(),
            "オラクル自己検証: include_limits=false で部分食域は None"
        );
        assert!(
            reference.samples.is_empty(),
            "オラクル自己検証: include_limits=false で samples は空（M9.7 lockstep）"
        );
        assert!(
            reference.center_line.is_some(),
            "オラクル自己検証: 中心食なので中心線は残る"
        );

        let out = run_path(&args).expect("--no-limits でも経路計算は成功する");
        for key in ["北限", "南限", "部分食域"] {
            let line = line_containing(&out, key);
            assert!(
                line.contains("なし"),
                "--no-limits では '{key}' の行が「なし」と明示される（空欄にしない）: {line}"
            );
        }
        let samples_line = line_containing(&out, "samples");
        assert!(
            samples_line.contains('0'),
            "samples 件数は 0 と明示される: {samples_line}"
        );
        let center_line = line_containing(&out, "中心線");
        assert!(
            !center_line.contains("なし"),
            "--no-limits でも中心線は残る（「なし」にならない）: {center_line}"
        );
    }

    /// 【SLOW】【ISSUE-048 §確定仕様 5/7・conventions §11】部分食（中心線を持たない日食）の text は
    /// `None` の要素を **行ごと残して「なし」と明示**する（行を落とす／空欄にするのを禁止）。
    /// 対象日 2025-03-29 は「中心線・北限・南限が None」の日として **エンジン自身で検証**してから
    /// 縛る（外部表は使わない）。万一この日が中心食だった場合はオラクル自己検証で落ちる。
    /// 殺す変異: None の要素の行を出力から省く（サイレントドロップ）、`Option` を空文字で描く、
    ///   None を 0 件の偽データ（空の LineString 等）として描く。
    // SLOW
    #[test]
    fn run_path_partial_eclipse_text_names_absent_elements_as_nashi() {
        let args = path_args("2025-03-29", PathFormatArg::Text);
        let options = PathOptions {
            sample_interval_seconds: args.interval,
            include_limits: true,
            ..PathOptions::default()
        };
        let reference = engine_path("2025-03-29", options).expect("2025-03-29 には日食がある");
        assert!(
            reference.center_line.is_none(),
            "オラクル自己検証: 2025-03-29 は部分食で中心線 None"
        );

        let out = run_path(&args).expect("部分食でも成功する（中心線なしはエラーでない）");
        // None の要素はすべて「なし」と明示される（行の存在と内容を分けて確認）。
        for (key, absent) in [
            ("中心線", reference.center_line.is_none()),
            ("北限", reference.northern_limit.is_none()),
            ("南限", reference.southern_limit.is_none()),
            ("部分食域", reference.partial_limit.is_none()),
        ] {
            let line = line_containing(&out, key);
            if absent {
                assert!(
                    line.contains("なし"),
                    "None の要素 '{key}' は「なし」と明示される: {line}"
                );
            } else {
                assert!(
                    !line.contains("なし"),
                    "Some の要素 '{key}' を「なし」と書かない: {line}"
                );
            }
        }
        // 空の成功出力を作らない（最低限 最大食時刻の日付は出る）。
        assert!(
            out.contains("2025-03-29"),
            "部分食でも最大食時刻（日付）が出る＝空の成功出力にしない: {out}"
        );
    }

    /// 【SLOW】【ISSUE-048 §確定仕様 5・§受け入れテスト戦略 SLOW】実日食 2024-04-08（皆既）の text は
    /// 種別（`Total`）・最大食時刻（当日日付）・**中心線の点数 > 0** を出す。点数は外部表でなく
    /// 同一 options の `engine.path` から導出した値と一致することで縛る（オラクルのハードコード禁止）。
    /// 殺す変異: 種別を出さない、最大食時刻を出さない、中心線点数を常に 0 と書く、
    ///   `--interval` を `PathOptions` に配線せず既定値で計算する（点数が参照と食い違う）。
    // SLOW
    #[test]
    fn run_path_2024_total_text_shows_kind_greatest_and_center_points() {
        let args = path_args("2024-04-08", PathFormatArg::Text);
        let options = PathOptions {
            sample_interval_seconds: args.interval,
            include_limits: true,
            ..PathOptions::default()
        };
        let reference = engine_path("2024-04-08", options).expect("2024-04-08 には日食がある");
        let n = reference
            .center_line
            .as_ref()
            .expect("2024-04-08 は中心食＝中心線あり")
            .points
            .len();
        assert!(n > 0, "オラクル自己検証: 中心線の点数 > 0（got {n}）");

        let out = run_path(&args).expect("2024-04-08 の経路計算は成功する");
        assert!(out.contains("Total"), "種別 Total が出力に含まれる: {out}");
        assert!(
            out.contains("2024-04-08"),
            "最大食時刻（当日日付）が出力に含まれる: {out}"
        );
        let center_line = line_containing(&out, "中心線");
        assert!(
            center_line.contains(&n.to_string()),
            "中心線の点数（参照計算 {n} 点・>0）が行に出る: {center_line}"
        );
    }

    /// 【SLOW】【ISSUE-048 §確定仕様 4・`PathOptions` 構築】`--interval` が **実際に**
    /// `PathOptions::sample_interval_seconds` として engine へ渡ることを、出力そのもので縛る。
    /// 同一日食（2024-04-08）を 2 つの間隔で解き、**細かい方が中心線の点数が厳密に多い**ことを
    /// assert する。点数の絶対値は外部表から取らず、2 回の実行を相互比較するだけ（自己参照）。
    ///
    /// 殺す実装: 「要求された `--interval` を捨てて既定値を黙って使う」
    /// （`PathOptions { .. }` から `sample_interval_seconds` の指定が落ち、構造体更新構文で
    /// `PathOptions::default()` の 60 s に戻る）。
    ///
    /// **間隔の選び方（この選択でなければ変異を殺せない理由）**: 2 値とも
    /// `PathOptions::default().sample_interval_seconds` と**異なる**値（既定の 2 倍と 6 倍）を選ぶ。
    /// こうすると、この変異下では 2 回の実行がどちらも既定 60 s で走り、
    /// 同一入力・同一 options ＝**点数が等しく**なるため、厳密不等号 `fine > coarse` が破れて
    /// テストが落ちる。正しい実装では 120 s と 360 s で走り、走査点数が 3 倍差になるので通る。
    /// （どちらか一方を既定値に選ぶと、変異下でも「既定 vs 既定」でなく偶然差が出るか等号になるかが
    /// 入力依存になり、殺せる保証が無い。）
    /// 実エンジンを 2 回実走するため SLOW（2 回に限定する）。
    // SLOW
    #[test]
    fn run_path_interval_is_wired_into_path_options_finer_gives_more_points() {
        /// text 出力の `  中心線: {n} 点 ...` 行から点数 `n` を取り出す（列位置に依存しない）。
        fn center_point_count(out: &str) -> usize {
            let line = out
                .lines()
                .find(|l| l.contains("中心線"))
                .unwrap_or_else(|| panic!("出力に「中心線」行が無い: {out}"));
            let rest = line
                .split_once("中心線:")
                .unwrap_or_else(|| panic!("「中心線:」の区切りが無い: {line}"))
                .1;
            let token = rest
                .split_whitespace()
                .next()
                .unwrap_or_else(|| panic!("「中心線:」の後に値が無い: {line}"));
            token.parse::<usize>().unwrap_or_else(|e| {
                panic!("中心線の点数をパースできない（{token:?}: {e}）: {line}")
            })
        }

        let default_interval = PathOptions::default().sample_interval_seconds;
        let fine = default_interval * 2.0;
        let coarse = default_interval * 6.0;
        assert!(
            fine != default_interval && coarse != default_interval && fine < coarse,
            "オラクル自己検証: 2 値とも既定 {default_interval} s と異なり fine < coarse \
             （fine={fine}, coarse={coarse}）"
        );

        let mut args = path_args("2024-04-08", PathFormatArg::Text);
        args.interval = fine;
        let fine_out = run_path(&args).expect("2024-04-08・細かい interval の経路計算は成功する");
        args.interval = coarse;
        let coarse_out = run_path(&args).expect("2024-04-08・粗い interval の経路計算は成功する");

        let fine_n = center_point_count(&fine_out);
        let coarse_n = center_point_count(&coarse_out);
        assert!(
            coarse_n > 0,
            "オラクル自己検証: 粗い方も中心線は存在する（got {coarse_n}）"
        );
        assert!(
            fine_n > coarse_n,
            "細かい interval（{fine} s）の中心線点数 {fine_n} は粗い interval（{coarse} s）の \
             {coarse_n} より厳密に多い＝--interval が PathOptions へ配線されている\n\
             fine:\n{fine_out}\ncoarse:\n{coarse_out}"
        );
    }

    // ------------------------------------------------------------------
    // 5. サンプル数上限（ISSUE-048 §確定仕様 4 追補・ハング防止）
    // ------------------------------------------------------------------
    // 正の有限値でも極端に小さい `--interval` は `path` の走査回数を事実上無限にし、CLI が無反応
    // のままハングする。日食特定後・`path` 呼び出し前に推定サンプル数 `span(P1..P4)/interval` が
    // `MAX_PATH_SAMPLES = 100_000` を超えたら `CliError::IntervalTooSmall` を返す契約を縛る。
    //
    // ## オラクル戦略
    // 閾値近傍の interval は **テスト内で P1〜P4 の span を公開 API（`eclipse.global`）から計算**して
    // 導出する（100_000 という定数以外の数値をハードコードしない）。

    /// 本テスト節が縛る上限定数（確定仕様 §4 追補）。
    use umbra_eclipse::MAX_PATH_SAMPLES as MAX_PATH_SAMPLES_ORACLE;

    /// 当日最初の日食を公開 API で取得する（上限検査の span オラクル用）。
    fn day_first_eclipse(date: &str) -> Option<SolarEclipse> {
        let start = parse_date(date).expect("妥当日付");
        let end = UtcInstant::from_jd2(JulianDate2::from_jd(start.jd2().jd() + 1.0));
        standard_engine(bundled_time_data())
            .search(UtcRange { start, end })
            .expect("1 日範囲の探索は成功する")
            .into_iter()
            .next()
    }

    /// P1〜P4（`partial_begin`/`partial_end`）の秒数（実装と同じ公開フィールドから導出）。
    fn partial_span_seconds(eclipse: &SolarEclipse) -> f64 {
        let p1 = eclipse
            .global
            .partial_begin
            .as_ref()
            .expect("実日食は P1 を持つ");
        let p4 = eclipse
            .global
            .partial_end
            .as_ref()
            .expect("実日食は P4 を持つ");
        (p4.time_utc.jd2().jd() - p1.time_utc.jd2().jd()) * 86_400.0
    }

    /// 【SLOW】【ISSUE-048 §確定仕様 4 追補】極端に小さい `--interval`（`1e-6` s）は、日食特定後・
    /// `path` 呼び出し**前**に `Err(CliError::IntervalTooSmall { interval, estimated_samples })`。
    /// `interval` には違反入力がそのまま載り、`estimated_samples` は有限かつ上限 100_000 を超える。
    /// 本テストが（分オーダーでなく）通常の探索時間で終わること自体が「`path` に入っていない」証拠
    /// （1e-6 s 間隔で `path` を実走すれば数十億点＝事実上ハングし、テストは終わらない）。
    /// 殺す変異: 上限検査を入れない（ハング）、検査を `path` の後に置く（ハング）、
    ///   `InvalidInterval` など別 variant を返す、`interval`/`estimated_samples` に別値を載せる。
    // SLOW
    #[test]
    fn run_path_tiny_interval_is_interval_too_small() {
        let mut args = path_args("2024-04-08", PathFormatArg::Text);
        args.interval = 1e-6;
        let r = run_path(&args);
        match r {
            Err(CliError::IntervalTooSmall {
                interval,
                estimated_samples,
            }) => {
                assert_eq!(interval, 1e-6, "IntervalTooSmall は違反入力 1e-6 を保持");
                assert!(
                    estimated_samples.is_finite(),
                    "estimated_samples は有限（got {estimated_samples}）"
                );
                assert!(
                    estimated_samples > MAX_PATH_SAMPLES_ORACLE,
                    "estimated_samples は上限 100_000 を超える（got {estimated_samples}）"
                );
            }
            other => panic!("expected Err(IntervalTooSmall {{..}}), got {other:?}"),
        }
    }

    /// 【SLOW】【ISSUE-048 §確定仕様 4 追補】上限は**通常利用を制限しない**: 閾値
    /// `span(P1..P4) / 100_000` から導出した合法な interval（閾値の 100 倍・既定 60 s より細かい）で
    /// `Ok` を返す。閾値そのものはテスト内で公開 API（`eclipse.global` の P1/P4）から計算し、
    /// 数値をハードコードしない。
    /// **なぜ閾値ぎりぎり（×1.01）でないか**: 閾値直上は約 100_000 点の経路追跡を実走することになり、
    /// テスト時間が実用外になる。ここでは「閾値より十分小さい interval でも弾かれない」ことを、
    /// 実行可能な点数（約 1,000 点）で縛る。
    /// 殺す変異: 上限検査を過剰に厳しくする（不等号の向き/等号ずれ、閾値を 1/100 等に取り違える）、
    ///   span を P1..P4 でなく別区間（中心食区間・1 日）で取って推定を過大にする。
    // SLOW
    #[test]
    fn run_path_interval_above_guard_threshold_succeeds() {
        let eclipse = day_first_eclipse("2024-04-08").expect("2024-04-08 には日食がある");
        let span = partial_span_seconds(&eclipse);
        assert!(
            span > 0.0,
            "オラクル自己検証: P1..P4 の span > 0（got {span}）"
        );
        let threshold = span / MAX_PATH_SAMPLES_ORACLE;

        let mut args = path_args("2024-04-08", PathFormatArg::Text);
        // 閾値の 100 倍＝推定 1,000 点（合法・実行可能）。
        args.interval = threshold * 100.0;
        assert!(
            args.interval < PathOptions::default().sample_interval_seconds,
            "オラクル自己検証: 既定 60 s より細かい interval を試している（got {})",
            args.interval
        );
        let out = run_path(&args).expect("閾値より粗い interval は上限に弾かれず成功する");
        assert!(
            !out.is_empty(),
            "成功出力は空でない（空の成功出力を作らない, conventions §11）"
        );
    }

    /// 【SLOW】【ISSUE-048 §確定仕様 4 追補】**既定 interval では上限が発火しない**
    /// （`PathOptions::default().sample_interval_seconds` で `Ok`）。実日食の推定点数が
    /// 上限を大きく下回ることをテスト内で確認してから縛る。
    /// 殺す変異: 上限を既定利用にかかる値（例 100）に取り違える、比較の向きを反転して常に弾く。
    // SLOW
    #[test]
    fn run_path_default_interval_is_not_rejected_by_guard() {
        let default_interval = PathOptions::default().sample_interval_seconds;
        let eclipse = day_first_eclipse("2024-04-08").expect("2024-04-08 には日食がある");
        let estimated = partial_span_seconds(&eclipse) / default_interval;
        assert!(
            estimated < MAX_PATH_SAMPLES_ORACLE,
            "オラクル自己検証: 既定 interval の推定点数は上限を下回る（got {estimated}）"
        );

        let args = path_args("2024-04-08", PathFormatArg::Text);
        assert_eq!(
            args.interval, default_interval,
            "ヘルパは既定 interval を使っている"
        );
        run_path(&args).expect("既定 interval は上限に弾かれない");
    }

    /// 【SLOW】【ISSUE-048 §確定仕様 2 と 4 追補の順序契約】上限検査は **日食を特定した後**に行う
    /// ため、日食の無い日（2024-06-15）＋極小 interval は `IntervalTooSmall` ではなく
    /// **日食なしの成功**を返す（非正・非有限 interval が `search` 前に弾かれるのとは対照的）。
    /// 殺す変異: 上限検査を `search` より前に置いて span を捏造する（日食が無いのに区間を仮定する）、
    ///   日食なし日で Err を返す、P1/P4 が無い場合に 0 や 1 日を span として使う。
    // SLOW
    #[test]
    fn run_path_no_eclipse_date_with_tiny_interval_still_reports_no_eclipse() {
        let mut args = path_args("2024-06-15", PathFormatArg::Text);
        args.interval = 1e-6;
        let out =
            run_path(&args).expect("日食なし日は極小 interval でも成功（上限検査は日食特定後）");
        assert!(
            out.contains("No solar eclipse on 2024-06-15."),
            "日食なしの告知 1 行が返る: {out}"
        );
    }

    // ------------------------------------------------------------------
    // 6. check_path_sample_count 単体（FAST・境界と分岐を厳密に縛る）
    // ------------------------------------------------------------------
    // 上限検査が独立関数になったので、CLI 表面からは到達できない境界・分岐を直接縛る。
    // 数値は `MAX_PATH_SAMPLES` から導出し、除算が**二進で厳密**になる組（interval は 2 の冪）を
    // 選ぶ。これにより「丸め次第でどちらにも転ぶ」曖昧さ無しに `>` と `>=` を判別できる。

    /// 【ISSUE-048 §確定仕様 4 追補】**境界は strict `>`**: 推定サンプル数がちょうど
    /// `MAX_PATH_SAMPLES` に等しいときは **受理**（`Ok`）。`interval = 1.0`（2 の冪）・
    /// `span = MAX_PATH_SAMPLES * 1.0` なので `span / interval` は二進で厳密に
    /// `MAX_PATH_SAMPLES` と一致し、丸めに依存しない。
    /// 殺す変異: 比較を `>=` にする（境界ちょうどを弾く）、上限に ±1 のオフセットを入れる。
    #[test]
    fn check_path_sample_count_accepts_exact_boundary() {
        let interval = 1.0;
        let span = MAX_PATH_SAMPLES * interval;
        assert_eq!(
            span / interval,
            MAX_PATH_SAMPLES,
            "オラクル自己検証: 除算が二進で厳密（境界ちょうど）"
        );
        let r = check_path_sample_count(Some(span), interval);
        assert!(
            r.is_ok(),
            "推定サンプル数が上限ちょうど（{MAX_PATH_SAMPLES}）なら受理する（strict >）: {r:?}"
        );
    }

    /// 【ISSUE-048 §確定仕様 4 追補】境界の**直上**は拒否: (a) 1 サンプル分だけ長い span、
    /// (b) 境界 span の次の表現可能値（1 ULP 上）。いずれも `IntervalTooSmall`。
    /// (b) により「上限に緩衝を足して実質 `>=` にする」実装も落ちる。
    /// 殺す変異: 比較を `>=`/`>` 以外（`>= MAX + 1` 等）に緩める、判定を span/interval でなく
    ///   別式にする、境界直上を受理する。
    #[test]
    fn check_path_sample_count_rejects_just_above_boundary() {
        let interval = 1.0;
        let boundary_span = MAX_PATH_SAMPLES * interval;

        // (a) ちょうど 1 サンプル分だけ超過。
        let span_a = boundary_span + interval;
        assert!(
            matches!(
                check_path_sample_count(Some(span_a), interval),
                Err(CliError::IntervalTooSmall { .. })
            ),
            "上限 +1 サンプルは IntervalTooSmall"
        );

        // (b) 境界 span の 1 ULP 上（最小の超過）。
        let span_b = f64::from_bits(boundary_span.to_bits() + 1);
        assert!(
            span_b / interval > MAX_PATH_SAMPLES,
            "オラクル自己検証: 1 ULP 上は上限を厳密に超える"
        );
        assert!(
            matches!(
                check_path_sample_count(Some(span_b), interval),
                Err(CliError::IntervalTooSmall { .. })
            ),
            "上限の 1 ULP 上も IntervalTooSmall（緩衝付き実装を撃破）"
        );
    }

    /// 【ISSUE-048 §確定仕様 4 追補】`span_seconds = None`（P1/P4 が無い）のときは **検査しない**＝
    /// `Ok(())`。極小 interval（`1e-300`）でも span を捏造してエラーにしない。
    /// CLI 表面（`run_path`）からは到達できない分岐をここで直接縛る。
    /// 殺す変異: `None` を 0 や 1 日（86_400 s）などの既定 span に差し替えて判定する、
    ///   `None` で無条件に Err にする、`unwrap()` で panic する。
    #[test]
    fn check_path_sample_count_none_span_is_ok() {
        let r = check_path_sample_count(None, 1e-300);
        assert!(
            r.is_ok(),
            "span 不明（None）なら捏造せず検査もしない（Ok）: {r:?}"
        );
    }

    /// 【ISSUE-048 §確定仕様 4 追補】拒否時のペイロード 2 フィールドを同時に縛る:
    /// `interval` は**渡した値そのもの**、`estimated_samples` は **`span / interval`**
    /// （上限値でも span でもない）。除算が二進で厳密になる組
    /// （`span = 1024.0`, `interval = 2^-10`）を使い、期待値をテスト内で独立計算する。
    /// 殺す変異: `estimated_samples` に `MAX_PATH_SAMPLES` や span をそのまま載せる、
    ///   2 フィールドを入れ替える、interval を正規化した値に差し替える。
    #[test]
    fn check_path_sample_count_error_payload_is_interval_and_estimate() {
        let span = 1024.0_f64;
        let interval = 1.0 / 1024.0; // 2^-10（二進で厳密）。
        let expected = span / interval; // 1_048_576（> MAX_PATH_SAMPLES）。
        assert!(
            expected > MAX_PATH_SAMPLES,
            "オラクル自己検証: この組は上限を超える（got {expected}）"
        );

        match check_path_sample_count(Some(span), interval) {
            Err(CliError::IntervalTooSmall {
                interval: got_interval,
                estimated_samples,
            }) => {
                assert_eq!(got_interval, interval, "interval は渡した値そのもの");
                assert_eq!(
                    estimated_samples, expected,
                    "estimated_samples は span / interval（上限値でも span でもない）"
                );
                assert_ne!(
                    estimated_samples, MAX_PATH_SAMPLES,
                    "estimated_samples に上限値を載せない"
                );
                assert_ne!(
                    estimated_samples, span,
                    "estimated_samples に span を載せない"
                );
            }
            other => panic!("expected Err(IntervalTooSmall {{..}}), got {other:?}"),
        }
    }

    /// 【SLOW】【ISSUE-048 §確定仕様 4 追補】`run_path` が上限検査へ渡す span が **P1〜P4**
    /// （`partial_begin`〜`partial_end`）であることを、返ってきた `estimated_samples` が
    /// テスト側で独立計算した `span(P1..P4) / interval` と**厳密一致**することで縛る。
    /// `estimated_samples > MAX` だけでは、U1〜U4・P1〜最大食・丸一日など「誤っているが十分大きい」
    /// span が素通りするため、span の出所そのものを固定する。
    /// 殺す変異: span を U1〜U4（central_begin/central_end）にする、P1〜最大食の片側だけにする、
    ///   1 日（86_400 s）や探索範囲を span に使う、秒換算係数（86_400）を取り違える。
    // SLOW
    #[test]
    fn run_path_guard_span_is_p1_to_p4() {
        let eclipse = day_first_eclipse("2024-04-08").expect("2024-04-08 には日食がある");
        let span = partial_span_seconds(&eclipse);

        let mut args = path_args("2024-04-08", PathFormatArg::Text);
        args.interval = 1e-6;
        let expected = span / args.interval;

        match run_path(&args) {
            Err(CliError::IntervalTooSmall {
                estimated_samples, ..
            }) => {
                assert_eq!(
                    estimated_samples, expected,
                    "estimated_samples は P1〜P4 span（{span} s）/ interval と厳密一致する"
                );
            }
            other => panic!("expected Err(IntervalTooSmall {{..}}), got {other:?}"),
        }
    }

    // ==================================================================
    // 6. text 整形の純関数単体（FAST・合成データのみ／ISSUE-048 §5）
    // ==================================================================
    // ## なぜ上の `run_path` テストと重複して見えるか（意図的）
    // 上節（区分「出力契約」）は同じ §5 契約を **実エンジン経由**で縛るため 1 件あたり数十秒かかり、
    // `docs/reviews/mutation-cli-path.md` のとおり text 整形系 32 変異の判別力が **未証明**のまま
    // 残った（完走に 3〜4 時間）。本節は `format_path_text` / `format_limit_line` を
    // **合成 `SolarEclipse` + `EclipsePath` で直接呼ぶ**（エンジン・探索・日付解決なし＝ミリ秒）。
    // 役割分担: 上節＝「実データでも契約が成り立つ」、本節＝「整形分岐の全経路を高速に殺す」。
    // どちらも消さない（上を消すと配線が、下を消すと分岐が無検証になる）。
    //
    // ## オラクル戦略
    // - 値はすべて **相互に異なる非対称値**（点数 5/7/3、環数 2・外環 9・穴 4、samples 6、
    //   lat≠lon・始点≠終点・符号違い）。取り違え（lat↔lon、始↔終、環数↔頂点数、北↔南）が
    //   必ず出力差になるよう選ぶ。
    // - `format_limit_line` は crate 内の私有関数だが、本テストモジュールは同一 crate 内なので
    //   `super::*` 経由で **直接呼べる**（境界 `n > 0` を単体で縛る）。

    use umbra_eclipse::PathSample;
    use umbra_geo::{GeoLine, GeoPolygon};

    /// 度の組から折れ線を組む。
    fn synth_line(points: &[(f64, f64)]) -> GeoLine {
        GeoLine::new(points.iter().map(|&(la, lo)| geo(la, lo)).collect())
    }

    /// 度の組の列から環（`Vec<GeoPoint>`）を組む。
    fn synth_ring(points: &[(f64, f64)]) -> Vec<GeoPoint> {
        points.iter().map(|&(la, lo)| geo(la, lo)).collect()
    }

    /// `n` 個のサンプル（内容は整形に使われないので最小・件数だけが契約）。
    fn synth_samples(n: usize) -> Vec<PathSample> {
        std::iter::repeat_with(|| PathSample {
            time_utc: utc(2024, 4, 8, 18, 0, 0.0),
            center: geo(1.5, -2.75),
            duration_seconds: 201.0,
            sun_altitude: Degrees(60.0),
            path_width: Kilometers(190.0),
            kind: SolarEclipseKind::Total,
        })
        .take(n)
        .collect()
    }

    /// 全要素 `None`・samples 空の経路（最大食点だけを持つ骨組み）。
    fn synth_empty_path() -> EclipsePath {
        EclipsePath {
            center_line: None,
            northern_limit: None,
            southern_limit: None,
            partial_limit: None,
            partial_limit_pole: None,
            // lat と lon は絶対値・符号ともに異なる（取り違えが必ず出力差になる）。
            greatest_point: geo(17.5, -66.25),
            samples: vec![],
        }
    }

    /// 全要素 `Some` かつ非空・各件数が相互に異なる経路。
    /// 中心線 5 点（始点≠終点）／北限 7 点／南限 3 点／環数 2・外環 9 頂点・穴 4 頂点／samples 6 件。
    fn synth_full_path() -> EclipsePath {
        EclipsePath {
            center_line: Some(synth_line(&[
                (11.25, -33.5),
                (12.5, -20.0),
                (13.75, -10.0),
                (-5.0, 60.0),
                (-44.75, 122.25),
            ])),
            northern_limit: Some(synth_line(&[
                (1.0, 2.0),
                (3.0, 4.0),
                (5.0, 6.0),
                (7.0, 8.0),
                (9.0, 10.0),
                (11.0, 12.0),
                (13.0, 14.0),
            ])),
            southern_limit: Some(synth_line(&[(-1.0, -2.0), (-3.0, -4.0), (-5.0, -6.0)])),
            partial_limit: Some(GeoPolygon::new(vec![
                // 外環 9 頂点（環数 2 とも穴 4 とも異なる＝取り違え検出可能）。
                synth_ring(&[
                    (0.0, 0.0),
                    (10.0, 0.0),
                    (20.0, 5.0),
                    (30.0, 10.0),
                    (35.0, 20.0),
                    (30.0, 30.0),
                    (20.0, 35.0),
                    (10.0, 30.0),
                    (0.0, 20.0),
                ]),
                // 穴 4 頂点。
                synth_ring(&[(15.0, 15.0), (16.0, 15.0), (16.0, 16.0), (15.0, 16.0)]),
            ])),
            partial_limit_pole: None,
            greatest_point: geo(17.5, -66.25),
            samples: synth_samples(6),
        }
    }

    /// §5 のラベル表の全ラベルが **過不足なく 1 回ずつ・表の順序で** 出る
    /// （種別 → 最大食 → 中心線 → 北限 → 南限 → 部分食域 → samples）。
    /// 殺す変異: いずれかの行を落とす（サイレントドロップ）、行を二重に出す、順序を入れ替える、
    ///   ラベル文字列を別語に変える。
    #[test]
    fn format_path_text_prints_all_spec_labels_once_in_order() {
        let out = format_path_text(&total_eclipse(), &synth_full_path());
        let labels = [
            "種別",
            "最大食",
            "中心線",
            "北限",
            "南限",
            "部分食域",
            "samples",
        ];
        let mut previous = 0usize;
        for label in labels {
            assert_eq!(
                out.matches(label).count(),
                1,
                "ラベル '{label}' はちょうど 1 回だけ出る: {out}"
            );
            let at = out.find(label).expect("ラベルが存在する");
            assert!(
                at > previous,
                "ラベル '{label}' は §5 表の順序（前のラベルより後ろ）に出る: {out}"
            );
            previous = at;
        }
    }

    /// §5「`None` の要素は行を省略せず『なし』」— `center_line` のみ `None`。
    /// 他要素は `Some` のまま残すので、「全部なしにする」実装では落ちる。
    /// 殺す変異: 中心線が `None` のとき行ごと省く、空文字にする、0 点と偽って書く、
    ///   `None` 判定を他フィールドへ流用する。
    #[test]
    fn format_path_text_center_line_none_reads_nashi() {
        let mut path = synth_full_path();
        path.center_line = None;
        let out = format_path_text(&total_eclipse(), &path);
        assert_eq!(
            line_containing(&out, "中心線"),
            "  中心線: なし",
            "center_line=None は「中心線: なし」の 1 行: {out}"
        );
        assert!(
            !line_containing(&out, "北限").contains("なし"),
            "他要素（北限）は Some のまま件数を出す: {out}"
        );
    }

    /// §5「`None` の要素は行を省略せず『なし』」— `northern_limit` のみ `None`。
    /// 南限を非空のまま残すので、北限と南限を取り違える実装も落ちる。
    /// 殺す変異: 北限行の省略、北限に南限の点数を書く、`None` を 0 点として書く。
    #[test]
    fn format_path_text_northern_limit_none_reads_nashi() {
        let mut path = synth_full_path();
        path.northern_limit = None;
        let out = format_path_text(&total_eclipse(), &path);
        assert_eq!(
            line_containing(&out, "北限"),
            "  北限: なし",
            "northern_limit=None は「北限: なし」: {out}"
        );
        assert_eq!(
            line_containing(&out, "南限"),
            "  南限: 3 点",
            "南限は Some のまま 3 点（北↔南の取り違え検出）: {out}"
        );
    }

    /// §5「`None` の要素は行を省略せず『なし』」— `southern_limit` のみ `None`。
    /// 北限を非空のまま残すので、南限行に北限の値を流用する実装も落ちる。
    /// 殺す変異: 南限行の省略、南限に北限の点数を書く、`None` を 0 点として書く。
    #[test]
    fn format_path_text_southern_limit_none_reads_nashi() {
        let mut path = synth_full_path();
        path.southern_limit = None;
        let out = format_path_text(&total_eclipse(), &path);
        assert_eq!(
            line_containing(&out, "南限"),
            "  南限: なし",
            "southern_limit=None は「南限: なし」: {out}"
        );
        assert_eq!(
            line_containing(&out, "北限"),
            "  北限: 7 点",
            "北限は Some のまま 7 点（南↔北の取り違え検出）: {out}"
        );
    }

    /// §5「`None` の要素は行を省略せず『なし』」— `partial_limit` のみ `None`。
    /// 殺す変異: 部分食域行の省略、`None` を「0 環 外環 0 頂点」と偽る、`Option` 分岐の反転。
    #[test]
    fn format_path_text_partial_limit_none_reads_nashi() {
        let mut path = synth_full_path();
        path.partial_limit = None;
        let out = format_path_text(&total_eclipse(), &path);
        assert_eq!(
            line_containing(&out, "部分食域"),
            "  部分食域: なし",
            "partial_limit=None は「部分食域: なし」: {out}"
        );
        assert!(
            line_containing(&out, "中心線").contains("5 点"),
            "他要素（中心線）は Some のまま: {out}"
        );
    }

    /// §5「中心線＝点数と始点・終点／北限・南限＝点数／部分食域＝環数と外環頂点数」を
    /// **相互に異なる件数**（5 / 7 / 3 / 2 環・外環 9）で同時に縛る。
    /// 殺す変異: `len()` を `len()±1` にする、北限と南限の値を入れ替える、環数と外環頂点数を
    ///   入れ替える、穴の頂点数（4）を外環として書く、いずれかを定数にする。
    #[test]
    fn format_path_text_reports_distinct_counts_for_each_element() {
        let out = format_path_text(&total_eclipse(), &synth_full_path());
        assert!(
            line_containing(&out, "中心線").starts_with("  中心線: 5 点"),
            "中心線は 5 点: {out}"
        );
        assert_eq!(
            line_containing(&out, "北限"),
            "  北限: 7 点",
            "北限は 7 点: {out}"
        );
        assert_eq!(
            line_containing(&out, "南限"),
            "  南限: 3 点",
            "南限は 3 点: {out}"
        );
        assert_eq!(
            line_containing(&out, "部分食域"),
            "  部分食域: 2 環  外環 9 頂点",
            "部分食域は 2 環・外環 9 頂点（穴の 4 ではない）: {out}"
        );
    }

    /// §5「中心線は始点・終点を出す」— 始点 `(11.25, -33.5)`・終点 `(-44.75, 122.25)` は
    /// 緯度経度も符号も互いに異なるので、**始↔終の転置**と **lat↔lon の転置**の双方が検出される。
    /// 殺す変異: `points[0]` と `points[len-1]` を入れ替える、`points[1]`/`points[len-2]` を使う、
    ///   lat と lon を入れ替えて書く、終点を始点で代用する。
    #[test]
    fn format_path_text_center_line_prints_first_then_last_point() {
        let out = format_path_text(&total_eclipse(), &synth_full_path());
        let line = line_containing(&out, "中心線");
        assert_eq!(
            line, "  中心線: 5 点  始 lat 11.2500° lon -33.5000°  終 lat -44.7500° lon 122.2500°",
            "始点＝先頭要素・終点＝末尾要素を lat, lon の順で出す: {out}"
        );
    }

    /// §5・conventions §11「空の成功出力を作らない」— `Some` でも **点列が空**の中心線は
    /// 「0 点」ではなく「なし」（現行の確定挙動の固定）。
    /// 殺す変異: 空点列で `points[0]` を触って panic する、`Some(空)` を「0 点」と書く、
    ///   空判定 `!is_empty()` を落とす。
    #[test]
    fn format_path_text_empty_center_line_reads_nashi() {
        let mut path = synth_full_path();
        path.center_line = Some(GeoLine::new(vec![]));
        let out = format_path_text(&total_eclipse(), &path);
        assert_eq!(
            line_containing(&out, "中心線"),
            "  中心線: なし",
            "空点列の Some も「なし」: {out}"
        );
    }

    /// §5・conventions §11 — `Some` でも **点列が空**の北限・南限は「0 点」ではなく「なし」。
    /// 北限のみ空・南限は非空にして、片側だけ壊す変異も捕まえる。
    /// 殺す変異: 空点列を「0 点」と書く、`n > 0` の境界を `n >= 0` にする、空判定を落とす。
    #[test]
    fn format_path_text_empty_limit_lines_read_nashi() {
        let mut path = synth_full_path();
        path.northern_limit = Some(GeoLine::new(vec![]));
        let out = format_path_text(&total_eclipse(), &path);
        assert_eq!(
            line_containing(&out, "北限"),
            "  北限: なし",
            "空点列の北限は「なし」: {out}"
        );
        assert_eq!(
            line_containing(&out, "南限"),
            "  南限: 3 点",
            "南限は非空のまま 3 点: {out}"
        );

        let mut path = synth_full_path();
        path.southern_limit = Some(GeoLine::new(vec![]));
        let out = format_path_text(&total_eclipse(), &path);
        assert_eq!(
            line_containing(&out, "南限"),
            "  南限: なし",
            "空点列の南限は「なし」: {out}"
        );
        assert_eq!(
            line_containing(&out, "北限"),
            "  北限: 7 点",
            "北限は非空のまま 7 点: {out}"
        );
    }

    /// §5・conventions §11 — `Some` でも **環が 1 本も無い**部分食域は「0 環」ではなく「なし」。
    /// 殺す変異: 空 `rings` で `rings[0]` に触って panic する、「0 環  外環 0 頂点」と書く、
    ///   `!rings.is_empty()` の判定を落とす／反転する。
    #[test]
    fn format_path_text_empty_partial_limit_polygon_reads_nashi() {
        let mut path = synth_full_path();
        path.partial_limit = Some(GeoPolygon::new(vec![]));
        let out = format_path_text(&total_eclipse(), &path);
        assert_eq!(
            line_containing(&out, "部分食域"),
            "  部分食域: なし",
            "環の無い Some も「なし」: {out}"
        );
    }

    /// §5「部分食域＝環数と外環頂点数」— 穴つき多角形（環数 2・外環 9 頂点・穴 4 頂点）で、
    /// **報告されるのは外環（`rings[0]`）の頂点数**であることを縛る。3 数が相互に異なるので
    /// どの取り違えも検出できる。
    /// 殺す変異: `rings[1]`（穴）の頂点数を外環として書く、環数と頂点数を入れ替える、
    ///   全環の頂点合計（13）を書く、環数を定数 1 にする。
    #[test]
    fn format_path_text_partial_limit_reports_rings_and_outer_ring_vertices() {
        let out = format_path_text(&total_eclipse(), &synth_full_path());
        let line = line_containing(&out, "部分食域");
        assert_eq!(
            line, "  部分食域: 2 環  外環 9 頂点",
            "環数 2・外環頂点数 9（穴の 4 でも合計 13 でもない）: {out}"
        );
    }

    /// §5「`samples` は `None` ではなく 0 件で表現」— 空と 6 件の双方を縛る。
    /// 殺す変異: samples 行を落とす、件数を定数にする、0 件を「なし」と書く（§5 表は `0`）、
    ///   `len()±1`。
    #[test]
    fn format_path_text_samples_count_includes_zero() {
        let out = format_path_text(&total_eclipse(), &synth_empty_path());
        assert_eq!(
            line_containing(&out, "samples"),
            "  samples: 0",
            "空 samples は 0 件と書く（「なし」でも行省略でもない）: {out}"
        );

        let out = format_path_text(&total_eclipse(), &synth_full_path());
        assert_eq!(
            line_containing(&out, "samples"),
            "  samples: 6",
            "6 件の samples は 6 と書く: {out}"
        );
    }

    /// §5「種別」「最大食＝UTC 時刻と地点」— fixture の最大食 UTC は 2024-04-08 18:17:00.0、
    /// TT は同日 12:00 + 0.123 日（≒14:57）なので **UTC/TT 取り違えは時刻差として現れる**。
    /// 地点は lat 17.5 / lon -66.25 と符号も絶対値も異なるので **lat↔lon 転置も検出**される。
    /// 座標は `SolarEclipse` 側の `greatest.position` ではなく `EclipsePath::greatest_point` から
    /// 取る契約（fixture 側は lat 25 / lon -104 と別値にしてある）。
    /// 殺す変異: `time_tt` を印字する、地点の lat と lon を入れ替える、
    ///   `eclipse.global.greatest.position` を使う、種別行を落とす／固定値にする。
    #[test]
    fn format_path_text_prints_kind_and_greatest_time_and_position() {
        let out = format_path_text(&total_eclipse(), &synth_full_path());
        assert_eq!(
            line_containing(&out, "種別"),
            "  種別: Total",
            "種別は SolarEclipseKind の表示: {out}"
        );
        let line = line_containing(&out, "最大食");
        assert!(
            line.contains("2024-04-08 18:17:00.0 UTC"),
            "最大食は UTC 時刻（TT の 14:57 ではない）: {out}"
        );
        assert!(
            line.contains("lat 17.5000°") && line.contains("lon -66.2500°"),
            "最大食地点は EclipsePath::greatest_point の lat/lon（転置なし）: {out}"
        );
        assert!(
            !line.contains("25.0000") && !line.contains("-104.0000"),
            "SolarEclipse 側の greatest.position（25, -104）ではない: {out}"
        );

        // 種別は fixture 依存（固定文字列でない）ことを別 fixture で確認。
        let out = format_path_text(&partial_eclipse(), &synth_empty_path());
        assert_eq!(
            line_containing(&out, "種別"),
            "  種別: Partial",
            "部分食 fixture では Partial: {out}"
        );
        assert!(
            line_containing(&out, "最大食").contains("2025-03-29 10:47:00.0 UTC"),
            "最大食時刻も fixture 依存: {out}"
        );
    }

    /// §5 限界線行の整形単体（私有関数 `format_limit_line` を同一 crate 内から直接呼ぶ）。
    /// `None` と `Some(0)` はいずれも「なし」、`Some(n>0)` は「n 点」。境界は `n > 0`。
    /// 殺す変異: `n > 0` を `n >= 0` にする（`Some(0)` が「0 点」になる）、`None` 分岐を消す、
    ///   ラベルを引数でなく定数にする、改行を落とす。
    #[test]
    fn format_limit_line_none_and_zero_are_nashi_positive_is_count() {
        assert_eq!(format_limit_line("北限", None), "  北限: なし\n");
        assert_eq!(format_limit_line("南限", Some(0)), "  南限: なし\n");
        assert_eq!(format_limit_line("北限", Some(1)), "  北限: 1 点\n");
        assert_eq!(format_limit_line("南限", Some(7)), "  南限: 7 点\n");
        // ラベルは引数どおり（定数化・入れ替えを撃破）。
        assert_eq!(format_limit_line("南限", Some(3)), "  南限: 3 点\n");
    }
}
