//! `xtask sweep` — 全食スイープの自己カタログ集計＋完備性突合（ISSUE-030・accuracy.md §3.4）。
//!
//! 解析暦エンジン（[`standard_engine`]）で指定範囲（既定 1900-2100）を `search` し、検出食群を
//! [`summarize_sweep`] でカタログ集計（種別件数・`|gamma|`/食分の範囲統計）し、NASA 4 区分の
//! **完備性突合**（偽陰性＝取りこぼし検知）を行う。完備性の期待件数（オラクル）は **CLI フラグで
//! 供給**する（`--expected-total/annular/hybrid/partial`）＝オラクル件数は利用者が出典付きで指定
//! （ハードコードしない・provenance 規律, data-sources §0/§4。参考: NASA 5MCSE 2001-2100=全 224 件）。
//! 既定（フラグ無し）は期待件数 0＝カタログ集計のみ（完備性は vacuous）。
//!
//! 純粋寄りの結線（[`parse_range`]/[`parse_expected_counts`]/[`sweep_report`]）は高速ユニットで
//! 検証し、実エンジン `search`（[`run_sweep`]・範囲が広いと**非常に低速**）は SLOW 統合で担保する。
//!
//! 注: search の偽陰性ゼロ・マージンの実余裕統計（D6・accuracy.md §3.4）は coarse-scan 内部の計装を
//! 要し本スコープ外（後続課題）。本ランナーは**集計レベルの完備性突合**で取りこぼしを検知する。

use umbra_core::UtcInstant;
use umbra_eclipse::{standard_engine, SolarEclipse, UtcRange};
use umbra_ephemeris::bundled_time_data;
use umbra_fixtures::{render_sweep_json, render_sweep_text, summarize_sweep, CatalogCounts};

use crate::error::XtaskError;
use crate::validate::{flag_value, parse_format, ReportFormat};

/// 既定の探索開始年（解析暦の対応域下端・accuracy.md §1）。
const DEFAULT_FROM_YEAR: i32 = 1900;
/// 既定の探索終了年（解析暦の対応域上端）。
const DEFAULT_TO_YEAR: i32 = 2100;

/// `flag` の整数値を取り出す（無指定→`default`・非整数→[`XtaskError::InvalidArgument`]）。
fn parse_int_flag(args: &[String], flag: &str, default: i32) -> Result<i32, XtaskError> {
    match flag_value(args, flag)? {
        None => Ok(default),
        Some(value) => value
            .parse::<i32>()
            .map_err(|_| XtaskError::InvalidArgument {
                flag: flag.to_string(),
                value: value.to_string(),
            }),
    }
}

/// `flag` の件数（非負整数）を取り出す（無指定→0・非整数→[`XtaskError::InvalidArgument`]）。
fn parse_count_flag(args: &[String], flag: &str) -> Result<usize, XtaskError> {
    match flag_value(args, flag)? {
        None => Ok(0),
        Some(value) => value
            .parse::<usize>()
            .map_err(|_| XtaskError::InvalidArgument {
                flag: flag.to_string(),
                value: value.to_string(),
            }),
    }
}

/// 指定年の 1/1 0h UTC を組む（暦変換失敗は [`XtaskError::InvalidArgument`]）。
fn year_start(year: i32, flag: &str) -> Result<UtcInstant, XtaskError> {
    UtcInstant::from_gregorian(year, 1, 1, 0, 0, 0.0).map_err(|e| XtaskError::InvalidArgument {
        flag: flag.to_string(),
        value: format!("year {year}: {e:?}"),
    })
}

/// `--from <year>` / `--to <year>` を解釈し `[start, end)` の UTC 対を返す（各年 1/1 0h UTC）。
///
/// 既定 from=1900・to=2100。非整数は [`XtaskError::InvalidArgument`]。`from > to`（範囲逆転）も
/// [`XtaskError::InvalidArgument`]（`from == to` は空区間として `Ok`）。
pub fn parse_range(args: &[String]) -> Result<(UtcInstant, UtcInstant), XtaskError> {
    let from_year = parse_int_flag(args, "--from", DEFAULT_FROM_YEAR)?;
    let to_year = parse_int_flag(args, "--to", DEFAULT_TO_YEAR)?;
    if from_year > to_year {
        return Err(XtaskError::InvalidArgument {
            flag: "--from".to_string(),
            value: format!("{from_year} > --to {to_year}（範囲逆転）"),
        });
    }
    Ok((
        year_start(from_year, "--from")?,
        year_start(to_year, "--to")?,
    ))
}

/// `--expected-total/annular/hybrid/partial <N>` を解釈し [`CatalogCounts`] を返す（既定すべて 0）。
///
/// 非整数は [`XtaskError::InvalidArgument`]。一部のみ指定なら他区分は 0（部分指定は独立）。
pub fn parse_expected_counts(args: &[String]) -> Result<CatalogCounts, XtaskError> {
    Ok(CatalogCounts {
        total: parse_count_flag(args, "--expected-total")?,
        annular: parse_count_flag(args, "--expected-annular")?,
        hybrid: parse_count_flag(args, "--expected-hybrid")?,
        partial: parse_count_flag(args, "--expected-partial")?,
    })
}

/// 検出食群を [`summarize_sweep`] し、`format` に応じ整形済み文字列を返す（エンジン非依存・純粋寄り）。
///
/// [`render_sweep_text`]（Text）/ [`render_sweep_json`]（Json・末尾改行）の結線。JSON 整形失敗は
/// [`XtaskError::Json`] へ透過する。
pub fn sweep_report(
    eclipses: &[SolarEclipse],
    expected: CatalogCounts,
    format: ReportFormat,
) -> Result<String, XtaskError> {
    let summary = summarize_sweep(eclipses, expected);
    match format {
        ReportFormat::Text => Ok(render_sweep_text(&summary)),
        ReportFormat::Json => Ok(render_sweep_json(&summary)?),
    }
}

/// `sweep` サブコマンド実体。`--from`/`--to`/`--format`/`--expected-*` を解釈し、解析暦エンジンで
/// 範囲を `search` → [`sweep_report`] → 標準出力へ印字する（**実エンジン実走＝範囲が広いと非常に低速**）。
pub fn run_sweep(args: &[String]) -> Result<(), XtaskError> {
    let format = parse_format(args)?;
    let (start, end) = parse_range(args)?;
    let expected = parse_expected_counts(args)?;
    let engine = standard_engine(bundled_time_data());
    eprintln!("[sweep] search {start:?} .. {end:?}（解析暦・実エンジン実走＝低速）...");
    let eclipses = engine.search(UtcRange { start, end })?;
    eprintln!("[sweep] detected {} eclipse(s)", eclipses.len());
    let output = sweep_report(&eclipses, expected, format)?;
    print!("{output}");
    Ok(())
}
