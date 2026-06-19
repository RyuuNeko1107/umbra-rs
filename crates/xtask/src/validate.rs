//! `xtask validate` — ゴールデン照合の実行とレポート出力（ISSUE-030 S30f）。
//!
//! `umbra-fixtures` の検証ハーネス（[`report_against_golden`]）を**実エンジン**で同梱ゴールデンに
//! 適用し、誤差統計レポートを人間可読テキスト／機械可読 JSON で出力する。`umbra-fixtures` は
//! 検証専用のため本番 `umbra` CLI には含めず、dev ツールである xtask に結線する（architecture §1）。
//!
//! 純粋寄りの結線（[`validate_report`]・引数解釈）は高速ユニットで検証し、実エンジン経路
//! ([`EngineGoldenComputer`]) は SLOW 統合で担保する（負荷配分）。

use umbra_core::{EspenakMeeusDeltaT, JulianDate2, Observer, UtcInstant};
use umbra_eclipse::{
    EclipseEngine, EclipseError, EngineConfig, LocalCircumstances, LunarRadiusModel, SolarEclipse,
    StandardEngine, UtcRange,
};
use umbra_ephemeris::{bundled_time_data, AnalyticalEphemeris};
use umbra_fixtures::{
    golden_eclipses, render_json, render_text, report_against_golden, GoldenComputer,
    GoldenEclipse, GoldenLocation, ToleranceProfile,
};

use crate::error::XtaskError;

/// 出力形式（`--format <text|json>`・既定 text）。
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ReportFormat {
    /// 人間可読テキストサマリ。
    Text,
    /// 機械可読 JSON（CI/履歴比較）。
    Json,
}

/// 精度プロファイル引数（`--accuracy <standard|reference>`・既定 standard）。
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum AccuracyArg {
    /// 本番標準。
    Standard,
    /// 高精度参照。
    Reference,
}

/// `args` から `flag` の値を取り出す（無指定→`None`、値欠落→[`XtaskError::MissingArgument`]）。
fn flag_value<'a>(args: &'a [String], flag: &str) -> Result<Option<&'a str>, XtaskError> {
    match args.iter().position(|arg| arg == flag) {
        None => Ok(None),
        Some(index) => args
            .get(index + 1)
            .map(|value| Some(value.as_str()))
            .ok_or_else(|| XtaskError::MissingArgument(flag.to_string())),
    }
}

/// `--format <text|json>` を解釈する（既定 [`ReportFormat::Text`]・未知値は [`XtaskError::InvalidArgument`]）。
pub fn parse_format(args: &[String]) -> Result<ReportFormat, XtaskError> {
    match flag_value(args, "--format")? {
        None | Some("text") => Ok(ReportFormat::Text),
        Some("json") => Ok(ReportFormat::Json),
        Some(other) => Err(XtaskError::InvalidArgument {
            flag: "--format".to_string(),
            value: other.to_string(),
        }),
    }
}

/// `--accuracy <standard|reference>` を解釈する（既定 [`AccuracyArg::Standard`]・未知値は
/// [`XtaskError::InvalidArgument`]）。
pub fn parse_accuracy(args: &[String]) -> Result<AccuracyArg, XtaskError> {
    match flag_value(args, "--accuracy")? {
        None | Some("standard") => Ok(AccuracyArg::Standard),
        Some("reference") => Ok(AccuracyArg::Reference),
        Some(other) => Err(XtaskError::InvalidArgument {
            flag: "--accuracy".to_string(),
            value: other.to_string(),
        }),
    }
}

/// `AccuracyArg` に対応する許容プロファイル。
fn tolerance_profile(accuracy: AccuracyArg) -> ToleranceProfile {
    match accuracy {
        AccuracyArg::Standard => ToleranceProfile::standard(),
        AccuracyArg::Reference => ToleranceProfile::reference(),
    }
}

/// 実エンジン版 [`GoldenComputer`]。ゴールデンの最大食 UTC を中心に ±0.5 日窓を `search` し最初の食を
/// 採用、地点は [`Observer`] を構築して `local_circumstances` を解く。
#[derive(Debug)]
pub struct EngineGoldenComputer {
    engine: StandardEngine,
}

impl EngineGoldenComputer {
    /// 精度プロファイルで同梱データ（[`bundled_time_data`]）からエンジンを構築する。
    pub fn new(accuracy: AccuracyArg) -> Self {
        let time = bundled_time_data();
        let mut config = match accuracy {
            AccuracyArg::Standard => EngineConfig::standard(),
            AccuracyArg::Reference => EngineConfig::reference(),
        };
        // 月半径 k 慣習の切替（オラクル整合の裏取り用）。NASA/USNO は Espenak 2値 k を採用、
        // エンジン既定は IauMean。`UMBRA_VALIDATE_K=espenak-umbral|espenak-penumbral|iau-mean`。
        if let Ok(k) = std::env::var("UMBRA_VALIDATE_K") {
            config.lunar_radius_model = match k.as_str() {
                "espenak-umbral" => LunarRadiusModel::EspenakUmbral,
                "espenak-penumbral" => LunarRadiusModel::EspenakPenumbral,
                "iau-mean" => LunarRadiusModel::IauMean,
                other => {
                    eprintln!("[validate] unknown UMBRA_VALIDATE_K={other}; keeping default");
                    config.lunar_radius_model
                }
            };
            eprintln!(
                "[validate] lunar_radius_model = {}",
                config.lunar_radius_model.name()
            );
        }
        let earth_orientation = time.eop().clone();
        let engine = EclipseEngine::new(
            AnalyticalEphemeris::new(),
            EspenakMeeusDeltaT,
            earth_orientation,
            time,
            config,
        );
        Self { engine }
    }
}

impl GoldenComputer for EngineGoldenComputer {
    fn eclipse_on(&self, golden: &GoldenEclipse) -> Result<Option<SolarEclipse>, EclipseError> {
        // ゴールデン最大食 UTC を中心に ±0.5 日窓で探索（当日の食を確実に括る）。
        let center_jd = golden.greatest_time_utc.jd2().jd();
        let start = UtcInstant::from_jd2(JulianDate2::from_jd(center_jd - 0.5));
        let end = UtcInstant::from_jd2(JulianDate2::from_jd(center_jd + 0.5));
        let t0 = std::time::Instant::now();
        eprintln!("[validate] search {} ...", golden.event_key);
        let found = self
            .engine
            .search(UtcRange { start, end })?
            .into_iter()
            .next();
        // 全球最大食時刻の誤差（engine − golden, 秒）を per-eclipse 出力（異常食の特定用）。
        let greatest_err_s = found.as_ref().map(|e| {
            (e.global.greatest.time_utc.jd2().jd() - golden.greatest_time_utc.jd2().jd()) * 86_400.0
        });
        eprintln!(
            "[validate] search {} done in {:.1}s (found={}) greatest_err={:?}s",
            golden.event_key,
            t0.elapsed().as_secs_f64(),
            found.is_some(),
            greatest_err_s.map(|s| (s * 10.0).round() / 10.0)
        );
        Ok(found)
    }

    fn local_at(
        &self,
        eclipse: &SolarEclipse,
        location: &GoldenLocation,
    ) -> Result<LocalCircumstances, EclipseError> {
        // 西経入力は Observer::from_degrees が東経へ吸収。DomainError は EclipseError へ透過。
        let observer = Observer::from_degrees(
            location.latitude_deg,
            location.east_longitude_deg,
            location.elevation_m,
        )?;
        let t0 = std::time::Instant::now();
        eprintln!("[validate]   local {} ...", location.name);
        let result = self.engine.local_circumstances(eclipse, observer);
        // per-location 診断: 最大食時刻誤差(engine − golden UTC, 秒)と可視性一致。
        if let Ok(local) = result.as_ref() {
            let max_err = (local.contacts.maximum.time_utc.jd2().jd()
                - location.maximum.time_utc.jd2().jd())
                * 86_400.0;
            let vis_ok = local.visibility == location.visibility_expected;
            eprintln!(
                "[validate]   local {} done in {:.1}s max_err={:.1}s vis={}{}",
                location.name,
                t0.elapsed().as_secs_f64(),
                max_err,
                if vis_ok { "OK" } else { "MISMATCH" },
                if vis_ok {
                    String::new()
                } else {
                    format!(
                        " (eng={:?} golden={:?})",
                        local.visibility, location.visibility_expected
                    )
                },
            );
        } else {
            eprintln!(
                "[validate]   local {} done in {:.1}s (err)",
                location.name,
                t0.elapsed().as_secs_f64()
            );
        }
        result
    }
}

/// computer でゴールデン照合し、`format` に応じて整形済みレポート文字列を返す（エンジン非依存・
/// computer 注入）。[`report_against_golden`] → [`render_text`]/[`render_json`] の結線。
pub fn validate_report<C: GoldenComputer>(
    computer: &C,
    golden: &[GoldenEclipse],
    profile: &ToleranceProfile,
    format: ReportFormat,
) -> Result<String, XtaskError> {
    let report = report_against_golden(computer, golden, profile)?;
    match format {
        ReportFormat::Text => Ok(render_text(&report)),
        ReportFormat::Json => Ok(render_json(&report)?),
    }
}

/// `validate` サブコマンド実体。`--accuracy`/`--format` を解釈し、実エンジンで同梱ゴールデン
/// （[`golden_eclipses`]）を照合してレポートを標準出力へ印字する（実エンジン実走＝低速）。
pub fn run_validate(args: &[String]) -> Result<(), XtaskError> {
    let accuracy = parse_accuracy(args)?;
    let format = parse_format(args)?;
    let profile = tolerance_profile(accuracy);
    let computer = EngineGoldenComputer::new(accuracy);
    let mut golden = golden_eclipses();
    // 開発時の部分実行: UMBRA_VALIDATE_ONLY=key1,key2 で event_key 部分一致のみ照合（実エンジンが
    // 重く全20件は ~2h かかるため、代表食での素早い精度確認に使う。未設定なら全件）。
    if let Ok(filter) = std::env::var("UMBRA_VALIDATE_ONLY") {
        let keys: Vec<&str> = filter
            .split(',')
            .map(|s| s.trim())
            .filter(|s| !s.is_empty())
            .collect();
        golden.retain(|g| keys.iter().any(|k| g.event_key.contains(k)));
        eprintln!(
            "[validate] filtered to {} eclipse(s): {:?}",
            golden.len(),
            keys
        );
    }
    let output = validate_report(&computer, &golden, &profile, format)?;
    print!("{output}");
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::{tolerance_profile, AccuracyArg};
    use umbra_fixtures::ToleranceProfile;

    /// 非公開 `tolerance_profile` の写像: Standard→standard()・Reference→reference()。
    /// run_validate（実エンジン・SLOW のみ）経由でしか踏まれないため直接縛る。
    /// 殺す変異: standard↔reference 取り違え（誤った許容で合否判定が狂う）。
    #[test]
    fn tolerance_profile_maps_accuracy_to_profile() {
        assert_eq!(
            tolerance_profile(AccuracyArg::Standard),
            ToleranceProfile::standard(),
            "Standard は standard() 許容"
        );
        assert_eq!(
            tolerance_profile(AccuracyArg::Reference),
            ToleranceProfile::reference(),
            "Reference は reference() 許容"
        );
    }
}
