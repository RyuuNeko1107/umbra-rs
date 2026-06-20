//! `xtask differential` — DE 差分・誤差層分解の実行とレポート出力（ISSUE-030・accuracy.md §4.1）。
//!
//! `umbra-fixtures` の層分解ハーネス（[`report_differential`]）を、**解析暦エンジン**（[`EngineGoldenComputer`]）
//! と **JPL DE Reference エンジン**（[`JplGoldenComputer`]・DE440s）の 2 系統で同梱ゴールデンに適用し、
//! 各 metric の誤差を暦層（analytical−DE）/ 幾何・数値層（DE−golden）/ 総合（analytical−golden）へ
//! 分解する（accuracy.md §4.1 の 2 バケット層分解）。`umbra-fixtures`・DE は検証専用のため本番 `umbra`
//! CLI には含めず、dev ツールである xtask に結線する（architecture §1）。
//!
//! 純粋寄りの結線（[`differential_report`]・引数解釈）は高速ユニットで検証し、実 DE エンジン経路
//! （[`JplGoldenComputer`]）は SLOW 統合で担保する（負荷配分）。DE440s SPK は非同梱・実行時に解決する。

use std::path::Path;

use umbra_core::{EspenakMeeusDeltaT, IersEopData, JulianDate2, Observer, UtcInstant};
use umbra_eclipse::{
    EclipseEngine, EclipseError, EngineConfig, LocalCircumstances, LunarRadiusModel, SolarEclipse,
    UtcRange,
};
use umbra_ephemeris::{bundled_time_data, JplEphemeris};
use umbra_fixtures::{
    golden_eclipses, render_differential_json, render_differential_text, report_differential,
    GoldenComputer, GoldenEclipse, GoldenLocation,
};

use crate::error::XtaskError;
use crate::validate::{
    parse_accuracy, parse_format, AccuracyArg, EngineGoldenComputer, ReportFormat,
};

/// 既定の DE440s SPK パス（リポジトリ root の data/spk）。CARGO_MANIFEST_DIR は crates/xtask。
const DE440S_DEFAULT_PATH: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../data/spk/de440s.bsp");

/// DE Reference 版エンジン型（JPL DE440s バックエンド・解析暦版 [`StandardEngine`] と暦のみ差し替え）。
///
/// [`crate::validate::EngineGoldenComputer`] が包む `StandardEngine`
/// （= `EclipseEngine<AnalyticalEphemeris, …>`）の暦バックエンドだけを [`JplEphemeris`] に替えたもの。
/// ΔT・EOP・幾何/接触パイプラインは共通＝[`report_differential`] が暦層を清浄に切り出せる（§4.1）。
type ReferenceDeEngine = EclipseEngine<JplEphemeris, EspenakMeeusDeltaT, IersEopData>;

/// JPL DE440s Reference 版 [`GoldenComputer`]。ゴールデンの最大食 UTC を中心に ±0.5 日窓を `search` し
/// 最初の食を採用、地点は [`Observer`] を構築して `local_circumstances` を解く（解析暦版と同一挙動）。
#[derive(Debug)]
pub struct JplGoldenComputer {
    engine: ReferenceDeEngine,
}

impl JplGoldenComputer {
    /// SPK（de440s.bsp）パスと精度プロファイルから DE エンジンを構築する。
    ///
    /// SPK の IO/形式エラー（[`JplEphemeris::from_spk_path`] 由来）は descriptive な
    /// [`XtaskError::MalformedSource`] へ写像する（不在/不正パスは `Err`）。月半径 k 慣習は
    /// 解析暦版（[`EngineGoldenComputer`]）と同じく `UMBRA_VALIDATE_K` で切替可能（オラクル整合の裏取り）。
    pub fn from_spk_path(spk_path: &Path, accuracy: AccuracyArg) -> Result<Self, XtaskError> {
        let jpl = JplEphemeris::from_spk_path(spk_path).map_err(|e| {
            XtaskError::MalformedSource(format!("DE440s SPK {}: {e:?}", spk_path.display()))
        })?;
        let time = bundled_time_data();
        let mut config = match accuracy {
            AccuracyArg::Standard => EngineConfig::standard(),
            AccuracyArg::Reference => EngineConfig::reference(),
        };
        // 解析暦版と同じ k 慣習切替（NASA/USNO は Espenak 2 値・エンジン既定は IauMean）。
        if let Ok(k) = std::env::var("UMBRA_VALIDATE_K") {
            config.lunar_radius_model = match k.as_str() {
                "espenak-umbral" => LunarRadiusModel::EspenakUmbral,
                "espenak-penumbral" => LunarRadiusModel::EspenakPenumbral,
                "iau-mean" => LunarRadiusModel::IauMean,
                other => {
                    eprintln!("[differential] unknown UMBRA_VALIDATE_K={other}; keeping default");
                    config.lunar_radius_model
                }
            };
            eprintln!(
                "[differential] DE lunar_radius_model = {}",
                config.lunar_radius_model.name()
            );
        }
        let earth_orientation = time.eop().clone();
        let engine = EclipseEngine::new(jpl, EspenakMeeusDeltaT, earth_orientation, time, config);
        Ok(Self { engine })
    }
}

impl GoldenComputer for JplGoldenComputer {
    fn eclipse_on(&self, golden: &GoldenEclipse) -> Result<Option<SolarEclipse>, EclipseError> {
        // ゴールデン最大食 UTC を中心に ±0.5 日窓で探索（当日の食を確実に括る）。
        let center_jd = golden.greatest_time_utc.jd2().jd();
        let start = UtcInstant::from_jd2(JulianDate2::from_jd(center_jd - 0.5));
        let end = UtcInstant::from_jd2(JulianDate2::from_jd(center_jd + 0.5));
        let t0 = std::time::Instant::now();
        eprintln!("[differential] DE search {} ...", golden.event_key);
        let found = self
            .engine
            .search(UtcRange { start, end })?
            .into_iter()
            .next();
        // 実 DE は重いので進捗診断（解析暦版 EngineGoldenComputer と対称）。全球最大食時刻誤差
        // （DE − golden, 秒）も per-eclipse 出力（異常食の特定用）。
        let greatest_err_s = found.as_ref().map(|e| {
            (e.global.greatest.time_utc.jd2().jd() - golden.greatest_time_utc.jd2().jd()) * 86_400.0
        });
        eprintln!(
            "[differential] DE search {} done in {:.1}s (found={}) greatest_err={:?}s",
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
        self.engine.local_circumstances(eclipse, observer)
    }
}

/// 注入された解析暦・DE の 2 computer でゴールデンを層分解し、`format` に応じ整形済み文字列を返す。
///
/// [`report_differential`]（暦層/幾何層/総合の 3 層分解）→ [`render_differential_text`] /
/// [`render_differential_json`] の結線。エンジン非依存（computer 注入）で、`report_differential` の
/// [`EclipseError`] は [`XtaskError::Eclipse`]、JSON 整形失敗は [`XtaskError::Json`] へ透過する。
pub fn differential_report<A, D>(
    analytical: &A,
    de: &D,
    golden: &[GoldenEclipse],
    format: ReportFormat,
) -> Result<String, XtaskError>
where
    A: GoldenComputer,
    D: GoldenComputer,
{
    let report = report_differential(analytical, de, golden)?;
    match format {
        ReportFormat::Text => Ok(render_differential_text(&report)),
        ReportFormat::Json => Ok(render_differential_json(&report)?),
    }
}

/// `differential` サブコマンド実体。`--accuracy`/`--format` を解釈し、解析暦（[`EngineGoldenComputer`]）と
/// DE（[`JplGoldenComputer`]・既定 [`DE440S_DEFAULT_PATH`]）で同梱ゴールデン（[`golden_eclipses`]）を
/// 層分解し、レポートを標準出力へ印字する（**実 DE×解析暦の 2 エンジン実走＝非常に低速**）。
///
/// `UMBRA_VALIDATE_ONLY=key1,key2` で event_key 部分一致のみ照合（代表食での素早い確認用）。
pub fn run_differential(args: &[String]) -> Result<(), XtaskError> {
    let accuracy = parse_accuracy(args)?;
    let format = parse_format(args)?;
    let de = JplGoldenComputer::from_spk_path(Path::new(DE440S_DEFAULT_PATH), accuracy)?;
    let analytical = EngineGoldenComputer::new(accuracy);
    let mut golden = golden_eclipses();
    // 開発時の部分実行: UMBRA_VALIDATE_ONLY=key1,key2 で event_key 部分一致のみ照合（2 エンジン実走は
    // 非常に重く全件は長時間かかるため、代表食での素早い精度・層分解確認に使う。未設定なら全件）。
    if let Ok(filter) = std::env::var("UMBRA_VALIDATE_ONLY") {
        let keys: Vec<&str> = filter
            .split(',')
            .map(|s| s.trim())
            .filter(|s| !s.is_empty())
            .collect();
        golden.retain(|g| keys.iter().any(|k| g.event_key.contains(k)));
        eprintln!(
            "[differential] filtered to {} eclipse(s): {:?}",
            golden.len(),
            keys
        );
    }
    let output = differential_report(&analytical, &de, &golden, format)?;
    print!("{output}");
    Ok(())
}
