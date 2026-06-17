//! TOML フィクスチャのデシリアライズと公開型への変換（ISSUE-029 §アルゴリズム概要）。
//!
//! データは `data/golden/*.toml`（手書き Rust 配列直書き禁止・architecture §11）。serde の生 DTO で
//! 読み、境界で度・時刻文字列を公開型（`UtcInstant`/`TtInstant`・enum）へ変換する。
//! バンドル済みの自データを対象とするため、フォーマット不正は**即パニック**で露見させる
//! （テスト専用 crate。データ作成ミスを沈黙させない）。

use serde::Deserialize;

use umbra_core::{gregorian_to_jd2, TtInstant, UtcInstant};
use umbra_eclipse::{SolarEclipseKind, Visibility};

use crate::types::{GoldenContact, GoldenEclipse, GoldenLocation, LocationClass, OracleSource};
use crate::FIXTURE_FILES;

/// 現在同梱されている全ゴールデン日食を返す（固定回帰の正本・ISSUE-029）。
///
/// v0.1 最終目標は 20 日食。本 seed では実在オラクルから完全転記した少数を返す。
pub fn golden_eclipses() -> Vec<GoldenEclipse> {
    FIXTURE_FILES
        .iter()
        .map(|toml_src| {
            let raw: RawEclipse =
                toml::from_str(toml_src).expect("bundled golden fixture must be valid TOML");
            convert(raw)
        })
        .collect()
}

// ------------------------------------------------------------------
// 生 DTO（TOML スキーマ）
// ------------------------------------------------------------------

#[derive(Deserialize)]
struct RawEclipse {
    event_key: String,
    kind: String,
    greatest_utc: String,
    greatest_tt: Option<String>,
    gamma: f64,
    magnitude: f64,
    delta_t_seconds: Option<f64>,
    source: RawSource,
    #[serde(default)]
    location: Vec<RawLocation>,
}

#[derive(Deserialize)]
struct RawSource {
    name: String,
    url: String,
    retrieved: String,
    delta_t_convention: String,
    k_convention: String,
    license_note: String,
}

#[derive(Deserialize)]
struct RawLocation {
    name: String,
    latitude_deg: f64,
    east_longitude_deg: f64,
    elevation_m: f64,
    class: String,
    visibility: String,
    magnitude: f64,
    obscuration: f64,
    max_altitude_deg: f64,
    max_azimuth_deg: f64,
    maximum: RawContact,
    c1: Option<RawContact>,
    c2: Option<RawContact>,
    c3: Option<RawContact>,
    c4: Option<RawContact>,
}

#[derive(Deserialize)]
struct RawContact {
    /// `YYYY-MM-DD HH:MM:SS.s`（UTC）。
    utc: String,
    /// 太陽高度（度）。
    altitude_deg: f64,
    /// TT(=TD) 文字列（USNO は UT のみのため通常省略）。
    tt: Option<String>,
}

// ------------------------------------------------------------------
// 変換（生 DTO → 公開型）
// ------------------------------------------------------------------

fn convert(raw: RawEclipse) -> GoldenEclipse {
    GoldenEclipse {
        event_key: raw.event_key,
        kind_expected: parse_kind(&raw.kind),
        greatest_time_utc: parse_utc(&raw.greatest_utc),
        greatest_time_tt: raw.greatest_tt.as_deref().map(parse_tt),
        gamma: raw.gamma,
        magnitude: raw.magnitude,
        delta_t_seconds: raw.delta_t_seconds,
        locations: raw.location.into_iter().map(convert_location).collect(),
        source: OracleSource {
            name: raw.source.name,
            url: raw.source.url,
            retrieved: raw.source.retrieved,
            delta_t_convention: raw.source.delta_t_convention,
            k_convention: raw.source.k_convention,
            license_note: raw.source.license_note,
        },
    }
}

fn convert_location(raw: RawLocation) -> GoldenLocation {
    GoldenLocation {
        name: raw.name,
        latitude_deg: raw.latitude_deg,
        east_longitude_deg: raw.east_longitude_deg,
        elevation_m: raw.elevation_m,
        location_class: parse_class(&raw.class),
        c1: raw.c1.as_ref().map(convert_contact),
        c2: raw.c2.as_ref().map(convert_contact),
        maximum: convert_contact(&raw.maximum),
        c3: raw.c3.as_ref().map(convert_contact),
        c4: raw.c4.as_ref().map(convert_contact),
        magnitude: raw.magnitude,
        obscuration: raw.obscuration,
        max_altitude_deg: raw.max_altitude_deg,
        max_azimuth_deg: raw.max_azimuth_deg,
        visibility_expected: parse_visibility(&raw.visibility),
    }
}

fn convert_contact(raw: &RawContact) -> GoldenContact {
    GoldenContact {
        time_utc: parse_utc(&raw.utc),
        time_tt: raw.tt.as_deref().map(parse_tt),
        altitude_deg: raw.altitude_deg,
    }
}

// ------------------------------------------------------------------
// 文字列パーサ（バンドル自データ前提・不正は panic）
// ------------------------------------------------------------------

/// `YYYY-MM-DD HH:MM:SS.s` を (年, 月, 日, 時, 分, 秒) に分解する。
fn parse_datetime(s: &str) -> (i32, u8, u8, u8, u8, f64) {
    let (date, time) = s
        .split_once(' ')
        .unwrap_or_else(|| panic!("fixture datetime '{s}' must be 'YYYY-MM-DD HH:MM:SS'"));
    let mut d = date.split('-');
    let year: i32 = next_field(&mut d, s, "year");
    let month: u8 = next_field(&mut d, s, "month");
    let day: u8 = next_field(&mut d, s, "day");
    let mut t = time.split(':');
    let hour: u8 = next_field(&mut t, s, "hour");
    let minute: u8 = next_field(&mut t, s, "minute");
    let second: f64 = next_field(&mut t, s, "second");
    (year, month, day, hour, minute, second)
}

fn next_field<'a, T, I>(it: &mut I, whole: &str, field: &str) -> T
where
    T: std::str::FromStr,
    I: Iterator<Item = &'a str>,
{
    it.next()
        .unwrap_or_else(|| panic!("fixture datetime '{whole}' missing {field}"))
        .parse()
        .unwrap_or_else(|_| panic!("fixture datetime '{whole}' has invalid {field}"))
}

fn parse_utc(s: &str) -> UtcInstant {
    let (y, mo, d, h, mi, sec) = parse_datetime(s);
    UtcInstant::from_gregorian(y, mo, d, h, mi, sec)
        .unwrap_or_else(|_| panic!("fixture UTC '{s}' is not a valid calendar instant"))
}

/// 入力文字列は**既に TT(=TD) スケール**の暦表現であることを前提とする（オラクル NASA が TD で
/// 与える値をそのまま取り込む。UTC↔TT のオフセットは適用しない・accuracy.md §0）。USNO 局地の
/// UT 値に対しては呼ばない（TT を捏造しない＝接触の `time_tt` は `None`）。
fn parse_tt(s: &str) -> TtInstant {
    let (y, mo, d, h, mi, sec) = parse_datetime(s);
    let jd = gregorian_to_jd2(y, mo, d, h, mi, sec)
        .unwrap_or_else(|_| panic!("fixture TT '{s}' is not a valid calendar instant"));
    TtInstant::from_jd2(jd)
}

fn parse_kind(s: &str) -> SolarEclipseKind {
    match s {
        "Total" => SolarEclipseKind::Total,
        "Annular" => SolarEclipseKind::Annular,
        "Partial" => SolarEclipseKind::Partial,
        "Hybrid" => SolarEclipseKind::Hybrid,
        "NonCentralTotal" => SolarEclipseKind::NonCentralTotal,
        "NonCentralAnnular" => SolarEclipseKind::NonCentralAnnular,
        other => panic!("unknown SolarEclipseKind '{other}' in fixture"),
    }
}

fn parse_class(s: &str) -> LocationClass {
    match s {
        "Centerline" => LocationClass::Centerline,
        "NearLimit" => LocationClass::NearLimit,
        "PartialZone" => LocationClass::PartialZone,
        "Sunrise" => LocationClass::Sunrise,
        "Sunset" => LocationClass::Sunset,
        "BelowHorizon" => LocationClass::BelowHorizon,
        "HighElevation" => LocationClass::HighElevation,
        other => panic!("unknown LocationClass '{other}' in fixture"),
    }
}

fn parse_visibility(s: &str) -> Visibility {
    match s {
        "NotVisible" => Visibility::NotVisible,
        "BelowHorizon" => Visibility::BelowHorizon,
        "SunriseEclipse" => Visibility::SunriseEclipse,
        "SunsetEclipse" => Visibility::SunsetEclipse,
        "PartialVisible" => Visibility::PartialVisible,
        "FullyVisible" => Visibility::FullyVisible,
        other => panic!("unknown Visibility '{other}' in fixture"),
    }
}
