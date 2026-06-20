//! M9.3 経路生成スライス（南北限界線・geometric 近似）の統合テスト（確定 path）。
//!
//! `umbra-eclipse` の**公開 API のみ**を対象とした統合テスト（tests/ 配下・別クレート境界）。
//! 対象は `EclipseEngine::path(&SolarEclipse, PathOptions) -> Result<EclipsePath, EclipseError>`。
//! M9.1（中心線）・M9.2（GeoJSON）は実装済み。本スライスは中心食の本影帯の北/南縁
//! （`northern_limit` / `southern_limit`）の生成を縛る。
//!
//! ## 確定セマンティクス（geometric 近似・テストで縛る）
//! 1. `greatest_point`・`center_line`・`samples`・`partial_limit` は M9.1/現状どおり
//!    （center_line は中心食で Some・samples 空・partial_limit None）。
//! 2. **中心食（central_begin/central_end 両方 Some）かつ `include_limits == true`**:
//!    `northern_limit = Some(GeoLine)` かつ `southern_limit = Some(GeoLine)`。中心線と同じサンプル
//!    時刻列で、各時刻に本影帯の北縁/南縁の地表点を結ぶ。**高緯度側が北限・低緯度側が南限**。
//! 3. **`include_limits == false`**: 限界線は両方 None（center_line は include_limits に依らず
//!    中心食なら Some）。
//! 4. **非中心**（central_begin か central_end が None）: 限界線も center_line も None。
//! 5. 限界線は GeoJSON にはまだ出さない（本スライスでは to_geojson を変えない）。
//!
//! ## テスト戦略（strict / mutation-resistant / 負荷配分）
//! - FAST（実 search 非実走・合成 SolarEclipse）: path_center_line.rs の合成中心食ヘルパを流用。
//!   limits=Some・点列非空・各点妥当・北限緯度 ≥ 中心線緯度 ≥ 南限緯度（対応サンプルで）・
//!   北限と南限が分離（緯度差が正で過大でない）を縛る。include_limits=false / 非中心は None。
//! - SLOW（実エンジン 1 件・2024-04-08 皆既）: search → path()。限界線 Some・中心線が南北の間・
//!   最大食付近の帯幅が NASA 公表 ~197 km の妥当域（geometric 近似ゆえ広め 100–350 km）。
//!
//! ## 期待される RED（実装前）
//! 現状 path() は `northern_limit = None` / `southern_limit = None` を返すため、limits=Some を
//! 要求する FAST/SLOW assert が落ちる。コンパイルは通る（新シンボル無し）。

use umbra_core::{JulianDate2, TimeInterval, TtInstant, UtcInstant};
use umbra_eclipse::{
    standard_engine, AccuracyProfile, BesselFitError, BesselianPolynomial, CalculationMetadata,
    EclipseMagnitude, GlobalCircumstances, GlobalContact, GreatestEclipse, Obscuration,
    PathOptions, Polynomial, SolarEclipse, SolarEclipseKind,
};
use umbra_ephemeris::bundled_time_data;

// ============================================================
// 時刻 / 地理ヘルパ（path_center_line.rs ミラー）
// ============================================================

/// TT 時刻を 2 要素 JD から構築。
fn tt(jd1: f64, jd2: f64) -> TtInstant {
    TtInstant::from_jd2(JulianDate2::new(jd1, jd2))
}

/// epoch から経過時間 `hours`[hour] だけ進んだ TT 時刻。
fn tt_at_hours(epoch: TtInstant, hours: f64) -> TtInstant {
    TtInstant::from_jd2(epoch.jd2().add_days(hours / 24.0))
}

/// UTC 瞬時（合成日食の時刻ラベル用・幾何には無関係）。
fn utc(year: i32, month: u8, day: u8, hour: u8, minute: u8, second: f64) -> UtcInstant {
    UtcInstant::from_gregorian(year, month, day, hour, minute, second).expect("有効な UTC 日時")
}

/// 地表点（度）。
fn geo(lat: f64, lon: f64) -> umbra_geo::GeoPoint {
    umbra_geo::GeoPoint::from_degrees(lat, lon).expect("有効な地表点")
}

/// 緯度（度）取り出し。
fn lat_deg(p: &umbra_geo::GeoPoint) -> f64 {
    p.lat.degrees().0
}

/// 経度（度）取り出し。
fn lon_deg(p: &umbra_geo::GeoPoint) -> f64 {
    p.lon.degrees().0
}

/// 緯度・経度が妥当な範囲か（限界線各点のサニティ）。
fn lat_lon_in_range(p: &umbra_geo::GeoPoint) -> bool {
    let lat = lat_deg(p);
    let lon = lon_deg(p);
    (-90.0..=90.0).contains(&lat) && (-180.0..=180.0).contains(&lon)
}

/// 2 点間の概算大円距離 [km]（haversine・WGS84 平均半径 6371 km）。
fn great_circle_km(a: &umbra_geo::GeoPoint, b: &umbra_geo::GeoPoint) -> f64 {
    let r = 6371.0_f64;
    let (lat1, lon1) = (lat_deg(a).to_radians(), lon_deg(a).to_radians());
    let (lat2, lon2) = (lat_deg(b).to_radians(), lon_deg(b).to_radians());
    let dlat = lat2 - lat1;
    let dlon = lon2 - lon1;
    let h = (dlat / 2.0).sin().powi(2) + lat1.cos() * lat2.cos() * (dlon / 2.0).sin().powi(2);
    2.0 * r * h.sqrt().asin()
}

// ============================================================
// 合成 BesselianPolynomial / SolarEclipse（path_center_line.rs ミラー）
// ============================================================

/// 合成日食の epoch（窓の中心 TT）。J2000 を借用（解析暦の評価可能域内・幾何のみ使用）。
fn synth_epoch() -> TtInstant {
    tt(2_451_545.0, 0.0)
}

/// 中心食を模す合成 bessel poly（影軸が地表に当たる小さい x,y）。
/// path_center_line.rs の central_bessel と同値。
fn central_bessel() -> BesselianPolynomial {
    let epoch = synth_epoch();
    let c = |v: f64| Polynomial {
        coefficients: vec![v],
    };
    BesselianPolynomial {
        epoch_tt: epoch,
        x: Polynomial {
            coefficients: vec![0.0, 0.05],
        },
        y: c(0.10),
        d: c(0.20),
        mu: c(1.2),
        l1: c(0.54),
        l2: c(-0.009),
        tan_f1: 0.004_65,
        tan_f2: 0.004_63,
        fit_interval: TimeInterval {
            start: tt_at_hours(epoch, -2.0),
            end: tt_at_hours(epoch, 2.0),
        },
        fit_error: BesselFitError {
            max_x: 1.0e-7,
            max_y: 1.0e-7,
            max_l1: 1.0e-7,
            max_l2: 1.0e-7,
        },
    }
}

/// 代表メタデータ（幾何に無関係）。
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

/// U1/U4 の全球接触点（時刻のみ幾何に効く）。
fn contact(time_tt: TtInstant) -> GlobalContact {
    GlobalContact {
        time_utc: utc(2024, 4, 8, 18, 0, 0.0),
        time_tt,
        position: geo(20.0, -100.0),
    }
}

/// 中心食の最大食。
fn greatest_at(position: umbra_geo::GeoPoint) -> GreatestEclipse {
    GreatestEclipse {
        time_utc: utc(2024, 4, 8, 18, 17, 0.0),
        time_tt: synth_epoch(),
        position,
        magnitude: EclipseMagnitude(1.05),
        obscuration: Obscuration(1.0),
        path_width: Some(umbra_core::Kilometers(180.0)),
        central_duration: Some(200.0),
        sun_altitude: umbra_core::Degrees(70.0),
    }
}

/// 合成「中心食」SolarEclipse（central_begin/central_end=Some・kind=Total）。
fn central_eclipse(greatest_position: umbra_geo::GeoPoint) -> SolarEclipse {
    let epoch = synth_epoch();
    let u1 = contact(tt_at_hours(epoch, -1.0));
    let u4 = contact(tt_at_hours(epoch, 1.0));
    let global = GlobalCircumstances {
        kind: SolarEclipseKind::Total,
        partial_begin: None,
        central_begin: Some(u1),
        greatest: greatest_at(greatest_position),
        central_end: Some(u4),
        partial_end: None,
        gamma: 0.10,
    };
    SolarEclipse {
        event_key: "synthetic-central#0".to_string(),
        kind: SolarEclipseKind::Total,
        global,
        bessel: central_bessel(),
        metadata: metadata(),
    }
}

/// 合成「非中心/部分食」SolarEclipse（central_begin か central_end が None）。
fn noncentral_eclipse(
    greatest_position: umbra_geo::GeoPoint,
    with_begin: bool,
    with_end: bool,
) -> SolarEclipse {
    let epoch = synth_epoch();
    let central_begin = with_begin.then(|| contact(tt_at_hours(epoch, -1.0)));
    let central_end = with_end.then(|| contact(tt_at_hours(epoch, 1.0)));
    let global = GlobalCircumstances {
        kind: SolarEclipseKind::Partial,
        partial_begin: Some(contact(tt_at_hours(epoch, -1.5))),
        central_begin,
        greatest: greatest_at(greatest_position),
        central_end,
        partial_end: Some(contact(tt_at_hours(epoch, 1.5))),
        gamma: 0.10,
    };
    SolarEclipse {
        event_key: "synthetic-noncentral#0".to_string(),
        kind: SolarEclipseKind::Partial,
        global,
        bessel: central_bessel(),
        metadata: metadata(),
    }
}

// ============================================================
// FAST: 中心食＋include_limits=true → 北/南限界線 Some・非空・各点妥当
// ============================================================

/// FAST / 新規: 中心食＋include_limits 既定(true) で northern/southern_limit=Some・点列非空（≥2）・
/// 各点が妥当な緯度経度。center_line も Some（限界線と取り違えていない＝3 本とも独立に存在）。
///
/// 殺す変異: 限界線を常に None にする・空点列を返す・限界線を center_line の別名にして 1 本しか作らない・
///   include_limits を無視して常に None にする。
#[test]
fn central_eclipse_with_limits_produces_nonempty_north_and_south() {
    let engine = standard_engine(bundled_time_data());
    let eclipse = central_eclipse(geo(12.5, 77.5));

    let path = engine
        .path(&eclipse, PathOptions::default())
        .expect("中心食の path() は成功する");

    let north = path
        .northern_limit
        .as_ref()
        .expect("中心食＋include_limits=true では northern_limit=Some");
    let south = path
        .southern_limit
        .as_ref()
        .expect("中心食＋include_limits=true では southern_limit=Some");
    assert!(
        path.center_line.is_some(),
        "中心食では center_line も Some（限界線と独立）"
    );

    assert!(
        north.points.len() >= 2,
        "北限界線は ≥2 点, got {}",
        north.points.len()
    );
    assert!(
        south.points.len() >= 2,
        "南限界線は ≥2 点, got {}",
        south.points.len()
    );

    for p in north.points.iter().chain(south.points.iter()) {
        assert!(
            lat_lon_in_range(p),
            "限界線の点が妥当な緯度経度域にない: lat={} lon={}",
            lat_deg(p),
            lon_deg(p)
        );
    }
}

/// FAST / 新規: 対応サンプルで **北限緯度 ≥ 中心線緯度 ≥ 南限緯度**（高緯度側=北限・低緯度側=南限）。
/// 近似ゆえ等号許容（小マージン）。3 本が同じサンプル時刻列＝同点数で並ぶことも縛る。
///
/// 殺す変異: 北限と南限を入れ替える・中心線の外に両方とも同じ側へずらす・限界線を中心線のコピーにする
///   （緯度差ゼロ）・サンプル列を北/南で食い違わせる。
#[test]
fn northern_limit_is_north_of_center_is_north_of_southern() {
    let engine = standard_engine(bundled_time_data());
    let eclipse = central_eclipse(geo(0.0, 0.0));

    let path = engine
        .path(&eclipse, PathOptions::default())
        .expect("中心食の path() は成功する");

    let center = path.center_line.as_ref().expect("center_line=Some");
    let north = path.northern_limit.as_ref().expect("northern_limit=Some");
    let south = path.southern_limit.as_ref().expect("southern_limit=Some");

    // 同一サンプル時刻列＝3 本が同点数（取り違え・食い違いを撃破）。
    assert_eq!(
        north.points.len(),
        center.points.len(),
        "北限と中心線は同点数（同サンプル列）"
    );
    assert_eq!(
        south.points.len(),
        center.points.len(),
        "南限と中心線は同点数（同サンプル列）"
    );

    // 各対応サンプルで 北 ≥ 中心 ≥ 南（近似ゆえ微小マージン許容）。
    const EPS: f64 = 1.0e-6;
    for i in 0..center.points.len() {
        let n = lat_deg(&north.points[i]);
        let c = lat_deg(&center.points[i]);
        let s = lat_deg(&south.points[i]);
        assert!(
            n >= c - EPS,
            "サンプル{i}: 北限緯度 {n} ≥ 中心線緯度 {c}（高緯度側=北限）"
        );
        assert!(
            c >= s - EPS,
            "サンプル{i}: 中心線緯度 {c} ≥ 南限緯度 {s}（低緯度側=南限）"
        );
    }
}

/// FAST / 新規: 北限と南限が**分離**している（帯幅が正で過大でない）。代表サンプルで北限・南限の
/// 緯度差が正、かつ各対応点の概算距離が 0 でなく数百 km オーダー（数度未満＝過大でない）に収まる。
/// 脆い絶対値固定は避け、下限>0・上限を緩い帯（< ~5°/~600 km）で縛る。
///
/// 殺す変異: 北限=南限（幅ゼロ）にする・帯幅を桁違いに大きく（地球規模）/小さく（数値誤差）する・
///   ±オフセットの符号を片側に倒して幅を消す。
#[test]
fn limits_are_separated_with_plausible_band_width() {
    let engine = standard_engine(bundled_time_data());
    let eclipse = central_eclipse(geo(0.0, 0.0));

    let path = engine
        .path(&eclipse, PathOptions::default())
        .expect("中心食の path() は成功する");
    let north = path.northern_limit.as_ref().expect("northern_limit=Some");
    let south = path.southern_limit.as_ref().expect("southern_limit=Some");
    let center = path.center_line.as_ref().expect("center_line=Some");

    // 代表サンプル（中央付近）で北限・南限が分離している。
    let mid = center.points.len() / 2;
    let n = &north.points[mid];
    let s = &south.points[mid];

    let dlat = lat_deg(n) - lat_deg(s);
    assert!(
        dlat > 0.0,
        "中央サンプルで北限緯度 - 南限緯度 = {dlat} > 0（帯が分離）"
    );
    assert!(
        dlat < 5.0,
        "帯幅（緯度差 {dlat}°）が過大でない（< 5°＝数百 km オーダー）"
    );

    let band_km = great_circle_km(n, s);
    assert!(
        band_km > 1.0,
        "北限-南限の概算距離 {band_km} km が 0 でない（>1 km）"
    );
    assert!(
        band_km < 600.0,
        "北限-南限の概算距離 {band_km} km が過大でない（< 600 km）"
    );

    // 中心点が北限と南限の緯度で挟まれる（代表サンプル）。
    let c_lat = lat_deg(&center.points[mid]);
    assert!(
        lat_deg(s) <= c_lat && c_lat <= lat_deg(n),
        "中心点緯度 {c_lat} が南限 {} ≤ ・北限 {} ≥ で挟まれる",
        lat_deg(s),
        lat_deg(n)
    );
}

// ============================================================
// FAST: include_limits=false → 限界線 None（center_line は Some のまま）
// ============================================================

/// FAST / 新規: 中心食でも include_limits=false なら northern/southern_limit=None。center_line は
/// include_limits に依らず中心食なら Some（限界線フラグが中心線生成を巻き込まない）。
///
/// 殺す変異: include_limits を無視して常に限界線を作る・include_limits=false で center_line まで
///   None にする・フラグを反転して解釈する。
#[test]
fn include_limits_false_yields_no_limits_but_keeps_center_line() {
    let engine = standard_engine(bundled_time_data());
    let eclipse = central_eclipse(geo(0.0, 0.0));

    let path = engine
        .path(
            &eclipse,
            PathOptions {
                include_limits: false,
                ..PathOptions::default()
            },
        )
        .expect("中心食の path() は成功する");

    assert!(
        path.northern_limit.is_none(),
        "include_limits=false では northern_limit=None"
    );
    assert!(
        path.southern_limit.is_none(),
        "include_limits=false では southern_limit=None"
    );
    assert!(
        path.center_line.is_some(),
        "include_limits=false でも中心食なら center_line=Some"
    );
    // partial_limit / samples は本スライスでも未生成。
    assert!(path.partial_limit.is_none(), "partial_limit None");
    assert!(
        path.samples.is_empty(),
        "samples 空, got {}",
        path.samples.len()
    );
}

// ============================================================
// FAST: 非中心 → 限界線も中心線も None（include_limits=true でも）
// ============================================================

/// FAST / 新規: 非中心（central_begin か central_end が None）では include_limits=true でも
/// northern/southern_limit=None・center_line=None。「両方 Some」のときだけ限界線を作る（&& 条件）。
///
/// 殺す変異: 非中心でも限界線を作る・central_begin/end の片方だけ見て限界線を作る（|| 化）・
///   include_limits=true なら無条件で限界線を出す。
#[test]
fn noncentral_eclipse_has_no_limits() {
    let engine = standard_engine(bundled_time_data());
    let greatest_position = geo(-33.0, 151.0);

    for (with_begin, with_end) in [(false, false), (true, false), (false, true)] {
        let eclipse = noncentral_eclipse(greatest_position, with_begin, with_end);
        let path = engine
            .path(&eclipse, PathOptions::default())
            .expect("非中心の path() も成功する");

        assert!(
            path.northern_limit.is_none(),
            "非中心(begin={with_begin}, end={with_end}) では northern_limit=None"
        );
        assert!(
            path.southern_limit.is_none(),
            "非中心(begin={with_begin}, end={with_end}) では southern_limit=None"
        );
        assert!(
            path.center_line.is_none(),
            "非中心(begin={with_begin}, end={with_end}) では center_line=None"
        );
    }
}

// ============================================================
// SLOW: 実 2024-04-08 皆既を search → path() の限界線が NASA 帯幅域
// ============================================================

/// SLOW / 新規: 実エンジンで 2024-04-08 皆既を search → path()。北/南限界線が Some・各点妥当・
/// 中心線が（各点近傍で）南北限界の間にあり、最大食付近の北限と南限の概算距離が NASA 公表帯幅
/// ~197 km の妥当域（geometric 近似ゆえ広め 100–350 km）にあることを縛る。de440s 不要（解析暦）。
///
/// 殺す変異: 限界線を捏造/空にする・南北を取り違える・帯幅を桁違いにする・中心線が帯の外に出る・
///   実日食で限界線を生成しない。
#[test]
fn real_2024_eclipse_limits_match_nasa_band_width() {
    let engine = standard_engine(bundled_time_data());
    let range = umbra_core::TimeRange {
        start: utc(2024, 4, 8, 0, 0, 0.0),
        end: utc(2024, 4, 9, 0, 0, 0.0),
    };
    let eclipses = engine
        .search(range)
        .expect("2024-04-08 範囲の search は成功する");
    let eclipse = eclipses
        .iter()
        .find(|e| matches!(e.kind, SolarEclipseKind::Total))
        .expect("2024-04-08 皆既が見つかる");

    let path = engine
        .path(eclipse, PathOptions::default())
        .expect("実皆既の path() は成功する");

    let center = path
        .center_line
        .as_ref()
        .expect("皆既なので center_line=Some");
    let north = path
        .northern_limit
        .as_ref()
        .expect("皆既なので northern_limit=Some");
    let south = path
        .southern_limit
        .as_ref()
        .expect("皆既なので southern_limit=Some");

    assert!(
        north.points.len() >= 2,
        "実北限界線は ≥2 点, got {}",
        north.points.len()
    );
    assert_eq!(
        north.points.len(),
        center.points.len(),
        "北限と中心線は同サンプル列"
    );
    assert_eq!(
        south.points.len(),
        center.points.len(),
        "南限と中心線は同サンプル列"
    );

    for p in north.points.iter().chain(south.points.iter()) {
        assert!(
            lat_lon_in_range(p),
            "実限界線の点が妥当な緯度経度域にない: lat={} lon={}",
            lat_deg(p),
            lon_deg(p)
        );
    }

    // 各対応サンプルで 北 ≥ 中心 ≥ 南（近似ゆえ微小マージン許容）。
    const EPS: f64 = 1.0e-3;
    for i in 0..center.points.len() {
        let n = lat_deg(&north.points[i]);
        let c = lat_deg(&center.points[i]);
        let s = lat_deg(&south.points[i]);
        assert!(n >= c - EPS, "サンプル{i}: 北限 {n} ≥ 中心 {c}");
        assert!(c >= s - EPS, "サンプル{i}: 中心 {c} ≥ 南限 {s}");
    }

    // 最大食点に最も近いサンプル付近で帯幅 ~197 km の妥当域（geometric 近似ゆえ広め 100–350 km）。
    let g_lat = lat_deg(&path.greatest_point);
    let g_lon = lon_deg(&path.greatest_point);
    let mid = (0..center.points.len())
        .min_by(|&a, &b| {
            let da = {
                let dlat = lat_deg(&center.points[a]) - g_lat;
                let dlon = lon_deg(&center.points[a]) - g_lon;
                dlat * dlat + dlon * dlon
            };
            let db = {
                let dlat = lat_deg(&center.points[b]) - g_lat;
                let dlon = lon_deg(&center.points[b]) - g_lon;
                dlat * dlat + dlon * dlon
            };
            da.partial_cmp(&db).expect("有限距離")
        })
        .expect("中心線は非空");

    let band_km = great_circle_km(&north.points[mid], &south.points[mid]);
    assert!(
        (100.0..=350.0).contains(&band_km),
        "最大食付近の帯幅 {band_km} km が NASA ~197 km の妥当域 [100, 350] に入る（geometric 近似）"
    );
}
