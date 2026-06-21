//! M9 経路生成 第1スライス（中心線トラック）の統合テスト（ISSUE-043 S9 / 確定 path）。
//!
//! `umbra-eclipse` の**公開 API のみ**を対象とした統合テスト（tests/ 配下・別クレート境界）。
//! 対象は `EclipseEngine::path(&SolarEclipse, PathOptions) -> Result<EclipsePath, EclipseError>`。
//! 本ファイルは中心線 `center_line: Option<GeoLine>` と `greatest_point` の生成、および
//! 中心食での samples 充足（中心線との lockstep）を縛る。部分食域は未実装（常に None）。
//!
//! ## 確定セマンティクス（テストで縛る）
//! 1. `greatest_point == eclipse.global.greatest.position`（passthrough・中心/非中心問わず）。
//! 2. 中心食（`global.central_begin` と `global.central_end` がともに Some）:
//!    `center_line = Some(GeoLine)`。GeoLine は [central_begin.time_tt, central_end.time_tt] を
//!    `options.sample_interval_seconds` 刻みでサンプルし、各時刻で `bessel.at(t)` の瞬時要素から
//!    影軸地表貫通点を求めて結んだ点列。軸が地球を外す時刻はスキップ。中心食では非空（≥2 点）。
//! 3. 非中心（central_begin か central_end のいずれかが None）: `center_line = None`。
//! 4. `partial_limit` は常に None。`samples` は M9.7 以降、中心食かつ `include_limits` 既定(true) で
//!    中心線と lockstep 充足（非中心・`include_limits=false` では空）。`northern_limit`/`southern_limit`
//!    は M9.3 以降、中心食かつ `include_limits` 既定(true) のとき Some（本ファイルは存在・lockstep のみ
//!    縛り、限界線の幾何性質・samples 各フィールドは `path_limits.rs` で縛る）。非中心では限界線も None。
//! 5. `bessel.at`/影軸貫通の `RootNotBracketed` 以外の Err は伝播。
//!
//! ## テスト戦略（strict / mutation-resistant / 負荷配分）
//! - FAST（実 search 非実走）: `standard_engine(bundled_time_data())` でエンジンを作り（search しない
//!   ＝速い）、**合成 SolarEclipse**（中心食/非中心）を渡して path() の確定セマンティクスを縛る。
//!   合成 bessel poly は影軸が WGS84 に当たる小さい x,y を持ち、center_line が実際に非空になる
//!   ことを `bessel.at(t).gamma() < 1` で独立に裏取りする。
//! - SLOW（実エンジン・1 件）: `standard_engine` で 2017-08-21 皆既を search → path() で得た中心線が
//!   北米を横断する地理範囲に収まり、最大食点近傍を含むことを範囲・包含で縛る（脆い座標固定は避ける）。
//!
//! ## 期待される RED（実装前）
//! 現状 `path()` は常に `Err(EclipseError::NotImplemented)` を返すため、`.expect()` で panic する／
//! center_line を読む assert が落ちる。これが想定どおりの赤（新シンボル追加なし・既存 path() の挙動変更）。

use umbra_core::{JulianDate2, TimeInterval, TtInstant, UtcInstant};
use umbra_eclipse::{
    standard_engine, AccuracyProfile, BesselFitError, BesselianPolynomial, BesselianSource,
    CalculationMetadata, EclipseError, EclipseMagnitude, GlobalCircumstances, GlobalContact,
    GreatestEclipse, Obscuration, PathOptions, Polynomial, SolarEclipse, SolarEclipseKind,
};
use umbra_ephemeris::bundled_time_data;

// ============================================================
// 時刻 / 地理ヘルパ
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

// ============================================================
// 合成 BesselianPolynomial / SolarEclipse 構築
//
// 合成 bessel poly は影軸が WGS84 に当たる小さい x,y を持つ（gamma=√(x²+y²)≪1）。
// epoch を窓の中心に置き、x(t) を緩やかな線形ランプ、y は小さい定数にして、
// [U1,U4] 全域で gamma<1（軸が地表に当たる）かつ点が時間とともに動くようにする。
// fit_interval は [U1,U4] を内包する余裕を持たせる（at() の区間チェックを通す）。
// ============================================================

/// 合成日食の epoch（窓の中心 TT）。J2000 を借用（解析暦の評価可能域内・幾何のみ使用）。
fn synth_epoch() -> TtInstant {
    tt(2_451_545.0, 0.0)
}

/// 中心食を模す合成 bessel poly（影軸が地表に当たる小さい x,y）。
///
/// x(t)=0.05·t[hour]（±1h で ±0.05）, y(t)=0.10（定数）, d≈0.20 rad, μ=1.2,
/// l1=0.54, l2=−0.009（皆既）。gamma=√(x²+y²) は窓全域で ≲0.11 ≪ 1 ＝軸が必ず地表に当たる。
/// fit_interval は epoch±2h（[U1,U4]=epoch±1h を内包）。
fn central_bessel() -> BesselianPolynomial {
    let epoch = synth_epoch();
    let c = |v: f64| Polynomial {
        coefficients: vec![v],
    };
    BesselianPolynomial {
        epoch_tt: epoch,
        // x(t) は t[hour] の線形ランプ（点が時間とともに東西に動く）。
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
        // [U1,U4]=epoch±1h を内包する余裕（epoch±2h）。
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

/// 代表メタデータ（results.rs テストと同パターン・幾何に無関係）。
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

/// U1/U4 の全球接触点（時刻のみ幾何に効く・position は妥当な点でラベル）。
fn contact(time_tt: TtInstant) -> GlobalContact {
    GlobalContact {
        time_utc: utc(2024, 4, 8, 18, 0, 0.0),
        time_tt,
        position: geo(20.0, -100.0),
    }
}

/// 中心食の最大食（greatest_point passthrough を厳密一致で縛るため既知の固有点を仕込む）。
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
/// U1=epoch−1h, U4=epoch+1h（fit_interval epoch±2h に内包）。
/// `greatest_position` は global.greatest.position に入る既知点（passthrough 検証用）。
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
/// `with_begin`/`with_end` で 4 通り（両 None / begin だけ Some / end だけ Some）を作れる。
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

/// 緯度・経度が妥当な範囲か（合成中心線の各点のサニティ）。
fn lat_lon_in_range(p: &umbra_geo::GeoPoint) -> bool {
    let lat = p.lat.degrees().0;
    let lon = p.lon.degrees().0;
    (-90.0..=90.0).contains(&lat) && (-180.0..=180.0).contains(&lon)
}

// ============================================================
// FAST: 合成中心食の自己検証（中心線が実際に非空になる前提を裏取り）
// ============================================================

/// 合成 bessel poly が [U1,U4] 全域で gamma<1（影軸が WGS84 に当たる）ことを公開 IF
/// （`BesselianSource::at`）で独立に確認する。これが満たされなければ、中心食でも center_line が
/// 空になりうる（仕様の留意点）。本テストが落ちたら合成 x,y を縮める必要がある（テスト設計の前提固定）。
#[test]
fn synthetic_central_axis_hits_earth_across_window() {
    let bessel = central_bessel();
    let epoch = synth_epoch();
    // U1, 中央, U4（path のサンプル両端と中央）で gamma<1 を確認。
    for hours in [-1.0_f64, 0.0, 1.0] {
        let t = tt_at_hours(epoch, hours);
        let e = bessel.at(t).expect("合成区間内の at() は成功する");
        assert!(
            e.gamma() < 1.0,
            "gamma={} at t={hours}h must be <1（軸が地表に当たる前提）",
            e.gamma()
        );
    }
}

// ============================================================
// FAST: 中心食 → center_line=Some・非空・各点妥当・passthrough・限界線 Some・samples lockstep 充足
// ============================================================

/// 中心食合成: path() は center_line=Some・点列非空（≥2）・各点が妥当な緯度経度を返す。
/// greatest_point は global.greatest.position と厳密一致（passthrough）。
/// M9.3 以降は include_limits 既定(true) で northern/southern_limit=Some（本テストは存在のみ縛り、
/// 限界線の幾何性質は path_limits.rs で縛る）。partial_limit は None。M9.7 以降は samples が中心線と
/// lockstep で充足（同点数。各フィールドオラクルは path_limits.rs）。
///
/// 殺す変異: center_line を常に None にする・空点列を返す・サンプル間隔/区間の取り違えで <2 点になる・
///   greatest_point を別ソースから取る/捏造する・partial_limit や samples を捏造する・
///   中心食でも限界線を None にする。
#[test]
fn central_eclipse_produces_nonempty_center_line() {
    let engine = standard_engine(bundled_time_data());
    // global.greatest.position に入れる既知の固有点（passthrough 検証用・北米外の任意点）。
    let greatest_position = geo(12.5, 77.5);
    let eclipse = central_eclipse(greatest_position);

    let path = engine
        .path(&eclipse, PathOptions::default())
        .expect("中心食の path() は成功する（NotImplemented であってはならない）");

    // center_line=Some・非空（≥2 点）。
    let line = path
        .center_line
        .as_ref()
        .expect("中心食では center_line=Some");
    assert!(
        line.points.len() >= 2,
        "中心食の中心線は ≥2 点（軸が [U1,U4] で地表に当たる）, got {}",
        line.points.len()
    );
    // 各点が妥当な緯度経度。
    for p in &line.points {
        assert!(
            lat_lon_in_range(p),
            "中心線の点が妥当な緯度経度域にない: lat={} lon={}",
            p.lat.degrees().0,
            p.lon.degrees().0
        );
    }

    // greatest_point は global.greatest.position の passthrough（厳密一致）。
    assert_eq!(
        path.greatest_point, greatest_position,
        "greatest_point は global.greatest.position の passthrough"
    );

    // M9.3 以降: 中心食＋include_limits 既定(true) では北/南限界線は Some（中心線同様に生成）。
    // M9.7 以降: samples も中心線と lockstep で充足。`central_eclipse` は部分 phase（P1/P4）を持たない
    // ので partial_limit は None（3c-ii リボン法は部分 phase があるときのみ部分食域を組む）。詳細フィールド
    // オラクルは path_limits.rs が担保。ここでは充足と lockstep（中心線と同点数）のみ縛る。
    assert!(
        path.northern_limit.is_some(),
        "中心食＋include_limits 既定(true) では northern_limit=Some（M9.3）"
    );
    assert!(
        path.southern_limit.is_some(),
        "中心食＋include_limits 既定(true) では southern_limit=Some（M9.3）"
    );
    assert!(
        path.partial_limit.is_none(),
        "central_eclipse は部分 phase 無し（P1/P4=None）ゆえ partial_limit=None"
    );
    assert_eq!(
        path.samples.len(),
        path.center_line.as_ref().unwrap().points.len(),
        "samples は中心線と lockstep で充足（M9.7）, got {} 点",
        path.samples.len()
    );
}

/// サンプル間隔を粗くすると中心線の点数が減る（区間を sample_interval_seconds 刻みで配線している証拠）。
/// 既定 60s と粗い 1800s（30 分）を比較。U1→U4=2h の窓で、60s なら多数点・1800s なら数点に減る。
///
/// 殺す変異: sample_interval_seconds を無視して固定点数を返す・区間を間隔で割らずに常に同数を返す。
#[test]
fn coarser_interval_yields_fewer_center_line_points() {
    let engine = standard_engine(bundled_time_data());
    let eclipse = central_eclipse(geo(0.0, 0.0));

    let fine = engine
        .path(
            &eclipse,
            PathOptions {
                sample_interval_seconds: 60.0,
                ..PathOptions::default()
            },
        )
        .expect("細かい間隔の path() は成功する");
    let coarse = engine
        .path(
            &eclipse,
            PathOptions {
                sample_interval_seconds: 1800.0,
                ..PathOptions::default()
            },
        )
        .expect("粗い間隔の path() は成功する");

    let fine_n = fine
        .center_line
        .as_ref()
        .expect("中心食 center_line=Some")
        .points
        .len();
    let coarse_n = coarse
        .center_line
        .as_ref()
        .expect("中心食 center_line=Some")
        .points
        .len();

    assert!(coarse_n >= 2, "粗い間隔でも中心線は ≥2 点, got {coarse_n}");
    assert!(
        fine_n > coarse_n,
        "細かい間隔(60s)の点数 {fine_n} は粗い間隔(1800s)の点数 {coarse_n} より多い"
    );
}

// ============================================================
// FAST: 非中心/部分食 → center_line=None・greatest_point passthrough・samples 空
// ============================================================

/// 非中心/部分食合成（central_begin/central_end の少なくとも一方が None）: center_line=None。
/// greatest_point は passthrough、限界線 None・samples 空。
/// 「両 None」「begin だけ Some」「end だけ Some」の 3 通りで、**両方 Some のときだけ Some**を縛る
/// （片方でも None なら中心線を作らない = `&&` 条件）。
///
/// 殺す変異: center_begin/center_end のいずれか片方だけ見て center_line を作る（`||` 化）・
///   非中心でも center_line を Some にする・greatest_point を捏造する。
#[test]
fn noncentral_eclipse_has_no_center_line() {
    let engine = standard_engine(bundled_time_data());
    let greatest_position = geo(-33.0, 151.0);

    // (with_begin, with_end) = 両 None / begin だけ / end だけ。いずれも非中心 → center_line=None。
    for (with_begin, with_end) in [(false, false), (true, false), (false, true)] {
        let eclipse = noncentral_eclipse(greatest_position, with_begin, with_end);
        let path = engine
            .path(&eclipse, PathOptions::default())
            .expect("非中心の path() も成功する（NotImplemented であってはならない）");

        assert!(
            path.center_line.is_none(),
            "非中心(begin={with_begin}, end={with_end}) では center_line=None"
        );
        // greatest_point は passthrough。
        assert_eq!(
            path.greatest_point, greatest_position,
            "非中心でも greatest_point は passthrough"
        );
        // 本影南北限界線 None・samples 空（中心食でない）。
        assert!(path.northern_limit.is_none(), "northern_limit None");
        assert!(path.southern_limit.is_none(), "southern_limit None");
        // 部分食域 partial_limit は中心食と独立（M9 残(3) 3c-ii リボン法）: この fixture は
        // partial_begin/partial_end=Some（部分 phase あり）＋include_limits 既定(true) なので
        // 非中心でも partial_limit=Some（部分食 only でも部分食域は存在しうる・center_line とは別経路）。
        assert!(
            path.partial_limit.is_some(),
            "非中心でも部分 phase（P1/P4=Some）＋include_limits=true なら partial_limit=Some（中心食と独立）"
        );
        assert!(
            path.samples.is_empty(),
            "samples は空, got {} 点",
            path.samples.len()
        );
    }
}

// ============================================================
// SLOW: 実 2017-08-21 皆既を search → path() の中心線が北米を横断
// ============================================================

/// 実エンジンで 2017-08-21 皆既を search → 得た SolarEclipse で path() を実行し、中心線が
/// 北米を横断する地理範囲（緯度 ~30–48°N・経度 ~−125〜−78°E）に収まり、最大食点近傍を含むことを
/// 範囲・包含で縛る（脆い厳密座標固定は避ける）。de440s 不要（解析暦）。
///
/// 殺す変異: 中心線を捏造/空にする・南北/東西を取り違える・最大食点を通らない線を返す・
///   greatest_point passthrough の破れ。
#[test]
fn real_2017_eclipse_center_line_crosses_north_america() {
    let engine = standard_engine(bundled_time_data());
    // 2017-08-21 を内包する控えめな探索窓（前後 0.5 日）。
    let range = umbra_core::TimeRange {
        start: utc(2017, 8, 20, 12, 0, 0.0),
        end: utc(2017, 8, 21, 12, 0, 0.0),
    };
    let eclipses = engine
        .search(range)
        .expect("2017-08 範囲の search は成功する");
    // 範囲内に皆既日食が 1 件（2017-08-21）。
    let eclipse = eclipses
        .iter()
        .find(|e| matches!(e.kind, SolarEclipseKind::Total))
        .expect("2017-08-21 皆既が見つかる");

    let path = engine
        .path(eclipse, PathOptions::default())
        .expect("実皆既の path() は成功する");

    let line = path
        .center_line
        .as_ref()
        .expect("皆既なので center_line=Some");
    assert!(
        line.points.len() >= 2,
        "実中心線は ≥2 点, got {}",
        line.points.len()
    );

    // 全点が大域的に妥当な緯度経度。2017 皆既の本影中心線は U1（太平洋・日の出本影接触, lon≈−166°）から
    // U4（大西洋・日没本影接触）まで延びるため、北米帯に限定されない（米国横断は経路の中間部のみ）。
    for p in &line.points {
        assert!(
            lat_lon_in_range(p),
            "中心線の点が妥当な緯度経度域にない: lat={} lon={}",
            p.lat.degrees().0,
            p.lon.degrees().0
        );
    }
    // 中心線は北米を横断する（米国本土の箱内に少なくとも 1 点を含む）。全点が箱内ではなく「箱を通る」を縛る。
    let crosses_us = line.points.iter().any(|p| {
        let lat = p.lat.degrees().0;
        let lon = p.lon.degrees().0;
        (30.0..=48.0).contains(&lat) && (-125.0..=-78.0).contains(&lon)
    });
    assert!(
        crosses_us,
        "中心線が米国本土（lat30–48°N・lon−125〜−78°E）を横断する点を含む（2017 皆既は北米横断）"
    );

    // greatest_point は passthrough（search 由来 global.greatest.position と厳密一致）。
    assert_eq!(
        path.greatest_point, eclipse.global.greatest.position,
        "greatest_point は search 由来 global.greatest.position の passthrough"
    );

    // 中心線は最大食点の近傍を通る（最大食点に最も近い中心線点が ~数百 km 内）。
    let g_lat = path.greatest_point.lat.degrees().0;
    let g_lon = path.greatest_point.lon.degrees().0;
    let min_deg = line
        .points
        .iter()
        .map(|p| {
            let dlat = p.lat.degrees().0 - g_lat;
            let dlon = p.lon.degrees().0 - g_lon;
            (dlat * dlat + dlon * dlon).sqrt()
        })
        .fold(f64::INFINITY, f64::min);
    assert!(
        min_deg < 3.0,
        "中心線は最大食点の近傍を通る（最近点 {min_deg}° < 3°）"
    );

    // 実皆既＋include_limits 既定(true) では北/南限界線も Some（M9.3）。実皆既は部分 phase（P1/P4）も
    // 持つので partial_limit も Some（M9 残(3) 3c-ii リボン法・部分食域は中心食と独立に存在）。
    // M9.7 以降: samples も中心線と lockstep で充足（詳細フィールドは path_limits.rs で縛る）。
    assert!(
        path.northern_limit.is_some(),
        "皆既＋既定 で northern_limit=Some（M9.3）"
    );
    assert!(
        path.southern_limit.is_some(),
        "皆既＋既定 で southern_limit=Some（M9.3）"
    );
    assert!(
        path.partial_limit.is_some(),
        "実皆既は部分 phase を持つので partial_limit=Some（3c-ii リボン法）"
    );
    assert_eq!(
        path.samples.len(),
        path.center_line.as_ref().unwrap().points.len(),
        "samples は中心線と lockstep で充足（M9.7）, got {} 点",
        path.samples.len()
    );

    // エラー型が公開参照可能（path の戻り値型に使う）。
    fn _accepts_err(_e: EclipseError) {}
}
