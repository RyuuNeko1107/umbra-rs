//! M9.2 経路の GeoJSON 出力（FeatureCollection 組立）の統合テスト（umbra-eclipse 公開 API のみ）。
//!
//! 対象は `EclipsePath::to_geojson() -> Result<String, serde_json::Error>`。
//! `EclipsePath` を公開フィールドで直接構築し（実エンジン非実走＝FAST）、直列化結果を
//! `serde_json::from_str` で `Value` に戻して構造を検証する。
//!
//! ## 確定セマンティクス（テストで縛る）
//! - トップは `{"type":"FeatureCollection","features":[...]}`。pretty・**末尾改行**。
//! - features に **greatest_point の Point feature**
//!   （`{"type":"Feature","geometry":{Point...},"properties":{"role":"greatest"}}`）を含む。
//! - `center_line=Some` のとき **center_line の Feature**（geometry は LineString/MultiLineString・
//!   `properties.role="center_line"`）を含む。`center_line=None` なら center_line feature を含まない
//!   （greatest のみ）。
//! - **限界線 Feature（M9.5）**: `northern_limit=Some` → Feature（geometry=GeoLine の GeoJSON・
//!   `properties.role="northern_limit"`）、`southern_limit=Some` → 同（`role="southern_limit"`）。
//!   いずれも `None` なら当該 feature を出さない。
//! - **partial_limit は features に出さない**（常に None・GeoPolygon の GeoJSON 化は後続スライス）。
//! - feature 順序は **greatest → center_line → northern_limit → southern_limit**（決定的）。
//! - samples は features に出さない（本スライス）。
//! - GeoJSON 座標順は [経度, 緯度]（lon, lat）。
//!
//! ## テスト戦略（strict / mutation-resistant / FAST）
//! 全 FAST。Value 構造で縛り、脆い完全文字列一致は避ける（末尾改行のみ文字列で確認）。
//! greatest_point は lon≠lat の非対称値で取り違え・逆順変異を殺す。
//!
//! ## 期待される RED（実装前）
//! `EclipsePath::to_geojson` が未定義のためメソッド解決不能（E0599）でコンパイルできない。
//! これが想定どおりの赤。

use serde_json::Value;
use umbra_eclipse::EclipsePath;
use umbra_geo::{GeoLine, GeoPoint};

/// 数値比較の許容。
const EPS: f64 = 1e-9;

/// 度の GeoPoint。
fn pt(lat: f64, lon: f64) -> GeoPoint {
    GeoPoint::from_degrees(lat, lon).expect("有効な緯度経度")
}

/// 浮動小数の近接。
fn close(a: f64, b: f64) -> bool {
    (a - b).abs() < EPS
}

/// center_line=Some の中心食 EclipsePath を作る（北/南限界線 Some・部分食域 None・samples 空）。
/// greatest_point は lat=12.5・lon=77.5 の非対称値（逆順変異検出用）。
/// 北限界線と南限界線は **互いに異なる座標**（北=北寄り正緯度、南=南寄り負緯度）にして、
/// role↔座標 の取り違え（北南入替）を殺せるようにする。
fn central_path() -> EclipsePath {
    EclipsePath {
        center_line: Some(GeoLine::new(vec![
            pt(0.0, 0.0),
            pt(10.0, 30.0),
            pt(-5.0, 60.0),
        ])),
        // 北限界線: 正緯度・東経（中心線とも南限界線とも別座標）。
        northern_limit: Some(GeoLine::new(vec![pt(11.0, 31.0), pt(12.0, 61.0)])),
        // 南限界線: 負緯度・別経度（北限界線と区別可能）。
        southern_limit: Some(GeoLine::new(vec![pt(-11.0, 29.0), pt(-12.0, 59.0)])),
        partial_limit: None,
        greatest_point: pt(12.5, 77.5),
        samples: Vec::new(),
    }
}

/// feature の geometry.coordinates（LineString）を (lon, lat) ペア列で取り出す。
fn line_coords(feature: &Value) -> Vec<(f64, f64)> {
    feature["geometry"]["coordinates"]
        .as_array()
        .expect("LineString coordinates は配列")
        .iter()
        .map(|c| {
            let a = c.as_array().expect("ペアは配列");
            (a[0].as_f64().expect("lon"), a[1].as_f64().expect("lat"))
        })
        .collect()
}

/// features 配列の role 順序（出現順）を取り出す。
fn role_order(root: &Value) -> Vec<String> {
    root["features"]
        .as_array()
        .expect("features は配列")
        .iter()
        .map(|f| {
            f["properties"]["role"]
                .as_str()
                .expect("role は文字列")
                .to_string()
        })
        .collect()
}

/// center_line=None の非中心 EclipsePath（greatest_point のみが出る想定）。
fn noncentral_path() -> EclipsePath {
    EclipsePath {
        center_line: None,
        northern_limit: None,
        southern_limit: None,
        partial_limit: None,
        greatest_point: pt(-33.0, 151.0),
        samples: Vec::new(),
    }
}

/// FeatureCollection の features 配列から、指定 role を持つ feature を集める。
fn features_with_role<'a>(root: &'a Value, role: &str) -> Vec<&'a Value> {
    root["features"]
        .as_array()
        .expect("features は配列")
        .iter()
        .filter(|f| f["properties"]["role"] == Value::String(role.to_string()))
        .collect()
}

/// feature の geometry.coordinates の単一ペア [lon, lat] を取り出す（Point 用）。
fn point_coord(feature: &Value) -> (f64, f64) {
    let arr = feature["geometry"]["coordinates"]
        .as_array()
        .expect("Point coordinates は配列");
    assert_eq!(arr.len(), 2, "Point coordinates は [lon, lat]");
    (arr[0].as_f64().expect("lon"), arr[1].as_f64().expect("lat"))
}

// ============================================================
// to_geojson — FeatureCollection の骨格・末尾改行・valid JSON
// ============================================================

/// `to_geojson` は valid JSON で `type="FeatureCollection"`・`features` 配列・末尾改行を持つ。
///
/// 殺す変異: type 文字列の改変・features キーの欠落・末尾改行の脱落・pretty でない直列化。
#[test]
fn to_geojson_is_feature_collection_with_trailing_newline() {
    let s = central_path().to_geojson().expect("直列化は成功する");
    // 末尾改行（文字列レベルで縛る・他は Value で縛る）。
    assert!(s.ends_with('\n'), "出力は末尾改行で終わる");
    let root: Value = serde_json::from_str(&s).expect("出力は valid JSON");
    assert_eq!(
        root["type"],
        Value::String("FeatureCollection".to_string()),
        "type=FeatureCollection"
    );
    assert!(root["features"].is_array(), "features は配列");
}

// ============================================================
// to_geojson — 中心食: greatest + center_line + northern_limit + southern_limit
// ============================================================

/// 中心食（center_line/northern/southern すべて Some）の features は
/// greatest + center_line + northern_limit + southern_limit の**ちょうど 4 件**。
/// 各 role はちょうど 1 件・partial_limit は None なので feature を出さない。
///
/// 殺す変異: 限界線 feature を出さない（旧 2 件仕様）・partial_limit feature を捏造する・
///   いずれかの role を重複/欠落させる・greatest や center_line を出さない・samples を混ぜる。
#[test]
fn to_geojson_central_has_greatest_center_line_and_both_limits() {
    let s = central_path().to_geojson().expect("直列化は成功する");
    let root: Value = serde_json::from_str(&s).expect("valid JSON");

    // 各 role ちょうど 1 件。
    assert_eq!(
        features_with_role(&root, "greatest").len(),
        1,
        "greatest 1 件"
    );
    assert_eq!(
        features_with_role(&root, "center_line").len(),
        1,
        "center_line 1 件"
    );
    assert_eq!(
        features_with_role(&root, "northern_limit").len(),
        1,
        "northern_limit 1 件"
    );
    assert_eq!(
        features_with_role(&root, "southern_limit").len(),
        1,
        "southern_limit 1 件"
    );
    // partial_limit feature は出さない（None・後続スライス）。
    assert!(
        features_with_role(&root, "partial_limit").is_empty(),
        "partial_limit feature は出さない"
    );

    // greatest は Point geometry・Feature 型。
    let greatest = features_with_role(&root, "greatest");
    assert_eq!(
        greatest[0]["type"],
        Value::String("Feature".to_string()),
        "greatest は Feature"
    );
    assert_eq!(
        greatest[0]["geometry"]["type"],
        Value::String("Point".to_string()),
        "greatest の geometry は Point"
    );
    // center_line は LineString（跨ぎ無し中心線）。
    assert_eq!(
        features_with_role(&root, "center_line")[0]["geometry"]["type"],
        Value::String("LineString".to_string()),
        "center_line は LineString（跨ぎ無し）"
    );

    // features はちょうど 4 件。
    let all = root["features"].as_array().expect("features は配列");
    assert_eq!(
        all.len(),
        4,
        "features は greatest + center_line + northern_limit + southern_limit の 4 件, got {}",
        all.len()
    );
}

/// features の出現順は greatest → center_line → northern_limit → southern_limit（決定的）。
///
/// 殺す変異: feature 追加順を入れ替える・限界線を center_line より前に出す・
///   北南の追加順を逆にする。
#[test]
fn to_geojson_feature_order_is_deterministic() {
    let s = central_path().to_geojson().expect("直列化は成功する");
    let root: Value = serde_json::from_str(&s).expect("valid JSON");
    assert_eq!(
        role_order(&root),
        vec![
            "greatest".to_string(),
            "center_line".to_string(),
            "northern_limit".to_string(),
            "southern_limit".to_string(),
        ],
        "feature 順序は greatest→center_line→northern_limit→southern_limit"
    );
}

/// northern_limit / southern_limit の geometry 座標が**対応する** GeoLine と一致する。
/// 北と南は別座標（北=正緯度 [31,11],[61,12]、南=負緯度 [29,-11],[59,-12]・[lon,lat]）なので、
/// role と座標の対応が崩れる（北南取り違え）変異を殺す。
///
/// 殺す変異: northern_limit に southern の座標を入れる（role↔geometry 入替）・
///   [lat,lon] 逆順・点の脱落/順序入替・別の線を出す。
#[test]
fn to_geojson_limits_geometry_matches_corresponding_geoline() {
    let s = central_path().to_geojson().expect("直列化は成功する");
    let root: Value = serde_json::from_str(&s).expect("valid JSON");

    let north = features_with_role(&root, "northern_limit");
    assert_eq!(north.len(), 1, "northern_limit 1 件");
    assert_eq!(
        north[0]["geometry"]["type"],
        Value::String("LineString".to_string()),
        "northern_limit は LineString（跨ぎ無し）"
    );
    let n_coords = line_coords(north[0]);
    let n_expected = [(31.0, 11.0), (61.0, 12.0)]; // (lon, lat)
    assert_eq!(n_coords.len(), n_expected.len(), "北限界線の点数一致");
    for (i, exp) in n_expected.iter().enumerate() {
        assert!(close(n_coords[i].0, exp.0), "北 点{i} lon");
        assert!(close(n_coords[i].1, exp.1), "北 点{i} lat");
    }

    let south = features_with_role(&root, "southern_limit");
    assert_eq!(south.len(), 1, "southern_limit 1 件");
    let s_coords = line_coords(south[0]);
    let s_expected = [(29.0, -11.0), (59.0, -12.0)]; // (lon, lat)
    assert_eq!(s_coords.len(), s_expected.len(), "南限界線の点数一致");
    for (i, exp) in s_expected.iter().enumerate() {
        assert!(close(s_coords[i].0, exp.0), "南 点{i} lon");
        assert!(close(s_coords[i].1, exp.1), "南 点{i} lat");
    }

    // 念のため北と南が別座標であることを確認（取り違え検出の前提）。
    assert_ne!(n_coords, s_coords, "北限界線と南限界線は別座標");
}

/// greatest_point の Point 座標が [lon, lat] 順で正しい（lat=12.5・lon=77.5 → [77.5, 12.5]）。
///
/// 殺す変異: coordinates を [lat, lon] 逆順にする・別フィールドから座標を取る。
#[test]
fn to_geojson_greatest_point_uses_lon_lat_order() {
    let s = central_path().to_geojson().expect("直列化は成功する");
    let root: Value = serde_json::from_str(&s).expect("valid JSON");
    let greatest = features_with_role(&root, "greatest");
    assert_eq!(greatest.len(), 1, "greatest はちょうど 1 つ");
    let (lon, lat) = point_coord(greatest[0]);
    assert!(close(lon, 77.5), "greatest lon=77.5, got {lon}");
    assert!(close(lat, 12.5), "greatest lat=12.5, got {lat}");
}

/// center_line feature の座標列が元の GeoLine と一致する（順序・[lon,lat]）。
/// lon=[0,30,60]・lat=[0,10,-5] → coordinates=[[0,0],[30,10],[60,-5]]。
///
/// 殺す変異: 中心線の点を脱落/順序入れ替え・[lat,lon] 逆順・別の線を出す。
#[test]
fn to_geojson_center_line_coordinates_match_geoline() {
    let s = central_path().to_geojson().expect("直列化は成功する");
    let root: Value = serde_json::from_str(&s).expect("valid JSON");
    let center = features_with_role(&root, "center_line");
    assert_eq!(center.len(), 1, "center_line はちょうど 1 つ");
    let coords = center[0]["geometry"]["coordinates"]
        .as_array()
        .expect("LineString coordinates は配列");
    let expected = [(0.0, 0.0), (30.0, 10.0), (60.0, -5.0)]; // (lon, lat)
    assert_eq!(coords.len(), expected.len(), "中心線の点数一致");
    for (i, exp) in expected.iter().enumerate() {
        let arr = coords[i].as_array().expect("ペアは配列");
        let lon = arr[0].as_f64().expect("lon");
        let lat = arr[1].as_f64().expect("lat");
        assert!(close(lon, exp.0), "点{i} lon={lon} expected {}", exp.0);
        assert!(close(lat, exp.1), "点{i} lat={lat} expected {}", exp.1);
    }
}

// ============================================================
// to_geojson — 非中心（center_line=None）: greatest のみ
// ============================================================

/// 非中心（center_line=None）の features は role=greatest の Point **のみ**
/// （center_line feature を含まない）。
///
/// 殺す変異: center_line=None でも center_line feature を捏造する・greatest を出さない・
///   features を空にする。
#[test]
fn to_geojson_noncentral_has_only_greatest_feature() {
    let s = noncentral_path().to_geojson().expect("直列化は成功する");
    let root: Value = serde_json::from_str(&s).expect("valid JSON");

    let all = root["features"].as_array().expect("features は配列");
    assert_eq!(all.len(), 1, "非中心では feature は greatest のみの 1 件");

    let greatest = features_with_role(&root, "greatest");
    assert_eq!(greatest.len(), 1, "role=greatest はちょうど 1 つ");
    let center = features_with_role(&root, "center_line");
    assert!(
        center.is_empty(),
        "center_line=None では center_line feature を含まない"
    );
    // northern/southern_limit=None では限界線 feature も含まない（include_limits 相当）。
    assert!(
        features_with_role(&root, "northern_limit").is_empty(),
        "northern_limit=None では northern_limit feature を含まない"
    );
    assert!(
        features_with_role(&root, "southern_limit").is_empty(),
        "southern_limit=None では southern_limit feature を含まない"
    );

    // greatest の座標は [lon, lat]=[151.0, -33.0]。
    let (lon, lat) = point_coord(greatest[0]);
    assert!(close(lon, 151.0), "greatest lon=151.0, got {lon}");
    assert!(close(lat, -33.0), "greatest lat=-33.0, got {lat}");
}

// ============================================================
// to_geojson — 跨ぎ中心線が MultiLineString になる（umbra-geo 分割が経由でも効く）
// ============================================================

/// 中心線が日付変更線を跨ぐ（lon=[170, -170, -160]・lat=[1,2,3]）と center_line feature の
/// geometry は MultiLineString になり、±180 補間（M9.5）が EclipsePath 経由でも効く。
/// 東進跨ぎ 170→−170: t=(180−170)/(360−340)=0.5, lat_c=1+0.5·(2−1)=1.5。
/// seg0 末尾=[+180,1.5]・seg1 先頭=[−180,1.5]。
///
/// 殺す変異: 経路出力で GeoLine の分割を使わず LineString に固定する・跨ぎを無視する・
///   ±180 補間点を挿入しない（隙間旧仕様）・境界符号の取り違え。
#[test]
fn to_geojson_antimeridian_center_line_is_multilinestring() {
    let path = EclipsePath {
        center_line: Some(GeoLine::new(vec![
            pt(1.0, 170.0),
            pt(2.0, -170.0),
            pt(3.0, -160.0),
        ])),
        northern_limit: None,
        southern_limit: None,
        partial_limit: None,
        greatest_point: pt(0.0, 175.0),
        samples: Vec::new(),
    };
    let s = path.to_geojson().expect("直列化は成功する");
    let root: Value = serde_json::from_str(&s).expect("valid JSON");
    let center = features_with_role(&root, "center_line");
    assert_eq!(center.len(), 1, "center_line はちょうど 1 つ");
    assert_eq!(
        center[0]["geometry"]["type"],
        Value::String("MultiLineString".to_string()),
        "跨ぎ中心線の geometry は MultiLineString"
    );
    // 跨ぎ 1 箇所 → 2 セグメント。
    let segs = center[0]["geometry"]["coordinates"]
        .as_array()
        .expect("coordinates は配列");
    assert_eq!(segs.len(), 2, "跨ぎ 1 箇所 → 2 セグメント");

    // ±180 補間（オラクル独立計算）: 東進 t=0.5・lat_c=1.5。
    let t = (180.0 - 170.0) / (360.0 + (-170.0 - 170.0));
    let lat_c = 1.0 + t * (2.0 - 1.0);
    assert!(close(lat_c, 1.5), "オラクル lat_c=1.5, got {lat_c}");

    // seg0 末尾=[+180, lat_c]。
    let s0 = segs[0].as_array().expect("seg0");
    let last0 = s0.last().expect("seg0 末尾").as_array().expect("ペア");
    assert!(
        close(last0[0].as_f64().expect("lon"), 180.0),
        "東進: seg0 末尾の境界は +180"
    );
    assert!(
        close(last0[1].as_f64().expect("lat"), lat_c),
        "seg0 末尾 lat=lat_c"
    );

    // seg1 先頭=[−180, lat_c]。
    let s1 = segs[1].as_array().expect("seg1");
    let first1 = s1[0].as_array().expect("ペア");
    assert!(
        close(first1[0].as_f64().expect("lon"), -180.0),
        "東進: seg1 先頭の境界は −180"
    );
    assert!(
        close(first1[1].as_f64().expect("lat"), lat_c),
        "seg1 先頭 lat=lat_c"
    );
}
