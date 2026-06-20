//! M9.2 GeoJSON 出力（幾何プリミティブ）の統合テスト（umbra-geo 公開 API のみ）。
//!
//! 対象は `GeoPoint::geojson_geometry()`（Point）と `GeoLine::geojson_geometry()`
//! （LineString / 日付変更線跨ぎで MultiLineString に分割）。
//!
//! ## 確定セマンティクス（テストで縛る）
//! - Point: `{"type":"Point","coordinates":[lon_deg, lat_deg]}`。**座標順は [経度, 緯度]**
//!   （GeoJSON RFC 7946）。値は度。
//! - LineString（跨ぎ無し・全連続点で |Δlon| ≤ 180）:
//!   `{"type":"LineString","coordinates":[[lon,lat],...]}`（全点・順序保持・[lon,lat]）。
//! - MultiLineString（ある連続 2 点で |Δlon| > 180）: 跨ぎ位置で切り、各セグメントは連続
//!   [lon,lat] 列。閾値ちょうど（|Δlon| = 180）は**跨ぎでない**（LineString）。
//! - 退行（0/1 点）: LineString で coordinates 長 0/1（不正だがそのまま・panic しない）。
//!
//! ## テスト戦略（strict / mutation-resistant / FAST）
//! 全 FAST（実エンジン不要）。`serde_json::from_str` ではなく直接返る `Value` を構造で検証する。
//! lon と lat を別値にして座標取り違え変異を殺す。完全文字列一致は避け Value 構造で縛る。
//!
//! ## 期待される RED（実装前）
//! `GeoPoint::geojson_geometry` / `GeoLine::geojson_geometry` が未定義のためメソッド解決不能
//! （E0599）でコンパイルできない。これが想定どおりの赤。

use serde_json::Value;
use umbra_geo::{GeoLine, GeoPoint};

/// 数値比較の許容（度の桁落ち程度）。
const EPS: f64 = 1e-9;

/// 度の GeoPoint を作るヘルパ。
fn pt(lat: f64, lon: f64) -> GeoPoint {
    GeoPoint::from_degrees(lat, lon).expect("有効な緯度経度")
}

/// `Value` の座標ペア `[lon, lat]` を f64 タプルで取り出す（型・長さも検証）。
fn coord_pair(v: &Value) -> (f64, f64) {
    let arr = v.as_array().expect("座標ペアは配列");
    assert_eq!(arr.len(), 2, "座標ペアは [lon, lat] の長さ 2");
    let lon = arr[0].as_f64().expect("lon は数値");
    let lat = arr[1].as_f64().expect("lat は数値");
    (lon, lat)
}

/// 浮動小数の近接。
fn close(a: f64, b: f64) -> bool {
    (a - b).abs() < EPS
}

// ============================================================
// GeoPoint::geojson_geometry — Point・[lon, lat] 順
// ============================================================

/// `GeoPoint::geojson_geometry` は `type="Point"`・`coordinates=[lon, lat]` を返す。
/// lat=12.5・lon=77.5 という**非対称値**を与え、coordinates[0]=lon=77.5・
/// coordinates[1]=lat=12.5 を厳密に縛る。
///
/// 殺す変異: coordinates を [lat, lon] 逆順にする・type 文字列の改変・lon/lat の取り違え。
#[test]
fn geo_point_geojson_is_point_with_lon_lat_order() {
    let g = pt(12.5, 77.5).geojson_geometry();
    // type は "Point"。
    assert_eq!(g["type"], Value::String("Point".to_string()), "type=Point");
    // coordinates は [lon, lat]=[77.5, 12.5]（順序厳密）。
    let (lon, lat) = coord_pair(&g["coordinates"]);
    assert!(close(lon, 77.5), "coordinates[0]=lon=77.5, got {lon}");
    assert!(close(lat, 12.5), "coordinates[1]=lat=12.5, got {lat}");
}

/// 負の緯度・経度でも [lon, lat] 順が保たれる（符号の取り違え・逆順を殺す）。
/// lat=-33.0・lon=151.0（シドニー近傍）→ coordinates=[151.0, -33.0]。
#[test]
fn geo_point_geojson_preserves_sign_and_order() {
    let g = pt(-33.0, 151.0).geojson_geometry();
    let (lon, lat) = coord_pair(&g["coordinates"]);
    assert!(close(lon, 151.0), "lon=151.0, got {lon}");
    assert!(close(lat, -33.0), "lat=-33.0, got {lat}");
}

// ============================================================
// GeoLine::geojson_geometry — 跨ぎ無し LineString
// ============================================================

/// 跨ぎ無し（全連続点で |Δlon| ≤ 180）の折れ線は `type="LineString"`。
/// coordinates は全点・順序保持・各 [lon, lat]。
/// lon=[0, 30, 60]（緩やかに東進）・lat=[0, 10, -5]。
///
/// 殺す変異: 点の脱落・順序入れ替え・[lat,lon] 逆順・type を MultiLineString に固定。
#[test]
fn geo_line_geojson_no_crossing_is_linestring() {
    let line = GeoLine::new(vec![pt(0.0, 0.0), pt(10.0, 30.0), pt(-5.0, 60.0)]);
    let g = line.geojson_geometry();
    assert_eq!(
        g["type"],
        Value::String("LineString".to_string()),
        "跨ぎ無しは LineString"
    );
    let coords = g["coordinates"].as_array().expect("coordinates は配列");
    // 全 3 点が出る（順序保持）。
    assert_eq!(coords.len(), 3, "全 3 点が出る");
    let expected = [(0.0, 0.0), (30.0, 10.0), (60.0, -5.0)]; // (lon, lat)
    for (i, exp) in expected.iter().enumerate() {
        let (lon, lat) = coord_pair(&coords[i]);
        assert!(close(lon, exp.0), "点{i} lon={} expected {}", lon, exp.0);
        assert!(close(lat, exp.1), "点{i} lat={} expected {}", lat, exp.1);
    }
}

// ============================================================
// GeoLine::geojson_geometry — 跨ぎ MultiLineString
// ============================================================

/// 1 箇所跨ぎ（lon=[170, -170, -160]）は `type="MultiLineString"`・2 セグメント。
/// 170→−170 は |Δ|=340>180 ゆえ跨ぎ → seg0=[[170,..]]・seg1=[[-170,..],[-160,..]]。
///
/// 殺す変異: 跨ぎを検出せず LineString のまま・分割位置のずれ・セグメント数の誤り・
///   座標の [lat,lon] 逆順・点の脱落。
#[test]
fn geo_line_geojson_single_crossing_is_multilinestring_two_segments() {
    let line = GeoLine::new(vec![pt(1.0, 170.0), pt(2.0, -170.0), pt(3.0, -160.0)]);
    let g = line.geojson_geometry();
    assert_eq!(
        g["type"],
        Value::String("MultiLineString".to_string()),
        "跨ぎ有りは MultiLineString"
    );
    let segs = g["coordinates"].as_array().expect("coordinates は配列");
    assert_eq!(segs.len(), 2, "跨ぎ 1 箇所 → 2 セグメント");

    // seg0 は跨ぎ前の 1 点 [170, 1]。
    let s0 = segs[0].as_array().expect("seg0 は配列");
    assert_eq!(s0.len(), 1, "seg0 は 1 点");
    let (lon, lat) = coord_pair(&s0[0]);
    assert!(close(lon, 170.0) && close(lat, 1.0), "seg0[0]=[170,1]");

    // seg1 は跨ぎ後の 2 点 [-170, 2], [-160, 3]。
    let s1 = segs[1].as_array().expect("seg1 は配列");
    assert_eq!(s1.len(), 2, "seg1 は 2 点");
    let (lon0, lat0) = coord_pair(&s1[0]);
    assert!(close(lon0, -170.0) && close(lat0, 2.0), "seg1[0]=[-170,2]");
    let (lon1, lat1) = coord_pair(&s1[1]);
    assert!(close(lon1, -160.0) && close(lat1, 3.0), "seg1[1]=[-160,3]");
}

/// 複数跨ぎ（lon=[170, -170, 170]）は 3 セグメント（各跨ぎで切る）。
/// 170→−170（|Δ|=340>180）と −170→170（|Δ|=340>180）の 2 跨ぎ → 3 セグ（各 1 点）。
///
/// 殺す変異: 最初の跨ぎだけ切って 2 セグにする・全部 1 本に戻す・跨ぎ回数の数え誤り。
#[test]
fn geo_line_geojson_double_crossing_is_three_segments() {
    let line = GeoLine::new(vec![pt(1.0, 170.0), pt(2.0, -170.0), pt(3.0, 170.0)]);
    let g = line.geojson_geometry();
    assert_eq!(
        g["type"],
        Value::String("MultiLineString".to_string()),
        "2 跨ぎは MultiLineString"
    );
    let segs = g["coordinates"].as_array().expect("coordinates は配列");
    assert_eq!(segs.len(), 3, "跨ぎ 2 箇所 → 3 セグメント");
    // 各セグメントは 1 点ずつ（lon=170, -170, 170）。
    let expected_lon = [170.0, -170.0, 170.0];
    for (i, &elon) in expected_lon.iter().enumerate() {
        let seg = segs[i].as_array().expect("セグは配列");
        assert_eq!(seg.len(), 1, "seg{i} は 1 点");
        let (lon, _lat) = coord_pair(&seg[0]);
        assert!(close(lon, elon), "seg{i} lon={lon} expected {elon}");
    }
}

// ============================================================
// GeoLine::geojson_geometry — 跨ぎ閾値の境界（|Δlon|=180 は跨ぎでない）
// ============================================================

/// |Δlon| = 180 ちょうどは**跨ぎでない**＝ LineString（1 本）のまま。
/// lon=[0, -180]（正規化後 -180）。|Δ| = |-180 − 0| = 180 ≤ 180 → 切らない。
///
/// 殺す変異: 閾値を `>=` にして 180 ちょうどでも切る（オフバイワン）・常に MultiLineString。
#[test]
fn geo_line_geojson_exactly_180_is_not_crossing() {
    let line = GeoLine::new(vec![pt(0.0, 0.0), pt(0.0, -180.0)]);
    let g = line.geojson_geometry();
    assert_eq!(
        g["type"],
        Value::String("LineString".to_string()),
        "|Δlon|=180 ちょうどは跨ぎでない → LineString"
    );
    let coords = g["coordinates"].as_array().expect("coordinates は配列");
    assert_eq!(coords.len(), 2, "2 点とも 1 本の LineString に残る");
}

/// |Δlon| > 180（180 超）は跨ぎ＝ MultiLineString に分割。
/// lon=[10, -171]: |Δ| = |-171 − 10| = 181 > 180 → 2 セグメント（各 1 点）。
/// 上の「180 ちょうど」テストと対にして、閾値が「> 180」であることを両側から縛る。
///
/// 殺す変異: 閾値を `> 180` でなく `>= 181` 等にして 181 を跨ぎとしない・閾値方向の誤り。
#[test]
fn geo_line_geojson_just_over_180_is_crossing() {
    let line = GeoLine::new(vec![pt(0.0, 10.0), pt(0.0, -171.0)]);
    let g = line.geojson_geometry();
    assert_eq!(
        g["type"],
        Value::String("MultiLineString".to_string()),
        "|Δlon|=181 > 180 は跨ぎ → MultiLineString"
    );
    let segs = g["coordinates"].as_array().expect("coordinates は配列");
    assert_eq!(segs.len(), 2, "181° 差で 2 セグメント");
    assert_eq!(segs[0].as_array().expect("seg0").len(), 1, "seg0 は 1 点");
    assert_eq!(segs[1].as_array().expect("seg1").len(), 1, "seg1 は 1 点");
}

// ============================================================
// GeoLine::geojson_geometry — 退行（0/1 点）
// ============================================================

/// 0 点の折れ線は LineString で coordinates 長 0（panic しない）。
/// 不正な GeoJSON だが本スライスはそのまま出す（呼び出し側責務）。
///
/// 殺す変異: 退行で panic する・点を捏造して長さ非 0 にする・type を変える。
#[test]
fn geo_line_geojson_empty_is_linestring_with_no_coords() {
    let line = GeoLine::new(Vec::new());
    let g = line.geojson_geometry();
    assert_eq!(
        g["type"],
        Value::String("LineString".to_string()),
        "0 点でも LineString"
    );
    let coords = g["coordinates"].as_array().expect("coordinates は配列");
    assert_eq!(coords.len(), 0, "0 点 → coordinates 長 0");
}

/// 1 点の折れ線は LineString で coordinates 長 1（panic しない・跨ぎ判定の対象外）。
///
/// 殺す変異: 1 点を 0 点に落とす・捏造で増やす・1 点で MultiLineString にする。
#[test]
fn geo_line_geojson_single_point_is_linestring_with_one_coord() {
    let line = GeoLine::new(vec![pt(7.0, 8.0)]);
    let g = line.geojson_geometry();
    assert_eq!(
        g["type"],
        Value::String("LineString".to_string()),
        "1 点でも LineString"
    );
    let coords = g["coordinates"].as_array().expect("coordinates は配列");
    assert_eq!(coords.len(), 1, "1 点 → coordinates 長 1");
    let (lon, lat) = coord_pair(&coords[0]);
    assert!(
        close(lon, 8.0) && close(lat, 7.0),
        "唯一点=[8,7]（[lon,lat]）"
    );
}
