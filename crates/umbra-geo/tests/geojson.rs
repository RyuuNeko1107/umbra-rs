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
//!   **±180 補間（M9.5）**: 跨ぎ点 (lon1,lat1)→(lon2,lat2) では交点緯度 lat_c を ±180 子午線上に
//!   線形補間し、前セグメント末尾に `[±180, lat_c]`・次セグメント先頭に `[∓180, lat_c]` を追加する
//!   （隙間を埋める）。
//!   - 東進 Δlon<−180（例 170→−170）: 前末尾 `[+180, lat_c]`、次先頭 `[−180, lat_c]`。
//!     `t = (180 − lon1) / (360 + Δlon)`, `lat_c = lat1 + t·(lat2 − lat1)`。
//!   - 西進 Δlon>+180（例 −170→170）: 前末尾 `[−180, lat_c]`、次先頭 `[+180, lat_c]`。
//!     `t = (lon1 + 180) / (360 − Δlon)`, `lat_c = lat1 + t·(lat2 − lat1)`。
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
use umbra_geo::{GeoLine, GeoPoint, GeoPolygon};

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

/// 東進跨ぎ（Δlon<−180、例 170→−170）の交点緯度 lat_c を**独立に**計算する。
/// `t = (180 − lon1) / (360 + Δlon)`, `lat_c = lat1 + t·(lat2 − lat1)`。
/// （実装式を写経せず、オラクル根拠の式を別途展開して t・lat 補間を縛る。）
fn lat_c_east(lon1: f64, lat1: f64, lon2: f64, lat2: f64) -> f64 {
    let dlon = lon2 - lon1;
    let t = (180.0 - lon1) / (360.0 + dlon);
    lat1 + t * (lat2 - lat1)
}

/// 西進跨ぎ（Δlon>+180、例 −170→170）の交点緯度 lat_c を**独立に**計算する。
/// `t = (lon1 + 180) / (360 − Δlon)`, `lat_c = lat1 + t·(lat2 − lat1)`。
fn lat_c_west(lon1: f64, lat1: f64, lon2: f64, lat2: f64) -> f64 {
    let dlon = lon2 - lon1;
    let t = (lon1 + 180.0) / (360.0 - dlon);
    lat1 + t * (lat2 - lat1)
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

/// 東進 1 箇所跨ぎ（lon=[170, -170, -160]・lat=[1, 2, 3]）は `type="MultiLineString"`・2 セグメント。
/// 170→−170 は Δlon=−340<−180 ゆえ**東進**跨ぎ。±180 補間（M9.5）で:
/// t=(180−170)/(360−340)=10/20=0.5, lat_c=1+0.5·(2−1)=1.5。
/// → seg0=[[170,1],[+180,1.5]]、seg1=[[−180,1.5],[-170,2],[-160,3]]。
///
/// 殺す変異: 跨ぎを検出せず LineString のまま・±180 補間点を挿入しない（隙間が残る旧仕様）・
///   境界符号の取り違え（前末尾を −180／次先頭を +180 にする）・lat_c の補間誤り・
///   分割位置のずれ・座標の [lat,lon] 逆順・点の脱落。
#[test]
fn geo_line_geojson_east_crossing_inserts_pm180_interpolation() {
    let line = GeoLine::new(vec![pt(1.0, 170.0), pt(2.0, -170.0), pt(3.0, -160.0)]);
    let g = line.geojson_geometry();
    assert_eq!(
        g["type"],
        Value::String("MultiLineString".to_string()),
        "跨ぎ有りは MultiLineString"
    );
    let segs = g["coordinates"].as_array().expect("coordinates は配列");
    assert_eq!(segs.len(), 2, "跨ぎ 1 箇所 → 2 セグメント");

    // オラクル：交点緯度（独立計算）。t=0.5・lat_c=1.5 になるはず。
    let lat_c = lat_c_east(170.0, 1.0, -170.0, 2.0);
    assert!(close(lat_c, 1.5), "オラクル lat_c=1.5, got {lat_c}");

    // seg0 = [[170,1], [+180, lat_c]]（末尾に +180 子午線上の補間点）。
    let s0 = segs[0].as_array().expect("seg0 は配列");
    assert_eq!(s0.len(), 2, "seg0 は元 1 点 + 補間 1 点 = 2 点");
    let (lon, lat) = coord_pair(&s0[0]);
    assert!(close(lon, 170.0) && close(lat, 1.0), "seg0[0]=[170,1]");
    let (lon_b, lat_b) = coord_pair(&s0[1]);
    assert!(
        close(lon_b, 180.0),
        "東進: 前セグメント末尾の境界は +180, got {lon_b}"
    );
    assert!(
        close(lat_b, lat_c),
        "seg0 末尾 lat=lat_c={lat_c}, got {lat_b}"
    );

    // seg1 = [[−180, lat_c], [-170,2], [-160,3]]（先頭に −180 子午線上の補間点）。
    let s1 = segs[1].as_array().expect("seg1 は配列");
    assert_eq!(s1.len(), 3, "seg1 は 補間 1 点 + 元 2 点 = 3 点");
    let (lon_f, lat_f) = coord_pair(&s1[0]);
    assert!(
        close(lon_f, -180.0),
        "東進: 次セグメント先頭の境界は −180, got {lon_f}"
    );
    assert!(
        close(lat_f, lat_c),
        "seg1 先頭 lat=lat_c={lat_c}, got {lat_f}"
    );
    let (lon0, lat0) = coord_pair(&s1[1]);
    assert!(close(lon0, -170.0) && close(lat0, 2.0), "seg1[1]=[-170,2]");
    let (lon1, lat1) = coord_pair(&s1[2]);
    assert!(close(lon1, -160.0) && close(lat1, 3.0), "seg1[2]=[-160,3]");
}

/// 東進跨ぎの**非対称**ケース（t≠0.5・lat 非対称）で t の分母/分子・lat 補間を個別に縛る。
/// lon=[170, -175]（Δlon=−345<−180・東進）・lat=[2, 8]。
/// t=(180−170)/(360−345)=10/15=2/3≈0.6667, lat_c=2+(2/3)·(8−2)=2+4=6.0。
/// 緯度を 2→8（Δlat=6）の**非対称値**にして t·Δlat を縛り、t=0.5 固定・分母を 360 固定等の変異を殺す。
///
/// 殺す変異: t の分母を (360+Δlon) でなく定数/別式にする・分子 (180−lon1) を誤る・
///   lat_c を中点固定（0.5）にする・lat1/lat2 取り違え。
#[test]
fn geo_line_geojson_east_crossing_asymmetric_t_and_lat() {
    let line = GeoLine::new(vec![pt(2.0, 170.0), pt(8.0, -175.0)]);
    let g = line.geojson_geometry();
    let segs = g["coordinates"].as_array().expect("coordinates は配列");
    assert_eq!(segs.len(), 2, "1 跨ぎ → 2 セグメント");

    let lat_c = lat_c_east(170.0, 2.0, -175.0, 8.0);
    assert!(close(lat_c, 6.0), "オラクル lat_c=6.0, got {lat_c}");

    // seg0 末尾 = [+180, lat_c]。
    let s0 = segs[0].as_array().expect("seg0 は配列");
    let (lon_b, lat_b) = coord_pair(s0.last().expect("seg0 末尾"));
    assert!(close(lon_b, 180.0), "前末尾境界 +180, got {lon_b}");
    assert!(close(lat_b, lat_c), "前末尾 lat=lat_c={lat_c}, got {lat_b}");

    // seg1 先頭 = [−180, lat_c]。
    let s1 = segs[1].as_array().expect("seg1 は配列");
    let (lon_f, lat_f) = coord_pair(s1.first().expect("seg1 先頭"));
    assert!(close(lon_f, -180.0), "次先頭境界 −180, got {lon_f}");
    assert!(close(lat_f, lat_c), "次先頭 lat=lat_c={lat_c}, got {lat_f}");
}

/// 西進 1 箇所跨ぎ（lon=[−170, 175]・lat=[3, 9]）は MultiLineString・2 セグメント。
/// −170→175 は Δlon=+345>180 ゆえ**西進**跨ぎ。±180 補間（M9.5）で:
/// t=(−170+180)/(360−345)=10/15=2/3, lat_c=3+(2/3)·(9−3)=3+4=7.0。
/// → seg0 末尾=[−180, lat_c]、seg1 先頭=[+180, lat_c]（東進と境界符号が**逆**）。
///
/// 殺す変異: 西進で境界符号を東進と同じ（前+180/次−180）にする取り違え・
///   西進 t の分子を (lon1+180) でなく (180−lon1) にする・lat 補間誤り・
///   西進を跨ぎと認識しない（Δlon>+180 の判定脱落）。
#[test]
fn geo_line_geojson_west_crossing_inserts_pm180_interpolation() {
    let line = GeoLine::new(vec![pt(3.0, -170.0), pt(9.0, 175.0)]);
    let g = line.geojson_geometry();
    assert_eq!(
        g["type"],
        Value::String("MultiLineString".to_string()),
        "西進跨ぎも MultiLineString"
    );
    let segs = g["coordinates"].as_array().expect("coordinates は配列");
    assert_eq!(segs.len(), 2, "西進 1 跨ぎ → 2 セグメント");

    let lat_c = lat_c_west(-170.0, 3.0, 175.0, 9.0);
    assert!(close(lat_c, 7.0), "オラクル lat_c=7.0, got {lat_c}");

    // seg0 = [[-170,3], [−180, lat_c]]（西進: 前末尾は −180）。
    let s0 = segs[0].as_array().expect("seg0 は配列");
    let (lon0, lat0) = coord_pair(&s0[0]);
    assert!(close(lon0, -170.0) && close(lat0, 3.0), "seg0[0]=[-170,3]");
    let (lon_b, lat_b) = coord_pair(s0.last().expect("seg0 末尾"));
    assert!(
        close(lon_b, -180.0),
        "西進: 前セグメント末尾の境界は −180, got {lon_b}"
    );
    assert!(close(lat_b, lat_c), "前末尾 lat=lat_c={lat_c}, got {lat_b}");

    // seg1 = [[+180, lat_c], [175,9]]（西進: 次先頭は +180）。
    let s1 = segs[1].as_array().expect("seg1 は配列");
    let (lon_f, lat_f) = coord_pair(&s1[0]);
    assert!(
        close(lon_f, 180.0),
        "西進: 次セグメント先頭の境界は +180, got {lon_f}"
    );
    assert!(close(lat_f, lat_c), "次先頭 lat=lat_c={lat_c}, got {lat_f}");
    let (lon1, lat1) = coord_pair(s1.last().expect("seg1 末尾"));
    assert!(close(lon1, 175.0) && close(lat1, 9.0), "seg1 末尾=[175,9]");
}

/// 複数跨ぎ（lon=[170, -170, 170]・lat=[1, 4, 7]）は 3 セグメント（各跨ぎで切る）。
/// 170→−170（Δ=−340・東進）と −170→170（Δ=+340・西進）の 2 跨ぎ → 3 セグ。
/// 各セグメントは元 1 点 + 跨ぎごとの補間点を持つ:
///   seg0=[[170,1],[+180,lat_c0]]、seg1=[[−180,lat_c0],[-170,4],[−180,lat_c1]]、
///   seg2=[[+180,lat_c1],[170,7]]。
/// lat_c0=lat_c_east(170,1,-170,4)=2.5、lat_c1=lat_c_west(-170,4,170,7)=5.5（独立計算）。
///
/// 殺す変異: 最初の跨ぎだけ切って 2 セグにする・全部 1 本に戻す・跨ぎ回数の数え誤り・
///   2 つ目（西進）で境界符号を取り違える・補間点を入れない（隙間旧仕様）。
#[test]
fn geo_line_geojson_double_crossing_inserts_pm180_per_crossing() {
    let line = GeoLine::new(vec![pt(1.0, 170.0), pt(4.0, -170.0), pt(7.0, 170.0)]);
    let g = line.geojson_geometry();
    assert_eq!(
        g["type"],
        Value::String("MultiLineString".to_string()),
        "2 跨ぎは MultiLineString"
    );
    let segs = g["coordinates"].as_array().expect("coordinates は配列");
    assert_eq!(segs.len(), 3, "跨ぎ 2 箇所 → 3 セグメント");

    let lat_c0 = lat_c_east(170.0, 1.0, -170.0, 4.0);
    let lat_c1 = lat_c_west(-170.0, 4.0, 170.0, 7.0);
    assert!(close(lat_c0, 2.5), "lat_c0=2.5, got {lat_c0}");
    assert!(close(lat_c1, 5.5), "lat_c1=5.5, got {lat_c1}");

    // seg0: [[170,1], [+180, lat_c0]]。
    let s0 = segs[0].as_array().expect("seg0");
    assert_eq!(s0.len(), 2, "seg0 は 2 点");
    let (l, la) = coord_pair(&s0[0]);
    assert!(close(l, 170.0) && close(la, 1.0), "seg0[0]=[170,1]");
    let (l, la) = coord_pair(&s0[1]);
    assert!(
        close(l, 180.0) && close(la, lat_c0),
        "seg0 末尾=[+180,lat_c0]"
    );

    // seg1: [[−180, lat_c0], [-170,4], [−180, lat_c1]]。
    let s1 = segs[1].as_array().expect("seg1");
    assert_eq!(s1.len(), 3, "seg1 は 3 点（両端に補間）");
    let (l, la) = coord_pair(&s1[0]);
    assert!(
        close(l, -180.0) && close(la, lat_c0),
        "seg1 先頭=[−180,lat_c0]"
    );
    let (l, la) = coord_pair(&s1[1]);
    assert!(close(l, -170.0) && close(la, 4.0), "seg1[1]=[-170,4]");
    let (l, la) = coord_pair(&s1[2]);
    assert!(
        close(l, -180.0) && close(la, lat_c1),
        "seg1 末尾=[−180,lat_c1]"
    );

    // seg2: [[+180, lat_c1], [170,7]]。
    let s2 = segs[2].as_array().expect("seg2");
    assert_eq!(s2.len(), 2, "seg2 は 2 点");
    let (l, la) = coord_pair(&s2[0]);
    assert!(
        close(l, 180.0) && close(la, lat_c1),
        "seg2 先頭=[+180,lat_c1]"
    );
    let (l, la) = coord_pair(&s2[1]);
    assert!(close(l, 170.0) && close(la, 7.0), "seg2[1]=[170,7]");
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

/// 西進側の境界 Δlon = +180 ちょうど（lon1 < lon2）も**跨ぎでない**＝ LineString（1 本）。
/// lon=[−90, +90]（正規化不変・Δlon = +90 − (−90) = +180）。`from_degrees(_, 90.0)` /
/// `from_degrees(_, -90.0)` はそのまま 90 / −90。lat は非対称値（10, 20）にして座標一致も縛る。
/// 上の「東進ちょうど 180」テスト（exactly_180・Δlon=−180）と対称に、西進ちょうど +180 を縛る。
///
/// 殺す変異: 西進 `delta > 180.0` → `>=`（Δlon=+180 を誤跨ぎ分割し MultiLineString・補間点挿入）。
#[test]
fn geo_line_geojson_exactly_180_west_is_not_crossing() {
    let line = GeoLine::new(vec![pt(10.0, -90.0), pt(20.0, 90.0)]);
    let g = line.geojson_geometry();
    assert_eq!(
        g["type"],
        Value::String("LineString".to_string()),
        "西進 Δlon=+180 ちょうどは跨ぎでない → LineString"
    );
    let coords = g["coordinates"].as_array().expect("coordinates は配列");
    // 元の 2 点のみ（±180 補間点を挿入しない）。
    assert_eq!(coords.len(), 2, "元の 2 点のみ・補間点を挿入しない");
    let (lon0, lat0) = coord_pair(&coords[0]);
    assert!(
        close(lon0, -90.0) && close(lat0, 10.0),
        "coords[0]=[-90,10]"
    );
    let (lon1, lat1) = coord_pair(&coords[1]);
    assert!(close(lon1, 90.0) && close(lat1, 20.0), "coords[1]=[90,20]");
}

/// |Δlon| > 180（180 超）は跨ぎ＝ MultiLineString に分割。
/// lon=[10, -171]: Δlon = −181 < −180（東進）→ 2 セグメント。
/// ±180 補間（M9.5）で各セグメントは元 1 点 + 境界補間 1 点 = 2 点になる。
/// lat は両端 0 なので lat_c=0。前末尾=[+180,0]・次先頭=[−180,0]。
/// 上の「180 ちょうど」テストと対にして、閾値が「> 180」であることを両側から縛る。
///
/// 殺す変異: 閾値を `> 180` でなく `>= 181` 等にして 181 を跨ぎとしない・閾値方向の誤り・
///   補間点の挿入脱落（旧仕様で seg 長が 1 になる）。
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

    // 東進: seg0 末尾=[+180, 0]、seg1 先頭=[−180, 0]（lat_c=0）。
    let s0 = segs[0].as_array().expect("seg0");
    assert_eq!(s0.len(), 2, "seg0 は 元 1 点 + 補間 1 点");
    let (lon_b, lat_b) = coord_pair(s0.last().expect("seg0 末尾"));
    assert!(
        close(lon_b, 180.0) && close(lat_b, 0.0),
        "seg0 末尾=[+180,0]"
    );

    let s1 = segs[1].as_array().expect("seg1");
    assert_eq!(s1.len(), 2, "seg1 は 補間 1 点 + 元 1 点");
    let (lon_f, lat_f) = coord_pair(&s1[0]);
    assert!(
        close(lon_f, -180.0) && close(lat_f, 0.0),
        "seg1 先頭=[−180,0]"
    );
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

// ============================================================
// GeoPolygon::geojson_geometry — Polygon 構造・[lon,lat]順・閉リング・環向き
// ============================================================

/// 度の GeoPolygon リング（点列）を作るヘルパ。各点は (lat, lon)。
fn ring(points: &[(f64, f64)]) -> Vec<GeoPoint> {
    points.iter().map(|&(lat, lon)| pt(lat, lon)).collect()
}

/// `Value` のリング（座標ペア配列）を (lon, lat) タプル列で取り出す。
fn ring_coords(v: &Value) -> Vec<(f64, f64)> {
    v.as_array()
        .expect("リングは座標配列")
        .iter()
        .map(coord_pair)
        .collect()
}

/// (lon, lat) 平面の符号付き面積（shoelace）を**独立に**計算する（オラクル）。
/// CCW で正・CW で負（RFC 7946 §3.1.6 の右手則）。実装の判定式を写経せず別途展開して環向きを縛る。
fn signed_area(coords: &[(f64, f64)]) -> f64 {
    let n = coords.len();
    if n < 3 {
        return 0.0;
    }
    let mut s = 0.0;
    for i in 0..n {
        let (x1, y1) = coords[i];
        let (x2, y2) = coords[(i + 1) % n];
        s += x1 * y2 - x2 * y1;
    }
    s / 2.0
}

/// 単一外環の `GeoPolygon`（**非閉**入力・CCW 向き）→
/// `{"type":"Polygon","coordinates":[ring]}`。ring は閉（先頭==末尾複製）・[lon,lat] 順。
/// 外環を CCW（符号付き面積>0）で与え、出力も CCW のまま（規約に既に合致＝反転しない）。
/// 非対称な三角形（lon/lat が各点で別値）で [lat,lon] 逆順変異を殺す。
///
/// 殺す変異: type を Polygon 以外に固定・coordinates のネスト段数誤り（リング配列でなく座標直置き）・
///   [lat,lon] 逆順・閉リングを作らない（先頭複製の脱落）・既に CCW なのに反転する。
#[test]
fn geo_polygon_geojson_single_ring_closes_and_keeps_lon_lat() {
    // CCW 三角形（(lon,lat)=(0,0)->(40,0)->(20,30)）。shoelace>0。
    let outer = ring(&[(0.0, 0.0), (0.0, 40.0), (30.0, 20.0)]);
    let g = GeoPolygon::new(vec![outer]).geojson_geometry();

    assert_eq!(
        g["type"],
        Value::String("Polygon".to_string()),
        "type=Polygon"
    );
    let rings = g["coordinates"]
        .as_array()
        .expect("coordinates はリング配列");
    assert_eq!(rings.len(), 1, "外環 1 つだけ → リング 1 つ");

    let r0 = ring_coords(&rings[0]);
    // 閉リング: 3 点 + 先頭複製 = 4 点。
    assert_eq!(r0.len(), 4, "非閉 3 点入力 → 末尾に先頭複製で 4 点");
    assert_eq!(r0.first(), r0.last(), "先頭==末尾（閉リング）");

    // [lon,lat] 順で元の点が並ぶ（最初の 3 点）。
    let expected = [(0.0, 0.0), (40.0, 0.0), (20.0, 30.0)]; // (lon, lat)
    for (i, exp) in expected.iter().enumerate() {
        assert!(
            close(r0[i].0, exp.0),
            "点{i} lon={} expected {}",
            r0[i].0,
            exp.0
        );
        assert!(
            close(r0[i].1, exp.1),
            "点{i} lat={} expected {}",
            r0[i].1,
            exp.1
        );
    }
    // 末尾は先頭の複製。
    assert!(
        close(r0[3].0, 0.0) && close(r0[3].1, 0.0),
        "末尾=先頭=[0,0]"
    );

    // 外環は CCW（符号付き面積>0・独立計算オラクル）。
    assert!(
        signed_area(&r0) > 0.0,
        "外環は CCW（面積>0）, got {}",
        signed_area(&r0)
    );
}

/// 既に閉じている入力リング（先頭==末尾）は**二重化しない**（末尾に更に複製を足さない）。
/// CCW 閉四角形を与え、出力リング長が入力と同じであることを縛る。
///
/// 殺す変異: 閉入力でも無条件に先頭を足して二重化する・閉判定を逆にする。
#[test]
fn geo_polygon_geojson_already_closed_ring_not_duplicated() {
    // CCW 閉四角形（先頭==末尾）。(lon,lat): (0,0)->(40,0)->(40,20)->(0,20)->(0,0)。
    let closed = ring(&[
        (0.0, 0.0),
        (0.0, 40.0),
        (20.0, 40.0),
        (20.0, 0.0),
        (0.0, 0.0),
    ]);
    let g = GeoPolygon::new(vec![closed]).geojson_geometry();
    let rings = g["coordinates"].as_array().expect("リング配列");
    let r0 = ring_coords(&rings[0]);
    // 入力 5 点（既に閉）→ 二重化せず 5 点のまま。
    assert_eq!(r0.len(), 5, "既に閉なら二重化しない（5 点のまま）");
    assert_eq!(r0.first(), r0.last(), "先頭==末尾は保持");
    // 末尾の次（=4 番目）が先頭と一致し、6 点目を作っていない。
    assert!(close(r0[4].0, 0.0) && close(r0[4].1, 0.0), "末尾=[0,0]");
}

/// CW で与えた外環は **CCW に正規化**して出力（点列を反転）。
/// CW 三角形（符号付き面積<0）を入力し、出力リングの符号付き面積>0 を縛る。
/// 反転後も閉リング（先頭==末尾）であること・元の点集合が保たれること（捏造なし）を確認。
///
/// 殺す変異: 環向き判定の符号反転（CW を CCW と誤判定し反転しない）・反転処理の欠落・
///   外環を CW のまま出す。
#[test]
fn geo_polygon_geojson_cw_outer_ring_normalized_to_ccw() {
    // CW 三角形（(lon,lat)=(0,0)->(20,30)->(40,0)）。shoelace<0。
    let outer = ring(&[(0.0, 0.0), (30.0, 20.0), (0.0, 40.0)]);
    // 入力が CW であることをオラクルで確認（前提）。
    let input_coords: Vec<(f64, f64)> = vec![(0.0, 0.0), (20.0, 30.0), (40.0, 0.0)];
    assert!(signed_area(&input_coords) < 0.0, "入力外環は CW（前提）");

    let g = GeoPolygon::new(vec![outer]).geojson_geometry();
    let rings = g["coordinates"].as_array().expect("リング配列");
    let r0 = ring_coords(&rings[0]);
    assert_eq!(r0.first(), r0.last(), "出力も閉リング");
    // CW 入力 → CCW 出力（面積>0）。
    assert!(
        signed_area(&r0) > 0.0,
        "CW 外環は CCW に正規化（面積>0）, got {}",
        signed_area(&r0)
    );
}

/// 穴（rings[1..]）は **CW に正規化**して出力。外環 CCW・穴を CCW で与え、
/// 出力で外環は CCW（面積>0）・穴は CW（面積<0）になることを縛る。
/// 外環と穴で逆の規約が適用されること（取り違え検出）を非対称座標で確認。
///
/// 殺す変異: 外環と穴に同じ向き規約を適用する・穴の規約を CCW にする・
///   外環/穴のインデックス取り違え（rings[0] を穴扱い）。
#[test]
fn geo_polygon_geojson_hole_normalized_to_cw_while_outer_ccw() {
    // 外環 CCW 大三角形・穴 CCW 小三角形（両方 CCW 入力）。
    let outer = ring(&[(0.0, 0.0), (0.0, 60.0), (30.0, 30.0)]); // CCW
    let hole = ring(&[(10.0, 20.0), (10.0, 40.0), (20.0, 30.0)]); // CCW
    let g = GeoPolygon::new(vec![outer, hole]).geojson_geometry();

    let rings = g["coordinates"].as_array().expect("リング配列");
    assert_eq!(rings.len(), 2, "外環 + 穴 = 2 リング");

    let r0 = ring_coords(&rings[0]);
    let r1 = ring_coords(&rings[1]);
    // 外環は CCW（面積>0）。
    assert!(signed_area(&r0) > 0.0, "外環 CCW, got {}", signed_area(&r0));
    // 穴は CW（面積<0）。
    assert!(signed_area(&r1) < 0.0, "穴 CW, got {}", signed_area(&r1));

    // 穴の座標集合は元の穴に由来（外環座標でない＝取り違えなし）。
    // 外環の経度は 0/60 系・穴の経度は 20/40 系（重ならない）。
    assert!(
        r1.iter().all(|&(lon, _)| (15.0..=45.0).contains(&lon)),
        "穴の経度は穴入力由来（20/30/40 系）"
    );
}

/// 空 `rings` → `coordinates: []`（リングを捏造しない）。type は Polygon のまま。
///
/// 殺す変異: 空でリングを捏造して非空にする・panic・type を変える。
#[test]
fn geo_polygon_geojson_empty_rings_is_empty_coordinates() {
    let g = GeoPolygon::new(Vec::new()).geojson_geometry();
    assert_eq!(
        g["type"],
        Value::String("Polygon".to_string()),
        "空でも type=Polygon"
    );
    let rings = g["coordinates"].as_array().expect("coordinates は配列");
    assert_eq!(rings.len(), 0, "空 rings → coordinates 長 0");
}

/// ±180 を跨ぐ単一リングでも v1 では **分割しない**＝ `type=="Polygon"`（MultiPolygon にしない）。
/// 経度がそのまま並ぶ単一 Polygon のまま（クリッピングは後続スライス）。
/// lon=[170, -170, -160] を含むリングで MultiPolygon 化しないこと・経度がそのまま並ぶことを縛る。
///
/// 殺す変異: 跨ぎを検出して MultiPolygon にする（v1 仕様逸脱）・跨ぎ経度を勝手に補間/分割する。
#[test]
fn geo_polygon_geojson_antimeridian_stays_single_polygon() {
    // ±180 を跨ぐ三角形（非閉入力・CCW になるよう順序選定は問わない）。
    let outer = ring(&[(1.0, 170.0), (2.0, -170.0), (3.0, -160.0)]);
    let g = GeoPolygon::new(vec![outer]).geojson_geometry();

    // v1: 跨いでも Polygon のまま（MultiPolygon にしない）。
    assert_eq!(
        g["type"],
        Value::String("Polygon".to_string()),
        "±180 跨ぎでも v1 は単一 Polygon（MultiPolygon にしない）"
    );
    let rings = g["coordinates"].as_array().expect("リング配列");
    assert_eq!(rings.len(), 1, "単一リング（分割しない）");

    // 元の経度（170, -170, -160）がそのまま含まれる（±180 補間点を挿入していない）。
    let r0 = ring_coords(&rings[0]);
    let lons: Vec<f64> = r0.iter().map(|&(lon, _)| lon).collect();
    assert!(lons.iter().any(|&l| close(l, 170.0)), "経度 170 がそのまま");
    assert!(
        lons.iter().any(|&l| close(l, -170.0)),
        "経度 -170 がそのまま"
    );
    assert!(
        lons.iter().any(|&l| close(l, -160.0)),
        "経度 -160 がそのまま"
    );
    // ±180 ちょうどの補間点を作っていない（GeoLine と違いポリゴンは v1 で補間しない）。
    assert!(
        !lons.iter().any(|&l| close(l.abs(), 180.0)),
        "v1 は ±180 補間点を挿入しない"
    );
}

// ============================================================
// GeoPolygon::geojson_geometry — 環向き正規化の境界（面積ゼロ・shoelace 符号）
//   mutation 工程: 生存 3 変異を撃つ。各テストは「reverse される/されない」を
//   出力座標の**順序**で弁別する（向き判定は独立 shoelace オラクルで手計算する）。
// ============================================================

/// 面積ちょうど 0 の**外環**（lon 軸上の共線 3 点）は反転しない（入力順のまま閉じる）。
///
/// 撃つ変異: `geometry.rs:164` 外環側 `area < 0.0` → `area <= 0.0`。
/// 共線リング `[lon,lat] = (0,0),(2,0),(1,0)` を閉じると closed=`[(0,0),(2,0),(1,0),(0,0)]`。
/// shoelace（独立計算）はちょうど 0（lat が全点 0 ゆえ Σ(x1·y2−x2·y1)=0）。
/// original の条件 `want_ccw && area < 0.0` は `0 < 0` が偽 → **reverse しない**ので出力は
/// 入力順 `lon=[0,2,1,0]`。変異 `<=` は `0 <= 0` が真 → reverse し `lon=[0,1,2,0]` になる。
/// 出力リングの 2 番目の経度（original=2 / 変異=1）で弁別する。
///
/// `ring()` は (lat, lon) を取るので共線（lat=0 固定・lon=0/2/1）は (0,0),(0,2),(0,1)。
#[test]
fn geo_polygon_geojson_zero_area_outer_not_reversed() {
    // lon 軸上の共線 3 点（lat=0 固定）。面積ちょうど 0。
    let outer = ring(&[(0.0, 0.0), (0.0, 2.0), (0.0, 1.0)]);
    let g = GeoPolygon::new(vec![outer]).geojson_geometry();
    let rings = g["coordinates"].as_array().expect("リング配列");
    assert_eq!(rings.len(), 1, "外環 1 つ");
    let r0 = ring_coords(&rings[0]);

    // 閉リング: 共線 3 点 + 先頭複製 = 4 点。
    assert_eq!(r0.len(), 4, "共線 3 点入力 → 末尾に先頭複製で 4 点");
    assert_eq!(r0.first(), r0.last(), "先頭==末尾（閉リング）");

    // 独立オラクル: shoelace（共線）は厳密に 0。面積 0 → original は reverse しない。
    assert!(
        signed_area(&r0).abs() < EPS,
        "共線外環の面積はちょうど 0, got {}",
        signed_area(&r0)
    );

    // 入力順のまま（reverse しない）: lon=[0,2,1,0]。
    // 変異 `<=` だと reverse され lon=[0,1,2,0] になる（2 番目の経度 2→1 で弁別）。
    let lons: Vec<f64> = r0.iter().map(|&(lon, _)| lon).collect();
    let expected_lons = [0.0, 2.0, 1.0, 0.0];
    for (i, exp) in expected_lons.iter().enumerate() {
        assert!(
            close(lons[i], *exp),
            "外環は面積 0 で反転しない: 点{i} lon={} expected {}（変異 `<=` は反転して {:?} になる）",
            lons[i],
            exp,
            [0.0, 1.0, 2.0, 0.0]
        );
    }
}

/// 面積ちょうど 0 の**穴**（rings[1]・lon 軸上の共線 3 点）は反転しない（入力順のまま）。
///
/// 撃つ変異: `geometry.rs:164` 穴側 `area > 0.0` → `area >= 0.0`。
/// 穴 `[lon,lat] = (0,0),(2,0),(1,0)` を閉じると closed の shoelace はちょうど 0。
/// original の条件 `!want_ccw && area > 0.0` は `0 > 0` が偽 → **reverse しない**ので穴は
/// 入力順 `lon=[0,2,1,0]`。変異 `>=` は `0 >= 0` が真 → reverse し `lon=[0,1,2,0]`。
/// 穴リング（rings[1]）の 2 番目の経度で弁別する。外環（rings[0]）は非退行 CCW 三角形を与える。
#[test]
fn geo_polygon_geojson_zero_area_hole_not_reversed() {
    // 外環は非退行 CCW 三角形（面積>0・正規化で不変）。共線穴と経度域が重ならないよう離す。
    let outer = ring(&[(0.0, 100.0), (0.0, 160.0), (30.0, 130.0)]);
    // 穴: lon 軸（lat=0 固定）上の共線 3 点・面積ちょうど 0。lon=10/12/11 系。
    let hole = ring(&[(0.0, 10.0), (0.0, 12.0), (0.0, 11.0)]);
    let g = GeoPolygon::new(vec![outer, hole]).geojson_geometry();
    let rings = g["coordinates"].as_array().expect("リング配列");
    assert_eq!(rings.len(), 2, "外環 + 穴 = 2 リング");

    let r1 = ring_coords(&rings[1]);
    assert_eq!(r1.len(), 4, "共線穴 3 点 → 閉じて 4 点");
    assert!(
        signed_area(&r1).abs() < EPS,
        "共線穴の面積はちょうど 0, got {}",
        signed_area(&r1)
    );

    // 穴は面積 0 で反転しない: lon=[10,12,11,10]。
    // 変異 `>=` だと reverse され lon=[10,11,12,10] になる（2 番目 12→11 で弁別）。
    let lons: Vec<f64> = r1.iter().map(|&(lon, _)| lon).collect();
    let expected_lons = [10.0, 12.0, 11.0, 10.0];
    for (i, exp) in expected_lons.iter().enumerate() {
        assert!(
            close(lons[i], *exp),
            "穴は面積 0 で反転しない: 点{i} lon={} expected {}（変異 `>=` は反転して {:?} になる）",
            lons[i],
            exp,
            [10.0, 11.0, 12.0, 10.0]
        );
    }
}

/// shoelace の `x1*y2 − x2*y1` の `*`→`+` 変異を環向き正規化の結果（出力座標順）で弁別する。
///
/// 撃つ変異: `geometry.rs:181` `x1 * y2 − x2 * y1` → `(x1 + y2) − x2 * y1`。
/// 外環 `[lon,lat] = (0,0),(0,1),(1,0)` を閉じると closed=`[(0,0),(0,1),(1,0),(0,0)]`。
/// **original 式（手計算）**: 各辺 `x1·y2 − x2·y1` の和
///   (0·1 − 0·0) + (0·0 − 1·1) + (1·0 − 0·0) + 0 = 0 − 1 + 0 = **−1.0**（CW・面積<0）。
/// **変異 `+` 式（手計算）**: 各辺 `(x1+y2) − x2·y1`
///   (0+1 − 0·0) + (0+0 − 1·1) + (1+0 − 0·0) + (0+0 − 0·0) = 1 − 1 + 1 + 0 = **+1.0**（>0）。
/// 外環（want_ccw）の条件 `area < 0.0`:
///   - original area=−1.0 < 0 → **reverse する** → 反転後の閉リングは
///     `[(0,0),(1,0),(0,1),(0,0)]` ＝ lon=[0,1,0,0]。
///   - 変異 area=+1.0 < 0 偽 → **reverse しない** → 入力順の閉リング lon=[0,0,1,0]。
///
/// 出力リングの経度列（2 番目: original=1 / 変異=0、3 番目: original=0 / 変異=1）で弁別する。
///
/// `ring()` は (lat, lon) を取るので入力点 `[lon,lat]=(0,0),(0,1),(1,0)` は (lat,lon)=(0,0),(1,0),(0,1)。
#[test]
fn geo_polygon_geojson_shoelace_product_sign_decides_winding() {
    // [lon,lat] = (0,0),(0,1),(1,0)。(lat,lon) で与える。
    let outer = ring(&[(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)]);
    let g = GeoPolygon::new(vec![outer]).geojson_geometry();
    let rings = g["coordinates"].as_array().expect("リング配列");
    assert_eq!(rings.len(), 1, "外環 1 つ");
    let r0 = ring_coords(&rings[0]);
    assert_eq!(r0.len(), 4, "非閉 3 点 → 閉じて 4 点");
    assert_eq!(r0.first(), r0.last(), "閉リング");

    // 独立オラクル: 入力閉リングの真の面積は CW（負）。original は reverse して CCW 出力にする。
    let input_closed = [(0.0, 0.0), (0.0, 1.0), (1.0, 0.0), (0.0, 0.0)];
    assert!(
        signed_area(&input_closed) < 0.0,
        "入力外環は CW（面積<0・前提）, got {}",
        signed_area(&input_closed)
    );
    // 出力は CCW に正規化されている（面積>0）。変異 `+` は面積を +1 と誤算し反転を省くため CW のまま。
    assert!(
        signed_area(&r0) > 0.0,
        "外環は CCW に正規化（面積>0）, got {}（変異 `*`→`+` は誤って反転せず CW のまま残す）",
        signed_area(&r0)
    );

    // 出力座標順でも弁別: original は反転後 lon=[0,1,0,0]。変異は入力順 lon=[0,0,1,0]。
    // 2・3 番目の経度（original=1,0 / 変異=0,1）で確実に分かれる。
    let lons: Vec<f64> = r0.iter().map(|&(lon, _)| lon).collect();
    let expected_lons = [0.0, 1.0, 0.0, 0.0];
    for (i, exp) in expected_lons.iter().enumerate() {
        assert!(
            close(lons[i], *exp),
            "original は CW 入力を反転: 点{i} lon={} expected {}（変異 `+` は反転せず {:?}）",
            lons[i],
            exp,
            [0.0, 0.0, 1.0, 0.0]
        );
    }
}
