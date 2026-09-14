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
//! - Polygon / **MultiPolygon（(3g)・§11.8 反子午線分割）**: `GeoPolygon::geojson_geometry` は
//!   リングが ±180 を跨ぐ場合に断片へ分割して MultiPolygon を返す。跨がないリングの出力は不変。
//!   詳細は後半の「(3g) 反子午線分割」節の見出しコメント参照。
//!
//! ## テスト戦略（strict / mutation-resistant / FAST）
//! 全 FAST（実エンジン不要）。`serde_json::from_str` ではなく直接返る `Value` を構造で検証する。
//! lon と lat を別値にして座標取り違え変異を殺す。完全文字列一致は避け Value 構造で縛る。
//!
//! ## 期待される RED（実装前）
//! `GeoPoint::geojson_geometry` / `GeoLine::geojson_geometry` が未定義のためメソッド解決不能
//! （E0599）でコンパイルできない。これが想定どおりの赤。

use serde_json::Value;
use umbra_geo::{EnclosedPole, GeoLine, GeoPoint, GeoPolygon};

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

// （旧 `geo_polygon_geojson_antimeridian_stays_single_polygon` は (3g) §11.8 で仕様が
//   反転したため削除。跨ぎリングの出力は下の「反子午線分割」節で MultiPolygon として縛る。）

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

// ============================================================
// GeoPolygon::geojson_geometry — (3g) 反子午線分割（§11.8 確定仕様）
//
// ## 縛る仕様（§11.8）
// - (b) 跨ぎ判定は閉リングの連続 2 点の `Δlon = lon2 − lon1`:
//     `Δlon < −180` = 東進（+180 → −180）、`Δlon > +180` = 西進（−180 → +180）。
//     `|Δlon| = 180` ちょうどは**跨ぎでない**（測度ゼロ境界・GeoLine (3d) と同一規約）。
//   跨ぎ点の緯度は子午線上で線形補間（GeoLine と同一式・ヘルパ `lat_c_east` / `lat_c_west`）。
// - (b) 各半球の開いた弧は、子午線上を **lat_c の緯度順**で対にして閉じる。
// - (c) 断片 1 枚 → `{"type":"Polygon"}`（跨がないリングの出力は**現行と不変**）、
//        2 枚以上 → `{"type":"MultiPolygon","coordinates":[polygon, …]}`。
//        穴は分割後、それを含む（最小面積の）外環断片へ割り当てる。
// - 既存規約は維持: 閉リング（先頭==末尾）・外環 CCW / 穴 CW・点の捏造なし・空 rings は `coordinates: []`。
//
// ## 仕様が定めていない自由度（テストでは固定しない）
// - MultiPolygon 内の**多角形の並び順**、および断片リングの**開始頂点**（どの頂点から書き出すか）。
//   → `sorted_polys`（外環の (min lon, min lat) でソート）と `ring_cycle`（巡回正規化）で吸収し、
//     「座標列・向き・閉性」だけを厳密に縛る。
// ============================================================

/// 閉リング（先頭==末尾）を検証し、末尾複製を落として**巡回正規化**した開リングを返す。
/// 辞書順最小の頂点が先頭に来るよう回転するだけで、**向き（周回方向）は変えない**。
/// これにより「開始頂点は仕様未定義」を吸収しつつ、座標値・順序・向きは厳密に比較できる。
fn ring_cycle(coords: &[(f64, f64)]) -> Vec<(f64, f64)> {
    assert!(
        coords.len() >= 4,
        "面のあるリングは閉じて 4 点以上, got {coords:?}"
    );
    assert_eq!(
        coords.first(),
        coords.last(),
        "リングは閉じる（先頭==末尾）"
    );
    let open = &coords[..coords.len() - 1];
    let mut best = 0usize;
    for (i, p) in open.iter().enumerate() {
        let b = open[best];
        if (p.0, p.1) < (b.0, b.1) {
            best = i;
        }
    }
    let mut out = Vec::with_capacity(open.len());
    for k in 0..open.len() {
        out.push(open[(best + k) % open.len()]);
    }
    out
}

/// 展開済みリング座標を「閉性 + 巡回正規化後の座標列」で期待値（開リング・(lon,lat)）と厳密比較する。
/// 開始頂点の自由度だけを吸収し、頂点数・座標値・周回方向は厳密に縛る。
fn assert_ring_coords_eq(got: &[(f64, f64)], expected_open: &[(f64, f64)], what: &str) {
    let cyc = ring_cycle(got);
    let closed: Vec<(f64, f64)> = expected_open
        .iter()
        .copied()
        .chain(std::iter::once(expected_open[0]))
        .collect();
    let exp = ring_cycle(&closed);
    assert_eq!(
        cyc.len(),
        exp.len(),
        "{what}: 頂点数 {} expected {}（got={got:?}）",
        cyc.len(),
        exp.len()
    );
    for (i, (g, e)) in cyc.iter().zip(exp.iter()).enumerate() {
        assert!(
            close(g.0, e.0) && close(g.1, e.1),
            "{what}: 頂点{i} = [{}, {}] expected [{}, {}]（巡回正規化 got={cyc:?} expected={exp:?}）",
            g.0,
            g.1,
            e.0,
            e.1
        );
    }
}

/// MultiPolygon の `coordinates` を「多角形 → リング → (lon,lat) 列」に展開し、
/// **外環の (min lon, min lat) 昇順**で並べ替えて返す（多角形の並び順は仕様未定義のため）。
fn sorted_polys(g: &Value) -> Vec<Vec<Vec<(f64, f64)>>> {
    let polys = g["coordinates"].as_array().expect("coordinates は配列");
    let mut out: Vec<Vec<Vec<(f64, f64)>>> = polys
        .iter()
        .map(|p| {
            p.as_array()
                .expect("多角形はリング配列")
                .iter()
                .map(ring_coords)
                .collect()
        })
        .collect();
    fn key(p: &[Vec<(f64, f64)>]) -> (f64, f64) {
        let lon = p[0].iter().map(|&(l, _)| l).fold(f64::INFINITY, f64::min);
        let lat = p[0].iter().map(|&(_, a)| a).fold(f64::INFINITY, f64::min);
        (lon, lat)
    }
    out.sort_by(|a, b| key(a).partial_cmp(&key(b)).expect("有限値"));
    out
}

// ------------------------------------------------------------
// 1. 跨がないリングは現行どおり Polygon（分割は共通ケースで透明）
// ------------------------------------------------------------

/// **§11.8(c)「跨がないリングの出力はバイト不変」**（+180 側に近いが跨がない）。
/// lon=[170, 179, 174]（最大 |Δlon| = 9）・lat=[5, 8, 21]（CCW）→ `type="Polygon"`・
/// リング 1 本・閉じて 4 点・入力順そのまま・**±180 の頂点を 1 つも作らない**。
///
/// 殺す変異: 経度が大きければ（例 |lon| > 160）無条件に分割する・跨ぎ判定を `|Δlon| > 0` 等に緩める・
///   跨がないのに子午線頂点 ±180 を挿入する・常に MultiPolygon を返す。
#[test]
fn geo_polygon_antimeridian_no_crossing_near_plus180_stays_polygon() {
    // (lat, lon) 順で与える。(lon,lat) = (170,5),(179,8),(174,21)。
    let outer = ring(&[(5.0, 170.0), (8.0, 179.0), (21.0, 174.0)]);
    let g = GeoPolygon::new(vec![outer]).geojson_geometry();
    assert_eq!(
        g["type"],
        Value::String("Polygon".to_string()),
        "跨がない（|Δlon| ≤ 9）→ Polygon のまま"
    );
    let rings = g["coordinates"].as_array().expect("リング配列");
    assert_eq!(rings.len(), 1, "リング 1 本（分割しない）");
    let r0 = ring_coords(&rings[0]);
    assert_eq!(r0.len(), 4, "非閉 3 点 → 閉じて 4 点");
    assert_eq!(r0.first(), r0.last(), "閉リング");
    // 入力は CCW（面積 +66）なので反転せず入力順のまま出る。
    let expected = [(170.0, 5.0), (179.0, 8.0), (174.0, 21.0), (170.0, 5.0)];
    for (i, e) in expected.iter().enumerate() {
        assert!(
            close(r0[i].0, e.0) && close(r0[i].1, e.1),
            "点{i} = [{}, {}] expected [{}, {}]",
            r0[i].0,
            r0[i].1,
            e.0,
            e.1
        );
    }
    assert!(
        !r0.iter().any(|&(lon, _)| close(lon.abs(), 180.0)),
        "跨がないので ±180 子午線頂点を作らない, got {r0:?}"
    );
}

/// **§11.8(c)** の鏡像（−180 側に近いが跨がない）。lon=[−179, −170, −175]・lat=[−7, −4, 9]（CCW）。
/// 東側テストと対にして「符号だけで分割を決める」実装を殺す。
///
/// 殺す変異: `lon < 0` 側だけ／`lon > 0` 側だけ跨ぎ扱いにする・判定に Δlon でなく lon の絶対値を使う・
///   負経度で ±180 頂点を挿入する。
#[test]
fn geo_polygon_antimeridian_no_crossing_near_minus180_stays_polygon() {
    // (lon,lat) = (-179,-7),(-170,-4),(-175,9)。面積 +66（CCW）。
    let outer = ring(&[(-7.0, -179.0), (-4.0, -170.0), (9.0, -175.0)]);
    let g = GeoPolygon::new(vec![outer]).geojson_geometry();
    assert_eq!(
        g["type"],
        Value::String("Polygon".to_string()),
        "Polygon のまま"
    );
    let rings = g["coordinates"].as_array().expect("リング配列");
    assert_eq!(rings.len(), 1, "リング 1 本");
    let r0 = ring_coords(&rings[0]);
    let expected = [
        (-179.0, -7.0),
        (-170.0, -4.0),
        (-175.0, 9.0),
        (-179.0, -7.0),
    ];
    assert_eq!(r0.len(), expected.len(), "閉じて 4 点");
    for (i, e) in expected.iter().enumerate() {
        assert!(
            close(r0[i].0, e.0) && close(r0[i].1, e.1),
            "点{i} = [{}, {}] expected [{}, {}]",
            r0[i].0,
            r0[i].1,
            e.0,
            e.1
        );
    }
    assert!(
        !r0.iter().any(|&(lon, _)| close(lon.abs(), 180.0)),
        "±180 子午線頂点を作らない, got {r0:?}"
    );
}

/// **§11.8(b) 測度ゼロ境界**: `|Δlon| = 180` ちょうどは跨ぎ**でない**。
/// (lon,lat) = (−90,5),(90,5),(90,25)：Δlon は順に +180・0・−180（閉じ辺）。
/// `< −180` / `> +180` の厳密不等号ゆえどれも跨ぎでなく、出力は Polygon 1 リング（入力順・CCW）。
///
/// 殺す変異: 判定を `Δlon <= −180` / `>= 180`（オフバイワン）にして 180 ちょうどを分割する・
///   `|Δlon| >= 180` にまとめる・閉じ辺（末尾→先頭）を跨ぎ判定から漏らす／誤検出する。
#[test]
fn geo_polygon_antimeridian_delta_exactly_180_is_not_crossing() {
    // (lat, lon)。(lon,lat) = (-90,5),(90,5),(90,25)。面積 +1800（CCW）。
    let outer = ring(&[(5.0, -90.0), (5.0, 90.0), (25.0, 90.0)]);
    let g = GeoPolygon::new(vec![outer]).geojson_geometry();
    assert_eq!(
        g["type"],
        Value::String("Polygon".to_string()),
        "|Δlon| = 180 ちょうどは跨ぎでない → Polygon"
    );
    let rings = g["coordinates"].as_array().expect("リング配列");
    assert_eq!(rings.len(), 1, "リング 1 本（分割しない）");
    let r0 = ring_coords(&rings[0]);
    let expected = [(-90.0, 5.0), (90.0, 5.0), (90.0, 25.0), (-90.0, 5.0)];
    assert_eq!(r0.len(), expected.len(), "閉じて 4 点");
    for (i, e) in expected.iter().enumerate() {
        assert!(
            close(r0[i].0, e.0) && close(r0[i].1, e.1),
            "点{i} = [{}, {}] expected [{}, {}]",
            r0[i].0,
            r0[i].1,
            e.0,
            e.1
        );
    }
    assert!(
        !r0.iter().any(|&(lon, _)| close(lon.abs(), 180.0)),
        "子午線頂点を挿入しない, got {r0:?}"
    );
}

// ------------------------------------------------------------
// 2. 東進跨ぎ（外環のみ）
// ------------------------------------------------------------

/// **§11.8(b)(c) 東進跨ぎ** の基本形。(lon,lat) の矩形
/// `(165,10) → (−172,10) → (−172,40) → (165,40)`:
/// - 辺 (165,10)→(−172,10): Δlon = −337 < −180 → **東進**。lat は両端 10 ゆえ lat_c = 10。
/// - 辺 (−172,40)→(165,40): Δlon = +337 > 180 → **西進**。lat_c = 40。
///
/// 期待（2 枚・東西で経度幅を 15 / 8 と**非対称**にし、左右取り違えを検出）:
/// - 東断片（+180 子午線で閉じる）: `[(180,40),(165,40),(165,10),(180,10)]`（CCW・面積 +450）
/// - 西断片（−180 子午線で閉じる）: `[(-180,10),(-172,10),(-172,40),(-180,40)]`（CCW・面積 +240）
///
/// 殺す変異: 跨ぎを検出せず 1 枚の Polygon のまま出す・東進で子午線を −180/+180 逆に付ける
///   （東側断片に −180 が現れる）・子午線頂点を挿入せず弧を開いたまま出す（閉リングでない）・
///   断片を CCW に正規化しない・経度を 0..360 に付け替えて誤魔化す・断片数を 1 や 3 にする。
#[test]
fn geo_polygon_antimeridian_east_crossing_splits_into_two_polygons() {
    let outer = ring(&[(10.0, 165.0), (10.0, -172.0), (40.0, -172.0), (40.0, 165.0)]);
    let g = GeoPolygon::new(vec![outer]).geojson_geometry();
    assert_eq!(
        g["type"],
        Value::String("MultiPolygon".to_string()),
        "東進跨ぎ → MultiPolygon"
    );
    let ps = sorted_polys(&g);
    assert_eq!(ps.len(), 2, "断片は 2 枚, got {}", ps.len());
    // ソート順は外環の min lon 昇順 → 西（−180 始まり）・東（165 始まり）。
    let west = &ps[0];
    let east = &ps[1];
    assert_eq!(west.len(), 1, "西断片は外環のみ（穴なし）");
    assert_eq!(east.len(), 1, "東断片は外環のみ（穴なし）");

    // 交点緯度は水平辺ゆえ端点と同値（独立オラクルでも確認）。
    assert!(
        close(lat_c_east(165.0, 10.0, -172.0, 10.0), 10.0),
        "東進 lat_c = 10"
    );
    assert!(
        close(lat_c_west(-172.0, 40.0, 165.0, 40.0), 40.0),
        "西進 lat_c = 40"
    );

    // 東断片: +180 子午線で閉じる。西断片: −180 子午線で閉じる（符号の取り違えを殺す）。
    assert_ring_coords_eq(
        &east[0],
        &[(180.0, 40.0), (165.0, 40.0), (165.0, 10.0), (180.0, 10.0)],
        "東断片の外環",
    );
    assert_ring_coords_eq(
        &west[0],
        &[
            (-180.0, 10.0),
            (-172.0, 10.0),
            (-172.0, 40.0),
            (-180.0, 40.0),
        ],
        "西断片の外環",
    );
    // 外環はいずれも CCW（面積 +450 / +240・独立 shoelace オラクル）。
    assert!(
        close(signed_area(&east[0]), 450.0),
        "東断片は CCW・面積 +450, got {}",
        signed_area(&east[0])
    );
    assert!(
        close(signed_area(&west[0]), 240.0),
        "西断片は CCW・面積 +240, got {}",
        signed_area(&west[0])
    );
    // 東断片に −180 は現れず、西断片に +180 は現れない（子午線の取り違え検出）。
    assert!(
        !east[0].iter().any(|&(lon, _)| close(lon, -180.0)),
        "東断片に −180 が混入, got {:?}",
        east[0]
    );
    assert!(
        !west[0].iter().any(|&(lon, _)| close(lon, 180.0)),
        "西断片に +180 が混入, got {:?}",
        west[0]
    );
}

// ------------------------------------------------------------
// 3. 西進跨ぎ（東進の鏡像・向きを逆に辿る）
// ------------------------------------------------------------

/// **§11.8(b) 西進跨ぎ**（`Δlon > +180`）。テスト 2 と**逆向き**に辿る矩形
/// `(−158,12) → (173,12) → (173,44) → (−158,44)`:
/// - 辺 (−158,12)→(173,12): Δlon = +331 > 180 → **西進**（−180 → +180）。lat_c = 12。
/// - 辺 (173,44)→(−158,44): Δlon = −331 < −180 → **東進**。lat_c = 44。
///
/// この向きだと各断片は (lon,lat) 平面で **CW** に出来上がるので、外環 CCW 正規化が必ず働く
/// （テスト 2 は正規化なしで CCW になる形＝両方で正規化の有無を分けて縛る）。
/// 緯度（12/44）・経度幅（東 7 / 西 22）をテスト 2 と別値かつ非対称にして、東西の取り違えを可視化する。
///
/// 期待:
/// - 東断片: `[(180,44),(173,44),(173,12),(180,12)]`（CCW・面積 7·32 = +224）
/// - 西断片: `[(-180,12),(-158,12),(-158,44),(-180,44)]`（CCW・面積 22·32 = +704）
///
/// 殺す変異: 西進（`Δlon > 180`）の判定を落として分割しない・西進でも東進と同じ子午線符号を使う
///   （西断片が +180 で閉じる）・西進 t の分子を `(180 − lon1)` にする・
///   逆向き入力の断片を CW のまま出す（外環 CCW 正規化の欠落）。
#[test]
fn geo_polygon_antimeridian_west_crossing_splits_and_normalizes_ccw() {
    let outer = ring(&[(12.0, -158.0), (12.0, 173.0), (44.0, 173.0), (44.0, -158.0)]);
    let g = GeoPolygon::new(vec![outer]).geojson_geometry();
    assert_eq!(
        g["type"],
        Value::String("MultiPolygon".to_string()),
        "西進跨ぎ → MultiPolygon"
    );
    let ps = sorted_polys(&g);
    assert_eq!(ps.len(), 2, "断片は 2 枚, got {}", ps.len());
    let west = &ps[0];
    let east = &ps[1];
    assert_eq!(west.len(), 1, "西断片は外環のみ");
    assert_eq!(east.len(), 1, "東断片は外環のみ");

    // 独立オラクル（水平辺ゆえ端点と同値）。
    assert!(
        close(lat_c_west(-158.0, 12.0, 173.0, 12.0), 12.0),
        "西進 lat_c = 12"
    );
    assert!(
        close(lat_c_east(173.0, 44.0, -158.0, 44.0), 44.0),
        "東進 lat_c = 44"
    );

    assert_ring_coords_eq(
        &east[0],
        &[(180.0, 44.0), (173.0, 44.0), (173.0, 12.0), (180.0, 12.0)],
        "東断片の外環",
    );
    assert_ring_coords_eq(
        &west[0],
        &[
            (-180.0, 12.0),
            (-158.0, 12.0),
            (-158.0, 44.0),
            (-180.0, 44.0),
        ],
        "西断片の外環",
    );
    assert!(
        close(signed_area(&east[0]), 224.0),
        "東断片 CCW・面積 +224（幅 7 × 高 32）, got {}",
        signed_area(&east[0])
    );
    assert!(
        close(signed_area(&west[0]), 704.0),
        "西断片 CCW・面積 +704（幅 22 × 高 32）, got {}",
        signed_area(&west[0])
    );
}

// ------------------------------------------------------------
// 4. 斜め辺＝真の線形補間
// ------------------------------------------------------------

/// **§11.8(b) 補間式 `lat_c = lat1 + t·(lat2 − lat1)`** を、跨ぎ辺が**水平でない**形で縛る。
/// リング（(lon,lat)）: `(168,10) → (−176,42) → (−176,50) → (168,60)`。
/// - 東進辺 (168,10)→(−176,42): Δlon = −344。t = (180−168)/(360−344) = 12/16 = **0.75**、
///   lat_c = 10 + 0.75·32 = **34**（両端 10・42 の**どちらとも異なる**真の内分点）。
/// - 西進辺 (−176,50)→(168,60): Δlon = +344。t = (−176+180)/(360−344) = 4/16 = **0.25**、
///   lat_c = 50 + 0.25·10 = **52.5**（こちらも両端と異なる）。
///   値はすべて 2 進で厳密（12/16・4/16 は 2 の冪分母）。
///
/// 期待:
/// - 東断片: `[(180,52.5),(168,60),(168,10),(180,34)]`（CCW・面積 +411）
/// - 西断片: `[(-180,34),(-176,42),(-176,50),(-180,52.5)]`（CCW・面積 +53）
///
/// 殺す変異: lat ではなく lon を補間する（子午線点が ±180 でなくなる）・t を 0.5 固定/端点代用にする
///   （lat_c が 34/52.5 でなく 10・42・50・60 のどれかになる）・東進と西進で t の式を取り違える
///   （34 ↔ 52.5 が入れ替わる）・分母を 360 固定にする・lat1 と lat2 を取り違える。
#[test]
fn geo_polygon_antimeridian_slanted_edges_interpolate_latitude() {
    let outer = ring(&[(10.0, 168.0), (42.0, -176.0), (50.0, -176.0), (60.0, 168.0)]);
    let g = GeoPolygon::new(vec![outer]).geojson_geometry();
    assert_eq!(
        g["type"],
        Value::String("MultiPolygon".to_string()),
        "斜め辺の跨ぎ → MultiPolygon"
    );

    // 独立オラクル（実装式を写経せずヘルパで再計算）。
    let lat_a = lat_c_east(168.0, 10.0, -176.0, 42.0);
    let lat_b = lat_c_west(-176.0, 50.0, 168.0, 60.0);
    assert!(close(lat_a, 34.0), "東進 lat_c = 34（t=0.75）, got {lat_a}");
    assert!(
        close(lat_b, 52.5),
        "西進 lat_c = 52.5（t=0.25）, got {lat_b}"
    );

    let ps = sorted_polys(&g);
    assert_eq!(ps.len(), 2, "断片は 2 枚, got {}", ps.len());
    let west = &ps[0];
    let east = &ps[1];
    assert_eq!(west.len(), 1, "西断片は外環のみ");
    assert_eq!(east.len(), 1, "東断片は外環のみ");

    assert_ring_coords_eq(
        &east[0],
        &[(180.0, lat_b), (168.0, 60.0), (168.0, 10.0), (180.0, lat_a)],
        "東断片の外環（斜め辺の補間点を含む）",
    );
    assert_ring_coords_eq(
        &west[0],
        &[
            (-180.0, lat_a),
            (-176.0, 42.0),
            (-176.0, 50.0),
            (-180.0, lat_b),
        ],
        "西断片の外環（斜め辺の補間点を含む）",
    );
    assert!(
        close(signed_area(&east[0]), 411.0),
        "東断片 CCW・面積 +411, got {}",
        signed_area(&east[0])
    );
    assert!(
        close(signed_area(&west[0]), 53.0),
        "西断片 CCW・面積 +53, got {}",
        signed_area(&west[0])
    );
    // 補間緯度は跨ぎ辺の端点緯度のいずれとも異なる（端点代用の変異を明示的に殺す）。
    for bad in [10.0_f64, 42.0, 50.0, 60.0] {
        assert!(!close(lat_a, bad), "lat_c=34 は端点 {bad} と別値のはず");
        assert!(!close(lat_b, bad), "lat_c=52.5 は端点 {bad} と別値のはず");
    }
}

// ------------------------------------------------------------
// 5. 外環と穴の**両方**が跨ぐ
// ------------------------------------------------------------

/// **§11.8(c) 穴の割当**: 外環も穴も跨ぐ場合、穴も同じ手順で分割し、各断片の外環へ
/// 「それを含むもの」に割り当てる。穴は CW・外環は CCW。
///
/// 入力（(lon,lat)）:
/// - 外環 `(160,5),(−160,5),(−160,55),(160,55)` → 東断片 lon 160..180 / 西断片 lon −180..−160、lat 5..55。
/// - 穴   `(168,18),(−176,18),(−176,44),(168,44)` → 東断片 lon 168..180 / 西断片 lon −180..−176、lat 18..44。
///   穴は東西で幅 12 / 4 と**非対称**なので、穴断片の左右取り違えが座標で露見する。
///
/// 期待: MultiPolygon 2 枚、各 2 リング（外環 + 穴 1 本）、穴は合計 2 本（消失も重複もなし）:
/// - 東: 外環 `[(180,55),(160,55),(160,5),(180,5)]`（CCW）／穴 `[(180,18),(168,18),(168,44),(180,44)]`（CW）
/// - 西: 外環 `[(-180,5),(-160,5),(-160,55),(-180,55)]`（CCW）／穴 `[(-180,44),(-176,44),(-176,18),(-180,18)]`（CW）
///
/// 殺す変異: 穴を分割せず跨いだまま片方（または両方）の断片に付ける・穴を両断片に重複して付ける・
///   穴を落とす（リング 1 本だけの断片になる）・穴を東西逆の断片に割り当てる・穴を CCW のまま出す・
///   穴を独立した外環（別の多角形）として出す（断片数が 4 になる）。
#[test]
fn geo_polygon_antimeridian_outer_and_hole_both_cross() {
    let outer = ring(&[(5.0, 160.0), (5.0, -160.0), (55.0, -160.0), (55.0, 160.0)]);
    let hole = ring(&[(18.0, 168.0), (18.0, -176.0), (44.0, -176.0), (44.0, 168.0)]);
    let g = GeoPolygon::new(vec![outer, hole]).geojson_geometry();
    assert_eq!(
        g["type"],
        Value::String("MultiPolygon".to_string()),
        "外環が跨ぐ → MultiPolygon"
    );
    let ps = sorted_polys(&g);
    assert_eq!(
        ps.len(),
        2,
        "断片は 2 枚（穴を別多角形にしない）, got {}",
        ps.len()
    );
    let west = &ps[0];
    let east = &ps[1];
    assert_eq!(west.len(), 2, "西断片は外環 + 穴 1 本, got {}", west.len());
    assert_eq!(east.len(), 2, "東断片は外環 + 穴 1 本, got {}", east.len());
    let total_holes: usize = ps.iter().map(|p| p.len() - 1).sum();
    assert_eq!(
        total_holes, 2,
        "穴は全体で 2 本（消失も重複もなし）, got {total_holes}"
    );

    assert_ring_coords_eq(
        &east[0],
        &[(180.0, 55.0), (160.0, 55.0), (160.0, 5.0), (180.0, 5.0)],
        "東断片の外環",
    );
    assert_ring_coords_eq(
        &east[1],
        &[(180.0, 18.0), (168.0, 18.0), (168.0, 44.0), (180.0, 44.0)],
        "東断片の穴",
    );
    assert_ring_coords_eq(
        &west[0],
        &[(-180.0, 5.0), (-160.0, 5.0), (-160.0, 55.0), (-180.0, 55.0)],
        "西断片の外環",
    );
    assert_ring_coords_eq(
        &west[1],
        &[
            (-180.0, 44.0),
            (-176.0, 44.0),
            (-176.0, 18.0),
            (-180.0, 18.0),
        ],
        "西断片の穴",
    );

    // 向き: 外環 CCW（面積>0）・穴 CW（面積<0）。面積値も厳密に縛る。
    assert!(
        close(signed_area(&east[0]), 1000.0),
        "東外環 CCW・面積 +1000（20×50）, got {}",
        signed_area(&east[0])
    );
    assert!(
        close(signed_area(&east[1]), -312.0),
        "東の穴は CW・面積 −312（12×26）, got {}",
        signed_area(&east[1])
    );
    assert!(
        close(signed_area(&west[0]), 1000.0),
        "西外環 CCW・面積 +1000, got {}",
        signed_area(&west[0])
    );
    assert!(
        close(signed_area(&west[1]), -104.0),
        "西の穴は CW・面積 −104（4×26）, got {}",
        signed_area(&west[1])
    );
}

// ------------------------------------------------------------
// 6. 外環だけが跨ぎ、穴は片半球に収まる
// ------------------------------------------------------------

/// **§11.8(c) 穴の割当（含む断片のみ）**: 外環だけが跨ぎ、穴は東半球に完全に収まる場合、
/// 穴は東断片にだけ付き、西断片は外環 1 本だけになる。
///
/// 入力: 外環はテスト 5 と同じ（lon 160 → −160・lat 5..55）。
/// 穴は `(166,14),(176,14),(176,30),(166,30)`（跨がない・東半球・CCW 入力 → 出力は CW 正規化）。
///
/// 期待: 2 枚。東 = 外環 + 穴（穴の座標は入力そのまま・向きだけ CW へ反転）、西 = 外環のみ（リング 1 本）。
///
/// 殺す変異: 穴を全断片に配る（西断片にもリング 2 本が出る）・穴を「最初の断片」へ固定で付ける
///   （並べ替え後の西に付く）・跨がない穴を無条件に分割して ±180 頂点を生やす・穴を落とす・
///   穴を CCW のまま出す。
#[test]
fn geo_polygon_antimeridian_hole_in_one_hemisphere_only() {
    let outer = ring(&[(5.0, 160.0), (5.0, -160.0), (55.0, -160.0), (55.0, 160.0)]);
    let hole = ring(&[(14.0, 166.0), (14.0, 176.0), (30.0, 176.0), (30.0, 166.0)]);
    let g = GeoPolygon::new(vec![outer, hole]).geojson_geometry();
    assert_eq!(
        g["type"],
        Value::String("MultiPolygon".to_string()),
        "MultiPolygon"
    );
    let ps = sorted_polys(&g);
    assert_eq!(ps.len(), 2, "断片は 2 枚, got {}", ps.len());
    let west = &ps[0];
    let east = &ps[1];
    assert_eq!(
        west.len(),
        1,
        "西断片は外環のみ（穴は含まれない）, got {} リング",
        west.len()
    );
    assert_eq!(
        east.len(),
        2,
        "東断片は外環 + 穴, got {} リング",
        east.len()
    );

    assert_ring_coords_eq(
        &west[0],
        &[(-180.0, 5.0), (-160.0, 5.0), (-160.0, 55.0), (-180.0, 55.0)],
        "西断片の外環",
    );
    assert_ring_coords_eq(
        &east[0],
        &[(180.0, 55.0), (160.0, 55.0), (160.0, 5.0), (180.0, 5.0)],
        "東断片の外環",
    );
    // 穴は跨がないので座標は入力のまま（±180 を生やさない）。向きだけ CW に正規化。
    assert_ring_coords_eq(
        &east[1],
        &[(166.0, 30.0), (176.0, 30.0), (176.0, 14.0), (166.0, 14.0)],
        "東断片の穴（跨がない・CW 正規化のみ）",
    );
    assert!(
        !east[1].iter().any(|&(lon, _)| close(lon.abs(), 180.0)),
        "跨がない穴に ±180 頂点を生やさない, got {:?}",
        east[1]
    );
    assert!(
        close(signed_area(&east[1]), -160.0),
        "穴は CW・面積 −160（10×16）, got {}",
        signed_area(&east[1])
    );
}

// ------------------------------------------------------------
// 7. 4 回跨ぎ（子午線上の対は緯度順）
// ------------------------------------------------------------

/// **§11.8(b)「子午線上を緯度順に対で結ぶ」**を、跨ぎ 4 回の櫛形で縛る。
/// リング（(lon,lat)・CCW）:
/// `(160,3) → (−152,3) → (−152,12) → (170,12) → (170,25) → (−152,25) → (−152,38) → (160,38)`
/// Δlon は順に −312（東進）・0・+322（西進）・0・−322（東進）・0・+312（西進）・0 ＝ **跨ぎ 4 回**。
/// 交点緯度は 3・12・25・38（水平辺ゆえ端点と同値）。
///
/// 形は「東側の帯（lon 160..180・lat 3..38）から、西側へ 2 本の指（lat 3..12 と lat 25..38）が伸びる」。
/// 東側は lat 12..25 で lon 170 まで凹む**1 つの連結領域**、西側は**2 つの矩形**。→ **断片は 3 枚**
/// （跨ぎ回数 4 / 2 = 2 ではない＝「交差数の半分＝断片数」と決め打つ実装を殺す）。
///
/// 期待:
/// - 東: `[(160,3),(180,3),(180,12),(170,12),(170,25),(180,25),(180,38),(160,38)]`（CCW・面積 20·35 − 10·13 = +570）
/// - 西下: `[(-180,3),(-152,3),(-152,12),(-180,12)]`（CCW・面積 28·9 = +252）
/// - 西上: `[(-180,25),(-152,25),(-152,38),(-180,38)]`（CCW・面積 28·13 = +364）
///
/// 殺す変異: 子午線上の点を**緯度順でなく出現順／最近傍**で対にする（西側が lat 12〜25 を跨いで
///   1 枚に繋がり断片 2 枚・面積 28·35 になる）・最初の跨ぎだけ処理して残りを落とす・
///   断片数を「跨ぎ数 / 2」と決め打つ・東側の凹み（lon 170 の指の壁）を子午線で埋めて矩形にする。
#[test]
fn geo_polygon_antimeridian_four_crossings_pairs_by_latitude() {
    let outer = ring(&[
        (3.0, 160.0),
        (3.0, -152.0),
        (12.0, -152.0),
        (12.0, 170.0),
        (25.0, 170.0),
        (25.0, -152.0),
        (38.0, -152.0),
        (38.0, 160.0),
    ]);
    let g = GeoPolygon::new(vec![outer]).geojson_geometry();
    assert_eq!(
        g["type"],
        Value::String("MultiPolygon".to_string()),
        "4 回跨ぎ → MultiPolygon"
    );
    let ps = sorted_polys(&g);
    assert_eq!(
        ps.len(),
        3,
        "断片は 3 枚（東 1 + 西 2）。跨ぎ 4 回でも 2 枚ではない, got {}",
        ps.len()
    );
    // ソートキー = 外環の (min lon, min lat)。西下(−180,3) < 西上(−180,25) < 東(160,3)。
    let west_low = &ps[0];
    let west_high = &ps[1];
    let east = &ps[2];
    for (i, p) in ps.iter().enumerate() {
        assert_eq!(p.len(), 1, "断片{i} は外環のみ（穴なし）, got {}", p.len());
    }

    assert_ring_coords_eq(
        &west_low[0],
        &[(-180.0, 3.0), (-152.0, 3.0), (-152.0, 12.0), (-180.0, 12.0)],
        "西の下側断片（子午線 3↔12 で閉じる）",
    );
    assert_ring_coords_eq(
        &west_high[0],
        &[
            (-180.0, 25.0),
            (-152.0, 25.0),
            (-152.0, 38.0),
            (-180.0, 38.0),
        ],
        "西の上側断片（子午線 25↔38 で閉じる）",
    );
    assert_ring_coords_eq(
        &east[0],
        &[
            (160.0, 3.0),
            (180.0, 3.0),
            (180.0, 12.0),
            (170.0, 12.0),
            (170.0, 25.0),
            (180.0, 25.0),
            (180.0, 38.0),
            (160.0, 38.0),
        ],
        "東断片（子午線 3↔12・25↔38 で閉じ、lat 12..25 は lon 170 の指の壁）",
    );

    assert!(
        close(signed_area(&west_low[0]), 252.0),
        "西下 CCW・面積 +252, got {}",
        signed_area(&west_low[0])
    );
    assert!(
        close(signed_area(&west_high[0]), 364.0),
        "西上 CCW・面積 +364, got {}",
        signed_area(&west_high[0])
    );
    assert!(
        close(signed_area(&east[0]), 570.0),
        "東 CCW・面積 +570（凹み 10×13 を差し引いた値）, got {}",
        signed_area(&east[0])
    );
}

// ------------------------------------------------------------
// 8. 退化ガード
// ------------------------------------------------------------

/// **空 rings は分割導入後も `coordinates: []`**（`type` は Polygon のまま・リングを捏造しない）。
/// 既存の `geo_polygon_geojson_empty_rings_is_empty_coordinates` の (3g) 後の再確認。
///
/// 殺す変異: 空入力で MultiPolygon（`coordinates: [[]]` 等）を返す・panic する・
///   跨ぎ判定でリング先頭要素を無条件参照して index out of bounds。
#[test]
fn geo_polygon_antimeridian_empty_rings_still_empty_coordinates() {
    let g = GeoPolygon::new(Vec::new()).geojson_geometry();
    assert_eq!(
        g["type"],
        Value::String("Polygon".to_string()),
        "空 rings は Polygon のまま（MultiPolygon にしない）"
    );
    let rings = g["coordinates"].as_array().expect("coordinates は配列");
    assert!(
        rings.is_empty(),
        "空 rings → coordinates 長 0, got {}",
        rings.len()
    );
}

/// **点を捏造しない**（§11.7「umbra-geo 側」規約の継続）: 2 点の退化リングが ±180 を跨いでも、
/// 面積のある多角形を作ってはならない。
/// 入力 `(170,10) → (−170,20)`（Δlon = −340・東進）は面積 0 の線分。
///
/// 仕様 §11.8 は退化リングの出力形（Polygon のまま / 分割して面積 0 の断片）を定めていないため、
/// ここでは**どちらでも通る**が「出力されたどのリングも符号付き面積 0」を縛る
/// （＝子午線の結線で面積を捏造したら落ちる）。
///
/// 殺す変異: 子午線上の点で退化リングを閉じて 2 つの三角形/矩形（面積≠0）を作る・
///   panic する・±180 まで引き延ばした面積のある多角形を返す。
#[test]
fn geo_polygon_antimeridian_two_point_degenerate_ring_fabricates_no_area() {
    let outer = ring(&[(10.0, 170.0), (20.0, -170.0)]);
    let g = GeoPolygon::new(vec![outer]).geojson_geometry();
    let ty = g["type"].as_str().expect("type は文字列");
    assert!(
        ty == "Polygon" || ty == "MultiPolygon",
        "退化でも Polygon/MultiPolygon のいずれか（panic せず）, got {ty}"
    );
    // すべてのリングを平坦に集めて面積 0 を確認する。
    let mut rings: Vec<Vec<(f64, f64)>> = Vec::new();
    if ty == "Polygon" {
        for r in g["coordinates"].as_array().expect("リング配列") {
            rings.push(ring_coords(r));
        }
    } else {
        for p in g["coordinates"].as_array().expect("多角形配列") {
            for r in p.as_array().expect("リング配列") {
                rings.push(ring_coords(r));
            }
        }
    }
    for (i, r) in rings.iter().enumerate() {
        let a = signed_area(r);
        assert!(
            a.abs() < EPS,
            "退化 2 点リング由来のリング{i} が面積を持つ（捏造）: area={a}, ring={r:?}"
        );
    }
}

// ------------------------------------------------------------
// 9. 極を囲むリング（跨ぎ回数が**奇数**）は分割しない（§11.8(d)「未対応」の明示的な契約）
//
// ## 事実（実データ実測・2026-09-14）
// 閉リングの経度は一周して戻るので、跨ぎ（|Δlon| > 180 の辺）を ±360 の補正として数えると
// 「跨ぎ回数が奇数」⇔「経度の巻き数が ±1」⇔ **リングが地球を一周している＝極を囲む**。
// 実日食での実測: 2016-03-09 の部分食域は **2 回**跨ぎ（真の跨ぎ・極を囲まない）だが、
// 2021-12-04（**1 回**）・2028-07-22（**3 回**）・2037-07-13（**3 回**）は南極を囲む南極日食で、
// いずれも奇数回跨ぐ。
//
// ## 縛る仕様（§11.8(d)）
// 極を囲む領域は **未対応**。よって §11.8(b) の「跨ぎは偶数回」という前提が崩れるこれらの入力では、
// 分割を**行わず**入力リングをそのまま（閉じて・従来どおり環向き正規化して）単一 `Polygon` で返す。
// これは「対応していないものを、黙って壊さずに素通しする」契約であり、
// 弧の消失・面積の捏造・panic のいずれも起こしてはならない。
// ------------------------------------------------------------

/// **§11.8(d)**: 南極を囲む（＝経度が単調東進で一周する）リングは **跨ぎ 1 回（奇数）** で、
/// 分割**せず** 単一 `Polygon` として原座標のまま返す。
///
/// 入力（(lon,lat)）: `(−150,−75) → (−60,−78) → (30,−75) → (120,−78) → (170,−72)`（閉リング）。
/// 各辺の Δlon は `+90, +90, +90, +50`、**閉じ辺** (170,−72)→(−150,−75) が `Δlon = −320 < −180`
/// ＝東進跨ぎ。よって跨ぎ回数は **1 回（奇数）**＝経度の巻き数 +1＝地球を一周＝南極を囲む。
/// （緯度を −72..−78 で振って平面 shoelace を非ゼロ +885 にし、退化リングと区別する。）
///
/// 期待: `type="Polygon"`・リング 1 本・入力 5 点 + 先頭複製 = 6 点・座標は入力そのまま
/// （shoelace +885 > 0 ＝ 既に CCW なので反転もしない）・**±180 の頂点を 1 つも生やさない**。
///
/// 殺す実装: 奇数回跨ぎで弧のペアリングに失敗し、**閉じられなかった弧を黙って捨てる**
///   （リングが消える／頂点が欠ける＝面積が失われる）・**未対応の弧で panic する**
///   （unwrap / index out of bounds）・片方の半球の弧だけを拾って
///   **1 枚だけの偽フラグメントを MultiPolygon で返す**（面積の捏造・領域の半分の消失）・
///   極側を子午線で無理に閉じて存在しない面積を作る。
#[test]
fn geo_polygon_antimeridian_pole_enclosing_odd_crossing_is_not_split() {
    // (lat, lon) で与える。(lon,lat) = (-150,-75),(-60,-78),(30,-75),(120,-78),(170,-72)。
    let outer = ring(&[
        (-75.0, -150.0),
        (-78.0, -60.0),
        (-75.0, 30.0),
        (-78.0, 120.0),
        (-72.0, 170.0),
    ]);

    // 独立オラクル(1): 跨ぎ回数は 1（奇数）＝極を囲む。閉じ辺を含めて数える。
    let lons = [-150.0_f64, -60.0, 30.0, 120.0, 170.0];
    let crossings = lons
        .iter()
        .enumerate()
        .filter(|(i, &l1)| {
            let l2 = lons[(i + 1) % lons.len()];
            (l2 - l1).abs() > 180.0
        })
        .count();
    assert_eq!(crossings, 1, "前提: 跨ぎは 1 回（奇数）＝極を囲む");

    // 独立オラクル(2): 入力閉リングの平面 shoelace は +885（CCW・非退化）→ 反転されない。
    let input_closed = [
        (-150.0, -75.0),
        (-60.0, -78.0),
        (30.0, -75.0),
        (120.0, -78.0),
        (170.0, -72.0),
        (-150.0, -75.0),
    ];
    assert!(
        close(signed_area(&input_closed), 885.0),
        "前提: 入力の shoelace は +885（CCW・面積非ゼロ）, got {}",
        signed_area(&input_closed)
    );

    let g = GeoPolygon::new(vec![outer]).geojson_geometry();
    assert_eq!(
        g["type"],
        Value::String("Polygon".to_string()),
        "極を囲む（跨ぎ奇数回）リングは分割しない → Polygon（§11.8(d) 未対応の素通し）"
    );
    let rings = g["coordinates"].as_array().expect("リング配列");
    assert_eq!(
        rings.len(),
        1,
        "リング 1 本（断片化しない）, got {}",
        rings.len()
    );
    let r0 = ring_coords(&rings[0]);
    assert_eq!(
        r0.len(),
        6,
        "入力 5 点 + 先頭複製 = 6 点（頂点の欠落なし）, got {}",
        r0.len()
    );
    assert_eq!(r0.first(), r0.last(), "閉リング（先頭==末尾）");
    assert_ring_coords_eq(
        &r0,
        &[
            (-150.0, -75.0),
            (-60.0, -78.0),
            (30.0, -75.0),
            (120.0, -78.0),
            (170.0, -72.0),
        ],
        "極を囲むリング（原座標保存・CCW のまま）",
    );
    assert!(
        !r0.iter().any(|&(lon, _)| close(lon.abs(), 180.0)),
        "分割しないので ±180 子午線頂点を生やさない, got {r0:?}"
    );
    // 面積は入力どおり（捏造も消失もしない）。
    assert!(
        close(signed_area(&r0), 885.0),
        "出力面積は入力どおり +885（捏造も消失もなし）, got {}",
        signed_area(&r0)
    );
}

/// **§11.8(d)** の 2 例目: 跨ぎ **3 回（奇数）** でも同様に分割しない。
/// 「1 回だけを特別扱いする」実装（奇数一般でなく `crossings == 1` を見る）を殺す。
///
/// 入力（(lon,lat)）:
/// `(−150,−75) → (−60,−78) → (30,−75) → (120,−78) → (170,−72) → (−175,−70) → (170,−68)`（閉リング）。
/// Δlon は順に `+90, +90, +90, +50, −345(東進跨ぎ①), +345(西進跨ぎ②)`、
/// **閉じ辺** (170,−68)→(−150,−75) が `−320`（東進跨ぎ③）＝**跨ぎ 3 回（奇数）**。
/// 跨ぎを ±360 で補正した経度の総変位は `+360`（巻き数 +1）＝やはり南極を囲む
/// （跨ぎ①②は打ち消し合い、③が正味の一周を与える）。
///
/// 期待: `type="Polygon"`・リング 1 本・入力 7 点 + 先頭複製 = 8 点・原座標保存
/// （shoelace +835 > 0 ＝ CCW なので反転なし）・±180 頂点なし。
///
/// 殺す実装: 「偶数回なら分割・それ以外は 1 回だけ素通し」と場当たりに書く・
///   打ち消し合う 2 回だけを処理して残り 1 回の弧を捨てる（頂点欠落＝面積損失）・
///   3 本の弧を無理に子午線で閉じて偽の断片（MultiPolygon）を返す・unpaired 弧で panic する。
#[test]
fn geo_polygon_antimeridian_pole_enclosing_three_crossings_is_not_split() {
    let outer = ring(&[
        (-75.0, -150.0),
        (-78.0, -60.0),
        (-75.0, 30.0),
        (-78.0, 120.0),
        (-72.0, 170.0),
        (-70.0, -175.0),
        (-68.0, 170.0),
    ]);

    // 独立オラクル(1): 跨ぎ 3 回（奇数）かつ ±360 補正後の総変位 +360（巻き数 +1＝極を囲む）。
    let lons = [-150.0_f64, -60.0, 30.0, 120.0, 170.0, -175.0, 170.0];
    let mut crossings = 0usize;
    let mut winding = 0.0_f64;
    for (i, &l1) in lons.iter().enumerate() {
        let l2 = lons[(i + 1) % lons.len()];
        let d = l2 - l1;
        if d < -180.0 {
            crossings += 1;
            winding += d + 360.0;
        } else if d > 180.0 {
            crossings += 1;
            winding += d - 360.0;
        } else {
            winding += d;
        }
    }
    assert_eq!(crossings, 3, "前提: 跨ぎは 3 回（奇数）");
    assert!(
        close(winding, 360.0),
        "前提: 経度の巻き数は +1（総変位 +360）＝極を囲む, got {winding}"
    );

    // 独立オラクル(2): 入力閉リングの平面 shoelace は +835（CCW・非退化）。
    let input_closed = [
        (-150.0, -75.0),
        (-60.0, -78.0),
        (30.0, -75.0),
        (120.0, -78.0),
        (170.0, -72.0),
        (-175.0, -70.0),
        (170.0, -68.0),
        (-150.0, -75.0),
    ];
    assert!(
        close(signed_area(&input_closed), 835.0),
        "前提: 入力の shoelace は +835（CCW）, got {}",
        signed_area(&input_closed)
    );

    let g = GeoPolygon::new(vec![outer]).geojson_geometry();
    assert_eq!(
        g["type"],
        Value::String("Polygon".to_string()),
        "跨ぎ 3 回（奇数・極を囲む）も分割しない → Polygon"
    );
    let rings = g["coordinates"].as_array().expect("リング配列");
    assert_eq!(rings.len(), 1, "リング 1 本, got {}", rings.len());
    let r0 = ring_coords(&rings[0]);
    assert_eq!(
        r0.len(),
        8,
        "入力 7 点 + 先頭複製 = 8 点（頂点の欠落なし）, got {}",
        r0.len()
    );
    assert_eq!(r0.first(), r0.last(), "閉リング");
    assert_ring_coords_eq(
        &r0,
        &[
            (-150.0, -75.0),
            (-60.0, -78.0),
            (30.0, -75.0),
            (120.0, -78.0),
            (170.0, -72.0),
            (-175.0, -70.0),
            (170.0, -68.0),
        ],
        "極を囲むリング（跨ぎ 3 回・原座標保存）",
    );
    assert!(
        !r0.iter().any(|&(lon, _)| close(lon.abs(), 180.0)),
        "±180 子午線頂点を生やさない, got {r0:?}"
    );
    assert!(
        close(signed_area(&r0), 835.0),
        "出力面積は入力どおり +835, got {}",
        signed_area(&r0)
    );
}

// 対比ガード（跨ぎ **偶数**回＝極を囲まないリングは分割する）は新規追加しない。
// 既存テストが同じことを縛っている:
//   - `geo_polygon_antimeridian_east_crossing_splits_into_two_polygons`（跨ぎ 2 回・東進 → MultiPolygon 2 枚）
//   - `geo_polygon_antimeridian_west_crossing_splits_and_normalizes_ccw`（跨ぎ 2 回・西進 → MultiPolygon 2 枚）
//   - `geo_polygon_antimeridian_four_crossings_pairs_by_latitude`（跨ぎ 4 回 → MultiPolygon 3 枚）

// ============================================================
// 10. mutation 工程（(3g) §11.8）生存変異の判別テスト
//   分割の**ガード条件・穴の割当規則・退行リングの扱い**を、出力で弁別できる配置で縛る。
// ============================================================

/// **§11.8(b)(c)**: 頂点 **ちょうど 3 個**の三角形でも、反子午線を 2 回跨ぐなら**分割する**。
///
/// 撃つ変異: 分割対象の最小頂点数ガード `if ring.len() < 3 { 分割しない }` → `<= 3`
///   （3 頂点リングを**跨いだまま単一 Polygon で素通し**し、地球を一周する不正な多角形を出す）。
///
/// 入力（(lon,lat)）: `(170,0) → (−170,10) → (−170,−10)`（3 頂点）。
/// - 辺 (170,0)→(−170,10): Δlon = −340 < −180 → **東進**。t = (180−170)/(360−340) = 10/20 = **0.5**
///   （2 進で厳密）、lat_c = 0 + 0.5·10 = **+5**。
/// - 辺 (−170,10)→(−170,−10): Δlon = 0（跨ぎなし）。
/// - 閉じ辺 (−170,−10)→(170,0): Δlon = +340 > 180 → **西進**。t = (−170+180)/(360−340) = **0.5**、
///   lat_c = −10 + 0.5·(0−(−10)) = **−5**（同じく厳密）。
///   跨ぎは **2 回（偶数）**＝極を囲まない（連続経度の巻き数 0）ので §11.8(d) の素通しにも当たらない。
///
/// 期待（2 枚・面積を東 50 / 西 150 と**非対称**にして左右取り違えも検出）:
/// - 東断片: `[(180,5),(170,0),(180,−5)]`（CCW・面積 +50）
/// - 西断片: `[(-180,−5),(-170,−10),(-170,10),(-180,5)]`（CCW・面積 +150）
#[test]
fn geo_polygon_antimeridian_three_vertex_triangle_is_split() {
    // (lat, lon) で与える。(lon,lat) = (170,0),(−170,10),(−170,−10)。
    let outer = ring(&[(0.0, 170.0), (10.0, -170.0), (-10.0, -170.0)]);
    assert_eq!(outer.len(), 3, "前提: 入力リングはちょうど 3 頂点");

    // 独立オラクル: 跨ぎは 2 回（偶数）・巻き数 0（極を囲まない）。
    let lons = [170.0_f64, -170.0, -170.0];
    let mut crossings = 0usize;
    let mut winding = 0.0_f64;
    for (i, &l1) in lons.iter().enumerate() {
        let d = lons[(i + 1) % lons.len()] - l1;
        if d < -180.0 {
            crossings += 1;
            winding += d + 360.0;
        } else if d > 180.0 {
            crossings += 1;
            winding += d - 360.0;
        } else {
            winding += d;
        }
    }
    assert_eq!(crossings, 2, "前提: 跨ぎは 2 回（偶数）");
    assert!(close(winding, 0.0), "前提: 巻き数 0（極を囲まない）");

    // 独立オラクル: 交点緯度は ±5（両端 0/10・0/−10 のいずれとも異なる真の内分点）。
    let lat_e = lat_c_east(170.0, 0.0, -170.0, 10.0);
    let lat_w = lat_c_west(-170.0, -10.0, 170.0, 0.0);
    assert!(close(lat_e, 5.0), "東進 lat_c = +5（t=0.5）, got {lat_e}");
    assert!(close(lat_w, -5.0), "西進 lat_c = −5（t=0.5）, got {lat_w}");

    let g = GeoPolygon::new(vec![outer]).geojson_geometry();
    assert_eq!(
        g["type"],
        Value::String("MultiPolygon".to_string()),
        "3 頂点でも 2 回跨ぐなら分割する → MultiPolygon（変異 `<= 3` は Polygon のまま素通しする）"
    );
    let ps = sorted_polys(&g);
    assert_eq!(ps.len(), 2, "断片は 2 枚, got {}", ps.len());
    let west = &ps[0];
    let east = &ps[1];
    assert_eq!(west.len(), 1, "西断片は外環のみ");
    assert_eq!(east.len(), 1, "東断片は外環のみ");

    assert_ring_coords_eq(
        &east[0],
        &[(180.0, lat_e), (170.0, 0.0), (180.0, lat_w)],
        "東断片の外環（+180 子午線で lat_c = +5 ↔ −5 を結ぶ）",
    );
    assert_ring_coords_eq(
        &west[0],
        &[
            (-180.0, lat_w),
            (-170.0, -10.0),
            (-170.0, 10.0),
            (-180.0, lat_e),
        ],
        "西断片の外環（−180 子午線で閉じる）",
    );
    assert!(
        close(signed_area(&east[0]), 50.0),
        "東断片 CCW・面積 +50, got {}",
        signed_area(&east[0])
    );
    assert!(
        close(signed_area(&west[0]), 150.0),
        "西断片 CCW・面積 +150, got {}",
        signed_area(&west[0])
    );
}

/// **§11.8(b) 弧の対応付けが成立しない入力では分割を「全部やらない」**（all-or-nothing）。
///
/// 撃つ変異: フォールバック条件 `if out.is_empty() || consumed != arc_count { 元のリングを返す }`
///   → `&&`（＝**閉じられた分だけの部分的な MultiPolygon を返し、対を成さなかった弧の領域を黙って失う**）。
///
/// 入力は**自己交差リング**（意図的）。跨ぎ 6 回（偶数＝§11.8(d) の極ガードは発火しない）で、
/// 東側の 3 弧のうち 2 弧は §11.8(b) の緯度規則で閉環を作れるが、残り 1 弧は結線先が無く**余る**。
/// 西側も同様に 2 弧が閉じ 1 弧が余る。よって `out` は非空・`consumed(4) != arc_count(6)` となり、
/// original は**分割を取り下げて入力リングをそのまま単一 Polygon で返す**。
/// 変異 `&&` は「閉じた 4 弧ぶんの断片」だけを返し、余った 2 弧の領域が消える（面積の損失）。
///
/// 入力（(lon,lat)・閉リング）:
/// `(−175,40),(−175,70),(160,70),(160,50),(−160,50),(−160,20),(170,20),(170,30),(−170,30),(−170,10),(165,10),(165,40)`
/// 跨ぎ辺はすべて水平なので交点緯度は端点緯度に一致し、東側の弧は
/// `(10→40)`・`(20→30)`・`(70→50)`、西側は `(40→70)`・`(50→20)`・`(30→10)` という
/// **入れ子でも素でもない（互いに交差する）緯度区間**になる。
///
/// 期待: `type="Polygon"`・リング 1 本・入力 12 点 + 先頭複製 = 13 点・座標は入力そのまま
/// （入力の平面 shoelace は +6400 ＝ CCW なので反転もしない）・**±180 頂点を 1 つも作らない**・
/// 面積は入力どおり +6400（消失も捏造もなし）。
#[test]
fn geo_polygon_antimeridian_unpairable_arcs_decline_split_without_area_loss() {
    // (lat, lon) で与える。
    let outer = ring(&[
        (40.0, -175.0),
        (70.0, -175.0),
        (70.0, 160.0),
        (50.0, 160.0),
        (50.0, -160.0),
        (20.0, -160.0),
        (20.0, 170.0),
        (30.0, 170.0),
        (30.0, -170.0),
        (10.0, -170.0),
        (10.0, 165.0),
        (40.0, 165.0),
    ]);

    // 独立オラクル(1): 跨ぎは 6 回（**偶数**＝極ガードではなく対応付け失敗で素通しすることを保証）・巻き数 0。
    let lons = [
        -175.0_f64, -175.0, 160.0, 160.0, -160.0, -160.0, 170.0, 170.0, -170.0, -170.0, 165.0,
        165.0,
    ];
    let mut crossings = 0usize;
    let mut winding = 0.0_f64;
    for (i, &l1) in lons.iter().enumerate() {
        let d = lons[(i + 1) % lons.len()] - l1;
        if d < -180.0 {
            crossings += 1;
            winding += d + 360.0;
        } else if d > 180.0 {
            crossings += 1;
            winding += d - 360.0;
        } else {
            winding += d;
        }
    }
    assert_eq!(crossings, 6, "前提: 跨ぎは 6 回（偶数）, got {crossings}");
    assert!(
        close(winding, 0.0),
        "前提: 巻き数 0（極を囲まない＝§11.8(d) の素通しには当たらない）, got {winding}"
    );

    // 独立オラクル(2): 入力閉リングの平面 shoelace は +6400（CCW・非退化）→ 反転されない。
    let input_open = [
        (-175.0, 40.0),
        (-175.0, 70.0),
        (160.0, 70.0),
        (160.0, 50.0),
        (-160.0, 50.0),
        (-160.0, 20.0),
        (170.0, 20.0),
        (170.0, 30.0),
        (-170.0, 30.0),
        (-170.0, 10.0),
        (165.0, 10.0),
        (165.0, 40.0),
    ];
    assert!(
        close(signed_area(&input_open), 6400.0),
        "前提: 入力の shoelace は +6400（CCW）, got {}",
        signed_area(&input_open)
    );

    let g = GeoPolygon::new(vec![outer]).geojson_geometry();
    assert_eq!(
        g["type"],
        Value::String("Polygon".to_string()),
        "弧が対を成さない入力は分割を取り下げて単一 Polygon（変異 `&&` は部分的な MultiPolygon を返す）"
    );
    let rings = g["coordinates"].as_array().expect("リング配列");
    assert_eq!(rings.len(), 1, "リング 1 本, got {}", rings.len());
    let r0 = ring_coords(&rings[0]);
    assert_eq!(
        r0.len(),
        13,
        "入力 12 点 + 先頭複製 = 13 点（頂点の欠落なし）, got {}",
        r0.len()
    );
    assert_eq!(r0.first(), r0.last(), "閉リング（先頭==末尾）");
    assert_ring_coords_eq(&r0, &input_open, "対応付け失敗で素通しされたリング");
    assert!(
        !r0.iter().any(|&(lon, _)| close(lon.abs(), 180.0)),
        "分割しないので ±180 子午線頂点を捏造しない, got {r0:?}"
    );
    assert!(
        close(signed_area(&r0), 6400.0),
        "出力面積は入力どおり +6400（部分的な分割による面積損失なし）, got {}",
        signed_area(&r0)
    );
}

// （穴プローブ頂点 `|lon| < 180` → `> 180` の判別テストは**追加しない**: 入力頂点は
//   `EastLongitude::from_degrees` が [−180, 180) に正規化するため |lon| > 180 の頂点は作れず、
//   ちょうど −180 でも `abs() > 180` は偽＝変異下ではプローブが**常に不成立で穴が全部捨てられる**。
//   これは既存の `..._outer_and_hole_both_cross`（穴 2 本）・`..._hole_in_one_hemisphere_only`
//   （東断片が外環+穴の 2 リング）が既に落とすので、新規テストは冗長。)

/// **§11.8(c) 手順 3「穴は、それを含む外環断片のうち*最小面積*のものへ割り当てる」**。
///
/// 撃つ変異: `signed_area2(o).abs() < signed_area2(outers[j]).abs()` → `>` / `==`
///   （`>` は**最大面積の断片**へ付ける＝穴が外側の大きな断片に開いて、内側の小さな断片が塞がったままになる。
///   `==` は比較が成立せず割当が破綻して穴を落とす/誤配する）。
///
/// 入力（自己交差する外環・意図的）（(lon,lat)）:
/// `(−170,50),(120,50),(120,10),(150,10),(150,40),(−170,40),(−170,20),(125,20),(125,10),(−170,10)`
/// 跨ぎは 4 回（偶数）。東側の 2 弧は緯度区間 `[40,50]` と `[10,20]` で**互いに素**なので
/// それぞれ単独で閉じ、**同じ東半球に 2 枚の外環断片**ができる:
/// - F1（大・面積 **+1500**）: `[(180,50),(120,50),(120,10),(150,10),(150,40),(180,40)]`（フック形）
/// - F2（小・面積 **+550**）: `[(180,20),(125,20),(125,10),(180,10)]`（矩形）
///   F2 の領域（lon 125..180・lat 10..20）は F1 の脚（lon 120..150・lat 10..40）と重なり、
///   重なり（lon 125..150・lat 10..20）に穴 `lon 130..145・lat 13..17`（面積 60）を置く。
///   この穴は **F1 と F2 の両方に含まれる**ので、「最小面積の断片」規則が観測可能になる。
///
/// 期待: 穴は**面積 +550 の F2 にだけ**付き（CW・面積 −60・座標そのまま）、
/// 面積 +1500 の F1 は**穴なし**（リング 1 本）。西断片も穴なし。穴は全体で 1 本。
///
/// 注: `<=` への変異は含有断片の面積が相異なる限り original と同値（tie でのみ差が出るが、
/// tie の解決順は仕様未定義）。本テストが撃つのは `>` と `==`。
#[test]
fn geo_polygon_antimeridian_hole_assigned_to_smallest_containing_fragment() {
    // (lat, lon) で与える。
    let outer = ring(&[
        (50.0, -170.0),
        (50.0, 120.0),
        (10.0, 120.0),
        (10.0, 150.0),
        (40.0, 150.0),
        (40.0, -170.0),
        (20.0, -170.0),
        (20.0, 125.0),
        (10.0, 125.0),
        (10.0, -170.0),
    ]);
    // 穴: lon 130..145・lat 13..17（入力で既に CW・面積 −60）。F1・F2 の**両方**の内部にある。
    let hole = ring(&[(13.0, 130.0), (17.0, 130.0), (17.0, 145.0), (13.0, 145.0)]);

    // 独立オラクル: 跨ぎ 4 回（偶数）・巻き数 0。
    let lons = [
        -170.0_f64, 120.0, 120.0, 150.0, 150.0, -170.0, -170.0, 125.0, 125.0, -170.0,
    ];
    let mut crossings = 0usize;
    let mut winding = 0.0_f64;
    for (i, &l1) in lons.iter().enumerate() {
        let d = lons[(i + 1) % lons.len()] - l1;
        if d < -180.0 {
            crossings += 1;
            winding += d + 360.0;
        } else if d > 180.0 {
            crossings += 1;
            winding += d - 360.0;
        } else {
            winding += d;
        }
    }
    assert_eq!(crossings, 4, "前提: 跨ぎは 4 回（偶数）, got {crossings}");
    assert!(close(winding, 0.0), "前提: 巻き数 0, got {winding}");

    // 独立オラクル: 期待断片の面積（F1 = 1500 > F2 = 550）と穴（−60）。
    let f1 = [
        (180.0, 50.0),
        (120.0, 50.0),
        (120.0, 10.0),
        (150.0, 10.0),
        (150.0, 40.0),
        (180.0, 40.0),
    ];
    let f2 = [(180.0, 20.0), (125.0, 20.0), (125.0, 10.0), (180.0, 10.0)];
    let hole_open = [(130.0, 13.0), (130.0, 17.0), (145.0, 17.0), (145.0, 13.0)];
    assert!(
        close(signed_area(&f1), 1500.0),
        "前提: F1 は CCW・面積 +1500, got {}",
        signed_area(&f1)
    );
    assert!(
        close(signed_area(&f2), 550.0),
        "前提: F2 は CCW・面積 +550（F1 より小さい）, got {}",
        signed_area(&f2)
    );
    assert!(
        close(signed_area(&hole_open), -60.0),
        "前提: 穴は CW・面積 −60, got {}",
        signed_area(&hole_open)
    );

    let g = GeoPolygon::new(vec![outer, hole]).geojson_geometry();
    assert_eq!(
        g["type"],
        Value::String("MultiPolygon".to_string()),
        "跨ぎ 4 回 → MultiPolygon"
    );
    let ps = sorted_polys(&g);
    assert_eq!(
        ps.len(),
        3,
        "断片は 3 枚（東 2 枚〔F1・F2〕+ 西 1 枚）, got {}",
        ps.len()
    );
    // ソートキー =(min lon, min lat): 西(−180) < F1(120) < F2(125)。
    let west = &ps[0];
    let big = &ps[1];
    let small = &ps[2];

    assert_ring_coords_eq(&big[0], &f1, "東の大断片 F1 の外環");
    assert_ring_coords_eq(&small[0], &f2, "東の小断片 F2 の外環");
    assert!(
        close(signed_area(&big[0]), 1500.0),
        "F1 の外環は面積 +1500, got {}",
        signed_area(&big[0])
    );
    assert!(
        close(signed_area(&small[0]), 550.0),
        "F2 の外環は面積 +550, got {}",
        signed_area(&small[0])
    );

    let total_holes: usize = ps.iter().map(|p| p.len() - 1).sum();
    assert_eq!(total_holes, 1, "穴は全体で 1 本, got {total_holes}");
    assert_eq!(
        big.len(),
        1,
        "面積 +1500 の F1 は**穴なし**（変異 `>` は最大面積側に穴を付ける）, got {} リング",
        big.len()
    );
    assert_eq!(
        west.len(),
        1,
        "西断片は外環のみ（穴を含まない）, got {} リング",
        west.len()
    );
    assert_eq!(
        small.len(),
        2,
        "面積 +550 の F2（最小の含有断片）が穴を持つ, got {} リング",
        small.len()
    );
    assert_ring_coords_eq(&small[1], &hole_open, "F2 の穴（座標そのまま・CW）");
    assert!(
        close(signed_area(&small[1]), -60.0),
        "穴は CW・面積 −60, got {}",
        signed_area(&small[1])
    );
}

/// **§11.8(d) 退行リング**: **1 点だけ**のリングは「捏造せずそのまま出す」＝点を**消してはならない**。
///
/// 撃つ変異: 閉じ重複頂点の除去条件 `coords.len() > 1 && coords[0] == coords[last]` → `>= 1`
///   （1 点リングでは `coords[0] == coords[0]` が真になるので唯一の点が**閉じ重複と誤認されて除去**され、
///   **空リング**（`coordinates: [[]]`）になる＝入力点の消失）。
///
/// 既存の退行契約に整合（`GeoLine` 0/1 点 → 長さ 0/1・空 rings → `coordinates: []`・
/// 2 点リング → 面積を捏造しない）: 1 点リングは**点を保持**し、その座標は入力の `[lon, lat]` だけ。
/// 閉じ処理で先頭複製が付くか否か（長さ 1 か 2 か）は §11.8 が定めていないので**どちらでも通す**が、
/// 「非空」「現れる座標はすべて入力点」「点を捏造しない（長さ ≤ 2）」を厳密に縛る。
#[test]
fn geo_polygon_geojson_single_point_ring_keeps_its_point() {
    // 1 点だけのリング。lat=7・lon=8 の非対称値で [lat,lon] 逆順も検出する。
    let g = GeoPolygon::new(vec![ring(&[(7.0, 8.0)])]).geojson_geometry();
    assert_eq!(
        g["type"],
        Value::String("Polygon".to_string()),
        "1 点リングでも type=Polygon（分割対象ではない）"
    );
    let rings = g["coordinates"].as_array().expect("coordinates は配列");
    assert_eq!(rings.len(), 1, "リング 1 本, got {}", rings.len());
    let r0 = ring_coords(&rings[0]);
    assert!(
        !r0.is_empty(),
        "1 点リングの点を消さない（変異 `>= 1` は唯一の点を閉じ重複と誤認して空リングにする）, got {r0:?}"
    );
    assert!(
        r0.len() <= 2,
        "点を捏造しない（そのまま、または閉じ重複 1 点まで）, got {r0:?}"
    );
    for (i, &(lon, lat)) in r0.iter().enumerate() {
        assert!(
            close(lon, 8.0) && close(lat, 7.0),
            "座標{i} = [8, 7]（[lon,lat] 順）, got [{lon}, {lat}]"
        );
    }
}

/// **§11.8(b)(c) 入力リング表現の不変性**: 反子午線を跨ぐリングを**閉表現**（先頭頂点を末尾に複製）
/// で与えても**開表現**（複製なし）で与えても、`geojson_geometry()` の出力は**完全に同一**。
/// 外環・穴の両方を 2 表現で与え、分割経路（跨ぎ有り）を通した上で同値性を縛る。
///
/// 実装は分割前に「末尾の閉じ重複頂点を剥がして開表現にする」正規化を行う。本テストはその契約を縛る。
///
/// 撃つ変異（誤実装）:
///   - 閉じ重複の剥がしを**やめる**（条件を常に偽にする）: 閉入力では末尾に `lon == 先頭 lon` の
///     重複頂点が残り、退化辺や余分な子午線交点・頂点数差を生んで開入力と異なる出力になる。
///   - 閉じ重複を**無条件に剥がす**（先頭==末尾の判定を落として常に末尾を捨てる）: 開入力で
///     最終頂点（(162,50) / (170,40)）が失われ、断片座標が期待値と食い違う。
///     さらに「両表現が同じように壊れる」変異に備え、具体構造（MultiPolygon・断片 2 枚・
///     外環座標と面積・穴の総数）も厳密に縛る。
///
/// 入力（(lon,lat)）— 経度はすべて (−180, 180) の内側（`from_degrees` は +180 を −180 に丸めるため）:
/// - 外環 `(162,8),(−164,8),(−164,50),(162,50)` → 東断片 lon 162..180（幅 18）/ 西断片 lon −180..−164（幅 16）、lat 8..50（高 42）。
/// - 穴   `(170,20),(−172,20),(−172,40),(170,40)` → 東断片 lon 170..180 / 西断片 lon −180..−172。
#[test]
fn geojson_is_invariant_to_closed_or_open_input_rings() {
    // (lat, lon) 順で与える開表現。
    let outer_open = ring(&[(8.0, 162.0), (8.0, -164.0), (50.0, -164.0), (50.0, 162.0)]);
    let hole_open = ring(&[(20.0, 170.0), (20.0, -172.0), (40.0, -172.0), (40.0, 170.0)]);
    // 閉表現（先頭頂点を末尾に複製）。
    let mut outer_closed = outer_open.clone();
    outer_closed.push(outer_open[0]);
    let mut hole_closed = hole_open.clone();
    hole_closed.push(hole_open[0]);

    let g_open = GeoPolygon::new(vec![outer_open, hole_open]).geojson_geometry();
    let g_closed = GeoPolygon::new(vec![outer_closed, hole_closed]).geojson_geometry();

    // (1) 2 表現の出力は完全一致（Value の厳密比較）。
    assert_eq!(
        g_open, g_closed,
        "閉表現と開表現で出力が異なる（閉じ重複頂点の剥がし漏れ／過剰剥がし）"
    );

    // (2) 具体構造（両表現が同じように壊れる変異を殺す）。
    for (label, g) in [("開表現", &g_open), ("閉表現", &g_closed)] {
        assert_eq!(
            g["type"],
            Value::String("MultiPolygon".to_string()),
            "{label}: 跨ぎリング → MultiPolygon"
        );
        let ps = sorted_polys(g);
        assert_eq!(ps.len(), 2, "{label}: 断片は 2 枚, got {}", ps.len());
        let west = &ps[0];
        let east = &ps[1];
        assert_eq!(west.len(), 2, "{label}: 西断片は外環 + 穴 1 本");
        assert_eq!(east.len(), 2, "{label}: 東断片は外環 + 穴 1 本");
        let total_holes: usize = ps.iter().map(|p| p.len() - 1).sum();
        assert_eq!(
            total_holes, 2,
            "{label}: 穴は全体で 2 本, got {total_holes}"
        );

        // 外環断片の座標を厳密比較（閉じ重複が残ると頂点数・座標が食い違う）。
        assert_ring_coords_eq(
            &east[0],
            &[(180.0, 50.0), (162.0, 50.0), (162.0, 8.0), (180.0, 8.0)],
            &format!("{label}: 東断片の外環"),
        );
        assert_ring_coords_eq(
            &west[0],
            &[(-180.0, 8.0), (-164.0, 8.0), (-164.0, 50.0), (-180.0, 50.0)],
            &format!("{label}: 西断片の外環"),
        );
        // 外環は CCW・面積は 18×42 = 756 / 16×42 = 672（独立 shoelace オラクル）。
        assert!(
            close(signed_area(&east[0]), 756.0),
            "{label}: 東断片 CCW・面積 +756, got {}",
            signed_area(&east[0])
        );
        assert!(
            close(signed_area(&west[0]), 672.0),
            "{label}: 西断片 CCW・面積 +672, got {}",
            signed_area(&west[0])
        );
        // 穴は CW（面積<0）。
        assert!(
            signed_area(&east[1]) < 0.0,
            "{label}: 東断片の穴は CW, got {}",
            signed_area(&east[1])
        );
        assert!(
            signed_area(&west[1]) < 0.0,
            "{label}: 西断片の穴は CW, got {}",
            signed_area(&west[1])
        );
    }
}

// ============================================================
// 11. ISSUE-051 極で閉じる（`geojson_geometry_with_pole`・§確定仕様 1〜4）
//
// ## 縛る仕様（ISSUE-051 §確定仕様）
// - (1) 外環の反子午線跨ぎ回数が**奇数**なら極を囲む。
// - (2) どちらの極かは **umbra-geo の外**（umbra-eclipse）が決める。geo 層は引数で受け取る。
//       `None` は「決められなかった」であり、**(3g) の素通し（現状維持）**にフォールバックする。
// - (3) 閉じ方: 弧の端点から子午線 `+180` / `−180` を極（`lat = ±90`）まで辿り、極上で渡り、戻る。
//       **挿入する頂点は既存の補間交点 `(±180, lat_c)` と極頂点 `(±180, ±90)` のみ**。
// - (4) 跨がない／偶数跨ぎの出力は極引数に**一切影響されない**（バイト不変）。
// - `geojson_geometry()` ≡ `geojson_geometry_with_pole(None)`。
//
// ## 期待される RED（実装前）
// `EnclosedPole` が未定義（E0432）・`geojson_geometry_with_pole` 未定義（E0599）でコンパイル不能。
//
// ## 仕様が定めていない自由度（テストでは固定しない）
// - 断片リングの開始頂点（`assert_ring_coords_eq` の巡回正規化で吸収）。
// - 極頂点が既存頂点と**完全一致**する退化入力での重複頂点の除去有無（テスト 7 参照）。
// ============================================================

/// 「地球を逆走する辺」の検出オラクル。連続 2 頂点の |Δlon| の最大値を返すが、
/// **両端がともに極（|lat| = 90）の辺は除外**する（極上の辺は球面上で長さ 0 の退化辺であり、
/// 経度差 360 を持つのが正しい＝§確定仕様 3 が明示的に導入する辺）。
/// 極で閉じられた正しいリングではこの値が 180 以下になり、
/// 今日の不正な出力（弧の末尾から先頭へ逆走する辺）では 180 を超える。
fn max_non_polar_lon_span(coords: &[(f64, f64)]) -> f64 {
    let mut worst = 0.0_f64;
    for w in coords.windows(2) {
        let (lon1, lat1) = w[0];
        let (lon2, lat2) = w[1];
        if close(lat1.abs(), 90.0) && close(lat2.abs(), 90.0) {
            continue; // 極上の渡り辺（球面では 1 点）。
        }
        worst = worst.max((lon2 - lon1).abs());
    }
    worst
}

/// リングの全頂点が `allowed` のいずれかと一致することを縛る（**点を捏造しない**・§確定仕様 3）。
fn assert_all_vertices_allowed(coords: &[(f64, f64)], allowed: &[(f64, f64)], what: &str) {
    for (i, &(lon, lat)) in coords.iter().enumerate() {
        assert!(
            allowed.iter().any(|&(l, a)| close(lon, l) && close(lat, a)),
            "{what}: 頂点{i} = [{lon}, {lat}] は許された点集合に無い（捏造）。allowed={allowed:?}"
        );
    }
}

/// ISSUE-051 テストの共通 fixture（跨ぎ **1 回**・極を囲む）。(lat, lon) で返す。
///
/// (lon,lat) = `(−150,−75) → (−60,−78) → (30,−75) → (120,−78) → (170,−72)`。
/// 辺の Δlon は `+90, +90, +90, +50`、**閉じ辺** (170,−72)→(−150,−75) が `Δlon = −320 < −180`
/// ＝東進跨ぎ **1 回（奇数）**。経度は単調に東進して一周するので巻き数 +1 ＝極を囲む。
/// 交点緯度は `lat_c_east(170,−72,−150,−75)`: t = (180−170)/(360−320) = 10/40 = **0.25**（2 進で厳密）、
/// `lat_c = −72 + 0.25·(−3) = −72.75`（両端 −72・−75 のどちらとも異なる真の内分点）。
fn pole_ring_one_crossing() -> Vec<GeoPoint> {
    ring(&[
        (-75.0, -150.0),
        (-78.0, -60.0),
        (-75.0, 30.0),
        (-78.0, 120.0),
        (-72.0, 170.0),
    ])
}

/// 閉リングの跨ぎ回数と（±360 補正後の）経度総変位を**独立に**数える（前提の明示）。
fn crossings_and_winding(lons: &[f64]) -> (usize, f64) {
    let mut crossings = 0usize;
    let mut winding = 0.0_f64;
    for (i, &l1) in lons.iter().enumerate() {
        let d = lons[(i + 1) % lons.len()] - l1;
        if d < -180.0 {
            crossings += 1;
            winding += d + 360.0;
        } else if d > 180.0 {
            crossings += 1;
            winding += d - 360.0;
        } else {
            winding += d;
        }
    }
    (crossings, winding)
}

/// **ISSUE-051 §確定仕様 1・3（南極で閉じる）**: 跨ぎ **1 回（奇数）**のリングを `Some(South)` で
/// 閉じると、弧は `+180` 子午線を `lat = −90` まで下り、極上で `−180` へ渡り、戻って閉じる。
///
/// fixture は `pole_ring_one_crossing`（交点緯度 −72.75）。弧は `(−180,−72.75)` から始まり
/// 入力 5 点を経て `(+180,−72.75)` で終わる（断片は 1 枚）。§確定仕様 3 の閉じ方で
/// `(180,−90)`・`(−180,−90)` の 2 頂点**だけ**が加わる。
///
/// 期待（開表現・CCW）:
/// `[(−180,−90),(180,−90),(180,−72.75),(170,−72),(120,−78),(30,−75),(−60,−78),(−150,−75),(−180,−72.75)]`
/// 内部は**南極冠**（帯より南）なので、(lon,lat) 平面ではこの向きが CCW（shoelace **+5055**・手計算）。
/// 逆向き（東進のまま辿る順）は −5055 ＝ CW であり RFC 7946 の外環規約に反する。
/// 断片は 1 枚なので §11.8(c) どおり `type="Polygon"`。
///
/// 殺す実装: 極引数を無視して (3g) の素通しを返す（**経度幅 320 の逆走辺が残る**＝今日の欠陥）・
///   極頂点を入れずに弧の両端を直結する・極頂点を `(0,−90)` 等の捏造経度で入れる・
///   `(±180, lat_c)` の補間点を落として弧を極へ直結する・±180 の符号を取り違えて
///   東側の弧末尾に `−180` を付ける・リングを閉じない・CW のまま出す・北極（+90）で閉じる。
#[test]
fn geojson_with_pole_south_closes_odd_crossing_ring_through_south_pole() {
    // 前提（独立オラクル）: 跨ぎ 1 回（奇数）・巻き数 +1。
    let lons = [-150.0_f64, -60.0, 30.0, 120.0, 170.0];
    let (crossings, winding) = crossings_and_winding(&lons);
    assert_eq!(crossings, 1, "前提: 跨ぎは 1 回（奇数）＝極を囲む");
    assert!(close(winding, 360.0), "前提: 巻き数 +1, got {winding}");
    // 前提（独立オラクル）: 交点緯度 −72.75。
    let lat_c = lat_c_east(170.0, -72.0, -150.0, -75.0);
    assert!(close(lat_c, -72.75), "前提: lat_c = −72.75, got {lat_c}");

    let g = GeoPolygon::new(vec![pole_ring_one_crossing()])
        .geojson_geometry_with_pole(Some(EnclosedPole::South));

    assert_eq!(
        g["type"],
        Value::String("Polygon".to_string()),
        "極で閉じた結果は断片 1 枚 → Polygon（§11.8(c)）"
    );
    let rings = g["coordinates"].as_array().expect("リング配列");
    assert_eq!(
        rings.len(),
        1,
        "外環 1 本のみ（穴なし）, got {}",
        rings.len()
    );
    let r0 = ring_coords(&rings[0]);

    // (1) 閉リング・頂点数（入力 5 + 交点 2 + 極 2 = 9、+ 先頭複製 = 10）。
    assert_eq!(r0.first(), r0.last(), "閉リング（先頭==末尾）");
    assert_eq!(
        r0.len(),
        10,
        "入力 5 + ±180 交点 2 + 極頂点 2 = 9 点 + 先頭複製, got {}: {r0:?}",
        r0.len()
    );

    // (2) 座標列・周回方向を厳密比較（開始頂点の自由度のみ吸収）。
    let expected = [
        (-180.0, -90.0),
        (180.0, -90.0),
        (180.0, lat_c),
        (170.0, -72.0),
        (120.0, -78.0),
        (30.0, -75.0),
        (-60.0, -78.0),
        (-150.0, -75.0),
        (-180.0, lat_c),
    ];
    assert_ring_coords_eq(&r0, &expected, "南極で閉じた外環");

    // (3) 環向き: CCW（面積 +5055・独立 shoelace オラクル）。逆向きなら −5055。
    assert!(
        close(signed_area(&r0), 5055.0),
        "外環は CCW・面積 +5055（内部＝南極冠）, got {}",
        signed_area(&r0)
    );

    // (4) 極頂点が両子午線に 1 つずつある。北極頂点は現れない（極の取り違え検出）。
    assert!(
        r0.iter()
            .any(|&(lon, lat)| close(lon, 180.0) && close(lat, -90.0)),
        "(+180, −90) の極頂点が無い, got {r0:?}"
    );
    assert!(
        r0.iter()
            .any(|&(lon, lat)| close(lon, -180.0) && close(lat, -90.0)),
        "(−180, −90) の極頂点が無い, got {r0:?}"
    );
    assert!(
        !r0.iter().any(|&(_, lat)| close(lat, 90.0)),
        "南極で閉じたのに lat=+90 の頂点がある, got {r0:?}"
    );

    // (5) **今日の欠陥の直撃**: 極上の渡り辺を除き、経度幅 180 を超える辺が無い。
    //     素通し出力では閉じ辺 (170,−72)→(−150,−75) が 320 になる。
    let span = max_non_polar_lon_span(&r0);
    assert!(
        span <= 180.0 + EPS,
        "極以外に経度幅 >180 の辺がある（地球を逆走する偽の辺）: max={span}, ring={r0:?}"
    );

    // (6) 捏造なし: 全頂点は「入力頂点」「(±180, lat_c)」「(±180, −90)」のいずれか。
    let allowed = [
        (-150.0, -75.0),
        (-60.0, -78.0),
        (30.0, -75.0),
        (120.0, -78.0),
        (170.0, -72.0),
        (180.0, lat_c),
        (-180.0, lat_c),
        (180.0, -90.0),
        (-180.0, -90.0),
    ];
    assert_all_vertices_allowed(&r0, &allowed, "南極で閉じた外環");
}

/// **ISSUE-051 §確定仕様 2・3（北極で閉じる＝極引数が効いている）**: 同一 fixture を
/// `Some(North)` で閉じると、弧は子午線を `lat = +90` まで**上り**、別の（同じく妥当な）多角形になる。
///
/// 期待（開表現・CCW）:
/// `[(−180,−72.75),(−150,−75),(−60,−78),(30,−75),(120,−78),(170,−72),(180,−72.75),(180,90),(−180,90)]`
/// 内部は**北側**（帯より北＝北極を含む側）なので、東進のまま辿る順が CCW（shoelace **+59745**・手計算）。
/// 南極版（+5055・lat −90）と**頂点も面積も異なる**ので、極引数を無視する実装は両方を通せない。
///
/// 殺す実装: 極引数を無視して常に南（または常に北）で閉じる・`EnclosedPole` の分岐を取り違える・
///   `lat = ±90` を `∓90` にする・南北で同じ Value を返す。
#[test]
fn geojson_with_pole_north_differs_from_south_for_same_ring() {
    let lat_c = lat_c_east(170.0, -72.0, -150.0, -75.0);
    let poly = GeoPolygon::new(vec![pole_ring_one_crossing()]);
    let g_north = poly.geojson_geometry_with_pole(Some(EnclosedPole::North));
    let g_south = poly.geojson_geometry_with_pole(Some(EnclosedPole::South));

    // 極引数が出力を変える（無視する実装は必ず落ちる）。
    assert_ne!(
        g_north, g_south,
        "North と South で出力が同一＝極引数が無視されている"
    );

    assert_eq!(
        g_north["type"],
        Value::String("Polygon".to_string()),
        "断片 1 枚 → Polygon"
    );
    let rings = g_north["coordinates"].as_array().expect("リング配列");
    assert_eq!(rings.len(), 1, "外環 1 本のみ, got {}", rings.len());
    let r0 = ring_coords(&rings[0]);

    assert_eq!(r0.first(), r0.last(), "閉リング");
    assert_eq!(
        r0.len(),
        10,
        "入力 5 + 交点 2 + 極 2 = 9 点 + 先頭複製, got {}: {r0:?}",
        r0.len()
    );
    let expected = [
        (-180.0, lat_c),
        (-150.0, -75.0),
        (-60.0, -78.0),
        (30.0, -75.0),
        (120.0, -78.0),
        (170.0, -72.0),
        (180.0, lat_c),
        (180.0, 90.0),
        (-180.0, 90.0),
    ];
    assert_ring_coords_eq(&r0, &expected, "北極で閉じた外環");
    assert!(
        close(signed_area(&r0), 59745.0),
        "外環は CCW・面積 +59745（内部＝北側）, got {}",
        signed_area(&r0)
    );
    assert!(
        r0.iter()
            .any(|&(lon, lat)| close(lon, 180.0) && close(lat, 90.0)),
        "(+180, +90) の極頂点が無い, got {r0:?}"
    );
    assert!(
        r0.iter()
            .any(|&(lon, lat)| close(lon, -180.0) && close(lat, 90.0)),
        "(−180, +90) の極頂点が無い, got {r0:?}"
    );
    assert!(
        !r0.iter().any(|&(_, lat)| close(lat, -90.0)),
        "北極で閉じたのに lat=−90 の頂点がある, got {r0:?}"
    );
    let span = max_non_polar_lon_span(&r0);
    assert!(
        span <= 180.0 + EPS,
        "極以外に経度幅 >180 の辺がある: max={span}, ring={r0:?}"
    );
    let allowed = [
        (-150.0, -75.0),
        (-60.0, -78.0),
        (30.0, -75.0),
        (120.0, -78.0),
        (170.0, -72.0),
        (180.0, lat_c),
        (-180.0, lat_c),
        (180.0, 90.0),
        (-180.0, 90.0),
    ];
    assert_all_vertices_allowed(&r0, &allowed, "北極で閉じた外環");
}

/// **ISSUE-051 §確定仕様 1・3（跨ぎ 3 回でも極で閉じる）**: 「奇数」一般で発火し、
/// **極で閉じるのは対を成さなかった 1 本の弧だけ**（対を成す弧は §11.8(b) どおり子午線で閉じる）。
///
/// 入力（(lon,lat)・8 頂点）:
/// `(−150,−70),(−60,−76),(30,−72),(120,−78),(176,−60),(−176,−64),(172,−64),(150,−70)`
/// Δlon は順に `+90,+90,+90,+56,−352(東進①),+348(西進②),−22,−300(東進③・閉じ辺)`
/// ＝**跨ぎ 3 回（奇数）**、±360 補正後の総変位 `+360`（巻き数 +1 ＝極を囲む）。
/// 交点緯度: ① t=(180−176)/(360−352)=4/8=**0.5** → `−60+0.5·(−4) = −62`（真の内分点）、
/// ②（水平辺）**−64**、③（水平辺）**−70**。
///
/// 弧は 3 本: α `(−180,−70)…(180,−62)`（両端が**逆**の子午線）、
/// β `(−180,−62)→(−176,−64)→(−180,−64)`（両端 −180）、
/// γ `(180,−64)→(172,−64)→(150,−70)→(180,−70)`（両端 +180）。
///
/// **結線規則（第一原理・本テストが縛る核心）**: 子午線上の結線は、その子午線が**領域の内部**を
/// 通る区間でなければならない。「弧の両端が同じ子午線にあるなら自己閉環する」は**誤り**で、
/// 内部区間はパリティで決まる。極を囲むので `lat = −90` は内部、境界は子午線 `±180` を
/// **−70・−64・−62** の 3 点で横切る ⇒ 極から北へ辿ると内部/外部が交互に入れ替わり、
/// 内部区間は **[−90,−70] と [−64,−62]**（[−70,−64] は**外部**）。したがって
/// - `+180`: 区間 [−64,−62] が α の終点(−62) と γ の始点(−64) を結ぶ。区間 [−90,−70] が
///   γ の終点(−70) を極へ落とす。
/// - `−180`: 区間 [−64,−62] が β を自己閉環させる。区間 [−90,−70] が α の始点(−70) を極へ落とす。
///
/// ⇒ **α+γ が 1 本の鎖に繋がり、その両端（±180 の −70）が極で閉じる。断片は 2 枚**。
/// 統一形で言えば「各子午線上の端点を緯度順に並べ、極を囲む場合は**極側の端に極頂点を挿入**してから
/// 先頭から 2 つずつ対にする」＝偶数跨ぎの既存規則（§11.8(b)・`..._four_crossings_pairs_by_latitude`
/// の 3↔12・25↔38）のそのままの拡張である。
///
/// γ を「両端が +180 だから」単独で閉じる実装は、**外部区間 [−70,−64] を内部として取り込み**
/// （偽の面積 114）、**内部区間 [−64,−62] の取り込みを落とす**ため、下の面積保存オラクル
/// （総面積 **6102**）に対して 6330 を返して落ちる。
///
/// 期待（開表現・CCW）:
/// - 大断片（面積 **+6098**）:
///   `[(−180,−90),(180,−90),(180,−70),(150,−70),(172,−64),(180,−64),(180,−62),(176,−60),(120,−78),(30,−72),(−60,−76),(−150,−70),(−180,−70)]`
/// - 小断片（面積 **+4**・幅 4 × 高 2 の直角三角形の 2 倍…幅 4・高 2 の三角形は面積 4）:
///   `[(−180,−64),(−176,−64),(−180,−62)]`
///
/// 殺す実装: `crossings == 1` だけを特別扱いして 3 回では素通しする（逆走辺が残る）・
///   奇数回で**全部の弧**を極で閉じる（β・γ に偽の極頂点が生え面積が爆発する）・
///   極で閉じる弧の選択を誤る（β を極で閉じて α+γ を落とす＝面積の消失）・
///   余った弧を捨てる（大断片が消える）・panic する。
#[test]
fn geojson_with_pole_three_crossings_closes_only_the_unpaired_arc() {
    // (lat, lon) で与える。
    let outer = ring(&[
        (-70.0, -150.0),
        (-76.0, -60.0),
        (-72.0, 30.0),
        (-78.0, 120.0),
        (-60.0, 176.0),
        (-64.0, -176.0),
        (-64.0, 172.0),
        (-70.0, 150.0),
    ]);

    // 前提（独立オラクル）: 跨ぎ 3 回（奇数）・巻き数 +1。
    let lons = [-150.0_f64, -60.0, 30.0, 120.0, 176.0, -176.0, 172.0, 150.0];
    let (crossings, winding) = crossings_and_winding(&lons);
    assert_eq!(crossings, 3, "前提: 跨ぎは 3 回（奇数）, got {crossings}");
    assert!(close(winding, 360.0), "前提: 巻き数 +1, got {winding}");
    // 前提: ① の交点緯度は −62（端点 −60・−64 のどちらとも異なる真の内分点）。
    let lat_c1 = lat_c_east(176.0, -60.0, -176.0, -64.0);
    assert!(close(lat_c1, -62.0), "前提: ① lat_c = −62, got {lat_c1}");

    // 前提（**面積保存オラクル・座標決定と独立**）: 経度を ±360 で連続化した境界と `lat = −90` の
    // 間の符号付き面積（台形則・厳密）。極冠の真の面積であり、正しい分解では
    // **全断片の符号付き面積の和**がこれに一致しなければならない（内部の取りこぼしも外部の
    // 取り込みも、この 1 つの数値で露見する）。
    // 連続化した経度列: −150 → −60 → 30 → 120 → 176 → 184(=−176) → 172 → 150 → 210(=−150+360)。
    let unwrapped = [
        (-150.0_f64, -70.0_f64),
        (-60.0, -76.0),
        (30.0, -72.0),
        (120.0, -78.0),
        (176.0, -60.0),
        (184.0, -64.0),
        (172.0, -64.0),
        (150.0, -70.0),
        (210.0, -70.0),
    ];
    let cap_area: f64 = unwrapped
        .windows(2)
        .map(|w| ((w[0].1 + w[1].1) / 2.0 + 90.0) * (w[1].0 - w[0].0))
        .sum();
    assert!(
        close(cap_area, 6102.0),
        "前提: 極冠の真の面積は 6102（台形則・手計算 1530+1440+1350+1176+224−312−506+1200）, got {cap_area}"
    );

    let g = GeoPolygon::new(vec![outer]).geojson_geometry_with_pole(Some(EnclosedPole::South));
    assert_eq!(
        g["type"],
        Value::String("MultiPolygon".to_string()),
        "断片 2 枚 → MultiPolygon"
    );
    let ps = sorted_polys(&g);
    assert_eq!(
        ps.len(),
        2,
        "断片は 2 枚（極で閉じた大断片 + 西の小三角）, got {}",
        ps.len()
    );
    // ソートキー =(min lon, min lat): 大(−180,−90) < 小(−180,−64)。
    let big = &ps[0];
    let small = &ps[1];
    assert_eq!(big.len(), 1, "大断片は外環のみ, got {}", big.len());
    assert_eq!(small.len(), 1, "小断片は外環のみ, got {}", small.len());

    let expected_big = [
        (-180.0, -90.0),
        (180.0, -90.0),
        (180.0, -70.0),
        (150.0, -70.0),
        (172.0, -64.0),
        (180.0, -64.0),
        (180.0, lat_c1),
        (176.0, -60.0),
        (120.0, -78.0),
        (30.0, -72.0),
        (-60.0, -76.0),
        (-150.0, -70.0),
        (-180.0, -70.0),
    ];
    assert_ring_coords_eq(&big[0], &expected_big, "極で閉じた大断片の外環");
    assert!(
        close(signed_area(&big[0]), 6098.0),
        "大断片は CCW・面積 +6098, got {}",
        signed_area(&big[0])
    );
    assert_ring_coords_eq(
        &small[0],
        &[(-180.0, -64.0), (-176.0, -64.0), (-180.0, lat_c1)],
        "西の小三角（子午線 −180 で自己閉環・極頂点を持たない）",
    );
    assert!(
        close(signed_area(&small[0]), 4.0),
        "小断片は CCW・面積 +4, got {}",
        signed_area(&small[0])
    );

    // 極頂点は大断片にのみ・各子午線に 1 つずつ（閉じ重複を除いて数える）。
    let big_open = &big[0][..big[0].len() - 1];
    let polar_big = big_open
        .iter()
        .filter(|&&(_, lat)| close(lat, -90.0))
        .count();
    assert_eq!(
        polar_big, 2,
        "大断片の極頂点は (±180,−90) の 2 つ, got {polar_big}: {:?}",
        big[0]
    );
    assert!(
        !small[0].iter().any(|&(_, lat)| close(lat.abs(), 90.0)),
        "対を成した弧の断片に極頂点を生やさない, got {:?}",
        small[0]
    );

    // 経度を逆走する辺が無い（両断片とも）。極上の渡り辺（(180,−90)→(−180,−90)・経度幅 360）は
    // 球面では 1 点なので `max_non_polar_lon_span` が除外する。断片が何枚になっても、
    // この除外は「両端が |lat|=90」という判定だけに依るので安全（対を成した断片は極頂点を持たない）。
    for (i, p) in ps.iter().enumerate() {
        let span = max_non_polar_lon_span(&p[0]);
        assert!(
            span <= 180.0 + EPS,
            "断片{i} に経度幅 >180 の辺がある: max={span}, ring={:?}",
            p[0]
        );
    }

    // **面積保存**: 全断片（全リング）の符号付き面積の和 = 極冠の真の面積 6102。
    // γ を単独で閉じる誤った分解は 6212 + 114 + 4 = 6330 になり、ここで落ちる。
    let total: f64 = ps
        .iter()
        .flat_map(|p| p.iter())
        .map(|r| signed_area(r))
        .sum();
    assert!(
        close(total, cap_area),
        "断片の総面積が極冠の真の面積と一致しない（外部の取り込み／内部の取りこぼし）: got {total}, expected {cap_area}"
    );

    // 捏造なし: 全頂点は「入力頂点」「(±180, lat_c)」「(±180, −90)」のいずれか。
    let allowed = [
        (-150.0, -70.0),
        (-60.0, -76.0),
        (30.0, -72.0),
        (120.0, -78.0),
        (176.0, -60.0),
        (-176.0, -64.0),
        (172.0, -64.0),
        (150.0, -70.0),
        (180.0, lat_c1),
        (-180.0, lat_c1),
        (180.0, -64.0),
        (-180.0, -64.0),
        (180.0, -70.0),
        (-180.0, -70.0),
        (180.0, -90.0),
        (-180.0, -90.0),
    ];
    for (i, p) in ps.iter().enumerate() {
        for r in p {
            assert_all_vertices_allowed(r, &allowed, &format!("断片{i}"));
        }
    }
}

/// **ISSUE-051 §確定仕様 2（`None` は現状維持）**: 極を決められない呼び出し側のために、
/// `None` は (3g) §11.8(d) の素通し（分割せず原座標の単一 Polygon）を**そのまま**返す。
/// 併せて `geojson_geometry()` ≡ `geojson_geometry_with_pole(None)` を Value の厳密比較で縛る。
///
/// 期待は既存テスト `geo_polygon_antimeridian_pole_enclosing_odd_crossing_is_not_split` と同一
/// （入力 5 点 + 先頭複製 = 6 点・原座標・shoelace +885・±180 頂点なし）。
///
/// 殺す実装: `None` でも勝手にどちらかの極で閉じる（極頂点が生える・面積が変わる）・
///   `None` で panic / 空を返す・`geojson_geometry()` を別経路にして両者を食い違わせる。
#[test]
fn geojson_with_pole_none_reproduces_todays_unsplit_output() {
    let poly = GeoPolygon::new(vec![pole_ring_one_crossing()]);
    let g_none = poly.geojson_geometry_with_pole(None);
    let g_legacy = poly.geojson_geometry();
    assert_eq!(
        g_none, g_legacy,
        "geojson_geometry() は geojson_geometry_with_pole(None) と同一でなければならない"
    );

    assert_eq!(
        g_none["type"],
        Value::String("Polygon".to_string()),
        "None は素通し → Polygon"
    );
    let rings = g_none["coordinates"].as_array().expect("リング配列");
    assert_eq!(rings.len(), 1, "リング 1 本, got {}", rings.len());
    let r0 = ring_coords(&rings[0]);
    assert_eq!(
        r0.len(),
        6,
        "入力 5 点 + 先頭複製 = 6 点（極頂点を生やさない）, got {}: {r0:?}",
        r0.len()
    );
    assert_ring_coords_eq(
        &r0,
        &[
            (-150.0, -75.0),
            (-60.0, -78.0),
            (30.0, -75.0),
            (120.0, -78.0),
            (170.0, -72.0),
        ],
        "None は原座標のまま素通し",
    );
    assert!(
        !r0.iter()
            .any(|&(lon, lat)| close(lon.abs(), 180.0) || close(lat.abs(), 90.0)),
        "None では ±180 / ±90 の頂点を 1 つも作らない, got {r0:?}"
    );
    assert!(
        close(signed_area(&r0), 885.0),
        "面積は入力どおり +885（捏造も消失もなし）, got {}",
        signed_area(&r0)
    );
}

/// **ISSUE-051 §確定仕様 4（偶数跨ぎは極引数に影響されない）**: 跨ぎ **2 回**のリングでは
/// `Some(North)` / `Some(South)` / `None` / 既存 `geojson_geometry()` の 4 出力が**完全に同一**。
/// 極のロジックが**奇数跨ぎでしか発火しない**ことを縛る。
///
/// 入力は既存テスト `geo_polygon_antimeridian_east_crossing_splits_into_two_polygons` と同じ矩形
/// `(165,10),(−172,10),(−172,40),(165,40)`（東進・跨ぎ 2 回・巻き数 0）。
/// 同一性に加えて具体構造（MultiPolygon 2 枚・東西断片の座標と面積 +450 / +240・
/// **±90 の頂点が 1 つも無い**）も縛り、「4 つとも同じように壊れる」変異を殺す。
///
/// 殺す実装: 極引数を偶数跨ぎにも適用して極頂点を生やす／断片を結合する・
///   極判定の偶奇を取り違える（偶数回を「極を囲む」と見る）・引数の有無で分割経路を変える。
#[test]
fn geojson_with_pole_even_crossing_ring_is_unaffected_by_pole_argument() {
    let lons = [165.0_f64, -172.0, -172.0, 165.0];
    let (crossings, winding) = crossings_and_winding(&lons);
    assert_eq!(crossings, 2, "前提: 跨ぎは 2 回（偶数）, got {crossings}");
    assert!(
        close(winding, 0.0),
        "前提: 巻き数 0（極を囲まない）, got {winding}"
    );

    let poly = GeoPolygon::new(vec![ring(&[
        (10.0, 165.0),
        (10.0, -172.0),
        (40.0, -172.0),
        (40.0, 165.0),
    ])]);
    let g_none = poly.geojson_geometry_with_pole(None);
    let g_north = poly.geojson_geometry_with_pole(Some(EnclosedPole::North));
    let g_south = poly.geojson_geometry_with_pole(Some(EnclosedPole::South));
    let g_legacy = poly.geojson_geometry();

    assert_eq!(g_north, g_none, "偶数跨ぎ: North と None で出力が違う");
    assert_eq!(g_south, g_none, "偶数跨ぎ: South と None で出力が違う");
    assert_eq!(g_legacy, g_none, "偶数跨ぎ: 既存 API と None で出力が違う");

    // 具体構造（(3g) の期待そのまま）。
    assert_eq!(
        g_none["type"],
        Value::String("MultiPolygon".to_string()),
        "偶数跨ぎ → MultiPolygon（従来どおり分割）"
    );
    let ps = sorted_polys(&g_none);
    assert_eq!(ps.len(), 2, "断片は 2 枚, got {}", ps.len());
    let west = &ps[0];
    let east = &ps[1];
    assert_ring_coords_eq(
        &east[0],
        &[(180.0, 40.0), (165.0, 40.0), (165.0, 10.0), (180.0, 10.0)],
        "東断片の外環",
    );
    assert_ring_coords_eq(
        &west[0],
        &[
            (-180.0, 10.0),
            (-172.0, 10.0),
            (-172.0, 40.0),
            (-180.0, 40.0),
        ],
        "西断片の外環",
    );
    assert!(
        close(signed_area(&east[0]), 450.0),
        "東断片 CCW・面積 +450, got {}",
        signed_area(&east[0])
    );
    assert!(
        close(signed_area(&west[0]), 240.0),
        "西断片 CCW・面積 +240, got {}",
        signed_area(&west[0])
    );
    for (i, p) in ps.iter().enumerate() {
        assert!(
            !p[0].iter().any(|&(_, lat)| close(lat.abs(), 90.0)),
            "偶数跨ぎの断片{i} に極頂点を生やさない, got {:?}",
            p[0]
        );
    }
}

/// **ISSUE-051 §確定仕様 4（非跨ぎも極引数に影響されない）**: 跨ぎ **0 回**のリングでは
/// 4 出力（None / North / South / 既存 API）が完全に同一で、`type="Polygon"`・原座標・
/// **±180 も ±90 も現れない**。
///
/// 入力は既存テストと同じ `(170,5),(179,8),(174,21)`（+180 に近いが跨がない・CCW）。
///
/// 殺す実装: 極引数があれば無条件に極で閉じる・跨ぎ 0 を「奇数」と誤判定する・
///   極引数の有無で環向き正規化の経路を変える。
#[test]
fn geojson_with_pole_non_crossing_ring_is_unaffected_by_pole_argument() {
    let lons = [170.0_f64, 179.0, 174.0];
    let (crossings, winding) = crossings_and_winding(&lons);
    assert_eq!(crossings, 0, "前提: 跨ぎは 0 回, got {crossings}");
    assert!(close(winding, 0.0), "前提: 巻き数 0, got {winding}");

    let poly = GeoPolygon::new(vec![ring(&[(5.0, 170.0), (8.0, 179.0), (21.0, 174.0)])]);
    let g_none = poly.geojson_geometry_with_pole(None);
    let g_north = poly.geojson_geometry_with_pole(Some(EnclosedPole::North));
    let g_south = poly.geojson_geometry_with_pole(Some(EnclosedPole::South));
    let g_legacy = poly.geojson_geometry();

    assert_eq!(g_north, g_none, "非跨ぎ: North と None で出力が違う");
    assert_eq!(g_south, g_none, "非跨ぎ: South と None で出力が違う");
    assert_eq!(g_legacy, g_none, "非跨ぎ: 既存 API と None で出力が違う");

    assert_eq!(
        g_none["type"],
        Value::String("Polygon".to_string()),
        "非跨ぎは Polygon のまま"
    );
    let rings = g_none["coordinates"].as_array().expect("リング配列");
    assert_eq!(rings.len(), 1, "リング 1 本, got {}", rings.len());
    let r0 = ring_coords(&rings[0]);
    assert_ring_coords_eq(
        &r0,
        &[(170.0, 5.0), (179.0, 8.0), (174.0, 21.0)],
        "非跨ぎリング（原座標・CCW のまま）",
    );
    assert!(
        !r0.iter()
            .any(|&(lon, lat)| close(lon.abs(), 180.0) || close(lat.abs(), 90.0)),
        "±180 / ±90 の頂点を作らない, got {r0:?}"
    );
}

/// **ISSUE-051 §確定仕様 3（退化: 交点緯度が既に極上）**: 交点緯度がちょうど `−90` になる
/// 極を囲むリングでは、極で閉じるために挿入すべき点 `(±180,−90)` が**補間交点そのもの**と一致し、
/// 「子午線を極まで辿る」区間が**長さ 0** になる。
///
/// 入力（(lon,lat)）: `(−150,−90),(−60,−78),(30,−75),(120,−78),(170,−90)`。
/// 閉じ辺 (170,−90)→(−150,−90) は `Δlon = −320` ＝東進跨ぎ **1 回（奇数）**・巻き数 +1。
/// 水平辺なので `lat_c = −90` ちょうど。`EastLongitude` の正規化は緯度に触れないので
/// `lat = ±90` の入力頂点は**構成可能**＝この退化は到達しうる（構成不能ではない）。
///
/// 仕様は重複頂点の除去を定めていないため**頂点数は縛らない**（`(±180,−90)` が 1 つでも
/// 2 つ連続でも可）。代わりに退化しても壊れない性質を縛る:
/// (a) panic せず断片 1 枚、(b) 閉リング、(c) 極以外に経度幅 >180 の辺が無い、
/// (d) 全頂点が「入力頂点 ∪ (±180,−90)」に属する（捏造なし）、(e) 両子午線の極頂点が存在する、
/// (f) 面積は **+3270**（重複頂点は面積に寄与しないので除去の有無に依らず一意・CCW）。
///
/// 殺す実装: 退化を検出して素通しに落とす（逆走辺 320 が残る）・長さ 0 の子午線区間で
///   0 除算 / unwrap して panic する・極頂点を捏造経度で足す・向きを CW のまま出す。
#[test]
fn geojson_with_pole_degenerate_crossing_at_pole_latitude() {
    let outer = ring(&[
        (-90.0, -150.0),
        (-78.0, -60.0),
        (-75.0, 30.0),
        (-78.0, 120.0),
        (-90.0, 170.0),
    ]);
    // 前提: 跨ぎ 1 回（奇数）・巻き数 +1・交点緯度ちょうど −90。
    let lons = [-150.0_f64, -60.0, 30.0, 120.0, 170.0];
    let (crossings, winding) = crossings_and_winding(&lons);
    assert_eq!(crossings, 1, "前提: 跨ぎ 1 回（奇数）, got {crossings}");
    assert!(close(winding, 360.0), "前提: 巻き数 +1, got {winding}");
    let lat_c = lat_c_east(170.0, -90.0, -150.0, -90.0);
    assert!(
        close(lat_c, -90.0),
        "前提: 交点緯度は極ちょうど, got {lat_c}"
    );

    let g = GeoPolygon::new(vec![outer]).geojson_geometry_with_pole(Some(EnclosedPole::South));

    // (a) panic せず Polygon（断片 1 枚）。
    assert_eq!(
        g["type"],
        Value::String("Polygon".to_string()),
        "断片 1 枚 → Polygon"
    );
    let rings = g["coordinates"].as_array().expect("リング配列");
    assert_eq!(rings.len(), 1, "リング 1 本, got {}", rings.len());
    let r0 = ring_coords(&rings[0]);

    // (b) 閉リング。
    assert_eq!(r0.first(), r0.last(), "閉リング（先頭==末尾）");

    // (c) 極以外に逆走辺が無い。
    let span = max_non_polar_lon_span(&r0);
    assert!(
        span <= 180.0 + EPS,
        "極以外に経度幅 >180 の辺がある（退化で素通しに落ちた）: max={span}, ring={r0:?}"
    );

    // (d) 捏造なし。
    let allowed = [
        (-150.0, -90.0),
        (-60.0, -78.0),
        (30.0, -75.0),
        (120.0, -78.0),
        (170.0, -90.0),
        (180.0, -90.0),
        (-180.0, -90.0),
    ];
    assert_all_vertices_allowed(&r0, &allowed, "退化（交点が極上）のリング");

    // (e) 両子午線に極頂点。
    assert!(
        r0.iter()
            .any(|&(lon, lat)| close(lon, 180.0) && close(lat, -90.0)),
        "(+180,−90) が無い, got {r0:?}"
    );
    assert!(
        r0.iter()
            .any(|&(lon, lat)| close(lon, -180.0) && close(lat, -90.0)),
        "(−180,−90) が無い, got {r0:?}"
    );

    // (f) 面積 +3270（CCW・重複頂点の有無に依らず一意）。
    assert!(
        close(signed_area(&r0), 3270.0),
        "CCW・面積 +3270（重複頂点は面積に寄与しない）, got {}",
        signed_area(&r0)
    );
}
