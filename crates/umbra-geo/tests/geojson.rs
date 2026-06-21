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
