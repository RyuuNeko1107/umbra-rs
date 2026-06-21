//! 地理座標の幾何プリミティブ（`docs/api-draft.md` §4）。
//!
//! 中心線・限界線・部分食域の表現に用いる点 [`GeoPoint`]・折れ線 [`GeoLine`]・
//! 多角形 [`GeoPolygon`] を提供する。緯度・経度の検証/正規化は
//! `umbra-core` の [`GeodeticLatitude`](umbra_core::GeodeticLatitude) /
//! [`EastLongitude`](umbra_core::EastLongitude) に委譲し、ここでは**伝播と保持**のみを担う。

use umbra_core::{DomainError, EastLongitude, GeodeticLatitude};

/// 地理座標点（測地緯度・東経）。
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GeoPoint {
    /// 測地緯度。
    pub lat: GeodeticLatitude,
    /// 東経（東を正）。
    pub lon: EastLongitude,
}

impl GeoPoint {
    /// 構成済みの緯度・経度から点を構築する。
    pub fn new(lat: GeodeticLatitude, lon: EastLongitude) -> Self {
        Self { lat, lon }
    }

    /// 度から構築する。緯度範囲外は [`DomainError::OutOfRange`]（[`GeodeticLatitude`] 検証を伝播）、
    /// 経度は正規化（[`EastLongitude`]・infallible）。
    pub fn from_degrees(lat_deg: f64, lon_deg: f64) -> Result<Self, DomainError> {
        Ok(Self::new(
            GeodeticLatitude::from_degrees(lat_deg)?,
            EastLongitude::from_degrees(lon_deg),
        ))
    }

    /// GeoJSON Point ジオメトリ（`{"type":"Point","coordinates":[経度, 緯度]}`・M9.2）。
    ///
    /// 座標順は GeoJSON RFC 7946 の **[経度, 緯度]**（lon, lat の順）。値は度（公開入出力は度・conventions §3）。
    pub fn geojson_geometry(&self) -> serde_json::Value {
        serde_json::json!({
            "type": "Point",
            "coordinates": [self.lon.degrees().0, self.lat.degrees().0],
        })
    }
}

/// `GeoPoint` の JSON 表現（ISSUE-031 S31b・A7: 数値の単位はフィールド名で明示）。
/// 公開入出力は測地緯度・東経の度（conventions §3）なので `{ "lat_deg", "lon_deg" }`。
/// 内部表現（ラジアン）でなく度で出力する。Serialize のみ（api-draft §0）。
#[cfg(feature = "serde")]
impl serde::Serialize for GeoPoint {
    fn serialize<S: serde::Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        use serde::ser::SerializeStruct;
        let mut st = serializer.serialize_struct("GeoPoint", 2)?;
        st.serialize_field("lat_deg", &self.lat.degrees().0)?;
        st.serialize_field("lon_deg", &self.lon.degrees().0)?;
        st.end()
    }
}

/// 折れ線（中心線・限界線）。
#[derive(Clone, Debug, PartialEq)]
pub struct GeoLine {
    /// 構成点（順序付き）。
    pub points: Vec<GeoPoint>,
}

impl GeoLine {
    /// 構成点から折れ線を構築する（空も可）。
    pub fn new(points: Vec<GeoPoint>) -> Self {
        Self { points }
    }

    /// GeoJSON 折れ線ジオメトリ（M9.2 / M9.5 ±180 補間）。日付変更線（±180°）を跨ぐ場合は
    /// `MultiLineString` に分割し、跨ぎが無ければ `LineString`。座標は [経度, 緯度] 順（RFC 7946）。
    ///
    /// 跨ぎ判定: 連続 2 点 (lon1,lat1)→(lon2,lat2) の `Δlon=lon2−lon1` が `Δlon < −180`（東進・+180 越え）
    /// または `Δlon > 180`（西進・−180 越え）。跨ぎ点では交点緯度を子午線上に**線形補間**し、前セグメントの
    /// 末尾と次セグメントの先頭に境界点を補う（RFC 7946 §3.1.9・隙間を残さない）。
    /// 東進: `t=(180−lon1)/(360+Δlon)`, 末尾 `[+180, lat_c]` / 先頭 `[−180, lat_c]`。
    /// 西進: `t=(lon1+180)/(360−Δlon)`, 末尾 `[−180, lat_c]` / 先頭 `[+180, lat_c]`。`lat_c=lat1+t·(lat2−lat1)`。
    /// `|Δlon|=180` ちょうどは跨ぎとしない（測度ゼロ境界）。点 0/1 個の退行は空/単一座標の `LineString`。
    pub fn geojson_geometry(&self) -> serde_json::Value {
        // [経度, 緯度] 列を跨ぎ位置で分割し、交点を ±180 子午線へ線形補間して両端に補う。
        let mut segments: Vec<Vec<[f64; 2]>> = Vec::new();
        let mut current: Vec<[f64; 2]> = Vec::new();
        let mut prev: Option<[f64; 2]> = None;
        for p in &self.points {
            let lon = p.lon.degrees().0;
            let lat = p.lat.degrees().0;
            if let Some([prev_lon, prev_lat]) = prev {
                let delta = lon - prev_lon;
                if delta < -180.0 {
                    // 東進（+180 越え）: prev_lon→+180→−180→lon。
                    let t = (180.0 - prev_lon) / (360.0 + delta);
                    let lat_c = prev_lat + t * (lat - prev_lat);
                    current.push([180.0, lat_c]);
                    segments.push(std::mem::take(&mut current));
                    current.push([-180.0, lat_c]);
                } else if delta > 180.0 {
                    // 西進（−180 越え）: prev_lon→−180→+180→lon。
                    let t = (prev_lon + 180.0) / (360.0 - delta);
                    let lat_c = prev_lat + t * (lat - prev_lat);
                    current.push([-180.0, lat_c]);
                    segments.push(std::mem::take(&mut current));
                    current.push([180.0, lat_c]);
                }
            }
            current.push([lon, lat]);
            prev = Some([lon, lat]);
        }
        segments.push(current);

        if segments.len() > 1 {
            serde_json::json!({ "type": "MultiLineString", "coordinates": segments })
        } else {
            // 跨ぎ無し（1 セグメント）。0/1 点なら空/単一座標の LineString。
            let coordinates = segments.into_iter().next().unwrap_or_default();
            serde_json::json!({ "type": "LineString", "coordinates": coordinates })
        }
    }
}

/// 多角形（部分食域。外周＋穴のリング列）。
#[derive(Clone, Debug, PartialEq)]
pub struct GeoPolygon {
    /// リング列（[0]=外周、以降=穴。各リングは点列）。
    pub rings: Vec<Vec<GeoPoint>>,
}

impl GeoPolygon {
    /// リング列から多角形を構築する（空も可）。
    pub fn new(rings: Vec<Vec<GeoPoint>>) -> Self {
        Self { rings }
    }
}

#[cfg(test)]
mod tests {
    use crate::geometry::{GeoLine, GeoPoint, GeoPolygon};
    use umbra_core::{DomainError, EastLongitude, GeodeticLatitude};

    /// 度の exact 比較許容（型委譲した値の桁落ち程度）。
    const TOL_DEG: f64 = 1e-9;

    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    // ============================================================
    // GeoPoint::new — フィールド保持
    // ============================================================

    /// `GeoPoint::new` は受け取った lat/lon をそのまま保持し、取り違えない。
    /// 緯度と経度を**非対称な値**（緯度 12.0°, 経度 34.0°）で与え、フィールドが
    /// 入れ替わる変異（lat↔lon スワップ）・別フィールド代入を撃破する。
    #[test]
    fn geo_point_new_preserves_lat_and_lon_without_swap() {
        let lat = GeodeticLatitude::from_degrees(12.0).expect("12° は有効");
        let lon = EastLongitude::from_degrees(34.0);
        let p = GeoPoint::new(lat, lon);
        // lat フィールドは渡した緯度そのもの（経度ではない）。
        assert_eq!(p.lat, lat);
        // lon フィールドは渡した経度そのもの（緯度ではない）。
        assert_eq!(p.lon, lon);
        // 取り違え検出の補助：度値でも区別される。
        assert!(
            close(p.lat.degrees().0, 12.0, TOL_DEG),
            "lat = {}",
            p.lat.degrees().0
        );
        assert!(
            close(p.lon.degrees().0, 34.0, TOL_DEG),
            "lon = {}",
            p.lon.degrees().0
        );
    }

    // ============================================================
    // GeoPoint::from_degrees — 委譲・伝播・正規化
    // ============================================================

    /// 有効な内部点（岡山 ≈ 北緯34.66°, 東経133.92°）→ Ok・値保持。
    /// `lat.degrees()≈34.66`, `lon.degrees()≈133.92`（範囲内なので不変）。
    /// 度→型構築の取り違え・lat/lon 引数順の取り違えを撃破する。
    #[test]
    fn geo_point_from_degrees_valid_preserves_values() {
        let p = GeoPoint::from_degrees(34.66, 133.92).expect("岡山は有効");
        assert!(
            close(p.lat.degrees().0, 34.66, TOL_DEG),
            "lat = {}",
            p.lat.degrees().0
        );
        assert!(
            close(p.lon.degrees().0, 133.92, TOL_DEG),
            "lon = {}",
            p.lon.degrees().0
        );
    }

    /// 緯度が範囲外（91°）→ `Err(OutOfRange{what:"geodetic latitude"})` を伝播する。
    /// 緯度検証の脱落（`?` 伝播の欠落・unwrap 等）・what 文字列改変を撃破する。
    #[test]
    fn geo_point_from_degrees_propagates_latitude_out_of_range() {
        let err = GeoPoint::from_degrees(91.0, 0.0).expect_err("緯度91°は範囲外");
        assert_eq!(
            err,
            DomainError::OutOfRange {
                what: "geodetic latitude"
            }
        );
    }

    /// 経度は正規化される（190° → -170°）。経度の正規化脱落（生値保持）を撃破する。
    /// なお緯度 0° は有効なので全体は Ok。
    #[test]
    fn geo_point_from_degrees_normalizes_longitude() {
        let p = GeoPoint::from_degrees(0.0, 190.0).expect("経度は正規化され有効");
        assert!(
            close(p.lon.degrees().0, -170.0, TOL_DEG),
            "lon = {}",
            p.lon.degrees().0
        );
    }

    /// 引数順を縛る非対称ケース：lat=10°(有効)・lon=200°(→-160° に正規化)。
    /// 仮に lat と lon を取り違えると lat=200° が緯度検証で Err になる/値がずれるため、
    /// 引数順の取り違え（lat_deg と lon_deg のスワップ）を撃破する。
    #[test]
    fn geo_point_from_degrees_argument_order_is_lat_then_lon() {
        let p = GeoPoint::from_degrees(10.0, 200.0).expect("lat=10°,lon=200° は有効");
        assert!(
            close(p.lat.degrees().0, 10.0, TOL_DEG),
            "lat = {}",
            p.lat.degrees().0
        );
        assert!(
            close(p.lon.degrees().0, -160.0, TOL_DEG),
            "lon = {}",
            p.lon.degrees().0
        );
    }

    // ============================================================
    // GeoLine::new — points 保持
    // ============================================================

    /// `GeoLine::new` は複数点をそのまま（順序込みで）保持する。
    /// 点の落とし込み・順序入れ替え・空化する変異を撃破する。
    #[test]
    fn geo_line_new_preserves_points_in_order() {
        let p0 = GeoPoint::from_degrees(0.0, 0.0).unwrap();
        let p1 = GeoPoint::from_degrees(10.0, 20.0).unwrap();
        let p2 = GeoPoint::from_degrees(-5.0, 30.0).unwrap();
        let line = GeoLine::new(vec![p0, p1, p2]);
        assert_eq!(line.points, vec![p0, p1, p2]);
        // 順序が保たれること（先頭/末尾を固定）。
        assert_eq!(line.points.first(), Some(&p0));
        assert_eq!(line.points.last(), Some(&p2));
    }

    /// `GeoLine::new` は空 Vec も許容しそのまま保持する（点を捏造しない）。
    #[test]
    fn geo_line_new_allows_empty() {
        let line = GeoLine::new(Vec::new());
        assert!(line.points.is_empty());
    }

    // ============================================================
    // GeoPolygon::new — rings 保持
    // ============================================================

    /// `GeoPolygon::new` は外周＋穴の2リングをそのまま（順序込みで）保持する。
    /// リングの落とし込み・順序入れ替え（外周と穴の取り違え）・空化を撃破する。
    #[test]
    fn geo_polygon_new_preserves_outer_and_hole_rings() {
        let outer = vec![
            GeoPoint::from_degrees(0.0, 0.0).unwrap(),
            GeoPoint::from_degrees(0.0, 10.0).unwrap(),
            GeoPoint::from_degrees(10.0, 10.0).unwrap(),
        ];
        let hole = vec![
            GeoPoint::from_degrees(2.0, 2.0).unwrap(),
            GeoPoint::from_degrees(2.0, 4.0).unwrap(),
        ];
        let poly = GeoPolygon::new(vec![outer.clone(), hole.clone()]);
        assert_eq!(poly.rings.len(), 2);
        // 第0リング＝外周、第1リング＝穴（順序が保たれる）。
        assert_eq!(poly.rings[0], outer);
        assert_eq!(poly.rings[1], hole);
    }

    /// `GeoPolygon::new` は空 rings も許容しそのまま保持する（リングを捏造しない）。
    #[test]
    fn geo_polygon_new_allows_empty() {
        let poly = GeoPolygon::new(Vec::new());
        assert!(poly.rings.is_empty());
    }

    // ============================================================
    // derive 特性（Copy / Clone / PartialEq）
    // ============================================================

    /// `GeoPoint` は `Copy`（move されず複製で渡る）。`#[derive(Copy)]` の脱落を
    /// コンパイル時に撃破する（消費後も元の束縛が有効）。
    #[test]
    fn geo_point_is_copy() {
        fn assert_copy<T: Copy>(_: T) {}
        let p = GeoPoint::from_degrees(34.66, 133.92).unwrap();
        let p2 = p; // Copy なら move されない
        assert_copy(p2);
        assert_eq!(p, p2); // 元の束縛 p も有効
    }

    /// `GeoPoint` の `PartialEq` は値の同一性を区別する（異値で ne）。
    /// 常時 true 化する変異・lat だけ/lon だけ比較する変異を撃破する。
    #[test]
    fn geo_point_partial_eq_distinguishes_values() {
        let base = GeoPoint::from_degrees(10.0, 20.0).unwrap();
        // 緯度だけ異なる。
        assert_ne!(base, GeoPoint::from_degrees(11.0, 20.0).unwrap());
        // 経度だけ異なる。
        assert_ne!(base, GeoPoint::from_degrees(10.0, 21.0).unwrap());
        // 同値は等しい。
        assert_eq!(base, GeoPoint::from_degrees(10.0, 20.0).unwrap());
    }

    /// `GeoLine` は `Clone`・`PartialEq`。複製は等しく、異なる points 列は不一致。
    /// Clone がフィールドを取りこぼす/PartialEq を常時 true 化する変異を撃破する。
    #[test]
    fn geo_line_clone_and_partial_eq() {
        let p0 = GeoPoint::from_degrees(0.0, 0.0).unwrap();
        let p1 = GeoPoint::from_degrees(10.0, 20.0).unwrap();
        let line = GeoLine::new(vec![p0, p1]);
        let cloned = line.clone();
        assert_eq!(line, cloned);
        // 点が1つ少ない別の線は不一致。
        assert_ne!(line, GeoLine::new(vec![p0]));
    }

    /// `GeoPolygon` は `Clone`・`PartialEq`。複製は等しく、異なる rings は不一致。
    #[test]
    fn geo_polygon_clone_and_partial_eq() {
        let ring = vec![
            GeoPoint::from_degrees(0.0, 0.0).unwrap(),
            GeoPoint::from_degrees(0.0, 10.0).unwrap(),
        ];
        let poly = GeoPolygon::new(vec![ring.clone()]);
        let cloned = poly.clone();
        assert_eq!(poly, cloned);
        // リング数が異なる多角形は不一致。
        assert_ne!(poly, GeoPolygon::new(vec![ring.clone(), ring]));
    }
}
