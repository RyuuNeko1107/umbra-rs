//! 経路型（中心線・限界線・部分食域・各サンプル点）。`docs/api-draft.md` §4。
//!
//! [`EclipsePath`] は [`EclipseEngine::path`](crate) が返す経路の表現で、中心線
//! [`EclipsePath::center_line`]・北/南限界線・部分食域・最大食点・各サンプル点
//! [`PathSample`] を保持する。[`PathOptions`] は経路生成のオプション。
//!
//! 本スライス（ISSUE-043 S5a）では**型定義（pub フィールドのデータコンテナ）のみ**。
//! v0.1 では `engine.path()` 本体は未実装（`Err(NotImplemented)` を返す, S5b）。
//!
//! **crate 配置（api-draft §4 からの意図的逸脱）**: api-draft §4 はこれらを umbra-geo に置くが、
//! [`PathSample::kind`] が [`SolarEclipseKind`]（umbra-eclipse）に依存し、umbra-eclipse は
//! GeoPoint 等で umbra-geo に依存するため、umbra-geo 側に置くと geo↔eclipse 循環になる。よって
//! 経路 **result** 型は engine 側（umbra-eclipse）に置く。umbra-geo は純幾何（GeoPoint/Line/Polygon）
//! と GeoJSON 出力のみ（M9 で経路生成本体を実装する際もこの配置を維持＝循環回避）。

use umbra_core::{Degrees, Kilometers, UtcInstant};
use umbra_geo::{GeoLine, GeoPoint, GeoPolygon};

use crate::global::SolarEclipseKind;

/// 経路（中心線・北/南限界線・部分食域・最大食点・各サンプル点）。
///
/// v0.1 では `EclipseEngine::path()` が `Err(EclipseError::NotImplemented)` を返すため、本型は
/// 経路実装（umbra-geo, M9）まで結果として生成されない（型のみ提供）。
#[derive(Clone, Debug)]
pub struct EclipsePath {
    /// 中心線（中心食のみ Some）。
    pub center_line: Option<GeoLine>,
    /// 北限界線（中心食のみ Some）。
    pub northern_limit: Option<GeoLine>,
    /// 南限界線（中心食のみ Some）。
    pub southern_limit: Option<GeoLine>,
    /// 部分食域（外周＋穴）。
    pub partial_limit: Option<GeoPolygon>,
    /// 最大食地点（常に存在）。
    pub greatest_point: GeoPoint,
    /// 経路サンプル点列。
    pub samples: Vec<PathSample>,
}

impl EclipsePath {
    /// 経路を GeoJSON `FeatureCollection`（pretty・末尾改行）に直列化する（M9.2 / M9.5 限界線）。
    ///
    /// feature を決定的順序で含む: `greatest_point`（Point・`role="greatest"`）→ `center_line`（Some 時・
    /// `role="center_line"`）→ `northern_limit`（Some 時・`role="northern_limit"`）→ `southern_limit`
    /// （Some 時・`role="southern_limit"`）→ `partial_limit`（Some 時・`role="partial_limit"`・M9 残(3) 3d）。折れ線は
    /// [`GeoLine::geojson_geometry`]（LineString/MultiLineString・±180 補間）、部分食域は
    /// [`GeoPolygon::geojson_geometry`]（Polygon・閉リング・環向き正規化・v1 は反子午線非分割）。**`samples` は未出力**。
    /// 座標順は [経度, 緯度]（RFC 7946）。
    pub fn to_geojson(&self) -> Result<String, serde_json::Error> {
        let mut features = vec![serde_json::json!({
            "type": "Feature",
            "geometry": self.greatest_point.geojson_geometry(),
            "properties": { "role": "greatest" },
        })];
        // 折れ線 feature を決定的順序（center_line → northern_limit → southern_limit）で追加。
        for (line, role) in [
            (&self.center_line, "center_line"),
            (&self.northern_limit, "northern_limit"),
            (&self.southern_limit, "southern_limit"),
        ] {
            if let Some(line) = line {
                features.push(serde_json::json!({
                    "type": "Feature",
                    "geometry": line.geojson_geometry(),
                    "properties": { "role": role },
                }));
            }
        }
        // 部分食域 feature（Some 時・southern_limit の後）。
        if let Some(polygon) = &self.partial_limit {
            features.push(serde_json::json!({
                "type": "Feature",
                "geometry": polygon.geojson_geometry(),
                "properties": { "role": "partial_limit" },
            }));
        }
        let collection = serde_json::json!({
            "type": "FeatureCollection",
            "features": features,
        });
        let mut out = serde_json::to_string_pretty(&collection)?;
        out.push('\n');
        Ok(out)
    }
}

/// 経路上の 1 サンプル点（時刻・中心点・継続・太陽高度・帯幅・種別）。
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PathSample {
    /// サンプル時刻（UTC）。
    pub time_utc: UtcInstant,
    /// 中心線上の地表点。
    pub center: GeoPoint,
    /// 中心食継続時間 \[s\]。
    pub duration_seconds: f64,
    /// 太陽高度。
    pub sun_altitude: Degrees,
    /// 帯幅 \[km\]。
    pub path_width: Kilometers,
    /// 日食種別。
    pub kind: SolarEclipseKind,
}

/// 経路生成オプション。
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PathOptions {
    /// サンプル間隔 \[s\]。
    pub sample_interval_seconds: f64,
    /// 北/南限界線・部分食域を含めるか。
    pub include_limits: bool,
    /// 日付変更線で経路を分割するか（GeoJSON 出力向け）。
    pub split_antimeridian: bool,
}

impl Default for PathOptions {
    /// 既定: サンプル間隔 60 s・限界線含む・日付変更線分割あり。
    fn default() -> Self {
        Self {
            sample_interval_seconds: 60.0,
            include_limits: true,
            split_antimeridian: true,
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::global::SolarEclipseKind;
    use crate::path::{EclipsePath, PathOptions, PathSample};
    use umbra_core::{Degrees, Kilometers, UtcInstant};
    use umbra_geo::{GeoLine, GeoPoint, GeoPolygon};

    /// テスト用 UtcInstant（適当な既知日時）。
    fn t(year: i32, month: u8, day: u8) -> UtcInstant {
        UtcInstant::from_gregorian(year, month, day, 0, 0, 0.0).expect("有効な暦日")
    }

    /// テスト用 GeoPoint。
    fn pt(lat: f64, lon: f64) -> GeoPoint {
        GeoPoint::from_degrees(lat, lon).expect("有効な緯度経度")
    }

    /// 既定値とは別の、非対称な値を持つ PathSample を作る。
    /// duration_seconds≠path_width、sun_altitude も区別できる値にしてフィールド
    /// 取り違え変異を撃破できるようにする。
    fn sample() -> PathSample {
        PathSample {
            time_utc: t(2035, 9, 2),
            center: pt(34.0, 135.0),
            duration_seconds: 123.0,
            sun_altitude: Degrees(45.0),
            path_width: Kilometers(67.0),
            kind: SolarEclipseKind::Total,
        }
    }

    // ============================================================
    // EclipsePath — 6 フィールド保持
    // ============================================================

    /// `EclipsePath` は 6 フィールドを struct リテラルで構築し、そのまま保持する。
    /// center_line=Some・northern/southern_limit=Some・partial_limit=Some・
    /// greatest_point（非 Option）・samples 複数を与え、各フィールドが取り違えなく
    /// 保持されることを確認する（フィールド代入の取り違え変異を撃破）。
    #[test]
    fn eclipse_path_preserves_all_six_fields() {
        let center = GeoLine::new(vec![pt(0.0, 0.0), pt(1.0, 1.0)]);
        let north = GeoLine::new(vec![pt(2.0, 2.0)]);
        let south = GeoLine::new(vec![pt(3.0, 3.0)]);
        let partial = GeoPolygon::new(vec![vec![pt(4.0, 4.0), pt(5.0, 5.0)]]);
        let greatest = pt(6.0, 6.0);
        let s0 = sample();
        let mut s1 = sample();
        s1.duration_seconds = 200.0;

        let path = EclipsePath {
            center_line: Some(center.clone()),
            northern_limit: Some(north.clone()),
            southern_limit: Some(south.clone()),
            partial_limit: Some(partial.clone()),
            greatest_point: greatest,
            samples: vec![s0, s1],
        };

        // center_line と northern/southern_limit は GeoLine（取り違えないこと）。
        assert_eq!(path.center_line, Some(center));
        assert_eq!(path.northern_limit, Some(north));
        assert_eq!(path.southern_limit, Some(south));
        // partial_limit は GeoPolygon。
        assert_eq!(path.partial_limit, Some(partial));
        // greatest_point は GeoPoint（非 Option）。
        assert_eq!(path.greatest_point, greatest);
        // samples は与えた 2 点を順序込みで保持。
        assert_eq!(path.samples, vec![s0, s1]);
    }

    /// `center_line`/各限界線/`partial_limit` は `Option`、`None` も保持できる。
    /// `samples` は空 Vec も保持できる。Option フィールドが `Some` に固定される変異・
    /// samples を捏造する変異を撃破する。
    #[test]
    fn eclipse_path_allows_none_options_and_empty_samples() {
        let path = EclipsePath {
            center_line: None,
            northern_limit: None,
            southern_limit: None,
            partial_limit: None,
            greatest_point: pt(10.0, 20.0),
            samples: Vec::new(),
        };
        // Option フィールドはすべて None を保持する。
        assert_eq!(path.center_line, None);
        assert_eq!(path.northern_limit, None);
        assert_eq!(path.southern_limit, None);
        assert_eq!(path.partial_limit, None);
        // greatest_point は非 Option なので常に値を持つ。
        assert_eq!(path.greatest_point, pt(10.0, 20.0));
        // samples は空のまま（点を捏造しない）。
        assert!(path.samples.is_empty());
    }

    /// `EclipsePath` は `Clone`（GeoLine/GeoPolygon/Vec を含むため Copy 不可）。
    /// Clone がフィールドを取りこぼさず一致することを確認する。
    /// PartialEq 派生は無いため、フィールドごとに同値を確認する。
    #[test]
    fn eclipse_path_clone_preserves_fields() {
        let center = GeoLine::new(vec![pt(0.0, 0.0), pt(1.0, 1.0)]);
        let partial = GeoPolygon::new(vec![vec![pt(4.0, 4.0)]]);
        let s = sample();
        let path = EclipsePath {
            center_line: Some(center.clone()),
            northern_limit: None,
            southern_limit: Some(GeoLine::new(vec![pt(3.0, 3.0)])),
            partial_limit: Some(partial.clone()),
            greatest_point: pt(6.0, 6.0),
            samples: vec![s],
        };

        let cloned = path.clone();
        // Clone は元と同じ内容（PartialEq 非派生のためフィールド単位で確認）。
        assert_eq!(cloned.center_line, path.center_line);
        assert_eq!(cloned.northern_limit, path.northern_limit);
        assert_eq!(cloned.southern_limit, path.southern_limit);
        assert_eq!(cloned.partial_limit, path.partial_limit);
        assert_eq!(cloned.greatest_point, path.greatest_point);
        assert_eq!(cloned.samples, path.samples);
        // 元の path もまだ使える（Clone であり move していない）。
        assert_eq!(path.samples.len(), 1);
    }

    // ============================================================
    // PathSample — 6 フィールド保持・derive
    // ============================================================

    /// `PathSample` は 6 フィールドを構築・保持し、取り違えない。
    /// duration_seconds=123.0 と path_width=67.0 は非対称値にしてあり、両者を
    /// 取り違える変異・sun_altitude を別フィールドから読む変異を撃破する。
    #[test]
    fn path_sample_preserves_all_six_fields() {
        let s = sample();
        assert_eq!(s.time_utc, t(2035, 9, 2));
        assert_eq!(s.center, pt(34.0, 135.0));
        // duration_seconds（秒）と path_width（km）は別物（取り違え検出）。
        assert_eq!(s.duration_seconds, 123.0);
        assert_eq!(s.path_width, Kilometers(67.0));
        // sun_altitude は Degrees。
        assert_eq!(s.sun_altitude, Degrees(45.0));
        assert_eq!(s.kind, SolarEclipseKind::Total);
    }

    /// `PathSample` は `Copy`（move されず複製で渡る）。`#[derive(Copy)]` の脱落を
    /// コンパイル時に撃破する（消費後も元の束縛が有効）。
    #[test]
    fn path_sample_is_copy() {
        fn assert_copy<T: Copy>(_: T) {}
        let s = sample();
        let s2 = s; // Copy なら move されない
        assert_copy(s2);
        assert_eq!(s, s2); // 元の束縛 s も有効
    }

    /// `PathSample` の `PartialEq` は値の同一性を区別する（同値で eq、異値で ne）。
    /// 各フィールドを 1 つずつ変えて不一致になることを確認し、PartialEq を常時 true 化
    /// する変異・一部フィールドのみ比較する変異を撃破する。
    #[test]
    fn path_sample_partial_eq_distinguishes_each_field() {
        let base = sample();
        // 同値は等しい。
        assert_eq!(base, sample());

        // time_utc だけ異なる。
        let mut a = sample();
        a.time_utc = t(2035, 9, 3);
        assert_ne!(base, a);

        // center だけ異なる。
        let mut b = sample();
        b.center = pt(35.0, 135.0);
        assert_ne!(base, b);

        // duration_seconds だけ異なる。
        let mut c = sample();
        c.duration_seconds = 999.0;
        assert_ne!(base, c);

        // sun_altitude だけ異なる。
        let mut d = sample();
        d.sun_altitude = Degrees(46.0);
        assert_ne!(base, d);

        // path_width だけ異なる。
        let mut e = sample();
        e.path_width = Kilometers(68.0);
        assert_ne!(base, e);

        // kind だけ異なる。
        let mut f = sample();
        f.kind = SolarEclipseKind::Annular;
        assert_ne!(base, f);
    }

    // ============================================================
    // PathOptions — 3 フィールド保持・Default・derive
    // ============================================================

    /// `PathOptions` は 3 フィールドを構築・保持する。include_limits と
    /// split_antimeridian を**異なる bool**（false / true）にしてあり、両 bool
    /// フィールドの取り違え変異を撃破する。
    #[test]
    fn path_options_preserves_three_fields() {
        let opts = PathOptions {
            sample_interval_seconds: 30.0,
            include_limits: false,
            split_antimeridian: true,
        };
        assert_eq!(opts.sample_interval_seconds, 30.0);
        // include_limits=false（split_antimeridian と取り違えていない）。
        assert!(!opts.include_limits);
        // split_antimeridian=true（include_limits と取り違えていない）。
        assert!(opts.split_antimeridian);
    }

    /// `PathOptions::default()` の exact 値: sample_interval_seconds=60.0,
    /// include_limits=true, split_antimeridian=true。
    /// Default 値の改変（60.0→他値・true→false・両 bool の取り違え）を撃破する。
    #[test]
    fn path_options_default_exact_values() {
        let d = PathOptions::default();
        // 既定のサンプル間隔は 60.0 秒（exact）。
        assert_eq!(d.sample_interval_seconds, 60.0);
        // 既定で限界線を含む。
        assert!(d.include_limits);
        // 既定で日付変更線分割を行う。
        assert!(d.split_antimeridian);
    }

    /// `PathOptions` は `Copy`（move されず複製で渡る）。
    #[test]
    fn path_options_is_copy() {
        fn assert_copy<T: Copy>(_: T) {}
        let opts = PathOptions::default();
        let opts2 = opts; // Copy なら move されない
        assert_copy(opts2);
        assert_eq!(opts, opts2); // 元の束縛 opts も有効
    }

    /// `PathOptions` の `PartialEq` は値の同一性を区別する（各フィールドを 1 つずつ
    /// 変えて ne）。PartialEq を常時 true 化する変異・一部フィールドのみ比較する変異を
    /// 撃破する。include_limits と split_antimeridian は別個に区別される。
    #[test]
    fn path_options_partial_eq_distinguishes_each_field() {
        let base = PathOptions::default();
        // 同値は等しい。
        assert_eq!(base, PathOptions::default());

        // sample_interval_seconds だけ異なる。
        let a = PathOptions {
            sample_interval_seconds: 30.0,
            ..PathOptions::default()
        };
        assert_ne!(base, a);

        // include_limits だけ異なる（split_antimeridian は既定の true のまま）。
        let b = PathOptions {
            include_limits: false,
            ..PathOptions::default()
        };
        assert_ne!(base, b);

        // split_antimeridian だけ異なる（include_limits は既定の true のまま）。
        let c = PathOptions {
            split_antimeridian: false,
            ..PathOptions::default()
        };
        assert_ne!(base, c);
    }
}
