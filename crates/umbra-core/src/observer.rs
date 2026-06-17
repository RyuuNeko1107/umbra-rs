//! 観測者の入力型と距離・測地座標の newtype（`docs/api-draft.md` §1.5）。
//!
//! 距離 [`Meters`]/[`Kilometers`]、測地緯度 [`GeodeticLatitude`]（[-90°,90°] = [-π/2,π/2]）、
//! 東経 [`EastLongitude`]（東を正・[-π,π) に正規化）、観測地点 [`Observer`] を提供する。
//! 角度の内部表現はラジアン（`docs/conventions.md` §2）。範囲外の緯度は
//! [`DomainError::OutOfRange`]、経度は正規化（エラーなし）。
//!
//! 緯度の範囲検証は**度域**で行い（`90.0 <= 90.0`）、度→ラジアン変換の ULP で境界が
//! 弾かれないようにする。

use core::f64::consts::FRAC_PI_2;

use crate::angle::{Degrees, Radians};
use crate::error::DomainError;

/// 距離 \[m\]。順序付き（高度比較など, api-draft §1.1）。
#[derive(Clone, Copy, Debug, PartialEq, PartialOrd)]
pub struct Meters(pub f64);

/// 距離 \[km\]。順序付き（経路幅・高度比較など, api-draft §1.1）。
#[derive(Clone, Copy, Debug, PartialEq, PartialOrd)]
pub struct Kilometers(pub f64);

impl Meters {
    /// km へ変換（/1000）。
    pub fn to_kilometers(self) -> Kilometers {
        Kilometers(self.0 / 1000.0)
    }
}

impl Kilometers {
    /// m へ変換（×1000）。
    pub fn to_meters(self) -> Meters {
        Meters(self.0 * 1000.0)
    }
}

/// 測地緯度（[-90°, 90°] = [-π/2, π/2]）。内部表現はラジアン。
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GeodeticLatitude(Radians);

impl GeodeticLatitude {
    /// 度から構築。範囲 [-90, 90]（境界包含）外は [`DomainError::OutOfRange`]。
    pub fn from_degrees(deg: f64) -> Result<Self, DomainError> {
        if !(-90.0..=90.0).contains(&deg) {
            return Err(DomainError::OutOfRange {
                what: "geodetic latitude",
            });
        }
        Ok(Self(Degrees::new(deg).to_radians()))
    }

    /// ラジアンから構築。範囲 [-π/2, π/2]（境界包含）外は [`DomainError::OutOfRange`]。
    pub fn from_radians(rad: Radians) -> Result<Self, DomainError> {
        if !(-FRAC_PI_2..=FRAC_PI_2).contains(&rad.0) {
            return Err(DomainError::OutOfRange {
                what: "geodetic latitude",
            });
        }
        Ok(Self(rad))
    }

    /// ラジアン値。
    pub fn radians(&self) -> Radians {
        self.0
    }

    /// 度値。
    pub fn degrees(&self) -> Degrees {
        self.0.to_degrees()
    }
}

/// 東経（東を正）。内部は [-π, π) に正規化（conventions §2 `normalized_signed`）。
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EastLongitude(Radians);

impl EastLongitude {
    /// 度から構築（[-180, 180) へ正規化。例 190°→-170°, 360°→0°, 180°→-180°）。
    ///
    /// 任意の度は正規化で必ず [-π,π) に収まるため **infallible**（api-draft §1.2 の旧 `Result`
    /// 形・`from_signed_degrees` は本正規化に集約・改訂。西経入力は負の東経として吸収, conventions §3）。
    pub fn from_degrees(deg: f64) -> Self {
        Self::from_radians(Degrees::new(deg).to_radians())
    }

    /// ラジアンから構築（[-π, π) へ正規化）。
    pub fn from_radians(rad: Radians) -> Self {
        Self(rad.normalized_signed())
    }

    /// ラジアン値（[-π, π)）。
    pub fn radians(&self) -> Radians {
        self.0
    }

    /// 度値（[-180, 180)）。
    pub fn degrees(&self) -> Degrees {
        self.0.to_degrees()
    }
}

/// 観測地点（測地緯度・東経・楕円体高）。
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Observer {
    /// 測地緯度。
    pub latitude: GeodeticLatitude,
    /// 東経（東を正）。
    pub longitude: EastLongitude,
    /// 楕円体高 \[m\]（負も可＝海面下）。
    pub elevation: Meters,
}

impl Observer {
    /// 構成済みの newtype から観測地点を構築する。
    pub fn new(latitude: GeodeticLatitude, longitude: EastLongitude, elevation: Meters) -> Self {
        Self {
            latitude,
            longitude,
            elevation,
        }
    }

    /// 度・メートルから構築する。緯度範囲外は [`DomainError::OutOfRange`]。
    /// 経度は正規化、高さは検証なし（負も可）。
    pub fn from_degrees(lat_deg: f64, lon_deg: f64, elevation_m: f64) -> Result<Self, DomainError> {
        Ok(Self::new(
            GeodeticLatitude::from_degrees(lat_deg)?,
            EastLongitude::from_degrees(lon_deg),
            Meters(elevation_m),
        ))
    }
}

#[cfg(test)]
mod tests {
    use crate::angle::{Degrees, Radians};
    use crate::error::DomainError;
    use crate::observer::{EastLongitude, GeodeticLatitude, Kilometers, Meters, Observer};
    use core::f64::consts::PI;

    /// 度の exact 比較許容（正規化結果の桁落ち程度）。
    const TOL_DEG: f64 = 1e-9;
    /// ラジアンの厳密関係許容。
    const TOL_RAD: f64 = 1e-12;

    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    // ============================================================
    // Meters / Kilometers 相互変換
    // ============================================================

    /// `Meters::to_kilometers` は m/1000。例 1500 m → 1.5 km。
    /// 係数（1000 vs 他）と演算方向（/ vs *）の取り違えを撃破する exact オラクル。
    #[test]
    fn meters_to_kilometers_divides_by_thousand() {
        assert_eq!(Meters(1500.0).to_kilometers(), Kilometers(1.5));
        // ゼロは符号・係数に依らず 0（恒等の取り違え検出の補助）。
        assert_eq!(Meters(0.0).to_kilometers(), Kilometers(0.0));
    }

    /// `Kilometers::to_meters` は km*1000。例 2.5 km → 2500 m。
    /// `*1000` を `/1000` へ取り違える変異を撃破する exact オラクル。
    #[test]
    fn kilometers_to_meters_multiplies_by_thousand() {
        assert_eq!(Kilometers(2.5).to_meters(), Meters(2500.0));
    }

    /// 往復で値が一致（m→km→m）。係数の片側だけ誤る変異（例 一方が ×100）を撃破。
    #[test]
    fn meters_kilometers_round_trip_is_identity() {
        let m = Meters(1234.5);
        assert_eq!(m.to_kilometers().to_meters(), m);
    }

    // ============================================================
    // GeodeticLatitude
    // ============================================================

    /// 内部点（岡山 ≈ 34.66°）で from_degrees が Ok・radians/degrees round-trip。
    /// `radians()` は度→ラジアン換算と一致、`degrees()` は入力度に戻る。
    /// 度⇄ラジアン換算係数（180/π）の取り違えを撃破する。
    #[test]
    fn geodetic_latitude_from_degrees_round_trips_interior_point() {
        let lat = GeodeticLatitude::from_degrees(34.66).expect("34.66° は有効");
        assert!(
            close(lat.radians().0, 34.66_f64.to_radians(), TOL_RAD),
            "radians = {}",
            lat.radians().0
        );
        assert!(
            close(lat.degrees().0, 34.66, TOL_DEG),
            "degrees = {}",
            lat.degrees().0
        );
    }

    /// 負の内部点でも対称に動作（南緯 -34.66°）。符号の取り違え（abs 化など）を撃破。
    #[test]
    fn geodetic_latitude_handles_negative_interior_point() {
        let lat = GeodeticLatitude::from_degrees(-34.66).expect("-34.66° は有効");
        assert!(
            close(lat.degrees().0, -34.66, TOL_DEG),
            "degrees = {}",
            lat.degrees().0
        );
        assert!(close(lat.radians().0, (-34.66_f64).to_radians(), TOL_RAD));
    }

    /// 境界 +90°/-90° は**包含**（Ok）。`degrees()` ≈ ±90、`radians()` ≈ ±π/2。
    /// 境界判定を `<` にする（端を除外する）変異を撃破する。
    #[test]
    fn geodetic_latitude_includes_degree_boundaries() {
        let north = GeodeticLatitude::from_degrees(90.0).expect("+90° は包含");
        assert!(
            close(north.degrees().0, 90.0, TOL_DEG),
            "degrees = {}",
            north.degrees().0
        );
        assert!(close(north.radians().0, PI / 2.0, TOL_RAD));

        let south = GeodeticLatitude::from_degrees(-90.0).expect("-90° は包含");
        assert!(
            close(south.degrees().0, -90.0, TOL_DEG),
            "degrees = {}",
            south.degrees().0
        );
        assert!(close(south.radians().0, -PI / 2.0, TOL_RAD));
    }

    /// 範囲外の度はちょうど境界を超えた値（90.0001）で OutOfRange。
    /// variant と `what` 文字列も検証し、境界を `<=` から外す方向の変異・what 文字列改変を撃破。
    #[test]
    fn geodetic_latitude_rejects_just_above_north_boundary() {
        let err = GeodeticLatitude::from_degrees(90.0001).expect_err("90.0001° は範囲外");
        assert_eq!(
            err,
            DomainError::OutOfRange {
                what: "geodetic latitude"
            }
        );
    }

    /// 範囲外の度（91, -91）で OutOfRange。両端の範囲チェック脱落（片側のみ検査）を撃破。
    #[test]
    fn geodetic_latitude_rejects_out_of_range_degrees() {
        assert_eq!(
            GeodeticLatitude::from_degrees(91.0).expect_err("91° は範囲外"),
            DomainError::OutOfRange {
                what: "geodetic latitude"
            }
        );
        assert_eq!(
            GeodeticLatitude::from_degrees(-91.0).expect_err("-91° は範囲外"),
            DomainError::OutOfRange {
                what: "geodetic latitude"
            }
        );
    }

    /// NaN 入力は範囲チェックで弾かれる（`(-90..=90).contains(NaN)==false`）→ OutOfRange。
    /// 浮動小数の連鎖で NaN が混入した場合にサイレント伝播せず Err になる設計をコードで固定する。
    /// from_radians も同様（NaN ラジアン → Err）。
    #[test]
    fn geodetic_latitude_rejects_nan() {
        assert_eq!(
            GeodeticLatitude::from_degrees(f64::NAN).expect_err("NaN 度は範囲外"),
            DomainError::OutOfRange {
                what: "geodetic latitude"
            }
        );
        assert_eq!(
            GeodeticLatitude::from_radians(Radians(f64::NAN)).expect_err("NaN rad は範囲外"),
            DomainError::OutOfRange {
                what: "geodetic latitude"
            }
        );
    }

    /// `from_radians` の内部点は round-trip し、度換算も一致する。
    /// from_radians 経路でも換算係数（180/π）が正しいことを縛る。
    #[test]
    fn geodetic_latitude_from_radians_round_trips() {
        let r = Radians(0.5);
        let lat = GeodeticLatitude::from_radians(r).expect("0.5 rad は有効");
        assert!(
            close(lat.radians().0, 0.5, TOL_RAD),
            "radians = {}",
            lat.radians().0
        );
        assert!(close(
            lat.degrees().0,
            Degrees(0.5_f64 * 180.0 / PI).0,
            TOL_DEG
        ));
    }

    /// `from_radians` の境界 ±π/2 は包含（Ok）。ラジアン経路での境界除外変異を撃破。
    #[test]
    fn geodetic_latitude_from_radians_includes_boundaries() {
        let north = GeodeticLatitude::from_radians(Radians(PI / 2.0)).expect("+π/2 は包含");
        assert!(close(north.radians().0, PI / 2.0, TOL_RAD));
        let south = GeodeticLatitude::from_radians(Radians(-PI / 2.0)).expect("-π/2 は包含");
        assert!(close(south.radians().0, -PI / 2.0, TOL_RAD));
    }

    /// `from_radians` の範囲外（±π/2 を僅かに超える）で OutOfRange。
    /// ラジアン境界の向き・what 文字列を縛る。
    #[test]
    fn geodetic_latitude_from_radians_rejects_out_of_range() {
        assert_eq!(
            GeodeticLatitude::from_radians(Radians(PI / 2.0 + 1e-6)).expect_err("π/2 超過は範囲外"),
            DomainError::OutOfRange {
                what: "geodetic latitude"
            }
        );
        assert_eq!(
            GeodeticLatitude::from_radians(Radians(-PI / 2.0 - 1e-6))
                .expect_err("-π/2 未満は範囲外"),
            DomainError::OutOfRange {
                what: "geodetic latitude"
            }
        );
    }

    // ============================================================
    // EastLongitude（[-180,180) / [-π,π) 正規化）
    // ============================================================

    /// 正規化の代表ケースを degrees() で exact 確認。
    /// 190→-170, -190→170, 360→0, 180→-180（上端は下端へ折り返す半開区間）, 0→0, 270→-90。
    /// 正規化の有無・区間の半開性（180 を含めてしまう変異）を撃破する。
    #[test]
    fn east_longitude_normalizes_degrees_to_half_open_interval() {
        let cases = [
            (190.0, -170.0),
            (-190.0, 170.0),
            (360.0, 0.0),
            (180.0, -180.0),
            (0.0, 0.0),
            (270.0, -90.0),
        ];
        for (input, expected) in cases {
            let lon = EastLongitude::from_degrees(input);
            assert!(
                close(lon.degrees().0, expected, TOL_DEG),
                "from_degrees({input}).degrees() = {} (期待 {expected})",
                lon.degrees().0
            );
        }
    }

    /// 内部点（東経 133.92° 岡山）は正規化で不変。範囲内の値を誤って動かす変異を撃破。
    #[test]
    fn east_longitude_keeps_interior_value() {
        let lon = EastLongitude::from_degrees(133.92);
        assert!(
            close(lon.degrees().0, 133.92, TOL_DEG),
            "degrees = {}",
            lon.degrees().0
        );
    }

    /// `radians()` は常に [-π,π)：下端 -π を含み上端 π を含まない。
    /// 180°→-π（含む）/ 360°→0 で半開区間を縛り、`degrees()` は radians 由来で一致する。
    #[test]
    fn east_longitude_radians_in_half_open_interval() {
        let at_lower = EastLongitude::from_degrees(180.0);
        assert!(
            close(at_lower.radians().0, -PI, TOL_RAD),
            "180° の radians = {}",
            at_lower.radians().0
        );
        // 任意ケースで radians ∈ [-π, π)。
        for input in [190.0_f64, -190.0, 360.0, 270.0, 0.0, 133.92] {
            let r = EastLongitude::from_degrees(input).radians().0;
            assert!(
                (-PI..PI).contains(&r),
                "from_degrees({input}).radians() = {r}"
            );
        }
        // degrees() は radians 由来（換算係数 180/π の一致）。
        let lon = EastLongitude::from_degrees(270.0);
        assert!(close(
            lon.degrees().0,
            lon.radians().to_degrees().0,
            TOL_DEG
        ));
    }

    /// `from_radians` も [-π,π) へ正規化。3π/2 → -π/2、ちょうど π → -π（下端）。
    /// ラジアン経路の正規化脱落・上端の取り扱いを撃破する。
    #[test]
    fn east_longitude_from_radians_normalizes() {
        let a = EastLongitude::from_radians(Radians(3.0 * PI / 2.0));
        assert!(
            close(a.radians().0, -PI / 2.0, TOL_RAD),
            "radians = {}",
            a.radians().0
        );

        let b = EastLongitude::from_radians(Radians(PI));
        assert!(
            close(b.radians().0, -PI, TOL_RAD),
            "radians = {}",
            b.radians().0
        );

        let c = EastLongitude::from_radians(Radians(-PI / 4.0));
        assert!(
            close(c.radians().0, -PI / 4.0, TOL_RAD),
            "radians = {}",
            c.radians().0
        );
    }

    // ============================================================
    // Observer
    // ============================================================

    /// `Observer::new` はフィールドをそのまま保持（変換・正規化を挟まない）。
    /// フィールドの取り違え（lat/lon 入れ替え等）を撃破する。
    #[test]
    fn observer_new_preserves_fields() {
        let lat = GeodeticLatitude::from_degrees(34.66).unwrap();
        let lon = EastLongitude::from_degrees(133.92);
        let elev = Meters(50.0);
        let obs = Observer::new(lat, lon, elev);
        assert_eq!(obs.latitude, lat);
        assert_eq!(obs.longitude, lon);
        assert_eq!(obs.elevation, elev);
    }

    /// `from_degrees` 有効入力（岡山 ≈ 北緯34.66°,東経133.92°,高 50 m）→ Ok・値保持。
    /// latitude.degrees()≈lat、longitude≈lon（範囲内なので不変）、elevation==Meters(elev)。
    /// 度→型構築の取り違え・elevation の係数誤りを撃破する。
    #[test]
    fn observer_from_degrees_valid_preserves_values() {
        let obs = Observer::from_degrees(34.66, 133.92, 50.0).expect("岡山は有効");
        assert!(close(obs.latitude.degrees().0, 34.66, TOL_DEG));
        assert!(close(obs.longitude.degrees().0, 133.92, TOL_DEG));
        assert_eq!(obs.elevation, Meters(50.0));
    }

    /// `from_degrees` の緯度範囲外（91°）→ Err(OutOfRange{what:"geodetic latitude"})。
    /// 緯度検証の脱落・what 文字列改変を撃破する。
    #[test]
    fn observer_from_degrees_rejects_invalid_latitude() {
        let err = Observer::from_degrees(91.0, 0.0, 0.0).expect_err("緯度91°は範囲外");
        assert_eq!(
            err,
            DomainError::OutOfRange {
                what: "geodetic latitude"
            }
        );
    }

    /// `from_degrees` は経度を正規化する（190° → -170°）。経度の正規化脱落を撃破する。
    #[test]
    fn observer_from_degrees_normalizes_longitude() {
        let obs = Observer::from_degrees(0.0, 190.0, 0.0).expect("経度は正規化され有効");
        assert!(
            close(obs.longitude.degrees().0, -170.0, TOL_DEG),
            "longitude = {}",
            obs.longitude.degrees().0
        );
    }

    /// `from_degrees` の elevation は検証されず**負も可**（-100 m）、そのまま保持。
    /// elevation に範囲チェックを入れる/符号を落とす変異を撃破する。
    #[test]
    fn observer_from_degrees_allows_negative_elevation() {
        let obs = Observer::from_degrees(0.0, 0.0, -100.0).expect("負の高さも可");
        assert_eq!(obs.elevation, Meters(-100.0));
    }

    // ============================================================
    // newtype の derive 特性（Copy / PartialEq）
    // ============================================================

    /// 各 newtype は `Copy`（move されず複製で渡る）。`#[derive(Copy)]` の脱落を
    /// コンパイル時に撃破する（消費後も元の束縛が有効）。
    #[test]
    fn newtypes_are_copy() {
        fn assert_copy<T: Copy>(_: T) {}

        let m = Meters(1.0);
        let m2 = m;
        assert_copy(m2);
        assert_eq!(m, m2);

        let k = Kilometers(1.0);
        let k2 = k;
        assert_copy(k2);
        assert_eq!(k, k2);

        let lat = GeodeticLatitude::from_degrees(10.0).unwrap();
        let lat2 = lat;
        assert_copy(lat2);
        assert_eq!(lat, lat2);

        let lon = EastLongitude::from_degrees(10.0);
        let lon2 = lon;
        assert_copy(lon2);
        assert_eq!(lon, lon2);

        let obs = Observer::from_degrees(10.0, 20.0, 30.0).unwrap();
        let obs2 = obs;
        assert_copy(obs2);
        assert_eq!(obs, obs2);
    }

    /// `PartialEq` は値の同一性を区別する（異なる値は不一致）。
    /// `#[derive(PartialEq)]` を常時 true 化する変異を撃破する。
    #[test]
    fn newtypes_partial_eq_distinguishes_values() {
        assert_ne!(Meters(1.0), Meters(2.0));
        assert_ne!(Kilometers(1.0), Kilometers(2.0));
        assert_ne!(
            GeodeticLatitude::from_degrees(10.0).unwrap(),
            GeodeticLatitude::from_degrees(20.0).unwrap()
        );
        assert_ne!(
            EastLongitude::from_degrees(10.0),
            EastLongitude::from_degrees(20.0)
        );
        // Observer は全フィールド一致でのみ等しい（経度だけ違うペア）。
        let a = Observer::from_degrees(10.0, 20.0, 30.0).unwrap();
        let b = Observer::from_degrees(10.0, 21.0, 30.0).unwrap();
        assert_ne!(a, b);
        assert_eq!(a, a);
    }
}
