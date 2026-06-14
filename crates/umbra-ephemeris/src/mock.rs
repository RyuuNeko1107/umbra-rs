//! 人工配置の天体暦（`docs/issues/ISSUE-038`、accuracy.md §3.1）。
//!
//! DE 無し・純解析配置で日食幾何ロジック（影円錐・ベッセル要素・分類・局地）を検証する足場。
//! 各コンストラクタは**最大食付近の静的な太陽・月の地心位置**を返す。意図する (gamma, 種別) は
//! 各コンストラクタに記す。**l2 の符号・gamma の値そのものの検証は、幾何層（umbra-eclipse,
//! ISSUE-019/020/021/023）が構築されてから**そちらのテストで行う（本 Issue は位置の契約を検証）。

use crate::ephemeris::{
    Body, Ephemeris, EphemerisError, EphemerisFrame, EphemerisMetadata, Origin, StateVector,
};
use umbra_core::{JulianDate2, TdbInstant, TimeRange, Vector3};

const AU_KM: f64 = 149_597_870.7;
const MOON_PERIGEE_KM: f64 = 357_000.0;
const MOON_APOGEE_KM: f64 = 406_000.0;
const MOON_MEAN_KM: f64 = 384_400.0;
/// 地球質量 / 月質量。EMB 位置 = 月位置 / (1 + 本比)。
const EARTH_MOON_MASS_RATIO: f64 = 81.300_56;

/// 人工配置の天体暦（テスト専用）。太陽・月の地心位置を静的に保持する。
#[derive(Debug, Clone, Copy)]
pub struct MockEphemeris {
    sun: Vector3,
    moon: Vector3,
}

impl MockEphemeris {
    /// 完全皆既（影軸が地心を貫く中心配置 gamma≈0、近地点で月が大きく total ⇒ 意図 l2<0）。
    pub fn central_total() -> Self {
        MockEphemeris {
            sun: Vector3::new(AU_KM, 0.0, 0.0),
            moon: Vector3::new(MOON_PERIGEE_KM, 0.0, 0.0),
        }
    }

    /// 明確な金環（中心配置 gamma≈0、遠地点で月が小さく annular ⇒ 意図 l2>0）。
    pub fn clear_annular() -> Self {
        MockEphemeris {
            sun: Vector3::new(AU_KM, 0.0, 0.0),
            moon: Vector3::new(MOON_APOGEE_KM, 0.0, 0.0),
        }
    }

    /// 明確な部分食（影軸が地球縁を外す ⇒ |gamma|≈1.1 で中心食不成立・部分食）。
    pub fn clear_partial() -> Self {
        MockEphemeris {
            sun: Vector3::new(AU_KM, 0.0, 0.0),
            moon: Vector3::new(MOON_MEAN_KM, 7_000.0, 0.0),
        }
    }

    /// 非中心の皆既（近地点＋地球縁すれすれのオフセット ⇒ |gamma|≈1.002、本影が縁に接触）。
    pub fn non_central_total() -> Self {
        MockEphemeris {
            sun: Vector3::new(AU_KM, 0.0, 0.0),
            moon: Vector3::new(MOON_PERIGEE_KM, 6_390.0, 0.0),
        }
    }

    /// 日食なし（影軸を地球から大きく外す ⇒ 意図 |gamma| > 1 + l1）。
    pub fn shadow_misses_earth() -> Self {
        MockEphemeris {
            sun: Vector3::new(AU_KM, 0.0, 0.0),
            moon: Vector3::new(MOON_MEAN_KM, 40_000.0, 0.0),
        }
    }
}

impl Ephemeris for MockEphemeris {
    /// 地心位置を返す（`origin`/`frame` は無視する単純配置。`time` 非依存の静的配置）。
    fn state(
        &self,
        body: Body,
        _time: TdbInstant,
        _origin: Origin,
        _frame: EphemerisFrame,
    ) -> Result<StateVector, EphemerisError> {
        let position = match body {
            Body::Earth => Vector3::ZERO,
            Body::Sun => self.sun,
            Body::Moon => self.moon,
            Body::EarthMoonBarycenter => self.moon.scale(1.0 / (1.0 + EARTH_MOON_MASS_RATIO)),
        };
        Ok(StateVector {
            position,
            velocity: None,
        })
    }

    fn supported_range(&self) -> TimeRange<TdbInstant> {
        TimeRange {
            start: TdbInstant::from_jd2(JulianDate2::from_jd(2_415_020.0)), // ≈1900
            end: TdbInstant::from_jd2(JulianDate2::from_jd(2_488_070.0)),   // ≈2100
        }
    }

    fn metadata(&self) -> EphemerisMetadata {
        EphemerisMetadata {
            model: "MockEphemeris".to_string(),
            version: "artificial-configuration".to_string(),
            source: "synthetic (test-only)".to_string(),
            license: "test-only, not for production".to_string(),
            supported: self.supported_range(),
            max_residual_arcsec: f64::NAN, // 精度申告対象外（幾何ロジック検証用）
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tdb() -> TdbInstant {
        TdbInstant::from_jd2(JulianDate2::from_jd(2_451_545.0))
    }

    fn state_of(m: &MockEphemeris, body: Body) -> Vector3 {
        m.state(body, tdb(), Origin::Geocenter, EphemerisFrame::Icrs)
            .unwrap()
            .position
    }

    #[test]
    fn earth_is_geocentric_origin() {
        let m = MockEphemeris::central_total();
        assert_eq!(state_of(&m, Body::Earth), Vector3::ZERO);
    }

    #[test]
    fn sun_and_moon_positions_match_configuration() {
        let m = MockEphemeris::central_total();
        assert_eq!(state_of(&m, Body::Sun), Vector3::new(AU_KM, 0.0, 0.0));
        assert_eq!(
            state_of(&m, Body::Moon),
            Vector3::new(MOON_PERIGEE_KM, 0.0, 0.0)
        );
    }

    #[test]
    fn configurations_differ_as_intended() {
        // total は近地点（小距離）、annular は遠地点（大距離）。
        let total_moon = state_of(&MockEphemeris::central_total(), Body::Moon).x;
        let annular_moon = state_of(&MockEphemeris::clear_annular(), Body::Moon).x;
        assert!(total_moon < annular_moon, "{total_moon} < {annular_moon}");
        // partial は小オフセット、miss は大オフセット（軸ずれ）。
        let partial_off = state_of(&MockEphemeris::clear_partial(), Body::Moon).y;
        let miss_off = state_of(&MockEphemeris::shadow_misses_earth(), Body::Moon).y;
        assert!(partial_off > 0.0 && miss_off > partial_off);
    }

    #[test]
    fn emb_lies_between_earth_and_moon() {
        let m = MockEphemeris::central_total();
        let emb = state_of(&m, Body::EarthMoonBarycenter);
        let expected = MOON_PERIGEE_KM / (1.0 + EARTH_MOON_MASS_RATIO);
        assert!((emb.x - expected).abs() < 1e-9, "emb.x = {}", emb.x);
        // 地心と月の間（0 < emb < moon）。
        assert!(emb.x > 0.0 && emb.x < MOON_PERIGEE_KM);
    }

    #[test]
    fn velocity_is_none_for_static_mock() {
        let m = MockEphemeris::central_total();
        let s = m
            .state(Body::Sun, tdb(), Origin::Geocenter, EphemerisFrame::Icrs)
            .unwrap();
        assert!(s.velocity.is_none());
    }

    #[test]
    fn metadata_marks_test_only() {
        let m = MockEphemeris::central_total();
        let meta = m.metadata();
        assert_eq!(meta.model, "MockEphemeris");
        assert!(meta.max_residual_arcsec.is_nan());
    }

    #[test]
    fn supported_range_covers_2020() {
        let m = MockEphemeris::central_total();
        let r = m.supported_range();
        let y2020 = JulianDate2::from_jd(2_458_850.0).jd();
        assert!(r.start.jd2().jd() < y2020 && y2020 < r.end.jd2().jd());
    }

    #[test]
    fn is_send_and_sync() {
        fn assert_send_sync<T: Send + Sync>() {}
        assert_send_sync::<MockEphemeris>();
    }
}
