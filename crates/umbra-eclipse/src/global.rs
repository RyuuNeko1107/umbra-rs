//! 全球の日食種別判定（`docs/issues/ISSUE-023`、`docs/physical-models.md` §C3）。
//!
//! 影軸の地心最小距離 `gamma`（Re）と本影半径 `l2`（符号付き）で分類する（Meeus Ch.54 基準）:
//! ```text
//! |gamma| < 0.9972                    → 中心食（l2<0 皆既 / l2>0 金環）
//! 0.9972 ≤ |gamma| < 0.9972 + |l2|    → 非中心 皆既/金環
//! 0.9972 + |l2| ≤ |gamma| < 1.5433+l2 → 部分食
//! |gamma| ≥ 1.5433 + l2               → 日食なし
//! ```
//! 注: ハイブリッド（中心線上で l2 が符号反転）は単一時刻では判定不能。全球パス（時系列）で
//! 判定し本関数は瞬時の Total/Annular を返す（要確認: 0.9972/1.5433 の式番号・有効桁。§C3）。

use crate::besselian::BesselianElements;

/// 中心食境界（≈1 − 扁平縮約。Meeus）。
const CENTRAL_LIMIT: f64 = 0.9972;
/// 半影限界（Meeus）。
const PENUMBRA_LIMIT: f64 = 1.5433;

/// 太陽食の種別。
#[non_exhaustive]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SolarEclipseKind {
    /// 部分食。
    Partial,
    /// 金環食。
    Annular,
    /// 皆既食。
    Total,
    /// ハイブリッド（金環↔皆既。全球パスで判定）。
    Hybrid,
    /// 非中心の皆既。
    NonCentralTotal,
    /// 非中心の金環。
    NonCentralAnnular,
}

/// 瞬時ベッセル要素から日食種別を判定する。`None` は（その時刻に全球で）日食なし。
///
/// ハイブリッドは返さない（時系列が必要。上記注）。中心/非中心は l2 符号で皆既/金環を分ける。
pub fn classify(elements: &BesselianElements) -> Option<SolarEclipseKind> {
    let g = elements.gamma(); // ≥ 0
    let l2 = elements.l2;
    let total = l2 < 0.0; // l2<0 皆既 / l2>0 金環（正本 B1）

    if g < CENTRAL_LIMIT {
        Some(if total {
            SolarEclipseKind::Total
        } else {
            SolarEclipseKind::Annular
        })
    } else if g < CENTRAL_LIMIT + l2.abs() {
        Some(if total {
            SolarEclipseKind::NonCentralTotal
        } else {
            SolarEclipseKind::NonCentralAnnular
        })
    } else if g < PENUMBRA_LIMIT + l2 {
        Some(SolarEclipseKind::Partial)
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::SolarEclipseKind::{Annular, NonCentralAnnular, NonCentralTotal, Partial, Total};
    use super::*;
    use umbra_core::constants::{EARTH_EQUATORIAL_RADIUS_M, SOLAR_RADIUS_KM};
    use umbra_core::{JulianDate2, Radians, TdbInstant};
    use umbra_ephemeris::{Body, Ephemeris, EphemerisFrame, MockEphemeris, Origin};

    /// gamma=`g`（x=g, y=0）, 本影半径 `l2` のベッセル要素を作る。
    fn elem(g: f64, l2: f64) -> BesselianElements {
        BesselianElements {
            x: g,
            y: 0.0,
            declination: Radians(0.0),
            l1: 0.53,
            l2,
            tan_f1: 0.0047,
            tan_f2: 0.0046,
        }
    }

    #[test]
    fn central_total_and_annular_by_l2_sign() {
        assert_eq!(classify(&elem(0.5, -0.02)), Some(Total));
        assert_eq!(classify(&elem(0.5, 0.02)), Some(Annular));
    }

    #[test]
    fn non_central_band_by_l2_sign() {
        // 0.9972 ≤ g < 0.9972+|l2|（|l2|=0.02 → 上限 1.0172）。
        assert_eq!(classify(&elem(1.0, -0.02)), Some(NonCentralTotal));
        assert_eq!(classify(&elem(1.0, 0.02)), Some(NonCentralAnnular));
    }

    #[test]
    fn partial_band() {
        assert_eq!(classify(&elem(1.2, 0.01)), Some(Partial));
    }

    #[test]
    fn no_eclipse_when_gamma_too_large() {
        assert_eq!(classify(&elem(2.0, 0.01)), None);
    }

    #[test]
    fn central_to_noncentral_boundary() {
        // g=0.9972 ちょうどは中心食でない（< 厳密）→ 非中心。直下は中心。
        assert_eq!(classify(&elem(0.9972, -0.02)), Some(NonCentralTotal));
        assert_eq!(classify(&elem(0.9971, -0.02)), Some(Total));
    }

    #[test]
    fn noncentral_to_partial_boundary() {
        // 境界は実装と同じ計算 CENTRAL_LIMIT+|l2| で踏む（リテラル 1.0172 では f64 が
        // ぴったり一致せず < / <= を区別できない）。境界ちょうどは非中心でない → 部分。
        let b = CENTRAL_LIMIT + 0.02;
        assert_eq!(classify(&elem(b, -0.02)), Some(Partial));
        assert_eq!(classify(&elem(b - 1e-6, -0.02)), Some(NonCentralTotal));
    }

    #[test]
    fn l2_exactly_zero_is_annular_not_total() {
        // l2==0（皆既/金環の連続境界）は total=(l2<0)=false → 金環側に倒す（< 厳密）。
        assert_eq!(classify(&elem(0.5, 0.0)), Some(Annular));
    }

    #[test]
    fn partial_to_none_boundary() {
        // g=1.5433+l2 ちょうどは日食なし → 直下は部分。
        assert_eq!(classify(&elem(1.5433 + 0.01, 0.01)), None);
        assert_eq!(classify(&elem(1.5433 + 0.01 - 1e-6, 0.01)), Some(Partial));
    }

    #[test]
    fn partial_to_none_boundary_negative_l2() {
        // 皆既側(l2<0)では上限が 1.5433+l2 < 1.5433 に縮む（符号付き）。
        // `+l2` を `+|l2|` や `+0.01` に取り違えるとこのテストで露見する（H1）。
        let l2 = -0.02;
        assert_eq!(classify(&elem(1.5433 + l2 - 1e-6, l2)), Some(Partial));
        assert_eq!(classify(&elem(1.5433 + l2 + 1e-6, l2)), None);
    }

    #[test]
    fn partial_and_none_bands_with_negative_l2() {
        assert_eq!(classify(&elem(1.2, -0.02)), Some(Partial));
        assert_eq!(classify(&elem(2.0, -0.02)), None);
    }

    #[test]
    fn matches_mock_configurations() {
        let t = TdbInstant::from_jd2(JulianDate2::from_jd(2_451_545.0));
        let r_sun = SOLAR_RADIUS_KM;
        let r_moon = 0.2725076 * (EARTH_EQUATORIAL_RADIUS_M / 1000.0);
        let kind = |m: &MockEphemeris| {
            let pos = |b| {
                m.state(b, t, Origin::Geocenter, EphemerisFrame::Icrs)
                    .unwrap()
                    .position
            };
            let e = crate::besselian::besselian_elements(
                pos(Body::Sun),
                pos(Body::Moon),
                r_sun,
                r_moon,
            )
            .unwrap();
            classify(&e)
        };
        assert_eq!(kind(&MockEphemeris::central_total()), Some(Total));
        assert_eq!(kind(&MockEphemeris::clear_annular()), Some(Annular));
        assert_eq!(kind(&MockEphemeris::clear_partial()), Some(Partial));
        assert_eq!(kind(&MockEphemeris::shadow_misses_earth()), None);
        // 非中心皆既（暦→ベッセル→分類の貫通で NonCentralTotal バンドを踏む, M2）。
        assert_eq!(
            kind(&MockEphemeris::non_central_total()),
            Some(NonCentralTotal)
        );
    }
}
