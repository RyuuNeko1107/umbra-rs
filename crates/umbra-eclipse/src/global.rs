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

use crate::axis_intercept::shadow_axis_surface_point;
use crate::besselian::BesselianElements;
use crate::config::EngineConfig;
use crate::conjunction::RootConfig;
use crate::error::EclipseError;
use crate::global_contacts::{global_contact_ground_point, solve_global_contacts};
use crate::horizontal::{sun_horizontal, RefractionModel};
use crate::local_maximum::solve_local_maximum;
use crate::magnitude::{eclipse_magnitude, eclipse_obscuration};
use crate::projection::project_observer_to_fundamental;
use crate::results::GreatestEclipse;
use crate::source::BesselianSource;
use umbra_core::deltat::DeltaTModel;
use umbra_core::ellipsoid::{observer_geocentric, Ellipsoid, GeocentricObserver};
use umbra_core::{JulianDate2, Radians, SolverError, TtInstant};

/// 1 日 = 86400 SI 秒（root_tolerance を日へ換算）。
const SECONDS_PER_DAY: f64 = 86_400.0;
/// 最大食時刻 Brent 求根の反復上限。
const GREATEST_ROOT_MAX_ITER: usize = 200;

/// 中心食境界（≈1 − 扁平縮約。Meeus）。
const CENTRAL_LIMIT: f64 = 0.9972;
/// 半影限界（Meeus）。
const PENUMBRA_LIMIT: f64 = 1.5433;

/// 太陽食の種別。
#[non_exhaustive]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize))]
#[cfg_attr(feature = "serde", serde(tag = "type"))]
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

/// 最大食（global greatest eclipse）の解（時刻・地表点・食分/食面積/太陽高度）と gamma。
///
/// `kind` と全球接触 P1/U1/U4/P4・帯幅・中心食継続は S6b、`SolarEclipse` 組立は S6c の責務。
/// `path_width`/`central_duration` は本スライスでは常に `None`（S6b で充足）。
#[derive(Clone, Copy, Debug, PartialEq)]
#[allow(dead_code)] // S6c（classify_global / search 結線）が消費するまで未使用。
pub(crate) struct GreatestEclipseSolution {
    /// 最大食の地表点・食分・食面積・太陽高度・時刻（path/duration は None）。
    pub greatest: GreatestEclipse,
    /// 影軸の地心最小距離 gamma（Re, 符号なし）。
    pub gamma: f64,
}

/// 中心食の最大食状況（時刻・gamma・地表点・食分/食面積・太陽高度）を解く（ISSUE-043 S6a-ii）。
///
/// 検証済みプリミティブの合成: 最大食時刻と gamma は地心観測者（ρ=0）に対する
/// [`solve_local_maximum`] で得る（m²=(0−x)²+(0−y)²=x²+y²=gamma², ISSUE-026）。地表点は
/// [`shadow_axis_surface_point`]（S6a-i）、太陽高度は [`sun_horizontal`]（ISSUE-028）。食分・食面積は
/// 観測者 ζ で補正した本影/半影半径 L1'=l1−ζ·tanf1, L2'=l2−ζ·tanf2 から
/// [`eclipse_magnitude`]/[`eclipse_obscuration`] で評価する（中心点では m=0・視半径比 ρ=(L1'−L2')/(L1'+L2')）。
///
/// 中心食でない（軸が地表を外す）と地表点が無く `Err(Solver(RootNotBracketed))`（S6b で部分/非中心を扱う）。
#[allow(dead_code)] // S6c が消費するまで未使用。
pub(crate) fn solve_greatest_eclipse<B, M>(
    source: &B,
    delta_t: &M,
    config: &EngineConfig,
) -> Result<GreatestEclipseSolution, EclipseError>
where
    B: BesselianSource,
    M: DeltaTModel,
{
    // earth_model は現状 WGS84 のみ（projection/horizontal/axis_intercept と同様に定数で扱う）。
    let _ = config.earth_model;
    let ellipsoid = Ellipsoid::WGS84;

    // 1. 最大食時刻・gamma: 地心観測者（ρsin=ρcos=0）の局地最大食で得る。射影は ξ=η=ζ=0 となり
    //    m² = (0−x)² + (0−y)² = x² + y² = gamma²（影軸の地心距離²）。経度は ρcos=0 ゆえ無関係。
    let geocenter = GeocentricObserver {
        rho_sin_phi_prime: 0.0,
        rho_cos_phi_prime: 0.0,
    };
    let root_config = RootConfig {
        x_tolerance_days: config.root_tolerance_seconds / SECONDS_PER_DAY,
        max_iterations: GREATEST_ROOT_MAX_ITER,
    };
    let max = solve_local_maximum(
        source,
        &geocenter,
        Radians::new(0.0),
        source.fit_interval(),
        root_config,
    )?;
    let gamma = max.min_separation; // = √(x²+y²) at t_max

    // 2. 最大食時刻の瞬時ベッセル要素。
    let elements = source.at(max.time_tt)?;

    // 3. 最大食地点・観測者の基本面距離 m・ζ を中心食/部分食で分岐して得る。
    //    中心食: 影軸が地表を貫く（`shadow_axis_surface_point` 成功）→ 観測者は軸上で m=0、ζ は
    //      地表点を検証済み前方射影し直して取得（高さ補正用）。
    //    部分/非中心: 軸が地表を外す（RootNotBracketed）→ 地球縁で軸に最も近い点（縁点・ζ=0）が
    //      最大食地点。観測者の基本面距離 m = ρ_axis − ρ_g = gamma − 1（縁点が軸に最も近い）。
    //      **m は max(0) でクランプ**: 非中心帯（扁平楕円体では軸が外れるが球近似 gamma<1）では
    //      gamma−1<0 となるが、観測者-軸距離は非負・かつ中心点(m=0)が最大食分の上限。負 m を許すと
    //      中心値を超える非物理 magnitude になるため 0 に頭打ちする（球/扁平モデル差の吸収・要確認帯）。
    let (position, m, zeta) = match shadow_axis_surface_point(&elements, &ellipsoid) {
        Ok(p) => {
            let obs = observer_geocentric(&ellipsoid, p.lat.radians().0, 0.0);
            let zeta = project_observer_to_fundamental(&obs, p.lon.radians(), &elements).zeta;
            (p, 0.0, zeta)
        }
        Err(EclipseError::Solver(SolverError::RootNotBracketed)) => {
            let p = global_contact_ground_point(&elements, &ellipsoid)?;
            (p, (gamma - 1.0).max(0.0), 0.0)
        }
        Err(other) => return Err(other),
    };

    // 4. 食分・食面積（観測者 ζ で補正した半径 L1'=l1−ζ·tanf1, L2'=l2−ζ·tanf2）。
    //    食分 magnitude = (L1'−m)/(L1'+L2')。視半径比 ρ=(L1'−L2')/(L1'+L2')、視半径平面の中心離隔
    //    separation = (1+ρ)·m/L1'（m=0→0=同心, m=L1'→1+ρ=外接）。中心食は m=0 で従来と一致。
    let l1p = elements.l1 - zeta * elements.tan_f1;
    let l2p = elements.l2 - zeta * elements.tan_f2;
    let magnitude = eclipse_magnitude(m, l1p, l2p);
    let radius_ratio = (l1p - l2p) / (l1p + l2p);
    let separation = (1.0 + radius_ratio) * m / l1p;
    let obscuration = eclipse_obscuration(separation, 1.0, radius_ratio);

    // 5. 太陽の幾何学的高度（大気差なし, conventions §7 既定）。
    let sun_altitude = sun_horizontal(
        position.lat.radians(),
        position.lon.radians(),
        max.time_tt,
        RefractionModel::None,
        delta_t,
    )
    .altitude_geometric;

    let greatest = GreatestEclipse {
        time_utc: max.time_utc,
        time_tt: max.time_tt,
        position,
        magnitude,
        obscuration,
        path_width: None,
        central_duration: None,
        sun_altitude,
    };
    Ok(GreatestEclipseSolution { greatest, gamma })
}

/// 全球の日食種別（Total/Annular/Hybrid/Partial/NonCentral）を時系列込みで判定する（ISSUE-043 S6b-iii）。
///
/// 瞬時 [`classify`] が gamma＋最大食時 l2 符号で基本種別（Total/Annular/Partial/NonCentral 系）を返すが、
/// **Hybrid（中心線上で金環⇄皆既が切替わる）は単一時刻では判定不能**。本関数は最大食
/// （地心観測者 ρ=0 の [`solve_local_maximum`]）で基本種別を取り、中心食（Total/Annular）なら全球
/// 中心食区間 [U1,U4]（[`solve_global_contacts`]）で l2 の符号反転を走査して Hybrid を上書きする。
/// `None` は探索窓に日食なし。
#[allow(dead_code)] // S6c（classify_global / search 結線）が消費するまで未使用。
pub(crate) fn classify_global_kind<B: BesselianSource>(
    source: &B,
    config: RootConfig,
) -> Result<Option<SolarEclipseKind>, EclipseError> {
    // 1. 最大食（gamma 最小）時刻: 地心観測者 ρ=0 の局地最大食（投影 ξ=η=0 ⇒ m²=x²+y²=gamma²）。
    let geocenter = GeocentricObserver {
        rho_sin_phi_prime: 0.0,
        rho_cos_phi_prime: 0.0,
    };
    let max = solve_local_maximum(
        source,
        &geocenter,
        Radians::new(0.0),
        source.fit_interval(),
        config,
    )?;

    // 2. 最大食時の瞬時要素から基本種別（classify は瞬時 gamma＋l2 符号）。日食なしは None。
    let e = source.at(max.time_tt)?;
    let base = classify(&BesselianElements {
        x: e.x,
        y: e.y,
        declination: e.declination,
        l1: e.l1,
        l2: e.l2,
        tan_f1: e.tan_f1,
        tan_f2: e.tan_f2,
    });
    let Some(base) = base else { return Ok(None) };

    // 3. Hybrid 上書き: 中心食（Total/Annular）で全球中心食区間 [U1,U4] の l2 が符号反転なら Hybrid。
    //    部分食・非中心は中心食区間が無く Hybrid 対象外（base をそのまま返す）。
    if matches!(base, SolarEclipseKind::Total | SolarEclipseKind::Annular) {
        let contacts = solve_global_contacts(source, config)?;
        if let (Some(u1), Some(u4)) = (contacts.u1, contacts.u4) {
            if l2_changes_sign(source, u1.time_tt, u4.time_tt)? {
                return Ok(Some(SolarEclipseKind::Hybrid));
            }
        }
    }
    Ok(Some(base))
}

/// 区間 `[start_tt, end_tt]` で本影半径 l2 が符号反転する（正と負の両方を取る）かを粗走査で判定する。
///
/// Hybrid（中心線上で金環⇄皆既が切替わる）の検出に用いる。**サンプル数は load-bearing**: 細かさが
/// 足りないと、両端が同符号で中央だけ反対符号（金環-皆既-金環など）の短い切替帯を取りこぼし偽陰性に
/// なる（`scan_point_count` 系の純解像度＝等価とは異なる）。ハイブリッドの中心線皆既/金環区間（数分〜
/// オーダー）を確実に捉えるため十分大きく取る。符号判定 `l2>0`/`l2<0` も load-bearing。
#[allow(dead_code)]
#[allow(clippy::cast_precision_loss)]
fn l2_changes_sign<B: BesselianSource>(
    source: &B,
    start_tt: TtInstant,
    end_tt: TtInstant,
) -> Result<bool, EclipseError> {
    /// 中心食区間の l2 符号走査の分割数。5 時間級の [U1,U4] でも刻み ~1 分で短い切替帯を捉える。
    const SAMPLES: usize = 256;
    let t0 = start_tt.jd2().jd();
    let t1 = end_tt.jd2().jd();
    let mut saw_positive = false;
    let mut saw_negative = false;
    for i in 0..=SAMPLES {
        let frac = i as f64 / SAMPLES as f64;
        let jd = t0 + (t1 - t0) * frac;
        let l2 = source.at(TtInstant::from_jd2(JulianDate2::from_jd(jd)))?.l2;
        if l2 > 0.0 {
            saw_positive = true;
        }
        if l2 < 0.0 {
            saw_negative = true;
        }
    }
    Ok(saw_positive && saw_negative)
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

    // ====================================================================
    // solve_greatest_eclipse（ISSUE-043 S6a-ii・全球最大食組立）
    // ====================================================================
    //
    // ## オラクル戦略（追認回避）
    // 主オラクルは **独立内部再計算**。solver の内部手法（local_maximum・逆射影・magnitude/
    // obscuration の合成）には依存せず、外部観測可能な振る舞いを **検証済プリミティブ**
    // `source.at(t)`（ISSUE-037）＋ `project_observer_to_fundamental`（ISSUE-024）＋
    // `observer_geocentric`（ISSUE-010/011）から **別経路** で縛る:
    //   - gamma 再計算: 返った time_tt で √(x²+y²) を `source.at` から直接組む（1e-6 一致）。
    //   - 局地最小性: time_tt±δ の gamma が time_tt 以上（局地最小）。300s 両側で狭義増加。
    //   - 地表点往復: 返った position を前方射影し直すと (ξ,η)=(e.x,e.y)・ζ>0（逆射影の独立検証）。
    //   - 食面積==1（皆既の強い縛り）: 中心点では太陽が完全に隠れる。
    // NASA 公表値（gamma≈0.4367・greatest≈18:25UTC・37.0N/87.7W・magnitude≈1.031・alt 61–64°）は
    // **ballpark のみ**（k/ΔT 慣習差で秒値は再現しない）。range check に限定して用いる。
    //
    // 探索窓は `source.fit_interval()`。2017-08-21 最大食（TT-JD≈2457986.768）を内部に括る
    // 区間 [2457986, 2457988] を渡す（時刻 solver がブラケットできるよう最小が窓内部に来る）。
    //
    // 注（追補・実装レビューの結線網羅）: 当初省略した 2 経路を後段で追加した（実装が確定し
    // 金環分岐・非中心分岐が end-to-end で踏めるようになったため）:
    //   - 金環（2023-10-14）: l2>0 ⇒ magnitude<1 / obscuration<1 / radius_ratio<1 の分岐を踏む
    //     （`greatest_annular_2023_*`）。NASA 値は ballpark のみ・solver が窓内の真の最小を探す。
    //   - 軸が地表を外す→RootNotBracketed: 時不変 `ConstantSource` は時不変ゆえ「軸ミス」ではなく
    //     `solve_local_maximum` のブラケット不成立で別理由 Err になる。代わりに **時変** 合成供給源
    //     （gamma に内部極小を持つが極小値が >1）を使い、`solve_local_maximum` は成功し
    //     `shadow_axis_surface_point` 側で RootNotBracketed が起きることを縛る（`greatest_*_axis_miss_*`）。

    use umbra_core::ellipsoid::{observer_geocentric, Ellipsoid};
    use umbra_core::{EspenakMeeusDeltaT, TimeInterval, TtInstant};

    use crate::besselian::InstantaneousBesselianElements;
    use crate::projection::project_observer_to_fundamental;
    use crate::source::{BesselianSource, DirectBesselianSource};

    /// 1 日 = 86400 SI 秒。
    const SECONDS_PER_DAY: f64 = 86_400.0;
    /// 太陽物理半径[km]（local_maximum.rs / axis_intercept.rs と同一）。
    const G_R_SUN: f64 = SOLAR_RADIUS_KM;
    /// 月半径[km]（k·Re, IAU 慣習 k=0.2725076・同上）。
    const G_R_MOON: f64 = 0.2725076 * (EARTH_EQUATORIAL_RADIUS_M / 1000.0);

    /// 2017-08-21 最大食を内部に括る探索窓（TT-JD 2457986〜2457988, axis_intercept.rs と同形）。
    /// 最大食 TT-JD≈2457986.768 が区間内部にあり、時刻 solver が最小をブラケットできる。
    fn solve_window_2017() -> TimeInterval<TtInstant> {
        TimeInterval {
            start: TtInstant::from_jd2(JulianDate2::new(2_457_986.0, 0.0)),
            end: TtInstant::from_jd2(JulianDate2::new(2_457_988.0, 0.0)),
        }
    }

    /// この供給源の TT-JD（単一 f64）から TtInstant。
    fn g_tt_jd(jd: f64) -> TtInstant {
        TtInstant::from_jd2(JulianDate2::from_jd(jd))
    }

    /// 2017-08-21 中心皆既を括る `DirectBesselianSource`（fit_interval=探索窓）を解く。
    fn solve_2017<'d>(
        dt: &'d EspenakMeeusDeltaT,
    ) -> (
        DirectBesselianSource<'d, EspenakMeeusDeltaT>,
        GreatestEclipseSolution,
    ) {
        let src = DirectBesselianSource::new(G_R_SUN, G_R_MOON, dt, solve_window_2017());
        let config = crate::config::EngineConfig::standard();
        let sol = solve_greatest_eclipse(&src, dt, &config)
            .expect("2017 central total eclipse should yield a greatest-eclipse solution");
        (src, sol)
    }

    /// 影軸交点 (x,y) から geocentric な gamma=√(x²+y²) を `source.at` 由来で独立再計算する。
    fn gamma_at<B: BesselianSource>(src: &B, t: TtInstant) -> f64 {
        let e = src
            .at(t)
            .expect("source.at should succeed near 2017 eclipse");
        (e.x * e.x + e.y * e.y).sqrt()
    }

    /// 観点1（gamma ballpark, 緩め・NASA range check）: gamma ∈ [0.40, 0.47]。
    /// NASA 公表 gamma≈0.4367 を k/ΔT 慣習差の余裕込みで括る（絶対基準にしない）。
    #[test]
    fn greatest_gamma_in_nasa_ballpark() {
        let dt = EspenakMeeusDeltaT;
        let (_src, sol) = solve_2017(&dt);
        assert!(
            (0.40..=0.47).contains(&sol.gamma),
            "gamma={} not in NASA ballpark [0.40,0.47] (NASA≈0.4367)",
            sol.gamma
        );
    }

    /// 観点2（gamma 独立再計算, tight）: 返った gamma が、返った time_tt で `source.at` から
    /// 別経路に組んだ √(x²+y²) と 1e-6 Re 一致する。gamma が「その時刻の geocentric 軸距離」で
    /// あることを直接縛る（追認回避の主オラクル）。
    #[test]
    fn greatest_gamma_matches_independent_recomputation() {
        let dt = EspenakMeeusDeltaT;
        let (src, sol) = solve_2017(&dt);
        let g = gamma_at(&src, sol.greatest.time_tt);
        assert!(
            (sol.gamma - g).abs() < 1e-6,
            "gamma={} must match independent √(x²+y²)={g} at time_tt (tol 1e-6 Re)",
            sol.gamma
        );
    }

    /// 観点3（局地最小性, tight・主オラクル）: 返った time_tt で gamma が局地最小である。
    /// time_tt±60s / ±300s の gamma（`source.at` から独立再計算）が time_tt の gamma 以上で、
    /// ±300s では両側で狭義増加する（平坦/最大取り違えを弾く）。
    #[test]
    fn greatest_time_is_local_minimum_of_geocentric_gamma() {
        let dt = EspenakMeeusDeltaT;
        let (src, sol) = solve_2017(&dt);
        let t0 = sol.greatest.time_tt.jd2().jd();
        let g0 = gamma_at(&src, sol.greatest.time_tt);

        for &delta in &[60.0_f64, 300.0] {
            let plus = gamma_at(&src, g_tt_jd(t0 + delta / SECONDS_PER_DAY));
            let minus = gamma_at(&src, g_tt_jd(t0 - delta / SECONDS_PER_DAY));
            assert!(
                plus >= g0,
                "gamma(t+{delta}s)={plus} must be ≥ gamma(t_max)={g0} (local min)"
            );
            assert!(
                minus >= g0,
                "gamma(t−{delta}s)={minus} must be ≥ gamma(t_max)={g0} (local min)"
            );
        }
        // 300s 両側で狭義増加（谷であること＝最小/最大取り違えを弾く）。
        let p300 = gamma_at(&src, g_tt_jd(t0 + 300.0 / SECONDS_PER_DAY));
        let m300 = gamma_at(&src, g_tt_jd(t0 - 300.0 / SECONDS_PER_DAY));
        assert!(
            p300 > g0 && m300 > g0,
            "300s away on both sides must strictly increase: +300s={p300}, −300s={m300}, center={g0}"
        );
    }

    /// 観点4（TT/UTC 整合, tight）: time_utc == tt_to_utc(time_tt)（2017 は post-1972 で変換可能）。
    #[test]
    fn greatest_utc_matches_tt_to_utc() {
        let dt = EspenakMeeusDeltaT;
        let (_src, sol) = solve_2017(&dt);
        let want = umbra_core::time::tt_to_utc(sol.greatest.time_tt)
            .expect("2017 is post-1972, tt_to_utc must succeed");
        assert_eq!(
            sol.greatest.time_utc, want,
            "time_utc must equal tt_to_utc(time_tt)"
        );
    }

    /// 観点5（UTC 時刻 ballpark, 緩め）: time_utc の時刻が 17.5〜18.8 時（北米中緯度の最大食は
    /// 18時台 UTC）。NASA 公表 ≈18:25 UTC を range check で括る。
    #[test]
    fn greatest_utc_hour_in_ballpark() {
        let dt = EspenakMeeusDeltaT;
        let (_src, sol) = solve_2017(&dt);
        let (_y, _mo, _d, hh, mm, _ss) = sol.greatest.time_utc.to_gregorian();
        let hour = f64::from(hh) + f64::from(mm) / 60.0;
        assert!(
            (17.5..=18.8).contains(&hour),
            "greatest UTC hour {hour} not ~18:xx (NASA≈18:25 UTC)"
        );
    }

    /// 観点6（地表点 ballpark, 緩め）: 最大食地点が lat∈[30,42]N・lon∈[−95,−80]E。
    /// NASA 公表 ≈37.0N/87.7W を range check で括る。
    #[test]
    fn greatest_position_in_nasa_ballpark() {
        let dt = EspenakMeeusDeltaT;
        let (_src, sol) = solve_2017(&dt);
        let lat = sol.greatest.position.lat.degrees().0;
        let lon = sol.greatest.position.lon.degrees().0;
        assert!(
            (30.0..=42.0).contains(&lat),
            "greatest lat {lat}° not in ballpark [30,42]N (NASA≈37.0N)"
        );
        assert!(
            (-95.0..=-80.0).contains(&lon),
            "greatest lon {lon}°E not in ballpark [−95,−80] (NASA≈−87.7)"
        );
    }

    /// 観点7（地表点 往復, tight・主オラクル）: 返った position を **検証済み前方射影** へ通すと
    /// (ξ,η)=(e.x,e.y)・ζ>0 に往復一致する（影軸が太陽側で地表を貫く点であることの独立検証）。
    /// 逆射影の内部式は再実装せず、`project_observer_to_fundamental`（ISSUE-024）＋
    /// `observer_geocentric`（ISSUE-010/011）の独立経路で縛る（axis_intercept.rs と同方針）。
    #[test]
    fn greatest_position_roundtrips_through_forward_projection() {
        let dt = EspenakMeeusDeltaT;
        let (src, sol) = solve_2017(&dt);
        let e = src
            .at(sol.greatest.time_tt)
            .expect("source.at should succeed at greatest-eclipse time");

        let phi = sol.greatest.position.lat.radians().0;
        let lam = sol.greatest.position.lon.radians().0;
        let obs = observer_geocentric(&Ellipsoid::WGS84, phi, 0.0);
        let r = project_observer_to_fundamental(&obs, Radians::new(lam), &e);
        assert!(
            (r.xi - e.x).abs() < 1e-9,
            "ξ={} must round-trip to e.x={} (tol 1e-9)",
            r.xi,
            e.x
        );
        assert!(
            (r.eta - e.y).abs() < 1e-9,
            "η={} must round-trip to e.y={} (tol 1e-9)",
            r.eta,
            e.y
        );
        assert!(r.zeta > 0.0, "ζ={} must be sunward (>0)", r.zeta);
    }

    /// 観点8（magnitude, total）: 中心点の食分 > 1（皆既）かつ ballpark < 1.08。
    /// NASA 公表 magnitude≈1.031 を range check で括る（皆既は厳密に 1 超でなければならない）。
    #[test]
    fn greatest_magnitude_exceeds_one_for_total() {
        let dt = EspenakMeeusDeltaT;
        let (_src, sol) = solve_2017(&dt);
        let mag = sol.greatest.magnitude.0;
        assert!(mag > 1.0, "total eclipse magnitude {mag} must exceed 1");
        assert!(
            mag < 1.08,
            "total eclipse magnitude {mag} above ballpark <1.08 (NASA≈1.031)"
        );
    }

    /// 観点9（obscuration==1, total・強い縛り）: 皆既の中心点では太陽面が完全に隠れる。
    /// obscuration ≈ 1.0 を 1e-9 で縛る（中心点で太陽が月に完全内包される幾何の直接帰結）。
    #[test]
    fn greatest_obscuration_is_one_for_total() {
        let dt = EspenakMeeusDeltaT;
        let (_src, sol) = solve_2017(&dt);
        let obsc = sol.greatest.obscuration.0;
        assert!(
            (obsc - 1.0).abs() < 1e-9,
            "total eclipse central obscuration {obsc} must be exactly 1.0 (tol 1e-9)"
        );
    }

    /// 観点9b（ζ補正半径の結線ピン, tight）: 返った magnitude/obscuration が、**検証済み前方射影**
    /// 由来の ζ で **正しい符号** の補正半径 L1'=l1−ζ·tanf1・L2'=l2−ζ·tanf2 から組んだ値に厳密一致する。
    /// ballpark（mag<1.08）・obscuration==1（皆既で ζ 補正に鈍感）では捕捉できない、ζ 補正項の
    /// 符号・演算の取り違え（`l1−ζ·tanf1` → `+`/`/` 等）を撃破する結線ピン（S5b の半径結線ピンと同方針）。
    /// minus 符号は錐頂点が観測者より太陽側（apex ζ > 観測者 ζ）で半径が観測者へ向け縮む幾何に由来
    /// （独立レビューで besselian.rs apex 定義から確認済）。ζ は impl と同じ前方射影法だが本テストは
    /// **正符号を明示** するため、impl が符号を誤れば不一致で落ちる。
    #[test]
    fn greatest_magnitude_obscuration_match_zeta_corrected_radii() {
        let dt = EspenakMeeusDeltaT;
        let (src, sol) = solve_2017(&dt);
        let e = src
            .at(sol.greatest.time_tt)
            .expect("source.at should succeed at greatest-eclipse time");
        // ζ を検証済み前方射影から独立に取得（地表点を射影し直す）。
        let phi = sol.greatest.position.lat.radians().0;
        let lam = sol.greatest.position.lon.radians().0;
        let obs = observer_geocentric(&Ellipsoid::WGS84, phi, 0.0);
        let zeta = project_observer_to_fundamental(&obs, Radians::new(lam), &e).zeta;
        // 正しい符号（minus）の ζ 補正半径。
        let l1p = e.l1 - zeta * e.tan_f1;
        let l2p = e.l2 - zeta * e.tan_f2;
        let want_mag = crate::magnitude::eclipse_magnitude(0.0, l1p, l2p);
        let ratio = (l1p - l2p) / (l1p + l2p);
        let want_obsc = crate::magnitude::eclipse_obscuration(0.0, 1.0, ratio);
        assert_eq!(
            sol.greatest.magnitude, want_mag,
            "magnitude は ζ補正半径(minus)由来でなければならない"
        );
        assert_eq!(
            sol.greatest.obscuration, want_obsc,
            "obscuration は ζ補正半径(minus)由来でなければならない"
        );
    }

    /// 観点10（太陽高度 ballpark, 緩め）: 最大食地点の太陽幾何高度 ∈ [55,70]°。
    /// NASA 公表 ≈61–64° を range check で括る。
    #[test]
    fn greatest_sun_altitude_in_ballpark() {
        let dt = EspenakMeeusDeltaT;
        let (_src, sol) = solve_2017(&dt);
        let alt = sol.greatest.sun_altitude.0;
        assert!(
            (55.0..=70.0).contains(&alt),
            "greatest sun altitude {alt}° not in ballpark [55,70] (NASA≈61–64)"
        );
    }

    /// 観点11（path/duration は本スライス非責務）: path_width・central_duration はともに None。
    /// S6b（帯幅・中心食継続）が充足するまで本関数では常に None を返す契約を縛る。
    #[test]
    fn greatest_path_and_duration_are_none() {
        let dt = EspenakMeeusDeltaT;
        let (_src, sol) = solve_2017(&dt);
        assert_eq!(
            sol.greatest.path_width, None,
            "path_width must be None in this slice (S6b territory)"
        );
        assert_eq!(
            sol.greatest.central_duration, None,
            "central_duration must be None in this slice (S6b territory)"
        );
    }

    // ====================================================================
    // 追補A: 金環食 end-to-end（l2>0 ⇒ magnitude<1 / obscuration<1 / radius_ratio<1 分岐）
    // ====================================================================
    //
    // 2017 は皆既（l2<0・magnitude>1・obscuration==1）のみで、金環分岐
    // （L2'>0 ⇒ magnitude<1, 視半径比 ρ<1）が end-to-end で踏まれていなかった。
    // 2023-10-14 金環日食（greatest≈17:59:28 UTC, NASA ballpark: gamma≈0.375・
    // magnitude≈0.952・obscuration≈0.91・≈11.4°N/83.1°W）を実 epoch で踏む。
    // NASA 公表値は **ballpark のみ**（k/ΔT 慣習差で秒値・小数下位は再現しない）。

    /// 2023-10-14 金環食 greatest（TT-JD≈2_460_232.25）を内部に括る 1.5 日探索窓。
    /// gamma の極小が区間内部にあり、時刻 solver がブラケットできる。
    fn solve_window_2023() -> TimeInterval<TtInstant> {
        TimeInterval {
            start: TtInstant::from_jd2(JulianDate2::new(2_460_231.5, 0.0)),
            end: TtInstant::from_jd2(JulianDate2::new(2_460_233.0, 0.0)),
        }
    }

    /// 2023-10-14 金環食を括る `DirectBesselianSource`（半径・config は 2017 と同一）を解く。
    fn solve_2023<'d>(
        dt: &'d EspenakMeeusDeltaT,
    ) -> (
        DirectBesselianSource<'d, EspenakMeeusDeltaT>,
        GreatestEclipseSolution,
    ) {
        let src = DirectBesselianSource::new(G_R_SUN, G_R_MOON, dt, solve_window_2023());
        let config = crate::config::EngineConfig::standard();
        let sol = solve_greatest_eclipse(&src, dt, &config)
            .expect("2023 annular eclipse should yield a greatest-eclipse solution");
        (src, sol)
    }

    /// 追補A（金環 end-to-end）: 2023-10-14 金環食を解き、金環の契約と独立再計算を縛る。
    ///
    /// **本テストの主眼**は金環分岐（l2>0 ⇒ L2'>0）の網羅: 中心点でも太陽は完全には隠れず
    /// `magnitude<1`・`obscuration<1`・視半径比 ρ=(L1'−L2')/(L1'+L2')<1 を踏む（皆既と区別）。
    /// NASA 値（gamma≈0.375・≈11.4°N/83.1°W・≈17:59 UTC・mag≈0.952・obsc≈0.91）は **ballpark のみ**。
    /// 主オラクルは 2017 と同じ独立再計算（gamma=√(x²+y²)・前方射影往復・TT/UTC 整合）。
    #[test]
    fn greatest_annular_2023_contract_and_independent_checks() {
        let dt = EspenakMeeusDeltaT;
        let (src, sol) = solve_2023(&dt);

        // --- 金環の契約（本テストの主眼・l2>0 / radius_ratio<1 分岐を踏む）---
        let mag = sol.greatest.magnitude.0;
        let obsc = sol.greatest.obscuration.0;
        assert!(
            mag < 1.0,
            "annular central magnitude {mag} must be < 1 (金環は太陽を覆い切らない)"
        );
        assert!(
            mag > 0.85,
            "annular central magnitude {mag} below ballpark >0.85 (NASA≈0.952)"
        );
        assert!(
            obsc < 1.0,
            "annular central obscuration {obsc} must be < 1 (環が残る)"
        );
        assert!(
            obsc > 0.7,
            "annular central obscuration {obsc} below ballpark >0.7 (NASA≈0.91)"
        );

        // --- ballpark（NASA range check のみ・絶対基準にしない）---
        assert!(
            (0.30..=0.45).contains(&sol.gamma),
            "gamma={} not in NASA ballpark [0.30,0.45] (NASA≈0.375)",
            sol.gamma
        );
        let lat = sol.greatest.position.lat.degrees().0;
        let lon = sol.greatest.position.lon.degrees().0;
        assert!(
            (0.0..=22.0).contains(&lat),
            "annular lat {lat}° not in ballpark [0,22]N (NASA≈11.4N)"
        );
        assert!(
            (-92.0..=-74.0).contains(&lon),
            "annular lon {lon}°E not in ballpark [−92,−74] (NASA≈−83.1)"
        );
        let (_y, _mo, _d, hh, mm, _ss) = sol.greatest.time_utc.to_gregorian();
        let hour = f64::from(hh) + f64::from(mm) / 60.0;
        assert!(
            (17.0..=19.0).contains(&hour),
            "annular greatest UTC hour {hour} not ~18:xx (NASA≈17:59 UTC)"
        );

        // --- 独立再計算（tight・2017 と同じ主オラクル）---
        // gamma == √(x²+y²) at time_tt（source.at 由来の別経路, tol 1e-6 Re）。
        let g = gamma_at(&src, sol.greatest.time_tt);
        assert!(
            (sol.gamma - g).abs() < 1e-6,
            "gamma={} must match independent √(x²+y²)={g} at time_tt (tol 1e-6)",
            sol.gamma
        );

        // 地表点往復: 前方射影し直すと (ξ,η)=(e.x,e.y)・ζ>0。
        let e = src
            .at(sol.greatest.time_tt)
            .expect("source.at should succeed at 2023 greatest-eclipse time");
        let phi = sol.greatest.position.lat.radians().0;
        let lam = sol.greatest.position.lon.radians().0;
        let obs = observer_geocentric(&Ellipsoid::WGS84, phi, 0.0);
        let r = project_observer_to_fundamental(&obs, Radians::new(lam), &e);
        assert!(
            (r.xi - e.x).abs() < 1e-9,
            "ξ={} must round-trip to e.x={}",
            r.xi,
            e.x
        );
        assert!(
            (r.eta - e.y).abs() < 1e-9,
            "η={} must round-trip to e.y={}",
            r.eta,
            e.y
        );
        assert!(r.zeta > 0.0, "ζ={} must be sunward (>0)", r.zeta);

        // TT/UTC 整合。
        let want_utc = umbra_core::time::tt_to_utc(sol.greatest.time_tt)
            .expect("2023 is post-1972, tt_to_utc must succeed");
        assert_eq!(
            sol.greatest.time_utc, want_utc,
            "time_utc must equal tt_to_utc(time_tt)"
        );

        // path/duration は本スライス非責務（常に None）。
        assert_eq!(
            sol.greatest.path_width, None,
            "path_width must be None (S6b)"
        );
        assert_eq!(
            sol.greatest.central_duration, None,
            "central_duration must be None (S6b)"
        );
    }

    // ====================================================================
    // 追補B: 非中心/部分食（軸が地表を外す）→ S6c-i で Ok（縁点・部分食最大食）になった
    // ====================================================================
    //
    // 旧テスト `greatest_noncentral_axis_miss_propagates_root_not_bracketed`（軸ミス→
    // `Err(Solver(RootNotBracketed))`）は、S6c-i で `solve_greatest_eclipse` が部分食/非中心も
    // 扱う（軸ミス時は地球縁点・m=gamma−1）よう拡張されたため陳腐化し削除した。同じ「軸ミス・
    // 内部 gamma 極小 >1」構成は `greatest_partial_synthetic_*`（SyntheticGammaSource x_min=1.2）が
    // Ok（部分食最大食）として縛る。

    // ====================================================================
    // 追補C: 部分/非中心（軸が地表を外す）の最大食 → Ok(GreatestEclipseSolution)
    // ====================================================================
    //
    // S6c-i で `solve_greatest_eclipse` を **部分/非中心** へ拡張する。影軸が地表を外す
    // （`shadow_axis_surface_point` が解なし）場合、従来は `Err(Solver(RootNotBracketed))`
    // を返していた（追補B でその旧挙動を縛った）。新挙動は **Ok** を返し:
    //   - gamma  = 地心軸距離の最小値（部分は >1）。
    //   - greatest.position = 軸に最も近い地球リム点（sub-axis 終端点・地平線上で太陽が昇る点）。
    //   - magnitude  = eclipse_magnitude(m, l1, l2)（m=gamma−1, リムでは ζ=0 ⇒ 半径は無補正 l1,l2）。
    //   - obscuration = eclipse_obscuration((1+ρ)·m/l1, 1.0, ρ)（ρ=(l1−l2)/(l1+l2)）。
    //   - sun_altitude ≈ 0°（リム点）。
    //   - path_width = None / central_duration = None。
    //
    // ## オラクル戦略（追認回避）
    // 合成幾何そのものがオラクル: 既知の gamma/l1/l2 と文書化された式から magnitude/obscuration を
    // **手計算** して縛る。合成供給源は `SyntheticGammaSource`（x=x_min+K·(jd−center)², y=0, l1/l2 固定）。
    //   synthetic_source(1.2, 50.0, 0.54, −0.01, 0.15) ⇒ gamma_min=1.2（部分: 1<1.2<1.54）。
    //   m = gamma−1 = 0.2、magnitude = (0.54−0.2)/(0.54+(−0.01)) = 0.34/0.53 = 0.6415094…（0<mag<1 部分）。
    //   ρ = (0.54−(−0.01))/(0.54+(−0.01)) = 0.55/0.53、sep = (1+ρ)·m/0.54、obscuration = obsc(sep,1,ρ)∈(0,1)。
    // リム性は **物理的に独立** に縛る: 返った position を `observer_geocentric`＋
    // `project_observer_to_fundamental`（ISSUE-024）へ通し |ζ|<0.02（終端点・太陽が地平線上）。
    // 注: sun_altitude は **実太陽** から導かれ、合成供給源の d/μ は任意値ゆえ ≈0° を assert しない
    // （代わりにリム/終端性を縛る。spec の方針どおり）。

    /// gamma の中心極小値を `x_min`、曲率 `k`、本影/半影半径 `l1`/`l2`、半幅 `half_day`[day] に取る
    /// `SyntheticGammaSource`（部分テスト用・既存の struct 直書きと同形のコンストラクタ）。
    /// 最大食 TT-JD は 2017 中心（2457986.768）と同じ位置に置き、極小が窓内部に来る。
    fn synthetic_source(
        x_min: f64,
        k: f64,
        l1: f64,
        l2: f64,
        half_day: f64,
    ) -> SyntheticGammaSource {
        let center_jd = 2_457_986.768;
        SyntheticGammaSource {
            center_jd,
            x_min,
            k,
            l1,
            l2,
            window: TimeInterval {
                start: g_tt_jd(center_jd - half_day),
                end: g_tt_jd(center_jd + half_day),
            },
        }
    }

    /// 追補C-1（部分最大食・厳密ピン）: 部分食の合成供給源で `solve_greatest_eclipse` が **Ok** を返し、
    /// gamma/magnitude/obscuration/path/duration が部分食の契約と手計算オラクルに厳密一致する。
    ///
    /// gamma_min=1.2（部分: 1<1.2<1.54）。手計算: m=gamma−1=0.2、
    /// magnitude=(0.54−0.2)/(0.54+(−0.01))=0.34/0.53=0.6415094…（0<mag<1）。
    /// ρ=(0.54−(−0.01))/(0.54+(−0.01))、sep=(1+ρ)·m/0.54、obscuration=obsc(sep,1,ρ)∈(0,1)。
    /// 旧実装（軸ミスで Err を返す）に対しては Ok を取れず **red**（Err(Solver(RootNotBracketed))）になる。
    #[test]
    fn greatest_partial_synthetic_contract_and_hand_computed() {
        let dt = EspenakMeeusDeltaT;
        let src = synthetic_source(1.2, 50.0, 0.54, -0.01, 0.15);
        let config = crate::config::EngineConfig::standard();
        let sol = solve_greatest_eclipse(&src, &dt, &config)
            .expect("partial (axis-miss) eclipse must now yield Ok(GreatestEclipseSolution)");

        // gamma 独立再計算: 返った time_tt で √(x²+y²)（source.at 由来）と 1e-6 一致。
        let g = gamma_at(&src, sol.greatest.time_tt);
        assert!(
            (sol.gamma - g).abs() < 1e-6,
            "gamma={} must match independent √(x²+y²)={g} at time_tt (tol 1e-6)",
            sol.gamma
        );
        // gamma は部分食バンド（>1, ballpark ≈1.2）。
        assert!(
            (1.0..1.5433).contains(&sol.gamma),
            "partial gamma {} must be >1 and <penumbra limit (≈1.2)",
            sol.gamma
        );
        assert!(
            (sol.gamma - 1.2).abs() < 1e-3,
            "partial gamma {} should be ≈1.2 (synthetic min, solver tol)",
            sol.gamma
        );

        // --- magnitude: m=gamma−1（返った gamma 由来）で eclipse_magnitude に厳密一致 ---
        let m = sol.gamma - 1.0;
        let want_mag = eclipse_magnitude(m, 0.54, -0.01);
        assert_eq!(
            sol.greatest.magnitude, want_mag,
            "magnitude は m=gamma−1（リムで ζ=0 ⇒ 無補正 l1,l2）由来でなければならない"
        );
        // 手計算 ballpark: m≈0.2 ⇒ magnitude≈0.6415094…、かつ部分食 0<mag<1。
        let mag = sol.greatest.magnitude.0;
        assert!(
            (mag - 0.641_509_433_962_264).abs() < 1e-6,
            "partial magnitude {mag} must be ≈0.6415094 (=(0.54−0.2)/0.53), hand-computed"
        );
        assert!(
            (0.0..1.0).contains(&mag),
            "partial magnitude {mag} must be in (0,1) (太陽を覆い切らない)"
        );

        // --- obscuration: sep=(1+ρ)·m/l1（ρ=(l1−l2)/(l1+l2)）で eclipse_obscuration に厳密一致 ---
        let rho = (0.54 - (-0.01)) / (0.54 + (-0.01));
        let sep = (1.0 + rho) * m / 0.54;
        let want_obsc = eclipse_obscuration(sep, 1.0, rho);
        assert!(
            (sol.greatest.obscuration.0 - want_obsc.0).abs() < 1e-9,
            "obscuration {} must equal eclipse_obscuration((1+ρ)·m/l1,1,ρ)={} (tol 1e-9)",
            sol.greatest.obscuration.0,
            want_obsc.0
        );
        assert!(
            (0.0..1.0).contains(&sol.greatest.obscuration.0),
            "partial obscuration {} must be in (0,1)",
            sol.greatest.obscuration.0
        );

        // --- path/duration は本スライス非責務（常に None）---
        assert!(
            sol.greatest.path_width.is_none(),
            "path_width must be None for partial (S6b territory)"
        );
        assert!(
            sol.greatest.central_duration.is_none(),
            "central_duration must be None for partial (S6b territory)"
        );
    }

    /// 追補C-2（部分の地表点はリム上・物理的独立縛り）: 部分食の greatest.position を前方射影し直すと
    /// 終端点（|ζ|<0.02）に乗る ＝ 軸に最も近い sub-axis リム点（太陽が地平線上）であること。
    ///
    /// 逆射影の内部式は再実装せず、`observer_geocentric`（ISSUE-010/011）＋
    /// `project_observer_to_fundamental`（ISSUE-024）の独立経路で縛る。sun_altitude は実太陽由来で
    /// 合成供給源の d/μ は任意ゆえ ≈0° は assert せず、終端（リム）性のみを縛る（spec 方針）。
    #[test]
    fn greatest_partial_position_is_on_the_limb() {
        let dt = EspenakMeeusDeltaT;
        let src = synthetic_source(1.2, 50.0, 0.54, -0.01, 0.15);
        let config = crate::config::EngineConfig::standard();
        let sol = solve_greatest_eclipse(&src, &dt, &config)
            .expect("partial (axis-miss) eclipse must now yield Ok(GreatestEclipseSolution)");

        let e = src
            .at(sol.greatest.time_tt)
            .expect("source.at should succeed at partial greatest-eclipse time");
        let phi = sol.greatest.position.lat.radians().0;
        let lam = sol.greatest.position.lon.radians().0;
        let obs = observer_geocentric(&Ellipsoid::WGS84, phi, 0.0);
        let r = project_observer_to_fundamental(&obs, Radians::new(lam), &e);
        assert!(
            r.zeta.abs() < 0.02,
            "partial greatest position must lie on the terminator (limb): |ζ|={} < 0.02 \
             (sub-axis limb point, sun on the horizon)",
            r.zeta.abs()
        );
    }

    // 追補C-3（任意・実部分食 end-to-end）は **意図的に省略**する。clean な部分専用の実日食 epoch
    // ＋窓を解析暦/ΔT 慣習差に依らず安定に括るのは borderline（gamma>1 の余裕や greatest 時刻の
    // ばらつきで flaky になりうる）。部分食の堅牢なピンは合成テスト追補C-1（手計算オラクル）に委ねる。

    // ====================================================================
    // 追補C-4: 非中心帯（gamma<1 だが扁平楕円体で軸が外れる）の m クランプ
    // ====================================================================
    //
    // 部分分岐の `m = (gamma−1.0).max(0.0)` の `.max(0.0)` は、**軸が WGS84（扁平）地表を外すが
    // gamma<1** という非中心帯でのみ効く。極方向（扁平で縮んだ極半径 ~0.9966）へ向く軸では
    // gamma<1 でも軸が楕円体を外し（`shadow_axis_surface_point` が `RootNotBracketed`）部分分岐が
    // 取られる。そこでは gamma−1<0 ゆえクランプが m=0 に頭打ちする（観測者-軸距離は非負・中心点が
    // 食分の上限）。クランプなしなら m=gamma−1<0 となり magnitude が中心値を **超える** 非物理値に
    // なる。既存テスト（追補C-1 は gamma=1.2>1 で m>0）はこの帯を踏まず、クランプ（および `.max(0.0)`
    // の `remove`/`→.min` ミュータント）が gamma<1 で未検証だった。本テストがそれを縛る。

    /// 極方向（x=0, y 軸上）へ向く軸を持つ時変合成供給源。gamma=|y| が center で内部極小 Y_MIN を
    /// 取る（`y(jd)=Y_MIN+K·(jd−center)²`）。`shadow_axis_surface_point` の扁平楕円体半径
    /// `r(ζ)=ζ²+(y/(1−f))²−1` は、Y_MIN>1−f（≈0.99665）なら全 ζ で正 ⇒ 根なし ⇒ 部分分岐を踏む
    /// （gamma=Y_MIN<1 のまま）。Y_MIN=0.998 は (0.998/0.996647)²≈1.0027>1 を満たす。
    struct PolarAxisSource {
        center_jd: f64,
        y_min: f64,
        k: f64,
        l1: f64,
        l2: f64,
        window: TimeInterval<TtInstant>,
    }

    impl BesselianSource for PolarAxisSource {
        /// x=0（軸は y 軸上＝極方向）, y=Y_MIN+K·(jd−center)²（gamma=|y| が center で内部極小）, l1/l2 固定。
        fn at(&self, t: TtInstant) -> Result<InstantaneousBesselianElements, EclipseError> {
            let dj = t.jd2().jd() - self.center_jd;
            Ok(InstantaneousBesselianElements {
                x: 0.0,
                y: self.y_min + self.k * dj * dj,
                declination: Radians(0.0),
                mu: Radians(0.0),
                l1: self.l1,
                l2: self.l2,
                tan_f1: 0.0047,
                tan_f2: 0.0046,
                time_tt: t,
            })
        }

        fn fit_interval(&self) -> TimeInterval<TtInstant> {
            self.window
        }
    }

    /// 追補C-4（非中心帯の m クランプ・厳密ピン）: 極方向の軸（gamma<1 だが扁平楕円体で軸が外れる）の
    /// 非中心皆既で `m = (gamma−1).max(0.0)` が m=0 に頭打ちされ、magnitude が中心値（m=0）に一致する。
    ///
    /// 構成: Y_MIN=0.998（>1−f≈0.99665 ⇒ 扁平軸ミス）・K=50・l1=0.54・l2=−0.01（皆既）。center で
    /// gamma=0.998<1（非中心帯）。±0.05 day でも y=0.998+50·0.0025≈1.123 と窓内で内部極小は center に
    /// 立つ。
    ///
    /// 独立にレジームを縛る（追認回避）: (a) 返った greatest 時刻で gamma=√(x²+y²)≈0.998<1（1e-6）、
    /// (b) その時刻の `shadow_axis_surface_point(&elements, WGS84)` が `Err(Solver(RootNotBracketed))`
    /// ＝扁平軸ミスで部分分岐を踏む。クランプのピン: magnitude が **m=0** の `eclipse_magnitude(0,l1,l2)`
    /// ＝l1/(l1+l2)=0.54/0.53≈1.0189 に厳密一致（1e-9）。クランプなしなら m=gamma−1≈−0.002 で
    /// `eclipse_magnitude(−0.002,l1,l2)` となり中心値を **超える**（より大）ため、m=0 値との一致は
    /// クランプ除去 / `.max→.min` ミュータントを殺す。非中心皆既ゆえ magnitude>1・obscuration≈1。
    #[test]
    fn greatest_noncentral_total_clamps_m_to_zero() {
        let dt = EspenakMeeusDeltaT;
        let center_jd = 2_457_986.768;
        let half_day = 0.05;
        let l1 = 0.54;
        let l2 = -0.01;
        let src = PolarAxisSource {
            center_jd,
            y_min: 0.998,
            k: 50.0,
            l1,
            l2,
            window: TimeInterval {
                start: g_tt_jd(center_jd - half_day),
                end: g_tt_jd(center_jd + half_day),
            },
        };
        let config = crate::config::EngineConfig::standard();
        let sol = solve_greatest_eclipse(&src, &dt, &config)
            .expect("non-central (oblate axis-miss, gamma<1) eclipse must yield Ok");

        // (a) レジーム: gamma=√(x²+y²)≈0.998<1（非中心帯）を source.at 由来の別経路で縛る。
        let g = gamma_at(&src, sol.greatest.time_tt);
        assert!(
            (sol.gamma - g).abs() < 1e-6,
            "gamma={} must match independent √(x²+y²)={g} (tol 1e-6)",
            sol.gamma
        );
        assert!(
            (sol.gamma - 0.998).abs() < 1e-6,
            "gamma {} must be ≈0.998 (non-central band, <1)",
            sol.gamma
        );
        assert!(
            sol.gamma < 1.0,
            "gamma {} must be <1 (non-central band)",
            sol.gamma
        );

        // (b) レジーム: 扁平楕円体で軸が外れる（部分分岐＝RootNotBracketed）。gamma<1 でも軸ミス。
        let e = src
            .at(sol.greatest.time_tt)
            .expect("source.at should succeed at greatest-eclipse time");
        assert!(
            matches!(
                shadow_axis_surface_point(&e, &Ellipsoid::WGS84),
                Err(EclipseError::Solver(SolverError::RootNotBracketed))
            ),
            "oblate axis must miss the surface (RootNotBracketed) even though gamma<1"
        );

        // クランプのピン: magnitude は **m=0**（クランプ後）の eclipse_magnitude(0,l1,l2) に厳密一致。
        // l1/(l1+l2)=0.54/0.53≈1.0189。クランプなしなら m=gamma−1≈−0.002 ⇒ eclipse_magnitude(−0.002,..)
        // が **より大**（中心値超え）になるため、m=0 値との一致がクランプ除去 / `.max→.min` を殺す。
        let want_mag = eclipse_magnitude(0.0, l1, l2);
        assert_eq!(
            sol.greatest.magnitude, want_mag,
            "magnitude は m=0（.max(0.0) クランプ後）由来でなければならない"
        );
        assert!(
            (sol.greatest.magnitude.0 - 0.54 / 0.53).abs() < 1e-9,
            "clamped magnitude {} must be l1/(l1+l2)=0.54/0.53≈1.0189 (m=0)",
            sol.greatest.magnitude.0
        );
        assert!(
            sol.greatest.magnitude.0 > 1.0,
            "non-central total magnitude {} must exceed 1",
            sol.greatest.magnitude.0
        );

        // obscuration≈1: m=0 ⇒ separation=0、視半径比 ρ=(l1−l2)/(l1+l2)>1（皆既）ゆえ太陽が完全内包。
        assert!(
            (sol.greatest.obscuration.0 - 1.0).abs() < 1e-9,
            "non-central total central obscuration {} must be 1.0 (m=0, contained)",
            sol.greatest.obscuration.0
        );
    }

    // ====================================================================
    // classify_global_kind（ISSUE-043 S6b-iii・全球の日食種別判定）
    // ====================================================================
    //
    // ## オラクル戦略（追認回避）
    // - 実日食（2017/2023）の種別は **NASA 公表事実**（既知）をオラクルにする。実装内部の
    //   `classify` 出力と照合しない（追認回避）。同時に独立に l2 符号を `source.at(t).l2` で
    //   標本化し、純皆既なら全域で l2<0・純金環なら全域で l2>0 を確認して幾何を裏取りする。
    // - Hybrid 機構は **合成供給源** で証明する: gamma が中心内部極小（gamma_min<0.9972）を持ち、
    //   かつ l2 が中心食区間で 0 を跨ぐ（中心 −δ と +δ で逆符号）ことを **独立に標本化** して
    //   「真にハイブリッドな幾何」であることを縛る。単一時刻の classify は中心（l2≈0）では Hybrid を
    //   返せない（時系列 [U1,U4] 走査が必須）こともコメントで明示する。

    /// classify_global_kind 用の Brent 設定（接触/最大食 solver 共通・±2s 目標の 1/10 以下）。
    fn classify_config() -> RootConfig {
        RootConfig {
            x_tolerance_days: 1e-9,
            max_iterations: 200,
        }
    }

    /// 種別1（実 Total・2017-08-21）: 探索窓 [2457986,2457988] ⇒ Some(Total)（NASA 事実）。
    /// 独立裏取り: 最大食付近で l2<0 が **終始** 保たれる（符号反転なし＝純皆既で Hybrid ではない）。
    #[test]
    fn classify_global_total_real_2017() {
        let dt = EspenakMeeusDeltaT;
        let src = DirectBesselianSource::new(G_R_SUN, G_R_MOON, &dt, solve_window_2017());

        // 独立裏取り: 最大食 TT-JD≈2457986.768 の前後で l2<0（皆既・符号反転なし）。
        let center_jd = 2_457_986.768;
        for &off in &[-0.02_f64, -0.005, 0.0, 0.005, 0.02] {
            let e = src
                .at(g_tt_jd(center_jd + off))
                .expect("source.at should succeed near 2017 eclipse");
            assert!(
                e.l2 < 0.0,
                "2017 is a pure total: l2 must stay negative, got l2={} at off={off}",
                e.l2
            );
        }

        let kind =
            classify_global_kind(&src, classify_config()).expect("2017 eclipse should classify");
        assert_eq!(kind, Some(Total), "2017-08-21 is a total eclipse (NASA)");
    }

    /// 種別2（実 Annular・2023-10-14）: 探索窓 [2460231.5,2460233.0] ⇒ Some(Annular)（NASA 事実）。
    /// 独立裏取り: 最大食付近で l2>0 が **終始** 保たれる（金環・符号反転なし）。
    #[test]
    fn classify_global_annular_real_2023() {
        let dt = EspenakMeeusDeltaT;
        let src = DirectBesselianSource::new(G_R_SUN, G_R_MOON, &dt, solve_window_2023());

        // 独立裏取り: 最大食 TT-JD≈2460232.25 の前後で l2>0（金環・符号反転なし）。
        let center_jd = 2_460_232.25;
        for &off in &[-0.02_f64, -0.005, 0.0, 0.005, 0.02] {
            let e = src
                .at(g_tt_jd(center_jd + off))
                .expect("source.at should succeed near 2023 eclipse");
            assert!(
                e.l2 > 0.0,
                "2023 is a pure annular: l2 must stay positive, got l2={} at off={off}",
                e.l2
            );
        }

        let kind =
            classify_global_kind(&src, classify_config()).expect("2023 eclipse should classify");
        assert_eq!(
            kind,
            Some(Annular),
            "2023-10-14 is an annular eclipse (NASA)"
        );
    }

    /// Hybrid 機構を証明する時変合成供給源（global_contacts の `SyntheticGammaSource` と同形だが
    /// **l2 を時間の関数**にする）。
    ///   x(jd) = X_MIN + K·(jd−center)²（X_MIN=0.3<0.9972 ⇒ 中心食, K=50, y=0）。
    ///   l2(jd) = L2_SLOPE·(jd−center)（center で 0 を跨ぐ ⇒ 中心食区間で金環⇄皆既が切替わる）。
    /// gamma=√(x²+0)=x が center で内部極小。中心食区間 [U1,U4] で l2 が符号反転 ⇒ Hybrid。
    struct HybridSource {
        center_jd: f64,
        window: TimeInterval<TtInstant>,
    }

    /// gamma の中心極小値（=X_MIN）。0.9972 未満ゆえ中心食（Total/Annular）。
    const HYBRID_X_MIN: f64 = 0.3;
    /// 放物線曲率（小さく正・中心内部極小を作りつつ窓端で gamma を半影限界超へ持ち上げる）。
    const HYBRID_K: f64 = 50.0;
    /// l2 の時間勾配。l2=L2_SLOPE·(jd−center)。±0.02 day で |l2|≈0.01 に達する。
    const HYBRID_L2_SLOPE: f64 = 0.5;

    impl BesselianSource for HybridSource {
        /// x=X_MIN+K·(jd−center)²（中心極小 0.3）, y=0, l2=L2_SLOPE·(jd−center)（center で符号反転）。
        fn at(&self, t: TtInstant) -> Result<InstantaneousBesselianElements, EclipseError> {
            let dj = t.jd2().jd() - self.center_jd;
            Ok(InstantaneousBesselianElements {
                x: HYBRID_X_MIN + HYBRID_K * dj * dj,
                y: 0.0,
                declination: Radians(0.0),
                mu: Radians(0.0),
                l1: 0.54,
                l2: HYBRID_L2_SLOPE * dj,
                tan_f1: 0.0047,
                tan_f2: 0.0046,
                time_tt: t,
            })
        }

        fn fit_interval(&self) -> TimeInterval<TtInstant> {
            self.window
        }
    }

    /// 種別3（合成 Hybrid・厳密ピン）: gamma が中心内部極小（gamma_min=0.3<0.9972 ⇒ 中心食）を持ち、
    /// l2 が中心食区間で 0 を跨ぐ供給源 ⇒ `Some(Hybrid)`。
    ///
    /// 独立裏取り（呼出前）: (a) gamma_min=X_MIN=0.3<0.9972（中心食であること）、
    /// (b) center−δ と center+δ で l2 が逆符号（中心食区間で金環⇄皆既が切替わる＝真のハイブリッド幾何）。
    /// 注: center 単一時刻では l2≈0 で classify は Hybrid を返せない。時系列 [U1,U4] の符号走査が必須。
    #[test]
    fn classify_global_hybrid_synthetic() {
        let center_jd = 2_457_986.768;
        // 窓: center ± 0.2 day。端で x=0.3+50·0.04=2.3>半影限界 ⇒ P1/P4・U1/U4 が窓内に括られる。
        let half_day = 0.2;
        let src = HybridSource {
            center_jd,
            window: TimeInterval {
                start: g_tt_jd(center_jd - half_day),
                end: g_tt_jd(center_jd + half_day),
            },
        };

        // (a) 中心食であること: gamma_min=X_MIN=0.3 < 0.9972（CENTRAL_LIMIT）。
        let g_center = gamma_at(&src, g_tt_jd(center_jd));
        assert!(
            (g_center - HYBRID_X_MIN).abs() < 1e-12,
            "center gamma {g_center} should equal X_MIN={HYBRID_X_MIN}"
        );
        assert!(
            g_center < 0.9972,
            "center gamma {g_center} must be < 0.9972 (central ⇒ Total/Annular base kind)"
        );

        // (b) 中心食区間で l2 が両符号を取る（真のハイブリッド幾何）。
        let delta = 0.02; // 中心食区間 |dj|<0.118 day の内側。
        let l2_minus = src
            .at(g_tt_jd(center_jd - delta))
            .expect("source.at should succeed")
            .l2;
        let l2_plus = src
            .at(g_tt_jd(center_jd + delta))
            .expect("source.at should succeed")
            .l2;
        assert!(
            l2_minus < 0.0,
            "l2 at center−δ must be negative (total side), got {l2_minus}"
        );
        assert!(
            l2_plus > 0.0,
            "l2 at center+δ must be positive (annular side), got {l2_plus}"
        );
        // 中心（l2≈0）の単一時刻 classify は Hybrid を返せない。時系列 [U1,U4] 走査が必須。
        let e_center = src.at(g_tt_jd(center_jd)).expect("source.at");
        assert!(
            e_center.l2.abs() < 1e-9,
            "at center l2≈0 (single-instant classify cannot yield Hybrid), got {}",
            e_center.l2
        );

        let kind = classify_global_kind(&src, classify_config())
            .expect("hybrid synthetic source should classify");
        assert_eq!(
            kind,
            Some(SolarEclipseKind::Hybrid),
            "central + l2 sign change across [U1,U4] must be Hybrid"
        );
    }

    /// Hybrid 機構を **区間内部の符号反転** で証明する時変合成供給源（`HybridSource` と同形だが
    /// l2 を時間の**放物線**にする）。`classify_global_hybrid_synthetic` の線形 l2 は両端で逆符号
    /// なので端点 2 点だけでも検出でき、走査密度（SAMPLES）が load-bearing にならない。本源は
    ///   x(jd) = X_MIN + K·(jd−center)²（X_MIN=0.3<0.9972 ⇒ 中心食, K=50, y=0）。
    ///   l2(jd) = L2_A·((jd−center)²/HALF²) − L2_B（中心で −L2_B<0＝皆既, 端で正＝金環）。
    /// HALF を本影接触距離に取る（gamma=1+|l2|=1.004 となる dj）。これにより [U1,U4] の **両端**で
    /// l2>0（金環）・**中央**で l2<0（皆既）の「金環-皆既-金環」帯になり、端点だけ標本化すると両端
    /// 正で Annular に誤判定する。Hybrid 検出には区間内部の標本（SAMPLES）が必須＝SAMPLES が
    /// load-bearing になる（端点削減ミュータントを殺す）。
    struct HybridInteriorSource {
        center_jd: f64,
        window: TimeInterval<TtInstant>,
    }

    /// l2 放物線の中心値（−L2_B）。中心で l2<0（皆既）。
    const HYBRID_I_L2_B: f64 = 0.008;
    /// l2 放物線の振幅係数。l2(HALF)=L2_A−L2_B=+0.004（端で金環）。
    const HYBRID_I_L2_A: f64 = 0.012;
    /// 本影接触距離 HALF[day]: gamma(HALF)=X_MIN+K·HALF²=1+|l2(HALF)|=1.004 を満たす自己無撞着な値。
    /// HALF² = (1.004−X_MIN)/K。l2(HALF)=+0.004 ⇒ |l2|=0.004 ⇒ gamma=1.004 が閉じる。
    fn hybrid_i_half() -> f64 {
        ((1.0 + HYBRID_I_L2_A - HYBRID_I_L2_B - HYBRID_X_MIN) / HYBRID_K).sqrt()
    }

    /// HybridInteriorSource の l2(dj)（dj=jd−center）。放物線 L2_A·(dj²/HALF²)−L2_B。
    fn hybrid_i_l2(dj: f64) -> f64 {
        let half = hybrid_i_half();
        HYBRID_I_L2_A * (dj * dj / (half * half)) - HYBRID_I_L2_B
    }

    impl BesselianSource for HybridInteriorSource {
        /// x=X_MIN+K·(jd−center)²（中心極小 0.3）, y=0, l2 は中心負・端正の放物線（区間内部で符号反転）。
        fn at(&self, t: TtInstant) -> Result<InstantaneousBesselianElements, EclipseError> {
            let dj = t.jd2().jd() - self.center_jd;
            Ok(InstantaneousBesselianElements {
                x: HYBRID_X_MIN + HYBRID_K * dj * dj,
                y: 0.0,
                declination: Radians(0.0),
                mu: Radians(0.0),
                l1: 0.54,
                l2: hybrid_i_l2(dj),
                tan_f1: 0.0047,
                tan_f2: 0.0046,
                time_tt: t,
            })
        }

        fn fit_interval(&self) -> TimeInterval<TtInstant> {
            self.window
        }
    }

    /// 種別3b（合成 Hybrid・**区間内部**符号反転）: l2 が中央でのみ負・[U1,U4] の両端で正となる
    /// 放物線供給源 ⇒ `Some(Hybrid)`。`classify_global_hybrid_synthetic`（線形 l2・端点で逆符号）と
    /// 違い、**端点だけ標本化すると両端 l2>0 で Annular に誤判定**する。Hybrid を出すには区間内部の
    /// 標本（`l2_changes_sign` の SAMPLES）が必須 ⇒ 走査密度が load-bearing（SAMPLES 削減ミュータント
    /// を殺す）。
    ///
    /// 独立裏取り（呼出前・`src.at` 由来）: (a) gamma_min=X_MIN=0.3<0.9972（中心食 ⇒ base=Total/Annular）、
    /// (b) l2(center)<0（中央は皆既）、(c) l2 が **本影接触距離 ±HALF**（U1/U4 が立つ dj）で >0（端は金環）。
    /// すなわち符号反転は **区間内部**で起き、端点は同符号（両方正）。端点のみ標本化する実装はこれを
    /// 取りこぼし Annular に化けるため、SAMPLES の密度がここで load-bearing になる。
    #[test]
    fn classify_global_hybrid_interior_crossing() {
        let center_jd = 2_457_986.768;
        // 窓: center ± 0.2 day。端で x=0.3+50·0.04=2.3>1+l1=1.54 ⇒ P1/P4・U1/U4 が窓内に括られる。
        let half_day = 0.2;
        let src = HybridInteriorSource {
            center_jd,
            window: TimeInterval {
                start: g_tt_jd(center_jd - half_day),
                end: g_tt_jd(center_jd + half_day),
            },
        };

        // (a) 中心食であること: gamma_min=X_MIN=0.3 < 0.9972（CENTRAL_LIMIT）⇒ base=Total/Annular。
        let g_center = gamma_at(&src, g_tt_jd(center_jd));
        assert!(
            (g_center - HYBRID_X_MIN).abs() < 1e-12,
            "center gamma {g_center} should equal X_MIN={HYBRID_X_MIN}"
        );
        assert!(
            g_center < 0.9972,
            "center gamma {g_center} must be < 0.9972 (central ⇒ Total/Annular base kind)"
        );

        // (b) 中央は皆既: l2(center) < 0。本影接触距離 HALF。
        let half = hybrid_i_half();
        let l2_center = src.at(g_tt_jd(center_jd)).expect("source.at").l2;
        // (c) 本影接触距離 ±HALF（U1/U4 が立つ dj）では l2 > 0（金環）。符号反転は区間内部。
        let l2_u1 = src.at(g_tt_jd(center_jd - half)).expect("source.at").l2;
        let l2_u4 = src.at(g_tt_jd(center_jd + half)).expect("source.at").l2;
        assert!(
            l2_center < 0.0,
            "l2 at center must be negative (total in the middle), got {l2_center}"
        );
        assert!(
            l2_u1 > 0.0 && l2_u4 > 0.0,
            "l2 at the umbral-contact distance ±HALF must be positive (annular at the ends), \
             got u1={l2_u1}, u4={l2_u4}"
        );
        // 端点が同符号（両方正）ゆえ、[U1,U4] を **端点だけ** 標本化する実装は符号反転を取りこぼし
        // Annular に誤判定する。区間内部の標本（l2_changes_sign の SAMPLES）が必須＝SAMPLES が
        // load-bearing。下の確認で l2(center)<0 と l2(±HALF)>0 の符号反転が区間内部にあることを縛る。
        // HALF が本影接触距離（gamma=1+|l2|）であることの独立裏取り: gamma(HALF)≈1+|l2(HALF)|。
        let gamma_half = {
            let e = src.at(g_tt_jd(center_jd + half)).expect("source.at");
            (e.x * e.x + e.y * e.y).sqrt()
        };
        assert!(
            (gamma_half - (1.0 + l2_u4.abs())).abs() < 1e-7,
            "HALF must be the umbral-contact distance: gamma(HALF)={gamma_half} ≈ 1+|l2|={}",
            1.0 + l2_u4.abs()
        );

        let kind = classify_global_kind(&src, classify_config())
            .expect("interior-crossing hybrid synthetic source should classify");
        assert_eq!(
            kind,
            Some(SolarEclipseKind::Hybrid),
            "central with INTERIOR l2 sign change (endpoints both annular) must be Hybrid; \
             endpoint-only sampling would misclassify as Annular (SAMPLES is load-bearing)"
        );
    }

    /// 純皆既（pure total）を証明する時変合成供給源（`HybridSource` と同形だが l2 を**緩い線形**に
    /// し、[U1,U4] 内では終始 l2<0、零交差を [U1,U4] の **外**・1 日以内に置く）。
    ///   x(jd) = X_MIN + K·(jd−center)²（X_MIN=0.3<0.9972 ⇒ 中心食, K=50, y=0）。
    ///   l2(jd) = L2_C0 + L2_SLOPE·(jd−center)（中心 −0.008<0＝皆既）。
    /// [U1,U4]≈center±0.119 day（gamma=1+|l2|≈1.004）で l2∈[−0.0116,−0.0044]＝全域負（純皆既）。
    /// 零交差は center+L2_C0/|L2_SLOPE|=center+0.267 day で [U1,U4] の外。だが `l2_changes_sign` の
    /// 走査区間式 `(t1−t0)` を `(t1/t0)`（≈1.0 day）に取り違えるミュータントは ≈[U1,U1+1day] を走査し、
    /// そこ（dj≈0.119+1.0）で l2≈+0.026>0 に達して符号反転を検出 ⇒ Hybrid に誤判定する。
    struct PureTotalLinearSource {
        center_jd: f64,
        window: TimeInterval<TtInstant>,
    }

    /// l2 の中心値（−0.008<0＝皆既）。
    const PURE_TOTAL_L2_C0: f64 = -0.008;
    /// l2 の時間勾配。零交差は center+0.008/0.03≈center+0.267 day（[U1,U4] の外・1 日以内）。
    const PURE_TOTAL_L2_SLOPE: f64 = 0.03;

    /// PureTotalLinearSource の l2(dj)（dj=jd−center）。線形 L2_C0+L2_SLOPE·dj。
    fn pure_total_l2(dj: f64) -> f64 {
        PURE_TOTAL_L2_C0 + PURE_TOTAL_L2_SLOPE * dj
    }

    impl BesselianSource for PureTotalLinearSource {
        /// x=X_MIN+K·(jd−center)²（中心極小 0.3）, y=0, l2=L2_C0+L2_SLOPE·dj（中心負・緩い線形）。
        fn at(&self, t: TtInstant) -> Result<InstantaneousBesselianElements, EclipseError> {
            let dj = t.jd2().jd() - self.center_jd;
            Ok(InstantaneousBesselianElements {
                x: HYBRID_X_MIN + HYBRID_K * dj * dj,
                y: 0.0,
                declination: Radians(0.0),
                mu: Radians(0.0),
                l1: 0.54,
                l2: pure_total_l2(dj),
                tan_f1: 0.0047,
                tan_f2: 0.0046,
                time_tt: t,
            })
        }

        fn fit_interval(&self) -> TimeInterval<TtInstant> {
            self.window
        }
    }

    /// 種別3c（合成 純皆既・走査区間ピン）: [U1,U4] 内で l2<0 が終始保たれる純皆既供給源 ⇒
    /// `Some(Total)`（Hybrid に化けない）。`l2_changes_sign` の走査区間式 `t0+(t1−t0)·frac` を
    /// `t0+(t1/t0)·frac`（t1/t0≈1.0 ⇒ 走査が ≈[U1,U1+1day]）に取り違えるミュータントを殺す。
    ///
    /// なぜ殺せるか: 正しい [U1,U4]≈center±0.119 day では l2 が全域負ゆえ符号反転なし ⇒ Total。だが
    /// 零交差を center+0.267 day（[U1,U4] の外・1 日以内）に置いてあるので、ミュータントの広い
    /// ≈[U1,U1+1day] 区間は l2>0 の領域まで走査して両符号を検出 ⇒ Hybrid に誤判定する。よって
    /// `Some(Total)` を assert する本テストはミュータントを撃破する（純皆既が Hybrid に化けない縛り）。
    ///
    /// 独立裏取り（呼出前・`src.at` 由来）: (a) gamma_min=X_MIN=0.3<0.9972（中心食 ⇒ base=Total）、
    /// (b) l2(center)<0、(c) [U1,U4]≈center±0.119 day の **両端**で l2<0（純皆既が中心食区間で完結）。
    #[test]
    fn classify_global_total_not_hybrid_when_l2_negative_in_central_interval() {
        let center_jd = 2_457_986.768;
        // 窓: center ± 0.2 day。端で x=0.3+50·0.04=2.3>1.5433+l2 ⇒ P1/P4・U1/U4 が窓内に括られる。
        // かつ center+0.2 の l2=−0.008+0.03·0.2=−0.002<0 ゆえ fit 窓内には零交差(+0.267)を含めない
        // （greatest=center で l2<0 ⇒ base=Total、真の [U1,U4] も全域負を保つ）。
        let half_day = 0.2;
        let src = PureTotalLinearSource {
            center_jd,
            window: TimeInterval {
                start: g_tt_jd(center_jd - half_day),
                end: g_tt_jd(center_jd + half_day),
            },
        };

        // (a) 中心食であること: gamma_min=X_MIN=0.3 < 0.9972（CENTRAL_LIMIT）⇒ base=Total。
        let g_center = gamma_at(&src, g_tt_jd(center_jd));
        assert!(
            (g_center - HYBRID_X_MIN).abs() < 1e-12,
            "center gamma {g_center} should equal X_MIN={HYBRID_X_MIN}"
        );
        assert!(
            g_center < 0.9972,
            "center gamma {g_center} must be < 0.9972 (central ⇒ Total base kind)"
        );

        // (b)(c) 純皆既: l2(center)<0 かつ [U1,U4]≈center±0.119 day の両端でも l2<0（符号反転なし）。
        // U1/U4 は本影接触距離（gamma=1+|l2|≈1.004）で dj≈±0.119 day。
        let u_dj = 0.119;
        let l2_center = src.at(g_tt_jd(center_jd)).expect("source.at").l2;
        let l2_u1 = src.at(g_tt_jd(center_jd - u_dj)).expect("source.at").l2;
        let l2_u4 = src.at(g_tt_jd(center_jd + u_dj)).expect("source.at").l2;
        assert!(
            l2_center < 0.0,
            "l2 at center must be negative (total), got {l2_center}"
        );
        assert!(
            l2_u1 < 0.0 && l2_u4 < 0.0,
            "l2 at the umbral-contact distance ±0.119 day (≈U1/U4) must BOTH be negative \
             (pure total across [U1,U4]), got u1={l2_u1}, u4={l2_u4}"
        );
        // 純皆既ゆえ正しい [U1,U4] には符号反転がなく Total。だがミュータント `(t1−t0)→(t1/t0)` は
        // ≈[U1,U1+1day] を走査し、そこ（dj≈0.119+1.0）で l2 が正に転じて両符号を検出 ⇒ Hybrid に化ける。
        // よって正しい [U1,U4] 区間でのみ Total となり、走査区間をピン留めする。零交差が [U1,U4] の外・
        // 1 日以内にあることを独立に確認する。
        let l2_cross = src
            .at(g_tt_jd(
                center_jd + PURE_TOTAL_L2_C0.abs() / PURE_TOTAL_L2_SLOPE,
            ))
            .expect("source.at")
            .l2;
        assert!(
            l2_cross.abs() < 1e-9,
            "l2 zero crossing must sit at center+|C0|/SLOPE≈+0.267 day (outside [U1,U4]), got {l2_cross}"
        );
        let l2_mutant = src
            .at(g_tt_jd(center_jd + u_dj + 1.0))
            .expect("source.at")
            .l2;
        assert!(
            l2_mutant > 0.0,
            "within the mutant's wider ≈[U1,U1+1day] scan, l2 turns positive (would force Hybrid), \
             got {l2_mutant}"
        );

        let kind = classify_global_kind(&src, classify_config())
            .expect("pure total synthetic source should classify");
        assert_eq!(
            kind,
            Some(Total),
            "pure total (l2<0 across [U1,U4]) must be Total, not Hybrid; the scan interval must be \
             [U1,U4], not ≈[U1,U1+1day] (kills `(t1-t0)→(t1/t0)` mutant)"
        );
    }

    /// 種別4（合成 Partial）: gamma_min を中心限界と半影限界の **間** に置く合成供給源
    /// （X_MIN=1.2 ⇒ 0.9972<1.2<1.5433+l2）。l2 固定（−0.01）。⇒ `Some(Partial)`。
    /// 独立裏取り: gamma_min=1.2 は部分食バンド（中心食でない ⇒ U1/U4 は None・Hybrid 上書きなし）。
    #[test]
    fn classify_global_partial_synthetic() {
        let center_jd = 2_457_986.768;
        // global_contacts の非中心テストと同形: x_min=1.2, k=50, l1=0.54, l2=−0.01, 半幅 0.15。
        let half_day = 0.15;
        let src = SyntheticGammaSource {
            center_jd,
            x_min: 1.2,
            k: 50.0,
            l1: 0.54,
            l2: -0.01,
            window: TimeInterval {
                start: g_tt_jd(center_jd - half_day),
                end: g_tt_jd(center_jd + half_day),
            },
        };

        // 独立裏取り: gamma_min=1.2 が部分食バンド（0.9972 ≤ g < 1.5433+l2）にある。上限は符号付き
        // l2（−0.01）込みで 1.5433+(−0.01)=1.5323（PENUMBRA_LIMIT+l2, |l2| や +0.01 ではない）。
        let g_center = gamma_at(&src, g_tt_jd(center_jd));
        assert!(
            (0.9972..(1.5433 + (-0.01))).contains(&g_center),
            "gamma_min {g_center} must be in partial band (0.9972,1.5433+l2=1.5323)"
        );

        let kind = classify_global_kind(&src, classify_config())
            .expect("partial synthetic source should classify");
        assert_eq!(kind, Some(Partial), "gamma in partial band ⇒ Partial");
    }

    /// 種別5（合成 No eclipse）: gamma が全域で半影限界を超える供給源（X_MIN=2.0>1.5433）⇒ `Ok(None)`。
    /// 独立裏取り: gamma_min=2.0 は半影限界外（日食なし）。
    #[test]
    fn classify_global_no_eclipse_synthetic() {
        let center_jd = 2_457_986.768;
        let half_day = 0.05;
        let src = SyntheticGammaSource {
            center_jd,
            x_min: 2.0,
            k: 50.0,
            l1: 0.54,
            l2: -0.01,
            window: TimeInterval {
                start: g_tt_jd(center_jd - half_day),
                end: g_tt_jd(center_jd + half_day),
            },
        };

        // 独立裏取り: gamma_min=2.0 > 1.5433（半影限界）⇒ 日食なし。
        let g_center = gamma_at(&src, g_tt_jd(center_jd));
        assert!(
            g_center > 1.5433,
            "gamma_min {g_center} must exceed penumbra limit (no eclipse)"
        );

        let kind = classify_global_kind(&src, classify_config())
            .expect("no-eclipse window must be Ok(None), not Err");
        assert_eq!(kind, None, "gamma beyond penumbra limit ⇒ Ok(None)");
    }

    // 種別6（任意・実 Hybrid 2023-04-20）は **意図的に省略**する。2023-04-20 は実ハイブリッド食だが、
    // 解析暦が中心線上の l2 符号反転を再現するかは k/ΔT 慣習差に敏感で borderline になりうる。
    // flaky な実値を assert すると偽陽性/偽陰性を生むため、Hybrid の堅牢なピンは合成テスト
    // `classify_global_hybrid_synthetic`（種別3）に委ねる（spec の方針どおり）。

    /// 種別5 用の時変合成供給源（global_contacts の `SyntheticGammaSource` と同形・本モジュール内に複製）。
    /// gamma=|x| が center で内部極小を持つ。x=x_min+k·(jd−center)², y=0, l1/l2 固定。
    struct SyntheticGammaSource {
        center_jd: f64,
        x_min: f64,
        k: f64,
        l1: f64,
        l2: f64,
        window: TimeInterval<TtInstant>,
    }

    impl BesselianSource for SyntheticGammaSource {
        fn at(&self, t: TtInstant) -> Result<InstantaneousBesselianElements, EclipseError> {
            let dj = t.jd2().jd() - self.center_jd;
            Ok(InstantaneousBesselianElements {
                x: self.x_min + self.k * dj * dj,
                y: 0.0,
                declination: Radians(0.0),
                mu: Radians(0.0),
                l1: self.l1,
                l2: self.l2,
                tan_f1: 0.0047,
                tan_f2: 0.0046,
                time_tt: t,
            })
        }

        fn fit_interval(&self) -> TimeInterval<TtInstant> {
            self.window
        }
    }
}
