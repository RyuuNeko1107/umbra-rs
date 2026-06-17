//! 影軸の地表貫通点（逆射影, `docs/algorithms/08-global.md` 式 8.10、ISSUE-023 S6a-i,
//! Explanatory Supplement to the Astronomical Almanac §11 / Meeus *Astronomical Algorithms* Ch.54）。
//!
//! ある TT の瞬時ベッセル要素（影軸交点 (x,y)・赤緯 d・時角 μ）から、月影軸が WGS84 地表を
//! 貫く**太陽側（ζ>0）の点**の測地座標 (φ, λ_east) を求める。これは観測者→基本面の前方射影
//! [`project_observer_to_fundamental`](crate::projection) の**逆**であり、検証は前方射影の
//! 往復一致（前方射影し直すと (ξ,η)=(x,y)・ζ>0）で行う（ISSUE-043 S6a, ユーザ確定）。
//!
//! 中心食（影軸が地表に当たる）でのみ貫通点が存在する。影軸が地球を外れる（部分/非中心で
//! 軸が地表に届かない）場合は `Err(Solver(RootNotBracketed))` を返す（最大食地点の部分/非中心
//! 定義は S6b の責務）。

// shadow_axis_surface_point（pub(crate)）は ISSUE-043 S6a-ii（最大食組立）が消費するまで未使用。
#![allow(dead_code)]
// ζ 粗走査の添字 ⇔ f64 変換は小さな整数のみ（天文量ではない, local_maximum.rs と同方針）。
#![allow(clippy::cast_precision_loss)]

use umbra_core::ellipsoid::Ellipsoid;
use umbra_core::{brent_root, EastLongitude, GeodeticLatitude, Radians, SolverError};

use crate::besselian::InstantaneousBesselianElements;
use crate::error::EclipseError;

use umbra_geo::GeoPoint;

/// 影軸方向（ζ）の粗走査上限（Re）。地表点の ζ は地心動径 ρ ≤ 1 を超えないので 1.05 で十分な余裕。
const ZETA_SCAN_MAX: f64 = 1.05;
/// 粗走査の分割数（符号反転で太陽側根をブラケットするための解像度。掠め交点でも逃さない細かさ）。
const ZETA_SCAN_STEPS: usize = 256;
/// ζ 根の収束許容（Re）。中心線位置 sub-km（≈7.8e-5 Re, accuracy.md §2.1）を十分下回る。
const ZETA_ROOT_TOL: f64 = 1e-12;
/// Brent 反復上限。
const ZETA_ROOT_MAX_ITER: usize = 100;

/// 影軸の地表貫通点（太陽側 ζ>0）を測地座標で返す。中心食でなければ `Err(Solver(RootNotBracketed))`。
///
/// 検証済み前方射影 [`project_observer_to_fundamental`](crate::projection) の逆。前方射影の回転
/// `R_x(d)·R_z` の**転置**で観測者地心成分を影軸成分 ζ の関数として組み、海面（h=0）の子午線楕円
/// 拘束（`ellipsoid::observer_geocentric` 由来の不変量 `ρcos²+(ρsin/(1−f))²=1`）の零点を Brent で
/// 解く。零点の (ρsin, ρcos, 局地時角 H) から測地緯度 φ・東経 λ=H−μ を復元する。
///
/// 太陽側（ζ>0）の地表交点＝大きい方の ζ 根（+ζ が太陽向き, projection.rs）。軸が地表に届かない
/// （部分/非中心で軸が地球を外す）と零点が無く `Err(Solver(RootNotBracketed))`。
pub(crate) fn shadow_axis_surface_point(
    elements: &InstantaneousBesselianElements,
    ellipsoid: &Ellipsoid,
) -> Result<GeoPoint, EclipseError> {
    let x = elements.x;
    let y = elements.y;
    let (sin_d, cos_d) = elements.declination.0.sin_cos();
    let omf = 1.0 - ellipsoid.f; // b/a
    let inv_omf2 = 1.0 / (omf * omf); // 1/(1−f)²

    // 逆回転（前方 R_x(d)·R_z の転置）で観測者地心成分を ζ の関数に:
    //   px = ζ·cos d − y·sin d,  py = x,  pz = ζ·sin d + y·cos d
    //   ρcos² = px² + py²,  ρsin = pz
    // h=0 の子午線楕円拘束（ellipsoid.rs 不変量）: ρcos² + (ρsin/(1−f))² = 1。
    // 残差 r(ζ) = ρcos² + (ρsin/(1−f))² − 1。零点が軸の地表交点。
    let residual = |zeta: f64| -> f64 {
        let px = zeta * cos_d - y * sin_d;
        let pz = zeta * sin_d + y * cos_d;
        let rho_cos2 = px * px + x * x;
        rho_cos2 + pz * pz * inv_omf2 - 1.0
    };

    // 太陽側（最大 ζ>0）根のブラケットを降順粗走査で取り、Brent で精解する。
    // r は ζ² の係数 cos²d + sin²d/(1−f)² > 0 ゆえ下に凸（開口上向き）で、2 根の間で負・外で正
    // なので、上端から下る最初の符号反転が大きい方の根 ζ₊（＝太陽側）。軸が地表を外す（実根なし）と
    // 符号反転が無くブラケット不成立 → `RootNotBracketed`。
    let (a, b) = descending_sign_change_bracket(&residual, ZETA_SCAN_MAX, ZETA_SCAN_STEPS)
        .ok_or(EclipseError::Solver(SolverError::RootNotBracketed))?;
    let zeta = brent_root(residual, a, b, ZETA_ROOT_TOL, ZETA_ROOT_MAX_ITER)?;

    // 零点 ζ から地心成分・局地時角を復元。
    let px = zeta * cos_d - y * sin_d;
    let pz = zeta * sin_d + y * cos_d;
    let rho_cos = (px * px + x * x).sqrt();
    let rho_sin = pz;
    // 局地時角 H = atan2(py, px) = atan2(x, px)（前方の H=μ+λ_east と整合）。
    let h = x.atan2(px);
    // 東経 λ_east = H − μ（[-π,π) へ正規化）。
    let longitude = EastLongitude::from_radians(Radians::new(h - elements.mu.0));
    // 測地緯度: tan φ = ρsin / (ρcos·(1−f)²)（observer_geocentric(h=0) の逆）。
    // ρcos ≥ 0・(1−f)²>0 なので atan2 は [-π/2, π/2] を返し常に有効。
    let phi = rho_sin.atan2(rho_cos * omf * omf);
    let latitude = GeodeticLatitude::from_radians(Radians::new(phi))?;

    Ok(GeoPoint::new(latitude, longitude))
}

/// `[0, top]` を `top` から降順に粗走査し、`f` の符号反転を含む最初のブラケット `[lo, hi]` を返す。
///
/// **粗走査機構**（`conjunction::solve_zero_in_window` / `local_maximum::scan_point_count` と同カテゴリ）:
/// 区間内に根が一意なら、刻み・点配置によらず Brent が同じ真根へ収束するため、本関数内の算術
/// （`top/steps` の刻み・`top − step·i` の格子・`r_lo·r_hi` の符号判定）は結果に影響しない解像度要素。
/// `mutation.yml` で `--exclude-re 'in descending_sign_change_bracket'` 除外する
/// （`docs/reviews/mutation-axis-intercept.md`）。根が無ければ `None`。
fn descending_sign_change_bracket<F: Fn(f64) -> f64>(
    f: &F,
    top: f64,
    steps: usize,
) -> Option<(f64, f64)> {
    let step = top / steps as f64;
    let mut hi = top;
    let mut r_hi = f(hi);
    for i in 1..=steps {
        let lo = top - step * i as f64;
        let r_lo = f(lo);
        if r_lo * r_hi <= 0.0 {
            return Some((lo, hi));
        }
        hi = lo;
        r_hi = r_lo;
    }
    None
}

#[cfg(test)]
mod tests {
    //! ISSUE-023 S6a-i 受け入れテスト（strict・逆射影 = 前方射影の逆）。
    //!
    //! オラクル戦略（追認回避・ユーザ確定 ISSUE-043 S6a）:
    //! 本関数は**検証済み**前方射影 [`project_observer_to_fundamental`]（ISSUE-024）の逆。
    //! 逆の内部実装を再実装して突き合わせる（＝追認）ことはせず、戻り値 (φ,λ) を
    //! **前方射影し直して** (ξ,η)=(x,y)・ζ>0 の往復一致で縛る。
    //! `project_observer_to_fundamental` と `observer_geocentric`（ISSUE-010/011）は
    //! 独立に検証済みなので、これは独立オラクルになる。
    //! 具体的な数値（case 3 の φ=0・λ=−μ、case 4 の南北符号）は前方射影の式から手で導出する。
    #![allow(clippy::excessive_precision)]

    use super::*;
    use core::f64::consts::PI;
    use umbra_core::constants::{EARTH_EQUATORIAL_RADIUS_M, SOLAR_RADIUS_KM};
    use umbra_core::ellipsoid::{observer_geocentric, Ellipsoid};
    use umbra_core::{EspenakMeeusDeltaT, JulianDate2, Radians, TimeInterval, TtInstant};

    use crate::projection::project_observer_to_fundamental;
    use crate::source::{BesselianSource, DirectBesselianSource};

    const WGS84: Ellipsoid = Ellipsoid::WGS84;
    const R_SUN: f64 = SOLAR_RADIUS_KM;
    const R_MOON: f64 = 0.2725076 * (EARTH_EQUATORIAL_RADIUS_M / 1000.0);

    /// 往復一致の許容（前方射影し直したときの (ξ,η) 残差・Re）。
    /// 精度目標は中心線位置 sub-km（≲0.5 km ≈ 7.8e-5 Re, accuracy.md §2.1 / ISSUE-023）だが、
    /// 逆解は前方射影と数値的に厳密に閉じうるので往復は強い許容で締める。
    const TOL_ROUNDTRIP: f64 = 1e-9;

    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    /// テスト用の瞬時ベッセル要素を最小構成する（逆射影で使うのは x,y,d,μ）。
    /// 他フィールド（l1,l2,tan_f,time_tt）は貫通点計算に無関係なのでダミー（有限値）。
    fn elems(x: f64, y: f64, d: f64, mu: f64) -> InstantaneousBesselianElements {
        InstantaneousBesselianElements {
            x,
            y,
            declination: Radians(d),
            mu: Radians(mu),
            l1: 0.5,
            l2: -0.01,
            tan_f1: 0.0047,
            tan_f2: 0.0046,
            time_tt: TtInstant::from_jd2(JulianDate2::new(2_451_545.0, 0.0)),
        }
    }

    /// 返点 (φ,λ) を**検証済み前方射影**へ通し、(ξ,η)=(x,y)・ζ>0 を確認する独立オラクル。
    /// 逆の内部式は一切再実装しない（追認回避）。φ・λ は `GeoPoint` の radian アクセサ由来。
    fn assert_forward_roundtrip(p: &GeoPoint, e: &InstantaneousBesselianElements, tol: f64) -> f64 {
        let phi = p.lat.radians().0;
        let lam = p.lon.radians().0;
        let obs = observer_geocentric(&WGS84, phi, 0.0);
        let r = project_observer_to_fundamental(&obs, Radians::new(lam), e);
        assert!(
            close(r.xi, e.x, tol),
            "ξ={} expected x={} (φ={phi}, λ={lam})",
            r.xi,
            e.x
        );
        assert!(
            close(r.eta, e.y, tol),
            "η={} expected y={} (φ={phi}, λ={lam})",
            r.eta,
            e.y
        );
        // 太陽側選択（ζ>0）。逆の「太陽側 vs 反太陽側」分岐を縛る（要件6）。
        assert!(r.zeta > 0.0, "ζ={} must be sunward (>0)", r.zeta);
        r.zeta
    }

    // ---- 要件1: 前方往復（主検証） ----

    /// 軸が地球に当たる複数の合成要素（gamma=√(x²+y²)<1）で、返点を前方射影し直すと
    /// (ξ,η)=(x,y)・ζ>0 に往復一致する。
    /// 採用 config（gamma）:
    ///   (0.0, 0.0) gamma=0.0（影軸直下・直撃）
    ///   (0.3,−0.2) gamma=√0.13≈0.36056
    ///   (0.1, 0.5) gamma=√0.26≈0.50990
    /// d, μ は現実的な値（d≈太陽見かけ赤緯 ±0.2 rad, μ 任意）。
    #[test]
    fn forward_roundtrip_for_axis_hitting_earth() {
        let cases = [
            (0.0_f64, 0.0_f64, 0.20_f64, 1.0_f64),
            (0.3, -0.2, 0.20, 2.7),
            (0.1, 0.5, -0.15, -1.3),
        ];
        for (x, y, d, mu) in cases {
            let e = elems(x, y, d, mu);
            // gamma<1 を独立に確認（軸が地球に当たる前提）。
            let gamma = (x * x + y * y).sqrt();
            assert!(gamma < 1.0, "config gamma={gamma} must be <1");
            let p = shadow_axis_surface_point(&e, &WGS84)
                .expect("axis hits Earth (gamma<1) ⇒ surface point exists");
            assert_forward_roundtrip(&p, &e, TOL_ROUNDTRIP);
        }
    }

    // ---- 要件2: 実 2017-08-21（往復 + 大局的妥当性のみ） ----

    /// 2017-08-21 最大食付近の実瞬時要素（`DirectBesselianSource`）で逆射影し、
    /// (a) 厳密な前方往復一致、(b) 大局的妥当性のみ（北中緯度・米大陸経度）を確認する。
    /// NASA 秒のハードコードはしない（厳密オラクルではない）。
    #[test]
    fn real_2017_eclipse_roundtrips_and_is_plausible() {
        let dt = EspenakMeeusDeltaT;
        let iv = TimeInterval {
            start: TtInstant::from_jd2(JulianDate2::new(2_457_986.0, 0.0)),
            end: TtInstant::from_jd2(JulianDate2::new(2_457_988.0, 0.0)),
        };
        let src = DirectBesselianSource::new(R_SUN, R_MOON, &dt, iv);
        // 2017-08-21 最大食付近の TT（projection.rs / local_maximum.rs テストと同一エポック）。
        let tt = TtInstant::from_jd2(JulianDate2::new(2_457_986.5, 7.685_322_222_222_222e-1));
        let e = src
            .at(tt)
            .expect("real 2017 eclipse geometry should be valid");
        // 中心皆既なので gamma<1（軸は地球に当たる）。
        assert!(e.gamma() < 1.0, "2017 gamma={} should be <1", e.gamma());

        let p =
            shadow_axis_surface_point(&e, &WGS84).expect("2017 central eclipse ⇒ axis hits Earth");
        // (a) 厳密な前方往復。
        assert_forward_roundtrip(&p, &e, TOL_ROUNDTRIP);
        // (b) 大局的妥当性のみ（厳密オラクルではない）: 北中緯度・米大陸経度帯。
        let lat_deg = p.lat.degrees().0;
        let lon_deg = p.lon.degrees().0;
        assert!(
            (25.0..=50.0).contains(&lat_deg),
            "lat {lat_deg}° should be northern mid-latitude (Americas band)"
        );
        assert!(
            (-105.0..=-70.0).contains(&lon_deg),
            "lon {lon_deg}°E should be in the Americas"
        );
    }

    // ---- 要件3: 自明幾何（手導出） ----

    /// x=0, y=0, d=0, μ=1.0。前方射影式から (φ,λ) を**手で導出**する:
    ///   前方: H = μ + λ_east（signed 正規化）, 観測者高 0 で
    ///     ξ = ρcosφ′·sinH
    ///     η = ρsinφ′·cosd − ρcosφ′·sind·cosH
    ///     ζ = ρsinφ′·sind + ρcosφ′·cosd·cosH
    ///   d=0 を代入 → ξ = ρcosφ′·sinH, η = ρsinφ′, ζ = ρcosφ′·cosH。
    ///   要求 (ξ,η)=(0,0):
    ///     η = ρsinφ′ = 0 ⇒ φ′=0 ⇒ φ=0（赤道）。赤道海面では ρcosφ′=1。
    ///     ξ = sinH = 0 ⇒ H=0 または H=π。
    ///   太陽側 ζ>0: ζ = cosH > 0 ⇒ H=0（H=π は ζ=−1<0 で反太陽側ゆえ排除）。
    ///     H = μ + λ_east = 0 ⇒ λ_east = −μ = −1.0。−1.0 ∈ [−π,π) なので正規化後も −1.0。
    ///   ∴ 期待 φ≈0, λ≈−1.0。
    #[test]
    fn self_evident_geometry_phi_zero_lambda_minus_mu() {
        let mu = 1.0_f64;
        let e = elems(0.0, 0.0, 0.0, mu);
        let p = shadow_axis_surface_point(&e, &WGS84)
            .expect("axis through equator ⇒ surface point exists");
        let phi = p.lat.radians().0;
        let lam = p.lon.radians().0;
        // 手導出: φ=0（赤道）、λ_east = normalize(−μ) = −1.0。
        assert!(close(phi, 0.0, 1e-9), "φ={phi} expected 0 (equator)");
        let expected_lam = Radians::new(-mu).normalized_signed().0;
        assert!(
            close(lam, expected_lam, 1e-9),
            "λ={lam} expected normalize(−μ)={expected_lam}"
        );
        // 念のため往復・ζ>0 も確認（H=0 ⇒ ζ=+1）。
        let zeta = assert_forward_roundtrip(&p, &e, TOL_ROUNDTRIP);
        assert!(zeta > 0.0);
    }

    // ---- 要件4: 半球符号 ----

    /// y>0 ⇒ φ>0（北半球）、y<0 ⇒ φ<0（南半球）、同じ |y|。
    /// 独立幾何的期待: η は前方射影で ρsinφ′ 由来項（cosd>0）が主寄与で、同符号の緯度へ写る。
    /// （d 固定・x=0・同 μ で y の符号のみ反転。）
    #[test]
    fn hemisphere_sign_follows_y() {
        let d = 0.20_f64;
        let mu = 0.5_f64;
        let mag = 0.4_f64;
        let north = shadow_axis_surface_point(&elems(0.0, mag, d, mu), &WGS84)
            .expect("y>0 axis hits Earth");
        let south = shadow_axis_surface_point(&elems(0.0, -mag, d, mu), &WGS84)
            .expect("y<0 axis hits Earth");
        let phi_n = north.lat.radians().0;
        let phi_s = south.lat.radians().0;
        assert!(phi_n > 0.0, "y>0 should give northern φ, got {phi_n}");
        assert!(phi_s < 0.0, "y<0 should give southern φ, got {phi_s}");
        // 往復一致も維持（符号テストが偽の Ok を拾わないよう前方射影で裏取り）。
        assert_forward_roundtrip(&north, &elems(0.0, mag, d, mu), TOL_ROUNDTRIP);
        assert_forward_roundtrip(&south, &elems(0.0, -mag, d, mu), TOL_ROUNDTRIP);
    }

    // ---- 要件5: 軸が地球を外す ----

    /// gamma>1（x=1.5, y=0）⇒ 軸は地表に届かず `Err(Solver(RootNotBracketed))`。
    /// gamma>1 をテスト内で独立に確認する。
    #[test]
    fn axis_missing_earth_returns_root_not_bracketed() {
        let e = elems(1.5, 0.0, 0.20, 1.0);
        let gamma = (1.5_f64 * 1.5 + 0.0).sqrt();
        assert!(
            gamma > 1.0,
            "config gamma={gamma} must be >1 (axis misses Earth)"
        );
        let r = shadow_axis_surface_point(&e, &WGS84);
        assert!(
            matches!(
                r,
                Err(EclipseError::Solver(
                    umbra_core::SolverError::RootNotBracketed
                ))
            ),
            "expected Err(Solver(RootNotBracketed)), got {r:?}"
        );
    }

    // ---- 要件6: 太陽側選択 ----

    /// 太陽側（ζ>0）の交点を選ぶ。反太陽側を選ぶ変異を殺すため、往復で ζ>0 を明示確認。
    /// 反太陽側を選ぶと前方射影し直したとき ζ<0 になり `assert_forward_roundtrip` が落ちる。
    #[test]
    fn selects_sunward_intercept() {
        let e = elems(0.2, 0.1, 0.18, 3.0);
        let p = shadow_axis_surface_point(&e, &WGS84).expect("axis hits Earth");
        let zeta = assert_forward_roundtrip(&p, &e, TOL_ROUNDTRIP);
        // 太陽側選択を明示（assert_forward_roundtrip 内でも縛るが、要件6として独立に明示）。
        assert!(zeta > 0.0, "must select sunward intercept (ζ>0)");
    }

    // ---- 要件7: 経度正規化 ----

    /// 正規化前の素の λ_east（= H−μ 相当）が [−π,π) を外れうる構成で、返る λ_east ∈ [−π,π)。
    /// μ を大きく取り（μ=3.0, 軸直撃 x=y=0）、d=0 では λ_east=−μ=−3.0（これは [−π,π) 内）だが、
    /// 影軸直撃かつ μ=−3.0 にすると素の解は λ_east=+3.0（< π≈3.14159、ギリ内）。
    /// より確実に外側へ出すため μ=PI+0.5（H=0 解では λ=−(PI+0.5)=−3.6416 < −π）を使い、
    /// 正規化されて [−π,π) に収まることを縛る。
    #[test]
    fn longitude_is_normalized_into_half_open_interval() {
        // x=y=0, d=0 ⇒ 解は H=0 ⇒ λ_east=−μ。μ=PI+0.5 なら素の −μ=−3.6416<−π。
        let mu = PI + 0.5;
        let e = elems(0.0, 0.0, 0.0, mu);
        let p = shadow_axis_surface_point(&e, &WGS84).expect("axis through equator");
        let lam = p.lon.radians().0;
        assert!(
            (-PI..PI).contains(&lam),
            "λ_east={lam} must be normalized into [−π, π)"
        );
        // 正規化後も前方往復は成立しなければならない。
        assert_forward_roundtrip(&p, &e, TOL_ROUNDTRIP);
    }
}
