//! 観測者の基本面射影（`docs/issues/ISSUE-024`、`docs/conventions.md` §3/§5、
//! Explanatory Supplement to the Astronomical Almanac §11 / Meeus *Astronomical Algorithms* Ch.54）。
//!
//! 観測者（地心緯度成分 ρsinφ′/ρcosφ′ ＋ 東経 λ）を、ある TT における月影軸を z 軸とする
//! ベッセル基本面の座標 `(ξ, η, ζ)`（単位 Re）へ射影する。接触・最大食・食分 solver が
//! 時刻関数として呼ぶ最下層プリミティブ（本 issue は座標射影のみ）。
//!
//! `xi_rate, eta_rate` は「観測者が地球自転で運ばれることによる基本面速度」（影軸と d を固定し
//! 局地時角 H のみが地球自転角速度 μ′ で変化する成分, Re/SI秒）。影軸運動 (x′,y′)・d′ は本層に
//! 含めず、接触 solver が `d/dt(ξ−x)` として合成する（ISSUE-025）。

use umbra_core::constants::EARTH_ROTATION_RATE_RAD_PER_S;
use umbra_core::ellipsoid::GeocentricObserver;
use umbra_core::Radians;

use crate::besselian::InstantaneousBesselianElements;

/// 観測者の基本面座標（月影軸 z のベッセル系）。単位 Re。
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ObserverFundamental {
    /// 基本面 x̂（東）成分 ξ（Re）。
    pub xi: f64,
    /// 基本面 ŷ（天の北）成分 η（Re）。
    pub eta: f64,
    /// 影軸（太陽向き ẑ）成分 ζ（Re, 符号付き）。
    pub zeta: f64,
    /// dξ/dt（Re/SI秒, 地球自転による観測者速度のみ）。
    pub xi_rate: f64,
    /// dη/dt（Re/SI秒, 地球自転による観測者速度のみ）。
    pub eta_rate: f64,
}

/// 観測者を、瞬時ベッセル要素（d, μ）の基本面へ射影する。
///
/// - `observer`: ρsinφ′/ρcosφ′（WGS84 扁平込み, ISSUE-010/011・単位 Re）。
/// - `east_longitude`: 東経正（conventions §3）。
/// - `elements`: 瞬時ベッセル要素（赤緯 d・時角 μ を使用, ISSUE-021）。
pub fn project_observer_to_fundamental(
    observer: &GeocentricObserver,
    east_longitude: Radians,
    elements: &InstantaneousBesselianElements,
) -> ObserverFundamental {
    // 局地時角 H = μ − λ_east（東経正・conventions §3）。折返しで微分が壊れないよう
    // signed 正規化 [-π,π)（conventions §2）。三角関数は折返しに対し連続。
    let h = Radians::new(elements.mu.0 - east_longitude.0)
        .normalized_signed()
        .0;
    let (sin_h, cos_h) = h.sin_cos();
    let (sin_d, cos_d) = elements.declination.0.sin_cos();
    let rho_sin = observer.rho_sin_phi_prime;
    let rho_cos = observer.rho_cos_phi_prime;

    // Explanatory Supplement §11 / Meeus Ch.54 の局地予報式（単位 Re）。
    let xi = rho_cos * sin_h;
    let eta = rho_sin * cos_d - rho_cos * sin_d * cos_h;
    let zeta = rho_sin * sin_d + rho_cos * cos_d * cos_h;

    // 観測者が地球自転で運ばれる基本面速度（影軸と d を固定し H のみ μ′ で変化）。
    // μ′ = dERA/dt = 地球自転角速度（影軸赤経固定）。Re/SI秒。
    // 影軸運動 (x′,y′)・d′ は本層に含めない（接触 solver が d/dt(ξ−x) として合成・ISSUE-025）。
    let mu_rate = EARTH_ROTATION_RATE_RAD_PER_S;
    let xi_rate = rho_cos * cos_h * mu_rate;
    let eta_rate = rho_cos * sin_d * sin_h * mu_rate;

    ObserverFundamental {
        xi,
        eta,
        zeta,
        xi_rate,
        eta_rate,
    }
}

#[cfg(test)]
mod tests {
    //! ISSUE-024 受け入れテスト（strict）。
    //!
    //! オラクル戦略（追認回避）:
    //! - 不変量 ξ²+η²+ζ² = ρ²（純回転は長さ保存。成分式の内訳から完全独立）。
    //! - 既知配置の幾何的自明値（影軸直下で ξ=η≈0、赤道/極/H=±π/2 の極端配置）。
    //! - 既知地点（岡山）はテスト内に「契約式の独立再実装 (`oracle_xez`)」を置き、
    //!   入力（ρsinφ′/ρcosφ′ は `observer_geocentric` 由来, μ は手で与える）から
    //!   別経路で (ξ,η,ζ) を組む。`project_*` の戻り値と一致を見る。
    //!   ※「実装本体をコピーした追認」を避けるため、独立再実装は
    //!   回転行列形 R_x(d)·R_z(−H)·(ρ ベクトル) でも同じ値になることを別途確認する
    //!   （`okayama_matches_rotation_matrix_oracle`）。
    //! - 微分はベッセル要素 μ を ±δ 振って `project` を 2 回呼ぶ中心差分（d, ρ 固定）と、
    //!   解析 rate / μ′ の一致で検証。μ′ は IERS/SOFA era00 の dERA/dt を独立に導出。
    #![allow(clippy::excessive_precision)]

    use super::*;
    use core::f64::consts::{PI, TAU};
    use umbra_core::constants::EARTH_EQUATORIAL_RADIUS_M;
    use umbra_core::constants::SOLAR_RADIUS_KM;
    use umbra_core::ellipsoid::{observer_geocentric, Ellipsoid};
    use umbra_core::{JulianDate2, TtInstant};
    use umbra_ephemeris::{Body, Ephemeris, EphemerisFrame, MockEphemeris, Origin};

    use crate::besselian::{besselian_elements, besselian_elements_at};
    use umbra_core::EspenakMeeusDeltaT;

    const WGS84: Ellipsoid = Ellipsoid::WGS84;
    const R_SUN: f64 = SOLAR_RADIUS_KM;
    const R_MOON: f64 = 0.2725076 * (EARTH_EQUATORIAL_RADIUS_M / 1000.0);

    /// 地球自転角速度 μ′ = dERA/dt（rad/SI秒）。
    /// 出典: IERS/SOFA `iauEra00`。dERA/dt = 2π · 1.00273781191135448 / 86400。
    /// ベッセル要素の影軸を固定すれば μ′ = dERA/dt（α_axis は時間不変扱い）。独立導出。
    const D_ERA_DT: f64 = TAU * 1.002_737_811_911_354_48 / 86_400.0;

    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    /// テスト用に瞬時ベッセル要素を最小構成する（d=declination, μ=mu のみが射影で使われる）。
    /// 他フィールドは射影に無関係なのでダミー（有限値）を入れる。
    fn elems_with(d: f64, mu: f64) -> InstantaneousBesselianElements {
        InstantaneousBesselianElements {
            x: 0.0,
            y: 0.0,
            declination: Radians(d),
            mu: Radians(mu),
            l1: 0.5,
            l2: -0.01,
            tan_f1: 0.0047,
            tan_f2: 0.0046,
            time_tt: TtInstant::from_jd2(JulianDate2::new(2_451_545.0, 0.0)),
        }
    }

    /// 契約式の独立再実装（オラクル）。H=μ−λ を signed 正規化し (ξ,η,ζ) を返す。
    /// これは ISSUE-024 §座標定義の式を直接書き下したもの。実装本体とは別ファイル・別ロジック
    /// （実装は todo!()。本関数は本テストモジュール内の独立オラクル）。
    fn oracle_xez(rho_sin: f64, rho_cos: f64, d: f64, mu: f64, lambda: f64) -> (f64, f64, f64) {
        let h = Radians::new(mu - lambda).normalized_signed().0;
        let xi = rho_cos * h.sin();
        let eta = rho_sin * d.cos() - rho_cos * d.sin() * h.cos();
        let zeta = rho_sin * d.sin() + rho_cos * d.cos() * h.cos();
        (xi, eta, zeta)
    }

    // ---- 不変量（純回転は長さを保存。成分式から完全独立） ----

    /// ξ²+η²+ζ² = (ρsinφ′)²+(ρcosφ′)² = ρ²。任意 (lat,h,d,μ,λ) で成立すべき回転不変量。
    #[test]
    fn norm_squared_equals_rho_squared() {
        let cases = [
            (35.0_f64, 0.0, 0.20, 1.0, 0.5),
            (-42.0, 4000.0, -0.35, 2.7, -1.3),
            (0.0, 0.0, 0.0, 0.0, 0.0),
            (66.5, 1000.0, 0.45, -3.0, 3.0),
        ];
        for (lat_deg, h, d, mu, lam) in cases {
            let obs = observer_geocentric(&WGS84, lat_deg.to_radians(), h);
            let r = project_observer_to_fundamental(&obs, Radians::new(lam), &elems_with(d, mu));
            let rho2 = obs.rho_sin_phi_prime.powi(2) + obs.rho_cos_phi_prime.powi(2);
            let got = r.xi * r.xi + r.eta * r.eta + r.zeta * r.zeta;
            assert!(
                close(got, rho2, 1e-12),
                "‖(ξ,η,ζ)‖² = {got} expected ρ² = {rho2}"
            );
        }
    }

    // ---- 既知地点（岡山）独立オラクル一致 ----

    /// 岡山 34.507°N, 133.508°E, 10 m。独立オラクル `oracle_xez`（契約式の別実装）と一致。
    /// ρ 成分は `observer_geocentric`（WGS84 扁平込み）由来。μ は任意固定値を与える。
    #[test]
    fn okayama_matches_independent_oracle() {
        let lat = 34.507_f64.to_radians();
        let lon = 133.508_f64.to_radians();
        let h = 10.0;
        let (d, mu) = (0.2070_f64, 2.0_f64); // d≈+11.86°（2017 太陽見かけ赤緯相当）, μ 任意
        let obs = observer_geocentric(&WGS84, lat, h);
        let r = project_observer_to_fundamental(&obs, Radians::new(lon), &elems_with(d, mu));
        let (xi, eta, zeta) = oracle_xez(obs.rho_sin_phi_prime, obs.rho_cos_phi_prime, d, mu, lon);
        assert!(close(r.xi, xi, 1e-12), "ξ={} oracle={}", r.xi, xi);
        assert!(close(r.eta, eta, 1e-12), "η={} oracle={}", r.eta, eta);
        assert!(close(r.zeta, zeta, 1e-12), "ζ={} oracle={}", r.zeta, zeta);
    }

    /// 岡山の (ξ,η,ζ) を回転行列形 R_x(d)·R_z(−H)·(ρcosφ′,0,ρsinφ′) で組んだ第二オラクルと一致。
    /// 観測者地心ベクトルを「H 経度・φ′ 緯度」の球面ベクトルとして赤道系に置き、影軸系へ回す
    /// 経路。契約式の三角恒等式とは別表現なので追認にならない。
    #[test]
    fn okayama_matches_rotation_matrix_oracle() {
        let lat = 34.507_f64.to_radians();
        let lon = 133.508_f64.to_radians();
        let h = 10.0;
        let (d, mu) = (0.2070_f64, 2.0_f64);
        let obs = observer_geocentric(&WGS84, lat, h);
        let r = project_observer_to_fundamental(&obs, Radians::new(lon), &elems_with(d, mu));

        // 観測者地心ベクトル（赤道系, Re）: 時角 H からの赤経は (−H) 相当。
        // 影軸を z とする基本面は、赤道系を「赤経=影軸赤経」周りに回し、さらに赤緯 d だけ傾ける。
        // ここでは H = μ−λ を局地時角とし、赤道系で観測者を (ρcosφ′cosH', ρcosφ′ sinH', ρsinφ′)
        // のように置いた上で R_x(d) 回転で基本面 (ξ,η,ζ) を作る独立経路。
        let hh = Radians::new(mu - lon).normalized_signed().0;
        let rc = obs.rho_cos_phi_prime;
        let rs = obs.rho_sin_phi_prime;
        // 赤道系の中間座標（x 軸=影軸子午面方向, y 軸=東, z 軸=天の北）:
        //   p = (ρcosφ′·cosH, ρcosφ′·sinH, ρsinφ′)
        let p = (rc * hh.cos(), rc * hh.sin(), rs);
        // 影軸系へ: ξ=東成分=p.y、η・ζ は (p.x,p.z) を赤緯 d で回す。
        //   η =  −p.x·sin d + p.z·cos d
        //   ζ =   p.x·cos d + p.z·sin d
        let xi = p.1;
        let eta = -p.0 * d.sin() + p.2 * d.cos();
        let zeta = p.0 * d.cos() + p.2 * d.sin();
        assert!(close(r.xi, xi, 1e-12), "ξ={} rot={}", r.xi, xi);
        assert!(close(r.eta, eta, 1e-12), "η={} rot={}", r.eta, eta);
        assert!(close(r.zeta, zeta, 1e-12), "ζ={} rot={}", r.zeta, zeta);
    }

    // ---- 既知配置の幾何的自明値 ----

    /// 影軸直下の観測者（H=0, φ′ が d に一致する向き）で ξ=0, η=0, ζ=ρ。
    /// d=0 の赤道配置では、赤道上 (lat=0,h=0) かつ H=0（μ=λ）で観測者が影軸直下。
    #[test]
    fn observer_under_axis_has_zero_xi_eta() {
        let d = 0.0;
        let lambda = 0.7_f64;
        let mu = lambda; // H = μ − λ = 0
        let obs = observer_geocentric(&WGS84, 0.0, 0.0); // 赤道海面: ρcos=1, ρsin=0
        let r = project_observer_to_fundamental(&obs, Radians::new(lambda), &elems_with(d, mu));
        assert!(r.xi.abs() < 1e-12, "ξ={}", r.xi);
        assert!(r.eta.abs() < 1e-12, "η={}", r.eta);
        assert!(close(r.zeta, 1.0, 1e-12), "ζ={} (赤道海面 ρ=1)", r.zeta);
    }

    /// H=+π/2（観測者が影軸の東 90°）で ξ=+ρcosφ′（最大正）, η=−ρsinφ′·sin d。
    /// ξ の符号（東正）と sinH の効きを幾何的に固定する。
    #[test]
    fn east_quadrature_maximizes_xi() {
        let d = 0.3_f64;
        let lambda = 0.0_f64;
        let mu = PI / 2.0; // H = +π/2
        let obs = observer_geocentric(&WGS84, 20.0_f64.to_radians(), 0.0);
        let r = project_observer_to_fundamental(&obs, Radians::new(lambda), &elems_with(d, mu));
        assert!(close(r.xi, obs.rho_cos_phi_prime, 1e-12), "ξ={}", r.xi);
        // cos H = 0 なので η = ρsinφ′·cos d, ζ = ρsinφ′·sin d。
        assert!(
            close(r.eta, obs.rho_sin_phi_prime * d.cos(), 1e-12),
            "η={}",
            r.eta
        );
        assert!(
            close(r.zeta, obs.rho_sin_phi_prime * d.sin(), 1e-12),
            "ζ={}",
            r.zeta
        );
    }

    /// 北極の観測者（ρcosφ′≈0, ρsinφ′≈b/a）。ξ≈0, η≈ρsinφ′·cos d, ζ≈ρsinφ′·sin d（H 非依存）。
    #[test]
    fn north_pole_observer_is_hour_angle_independent() {
        let d = 0.25_f64;
        let obs = observer_geocentric(&WGS84, PI / 2.0, 0.0);
        let r1 = project_observer_to_fundamental(&obs, Radians::new(0.3), &elems_with(d, 1.0));
        let r2 = project_observer_to_fundamental(&obs, Radians::new(2.1), &elems_with(d, 4.5));
        // ρcosφ′≈0 → ξ≈0、H に依存しない。
        assert!(r1.xi.abs() < 1e-9, "ξ1={}", r1.xi);
        assert!(
            close(r1.eta, r2.eta, 1e-9),
            "η H依存: {} vs {}",
            r1.eta,
            r2.eta
        );
        assert!(
            close(r1.zeta, r2.zeta, 1e-9),
            "ζ H依存: {} vs {}",
            r1.zeta,
            r2.zeta
        );
        assert!(
            close(r1.eta, obs.rho_sin_phi_prime * d.cos(), 1e-9),
            "η={}",
            r1.eta
        );
        assert!(
            close(r1.zeta, obs.rho_sin_phi_prime * d.sin(), 1e-9),
            "ζ={}",
            r1.zeta
        );
    }

    // ---- MockEphemeris 貫通（暦→ベッセル要素→射影） ----

    /// central_total（影軸=+x, d≈0）で、赤道上・影軸直下経度（μ=λ）の観測者は ξ=η≈0。
    /// 暦→ベッセル要素の実チェーンを通す。μ は central_total では d=0, 影軸赤経=0 なので
    /// 観測者を影軸直下に置くには λ=μ（μ は時刻依存）。ここでは μ を直接与える等価構成にする。
    #[test]
    fn central_total_under_axis_gives_zero_transverse() {
        let m = MockEphemeris::central_total();
        let t = umbra_core::TdbInstant::from_jd2(JulianDate2::from_jd(2_451_545.0));
        let pos = |b| {
            m.state(b, t, Origin::Geocenter, EphemerisFrame::Icrs)
                .unwrap()
                .position
        };
        let geom = besselian_elements(pos(Body::Sun), pos(Body::Moon), R_SUN, R_MOON).unwrap();
        // central_total: 影軸 +x, d≈0。観測者を赤道海面・影軸直下（H=0）に置く。
        let d = geom.declination.0;
        assert!(d.abs() < 1e-6, "mock d should be ~0: {d}");
        let lambda = 1.234_f64;
        let elems = elems_with(d, lambda); // μ=λ → H=0
        let obs = observer_geocentric(&WGS84, 0.0, 0.0);
        let r = project_observer_to_fundamental(&obs, Radians::new(lambda), &elems);
        assert!(r.xi.abs() < 1e-9, "ξ={}", r.xi);
        assert!(r.eta.abs() < 1e-9, "η={}", r.eta);
        assert!(r.zeta > 0.0, "ζ={} (影軸前方=正)", r.zeta);
    }

    /// 影軸から既知横ずれ点で √(ξ²+η²) が独立計算と一致。
    /// 赤道海面 (ρ=1) で H をずらすと √(ξ²+η²) = √(sin²H + sin²d·cos²H)（d=0 なら |sinH|）。
    #[test]
    fn known_offset_transverse_distance_matches_oracle() {
        let d = 0.0_f64;
        let lambda = 0.0_f64;
        let h_angle = 0.4_f64; // 影軸から東へ 0.4 rad
        let obs = observer_geocentric(&WGS84, 0.0, 0.0); // ρ=1
        let r =
            project_observer_to_fundamental(&obs, Radians::new(lambda), &elems_with(d, h_angle));
        let m = (r.xi * r.xi + r.eta * r.eta).sqrt();
        // d=0, ρ=1: 独立オラクル m = |sin H|。
        assert!(
            close(m, h_angle.sin().abs(), 1e-12),
            "m={} expected {}",
            m,
            h_angle.sin()
        );
    }

    // ---- 標高差: ζ が標高分だけ変化 ----

    /// h=0 と h=4000 m で ζ が標高分変化（高い点ほど影軸前方＝ζ 増。中緯度・影軸近傍）。
    /// 独立オラクル: Δζ = (Δρsin·sin d + Δρcos·cos d·cos H)。ρ 増分は observer_geocentric 由来。
    #[test]
    fn altitude_changes_zeta_by_height_term() {
        let lat = 40.0_f64.to_radians();
        let lambda = 0.0_f64;
        let (d, mu) = (0.2_f64, 0.1_f64); // 影軸近傍（小 H）
        let o0 = observer_geocentric(&WGS84, lat, 0.0);
        let o4 = observer_geocentric(&WGS84, lat, 4000.0);
        let r0 = project_observer_to_fundamental(&o0, Radians::new(lambda), &elems_with(d, mu));
        let r4 = project_observer_to_fundamental(&o4, Radians::new(lambda), &elems_with(d, mu));
        let h = Radians::new(mu - lambda).normalized_signed().0;
        let d_rs = o4.rho_sin_phi_prime - o0.rho_sin_phi_prime;
        let d_rc = o4.rho_cos_phi_prime - o0.rho_cos_phi_prime;
        let expected_dz = d_rs * d.sin() + d_rc * d.cos() * h.cos();
        assert!(
            close(r4.zeta - r0.zeta, expected_dz, 1e-12),
            "Δζ={} expected {}",
            r4.zeta - r0.zeta,
            expected_dz
        );
        // 標高があるほど ζ は前方へ（影軸近傍 cos d·cos H>0, sin d>0 なので正）。
        assert!(r4.zeta > r0.zeta, "ζ4={} ζ0={}", r4.zeta, r0.zeta);
    }

    // ---- 西経入力の吸収 ----

    /// 東経 260° ≡ 西経 −100°。同じ μ・観測者で (ξ,η,ζ) 一致。
    #[test]
    fn west_longitude_absorbs_to_east() {
        let lat = 30.0_f64.to_radians();
        let (d, mu) = (0.1_f64, 1.0_f64);
        let obs = observer_geocentric(&WGS84, lat, 0.0);
        let east = 260.0_f64.to_radians();
        let west = (-100.0_f64).to_radians();
        let re = project_observer_to_fundamental(&obs, Radians::new(east), &elems_with(d, mu));
        let rw = project_observer_to_fundamental(&obs, Radians::new(west), &elems_with(d, mu));
        assert!(close(re.xi, rw.xi, 1e-12), "ξ {} vs {}", re.xi, rw.xi);
        assert!(close(re.eta, rw.eta, 1e-12), "η {} vs {}", re.eta, rw.eta);
        assert!(
            close(re.zeta, rw.zeta, 1e-12),
            "ζ {} vs {}",
            re.zeta,
            rw.zeta
        );
        assert!(close(re.xi_rate, rw.xi_rate, 1e-12));
        assert!(close(re.eta_rate, rw.eta_rate, 1e-12));
    }

    // ---- プロパティ: 周期 2π 不変 ----

    /// λ→λ+2π で (ξ,η,ζ,rate) 不変。
    #[test]
    fn invariant_under_lambda_plus_two_pi() {
        let lat = 25.0_f64.to_radians();
        let (d, mu) = (0.15_f64, 2.3_f64);
        let obs = observer_geocentric(&WGS84, lat, 100.0);
        let lam = 0.9_f64;
        let a = project_observer_to_fundamental(&obs, Radians::new(lam), &elems_with(d, mu));
        let b = project_observer_to_fundamental(&obs, Radians::new(lam + TAU), &elems_with(d, mu));
        assert!(
            close(a.xi, b.xi, 1e-12) && close(a.eta, b.eta, 1e-12) && close(a.zeta, b.zeta, 1e-12)
        );
        assert!(close(a.xi_rate, b.xi_rate, 1e-12) && close(a.eta_rate, b.eta_rate, 1e-12));
    }

    /// μ→μ+2π で (ξ,η,ζ,rate) 不変。
    #[test]
    fn invariant_under_mu_plus_two_pi() {
        let lat = 25.0_f64.to_radians();
        let (d, mu) = (0.15_f64, 2.3_f64);
        let obs = observer_geocentric(&WGS84, lat, 100.0);
        let lam = 0.9_f64;
        let a = project_observer_to_fundamental(&obs, Radians::new(lam), &elems_with(d, mu));
        let b = project_observer_to_fundamental(&obs, Radians::new(lam), &elems_with(d, mu + TAU));
        assert!(
            close(a.xi, b.xi, 1e-12) && close(a.eta, b.eta, 1e-12) && close(a.zeta, b.zeta, 1e-12)
        );
        assert!(close(a.xi_rate, b.xi_rate, 1e-12) && close(a.eta_rate, b.eta_rate, 1e-12));
    }

    // ---- 微分セマンティクス ----

    /// xi_rate / μ′ が ∂ξ/∂H の中心差分に一致（μ を ±δ 振り d, ρ 固定）。
    /// μ′ = D_ERA_DT（独立導出）。`project` の rate は μ′ を内包するので /μ′ で ∂ξ/∂H を取り出す。
    #[test]
    fn xi_rate_matches_central_difference_of_h() {
        let lat = 37.0_f64.to_radians();
        let (d, mu) = (0.22_f64, 1.7_f64);
        let lam = 0.5_f64;
        let obs = observer_geocentric(&WGS84, lat, 0.0);
        let delta = 1e-6_f64; // μ の微小振り（rad）
        let r_plus =
            project_observer_to_fundamental(&obs, Radians::new(lam), &elems_with(d, mu + delta));
        let r_minus =
            project_observer_to_fundamental(&obs, Radians::new(lam), &elems_with(d, mu - delta));
        let dxi_dh = (r_plus.xi - r_minus.xi) / (2.0 * delta);
        let r0 = project_observer_to_fundamental(&obs, Radians::new(lam), &elems_with(d, mu));
        // rate = (∂ξ/∂H)·μ′ → ∂ξ/∂H = rate / μ′。
        let dxi_dh_from_rate = r0.xi_rate / D_ERA_DT;
        let rel = (dxi_dh - dxi_dh_from_rate).abs() / dxi_dh.abs().max(1e-12);
        assert!(
            rel <= 1e-6,
            "∂ξ/∂H num={} rate/μ′={} rel={}",
            dxi_dh,
            dxi_dh_from_rate,
            rel
        );
    }

    /// eta_rate / μ′ が ∂η/∂H の中心差分に一致。
    #[test]
    fn eta_rate_matches_central_difference_of_h() {
        let lat = 37.0_f64.to_radians();
        let (d, mu) = (0.22_f64, 1.7_f64);
        let lam = 0.5_f64;
        let obs = observer_geocentric(&WGS84, lat, 0.0);
        let delta = 1e-6_f64;
        let r_plus =
            project_observer_to_fundamental(&obs, Radians::new(lam), &elems_with(d, mu + delta));
        let r_minus =
            project_observer_to_fundamental(&obs, Radians::new(lam), &elems_with(d, mu - delta));
        let deta_dh = (r_plus.eta - r_minus.eta) / (2.0 * delta);
        let r0 = project_observer_to_fundamental(&obs, Radians::new(lam), &elems_with(d, mu));
        let deta_dh_from_rate = r0.eta_rate / D_ERA_DT;
        let rel = (deta_dh - deta_dh_from_rate).abs() / deta_dh.abs().max(1e-12);
        assert!(
            rel <= 1e-6,
            "∂η/∂H num={} rate/μ′={} rel={}",
            deta_dh,
            deta_dh_from_rate,
            rel
        );
    }

    /// rate の解析オラクル一致: xi_rate = ρcosφ′·cosH·μ′, eta_rate = ρcosφ′·sin d·sinH·μ′。
    /// 契約に明記された解析式（中心差分とは別の閉形式）で符号・係数を固定する。
    #[test]
    fn rates_match_analytic_oracle() {
        let lat = 48.0_f64.to_radians();
        let (d, mu) = (-0.3_f64, 4.2_f64);
        let lam = 2.0_f64;
        let obs = observer_geocentric(&WGS84, lat, 500.0);
        let r = project_observer_to_fundamental(&obs, Radians::new(lam), &elems_with(d, mu));
        let h = Radians::new(mu - lam).normalized_signed().0;
        let xi_rate = obs.rho_cos_phi_prime * h.cos() * D_ERA_DT;
        let eta_rate = obs.rho_cos_phi_prime * d.sin() * h.sin() * D_ERA_DT;
        assert!(
            close(r.xi_rate, xi_rate, 1e-15),
            "ξ′={} expected {}",
            r.xi_rate,
            xi_rate
        );
        assert!(
            close(r.eta_rate, eta_rate, 1e-15),
            "η′={} expected {}",
            r.eta_rate,
            eta_rate
        );
    }

    // ---- 退化境界の有限性 ----

    /// 日付変更線付近（H が ±π 折返し近傍）でも全出力が有限・連続（NaN/Inf を出さない）。
    #[test]
    fn outputs_finite_near_antimeridian() {
        let lat = 10.0_f64.to_radians();
        let obs = observer_geocentric(&WGS84, lat, 0.0);
        for mu in [PI - 1e-9, PI + 1e-9, -PI + 1e-9] {
            let r = project_observer_to_fundamental(&obs, Radians::new(0.0), &elems_with(0.1, mu));
            for v in [r.xi, r.eta, r.zeta, r.xi_rate, r.eta_rate] {
                assert!(v.is_finite(), "non-finite at mu={mu}: {v}");
            }
        }
    }

    // ---- 実チェーン（暦→ΔT→μ）由来要素での μ 貫通 ----

    /// 合成 μ ではなく `besselian_elements_at`（apparent 位置→影軸→ΔT→μ の実チェーン）由来の
    /// 瞬時要素で射影し、独立オラクル `oracle_xez` と一致することを確認する。
    /// 2017-08-21 最大食付近の実エポックで評価（実太陽・月位置が日食配置で有効）。
    /// μ は実 ERA(UT1)−α_axis 由来の値が射影へ正しく貫通するかを縛る（`elems_with` の手与え μ では
    /// 通らない経路）。観測点は岡山経度。`EspenakMeeusDeltaT` を ΔT モデルに用いる。
    #[test]
    fn projection_uses_real_chain_mu_from_elements_at() {
        // 2017-08-21 最大食付近の TT（besselian.rs の tt_2017_max と同一エポック）。
        let tt = TtInstant::from_jd2(JulianDate2::new(2_457_986.5, 7.685_322_222_222_221_72e-1));
        let elems = besselian_elements_at(tt, R_SUN, R_MOON, &EspenakMeeusDeltaT)
            .expect("real apparent positions should yield valid eclipse geometry");

        let lambda = 133.508_f64.to_radians(); // 岡山経度（任意の既知地点でよい）
        let obs = observer_geocentric(&WGS84, 34.507_f64.to_radians(), 10.0);
        let r = project_observer_to_fundamental(&obs, Radians::new(lambda), &elems);

        // 実チェーン由来の d, μ を取り出して独立オラクルへ。
        let (xi, eta, zeta) = oracle_xez(
            obs.rho_sin_phi_prime,
            obs.rho_cos_phi_prime,
            elems.declination.0,
            elems.mu.0,
            lambda,
        );
        assert!(close(r.xi, xi, 1e-12), "ξ={} oracle={}", r.xi, xi);
        assert!(close(r.eta, eta, 1e-12), "η={} oracle={}", r.eta, eta);
        assert!(close(r.zeta, zeta, 1e-12), "ζ={} oracle={}", r.zeta, zeta);
    }

    // ---- H=±π 折返しでの rate 連続性 ----

    /// μ=π−ε と μ=π+ε（H が ±π 折返しをまたぐ）で xi_rate・eta_rate がそれぞれ連続。
    /// `normalized_signed` を rate 計算に使っても、cos/sin は H=±π で連続なので微分は壊れない。
    /// d, ρ, λ は固定。折返しで rate が不連続にジャンプ（例: H を [0,2π) に取った実装の取りこぼし）
    /// しないことを縛る。
    #[test]
    fn rates_continuous_across_antimeridian_wrap() {
        let eps = 1e-7_f64;
        let lat = 30.0_f64.to_radians();
        let lambda = 0.0_f64;
        let d = 0.2_f64;
        let obs = observer_geocentric(&WGS84, lat, 0.0);
        // λ=0 なので H=μ。μ=π∓ε で H が +π / −π 側の折返しをまたぐ。
        let r_minus =
            project_observer_to_fundamental(&obs, Radians::new(lambda), &elems_with(d, PI - eps));
        let r_plus =
            project_observer_to_fundamental(&obs, Radians::new(lambda), &elems_with(d, PI + eps));
        assert!(
            close(r_minus.xi_rate, r_plus.xi_rate, 1e-6),
            "xi_rate 不連続: {} vs {}",
            r_minus.xi_rate,
            r_plus.xi_rate
        );
        assert!(
            close(r_minus.eta_rate, r_plus.eta_rate, 1e-6),
            "eta_rate 不連続: {} vs {}",
            r_minus.eta_rate,
            r_plus.eta_rate
        );
    }

    // ---- d≠0 の横ずれ距離オラクル ----

    /// d≠0・赤道海面 ρ=1 で √(ξ²+η²) が独立式 √(sin²H + sin²d·cos²H) と一致。
    /// `known_offset_transverse_distance_matches_oracle`（d=0 限定）の補完。d≠0 では
    /// η = ρsinφ′·cosd − ρcosφ′·sind·cosH の sind·cosH 項が m に効く（ρ=1, φ′=0 なので
    /// ρsinφ′=0, ρcosφ′=1 → η = −sind·cosH, ξ = sinH）。この sind 寄与の取りこぼしを縛る。
    #[test]
    fn known_offset_transverse_distance_with_nonzero_d() {
        let d = 0.2_f64;
        let lambda = 0.0_f64;
        let h_angle = 0.4_f64; // 影軸から東へ 0.4 rad（H=μ−λ=μ）
        let obs = observer_geocentric(&WGS84, 0.0, 0.0); // 赤道海面 ρ=1
        let r =
            project_observer_to_fundamental(&obs, Radians::new(lambda), &elems_with(d, h_angle));
        let m = (r.xi * r.xi + r.eta * r.eta).sqrt();
        // ρ=1, φ′=0: ξ=sinH, η=−sind·cosH → m = √(sin²H + sin²d·cos²H)。
        let expected = (h_angle.sin().powi(2) + d.sin().powi(2) * h_angle.cos().powi(2)).sqrt();
        assert!(close(m, expected, 1e-12), "m={m} expected {expected}");
    }
}
