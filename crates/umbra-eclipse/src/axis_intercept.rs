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

use core::f64::consts::TAU;

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

/// 限界点の自己整合（相対速度包絡）不動点反復の収束許容（基本面 Re）。中心線位置 sub-km
/// （≈7.8e-5 Re）を大きく下回り、`L2'=l2−ζ·tan f2` の ζ 依存残差（tan f2·Δζ）を 1e-7 以下にする。
const LIMIT_FIXED_POINT_TOL: f64 = 1e-12;
/// 限界点 不動点反復の上限。半径 |L2'|≈0.01 Re は微小で ζ 依存が弱く、通常 3–4 回で収束。
const LIMIT_FIXED_POINT_MAX_ITER: usize = 16;
/// WGS84 平均半径 \[km\]（IUGG・帯幅の大圏距離換算用。±200 km 規模に対し球近似誤差 <0.5%）。
const EARTH_MEAN_RADIUS_KM: f64 = 6371.0;
/// rise/set limb 点（錐縁∩terminator 楕円）の媒介角 θ∈[0,2π) 粗走査の分割数（≤4 根を分離する解像度）。
const TERMINATOR_SCAN_STEPS: usize = 720;
/// terminator 交点の θ 根の収束許容（rad）。中心線位置 sub-km を十分下回る。
const TERMINATOR_ROOT_TOL: f64 = 1e-12;
/// terminator 交点 θ 求根の Brent 反復上限（円∩楕円の滑らかな残差ゆえ十分な余裕）。
const TERMINATOR_ROOT_MAX_ITER: usize = 100;

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
    surface_point_for_fundamental(
        elements.x,
        elements.y,
        elements.declination,
        elements.mu,
        ellipsoid,
    )
    .map(|(point, _zeta)| point)
}

/// 基本面交点 (ξ, η) を通る影軸方向の直線が WGS84 地表（太陽側 ζ>0）を貫く点と、その ζ を返す。
///
/// 影軸交点 (x,y) を渡せば中心線（[`shadow_axis_surface_point`]）、本影縁の点 (ξ,η) を渡せば限界線
/// （M9.3・`EclipseEngine::path`）の地表点になる。ζ は呼び出し側が ζ補正本影半径 `L2'=l2−ζ·tan f2` 等に
/// 使う。軸が地表を外す（実根なし）と `Err(Solver(RootNotBracketed))`。
pub(crate) fn surface_point_for_fundamental(
    xi: f64,
    eta: f64,
    declination: Radians,
    mu: Radians,
    ellipsoid: &Ellipsoid,
) -> Result<(GeoPoint, f64), EclipseError> {
    let (sin_d, cos_d) = declination.0.sin_cos();
    let omf = 1.0 - ellipsoid.f; // b/a
    let inv_omf2 = 1.0 / (omf * omf); // 1/(1−f)²

    // 逆回転（前方 R_x(d)·R_z の転置）で観測者地心成分を ζ の関数に:
    //   px = ζ·cos d − η·sin d,  py = ξ,  pz = ζ·sin d + η·cos d
    //   ρcos² = px² + py²,  ρsin = pz
    // h=0 の子午線楕円拘束（ellipsoid.rs 不変量）: ρcos² + (ρsin/(1−f))² = 1。
    // 残差 r(ζ) = ρcos² + (ρsin/(1−f))² − 1。零点が直線の地表交点。
    let residual = |zeta: f64| -> f64 {
        let px = zeta * cos_d - eta * sin_d;
        let pz = zeta * sin_d + eta * cos_d;
        let rho_cos2 = px * px + xi * xi;
        rho_cos2 + pz * pz * inv_omf2 - 1.0
    };

    // 太陽側（最大 ζ>0）根のブラケットを降順粗走査で取り、Brent で精解する。
    // r は ζ² の係数 cos²d + sin²d/(1−f)² > 0 ゆえ下に凸（開口上向き）で、2 根の間で負・外で正
    // なので、上端から下る最初の符号反転が大きい方の根 ζ₊（＝太陽側）。直線が地表を外す（実根なし）と
    // 符号反転が無くブラケット不成立 → `RootNotBracketed`。
    let (a, b) = descending_sign_change_bracket(&residual, ZETA_SCAN_MAX, ZETA_SCAN_STEPS)
        .ok_or(EclipseError::Solver(SolverError::RootNotBracketed))?;
    let zeta = brent_root(residual, a, b, ZETA_ROOT_TOL, ZETA_ROOT_MAX_ITER)?;

    // 零点 ζ・基本面交点 (ξ,η) から測地座標を復元。
    let point = fundamental_to_geodetic(xi, eta, zeta, declination, mu, ellipsoid)?;
    Ok((point, zeta))
}

/// 基本面座標 (ξ, η, ζ) と影軸の赤緯 d・時角 μ から測地座標（φ, λ_east）を復元する。
///
/// 前方射影 [`project_observer_to_fundamental`](crate::projection) の回転 `R_x(d)·R_z` の**転置**で
/// 観測者地心成分を組み（px=ζcosd−η sind, py=ξ, pz=ζsind+η cosd）、ρcos=√(px²+py²)・ρsin=pz・
/// 局地時角 H=atan2(py,px) を得て、東経 λ=H−μ（[-π,π) 正規化）・測地緯度
/// φ=atan2(ρsin, ρcos·(1−f)²)（`observer_geocentric(h=0)` の逆）を返す。ζ を与える側が決める
/// （影軸の地表貫通点 S6a-i は ζ 根、全球接触の地球縁点 S6b-ii は ζ=0）。
///
/// ρcos ≥ 0・(1−f)²>0 なので φ の atan2 は [-π/2, π/2] を返し `GeodeticLatitude` 検証は常に通る。
pub(crate) fn fundamental_to_geodetic(
    xi: f64,
    eta: f64,
    zeta: f64,
    declination: Radians,
    mu: Radians,
    ellipsoid: &Ellipsoid,
) -> Result<GeoPoint, EclipseError> {
    let (sin_d, cos_d) = declination.0.sin_cos();
    let omf = 1.0 - ellipsoid.f; // b/a
    let px = zeta * cos_d - eta * sin_d;
    let py = xi;
    let pz = zeta * sin_d + eta * cos_d;
    let rho_cos = (px * px + py * py).sqrt();
    let rho_sin = pz;
    let h = py.atan2(px);
    let longitude = EastLongitude::from_radians(Radians::new(h - mu.0));
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

/// 移動する影錐縁の片側（`sign=±1`）の縁点を相対速度包絡の不動点反復で厳密に解く（南北限界線・帯幅・
/// 部分食域）。錐半径は `cone_l`（ζ=0 での錐半径）と `cone_tan_f`（tan 半頂角）で**引数化**され、
/// **本影**縁は `(l2, tan f2)`・**半影**縁は `(l1, tan f1)` を渡して同一機構で解く（M9 残(3) 3a）。
///
/// 各反復で現在の縁推定 (ξ,η,ζ) から相対速度 rel=(x'−μ'(ζcos d−η sin d), y'−μ'ξ sin d) を作り、
/// それに直交する単位法線 n̂=(−rel_y, rel_x)/|rel| の `sign` 側へ ζ補正錐半径 `|cone_l−ζ·cone_tan_f|` だけ
/// 軸 (x,y) からオフセットし、地表へ射影して新しい (ξ,η,ζ) を得る。(ξ,η,ζ) が
/// [`LIMIT_FIXED_POINT_TOL`] 以内で安定（＝rel・半径とも自己整合）した縁点を返す。射影が
/// 地表を外す（`RootNotBracketed`）/ 相対速度ゼロ/ **[`LIMIT_FIXED_POINT_MAX_ITER`] 内に未収束**は
/// `Ok(None)`（lockstep スキップ）＝**未収束の近似を無警告で返さない**（誤差を隠さない・conventions §11）。
/// 他 `Err` は伝播。`vx`/`vy`/`mu_rate` は影軸運動 x'/y' と地球自転位相 μ'（基本面・時間単位は呼び側統一）。
///
/// `zeta0` は初期推定（中心軸 ζ）。半径（本影 |L2'|≈0.01 Re / 半影 |L1'|≈0.5 Re）でも写像は縮小で
/// **収束先は初期値に依らない**（良い初期値は反復数を減らすのみ）。よって `zeta0` を別値に変える/捨てる
/// 変異は結果を変えず等価（工程7 で生存列挙・許容）。
#[allow(clippy::too_many_arguments)]
pub(crate) fn solve_limit_edge(
    elements: &InstantaneousBesselianElements,
    cone_l: f64,
    cone_tan_f: f64,
    zeta0: f64,
    vx: f64,
    vy: f64,
    mu_rate: f64,
    sign: f64,
    ellipsoid: &Ellipsoid,
) -> Result<Option<GeoPoint>, EclipseError> {
    let (sin_d, cos_d) = elements.declination.0.sin_cos();
    // 縁推定（初期＝中心軸）。ξ,η は基本面オフセット点、ζ はその地表射影での影軸方向成分。
    let mut xi = elements.x;
    let mut eta = elements.y;
    let mut zeta = zeta0;
    let mut point: Option<GeoPoint> = None;
    let mut converged = false;
    for _ in 0..LIMIT_FIXED_POINT_MAX_ITER {
        // 影の地表に対する相対速度（ω×P の基本面 (ξ,η) 成分を影軸運動から差し引く）。
        let rel_x = vx - mu_rate * (zeta * cos_d - eta * sin_d);
        let rel_y = vy - mu_rate * xi * sin_d;
        let rel_speed = rel_x.hypot(rel_y);
        if rel_speed == 0.0 {
            return Ok(None); // 相対速度ゼロ（影が地表に静止）はスキップ。
        }
        // rel に直交する単位法線（路限界方向）。
        let nx = -rel_y / rel_speed;
        let ny = rel_x / rel_speed;
        // ζ補正錐半径（縁点自身の ζ で評価＝自己整合・本影 l2/半影 l1 は呼び側が cone_l で選ぶ）。
        let radius = (cone_l - zeta * cone_tan_f).abs();
        let xi_new = elements.x + sign * radius * nx;
        let eta_new = elements.y + sign * radius * ny;
        let (p, zeta_new) = match surface_point_for_fundamental(
            xi_new,
            eta_new,
            elements.declination,
            elements.mu,
            ellipsoid,
        ) {
            Ok(v) => v,
            Err(EclipseError::Solver(SolverError::RootNotBracketed)) => return Ok(None),
            Err(e) => return Err(e),
        };
        // 収束判定（過収束機構）: 反復は機械精度まで収束し ξ,η,ζ は同率で <TOL に達するため、境界
        // `< → <=` や `&& → ||` は真解に影響しない等価変異（`< → >`＝発散は caught）。mutation.yml で
        // `with <= in.*solve_limit_edge` / `with \|\| in.*solve_limit_edge` を除外（`||` は正規表現の空
        // alternation で全マッチするため必ず `\|\|` とエスケープ。docs/reviews/mutation-limit-line.md）。
        let step_converged = (xi_new - xi).abs() < LIMIT_FIXED_POINT_TOL
            && (eta_new - eta).abs() < LIMIT_FIXED_POINT_TOL
            && (zeta_new - zeta).abs() < LIMIT_FIXED_POINT_TOL;
        xi = xi_new;
        eta = eta_new;
        zeta = zeta_new;
        point = Some(p);
        if step_converged {
            converged = true;
            break;
        }
    }
    // 未収束（rel 悪条件等・実太陽日食では rel≈0.2 Re/h で到達しない）は近似を返さずスキップ。
    Ok(if converged { point } else { None })
}

/// 影錐縁（ζ=0・半径 `cone_l`）と WGS84 terminator（ζ=0 の地表＝太陽が地平）の交点を測地座標で返す
/// （rise/set limb 点・部分食域 §11.2/11.3・M9 残(3) 3b）。
///
/// 基本面で円 `(ξ−x)²+(η−y)²=cone_l²`（中心は影軸 (x,y)）と terminator 楕円 `ξ²+k·η²=1`
/// （`k=sin²d+cos²d/(1−f)²`・f は扁平）の交点を求める。楕円を θ で `ξ=cos θ, η=sin θ/√k` と媒介し、円残差
/// `r(θ)=(ξ−x)²+(η−y)²−cone_l²` の符号反転を [0,2π) の粗走査（[`TERMINATOR_SCAN_STEPS`]）でブラケットして
/// Brent で全根を求め、各根 (ξ,η) を [`fundamental_to_geodetic`]（ζ=0）で測地座標化し **θ 昇順**で返す。
/// 交点数は 0/1/2/(最大4)。錐が terminator に届かない/完全内包なら空 `Vec`。**接する端（P1/P4 の 1 点接触）は
/// 符号反転が無く拾わない**（交点数 0↔2 の遷移点＝外周端の扱いは後続 (3c) の責務）。
///
/// `cone_l` は ζ=0 の錐半径（半影は `elements.l1`）。terminator を WGS84 楕円で厳密化しているため、返点は
/// WGS84 前方射影で **ζ≈0（terminator 上）・軸からの面内距離≈cone_l（錐縁上）** に往復一致する（球近似 k=1 の
/// ~Re·f 残差を排除＝中心線/限界線と同じ WGS84-exact オラクルで検証可能）。`fundamental_to_geodetic` の緯度
/// 検証等以外の `Err` は伝播。
pub(crate) fn cone_terminator_intersections(
    elements: &InstantaneousBesselianElements,
    cone_l: f64,
    ellipsoid: &Ellipsoid,
) -> Result<Vec<GeoPoint>, EclipseError> {
    let (sin_d, cos_d) = elements.declination.0.sin_cos();
    let omf = 1.0 - ellipsoid.f; // b/a
                                 // terminator 楕円係数 k（球 f=0 で 1＝単位円・扁平で k≥1＝η 方向に縮む）。
    let k = sin_d * sin_d + cos_d * cos_d / (omf * omf);
    let sqrt_k = k.sqrt();
    let x = elements.x;
    let y = elements.y;

    // 媒介点 θ（terminator 楕円上 ξ=cos θ, η=sin θ/√k）に対する錐縁の円残差。
    let residual = |theta: f64| -> f64 {
        let (s, c) = theta.sin_cos();
        let xi = c;
        let eta = s / sqrt_k;
        (xi - x) * (xi - x) + (eta - y) * (eta - y) - cone_l * cone_l
    };

    // 符号反転の粗走査＋Brent（resolution 機構＝[`scan_periodic_sign_change_roots`]）で全根 θ を得て、
    // 各根を測地座標化（ζ=0）。点が無ければ空 Vec（錐が terminator に届かない/内包）。
    let mut points = Vec::new();
    for root in scan_periodic_sign_change_roots(&residual)? {
        let (s, c) = root.sin_cos();
        points.push(fundamental_to_geodetic(
            c,
            s / sqrt_k,
            0.0,
            elements.declination,
            elements.mu,
            ellipsoid,
        )?);
    }
    Ok(points)
}

/// 周期残差 `r(θ)`（θ∈[0,2π)）の符号反転を粗走査でブラケットし Brent で全根 θ を昇順に返す。
///
/// **resolution 機構**（[`descending_sign_change_bracket`] / `solve_zero_in_window` / `scan_point_count` と
/// 同カテゴリ）: 区間内で根が分離していれば、刻み `TAU/STEPS`・符号積 `prev_r·r`・境界 `< 0` の算術は
/// Brent の収束先（真根）に影響しない解像度要素。符号積の `* → /`（積と商は符号判定が同値）・`< → <=`
/// （厳密ゼロ＝グリッドが根に一致する測度ゼロのみ差）は真根を変えない等価変異ゆえ `mutation.yml` で
/// `--exclude-re 'in scan_periodic_sign_change_roots'` 除外（docs/reviews/mutation-rise-set.md）。閉ループ端
/// （θ=0 と 2π は同残差）で根を二重計上しないのは strict `< 0`（厳密 0 は測度ゼロ）による。
fn scan_periodic_sign_change_roots<F: Fn(f64) -> f64>(
    residual: &F,
) -> Result<Vec<f64>, EclipseError> {
    let step = TAU / TERMINATOR_SCAN_STEPS as f64;
    let mut roots = Vec::new();
    let mut prev_theta = 0.0_f64;
    let mut prev_r = residual(prev_theta);
    for i in 1..=TERMINATOR_SCAN_STEPS {
        let theta = step * i as f64; // i=STEPS で 2π（=0 と同点・閉ループの終端）
        let r = residual(theta);
        if prev_r * r < 0.0 {
            roots.push(brent_root(
                residual,
                prev_theta,
                theta,
                TERMINATOR_ROOT_TOL,
                TERMINATOR_ROOT_MAX_ITER,
            )?);
        }
        prev_theta = theta;
        prev_r = r;
    }
    Ok(roots)
}

/// 2 つの地表点間の概算大圏距離 \[km\]（haversine・[`EARTH_MEAN_RADIUS_KM`]）。帯幅（北縁-南縁距離）に使う。
pub(crate) fn great_circle_distance_km(a: &GeoPoint, b: &GeoPoint) -> f64 {
    let lat1 = a.lat.radians().0;
    let lat2 = b.lat.radians().0;
    let dlat = lat2 - lat1;
    let dlon = b.lon.radians().0 - a.lon.radians().0;
    let h = (dlat / 2.0).sin().powi(2) + lat1.cos() * lat2.cos() * (dlon / 2.0).sin().powi(2);
    2.0 * EARTH_MEAN_RADIUS_KM * h.sqrt().asin()
}

/// `from` から `to` への大圏**初期方位**（bearing, rad・`[0, 2π)`・北=0・東=π/2・南=π・西=3π/2）。
///
/// 部分食域の外周組立（§11.4・(3c-ii)）で、全境界点を最大食点まわりの方位で整列するのに使う。
/// 標準式 `θ = atan2(sinΔλ·cosφ2, cosφ1·sinφ2 − sinφ1·cosφ2·cosΔλ)`（Δλ=λ2−λ1）を `[0,2π)` へ正規化する。
/// 経度ラップに非依存（Δλ の三角関数経由）。
pub(crate) fn initial_bearing(from: &GeoPoint, to: &GeoPoint) -> f64 {
    let phi1 = from.lat.radians().0;
    let phi2 = to.lat.radians().0;
    let dlon = to.lon.radians().0 - from.lon.radians().0;
    let y = dlon.sin() * phi2.cos();
    let x = phi1.cos() * phi2.sin() - phi1.sin() * phi2.cos() * dlon.cos();
    y.atan2(x).rem_euclid(TAU)
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

    // ---- 要件8: 二次方程式オラクルによる ζ₊ と測地座標の絶対値ピン（strict, 独立算術） ----
    //
    // 戦略: 残差 r(ζ)=(px²+ξ²)+pz²/omf²−1（px=ζcosd−η sind, pz=ζsind+η cosd）は ζ の二次式
    //   A·ζ²+B·ζ+C で、A=cos²d+sin²d/omf², B=2η·sind·cosd·(1/omf²−1),
    //   C=(η²sin²d+ξ²)+η²cos²d/omf²−1。`surface_point_for_fundamental` の戻り ζ は太陽側（最大）根
    //   ζ₊=(−B+√(B²−4AC))/(2A)。テスト内でこの**閉形式**から ζ₊ と (φ,λ) を独立に組み、
    //   関数の戻りと突き合わせる。前方往復オラクル（要件1-7）とは別系統で、px/pz/ρcos²/pz²/omf² の
    //   各積・和の演算子が 1 つでも壊れれば ζ₊ または φ,λ がずれて落ちる。

    /// ξ=0.30, η=0.20, d=0.40, μ=1.10, WGS84。全項（cosd・sind の積、1/omf² 重み、ξ² 項、
    /// η の交差項 B）が非自明に効く構成。閉形式の ζ₊ と、それを使った px/pz/ρcos/ρsin/H から
    /// 復元した φ,λ（度）を独立に組み、関数戻りと ζ≈1e-9・lat/lon≈1e-7 deg で一致させる。
    #[test]
    fn quadratic_oracle_pins_root_and_geodetic() {
        let xi = 0.30_f64;
        let eta = 0.20_f64;
        let d = 0.40_f64;
        let mu = 1.10_f64;

        let (sin_d, cos_d) = d.sin_cos();
        let omf = 1.0 - WGS84.f;
        let inv_omf2 = 1.0 / (omf * omf);

        // 二次係数（独立算術。contract の A,B,C をそのまま組む）。
        let a = cos_d * cos_d + sin_d * sin_d * inv_omf2;
        let b = 2.0 * eta * sin_d * cos_d * (inv_omf2 - 1.0);
        let c =
            (eta * eta * sin_d * sin_d + xi * xi) + (eta * eta * cos_d * cos_d) * inv_omf2 - 1.0;

        let disc = b * b - 4.0 * a * c;
        assert!(
            disc > 0.0,
            "選定 (ξ,η,d) は実根を持つ必要がある（軸が地表に当たる）: disc={disc}"
        );
        // 太陽側＝大きい方の根 ζ₊。A>0 なので (−B+√disc)/(2A) が大きい根。
        let zeta_plus = (-b + disc.sqrt()) / (2.0 * a);
        assert!(zeta_plus > 0.0, "ζ₊={zeta_plus} は太陽側 (>0) のはず");

        // ζ₊ で残差が 0（点が確かに地表上）であることを独立に確認。
        let px_chk = zeta_plus * cos_d - eta * sin_d;
        let pz_chk = zeta_plus * sin_d + eta * cos_d;
        let residual = (px_chk * px_chk + xi * xi) + pz_chk * pz_chk * inv_omf2 - 1.0;
        assert!(
            residual.abs() < 1e-12,
            "ζ₊ で楕円体残差≈0 のはず: residual={residual}"
        );

        // 独立に測地座標を組む（contract の復元式）:
        //   px=ζcosd−η sind, py=ξ, pz=ζsind+η cosd
        //   ρcos=√(px²+py²), ρsin=pz, H=atan2(py,px)
        //   λ_east=normalize_[-π,π)(H−μ), φ=atan2(ρsin, ρcos·omf²)
        let px = zeta_plus * cos_d - eta * sin_d;
        let py = xi;
        let pz = zeta_plus * sin_d + eta * cos_d;
        let rho_cos = (px * px + py * py).sqrt();
        let rho_sin = pz;
        let h = py.atan2(px);
        let expected_lam = Radians::new(h - mu).normalized_signed().0;
        let expected_phi = rho_sin.atan2(rho_cos * omf * omf);

        // 関数を呼び、ζ と (φ,λ) を絶対値で突き合わせる。
        let (point, zeta) = surface_point_for_fundamental(xi, eta, Radians(d), Radians(mu), &WGS84)
            .expect("選定構成は実根を持つので Ok");

        // ζ: 二次閉形式の太陽側根に一致（A·ζ²+B·ζ+C 各係数の演算子を縛る）。
        assert!(
            close(zeta, zeta_plus, 1e-9),
            "ζ={zeta} expected ζ₊={zeta_plus}（二次オラクル）"
        );
        // φ: atan2(ρsin, ρcos·omf²) の独立復元に一致（pz・ρcos・omf² の積を縛る）。
        let got_phi = point.lat.radians().0;
        assert!(
            close(got_phi, expected_phi, 1e-7),
            "φ={got_phi} expected {expected_phi}（測地緯度復元）"
        );
        // λ_east: normalize(H−μ) の独立復元に一致（H=atan2(py,px) と −μ・正規化を縛る）。
        let got_lam = point.lon.radians().0;
        assert!(
            close(got_lam, expected_lam, 1e-7),
            "λ={got_lam} expected {expected_lam}（東経復元）"
        );
        // 度でも縛る（degrees アクセサ経路）。
        assert!(
            close(point.lat.degrees().0, expected_phi.to_degrees(), 1e-7),
            "φ_deg mismatch"
        );
        assert!(
            close(point.lon.degrees().0, expected_lam.to_degrees(), 1e-7),
            "λ_deg mismatch"
        );
    }

    /// 軸が地表を明確に外す構成（ξ=2.0, η=2.0 ⇒ B²−4AC<0, 実根なし）⇒
    /// `Err(Solver(RootNotBracketed))`。二次オラクルで判別式が負であることも独立確認。
    #[test]
    fn quadratic_oracle_missing_line_returns_root_not_bracketed() {
        let xi = 2.0_f64;
        let eta = 2.0_f64;
        let d = 0.40_f64;
        let mu = 1.10_f64;

        let (sin_d, cos_d) = d.sin_cos();
        let omf = 1.0 - WGS84.f;
        let inv_omf2 = 1.0 / (omf * omf);
        let a = cos_d * cos_d + sin_d * sin_d * inv_omf2;
        let b = 2.0 * eta * sin_d * cos_d * (inv_omf2 - 1.0);
        let c =
            (eta * eta * sin_d * sin_d + xi * xi) + (eta * eta * cos_d * cos_d) * inv_omf2 - 1.0;
        let disc = b * b - 4.0 * a * c;
        assert!(
            disc < 0.0,
            "選定 (ξ,η) は実根を持たない（軸が地表を外す）はず: disc={disc}"
        );

        let r = surface_point_for_fundamental(xi, eta, Radians(d), Radians(mu), &WGS84);
        assert!(
            matches!(r, Err(EclipseError::Solver(SolverError::RootNotBracketed))),
            "expected Err(Solver(RootNotBracketed)), got {r:?}"
        );
    }

    // ====================================================================
    // solve_limit_edge（錐半径引数化）FAST 単体テスト（M9 残(3) 3a）
    //
    // 戦略（追認回避）: 各縁点 P を**検証済み前方射影** project_observer_to_fundamental へ通して
    // 自身の (ξ,η,ζ) を独立復元し（path_limits.rs `assert_exact_limit_conditions` と同パターン）、
    // 次の 2 条件を絶対値で表明する:
    //   条件1（錐exact・自己整合ζ）: hypot(ξ−x, η−y) = |cone_l − ζ·cone_tan_f|（ζ は点自身の値）。
    //   条件2（包絡）: (ξ−x)·rel_vx + (η−y)·rel_vy ≈ 0,
    //     rel_vx = vx − μ'·(ζ·cos d − η·sin d), rel_vy = vy − μ'·ξ·sin d。
    // 半径以外（軸位置・rel・μ'・sign）は本影/半影で同形なので、半影版は半径を |l1−ζ·tan f1| に
    // するだけ。被テスト関数の戻りを期待値生成に流用しない（前方射影は独立に検証済み）。
    // ====================================================================

    /// 縁点 P を**検証済み前方射影**へ通し、自身の (ξ,η,ζ) を `ObserverFundamental` で返す独立オラクル。
    /// `assert_forward_roundtrip` は ζ だけ返すため、2 条件用に ξ,η,ζ を全部使えるヘルパを別途用意する。
    fn forward_of(
        p: &GeoPoint,
        e: &InstantaneousBesselianElements,
    ) -> crate::projection::ObserverFundamental {
        let phi = p.lat.radians().0;
        let lam = p.lon.radians().0;
        let obs = observer_geocentric(&WGS84, phi, 0.0);
        project_observer_to_fundamental(&obs, Radians::new(lam), e)
    }

    /// 縁点 P が錐半径 `(cone_l, cone_tan_f)` について 2 条件（自己整合ζの錐exact＋相対速度包絡⊥）を
    /// 満たすことを表明する共有チェック。半径は P 自身の ζ で評価する（自己整合）。
    /// `cone_tol`/`dot_tol` は呼び側が与える（合成＝機械精度）。
    #[allow(clippy::too_many_arguments)]
    fn assert_edge_conditions(
        p: &GeoPoint,
        e: &InstantaneousBesselianElements,
        cone_l: f64,
        cone_tan_f: f64,
        vx: f64,
        vy: f64,
        mu_rate: f64,
        cone_tol: f64,
        dot_tol: f64,
    ) {
        let (sin_d, cos_d) = e.declination.0.sin_cos();
        let of = forward_of(p, e);
        let off_x = of.xi - e.x;
        let off_y = of.eta - e.y;
        // 条件1: 面内距離 = |cone_l − ζ·cone_tan_f|（自己整合ζ）。
        let in_plane = off_x.hypot(off_y);
        let radius = (cone_l - of.zeta * cone_tan_f).abs();
        assert!(
            (in_plane - radius).abs() < cone_tol,
            "面内距離 {in_plane} = |cone_l − ζ·cone_tan_f| {radius}（自己整合ζ={}）でない",
            of.zeta
        );
        // 条件2: offset ⊥ rel 速度（μ' 項込み）。
        let rel_vx = vx - mu_rate * (of.zeta * cos_d - of.eta * sin_d);
        let rel_vy = vy - mu_rate * of.xi * sin_d;
        let dot = off_x * rel_vx + off_y * rel_vy;
        assert!(
            dot.abs() < dot_tol,
            "offset·rel = {dot}（≈0 でない＝包絡条件違反）"
        );
        // 面内距離の独立サニティ（半径が正＝縁が軸からオフセットしている）。
        assert!(
            in_plane > 0.0,
            "縁は軸からオフセットしているはず（面内距離>0）"
        );
    }

    /// 軸が地表に当たる合成中心食の代表構成（gamma≪1・μ'≠0 で包絡条件に分離力）。
    /// vx,vy,mu_rate は合成の既知値（多項式微分不要）。本影も半影も縁が地表に当たる小さい x,y。
    fn edge_config() -> (InstantaneousBesselianElements, f64, f64, f64, f64) {
        // x,y は影軸を地表中央付近に置く（gamma=√(0.1²+0.05²)≈0.112≪1）。
        // d は現実的な太陽見かけ赤緯, μ は任意。l1/l2/tan_f は elems() の既定。
        let e = elems(0.10, 0.05, 0.20, 1.10);
        // 影軸運動 vx,vy（Re/h 級）と地球自転位相 μ'（≠0 で包絡に効く）。
        let (vx, vy, mu_rate) = (0.45_f64, 0.06_f64, 0.26_f64);
        let zeta0 = 0.9_f64; // 初期推定（縮小写像ゆえ収束先に無関係）。
        (e, vx, vy, mu_rate, zeta0)
    }

    /// 半影縁の 2 条件: (l1, tan_f1) で呼ぶと、返点が条件1=面内距離 |l1−ζ·tan f1|・条件2=offset⊥rel を
    /// 機械精度で満たす。μ'≠0 で包絡条件に分離力。sign=±1 の両側を縛る。
    ///
    /// 殺す変異: cone_l↔cone_tan_f 取り違え（半径式が ζ·l1 等になり面内距離が大きく外れる）・
    ///   半径計算 `cone_l − ζ·cone_tan_f` の演算子（`−`→`+`/`*`/`/`）・rel の μ' 項脱落（条件2 が dot≠0）・
    ///   cos d↔sin d・ξ↔η・sign 反映漏れ。
    #[test]
    fn penumbral_edge_satisfies_cone_and_envelope_conditions() {
        let (e, vx, vy, mu_rate, zeta0) = edge_config();
        // μ'≠0 を独立確認（包絡条件の分離力の前提）。
        assert!(mu_rate.abs() > 1e-6, "μ'≠0 構成（rel 速度に μ' が効く）");
        for sign in [1.0_f64, -1.0_f64] {
            let p = solve_limit_edge(&e, e.l1, e.tan_f1, zeta0, vx, vy, mu_rate, sign, &WGS84)
                .expect("半影縁の解は Ok")
                .unwrap_or_else(|| panic!("sign={sign}: 半影縁点が存在するはず（gamma≪1）"));
            // 半影半径で 2 条件（自己整合ζ）を機械精度で。
            assert_edge_conditions(&p, &e, e.l1, e.tan_f1, vx, vy, mu_rate, 1e-7, 1e-9);
        }
    }

    /// 本影回帰: (l2, tan_f2) で呼ぶと、シグネチャ変更後も従来同様に 2 条件（半径=|l2−ζ·tan f2|）を
    /// 満たす（本影が壊れていない独立確認）。sign=±1 の両側。
    ///
    /// 殺す変異: 引数化で本影半径が l1 等に化ける・半径式の演算子・rel の μ' 項脱落・sign 反映漏れ。
    #[test]
    fn umbral_edge_regression_satisfies_cone_and_envelope_conditions() {
        let (e, vx, vy, mu_rate, zeta0) = edge_config();
        for sign in [1.0_f64, -1.0_f64] {
            let p = solve_limit_edge(&e, e.l2, e.tan_f2, zeta0, vx, vy, mu_rate, sign, &WGS84)
                .expect("本影縁の解は Ok")
                .unwrap_or_else(|| panic!("sign={sign}: 本影縁点が存在するはず（gamma≪1）"));
            // 本影半径で 2 条件（自己整合ζ）を機械精度で。
            assert_edge_conditions(&p, &e, e.l2, e.tan_f2, vx, vy, mu_rate, 1e-7, 1e-9);
        }
    }

    /// 半影 vs 本影の半径差: 同一 elements・同一 sign で半影縁の軸からの面内距離 > 本影縁のそれ。
    /// 半径 l1(~0.5)≫|l2|(~0.01) なので半影縁は本影縁より遥かに軸から遠い（l1↔l2 取り違えを撃つ）。
    ///
    /// 殺す変異: cone_l に l2 を使い続ける（半影呼び出しでも本影半径＝差が消える）・cone_l↔cone_tan_f
    ///   取り違え（両者とも壊れて差が偶然一致しうるが、上の 2 条件テストが別途撃つ）。
    #[test]
    fn penumbral_edge_is_farther_from_axis_than_umbral() {
        let (e, vx, vy, mu_rate, zeta0) = edge_config();
        let sign = 1.0_f64;
        let p_umbra = solve_limit_edge(&e, e.l2, e.tan_f2, zeta0, vx, vy, mu_rate, sign, &WGS84)
            .expect("本影縁 Ok")
            .expect("本影縁点が存在");
        let p_penumbra = solve_limit_edge(&e, e.l1, e.tan_f1, zeta0, vx, vy, mu_rate, sign, &WGS84)
            .expect("半影縁 Ok")
            .expect("半影縁点が存在");
        // 各縁点を前方射影し、軸 (x,y) からの基本面内距離を独立に測る。
        let of_u = forward_of(&p_umbra, &e);
        let of_p = forward_of(&p_penumbra, &e);
        let off_umbra = (of_u.xi - e.x).hypot(of_u.eta - e.y);
        let off_penumbra = (of_p.xi - e.x).hypot(of_p.eta - e.y);
        // |l2|≈0.01 vs l1≈0.5 ⇒ 半影は本影より桁違いに軸から遠い。
        assert!(
            off_penumbra > off_umbra,
            "半影縁の面内距離 {off_penumbra} > 本影縁 {off_umbra}（l1≫|l2|）"
        );
        // 桁差を独立に確認（取り違えで両者が近接する変異を撃つ）: 半影は本影の 10 倍以上遠い。
        assert!(
            off_penumbra > 10.0 * off_umbra,
            "半影縁は本影縁の 10 倍超遠いはず: penumbra {off_penumbra}, umbra {off_umbra}"
        );
    }

    /// sign ±1 が反対側の縁を与える（南北別側）: 同一 elements・同一半径で、sign=+1 と sign=−1 の
    /// 縁点は軸を挟んで反対側に出る。前方射影した offset ベクトルが互いに逆向き（内積<0）。
    ///
    /// 殺す変異: sign を無視して常に同じ側へ出す（offset が同一＝内積>0）・sign の符号反転を片側に倒す。
    #[test]
    fn sign_gives_opposite_sides_of_axis() {
        let (e, vx, vy, mu_rate, zeta0) = edge_config();
        // 半影半径で両側を取る（軸から遠く、対称性が見やすい）。
        let p_pos = solve_limit_edge(&e, e.l1, e.tan_f1, zeta0, vx, vy, mu_rate, 1.0, &WGS84)
            .expect("sign=+1 Ok")
            .expect("sign=+1 縁点が存在");
        let p_neg = solve_limit_edge(&e, e.l1, e.tan_f1, zeta0, vx, vy, mu_rate, -1.0, &WGS84)
            .expect("sign=−1 Ok")
            .expect("sign=−1 縁点が存在");
        let of_pos = forward_of(&p_pos, &e);
        let of_neg = forward_of(&p_neg, &e);
        let (ox_p, oy_p) = (of_pos.xi - e.x, of_pos.eta - e.y);
        let (ox_n, oy_n) = (of_neg.xi - e.x, of_neg.eta - e.y);
        // 反対側＝offset ベクトルが逆向き（内積<0）。
        let dot = ox_p * ox_n + oy_p * oy_n;
        assert!(
            dot < 0.0,
            "sign=+1 と sign=−1 の offset は反対側（内積<0）, got {dot}"
        );
    }

    /// 縁が地表を外す合成 ⇒ `Ok(None)`。軸を地表の縁ギリギリ（gamma<1）に置きつつ、半影半径
    /// l1≈0.5 で外側へ押すと縁点の射影が地表を外す（`RootNotBracketed`）→ `Ok(None)`。
    /// `Err` でも `Some` でもなく `Ok(None)`（未収束/外れは近似を返さずスキップ）。
    ///
    /// 殺す変異: `RootNotBracketed` を `Err` 伝播する・外れでも近似 `Some` を返す・None 分岐を潰す。
    #[test]
    fn edge_off_surface_returns_ok_none() {
        // x を限界近く（0.92）に置き gamma=0.92<1（軸は地表に当たる）。半影半径 l1≈0.5 で
        // sign=+1（外向き）に押すと縁点が gamma>1 域へ出て地表を外す。
        let e = elems(0.92, 0.0, 0.20, 1.10);
        let gamma = (e.x * e.x + e.y * e.y).sqrt();
        assert!(gamma < 1.0, "軸自身は地表に当たる前提 gamma={gamma}<1");
        let (vx, vy, mu_rate) = (0.45_f64, 0.06_f64, 0.26_f64);
        let r = solve_limit_edge(&e, e.l1, e.tan_f1, 0.4, vx, vy, mu_rate, 1.0, &WGS84);
        assert!(
            matches!(r, Ok(None)),
            "縁が地表を外す合成では Ok(None) のはず, got {r:?}"
        );
    }

    /// 相対速度ゼロ（vx=vy=0 かつ μ'=0）⇒ rel=0 ⇒ `Ok(None)`（影が地表に静止・スキップ）。
    ///
    /// 殺す変異: rel_speed==0 ガードを外す（0 除算で NaN を返す/panic）・Ok(None) を別値にする。
    #[test]
    fn zero_relative_speed_returns_ok_none() {
        let e = elems(0.10, 0.05, 0.20, 1.10);
        // vx=vy=0, μ'=0 ⇒ rel=(0,0)（最初の反復で rel_speed==0）。
        let r = solve_limit_edge(&e, e.l1, e.tan_f1, 0.9, 0.0, 0.0, 0.0, 1.0, &WGS84);
        assert!(
            matches!(r, Ok(None)),
            "相対速度ゼロでは Ok(None) のはず, got {r:?}"
        );
    }

    // ====================================================================
    // cone_terminator_intersections（錐縁 ∩ terminator 楕円）FAST 単体テスト（M9 残(3) 3b）
    //
    // 確定仕様（docs/algorithms/11-path-partial-domain.md §11.2/§11.3, WGS84 厳密）:
    //   terminator 楕円（基本面 ζ=0）: ξ² + k·η² = 1,  k = sin²d + cos²d/(1−f)²。
    //   錐縁の円（ζ=0）:               (ξ−x)² + (η−y)² = cone_l²（中心は影軸 (x,y)=(elements.x,y)）。
    //   交点 (ξ,η) を fundamental_to_geodetic(ξ,η,0,d,μ) で測地座標化して Vec<GeoPoint> で返す。
    //   交点数 0/1/2（…最大4）。届かない/完全内包なら空 Vec。
    //
    // 戦略（追認回避・独立オラクル）: 返る各点 P を**検証済み前方射影** project_observer_to_fundamental
    // （forward_of）へ通して自身の (ξ,η,ζ) を独立復元し、契約を絶対値で表明する:
    //   契約1（terminator 上）: ζ ≈ 0（日の出入りの最中）。球近似 k=1 なら扁平 d≠0 で ζ≠0 になり落ちる
    //     ＝WGS84 厳密 terminator を担保。
    //   契約2（錐縁上）: hypot(ξ−x, η−y) ≈ cone_l（面内距離が錐半径）。
    //   契約3: 妥当な緯度経度（GeodeticLatitude 検証通過＝Ok で返ること自体が担保）。
    // d=0（k=1）では二円交点の閉形式（根軸 ∩ 単位円）で (ξ,η) を手計算し、返点と突き合わせる。
    // ====================================================================

    /// 返点 P を前方射影し、契約1（ζ≈0）・契約2（面内距離≈cone_l）を機械精度で表明する共有チェック。
    /// 期待値生成に被テスト関数の戻りを流用しない（前方射影は独立に検証済み）。
    fn assert_terminator_intersection(
        p: &GeoPoint,
        e: &InstantaneousBesselianElements,
        cone_l: f64,
        zeta_tol: f64,
        cone_tol: f64,
    ) {
        let of = forward_of(p, e);
        // 契約1: terminator 上（ζ≈0・太陽が地平）。WGS84 厳密 ⇒ 機械精度で 0。
        assert!(
            of.zeta.abs() < zeta_tol,
            "ζ={} must be ≈0 (on terminator, WGS84-exact)",
            of.zeta
        );
        // 契約2: 錐縁上（軸からの基本面内距離 = cone_l）。
        let in_plane = (of.xi - e.x).hypot(of.eta - e.y);
        assert!(
            (in_plane - cone_l).abs() < cone_tol,
            "面内距離 {in_plane} = cone_l {cone_l} でない（錐縁上でない）"
        );
    }

    /// 2 交点・契約検証: 軸を端（x≈0.8）に置いた合成（d≠0・μ≠0）で 2 点返り、各点が前方射影往復で
    /// ζ≈0（契約1・WGS84 厳密 terminator）・面内距離≈cone_l（契約2）を機械精度で満たす。
    /// 円が terminator 楕円を跨ぐよう cone_l を適度に取る（軸が端なので円縁が楕円外周をまたぐ）。
    ///
    /// 殺す変異: 中心 (x,y) の取り違え（軸でない点を中心にすると面内距離≠cone_l）・cone_l の係数・
    ///   ζ=0 引数（fundamental_to_geodetic に 0 以外を渡すと ζ≠0 で契約1 違反）。
    #[test]
    fn two_intersections_satisfy_contracts() {
        // 軸を端に（x=0.80, y=0.10 ⇒ 円中心は楕円縁付近）。d≠0・μ≠0 で扁平と回転が効く。
        let e = elems(0.80, 0.10, 0.20, 1.10);
        let cone_l = 0.40_f64; // 円が terminator 楕円を跨ぐ適度な半径。
        let pts = cone_terminator_intersections(&e, cone_l, &WGS84)
            .expect("錐縁∩terminator は Ok のはず");
        assert_eq!(
            pts.len(),
            2,
            "軸を端＋適度な cone_l で 2 交点のはず, got {pts:?}"
        );
        for p in &pts {
            assert_terminator_intersection(p, &e, cone_l, 1e-7, 1e-7);
        }
    }

    /// k=1（d=π/2）で二円交点の閉形式一致: WGS84 terminator は基本面で楕円 ξ²+k·η²=1,
    /// k=sin²d+cos²d/(1−f)²。二円（円∩**単位円**）の閉形式が成立するのは k=1 のときのみ
    /// ＝cos d=0 ⇔ d=±π/2。そこで d=π/2 を採り真に単位円 ξ²+η²=1 にして既存の閉形式を有効化する。
    /// 根軸 x·ξ+y·η=(x²+y²+1−cone_l²)/2 ∩ 単位円の閉形式で 2 点を独立に組み、前方射影した
    /// 返点の (ξ,η) と集合一致（順序非依存）で突き合わせる。ζ≈0（面内距離＝cone_l）も確認。
    ///
    /// 殺す変異: 円側の中心 (x,y)・cone_l を撃つ。閉形式は厳密ゆえ機械精度で縛れる。
    #[test]
    fn d_zero_matches_two_circle_closed_form() {
        use std::f64::consts::FRAC_PI_2;
        let x = 0.30_f64;
        let y = 0.40_f64;
        let cone_l = 0.80_f64;
        // d=π/2 ⇒ k=sin²d+cos²d/(1−f)²=1 ⇒ terminator は単位円 ξ²+η²=1。円 (ξ−x)²+(η−y)²=cone_l²。
        // 根軸: x·ξ + y·η = c,  c = (x²+y²+1−cone_l²)/2。
        let c = (x * x + y * y + 1.0 - cone_l * cone_l) / 2.0;
        // 根軸 ∩ 単位円の閉形式（d2 = x²+y²>0）:
        //   ξ = (x·c ± y·√(d2−c²))/d2,  η = (y·c ∓ x·√(d2−c²))/d2。
        let d2 = x * x + y * y;
        let disc = d2 - c * c;
        assert!(disc > 0.0, "選定は 2 交点を持つ（判別式>0）: disc={disc}");
        let s = disc.sqrt();
        let sol = |sgn: f64| -> (f64, f64) {
            let xi = (x * c + sgn * y * s) / d2;
            let eta = (y * c - sgn * x * s) / d2;
            (xi, eta)
        };
        let expected = [sol(1.0), sol(-1.0)];
        // 閉形式が確かに両曲線上にあることを独立サニティ（単位円上・円上）。
        for (xi, eta) in expected {
            assert!(
                (xi * xi + eta * eta - 1.0).abs() < 1e-12,
                "閉形式は単位円上のはず"
            );
            assert!(
                ((xi - x).powi(2) + (eta - y).powi(2) - cone_l * cone_l).abs() < 1e-12,
                "閉形式は錐縁円上のはず"
            );
        }

        let e = elems(x, y, FRAC_PI_2, 0.7);
        let pts = cone_terminator_intersections(&e, cone_l, &WGS84).expect("Ok のはず");
        assert_eq!(pts.len(), 2, "k=1（d=π/2）二円交点で 2 点, got {pts:?}");
        // 順序非依存で集合一致: 各返点を前方射影し (ξ,η) を取り、期待 2 点のどちらかに一致。
        for p in &pts {
            let of = forward_of(p, &e);
            assert!(of.zeta.abs() < 1e-7, "ζ={} ≈0 (terminator)", of.zeta);
            let matched = expected
                .iter()
                .any(|(exi, eeta)| close(of.xi, *exi, 1e-7) && close(of.eta, *eeta, 1e-7));
            assert!(
                matched,
                "返点 (ξ={}, η={}) が閉形式 2 点 {expected:?} のどちらにも一致しない",
                of.xi, of.eta
            );
        }
        // 逆向きも: 期待 2 点それぞれに対応する返点が存在（全単射＝過不足なし）。
        for (exi, eeta) in expected {
            let found = pts.iter().any(|p| {
                let of = forward_of(p, &e);
                close(of.xi, exi, 1e-7) && close(of.eta, eeta, 1e-7)
            });
            assert!(found, "閉形式点 (ξ={exi}, η={eeta}) に対応する返点が無い");
        }
    }

    /// 空 Vec: 錐が terminator に届かない合成（軸 (x,y)≈0・cone_l 小 ⇒ 円が原点近傍で楕円
    /// |ξ²+kη²=1| に届かない）で交点 0。
    ///
    /// 殺す変異: 「届かない＝空」判定を潰して偽の交点を返す・空でなく Err にする。
    #[test]
    fn cone_not_reaching_terminator_returns_empty() {
        // 軸を原点付近、半径も小さく（円は原点近傍に完全に収まり、|ξ²+kη²|=1 楕円に届かない）。
        let e = elems(0.05, 0.0, 0.20, 1.10);
        let cone_l = 0.10_f64;
        let pts = cone_terminator_intersections(&e, cone_l, &WGS84)
            .expect("届かない合成でも Ok（空 Vec）のはず");
        assert!(
            pts.is_empty(),
            "錐が terminator に届かないので空のはず, got {pts:?}"
        );
    }

    /// cone_l 取り違え（l1 vs l2 相当）: cone_l を変えると交点位置/面内距離が連動する（半径の load-bearing 化）。
    /// 同一 elements で cone_l を 2 値に変え、各々で返点の面内距離が**渡した cone_l に一致**することを縛る。
    /// 半径を無視して固定値を使う変異は、片方で面内距離≠cone_l になり落ちる。
    ///
    /// 殺す変異: cone_l を定数で置換・l1/l2 取り違え・半径計算の脱落。
    #[test]
    fn cone_radius_is_load_bearing() {
        let e = elems(0.80, 0.10, 0.20, 1.10);
        for cone_l in [0.35_f64, 0.45_f64] {
            let pts = cone_terminator_intersections(&e, cone_l, &WGS84).expect("Ok");
            assert_eq!(pts.len(), 2, "cone_l={cone_l} で 2 交点のはず, got {pts:?}");
            for p in &pts {
                let of = forward_of(p, &e);
                let in_plane = (of.xi - e.x).hypot(of.eta - e.y);
                assert!(
                    (in_plane - cone_l).abs() < 1e-7,
                    "cone_l={cone_l} だが面内距離={in_plane}（半径が反映されていない）"
                );
            }
        }
        // 2 つの cone_l で交点集合が実際に違う（半径が位置を動かす）ことを独立に確認。
        let p_small = cone_terminator_intersections(&e, 0.35, &WGS84).expect("Ok");
        let p_large = cone_terminator_intersections(&e, 0.45, &WGS84).expect("Ok");
        let of_s = forward_of(&p_small[0], &e);
        let of_l = forward_of(&p_large[0], &e);
        let moved = !close(of_s.xi, of_l.xi, 1e-4) || !close(of_s.eta, of_l.eta, 1e-4);
        assert!(
            moved,
            "cone_l を変えても交点が動かない（半径が無視されている）"
        );
    }

    /// 決定性: 同一入力で同じ順序・同じ点を返す（実装は θ 昇順等で決定的）。
    /// 同じ呼び出しを 2 回行い、要素数・各 index の (φ,λ) が完全一致することを縛る。
    #[test]
    fn deterministic_order_and_points() {
        let e = elems(0.80, 0.10, 0.20, 1.10);
        let cone_l = 0.40_f64;
        let a = cone_terminator_intersections(&e, cone_l, &WGS84).expect("Ok");
        let b = cone_terminator_intersections(&e, cone_l, &WGS84).expect("Ok");
        assert_eq!(a.len(), b.len(), "要素数が決定的でない");
        for (pa, pb) in a.iter().zip(b.iter()) {
            // 決定性は厳密等価で確認（許容差付き close ではなく exact 一致）。
            assert!(
                pa == pb,
                "同一入力で index ごとの点が一致しない（非決定的）: a={pa:?} b={pb:?}"
            );
        }
    }

    /// terminator 楕円の扁平が効く: d≠0 で k=sin²d+cos²d/(1−f)²>1（η 方向が縮む）を反映した点になる。
    /// 球近似（k=1, 単位円）に terminator を置くと、その点を WGS84 前方射影で戻すと ζ≠0（~Re·f 残差）に
    /// なるが、本関数は WGS84 厳密 ⇒ 前方射影 ζ≈0 が機械精度で成立する。
    /// d を大きめ（0.6 rad）に取り k−1 を顕在化させ、各交点が確かに楕円 ξ²+kη²=1 上（球近似では外れる）
    /// にあることを独立に表明する。
    ///
    /// 殺す変異: k の cos²d/(1−f)² 項を cos²d（球近似）にする・(1−f)² を 1 にする・sin²d 項脱落
    ///   ⇒ terminator が単位円に化け、前方射影 ζ が機械精度 0 を外れて契約1 が落ちる。
    #[test]
    fn terminator_ellipse_flattening_matters() {
        let d = 0.60_f64; // 大きめの赤緯で扁平 k−1 を顕在化。
        let e = elems(0.70, 0.20, d, 0.9);
        let cone_l = 0.50_f64;
        // 独立に k を組み、k>1（η 方向に縮む楕円）を確認。
        let (sin_d, cos_d) = d.sin_cos();
        let omf = 1.0 - WGS84.f;
        let k = sin_d * sin_d + cos_d * cos_d / (omf * omf);
        assert!(k > 1.0, "扁平 d≠0 で k>1 のはず: k={k}");

        let pts = cone_terminator_intersections(&e, cone_l, &WGS84).expect("Ok");
        assert!(!pts.is_empty(), "交点が存在する合成のはず, got {pts:?}");
        for p in &pts {
            let of = forward_of(p, &e);
            // 契約1（WGS84 厳密）: ζ≈0 が機械精度で成立（球近似 k=1 なら ~Re·f≈3e-3 残差で落ちる）。
            assert!(
                of.zeta.abs() < 1e-7,
                "ζ={} ≈0（WGS84 厳密 terminator）。球近似なら ζ≠0 で落ちる",
                of.zeta
            );
            // 楕円 ξ²+k·η²=1 上にあることを独立に表明（球近似 ξ²+η²=1 とは k·η² 項で異なる）。
            let ellipse = of.xi * of.xi + k * of.eta * of.eta;
            assert!(
                (ellipse - 1.0).abs() < 1e-7,
                "ξ²+k·η²={ellipse} ≈1（WGS84 terminator 楕円上）でない（k={k}）"
            );
            // 錐縁上（契約2）も維持。
            let in_plane = (of.xi - e.x).hypot(of.eta - e.y);
            assert!(
                (in_plane - cone_l).abs() < 1e-7,
                "面内距離 {in_plane} = cone_l {cone_l} でない"
            );
        }
    }

    // ====================================================================
    // great_circle_distance_km の純関数ユニットテスト（高速・haversine 既知値オラクル）
    //
    // 既知 2 点→既知距離（球半径 EARTH_MEAN_RADIUS_KM=6371 km の大圏）で縛り、
    // `2.0 * EARTH_MEAN_RADIUS_KM` の `*`→`/`（半径が 2/6371 km に潰れ距離が桁違いに小さくなる）等を撃破する。
    // tol は ~0.5 km（haversine が球近似ゆえ厳密 0 ではない・WGS84 楕円体差は別途統合テスト）。
    // ====================================================================

    /// 赤道上の経度 1° 差の大圏距離 ≈ 6371·π/180 ≈ 111.195 km。
    /// `2.0 * R` の `*`→`/` は ~3.5e-3 km に潰れ tol 0.5 km を大きく外す。
    #[test]
    fn great_circle_distance_equator_one_degree_lon() {
        let a = GeoPoint::from_degrees(0.0, 0.0).expect("valid point");
        let b = GeoPoint::from_degrees(0.0, 1.0).expect("valid point");
        let expected = EARTH_MEAN_RADIUS_KM * (PI / 180.0); // 111.1949 km
        let got = great_circle_distance_km(&a, &b);
        assert!(
            (got - expected).abs() < 0.5,
            "equator 1° lon: got {got} km, expected ≈{expected} km"
        );
    }

    /// 子午線上の緯度 0°→1° の大圏距離 ≈ 6371·π/180 ≈ 111.195 km（lat 差経路も縛る）。
    #[test]
    fn great_circle_distance_one_degree_lat() {
        let a = GeoPoint::from_degrees(0.0, 0.0).expect("valid point");
        let b = GeoPoint::from_degrees(1.0, 0.0).expect("valid point");
        let expected = EARTH_MEAN_RADIUS_KM * (PI / 180.0); // 111.1949 km
        let got = great_circle_distance_km(&a, &b);
        assert!(
            (got - expected).abs() < 0.5,
            "1° lat: got {got} km, expected ≈{expected} km"
        );
    }

    /// 同一点の大圏距離 = 0（dlat=dlon=0 ⇒ h=0 ⇒ asin(0)=0）。
    #[test]
    fn great_circle_distance_same_point_is_zero() {
        let a = GeoPoint::from_degrees(35.0, 139.0).expect("valid point");
        let got = great_circle_distance_km(&a, &a);
        assert!(
            got.abs() < 1e-9,
            "same point distance must be 0, got {got} km"
        );
    }

    /// 赤道上の対蹠半周（経度 0°↔180°）= 半周 ≈ π·6371 ≈ 20015.09 km。
    /// この大スケールでも `2.0 * R` の係数 2 とスケールを縛る（`*`→`/` は ~6e-3 km）。
    #[test]
    fn great_circle_distance_antipodal_equator_is_half_circumference() {
        let a = GeoPoint::from_degrees(0.0, 0.0).expect("valid point");
        let b = GeoPoint::from_degrees(0.0, 180.0).expect("valid point");
        let expected = PI * EARTH_MEAN_RADIUS_KM; // 20015.086 km
        let got = great_circle_distance_km(&a, &b);
        assert!(
            (got - expected).abs() < 0.5,
            "antipodal equator: got {got} km, expected ≈{expected} km"
        );
    }

    /// 緯度・経度がともに非自明に異なる斜め 2 点（10°N,0° ↔ 50°N,40°E）の大圏距離 ≈ 5763.650 km。
    /// 両緯度の cos が異なり、かつ dlon≠0 で `lat1.cos() * lat2.cos()` 項が値を持つので、
    /// この積の `*`→`/`（cos(lat1)/cos(lat2)）を撃破する（赤道/同経度ケースでは cos が等しいか dlon=0 で
    /// 等価化していた）。期待値は haversine 既知式の独立手計算（R=6371 km）。
    #[test]
    fn great_circle_distance_diagonal_pair_known_value() {
        let a = GeoPoint::from_degrees(10.0, 0.0).expect("valid point");
        let b = GeoPoint::from_degrees(50.0, 40.0).expect("valid point");
        // haversine 手計算: h = sin²(Δφ/2) + cos φ1·cos φ2·sin²(Δλ/2), d = 2R·asin√h。
        let expected = 5_763.650_056_682_031_f64;
        let got = great_circle_distance_km(&a, &b);
        assert!(
            (got - expected).abs() < 0.5,
            "diagonal pair: got {got} km, expected ≈{expected} km"
        );
    }

    // ====================================================================
    // initial_bearing の純関数ユニットテスト（高速・既知方位の閉形式オラクル, M9 残(3) 3c-ii）
    //
    // 確定仕様（docs/algorithms/11-path-partial-domain.md §11.4）:
    //   θ = atan2(sinΔλ·cosφ2, cosφ1·sinφ2 − sinφ1·cosφ2·cosΔλ) を [0,2π) 正規化（Δλ=lon2−lon1）。
    //   北=0・東=π/2・南=π・西=3π/2。
    //
    // 戦略（mutation 意識・非対称オラクル）: from を非赤道・非自明な点（35°N,139°E）に置き、
    // 四方位（真北・真東・真南・真西）の既知値を絶対値で縛る。atan2 引数順の取り違え（x↔y）は
    // 90°回した値になり落ちる。sinΔλ↔cosΔλ・cosφ↔sinφ の取り違え・[0,2π) 正規化の脱落も撃つ。
    // ====================================================================

    /// 真東: 同緯度・東隣（赤道上で φ1=φ2=0・lon2>lon1）の初期方位 = π/2（東）。
    /// 赤道では大圏初期方位は厳密に真東。atan2 引数順を x↔y に取り違えると 0（北）になり落ちる。
    #[test]
    fn initial_bearing_due_east_on_equator() {
        let from = GeoPoint::from_degrees(0.0, 0.0).expect("valid");
        let to = GeoPoint::from_degrees(0.0, 10.0).expect("valid");
        let b = initial_bearing(&from, &to);
        assert!(
            close(b, PI / 2.0, 1e-9),
            "赤道で東隣の初期方位は π/2（東）, got {b}"
        );
    }

    /// 真西: 赤道上で西隣（lon2<lon1）の初期方位 = 3π/2（西・[0,2π) 正規化後）。
    /// 素の atan2 は −π/2 を返すが、[0,2π) 正規化で +3π/2。正規化脱落（−π/2 のまま）を撃つ。
    #[test]
    fn initial_bearing_due_west_normalized_to_three_half_pi() {
        let from = GeoPoint::from_degrees(0.0, 0.0).expect("valid");
        let to = GeoPoint::from_degrees(0.0, -10.0).expect("valid");
        let b = initial_bearing(&from, &to);
        assert!(
            close(b, 3.0 * PI / 2.0, 1e-9),
            "赤道で西隣の初期方位は 3π/2（西・正規化後）, got {b}"
        );
        // 正規化されて [0,2π) に収まる（負値のまま返さない）。
        assert!(
            (0.0..TAU).contains(&b),
            "方位は [0,2π) に正規化される, got {b}"
        );
    }

    /// 真北: 同経度・高緯度側（Δλ=0・φ2>φ1）の初期方位 = 0（北）。
    /// sinΔλ=0 ゆえ y 成分=0、x 成分=cosφ1·sinφ2−sinφ1·cosφ2·1=sin(φ2−φ1)>0 ⇒ atan2(0,+)=0。
    /// 非赤道 from（35°N）で sinφ1≠0・cosφ1≠0 の両項を効かせる（cosφ↔sinφ 取り違えを撃つ）。
    #[test]
    fn initial_bearing_due_north_same_meridian() {
        let from = GeoPoint::from_degrees(35.0, 139.0).expect("valid");
        let to = GeoPoint::from_degrees(60.0, 139.0).expect("valid");
        let b = initial_bearing(&from, &to);
        assert!(
            close(b, 0.0, 1e-9),
            "同経度で北側の初期方位は 0（北）, got {b}"
        );
    }

    /// 真南: 同経度・低緯度側（Δλ=0・φ2<φ1）の初期方位 = π（南）。
    /// y 成分=0、x 成分=sin(φ2−φ1)<0 ⇒ atan2(0,−)=π。非赤道 from で全項を効かせる。
    #[test]
    fn initial_bearing_due_south_same_meridian() {
        let from = GeoPoint::from_degrees(35.0, 139.0).expect("valid");
        let to = GeoPoint::from_degrees(10.0, 139.0).expect("valid");
        let b = initial_bearing(&from, &to);
        assert!(
            close(b, PI, 1e-9),
            "同経度で南側の初期方位は π（南）, got {b}"
        );
    }

    /// 非自明な斜め方位の絶対値ピン（北東向き）: from=35°N,139°E → to=45°N,150°E。
    /// 標準式から独立に手計算した θ（[0,2π) 正規化）と機械精度で一致。第1象限（0<θ<π/2）の北東。
    /// atan2 引数順（x↔y）・sinΔλ↔cosΔλ・cosφ1·sinφ2 と sinφ1·cosφ2·cosΔλ の項取り違え・符号を撃つ。
    #[test]
    fn initial_bearing_diagonal_matches_closed_form() {
        let from = GeoPoint::from_degrees(35.0, 139.0).expect("valid");
        let to = GeoPoint::from_degrees(45.0, 150.0).expect("valid");
        // 独立な閉形式（contract の式をそのまま組む。被テスト関数の戻りは使わない）。
        let phi1 = 35.0_f64.to_radians();
        let phi2 = 45.0_f64.to_radians();
        let dlon = (150.0_f64 - 139.0).to_radians();
        let y = dlon.sin() * phi2.cos();
        let x = phi1.cos() * phi2.sin() - phi1.sin() * phi2.cos() * dlon.cos();
        let expected = {
            let raw = y.atan2(x);
            raw.rem_euclid(TAU)
        };
        let b = initial_bearing(&from, &to);
        assert!(
            close(b, expected, 1e-9),
            "斜め方位 {b} が閉形式 {expected} に一致しない"
        );
        // 北東＝第1象限（0<θ<π/2）であることも独立に確認（東 π/2 と北 0 の取り違えを補強）。
        assert!(
            b > 0.0 && b < PI / 2.0,
            "北東向きは第1象限 (0,π/2), got {b}"
        );
    }

    /// 方位の向き非対称性: bearing(A→B) と bearing(B→A) は同じでない（東隣の往復で π/2 vs 3π/2）。
    /// Δλ=lon2−lon1 の引き算方向（lon1↔lon2 取り違え＝Δλ の符号反転）を撃つ。
    #[test]
    fn initial_bearing_is_directional_not_symmetric() {
        let a = GeoPoint::from_degrees(0.0, 0.0).expect("valid");
        let b = GeoPoint::from_degrees(0.0, 10.0).expect("valid");
        let ab = initial_bearing(&a, &b);
        let ba = initial_bearing(&b, &a);
        assert!(close(ab, PI / 2.0, 1e-9), "A→B は東 π/2, got {ab}");
        assert!(close(ba, 3.0 * PI / 2.0, 1e-9), "B→A は西 3π/2, got {ba}");
        // 明示的に非対称（Δλ 符号反転＝両者が一致する変異を撃つ）。
        assert!(
            !close(ab, ba, 1e-6),
            "方位は向き依存（A→B ≠ B→A）でなければならない: {ab} vs {ba}"
        );
    }
}
