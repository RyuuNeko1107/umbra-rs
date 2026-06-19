//! 太陽の地平座標（高度・方位）と可視性判定（`docs/issues/ISSUE-028`、`docs/conventions.md` §7,
//! `docs/physical-models.md` §C1/§C2, Meeus *Astronomical Algorithms* Ch.13/§16）。
//!
//! - 方位角: **北=0°、東回り**（北→東→南→西, conventions §7）。
//! - 高度: 既定**幾何学的高度**。大気差は `RefractionModel { None, Standard }` で分離し、補正前後を両方返す。
//!   `Standard` は Saemundsson（真→見かけ, physical-models §C1.1）標準条件。
//! - 可視性: 各接触/最大時点の太陽高度から 6 値（api-draft §3.4）。
//!
//! 赤道→地平変換: 局地時角 `H = (θ_ERA + λ_east) − α_sun`（CIO ベース・分点 GST 禁止, D4。
//! ISSUE-024 の時角符号 H=μ+λ_east と整合）。太陽見かけ赤経赤緯は CIRS（ISSUE-015）。
//!
//! 注: 接触点の位置角（position angle）は `LocalContact` 拡張（ISSUE-043 結線）に紐づくため本層では
//! 扱わない（ISSUE-025/026 で LocalContact を時刻のみに絞ったのと一貫）。

// pub(crate) 関数は ISSUE-043（EclipseEngine 結線）が消費するまで未使用。
#![allow(dead_code)]

use umbra_core::deltat::{tt_to_ut1, DeltaTModel};
use umbra_core::{Degrees, Radians, TtInstant};
use umbra_ephemeris::apparent::sun_apparent_cirs;
use umbra_ephemeris::frames::earth_rotation_angle;

/// 可視性判定の地平閾値（幾何学的高度, 度）。conventions §7 既定（physical-models §C2.1 の
/// −0.8333°＝太陽縁＋地平大気差は代替）。`< HORIZON_ALTITUDE_DEG` で地平下。
const HORIZON_ALTITUDE_DEG: f64 = 0.0;

/// 大気差モデル（physical-models §C1）。
#[non_exhaustive]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RefractionModel {
    /// 大気差なし（幾何学的高度のみ）。
    None,
    /// 標準大気差（Saemundsson 真→見かけ・標準条件 1010 hPa/10 ℃, physical-models §C1.1）。
    Standard,
}

/// 太陽の地平座標。高度は幾何学的（既定）と大気差補正後の両方、方位は北0東回り。
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Horizontal {
    /// 幾何学的高度（既定, conventions §7）。
    pub altitude_geometric: Degrees,
    /// 大気差補正後の見かけ高度（`RefractionModel::Standard` 時。`None` 時は幾何学的高度と同値）。
    pub altitude_apparent: Degrees,
    /// 方位角（北=0°、東回り, conventions §7）。範囲 `[0, 360)`。
    pub azimuth: Degrees,
}

/// 食の可視性（api-draft §3.4）。
#[non_exhaustive]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize))]
#[cfg_attr(feature = "serde", serde(tag = "type"))]
pub enum Visibility {
    /// その地点で食域外（接触なし）。
    NotVisible,
    /// 食全体が地平下（最大食も両接触も地平下＝観測不能。日の出/日没で端の相が地平上なら
    /// `SunriseEclipse`/`SunsetEclipse` となり本値にはならない）。
    BelowHorizon,
    /// 日の出中に食が進行（食開始 C1 が地平下、最大以降が地平上）。
    SunriseEclipse,
    /// 日没中に食が進行（食終了 C4 が地平下）。
    SunsetEclipse,
    /// 一部の接触のみ地平上（部分的に観測可能）。
    PartialVisible,
    /// C1〜C4 すべて地平上（全経過が観測可能）。
    FullyVisible,
}

/// 観測者（測地緯度・東経）と TT・大気差モデルから太陽の地平座標を求める。
///
/// - `geodetic_latitude`/`east_longitude`: 観測者（東経正, conventions §3）。
/// - `time_tt`: 評価時刻（TT）。`delta_t` で UT1 に変換し ERA（CIO 時角）を構成。
/// - `refraction`: `Standard` で見かけ高度に Saemundsson 大気差を加える（幾何学的高度も併せて返す）。
pub(crate) fn sun_horizontal<M: DeltaTModel>(
    geodetic_latitude: Radians,
    east_longitude: Radians,
    time_tt: TtInstant,
    refraction: RefractionModel,
    delta_t: &M,
) -> Horizontal {
    // 太陽見かけ地心位置（CIRS, ISSUE-015）→ 赤経 α・赤緯 δ（Meeus Ch.13）。
    let sun = sun_apparent_cirs(time_tt);
    let radius = (sun.x * sun.x + sun.y * sun.y + sun.z * sun.z).sqrt();
    let alpha = sun.y.atan2(sun.x);
    let delta = (sun.z / radius).asin();

    // 局地時角 H = (θ_ERA + λ_east) − α_sun（CIO ベース・D4。ISSUE-024 と整合）。
    let ut1 = tt_to_ut1(time_tt, delta_t);
    let era = earth_rotation_angle(ut1).0;
    let hour_angle = era + east_longitude.0 - alpha;

    let phi = geodetic_latitude.0;
    let (sin_phi, cos_phi) = phi.sin_cos();
    let (sin_d, cos_d) = delta.sin_cos();
    let (sin_h, cos_h) = hour_angle.sin_cos();

    // 高度（幾何学的）: sin a = sinφ sinδ + cosφ cosδ cosH。
    let altitude = (sin_phi * sin_d + cos_phi * cos_d * cos_h).asin();
    // 方位（北=0°、東回り, conventions §7）: A = atan2(−east, north)、east = cosδ sinH,
    // north = sinδ cosφ − cosδ cosH sinφ。これは南基準 A_south + π と恒等だが、+π の mod 2π
    // 等価変異を避けるため north/east 成分から直接 atan2 する（符号 − が load-bearing）。
    let east = cos_d * sin_h;
    let north = sin_d * cos_phi - cos_d * cos_h * sin_phi;
    let azimuth = Radians::new((-east).atan2(north)).normalized_two_pi();

    let altitude_geometric = Radians(altitude).to_degrees();
    let altitude_apparent = match refraction {
        RefractionModel::None => altitude_geometric,
        RefractionModel::Standard => {
            Degrees(altitude_geometric.0 + saemundsson_refraction_deg(altitude_geometric.0))
        }
    };

    Horizontal {
        altitude_geometric,
        altitude_apparent,
        azimuth: azimuth.to_degrees(),
    }
}

/// Saemundsson (1986) 大気差（真→見かけ, physical-models §C1.1, Meeus Ch.16）。標準条件 1010hPa/10℃。
/// 幾何学的高度 `h`[度] に対し `R[arcmin] = 1.02 / tan((h + 10.3/(h + 5.11))°)`、戻り値は度（R/60）。
/// 地平線付近（h≲5°）は非保証（accuracy.md §6）だが式はそのまま適用する。
fn saemundsson_refraction_deg(altitude_deg: f64) -> f64 {
    let arg_deg = altitude_deg + 10.3 / (altitude_deg + 5.11);
    let r_arcmin = 1.02 / arg_deg.to_radians().tan();
    r_arcmin / 60.0
}

/// 接触集合と最大食時点の太陽高度から可視性を分類する。
///
/// - `in_eclipse`: その地点が食域内か（C1/C4 が存在するか）。`false` なら `NotVisible`。
/// - `altitude_c1`/`altitude_c4`: C1/C4 時点の太陽幾何学的高度（部分食地点で内接が無くても外接はある）。
/// - `altitude_max`: 最大食時点の太陽幾何学的高度。
///
/// 地平閾値は幾何学的高度 0°（conventions §7 既定。physical-models §C2.1 の −0.8333° は代替）。
pub(crate) fn classify_visibility(
    in_eclipse: bool,
    altitude_c1: Option<Degrees>,
    altitude_max: Degrees,
    altitude_c4: Option<Degrees>,
) -> Visibility {
    if !in_eclipse {
        return Visibility::NotVisible;
    }
    let c1_above = altitude_c1.is_some_and(|a| a.0 >= HORIZON_ALTITUDE_DEG);
    let c4_above = altitude_c4.is_some_and(|a| a.0 >= HORIZON_ALTITUDE_DEG);
    let c1_below = altitude_c1.is_some_and(|a| a.0 < HORIZON_ALTITUDE_DEG);
    let c4_below = altitude_c4.is_some_and(|a| a.0 < HORIZON_ALTITUDE_DEG);

    if altitude_max.0 < HORIZON_ALTITUDE_DEG {
        // 幾何的最大食は地平下。ただし太陽が食の最中に地平を跨ぐ日の出/日没食では端の相が
        // 観測可能なため一律 BelowHorizon にはしない（観測可能最大はオラクルと定義が異なるが、
        // 可視性クラスは観測可能性で判定する。accuracy.md §2.1E）。
        // C4 地平上＝最大食後に太陽が昇る（日の出食・後半可視）、C1 地平上＝最大食前に沈む
        // （日没食・前半可視）。いずれも無ければ食全体が地平下＝観測不能。
        if c4_above {
            Visibility::SunriseEclipse
        } else if c1_above {
            Visibility::SunsetEclipse
        } else {
            Visibility::BelowHorizon
        }
    } else if c1_above && c4_above {
        // C1〜C4 すべて地平上。
        Visibility::FullyVisible
    } else if c1_below {
        // 食開始が地平下＝日の出中に進行（c1<0 を c4<0 より先に判定）。
        Visibility::SunriseEclipse
    } else if c4_below {
        // 食終了が地平下＝日没中に進行。
        Visibility::SunsetEclipse
    } else {
        // c1/c4 が窓外（None）等で一部のみ地平上。
        Visibility::PartialVisible
    }
}

#[cfg(test)]
mod tests {
    //! ISSUE-028 受け入れテスト（strict・太陽地平座標 / 大気差 / 可視性）。
    //!
    //! ## オラクル戦略（追認回避）
    //!
    //! 1. **alt/az 変換の独立再実装（主オラクル）**: テスト側に **ENU（East-North-Up）回転行列**
    //!    経由の equatorial→horizontal 変換を組み、`sun_horizontal` と一致を縛る。契約式は
    //!    高度 `sin a = sinφ sinδ + cosφ cosδ cosH`・方位 `A_south = atan2(cosδ sinH,
    //!    cosδ cosH sinφ − sinδ cosφ)`（南基準 → 北基準 +π）だが、テスト側は **観測者局所 ENU 単位
    //!    ベクトルを直接組む別表現**（同じ式の写しでない）:
    //!      天頂方向の地心赤道系単位ベクトル s_dir =(cosδ cosH', cosδ sinH'... ) ではなく、
    //!      時角系の east/north/up 基底に太陽方向ベクトルを射影して
    //!      `alt = asin(up)`, `az = atan2(−east, north)`（北0東回り）を得る。
    //!    α,δ,H は **公開プリミティブ** `sun_apparent_cirs`（ISSUE-015）＋
    //!    `earth_rotation_angle`（ISSUE-039）＋ `tt_to_ut1`（ISSUE-007）から取得する。これらは
    //!    solver 内部関数ではないため追認にならない（`sun_horizontal` 本体はこれらの合成・規約変換・
    //!    大気差適用を担う）。
    //!    出典: Meeus *Astronomical Algorithms* 2nd ed. Ch.13（式 13.5/13.6）/ Explanatory
    //!    Supplement §7。ENU 基底表現は球面天文の標準（例: SOFA `iauHd2ae` の幾何と等価だが
    //!    コードは移植せず独立に組む）。
    //!
    //! 2. **物理慣習チェック（独立）**: 任意実時刻で独立オラクルと一致を主軸にしつつ、
    //!    - azimuth ∈ [0,360)、altitude ≤ 90°。
    //!    - λ → λ+2π で alt/az 不変（H が 2π 周期）。
    //!    - H の符号と方位の整合（独立オラクル側の H から 真南正中 H=0 → A=180°,
    //!      真東 H<0 側 → 0<A<180°（東半分）, 真西 H>0 側 → 180<A<360°（西半分））。
    //!
    //! 3. **大気差（独立式・Saemundsson）**: 幾何高度 h[度] から
    //!    `R[arcmin] = 1.02 / tan((h + 10.3/(h+5.11))°)`、`apparent = geometric + R/60`[度]
    //!    （physical-models §C1.1, Meeus Ch.16）。実時刻の太陽幾何高度に対し `apparent−geometric`
    //!    が独立 Saemundsson 値に一致（≤1e-9°）。`None` では差 0。
    //!
    //! 4. **可視性 6 値**: `classify_visibility` は純ロジック（高度入力）。契約の判定木
    //!    （prompt / api-draft §3.4・physical-models §C2）を **6 値すべて＋境界**（altitude=0 ちょうど、
    //!    c1/c4=None 等）で網羅。閾値は **幾何学的高度 0°**（conventions §7 既定）。
    //!
    //! 注: 位置角（PA）は本 issue 非目的（LocalContact 拡張・ISSUE-043 結線）のため扱わない。
    //! pyerfa は使わず、テスト側 ENU 独立再実装を主オラクルとする（liberfa 不要・決定的）。

    // テスト doc の和文インデント（手順説明）が doc list lint と誤認される／可視性フィクスチャの
    // タプル型は test 専用で type_complexity を許容（ロジック非関与）。
    #![allow(
        clippy::doc_overindented_list_items,
        clippy::doc_lazy_continuation,
        clippy::type_complexity
    )]

    use super::*;

    use umbra_core::deltat::{tt_to_ut1, EspenakMeeusDeltaT};
    use umbra_core::{JulianDate2, Radians, TtInstant};
    use umbra_ephemeris::apparent::sun_apparent_cirs;
    use umbra_ephemeris::frames::earth_rotation_angle;

    /// 角度比較の許容（度）。alt/az は太陽位置(ISSUE-015)・ERA(ISSUE-039)由来の独立再計算と
    /// 完全に同じプリミティブを使うため、合成・規約変換のみの差。1e-9° で「式の同一性」を縛る。
    const TOL_DEG: f64 = 1e-9;

    /// TT を 2 要素 JD から構築（local_maximum.rs / local_contacts.rs テストと同形）。
    fn tt(jd1: f64, jd2: f64) -> TtInstant {
        TtInstant::from_jd2(JulianDate2::new(jd1, jd2))
    }

    /// 2017-08-21 皆既日食 最大食付近の TT エポック（local_maximum.rs と同一）。
    /// 北米中緯度で最大食（太陽が地平上・高度 ~60°級）の実時刻。
    fn tt_2017_max() -> TtInstant {
        tt(2_457_986.5, 7.685_322_222_222_222e-1)
    }

    /// もう一つの実時刻エポック（最大食 −3h ≒ 早朝・低高度側）。慣習・低高度大気差の確認に使う。
    fn tt_2017_minus_3h() -> TtInstant {
        let jd = tt_2017_max().jd2().jd() - 3.0 / 24.0;
        TtInstant::from_jd2(JulianDate2::from_jd(jd))
    }

    // ============================================================
    // 独立オラクル: 赤道→地平変換を ENU 回転行列で別経路再構成
    // ============================================================

    /// 太陽の見かけ CIRS 赤経・赤緯を公開プリミティブから取得（α,δ）。Meeus Ch.13:
    /// α = atan2(y, x), δ = asin(z/|v|)。`sun_apparent_cirs` は ISSUE-015 検証済。
    fn sun_alpha_delta(time_tt: TtInstant) -> (f64, f64) {
        let v = sun_apparent_cirs(time_tt);
        let r = (v.x * v.x + v.y * v.y + v.z * v.z).sqrt();
        let alpha = v.y.atan2(v.x); // [-π, π]
        let delta = (v.z / r).asin(); // [-π/2, π/2]
        (alpha, delta)
    }

    /// 局地時角 H = (θ_ERA + λ_east) − α_sun（CIO ベース・契約 D4）。
    /// θ_ERA は ut1=tt_to_ut1(time_tt) から、α_sun は上の独立取得から。
    fn hour_angle(east_longitude: f64, time_tt: TtInstant) -> f64 {
        let dt = EspenakMeeusDeltaT;
        let ut1 = tt_to_ut1(time_tt, &dt);
        let era = earth_rotation_angle(ut1).0;
        let (alpha, _delta) = sun_alpha_delta(time_tt);
        era + east_longitude - alpha
    }

    /// 赤道→地平の **独立再実装（ENU 回転行列・北0東回り）**。契約の atan2(A_south)+π 表現とは
    /// 別経路: 太陽の地心赤道(時角系)単位方向ベクトルを観測者局所 ENU 基底へ射影する。
    ///
    /// 時角系赤道直交基底で太陽方向（単位）。H は東向きを正に取る局地時角
    /// （H = θ_ERA + λ_east − α, 契約 D4）。x は子午線（真南）方向の水平射影に対応する極軸直交成分、
    /// y は **東方向**（H>0=午後で太陽が西へ動く＝y 成分が負）、z は極方向:
    ///   p = (cosδ cosH, cosδ sinH, sinδ)    （x:子午線方向, y:東方向, z:極方向）
    /// 観測者局所 ENU 基底（緯度 φ）を同じ系で:
    ///   Up    = ( cosφ, 0, sinφ )
    ///   North = (−sinφ, 0, cosφ )
    ///   East  = ( 0, 1, 0 )
    /// すると:
    ///   up    = p·Up    = sinφ sinδ + cosφ cosδ cosH        （= 契約の sin a。独立に導出）
    ///   north = p·North = cosφ sinδ − sinφ cosδ cosH
    ///   east  = p·East  = cosδ sinH
    /// 方位は **北0東回り**（契約 D4・conventions §7）。北0東回りでは方位は北から東へ向かって
    /// 増えるが、`atan2(east, north)` は北から東で減る（東西鏡像＝北0西回り）。よって east 射影の
    /// 符号を反転し:
    ///   alt = asin(up),  az = atan2(−east, north)  ∈ [0,2π)（北0東回り）
    ///
    /// これにより H<0（午前・太陽は地理的に東）で azimuth が東半分(0,180)、H=0（正中）で 180°
    /// （east=0 ゆえ符号反転の影響なし）、H>0（午後・西）で西半分(180,360) になり、契約
    /// `A_north = A_south + π`（A_south = atan2(cosδ sinH, cosδ cosH sinφ − sinδ cosφ)）と一致する。
    /// alt 式は契約と一致するのが物理的必然（同じ高度の定義）。
    /// 出典: Meeus Ch.13 / Explanatory Supplement §7（座標変換の標準）。
    fn horizontal_oracle(lat: f64, east_longitude: f64, time_tt: TtInstant) -> (f64, f64) {
        let (_alpha, delta) = sun_alpha_delta(time_tt);
        let h = hour_angle(east_longitude, time_tt);
        let (sinp, cosp) = (lat.sin(), lat.cos());
        let (sind, cosd) = (delta.sin(), delta.cos());
        let (sinh, cosh) = (h.sin(), h.cos());

        let up = sinp * sind + cosp * cosd * cosh;
        let north = cosp * sind - sinp * cosd * cosh;
        let east = cosd * sinh;

        let alt = up.clamp(-1.0, 1.0).asin();
        // 北0東回り: east 射影の符号を反転（atan2(east, north) は北0西回りになるため）。
        let az = (-east).atan2(north).rem_euclid(2.0 * std::f64::consts::PI);
        (alt.to_degrees(), az.to_degrees())
    }

    /// Saemundsson 大気差（真→見かけ, physical-models §C1.1, Meeus Ch.16・標準条件）を独立計算。
    /// R[arcmin] = 1.02 / tan((h + 10.3/(h+5.11))°)；apparent − geometric = R/60 [度]。
    fn saemundsson_delta_deg(h_geom_deg: f64) -> f64 {
        let arg_deg = h_geom_deg + 10.3 / (h_geom_deg + 5.11);
        let r_arcmin = 1.02 / arg_deg.to_radians().tan();
        r_arcmin / 60.0
    }

    fn deg_lat(d: f64) -> Radians {
        Radians::new(d.to_radians())
    }
    fn deg_lon(d: f64) -> Radians {
        Radians::new(d.to_radians())
    }

    fn dt() -> EspenakMeeusDeltaT {
        EspenakMeeusDeltaT
    }

    // ============================================================
    // alt/az: 独立 ENU オラクルとの一致（主オラクル）
    // ============================================================

    /// 主オラクル（最強）: 複数の実時刻・地点で `sun_horizontal`（None）の幾何高度・方位が
    /// 独立 ENU 再実装と 1e-9° 一致する。北0東回り変換・時角合成・幾何高度を一括で縛る。
    #[test]
    fn geometric_altaz_matches_independent_enu_oracle() {
        // (lat_deg, east_lon_deg, time)
        let cases: &[(f64, f64, TtInstant)] = &[
            (37.5, -89.2, tt_2017_max()),      // 北米中緯度・最大食（高高度）
            (37.5, -89.2, tt_2017_minus_3h()), // 同地点・早朝（低高度）
            (60.0, -100.0, tt_2017_max()),     // 高緯度
            (-33.9, 18.4, tt_2017_max()),      // 南半球（ケープタウン付近, 食域外でも幾何は定義可）
            (0.0, 120.0, tt_2017_max()),       // 赤道・東経
        ];
        for &(lat_d, lon_d, t) in cases {
            let h = sun_horizontal(
                deg_lat(lat_d),
                deg_lon(lon_d),
                t,
                RefractionModel::None,
                &dt(),
            );
            let (alt_o, az_o) = horizontal_oracle(lat_d.to_radians(), lon_d.to_radians(), t);
            assert!(
                (h.altitude_geometric.0 - alt_o).abs() < TOL_DEG,
                "alt mismatch lat={lat_d} lon={lon_d}: got {} want {alt_o}",
                h.altitude_geometric.0
            );
            // 方位は周期量。北0近傍の巻き戻りを跨ぐ差は ±360 で畳む。
            let mut daz = (h.azimuth.0 - az_o).rem_euclid(360.0);
            if daz > 180.0 {
                daz -= 360.0;
            }
            assert!(
                daz.abs() < TOL_DEG,
                "az mismatch lat={lat_d} lon={lon_d}: got {} want {az_o}",
                h.azimuth.0
            );
        }
    }

    /// プロパティ（L8）: azimuth ∈ [0,360)、altitude ≤ 90°、altitude ≥ −90°。複数地点で。
    #[test]
    fn altitude_and_azimuth_in_valid_ranges() {
        let cases: &[(f64, f64, TtInstant)] = &[
            (37.5, -89.2, tt_2017_max()),
            (60.0, -100.0, tt_2017_minus_3h()),
            (-33.9, 18.4, tt_2017_max()),
            (89.0, 0.0, tt_2017_max()),
            (0.0, 179.9, tt_2017_max()),
        ];
        for &(lat_d, lon_d, t) in cases {
            let h = sun_horizontal(
                deg_lat(lat_d),
                deg_lon(lon_d),
                t,
                RefractionModel::None,
                &dt(),
            );
            assert!(
                (0.0..360.0).contains(&h.azimuth.0),
                "azimuth {} out of [0,360) lat={lat_d} lon={lon_d}",
                h.azimuth.0
            );
            assert!(
                h.altitude_geometric.0 <= 90.0 + 1e-9 && h.altitude_geometric.0 >= -90.0 - 1e-9,
                "altitude {} out of [-90,90] lat={lat_d} lon={lon_d}",
                h.altitude_geometric.0
            );
        }
    }

    /// プロパティ: λ → λ + 2π で alt/az 不変（時角の 2π 周期。東経の循環）。
    #[test]
    fn longitude_plus_two_pi_is_invariant() {
        let lat = deg_lat(37.5);
        let lon = deg_lon(-89.2);
        let lon_wrapped = Radians::new(lon.0 + 2.0 * std::f64::consts::PI);
        let t = tt_2017_max();
        let a = sun_horizontal(lat, lon, t, RefractionModel::None, &dt());
        let b = sun_horizontal(lat, lon_wrapped, t, RefractionModel::None, &dt());
        assert!(
            (a.altitude_geometric.0 - b.altitude_geometric.0).abs() < TOL_DEG,
            "altitude must be invariant under λ→λ+2π: {} vs {}",
            a.altitude_geometric.0,
            b.altitude_geometric.0
        );
        let mut daz = (a.azimuth.0 - b.azimuth.0).rem_euclid(360.0);
        if daz > 180.0 {
            daz -= 360.0;
        }
        assert!(
            daz.abs() < TOL_DEG,
            "azimuth must be invariant under λ→λ+2π: {} vs {}",
            a.azimuth.0,
            b.azimuth.0
        );
    }

    // ============================================================
    // 方位規約: 北0東回り（時角符号 × 方位符号 交差検証, 契約 D4 必須）
    // ============================================================

    /// 時角符号 H=(θ_ERA+λ)−α と方位（北0東回り）の交差検証（独立オラクル側の H で判定）:
    ///  - 真南正中（|H|≈0, 北半球 δ<φ）→ azimuth ≈ 180°、かつ alt ≈ 90−(φ−δ)。
    /// 子午線通過時刻は α_sun と λ から逆算（H=θ_ERA+λ−α=0 を満たす λ を解析的に選ぶ）。
    /// 任意の固定時刻で α,θ_ERA を独立に得て、H=0 になる東経 λ* = α − θ_ERA を作れば
    /// その地点で必ず正中する（実装の内部 H に依存せず、独立に構成した正中条件）。
    #[test]
    fn meridian_transit_azimuth_is_south_180_north_hemisphere() {
        let t = tt_2017_max();
        let dt_m = dt();
        let ut1 = tt_to_ut1(t, &dt_m);
        let era = earth_rotation_angle(ut1).0;
        let (alpha, delta) = sun_alpha_delta(t);
        // H = era + λ − α = 0 → λ* = α − era（[-π,π) へ正規化）。
        let lon_star = Radians::new(alpha - era).normalized_signed();
        // 北半球で δ < φ となる緯度（夏至近傍 δ≈+11.8°、φ=40°）。真南正中で A=180°。
        let lat_deg = 40.0;
        let h = sun_horizontal(deg_lat(lat_deg), lon_star, t, RefractionModel::None, &dt_m);
        // 独立確認: この λ* で時角は ≈0。
        let hangle = (era + lon_star.0 - alpha).rem_euclid(2.0 * std::f64::consts::PI);
        let hangle_signed = if hangle > std::f64::consts::PI {
            hangle - 2.0 * std::f64::consts::PI
        } else {
            hangle
        };
        assert!(
            hangle_signed.abs() < 1e-9,
            "precondition: constructed longitude must put sun on meridian (H≈0), got H={hangle_signed} rad"
        );
        // 北半球・δ<φ → 真南 → 北0東回りで 180°。
        let mut daz = (h.azimuth.0 - 180.0).rem_euclid(360.0);
        if daz > 180.0 {
            daz -= 360.0;
        }
        assert!(
            daz.abs() < 1e-6,
            "northern-hemisphere meridian transit (δ<φ) azimuth must be ~180° (south), got {}",
            h.azimuth.0
        );
        // 正中高度 = 90 − (φ − δ)。
        let alt_expected = 90.0 - (lat_deg - delta.to_degrees());
        assert!(
            (h.altitude_geometric.0 - alt_expected).abs() < 1e-6,
            "meridian altitude must be 90−(φ−δ)={alt_expected}, got {}",
            h.altitude_geometric.0
        );
    }

    /// 北0東回りの東西非対称: 正中時刻からわずかに東へずらした地点（H<0, 太陽は東寄り＝午前側）で
    /// azimuth が東半分（0<A<180）、西へずらした地点（H>0, 午後側）で西半分（180<A<360）になる。
    /// H の符号を独立オラクルで決め、方位の象限が北0東回りと整合することを縛る（交差検証）。
    #[test]
    fn azimuth_east_west_quadrant_matches_hour_angle_sign() {
        let t = tt_2017_max();
        let dt_m = dt();
        let ut1 = tt_to_ut1(t, &dt_m);
        let era = earth_rotation_angle(ut1).0;
        let (alpha, _delta) = sun_alpha_delta(t);
        let lon_meridian = alpha - era; // H=0 の東経（rad, 未正規化でよい）
        let lat = deg_lat(40.0);

        // λ を西（−）にずらす → H=era+λ−α が小さく（負）→ 太陽は東側（午前）。
        let lon_west_shift = Radians::new(lon_meridian - 20.0_f64.to_radians());
        // λ を東（+）にずらす → H>0 → 太陽は西側（午後）。
        let lon_east_shift = Radians::new(lon_meridian + 20.0_f64.to_radians());

        let h_neg = sun_horizontal(lat, lon_west_shift, t, RefractionModel::None, &dt_m);
        let h_pos = sun_horizontal(lat, lon_east_shift, t, RefractionModel::None, &dt_m);

        // 独立 H 符号確認。
        let hsign = |lon: f64| {
            let x = (era + lon - alpha).rem_euclid(2.0 * std::f64::consts::PI);
            if x > std::f64::consts::PI {
                x - 2.0 * std::f64::consts::PI
            } else {
                x
            }
        };
        assert!(
            hsign(lon_west_shift.0) < 0.0,
            "precondition: west-shifted longitude must give H<0 (sun in the east)"
        );
        assert!(
            hsign(lon_east_shift.0) > 0.0,
            "precondition: east-shifted longitude must give H>0 (sun in the west)"
        );

        assert!(
            (0.0..180.0).contains(&h_neg.azimuth.0),
            "H<0 (morning, sun east) must give azimuth in eastern half (0,180), got {}",
            h_neg.azimuth.0
        );
        assert!(
            (180.0..360.0).contains(&h_pos.azimuth.0),
            "H>0 (afternoon, sun west) must give azimuth in western half (180,360), got {}",
            h_pos.azimuth.0
        );
    }

    // ============================================================
    // 大気差: Saemundsson 独立式との一致 / None は差0
    // ============================================================

    /// `RefractionModel::None` では altitude_apparent == altitude_geometric（補正なし）。
    #[test]
    fn refraction_none_apparent_equals_geometric() {
        let cases: &[(f64, f64, TtInstant)] = &[
            (37.5, -89.2, tt_2017_max()),
            (37.5, -89.2, tt_2017_minus_3h()),
            (60.0, -100.0, tt_2017_max()),
        ];
        for &(lat_d, lon_d, t) in cases {
            let h = sun_horizontal(
                deg_lat(lat_d),
                deg_lon(lon_d),
                t,
                RefractionModel::None,
                &dt(),
            );
            assert_eq!(
                h.altitude_apparent.0, h.altitude_geometric.0,
                "RefractionModel::None must leave apparent == geometric (lat={lat_d})"
            );
        }
    }

    /// `RefractionModel::Standard`: apparent − geometric が独立 Saemundsson 値に 1e-9° 一致。
    /// 幾何高度自体は None と同一（補正は apparent のみに乗る）。h>5° の妥当域（高高度・低高度
    /// 両方）の実時刻で縛る。標準条件 1010hPa/10℃（係数1）。
    #[test]
    fn standard_refraction_matches_saemundsson_independent_formula() {
        let cases: &[(f64, f64, TtInstant)] = &[
            (37.5, -89.2, tt_2017_max()),      // 高高度（~60°級）
            (37.5, -89.2, tt_2017_minus_3h()), // 低高度（朝・数十度以下）
            (20.0, -100.0, tt_2017_max()),     // 別緯度の昼側（亜太陽点 ≈11.9°N,96°W 付近で高高度）
        ];
        for &(lat_d, lon_d, t) in cases {
            let none = sun_horizontal(
                deg_lat(lat_d),
                deg_lon(lon_d),
                t,
                RefractionModel::None,
                &dt(),
            );
            let std = sun_horizontal(
                deg_lat(lat_d),
                deg_lon(lon_d),
                t,
                RefractionModel::Standard,
                &dt(),
            );

            // 幾何高度は屈折モデルに依存しない。
            assert!(
                (std.altitude_geometric.0 - none.altitude_geometric.0).abs() < TOL_DEG,
                "geometric altitude must be independent of refraction model (lat={lat_d})"
            );
            // 方位も屈折で不変。
            assert!(
                (std.azimuth.0 - none.azimuth.0).abs() < TOL_DEG,
                "azimuth must be independent of refraction model (lat={lat_d})"
            );

            let h_geom = std.altitude_geometric.0;
            // 妥当域（地平線付近は非保証）に限定して縛る。
            assert!(
                h_geom > 5.0,
                "test precondition: chosen geometric altitude {h_geom}° should be >5° (Saemundsson valid range, lat={lat_d})"
            );
            let expected = saemundsson_delta_deg(h_geom);
            let actual = std.altitude_apparent.0 - std.altitude_geometric.0;
            assert!(
                (actual - expected).abs() < TOL_DEG,
                "Standard refraction apparent−geometric={actual}° must match Saemundsson={expected}° \
                 (h_geom={h_geom}°, lat={lat_d})"
            );
            // 大気差は正（真→見かけで持ち上がる）。
            assert!(
                actual > 0.0,
                "Saemundsson refraction must lift apparent above geometric (got {actual}°, lat={lat_d})"
            );
        }
    }

    /// 大気差の単調性チェック（独立式の健全性, 観点補強）: 低高度ほど補正が大きい。
    /// 実装出力（高高度 vs 低高度の同地点）で apparent−geometric が低高度の方が大きい。
    #[test]
    fn refraction_is_larger_at_lower_altitude() {
        let lat = deg_lat(37.5);
        let lon = deg_lon(-89.2);
        let high = sun_horizontal(lat, lon, tt_2017_max(), RefractionModel::Standard, &dt());
        let low = sun_horizontal(
            lat,
            lon,
            tt_2017_minus_3h(),
            RefractionModel::Standard,
            &dt(),
        );
        // 前提: 朝の方が高度が低い。
        assert!(
            low.altitude_geometric.0 < high.altitude_geometric.0,
            "precondition: morning altitude {} should be below midday {}",
            low.altitude_geometric.0,
            high.altitude_geometric.0
        );
        let r_high = high.altitude_apparent.0 - high.altitude_geometric.0;
        let r_low = low.altitude_apparent.0 - low.altitude_geometric.0;
        assert!(
            r_low > r_high,
            "refraction must be larger at lower altitude: low={r_low}° high={r_high}°"
        );
    }

    /// 大気差の単調性（決定的・S1）: テスト側純関数 `saemundsson_delta_deg` を固定高度 5° と 45° で
    /// 直接比較する。実時刻2点に依存する `refraction_is_larger_at_lower_altitude` は太陽位置の変動で
    /// 将来 flaky 化しうるため、独立式そのものの単調性（h 小ほど補正大・正値）を決定的に縛る。
    #[test]
    fn refraction_monotonic_low_vs_high_fixed_altitudes() {
        let r_low = saemundsson_delta_deg(5.0);
        let r_high = saemundsson_delta_deg(45.0);
        assert!(
            r_low > 0.0 && r_high > 0.0,
            "Saemundsson refraction must be positive at both 5° and 45° (low={r_low}°, high={r_high}°)"
        );
        assert!(
            r_low > r_high,
            "Saemundsson refraction must be larger at lower altitude (5°={r_low}° must exceed 45°={r_high}°)"
        );
    }

    // ============================================================
    // classify_visibility: 6 値網羅 + 境界（閾値=幾何高度 0°）
    // ============================================================

    /// 食域外（in_eclipse=false）→ NotVisible。altitude_max が地平上でも NotVisible が優先。
    #[test]
    fn visibility_not_in_eclipse_is_not_visible() {
        assert_eq!(
            classify_visibility(
                false,
                Some(Degrees(30.0)),
                Degrees(45.0),
                Some(Degrees(20.0))
            ),
            Visibility::NotVisible,
            "in_eclipse=false must yield NotVisible regardless of altitudes"
        );
        // c1/c4 None でも NotVisible。
        assert_eq!(
            classify_visibility(false, None, Degrees(-10.0), None),
            Visibility::NotVisible
        );
    }

    /// 最大食も地平下（altitude_max < 0）→ BelowHorizon（全接触地平下）。
    #[test]
    fn visibility_max_below_horizon_is_below_horizon() {
        assert_eq!(
            classify_visibility(
                true,
                Some(Degrees(-20.0)),
                Degrees(-5.0),
                Some(Degrees(-1.0))
            ),
            Visibility::BelowHorizon,
            "altitude_max<0 must yield BelowHorizon"
        );
    }

    /// 受け入れ点（FIX 主題）: 最大食が地平直下だが C4 が地平上 → SunriseEclipse。
    /// Caribou ME 2025-03-29 型（幾何 max≈−0.05°, C4 地平上）。太陽が max〜C4 の間に昇り
    /// 後半相が観測可能なので SunriseEclipse。現行コード（max<0 で即 BelowHorizon）は
    /// BelowHorizon を返すため red。撃破する変異: 146 行「max<0→BelowHorizon early return」を
    /// 接触考慮分岐へ置換しないと落ちる（C4 地平上の救済漏れ）。
    #[test]
    fn visibility_max_below_with_c4_above_is_sunrise_eclipse() {
        assert_eq!(
            classify_visibility(true, Some(Degrees(-2.0)), Degrees(-0.05), Some(Degrees(3.0))),
            Visibility::SunriseEclipse,
            "max<0 but c4≥0 must yield SunriseEclipse (Sun rises between max and C4; later phase observable)"
        );
    }

    /// 受け入れ点（FIX 主題）: 最大食が地平直下だが C1 が地平上（C4 地平下）→ SunsetEclipse。
    /// 太陽が C1〜max の間に没し前半相のみ観測可能なので SunsetEclipse。現行コード
    /// （max<0 で即 BelowHorizon）は BelowHorizon を返すため red。撃破する変異: max<0 早期
    /// return（146 行）＋「C4 救済のみで C1 救済を欠く」部分修正を落とす（C1 地平上の救済）。
    #[test]
    fn visibility_max_below_with_c1_above_is_sunset_eclipse() {
        assert_eq!(
            classify_visibility(true, Some(Degrees(3.0)), Degrees(-0.05), Some(Degrees(-2.0))),
            Visibility::SunsetEclipse,
            "max<0 with c1≥0 and c4<0 must yield SunsetEclipse (Sun set between C1 and max; earlier phase observable)"
        );
    }

    /// 受け入れ点（FIX のガード）: max 地平下 ∧ 両接触とも地平下 → BelowHorizon（不変）。
    /// 食全体が地平下で観測不能。FIX が過剰救済しない（max<0 でも接触が両方下なら BelowHorizon を
    /// 維持する）ことを縛る。撃破する変異: 「max<0 ならば常に救済（無条件 Sunrise/Sunset）」へ
    /// 広げる過修正を落とす。
    #[test]
    fn visibility_max_below_with_both_contacts_below_is_below_horizon() {
        assert_eq!(
            classify_visibility(
                true,
                Some(Degrees(-3.0)),
                Degrees(-0.05),
                Some(Degrees(-2.0))
            ),
            Visibility::BelowHorizon,
            "max<0 with both contacts <0 must remain BelowHorizon (entire eclipse below horizon)"
        );
    }

    /// 受け入れ点（FIX のガード）: max 地平下 ∧ 両接触 None → BelowHorizon（真に観測不能）。
    /// 接触情報が無く max も地平下なら観測可能性を示す材料が無いので BelowHorizon。撃破する変異:
    /// 「max<0 で None を地平上扱いして救済」する過修正を落とす（`is_some_and` の None=false を縛る）。
    #[test]
    fn visibility_max_below_with_both_contacts_none_is_below_horizon() {
        assert_eq!(
            classify_visibility(true, None, Degrees(-5.0), None),
            Visibility::BelowHorizon,
            "max<0 with both contacts None must be BelowHorizon (no observable phase)"
        );
    }

    /// C1・C4 とも地平上（max も自動的に上）→ FullyVisible（全経過観測可能）。
    #[test]
    fn visibility_all_contacts_above_is_fully_visible() {
        assert_eq!(
            classify_visibility(
                true,
                Some(Degrees(10.0)),
                Degrees(40.0),
                Some(Degrees(15.0))
            ),
            Visibility::FullyVisible,
            "c1≥0 and c4≥0 (max≥0) must yield FullyVisible"
        );
    }

    /// 食開始 C1 が地平下（max≥0）→ SunriseEclipse（日の出中に食進行）。
    #[test]
    fn visibility_c1_below_is_sunrise_eclipse() {
        assert_eq!(
            classify_visibility(
                true,
                Some(Degrees(-3.0)),
                Degrees(20.0),
                Some(Degrees(35.0))
            ),
            Visibility::SunriseEclipse,
            "c1<0 with max≥0 and c4≥0 must yield SunriseEclipse"
        );
    }

    /// 食終了 C4 が地平下（max≥0, C1 地平上）→ SunsetEclipse（日没中に食終了）。
    #[test]
    fn visibility_c4_below_is_sunset_eclipse() {
        assert_eq!(
            classify_visibility(
                true,
                Some(Degrees(35.0)),
                Degrees(20.0),
                Some(Degrees(-3.0))
            ),
            Visibility::SunsetEclipse,
            "c4<0 with max≥0 and c1≥0 must yield SunsetEclipse"
        );
    }

    /// C1/C4 が窓外で None（max≥0, 一部の接触のみ地平上）→ PartialVisible。
    #[test]
    fn visibility_missing_contacts_is_partial_visible() {
        // c1=None, c4=Some≥0: 一部のみ可視 → PartialVisible。
        assert_eq!(
            classify_visibility(true, None, Degrees(10.0), Some(Degrees(5.0))),
            Visibility::PartialVisible,
            "c1=None with max≥0 must yield PartialVisible (partial coverage)"
        );
        // c4=None, c1=Some≥0。
        assert_eq!(
            classify_visibility(true, Some(Degrees(5.0)), Degrees(10.0), None),
            Visibility::PartialVisible,
            "c4=None with max≥0 must yield PartialVisible"
        );
        // 両方 None。
        assert_eq!(
            classify_visibility(true, None, Degrees(10.0), None),
            Visibility::PartialVisible,
            "both contacts None with max≥0 must yield PartialVisible"
        );
    }

    /// 境界（閾値=幾何高度 0° ちょうど・conventions §7 既定）: altitude=0 は「地平上（≥0）」側。
    /// max=0 ちょうどは BelowHorizon ではない（< 0 のみが地平下）。c1=0,c4=0 は FullyVisible。
    #[test]
    fn visibility_zero_altitude_boundary_is_on_horizon_side() {
        // max=0 ちょうど（c1,c4 ≥0）→ BelowHorizon ではなく FullyVisible（閾値は厳密 <0 で地平下）。
        assert_eq!(
            classify_visibility(true, Some(Degrees(0.0)), Degrees(0.0), Some(Degrees(0.0))),
            Visibility::FullyVisible,
            "altitude exactly 0° must count as above horizon (threshold is strict <0)"
        );
        // c1 = −ε（わずか地平下）→ SunriseEclipse（0 と −ε で分類が切り替わる境界）。
        assert_eq!(
            classify_visibility(true, Some(Degrees(-1e-9)), Degrees(5.0), Some(Degrees(5.0))),
            Visibility::SunriseEclipse,
            "c1 just below 0 must flip to SunriseEclipse (boundary at geometric 0°)"
        );
        // max = −ε（わずか地平下）だが C1・C4 とも地平上 → SunriseEclipse。
        // 修正前は「max<0 で即 BelowHorizon」だったが、接触が地平上なら食は観測可能。
        // 判定木は C4 地平上を先に見るため SunriseEclipse（後半相が観測可能）。
        assert_eq!(
            classify_visibility(true, Some(Degrees(1.0)), Degrees(-1e-9), Some(Degrees(1.0))),
            Visibility::SunriseEclipse,
            "max just below 0 with both contacts above must be SunriseEclipse (c4≥0 checked first; not BelowHorizon)"
        );
    }

    /// 判定木の優先順位（S2）: in_eclipse ∧ max≥0 ∧ c1<0 ∧ c4<0（日の出帯食と日没帯食が
    /// 同時成立する exotic ケース）。契約の判定木は in_eclipse→max<0→(c1≥0&c4≥0:Full)→
    /// (c1<0:Sunrise)→(c4<0:Sunset)→残り Partial の順で、c1<0 を c4<0 より先に判定するため
    /// SunriseEclipse になる。判定木の取りこぼし穴（c1/c4 同時地平下）を閉じる。
    #[test]
    fn visibility_both_c1_c4_below_with_max_above_is_sunrise_eclipse() {
        assert_eq!(
            classify_visibility(
                true,
                Some(Degrees(-3.0)),
                Degrees(10.0),
                Some(Degrees(-2.0))
            ),
            Visibility::SunriseEclipse,
            "c1<0 AND c4<0 with max≥0 must yield SunriseEclipse (c1<0 evaluated before c4<0)"
        );
    }

    /// Sunrise 境界の左端（S3）: c1 がちょうど 0°（地平上側）、max≥0、c4≥0 → FullyVisible。
    /// c1=0 は地平上（境界は厳密 <0 で Sunrise）であることを単独で明示し、Sunrise(c1<0) の
    /// 左端境界を固定する。
    #[test]
    fn visibility_c1_exactly_zero_is_fully_visible() {
        assert_eq!(
            classify_visibility(true, Some(Degrees(0.0)), Degrees(5.0), Some(Degrees(5.0))),
            Visibility::FullyVisible,
            "c1 exactly 0° (on horizon, not below) with max≥0 and c4≥0 must yield FullyVisible"
        );
    }

    /// 境界（S4・c1=0 ちょうど × c4 地平下）: c1=0°（地平上側）, max≥0, c4<0 → SunsetEclipse。
    /// 地平閾値 0° は地平上側（`< 0` のみ below）を境界で固定し、`< → <=` 変異を撃破。
    /// 正しい `<`: c1_below=(0<0)=false, c1_above=(0≥0)=true, c4_above=false → not Full,
    /// c1_below=false, c4_below=(−5<0)=true → SunsetEclipse。
    /// 変異 `<=`（148行）: c1_below=(0<=0)=true → SunriseEclipse になり assert が落ちる＝撃破。
    #[test]
    fn visibility_c1_exactly_zero_with_c4_below_is_sunset_eclipse() {
        assert_eq!(
            classify_visibility(true, Some(Degrees(0.0)), Degrees(5.0), Some(Degrees(-5.0))),
            Visibility::SunsetEclipse,
            "c1 exactly 0° (on horizon, not below) with c4<0 must yield SunsetEclipse"
        );
    }

    /// 境界（S5・c1=None × c4=0 ちょうど）: c1=None, max≥0, c4=0°（地平上側）→ PartialVisible。
    /// 地平閾値 0° は地平上側（`< 0` のみ below）を境界で固定し、`< → <=` 変異を撃破。
    /// 正しい `<`: c1_above=false(None), c4_above=(0≥0)=true → not Full(c1_above false),
    /// c1_below=false(None), c4_below=(0<0)=false → else → PartialVisible。
    /// 変異 `<=`（149行）: c4_below=(0<=0)=true → SunsetEclipse になり assert が落ちる＝撃破。
    #[test]
    fn visibility_c1_none_with_c4_exactly_zero_is_partial_visible() {
        assert_eq!(
            classify_visibility(true, None, Degrees(5.0), Some(Degrees(0.0))),
            Visibility::PartialVisible,
            "c1=None with c4 exactly 0° (on horizon, not below) must yield PartialVisible"
        );
    }

    /// 6 値すべてが少なくとも一度は生成されることのメタ確認（網羅性の自己点検）。
    /// 個別テストの期待値表をまとめ、Visibility の全バリアントが出ることを保証する。
    #[test]
    fn visibility_all_six_variants_are_reachable() {
        use Visibility::*;
        let table: &[(bool, Option<Degrees>, Degrees, Option<Degrees>, Visibility)] = &[
            (
                false,
                Some(Degrees(30.0)),
                Degrees(45.0),
                Some(Degrees(20.0)),
                NotVisible,
            ),
            (
                true,
                Some(Degrees(-20.0)),
                Degrees(-5.0),
                Some(Degrees(-1.0)),
                BelowHorizon,
            ),
            (
                true,
                Some(Degrees(10.0)),
                Degrees(40.0),
                Some(Degrees(15.0)),
                FullyVisible,
            ),
            (
                true,
                Some(Degrees(-3.0)),
                Degrees(20.0),
                Some(Degrees(35.0)),
                SunriseEclipse,
            ),
            (
                true,
                Some(Degrees(35.0)),
                Degrees(20.0),
                Some(Degrees(-3.0)),
                SunsetEclipse,
            ),
            (
                true,
                None,
                Degrees(10.0),
                Some(Degrees(5.0)),
                PartialVisible,
            ),
        ];
        // Visibility は Hash を導出しないため Vec で集める（公開型は変更しない）。
        let all = [
            NotVisible,
            BelowHorizon,
            FullyVisible,
            SunriseEclipse,
            SunsetEclipse,
            PartialVisible,
        ];
        let mut produced: Vec<Visibility> = Vec::new();
        for &(in_e, c1, mx, c4, want) in table {
            let got = classify_visibility(in_e, c1, mx, c4);
            assert_eq!(got, want, "table row mismatch: in_eclipse={in_e}");
            if !produced.contains(&got) {
                produced.push(got);
            }
        }
        for v in all {
            assert!(
                produced.contains(&v),
                "Visibility variant {v:?} was not produced by the fixture table"
            );
        }
        assert_eq!(
            produced.len(),
            6,
            "all six Visibility variants must be covered"
        );
    }
}
