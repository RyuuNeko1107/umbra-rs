//! 解析的太陽暦（VSOP87D・地球日心→太陽地心の反転）。ISSUE-013。
//!
//! ISSUE-033 が生成・コミットした packed 係数（`generated/vsop87/vsop87d_earth.bin`）を
//! `include_bytes!` で取り込み、地球の日心黄道座標 (L, B, R)（VSOP87D, 平均分点 of date）を
//! 評価する。太陽の地心位置はこれを符号反転して得る（補正前の幾何的地心位置）。光行時間・
//! 歳差章動・光行差は ISSUE-015、黄道 of date→ICRS は ISSUE-035。
//!
//! 評価式: `s = Σ_{α} T^α · Σ_k A_{α,k}·cos(B_{α,k} + C_{α,k}·T)`、T = ユリウス千年 from
//! J2000 TDB = `(JD_TDB − 2451545)/365250`（Bretagnon & Francou 1988）。

use std::sync::OnceLock;
use umbra_core::constants::{ASTRONOMICAL_UNIT_KM, JULIAN_MILLENNIUM_DAYS};
use umbra_core::{Radians, TdbInstant, Vector3};

/// packed VSOP87D 地球係数（ISSUE-033 生成物、flat little-endian f64）。
/// レイアウト: `[n_sections, <各セクション = [variable, power, n_terms, <各項 amp,phase,freq>]>...]`。
const PACKED: &[u8] = include_bytes!("../../../generated/vsop87/vsop87d_earth.bin");

/// VSOP87 級数 1 項: `amplitude·cos(phase + frequency·T)`。
struct Term {
    amplitude: f64,
    phase: f64,
    frequency: f64,
}

/// 1 セクション（変数 v × べき α）。
struct Section {
    /// 変数（1=L, 2=B, 3=R）。
    variable: u8,
    /// T のべき（0..=5）。
    power: u8,
    terms: Vec<Term>,
}

/// 検証済み f64 を非負カウントへ（自前生成・verify-generated ゲート済みのため信頼）。
#[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
fn nint(value: f64) -> usize {
    value.round() as usize
}

/// packed を 1 度だけ復号する。
fn model() -> &'static [Section] {
    static MODEL: OnceLock<Vec<Section>> = OnceLock::new();
    MODEL.get_or_init(|| {
        let values: Vec<f64> = PACKED
            .chunks_exact(8)
            .map(|chunk| f64::from_le_bytes(chunk.try_into().expect("packed length multiple of 8")))
            .collect();
        let n_sections = nint(values[0]);
        let mut idx = 1;
        let mut sections = Vec::with_capacity(n_sections);
        for _ in 0..n_sections {
            let variable = u8::try_from(nint(values[idx])).expect("variable fits u8");
            let power = u8::try_from(nint(values[idx + 1])).expect("power fits u8");
            let n_terms = nint(values[idx + 2]);
            idx += 3;
            let mut terms = Vec::with_capacity(n_terms);
            for _ in 0..n_terms {
                terms.push(Term {
                    amplitude: values[idx],
                    phase: values[idx + 1],
                    frequency: values[idx + 2],
                });
                idx += 3;
            }
            sections.push(Section {
                variable,
                power,
                terms,
            });
        }
        sections
    })
}

/// 地球の日心黄道座標 `(L, B, R)`（VSOP87D, 平均分点 of date）。
/// L は黄経 \[0,2π) rad、B は黄緯 rad、R は動径 AU。TDB 引数。
pub fn earth_heliocentric_lbr(time_tdb: TdbInstant) -> (Radians, Radians, f64) {
    let t = time_tdb.jd2().julian_millennia_since_j2000();
    // 変数 1=L, 2=B, 3=R を index 0,1,2 に集計。
    let mut lbr = [0.0_f64; 3];
    for section in model() {
        let series_sum: f64 = section
            .terms
            .iter()
            .map(|term| term.amplitude * (term.phase + term.frequency * t).cos())
            .sum();
        lbr[usize::from(section.variable - 1)] += t.powi(i32::from(section.power)) * series_sum;
    }
    (
        Radians::new(lbr[0]).normalized_two_pi(),
        Radians::new(lbr[1]),
        lbr[2],
    )
}

/// 太陽の地心位置（黄道・平均分点 of date の直交座標、AU）= −(地球日心直交座標)。
/// 補正前の幾何的地心位置（光行時間・光行差・章動は別 Issue）。TDB 引数。
pub fn sun_geocentric_ecliptic_of_date(time_tdb: TdbInstant) -> Vector3 {
    let (l, b, r) = earth_heliocentric_lbr(time_tdb);
    let (sin_l, cos_l) = l.0.sin_cos();
    let (sin_b, cos_b) = b.0.sin_cos();
    // 地球日心直交（黄道 of date）→ 符号反転で太陽地心。
    Vector3::new(-(r * cos_b * cos_l), -(r * cos_b * sin_l), -(r * sin_b))
}

/// 地球の日心速度（黄道 of date 直交, km/s）。VSOP87D 級数の**項別解析微分**
/// （numerical-policy §A2(1): `d/dT[A·cos(B+C·T)] = −A·C·sin(B+C·T)`、べき項は積の微分）。
/// 観測者（地心）速度 ≈ 地球日心速度として光行差・光行時間に用いる（ISSUE-015）。TDB 引数。
///
/// `s = Σ_α T^α·Σ_k A·cos(B+C·T)` の `ds/dT = Σ_α [α·T^(α−1)·Σ A cos + T^α·Σ A(−C)sin]`
/// を L,B,R で求め、球面→直交速度へ展開（AU/千年 → km/s に換算）。
pub fn earth_heliocentric_velocity_ecliptic_of_date(time_tdb: TdbInstant) -> Vector3 {
    let t = time_tdb.jd2().julian_millennia_since_j2000();
    // 変数 1=L,2=B,3=R を index 0,1,2 に、値と dT 微分（千年あたり）を集計。
    let mut val = [0.0_f64; 3];
    let mut dval = [0.0_f64; 3];
    for section in model() {
        // S = Σ A cos(B+C·t)、dS/dt = Σ A·(−C)·sin(B+C·t)。
        let mut s = 0.0;
        let mut ds = 0.0;
        for term in &section.terms {
            let arg = term.phase + term.frequency * t;
            s += term.amplitude * arg.cos();
            ds += term.amplitude * (-term.frequency) * arg.sin();
        }
        let idx = usize::from(section.variable - 1);
        let alpha = i32::from(section.power);
        val[idx] += t.powi(alpha) * s;
        // d/dT[t^α·S] = (α≥1 ? α·t^(α−1)·S : 0) + t^α·dS。
        let power_term = if alpha >= 1 {
            f64::from(alpha) * t.powi(alpha - 1) * s
        } else {
            0.0
        };
        dval[idx] += power_term + t.powi(alpha) * ds;
    }
    // 球面→黄道直交速度 [AU/千年] → km/s（1 千年 = 365250 日 × 86400 s）。
    let factor = ASTRONOMICAL_UNIT_KM / (JULIAN_MILLENNIUM_DAYS * 86_400.0);
    ecliptic_velocity_from_spherical(val[0], val[1], val[2], dval[0], dval[1], dval[2])
        .scale(factor)
}

/// 球面座標 `(l, b, r)` と各時間微分 `(dl, db, dr)` から黄道直交速度を構成する（純関数）。
/// 位置 `(r·cosb·cosl, r·cosb·sinl, r·sinb)` の時間微分。返り値の単位は r/時間（呼び出し側で換算）。
fn ecliptic_velocity_from_spherical(l: f64, b: f64, r: f64, dl: f64, db: f64, dr: f64) -> Vector3 {
    let (sin_l, cos_l) = l.sin_cos();
    let (sin_b, cos_b) = b.sin_cos();
    Vector3::new(
        dr * cos_b * cos_l - r * sin_b * cos_l * db - r * cos_b * sin_l * dl,
        dr * cos_b * sin_l - r * sin_b * sin_l * db + r * cos_b * cos_l * dl,
        dr * sin_b + r * cos_b * db,
    )
}

#[cfg(test)]
mod tests {
    // VSOP87 配布の検証値（vsop87.chk）は f64 表現可能桁数を超える桁で配布される。
    // provenance（一次ソースの逐語転記）を保つため桁を削らず、ここに限り過剰精度リント
    // （余剰桁は最近接 f64 へ丸められ値は不変）を許可する。nutation.rs と同様。
    #![allow(clippy::excessive_precision)]

    use super::*;
    use umbra_core::JulianDate2;

    // ------------------------------------------------------------------
    // 一次オラクル: VSOP87 配布の自己検証ファイル
    //   `data/coefficient-source/vsop87/vsop87.chk`（コミット済み）の
    //   `VSOP87D EARTH JD<...> ... l <L> rad  b <B> rad  r <R> au` 行。
    //   これはモデル（VSOP87D 級数）自身の独立参照値であり、本実装（packed 係数の
    //   f64 総和）とは別経路（配布元が公表した期待値）。DE 不要・M2 暫定ゲート。
    //
    //   chk の l は [0,2π) 正規化値。`earth_heliocentric_lbr` の L も [0,2π) 正規化で返す。
    //   各 JD（TDB, JD...0 = 当該日 12h TDB）と転記値:
    //     JD 2451545.0 (2000-01-01 12h): l=1.7519238681  b=-0.0000039656  r=0.9833276819
    //     JD 2415020.0 (1899-12-31 12h): l=1.7391225563  b=-0.0000005679  r=0.9832689778
    //     JD 2378495.0 (1799-12-30 12h): l=1.7262638916  b= 0.0000002083  r=0.9832274321
    //     JD 2341970.0 (1699-12-29 12h): l=1.7134419105  b= 0.0000025051  r=0.9831498441
    //   （chk ファイル該当行は vsop87.chk:1401-1414。値は逐語転記。）
    // ------------------------------------------------------------------

    /// chk テストエポックの `(JD_TDB, l_rad, b_rad, r_au)`。
    const CHK_CASES: [(f64, f64, f64, f64); 4] = [
        (2451545.0, 1.7519238681, -0.0000039656, 0.9833276819),
        (2415020.0, 1.7391225563, -0.0000005679, 0.9832689778),
        (2378495.0, 1.7262638916, 0.0000002083, 0.9832274321),
        (2341970.0, 1.7134419105, 0.0000025051, 0.9831498441),
    ];

    // 突合許容 [rad / au]。
    //
    // chk は 10 桁（~1e-10 rad/au 丸め）で配布される。本実装は VSOP87D 全項を f64 で
    // 総和するため評価誤差は ~1e-12。1e-9 は chk 丸め + 総和順差を吸収しつつ、係数読取り・
    // 評価式・時刻スケール（T = ユリウス千年）の実バグ（≫1e-6 rad/au）を確実に検出する。
    const TOL: f64 = 1e-9;

    /// TDB の単一要素 JD から `TdbInstant` を構築（part2=0）。
    fn tdb(jd: f64) -> TdbInstant {
        TdbInstant::from_jd2(JulianDate2::new(jd, 0.0))
    }

    /// 許容つきスカラ比較（nutation.rs の `close` 踏襲）。
    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    // ==================================================================
    // 1. earth_heliocentric_lbr vs vsop87.chk（多日付）
    //    J2000 / 1900 / 1800 / 1700 の 4 エポックで L,B,R を chk 値と TOL 以内照合。
    //    1700–1900 を含めることで T のべき（T^1..T^5）と periodic 項の時間依存を励起し、
    //    評価式・打切り・時刻スケール（千年 vs 世紀等）の誤りを検出する。
    // ==================================================================
    #[test]
    fn earth_lbr_matches_vsop87_chk() {
        for (jd, exp_l, exp_b, exp_r) in CHK_CASES {
            let (l, b, r) = earth_heliocentric_lbr(tdb(jd));
            assert!(
                close(l.0, exp_l, TOL),
                "L at JD{jd}: {} vs chk {exp_l} (|diff| = {:e} > TOL {:e})",
                l.0,
                (l.0 - exp_l).abs(),
                TOL
            );
            assert!(
                close(b.0, exp_b, TOL),
                "B at JD{jd}: {} vs chk {exp_b} (|diff| = {:e} > TOL {:e})",
                b.0,
                (b.0 - exp_b).abs(),
                TOL
            );
            assert!(
                close(r, exp_r, TOL),
                "R at JD{jd}: {} vs chk {exp_r} (|diff| = {:e} > TOL {:e})",
                r,
                (r - exp_r).abs(),
                TOL
            );
        }
    }

    // ==================================================================
    // 2. L の正規化: 返る黄経 L が [0,2π) に入る（複数エポックで）。
    //    正規化漏れ・符号誤りを検出する。
    // ==================================================================
    #[test]
    fn earth_longitude_is_normalized_to_two_pi() {
        use core::f64::consts::TAU;
        for (jd, ..) in CHK_CASES {
            let (l, ..) = earth_heliocentric_lbr(tdb(jd));
            assert!(
                (0.0..TAU).contains(&l.0),
                "L at JD{jd} not in [0,2π): {}",
                l.0
            );
        }
    }

    // ==================================================================
    // 3. オーダーサニティ: R は地球-太陽距離 ~0.98..1.02 au、黄緯 |B| は小（< 1e-4 rad）。
    //    定数返却・単位取り違え（千年 vs au 等）・スケール暴走を粗く検出する。
    // ==================================================================
    #[test]
    fn earth_radius_and_latitude_order_of_magnitude() {
        for (jd, ..) in CHK_CASES {
            let (_l, b, r) = earth_heliocentric_lbr(tdb(jd));
            assert!(
                (0.98..1.02).contains(&r),
                "R at JD{jd} out of Earth-Sun range: {r}"
            );
            assert!(b.0.abs() < 1e-4, "|B| at JD{jd} too large: {}", b.0.abs());
        }
    }

    // ==================================================================
    // 4. 太陽地心 = 地球日心直交の符号反転
    //    earth 直交: x=R·cosB·cosL, y=R·cosB·sinL, z=R·sinB（黄道 of date, au）。
    //    sun_geocentric_ecliptic_of_date(t) == (-x,-y,-z) を各成分 ~1e-12 au で確認。
    //    かつ earth_rect + sun_geo ≈ 0、norm(sun_geo) == R を ~1e-12 au で確認。
    //    同一実装由来のため厳密近く（1e-12）で一致するはず。
    // ==================================================================
    #[test]
    fn sun_geocentric_is_negated_earth_heliocentric_rectangular() {
        // 直交化は L,B,R から純粋に幾何で導く独立計算（評価器とは別ロジック）。
        const RECT_TOL: f64 = 1e-12;
        for (jd, ..) in CHK_CASES {
            let t = tdb(jd);
            let (l, b, r) = earth_heliocentric_lbr(t);
            let earth_rect = Vector3::new(
                r * b.0.cos() * l.0.cos(),
                r * b.0.cos() * l.0.sin(),
                r * b.0.sin(),
            );
            let sun = sun_geocentric_ecliptic_of_date(t);

            // 各成分が earth 直交の符号反転。
            assert!(
                close(sun.x, -earth_rect.x, RECT_TOL),
                "sun.x at JD{jd}: {} vs {} (|diff| = {:e})",
                sun.x,
                -earth_rect.x,
                (sun.x + earth_rect.x).abs()
            );
            assert!(
                close(sun.y, -earth_rect.y, RECT_TOL),
                "sun.y at JD{jd}: {} vs {} (|diff| = {:e})",
                sun.y,
                -earth_rect.y,
                (sun.y + earth_rect.y).abs()
            );
            assert!(
                close(sun.z, -earth_rect.z, RECT_TOL),
                "sun.z at JD{jd}: {} vs {} (|diff| = {:e})",
                sun.z,
                -earth_rect.z,
                (sun.z + earth_rect.z).abs()
            );

            // 反転の整合: earth_rect + sun ≈ 0。
            let sum = earth_rect + sun;
            assert!(
                sum.norm() < RECT_TOL,
                "earth_rect + sun_geo not ~0 at JD{jd}: |sum| = {:e}",
                sum.norm()
            );

            // sun_geo の大きさは R に等しい。
            assert!(
                close(sun.norm(), r, RECT_TOL),
                "|sun_geo| at JD{jd}: {} vs R {r} (|diff| = {:e})",
                sun.norm(),
                (sun.norm() - r).abs()
            );
        }
    }

    // ==================================================================
    // ISSUE-015 prereq: earth_heliocentric_velocity_ecliptic_of_date
    //   地球の日心速度（黄道 of date 直交, km/s）の解析微分を検証する。
    //   一次オラクル: sun_geocentric_ecliptic_of_date は太陽地心 = −(地球日心位置)
    //   なので、地球日心速度 = −d/dt[sun_geo]。これを sun_geo の**中心差分**で
    //   近似し、解析速度と照合する（評価器とは独立な数値微分経路）。
    // ==================================================================

    /// 速度テスト用エポック（TDB JD, part2=0）。J2000 / 近日点近傍 / 任意の近代日付。
    const VEL_JD: [f64; 3] = [2451545.0, 2444239.5, 2469807.0];
    /// 中心差分の半ステップ [日]。
    const H_DAYS: f64 = 0.5;
    /// 中心差分突合の許容 [km/s]。h=0.5 日の O(h²) 打切り誤差は (h²/6)·|x'''|≈3.7e-4 km/s
    /// （地球公転 v≈29.8, ω=2π/年, |x'''|≈v·ω²）。近日点近傍の局所躍度増を見て 1e-3 を採る。
    /// 項取りこぼし・符号・スケール誤りは各成分を ≫0.01 km/s ずらすため 1e-3 でも確実に検出。
    const VEL_TOL_KM_S: f64 = 1e-3;

    /// 要素ごと近接（clippy::float_cmp 回避）。
    fn vec_close(a: Vector3, b: Vector3, tol: f64) -> bool {
        (a.x - b.x).abs() < tol && (a.y - b.y).abs() < tol && (a.z - b.z).abs() < tol
    }

    /// `jd`（part1）に part2=`offset_days` を載せた TdbInstant。
    fn tdb_offset(jd: f64, offset_days: f64) -> TdbInstant {
        TdbInstant::from_jd2(JulianDate2::new(jd, offset_days))
    }

    /// sun_geo の中心差分から得た地球日心速度 [km/s]（= −d/dt[sun_geo]）。
    fn earth_velocity_central_difference(jd: f64, h_days: f64) -> Vector3 {
        let plus = sun_geocentric_ecliptic_of_date(tdb_offset(jd, h_days));
        let minus = sun_geocentric_ecliptic_of_date(tdb_offset(jd, -h_days));
        let d_sun_au_per_day = (plus - minus).scale(1.0 / (2.0 * h_days));
        d_sun_au_per_day.scale(-ASTRONOMICAL_UNIT_KM / 86_400.0)
    }

    /// 1. 解析速度 vs sun_geo の中心差分（複数エポック・各成分）。
    #[test]
    fn earth_velocity_matches_central_difference_of_sun_geocentric() {
        for jd in VEL_JD {
            let analytic = earth_heliocentric_velocity_ecliptic_of_date(tdb(jd));
            let numeric = earth_velocity_central_difference(jd, H_DAYS);
            assert!(
                vec_close(analytic, numeric, VEL_TOL_KM_S),
                "vel at JD{jd}: analytic ({},{},{}) vs central-diff ({},{},{})",
                analytic.x,
                analytic.y,
                analytic.z,
                numeric.x,
                numeric.y,
                numeric.z
            );
        }
    }

    /// 2. 速度の大きさ: 地球公転速度 遠日点~29.3 / 近日点~30.3 km/s。
    #[test]
    fn earth_speed_is_orbital_magnitude() {
        for jd in VEL_JD {
            let speed = earth_heliocentric_velocity_ecliptic_of_date(tdb(jd)).norm();
            assert!(
                (29.0..30.5).contains(&speed),
                "Earth speed at JD{jd} out of orbital range: {speed} km/s"
            );
        }
    }

    /// 3. 速度と日心位置の近直交性（軌道はほぼ円 e≈0.0167）。|cos∠(v,r)| < 0.05。
    #[test]
    fn earth_velocity_nearly_perpendicular_to_radius() {
        for jd in [VEL_JD[0], VEL_JD[1]] {
            let t = tdb(jd);
            let vel = earth_heliocentric_velocity_ecliptic_of_date(t);
            // 地球日心位置 [km] = 太陽地心の符号反転 ×AU。
            let pos = sun_geocentric_ecliptic_of_date(t).scale(-ASTRONOMICAL_UNIT_KM);
            let cos_angle = vel.dot(pos) / (vel.norm() * pos.norm());
            assert!(
                cos_angle.abs() < 0.05,
                "|cos∠(v,r)| at JD{jd} too large (orbit ~circular): {}",
                cos_angle.abs()
            );
        }
    }

    /// 4. 球面→直交速度の純関数 `ecliptic_velocity_from_spherical` を、**非特殊な合成軌道**
    ///    （緯度 b を 0 から離す＝全項を励起）の中心差分と照合。地球は b≈0 で緯度結合項が
    ///    休眠し XYZ テストでは隠れるため、一般入力で各項（sin_b, cos_b, db 結合）を直接検証する。
    #[test]
    fn ecliptic_velocity_from_spherical_matches_finite_difference() {
        // 全成分が distinct・非特殊（b0=0.4 で cos_b≈0.92, sin_b≈0.39, db≠0）。
        let (l0, b0, r0, dl, db, dr) = (0.7_f64, 0.4_f64, 2.3_f64, 1.1_f64, 0.6_f64, 0.3_f64);
        let pos = |t: f64| {
            let (l, b, r) = (l0 + dl * t, b0 + db * t, r0 + dr * t);
            Vector3::new(r * b.cos() * l.cos(), r * b.cos() * l.sin(), r * b.sin())
        };
        let h = 1e-6;
        let fd = (pos(h) - pos(-h)).scale(1.0 / (2.0 * h));
        let v = ecliptic_velocity_from_spherical(l0, b0, r0, dl, db, dr);
        assert!(
            vec_close(v, fd, 1e-6),
            "spherical velocity ({},{},{}) vs finite-diff ({},{},{})",
            v.x,
            v.y,
            v.z,
            fd.x,
            fd.y,
            fd.z
        );
    }
}
