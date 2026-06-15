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
}
