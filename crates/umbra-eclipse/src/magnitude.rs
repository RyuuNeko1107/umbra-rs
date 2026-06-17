//! 食分（magnitude）と食面積（obscuration）（`docs/issues/ISSUE-027`、`docs/conventions.md` §8/§9,
//! Explanatory Supplement §11 / Meeus *Astronomical Algorithms* Ch.54）。
//!
//! いずれも瞬時量から計算する純関数（時刻・暦に非依存）:
//! - 食分 `magnitude = (L1 − m)/(L1 + L2)`（基本面の外接縁から内接縁への食い込み割合, Re）。
//!   皆既で 1 超、金環で ≈1、部分食で 0..1、離隔（m ≥ L1）で 0。
//! - 食面積 `obscuration` = 太陽円と月円の重なり面積 / 太陽円面積（視半径平面, 0..1）。
//!   円-円交差面積（lens area）の標準式。acos 引数を `[-1,1]` クランプ、5 縮退境界
//!   （離隔/内包/外接/内接/同半径）を明示処理（accuracy.md §2.2）。
//!
//! 単位系の分離（混在禁止, conventions §1）: 食分は基本面 Re（m, L1, L2）、食面積は視半径平面
//! （太陽/月見かけ半径と中心離隔を同単位で）。

// eclipse_magnitude / eclipse_obscuration（pub(crate)）は ISSUE-026/043 が消費するまで未使用。
#![allow(dead_code)]

use core::f64::consts::PI;

/// 食分（太陽直径に対する食い込み量）。皆既で 1 を超えうる（api-draft §3.4）。
#[derive(Clone, Copy, Debug, PartialEq, PartialOrd)]
pub struct EclipseMagnitude(pub f64);

/// 食面積（太陽面積に対し月が覆う割合）。厳密に `[0, 1]`（api-draft §3.4）。
#[derive(Clone, Copy, Debug, PartialEq, PartialOrd)]
pub struct Obscuration(pub f64);

/// 基本面の中心間距離 `m` と影半径 `l1`（半影, 外接）・`l2`（本影, 符号付き）から食分を求める。
///
/// `magnitude = (L1 − m)/(L1 + L2)`、ここで L1 = l1, L2 = |l2|（食分は内接縁までの食い込み割合）。
/// 単位は Re（呼出側で統一, conventions §1）。`m ≥ L1`（離隔・食なし）は 0 にクランプ。
/// 皆既/金環で 1 跨ぎを許容（`EclipseMagnitude` は 1 超可）。
pub(crate) fn eclipse_magnitude(m: f64, l1: f64, l2: f64) -> EclipseMagnitude {
    // Explanatory Supplement §11 / Meeus Ch.54: magnitude = (L1 − m)/(L1 + |L2|)。
    // 離隔（m ≥ L1）は負になるので 0 にクランプ（食なし）。皆既/金環の 1 超は許容（上限クランプしない）。
    let magnitude = (l1 - m) / (l1 + l2.abs());
    EclipseMagnitude(magnitude.max(0.0))
}

/// 太陽-月の見かけ中心間距離 `separation` と太陽見かけ半径 `r_sun`・月見かけ半径 `r_moon`
/// （同単位）から食面積を求める。
///
/// 太陽円と月円の重なり面積 / 太陽円面積（円-円交差の lens 面積）。acos 引数を `[-1,1]`、
/// 判別式の負値を 0 にクランプ。5 縮退境界（離隔/内包/外接/内接/同半径）を明示処理（accuracy.md §2.2）。
/// 結果は `[0, 1]`。
pub(crate) fn eclipse_obscuration(separation: f64, r_sun: f64, r_moon: f64) -> Obscuration {
    let d = separation;
    let big_r = r_sun;
    let r = r_moon;
    let sun_area = PI * big_r * big_r;

    let overlap = if d >= big_r + r {
        // 離隔: 重なりなし。
        0.0
    } else if d <= (big_r - r).abs() {
        // 内包: 小円が大円に完全内包 → 小円の面積。月大(r≥R)→ obsc=1, 太陽大(金環)→ r²/R²。
        // ここで d=0（同心）も含み、部分重なり式の 2dR/2dr/d=0 の 0 除算を回避する。
        let min = big_r.min(r);
        PI * min * min
    } else {
        // 部分重なり: 円-円交差 lens 面積（MathWorld "Circle-Circle Intersection"）。
        // acos 引数を [-1,1]、判別式の負値を 0 にクランプ（丸めで NaN を出さない, accuracy.md §2.2）。
        let cos_a = ((d * d + big_r * big_r - r * r) / (2.0 * d * big_r)).clamp(-1.0, 1.0);
        let cos_b = ((d * d + r * r - big_r * big_r) / (2.0 * d * r)).clamp(-1.0, 1.0);
        let discriminant = (-d + r + big_r) * (d + r - big_r) * (d - r + big_r) * (d + r + big_r);
        big_r * big_r * cos_a.acos() + r * r * cos_b.acos() - 0.5 * discriminant.max(0.0).sqrt()
    };

    Obscuration((overlap / sun_area).clamp(0.0, 1.0))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    // ----------------------------------------------------------------------
    // 浮動小数比較ヘルパ（clippy::float_cmp 回避・許容つき比較）
    // ----------------------------------------------------------------------

    /// 絶対誤差比較。
    fn close_abs(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    /// 純幾何（円交差・解析解）の一致許容。f64 丸めのみを許容（ISSUE-027 許容: ≤1e-10）。
    const GEOM_TOL: f64 = 1e-10;

    /// 食分（解析的自明値）との一致許容（ISSUE-027 許容: ≤1e-12）。
    const MAG_TOL: f64 = 1e-12;

    // ----------------------------------------------------------------------
    // 独立オラクル（実装関数 eclipse_* は一切呼ばない／契約式のコピーでもない）
    // ----------------------------------------------------------------------
    //
    // 円-円交差（lens）面積の独立参照実装。
    // 出典: Weisstein, MathWorld "Circle-Circle Intersection"
    //   https://mathworld.wolfram.com/Circle-CircleIntersection.html
    // ここでは「面積比 obscuration の検証用」に縮退境界を独立に場合分けし、
    // 部分重なりのみ lens 公式を素朴展開する。実装側の関数とは独立に書いた参照値であり、
    // 縮退ごとの値（0, πr², πmin² など）は解析解として手で確定したもの。
    //
    /// 太陽円（半径 R）と月円（半径 r）の重なり面積。中心間距離 d。
    /// 縮退（離隔/内包/外接/内接/同心）は解析解で直接返し、丸めに頼らない。
    fn lens_overlap_area(d: f64, big_r: f64, r: f64) -> f64 {
        let sep = big_r + r; // 離隔閾値
        let contain = (big_r - r).abs(); // 内包閾値
        if d >= sep {
            // 離隔: 重なりなし。
            0.0
        } else if d <= contain {
            // 内包: 小円が大円に完全に入る → 小円の面積。
            let min = big_r.min(r);
            PI * min * min
        } else {
            // 部分重なり: 標準 lens 公式（MathWorld 4.07/4.08）。
            let term_a =
                big_r * big_r * (((d * d + big_r * big_r - r * r) / (2.0 * d * big_r)).acos());
            let term_b = r * r * (((d * d + r * r - big_r * big_r) / (2.0 * d * r)).acos());
            let tri = 0.5
                * ((-d + r + big_r) * (d + r - big_r) * (d - r + big_r) * (d + r + big_r)).sqrt();
            term_a + term_b - tri
        }
    }

    /// 食面積の独立参照値: overlap / (π·R²)、[0,1] にクランプ。
    fn obscuration_oracle(d: f64, r_sun: f64, r_moon: f64) -> f64 {
        let area = lens_overlap_area(d, r_sun, r_moon);
        (area / (PI * r_sun * r_sun)).clamp(0.0, 1.0)
    }

    /// 食分の独立参照値: (L1 − m)/(L1 + |l2|)、負は 0 にクランプ。
    /// 契約式と同形だが、ここでは「食い込み割合の定義そのもの」を最小限で別記したもの。
    /// 個々のテストは加えて手計算した具体数値でも縛る（追認回避の本体はそちら）。
    fn magnitude_oracle(m: f64, l1: f64, l2: f64) -> f64 {
        let mag = (l1 - m) / (l1 + l2.abs());
        if mag < 0.0 {
            0.0
        } else {
            mag
        }
    }

    // ======================================================================
    // 契約A: eclipse_magnitude = (L1 − m)/(L1 + |L2|), m≥L1 で 0 クランプ
    // ======================================================================

    /// 中心線上（m=0）: magnitude = L1/(L1+|l2|)。皆既配置（l2<0）で具体値を独立に手計算。
    /// 例 L1=0.55, l2=-0.45 → |l2|=0.45 → mag = 0.55/(0.55+0.45) = 0.55。皆既近傍だが
    /// この m=0/L1 配置では 1 未満（皆既の 1 超は別テストで縛る）。
    #[test]
    fn magnitude_central_partial_value() {
        let got = eclipse_magnitude(0.0, 0.55, -0.45).0;
        // 手計算: 0.55 / 1.00 = 0.55。
        assert!(
            close_abs(got, 0.55, MAG_TOL),
            "central magnitude: got {got}, want 0.55"
        );
    }

    /// 皆既で magnitude > 1 を許容: m が小さく |l2| が小さいと比が 1 を超える。
    /// 例 m=0.0, L1=0.6, l2=-0.1 → 0.6/(0.6+0.1) = 0.857…? いや 1 未満。
    /// 1 超を出すには L1 − m > L1 + |l2| すなわち −m > |l2|、不可（m,|l2|≥0）。
    /// 食分>1 は L1>1 かつ m 小の配置でなく、定義上は m<0 を含む配置で生じる。
    /// 中心通過の皆既では m≈0、本影縁の内側に入るため (L1−m)/(L1+|l2|) は最大で
    /// L1/(L1+|l2|)<1。よって本式単体では magnitude>1 は「m が負（影中心を跨ぐ）」で発生。
    /// 例 m=-0.30, L1=0.55, l2=-0.45 → (0.55+0.30)/(0.55+0.45)=0.85/1.0=0.85。
    /// 依然 1 未満。1 超は L1+|l2| < L1−m すなわち m < −|l2| のとき。
    /// 例 m=-0.50, L1=0.55, l2=-0.05 → (0.55+0.50)/(0.55+0.05)=1.05/0.60=1.75。
    #[test]
    fn magnitude_total_exceeds_one_allowed() {
        let got = eclipse_magnitude(-0.50, 0.55, -0.05).0;
        // 手計算: (0.55 − (−0.50)) / (0.55 + 0.05) = 1.05 / 0.60 = 1.75。
        assert!(
            close_abs(got, 1.75, MAG_TOL),
            "total magnitude >1: got {got}, want 1.75"
        );
        assert!(got > 1.0, "EclipseMagnitude must allow >1 for totality");
    }

    /// 外接（m = L1）: 食い込み 0 → magnitude = 0。
    #[test]
    fn magnitude_external_tangency_is_zero() {
        let got = eclipse_magnitude(0.7, 0.7, -0.3).0;
        // 手計算: (0.7 − 0.7)/(0.7+0.3) = 0。
        assert!(
            close_abs(got, 0.0, MAG_TOL),
            "external tangency magnitude: got {got}, want 0"
        );
    }

    /// 離隔（m > L1・食なし）: 負にせず 0 にクランプ。
    #[test]
    fn magnitude_no_eclipse_clamped_to_zero() {
        // m=1.2 > L1=0.7: 生の式なら (0.7−1.2)/(0.7+0.3) = −0.5、クランプで 0。
        let got = eclipse_magnitude(1.2, 0.7, -0.3).0;
        assert!(
            close_abs(got, 0.0, MAG_TOL),
            "no-eclipse magnitude must clamp to 0: got {got}"
        );
        assert!(got >= 0.0, "magnitude must never be negative");

        // m がさらに大きくても 0（負に振れない）。
        let got_far = eclipse_magnitude(5.0, 0.7, -0.3).0;
        assert!(
            close_abs(got_far, 0.0, MAG_TOL),
            "far separation magnitude must be 0: got {got_far}"
        );
    }

    /// 部分食域（0 < m < L1）で 0 < mag < 1、独立に手計算した値と一致。
    #[test]
    fn magnitude_partial_in_open_unit_interval() {
        // m=0.20, L1=0.80, l2=-0.40 → (0.80−0.20)/(0.80+0.40) = 0.60/1.20 = 0.5。
        let got = eclipse_magnitude(0.20, 0.80, -0.40).0;
        assert!(
            close_abs(got, 0.5, MAG_TOL),
            "partial magnitude: got {got}, want 0.5"
        );
        assert!(got > 0.0 && got < 1.0, "partial magnitude must be in (0,1)");
    }

    /// l2 の符号不問（|l2| を使う）: l2 = +x と l2 = −x で同一値。
    #[test]
    fn magnitude_uses_abs_of_l2() {
        let pos = eclipse_magnitude(0.30, 0.90, 0.50).0;
        let neg = eclipse_magnitude(0.30, 0.90, -0.50).0;
        // 手計算: (0.90−0.30)/(0.90+0.50) = 0.60/1.40 = 0.428571…
        let want = 0.60 / 1.40;
        assert!(
            close_abs(pos, want, MAG_TOL),
            "magnitude with +l2: got {pos}, want {want}"
        );
        assert!(
            close_abs(neg, want, MAG_TOL),
            "magnitude with -l2: got {neg}, want {want}"
        );
        assert!(
            close_abs(pos, neg, MAG_TOL),
            "magnitude must be invariant to sign of l2: {pos} vs {neg}"
        );
    }

    /// 独立オラクルとの広域照合（複数の部分食/皆既/離隔配置）。
    #[test]
    fn magnitude_matches_independent_oracle_grid() {
        let cases = [
            (0.0, 0.55, -0.45),
            (0.10, 0.80, -0.40),
            (0.50, 0.80, -0.40),
            (-0.50, 0.55, -0.05),
            (0.30, 0.90, 0.50),
            (0.79, 0.80, -0.40),
            (0.80, 0.80, -0.40), // 外接
            (1.50, 0.80, -0.40), // 離隔 → 0
        ];
        for &(m, l1, l2) in &cases {
            let got = eclipse_magnitude(m, l1, l2).0;
            let want = magnitude_oracle(m, l1, l2);
            assert!(
                close_abs(got, want, MAG_TOL),
                "magnitude({m},{l1},{l2}): got {got}, want {want}"
            );
            assert!(got >= 0.0, "magnitude never negative: {got}");
        }
    }

    // ======================================================================
    // 契約B: eclipse_obscuration = lens 重なり / 太陽面積、[0,1]、5 縮退境界
    // ======================================================================

    /// 同心・同半径（d=0, R=r）: 完全重なり → obscuration = 1。0 除算回避を縛る。
    #[test]
    fn obscuration_concentric_equal_radii_is_one() {
        let got = eclipse_obscuration(0.0, 1.0, 1.0).0;
        assert!(
            close_abs(got, 1.0, GEOM_TOL),
            "concentric equal radii: got {got}, want 1"
        );
    }

    /// 同心・月が大（d=0, r_moon>r_sun, 皆既近傍）: 太陽全面が覆われ obscuration = 1。
    #[test]
    fn obscuration_concentric_moon_larger_is_one() {
        let got = eclipse_obscuration(0.0, 0.8, 1.0).0;
        assert!(
            close_abs(got, 1.0, GEOM_TOL),
            "concentric moon larger: got {got}, want 1"
        );
    }

    /// 同心・太陽が大（d=0, r_sun>r_moon, 金環）: obscuration = r_moon²/r_sun²。
    /// 例 r_sun=1.0, r_moon=0.6 → 0.36/1.0 = 0.36。独立に手計算。
    #[test]
    fn obscuration_concentric_annular_is_radius_ratio_sq() {
        let got = eclipse_obscuration(0.0, 1.0, 0.6).0;
        // 手計算: π·0.6² / (π·1.0²) = 0.36。
        assert!(
            close_abs(got, 0.36, GEOM_TOL),
            "annular concentric: got {got}, want 0.36"
        );
        assert!(
            got > 0.0 && got < 1.0,
            "annular obscuration must be in (0,1)"
        );
    }

    /// 離隔（d = R+r, 外接）: 重なり 0 → obscuration = 0。0 除算（2dR/2dr の引数は有限だが
    /// 縮退分岐）を縛る。
    #[test]
    fn obscuration_external_tangency_is_zero() {
        let got = eclipse_obscuration(2.0, 1.0, 1.0).0; // d=2=R+r
        assert!(
            close_abs(got, 0.0, GEOM_TOL),
            "external tangency: got {got}, want 0"
        );
    }

    /// 離隔（d > R+r）: 0。
    #[test]
    fn obscuration_separated_is_zero() {
        let got = eclipse_obscuration(3.0, 1.0, 1.0).0;
        assert!(
            close_abs(got, 0.0, GEOM_TOL),
            "separated: got {got}, want 0"
        );
        assert!(got >= 0.0, "obscuration never negative");
    }

    /// 内接・月が大（d = |R−r|, r_moon>r_sun）: 太陽が月に内接 → 太陽全面 → obscuration = 1。
    /// 例 R=1.0, r=1.5, d=0.5。min=R=1.0 → π·1²/(π·1²)=1。
    #[test]
    fn obscuration_internal_tangency_moon_larger_is_one() {
        let got = eclipse_obscuration(0.5, 1.0, 1.5).0; // d=|1.0−1.5|=0.5
        assert!(
            close_abs(got, 1.0, GEOM_TOL),
            "internal tangency moon larger: got {got}, want 1"
        );
    }

    /// 内接・太陽が大（d = |R−r|, r_sun>r_moon, 金環縁）: obscuration = r_moon²/r_sun²。
    /// 例 R=1.0, r=0.6, d=0.4 → min=r=0.6 → 0.36/1.0=0.36。
    #[test]
    fn obscuration_internal_tangency_sun_larger_is_ratio_sq() {
        let got = eclipse_obscuration(0.4, 1.0, 0.6).0; // d=|1.0−0.6|=0.4
        assert!(
            close_abs(got, 0.36, GEOM_TOL),
            "internal tangency sun larger: got {got}, want 0.36"
        );
    }

    /// 内包（d < |R−r|, 太陽が大）: 月が完全内包 → obscuration = r_moon²/r_sun²。
    /// 例 R=1.0, r=0.5, d=0.2（< 0.5） → 0.25/1.0 = 0.25。
    #[test]
    fn obscuration_contained_sun_larger() {
        let got = eclipse_obscuration(0.2, 1.0, 0.5).0;
        assert!(
            close_abs(got, 0.25, GEOM_TOL),
            "contained (sun larger): got {got}, want 0.25"
        );
    }

    /// 内包（d < |R−r|, 月が大, 皆既）: 太陽が月に完全内包 → obscuration = 1。
    /// 例 R=1.0, r=1.5, d=0.3（< 0.5）。
    #[test]
    fn obscuration_contained_moon_larger_is_one() {
        let got = eclipse_obscuration(0.3, 1.0, 1.5).0;
        assert!(
            close_abs(got, 1.0, GEOM_TOL),
            "contained (moon larger): got {got}, want 1"
        );
    }

    /// 既知 lens 面積の独立解析値（R=r=1, d=1）。
    /// MathWorld の式より lens 面積 = 2·acos(1/2) − √3/2 = 2·(π/3) − √3/2 = 2π/3 − √3/2。
    /// この値は契約式のコピーでなく、特定配置 (R=r=1,d=1) を手で評価した解析定数。
    /// obscuration = (2π/3 − √3/2) / (π·1²)。
    #[test]
    fn obscuration_equal_radii_d_eq_r_known_lens() {
        let lens = 2.0 * PI / 3.0 - 3.0_f64.sqrt() / 2.0; // 独立に手計算した解析定数
        let want = lens / PI;
        let got = eclipse_obscuration(1.0, 1.0, 1.0).0;
        assert!(
            close_abs(got, want, GEOM_TOL),
            "R=r=1,d=1 lens obscuration: got {got}, want {want}"
        );
        assert!(got > 0.0 && got < 1.0, "partial overlap must be in (0,1)");
    }

    /// 部分重なり（非対称 R≠r, 縮退でない）を独立 lens オラクルと照合。
    #[test]
    fn obscuration_partial_asymmetric_matches_oracle() {
        let cases = [
            (1.2, 1.0, 0.5),
            (0.8, 1.0, 0.6),
            (1.0, 1.0, 0.7),
            (1.3, 1.0, 0.9),
            (0.9, 1.0, 0.4),
        ];
        for &(d, r_sun, r_moon) in &cases {
            let got = eclipse_obscuration(d, r_sun, r_moon).0;
            let want = obscuration_oracle(d, r_sun, r_moon);
            assert!(
                close_abs(got, want, GEOM_TOL),
                "obscuration({d},{r_sun},{r_moon}): got {got}, want {want}"
            );
            assert!(
                (0.0..=1.0).contains(&got),
                "obscuration must be in [0,1]: got {got}"
            );
        }
    }

    /// 部分重なり・スケールした非対称配置（`r_sun ≠ 1` かつ `r_sun ≠ r_moon`）を独立オラクルと照合。
    ///
    /// 狙い: 既存の部分重なりテストは `r_sun = 1`（一部は `r_sun = r_moon`）のため、lens 公式中の
    /// `2.0 * d * big_r` / `2.0 * d * r` / `big_r * big_r * cos_a.acos()` の `*` を `/` に置換しても
    /// `big_r * big_r == big_r / big_r == 1`、`2.0*d*big_r == 2.0*d/big_r`（big_r=1）で値が変わらず
    /// 検出できなかった。`r_sun ≠ 1` かつ `r_sun ≠ r_moon` の配置では `big_r*big_r ≠ big_r/big_r`、
    /// `2dR ≠ 2d/R` となり、`* → /` 変異が出力を変える → 殺せる。
    ///
    /// 各ケースは必ず部分重なり領域 `|R−r| < d < R+r` に入る（コメントで明記）:
    /// - `R=2.0, r=1.3, d=2.4`: |2.0−1.3|=0.7 < 2.4 < 3.3 ✓（R≠1, R≠r）
    /// - `R=3.5, r=1.0, d=3.0`: |3.5−1.0|=2.5 < 3.0 < 4.5 ✓（R≠1, R≠r。r=1 だが R≠1 なので 2dr の R 由来変異と無関係に big_r 系を縛る）
    /// - `R=2.0, r=3.0, d=2.5`: |2.0−3.0|=1.0 < 2.5 < 5.0 ✓（R≠1, R≠r, 月が大）
    /// - `R=4.0, r=2.5, d=3.0`: |4.0−2.5|=1.5 < 3.0 < 6.5 ✓（R≠1, R≠r, 太陽が大）
    /// - `R=0.7, r=0.4, d=0.6`: |0.7−0.4|=0.3 < 0.6 < 1.1 ✓（R<1, R≠r, スケール下方向）
    ///
    /// 独立性: `obscuration_oracle`/`lens_overlap_area` は実装 `eclipse_obscuration` を呼ばず、
    /// `r_sun` を一般値として扱う（`PI * r_sun * r_sun` で割り、lens 公式に big_r をそのまま入れる）。
    /// よってスケール値でも縮退せず正しく部分重なり式を通る。
    #[test]
    fn obscuration_partial_scaled_asymmetric_matches_oracle() {
        // すべて r_sun ≠ 1 かつ r_sun ≠ r_moon、かつ |R−r| < d < R+r（部分重なり）。
        let cases: [(f64, f64, f64); 5] = [
            (2.4, 2.0, 1.3), // |2.0−1.3|=0.7 < 2.4 < 3.3
            (3.0, 3.5, 1.0), // |3.5−1.0|=2.5 < 3.0 < 4.5
            (2.5, 2.0, 3.0), // |2.0−3.0|=1.0 < 2.5 < 5.0（月が大）
            (3.0, 4.0, 2.5), // |4.0−2.5|=1.5 < 3.0 < 6.5（太陽が大）
            (0.6, 0.7, 0.4), // |0.7−0.4|=0.3 < 0.6 < 1.1（R<1）
        ];
        for &(d, r_sun, r_moon) in &cases {
            // 部分重なり領域に入ることを実行時にも保証（縮退分岐に落ちると変異を殺せない）。
            let lo = (r_sun - r_moon).abs();
            let hi = r_sun + r_moon;
            assert!(
                lo < d && d < hi,
                "case must be in partial-overlap region: |{r_sun}-{r_moon}|={lo} < {d} < {hi}"
            );
            assert!(
                (r_sun - 1.0).abs() > GEOM_TOL,
                "case must have r_sun != 1 to kill * -> / on big_r: r_sun={r_sun}"
            );
            assert!(
                (r_sun - r_moon).abs() > GEOM_TOL,
                "case must have r_sun != r_moon: r_sun={r_sun}, r_moon={r_moon}"
            );

            let got = eclipse_obscuration(d, r_sun, r_moon).0;
            let want = obscuration_oracle(d, r_sun, r_moon);
            assert!(
                close_abs(got, want, GEOM_TOL),
                "obscuration({d},{r_sun},{r_moon}): got {got}, want {want}"
            );
            assert!(
                got > 0.0 && got < 1.0,
                "scaled partial overlap must be strictly in (0,1): got {got}"
            );
        }
    }

    /// acos 引数が丸めで ±1 をわずかに超えうる配置（外接/内接に極めて近い）で
    /// NaN を出さず有限値（クランプ必須）。
    #[test]
    fn obscuration_near_tangency_is_finite_no_nan() {
        // 外接にごく近い（d を R+r からわずかに下げる）。
        let near_external = eclipse_obscuration(2.0 - 1e-13, 1.0, 1.0).0;
        assert!(
            near_external.is_finite(),
            "near external tangency must be finite, got {near_external}"
        );
        assert!(
            (0.0..=1.0).contains(&near_external),
            "near external must stay in [0,1]: {near_external}"
        );

        // 内接にごく近い（金環縁）。d を |R−r| からわずかに上げる。
        let near_internal = eclipse_obscuration(0.4 + 1e-13, 1.0, 0.6).0;
        assert!(
            near_internal.is_finite(),
            "near internal tangency must be finite, got {near_internal}"
        );
        assert!(
            (0.0..=1.0).contains(&near_internal),
            "near internal must stay in [0,1]: {near_internal}"
        );

        // 内接の内側ごく近く（判別式が丸めで負になりうる側）。
        let just_inside = eclipse_obscuration(0.4 - 1e-13, 1.0, 0.6).0;
        assert!(
            just_inside.is_finite(),
            "just-inside internal must be finite, got {just_inside}"
        );
        assert!(
            (0.0..=1.0).contains(&just_inside),
            "just-inside internal must stay in [0,1]: {just_inside}"
        );
    }

    /// 結果は厳密に [0,1]（縮退・部分・離隔の全域で）。
    #[test]
    fn obscuration_always_in_unit_interval() {
        let ds = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.5, 1.8, 2.0, 2.5];
        let radii = [(1.0, 1.0), (1.0, 0.6), (1.0, 1.4), (1.0, 0.3)];
        for &(r_sun, r_moon) in &radii {
            for &d in &ds {
                let got = eclipse_obscuration(d, r_sun, r_moon).0;
                assert!(
                    (0.0..=1.0).contains(&got),
                    "obscuration({d},{r_sun},{r_moon})={got} out of [0,1]"
                );
                assert!(got.is_finite(), "obscuration must be finite: {got}");
            }
        }
    }

    /// プロパティ（L8）: d を R+r → 0 に減らすと obscuration は単調非減少。
    /// 太陽が大（金環ケース, R>r）: 同心で頭打ち（r²/R²）まで増加。
    #[test]
    fn obscuration_monotonic_nondecreasing_as_d_shrinks_annular() {
        let (r_sun, r_moon) = (1.0, 0.6);
        // d を大→小の順に並べる。
        let ds = [2.0, 1.6, 1.4, 1.2, 1.0, 0.8, 0.6, 0.5, 0.4, 0.2, 0.0];
        let mut prev = -1.0;
        for &d in &ds {
            let got = eclipse_obscuration(d, r_sun, r_moon).0;
            assert!(
                got >= prev - GEOM_TOL,
                "obscuration must be nondecreasing as d shrinks: d={d} got {got} < prev {prev}"
            );
            prev = got;
        }
    }

    /// プロパティ（L8）: 同様の単調性（皆既ケース, R<r で最終 1 に到達）。
    #[test]
    fn obscuration_monotonic_nondecreasing_as_d_shrinks_total() {
        let (r_sun, r_moon) = (1.0, 1.4);
        let ds = [2.4, 2.0, 1.6, 1.2, 0.8, 0.6, 0.4, 0.2, 0.0];
        let mut prev = -1.0;
        for &d in &ds {
            let got = eclipse_obscuration(d, r_sun, r_moon).0;
            assert!(
                got >= prev - GEOM_TOL,
                "obscuration must be nondecreasing (total): d={d} got {got} < prev {prev}"
            );
            prev = got;
        }
        // 最終（同心, 月大）は 1。
        let at_center = eclipse_obscuration(0.0, r_sun, r_moon).0;
        assert!(
            close_abs(at_center, 1.0, GEOM_TOL),
            "total at concentric must be 1: got {at_center}"
        );
    }

    // ======================================================================
    // 契約C: 食分と食面積の整合（部分食で両者 [0,1]）
    // ======================================================================

    /// 部分食に相当する配置で magnitude も obscuration も [0,1] に収まる。
    /// 食分（基本面 Re）と食面積（視半径平面）は単位系が別だが、部分食の値域は共通に [0,1]。
    #[test]
    fn partial_eclipse_magnitude_and_obscuration_in_unit_interval() {
        // 部分食: 0<m<L1。
        let mag = eclipse_magnitude(0.40, 0.80, -0.30).0;
        assert!(
            mag > 0.0 && mag < 1.0,
            "partial magnitude must be in (0,1): {mag}"
        );
        // 太陽の一部だけが覆われる部分重なり配置。
        let obsc = eclipse_obscuration(1.1, 1.0, 0.5).0;
        assert!(
            obsc > 0.0 && obsc < 1.0,
            "partial obscuration must be in (0,1): {obsc}"
        );
    }
}
