//! 局地接触点の位置角 PA（`docs/issues/ISSUE-043` S7a、`docs/conventions.md` §7,
//! `docs/algorithms/09-local.md` 式9.22, Explanatory Supplement §11 / Meeus Ch.54）。
//!
//! 接触点の位置角 = **太陽面の天の北点から反時計回り（東回り）に測った接触角**
//! （NASA/Espenak の定義 "P is the contact angle measured counter-clockwise from the north
//! point of the Sun's disk"、conventions §7「位置角: 天の北を 0°、東回り」と一致, 確定 S7a）。
//!
//! ## 幾何（本プロジェクトの基本面規約: ξ̂=天の東, η̂=天の北。`projection.rs`）
//! 影軸が基本面を貫く点 (x, y)・観測者射影 (ξ, η) から、観測者から見た**月中心の太陽中心に対する
//! 見かけずれ**（天球上の東・北成分）は基本面ベクトル `(x−ξ, y−η)` に比例する（月は観測者と太陽の
//! 間にあり、基本面 +東 のずれ → 天球上も +東）。PA（天の北 0・東回り）は東成分・北成分の atan2:
//! ```text
//! PA = atan2(σ·(x−ξ), σ·(y−η))   （[0, 360) に正規化）
//! ```
//! σ は接触点が月中心の**どちら側**かを表す符号（確定 S7a §3）:
//! - 外接 C1/C4・金環の内接 C2/C3: 接触点は月中心**方向** ⇒ σ = +1。
//! - 皆既の内接 C2/C3（本影 l2<0・月⊃太陽）: 接触点は月中心と**反対**側 ⇒ σ = −1。
//!
//! 物理整合（C1 は西縁≈270°系・C4 は東縁≈90°系。月は太陽の西から追いつき東へ抜ける）。

// contact_position_angle（pub(crate)）は ISSUE-043 S7b（local_circumstances 結線）が消費するまで
// 未使用。結線され次第この許容は外す（projection.rs / local_contacts.rs と同手順）。
#![allow(dead_code)]

use umbra_core::{Degrees, Radians};

use crate::besselian::InstantaneousBesselianElements;
use crate::projection::ObserverFundamental;

/// 接触点の位置角 PA（天の北 0°・東回り, `[0, 360)`）を返す。
///
/// - `elements`: 接触時刻の瞬時ベッセル要素（影軸交点 x, y を使用, ISSUE-021）。
/// - `observer`: 観測者の基本面射影 (ξ, η)（ISSUE-024 `project_observer_to_fundamental`）。
/// - `umbral_interior`: 皆既の内接接触（C2/C3 かつ本影 l2<0）なら `true`（σ=−1）、
///   外接 C1/C4・金環の内接 C2/C3 なら `false`（σ=+1）。
///
/// `PA = atan2(σ·(x−ξ), σ·(y−η))` を `[0, 360)` 度へ正規化（確定 S7a, conventions §7）。
pub(crate) fn contact_position_angle(
    elements: &InstantaneousBesselianElements,
    observer: &ObserverFundamental,
    umbral_interior: bool,
) -> Degrees {
    // σ: 接触点が月中心の方向（+1）か反対側（−1, 皆既の内接 C2/C3）か（確定 S7a §3）。
    let sigma = if umbral_interior { -1.0 } else { 1.0 };
    // 見かけずれの東・北成分（基本面 ξ̂=東, η̂=北）。月は観測者と太陽の間 → (x−ξ, y−η) が +東/+北。
    let east = sigma * (elements.x - observer.xi);
    let north = sigma * (elements.y - observer.eta);
    // 天の北 0・東回り: atan2(東, 北)（北=基準軸, 東へ増加）。horizontal.rs の方位と同じ正規化で
    // [0, 2π) → 度（[0, 360)）。
    Radians::new(east.atan2(north))
        .normalized_two_pi()
        .to_degrees()
}

#[cfg(test)]
mod tests {
    //! ISSUE-043 S7a 受け入れテスト（strict・接触点の位置角 PA）。
    //!
    //! ## オラクル戦略（追認回避 — `PA = atan2(σ·(x−ξ), σ·(y−η))` を**写経しない**）
    //!
    //! PA の契約は「天の北 0°・東回り（北→東→南→西 = 0→90→180→270）, `[0,360)`」。
    //! 基本面規約は ξ̂=天の東・η̂=天の北、見かけずれ (Δe, Δn) = (x−ξ, y−η)。
    //!
    //! 主オラクル（決定的）は **基数/対角方向の手計算値**を直接 assert する。北0東回りの定義
    //! だけから手で確定でき、実装式の代数とは独立:
    //!   真北   (Δe=0,  Δn>0) → 0°
    //!   真東   (Δe>0,  Δn=0) → 90°
    //!   真南   (Δe=0,  Δn<0) → 180°
    //!   真西   (Δe<0,  Δn=0) → 270°
    //!   北東   (Δe>0,  Δn>0, 等量) → 45°（他象限の対角も同様）
    //! これらは方位の定義（北を 0 とし東へ増える）から直接出る既知値であり、`atan2(Δe,Δn)` という
    //! 式を経由せずに「方向→角度」の対応表として独立に書き下している（写経ではない）。
    //!
    //! 補助オラクル:
    //! - σ（皆既内接=true）は同じ (x,y,ξ,η) に対し PA を 180° 反転させる（mod 360）。
    //! - スケール不変: ずれベクトルを正係数倍しても PA 不変（方向のみ依存）。
    //! - 値域 `[0,360)`、フィールド取り違え撃破（x,y,ξ,η を相異な非対称値に）、無関係フィールド
    //!   （declination/mu/zeta/l1/l2/tan_f1/tan_f2/rate/time）の不変。
    //! - 物理整合（第二義・緩い帯）: 2017-08-21 実日食の C1/C4 を独立検証済 solver から取り、
    //!   C1 が西半分(180,360)・C4 が東半分(0,180) になることを半球判定で縛る。
    //!
    //! 出典: NASA/Espenak「P measured counter-clockwise from the north point」, conventions §7。

    // 和文 doc 箇条書き（手順説明）が doc lint と誤認される場合に備える（既存テストの作法に倣う）。
    #![allow(clippy::doc_lazy_continuation)]

    use super::*;

    use umbra_core::{JulianDate2, Radians, TtInstant};

    /// PA 手計算オラクルの許容（度）。基数/対角方向は丸めのみなので極めて厳しく取る。
    const TOL_DEG: f64 = 1e-9;

    /// 周期量 PA の差を `[-180,180]` に畳んで比較（北0近傍の巻き戻りを跨いでも頑健）。
    fn angles_close(got: f64, want: f64, tol: f64) -> bool {
        let mut d = (got - want).rem_euclid(360.0);
        if d > 180.0 {
            d -= 360.0;
        }
        d.abs() < tol
    }

    /// テスト用に瞬時ベッセル要素を最小構成する。PA に効くのは影軸交点 x,y のみ。
    /// 他フィールド（declination/mu/l1/l2/tan_f1/tan_f2/time_tt）は PA に無関係なので
    /// 非対称なダミー有限値を入れる（無関係フィールド不変テストで実際に振って確かめる）。
    fn elems_xy(x: f64, y: f64) -> InstantaneousBesselianElements {
        InstantaneousBesselianElements {
            x,
            y,
            declination: Radians::new(0.207),
            mu: Radians::new(1.0),
            l1: 0.54,
            l2: -0.009,
            tan_f1: 0.004_65,
            tan_f2: 0.004_63,
            time_tt: TtInstant::from_jd2(JulianDate2::new(2_451_545.0, 0.0)),
        }
    }

    /// テスト用に観測者基本面射影を構成する。PA に効くのは ξ,η のみ。
    /// rate/zeta は非対称な有限ダミー（無関係フィールド不変テストで振って確かめる）。
    fn obs_xieta(xi: f64, eta: f64) -> ObserverFundamental {
        ObserverFundamental {
            xi,
            eta,
            zeta: 0.93,
            xi_rate: 1.1e-6,
            eta_rate: -2.3e-6,
        }
    }

    // ============================================================
    // 主オラクル: 基数/対角方向の手計算 PA（北0東回り）
    // ============================================================

    /// 真北のずれ（Δe=x−ξ=0, Δn=y−η>0）→ PA=0°。
    /// x↔y/ξ↔η/x↔ξ 取り違えで撃たれるよう 4 値すべて相異な非対称値にする。
    /// x=0.30, ξ=0.30 で Δe=0；y=0.50, η=0.20 で Δn=+0.30(>0)。手計算: 北 → 0°。
    #[test]
    fn cardinal_north_is_zero_degrees() {
        let e = elems_xy(0.30, 0.50);
        let o = obs_xieta(0.30, 0.20);
        let pa = contact_position_angle(&e, &o, false);
        assert!(
            angles_close(pa.0, 0.0, TOL_DEG),
            "Δe=0, Δn>0 (true north) must be PA=0°, got {}",
            pa.0
        );
    }

    /// 真東（Δe>0, Δn=0）→ PA=90°。
    /// x=0.40, ξ=0.10 → Δe=+0.30(>0)；y=0.25, η=0.25 → Δn=0。手計算: 東 → 90°。
    #[test]
    fn cardinal_east_is_ninety_degrees() {
        let e = elems_xy(0.40, 0.25);
        let o = obs_xieta(0.10, 0.25);
        let pa = contact_position_angle(&e, &o, false);
        assert!(
            angles_close(pa.0, 90.0, TOL_DEG),
            "Δe>0, Δn=0 (true east) must be PA=90°, got {}",
            pa.0
        );
    }

    /// 真南（Δe=0, Δn<0）→ PA=180°。
    /// x=0.15, ξ=0.15 → Δe=0；y=0.20, η=0.55 → Δn=−0.35(<0)。手計算: 南 → 180°。
    #[test]
    fn cardinal_south_is_one_eighty_degrees() {
        let e = elems_xy(0.15, 0.20);
        let o = obs_xieta(0.15, 0.55);
        let pa = contact_position_angle(&e, &o, false);
        assert!(
            angles_close(pa.0, 180.0, TOL_DEG),
            "Δe=0, Δn<0 (true south) must be PA=180°, got {}",
            pa.0
        );
    }

    /// 真西（Δe<0, Δn=0）→ PA=270°。
    /// x=0.05, ξ=0.45 → Δe=−0.40(<0)；y=0.33, η=0.33 → Δn=0。手計算: 西 → 270°。
    #[test]
    fn cardinal_west_is_two_seventy_degrees() {
        let e = elems_xy(0.05, 0.33);
        let o = obs_xieta(0.45, 0.33);
        let pa = contact_position_angle(&e, &o, false);
        assert!(
            angles_close(pa.0, 270.0, TOL_DEG),
            "Δe<0, Δn=0 (true west) must be PA=270°, got {}",
            pa.0
        );
    }

    /// 4 対角方向（北東45 / 南東135 / 南西225 / 北西315）。等量成分で手計算。
    /// 各ケースで x,y,ξ,η を相異な非対称値にし、Δe・Δn が ±等量になるよう設計。
    #[test]
    fn diagonal_directions_are_45_135_225_315() {
        // 北東: Δe>0, Δn>0 等量 → 45°。x=0.50,ξ=0.20→Δe=+0.30; y=0.60,η=0.30→Δn=+0.30。
        let pa_ne = contact_position_angle(&elems_xy(0.50, 0.60), &obs_xieta(0.20, 0.30), false);
        assert!(
            angles_close(pa_ne.0, 45.0, TOL_DEG),
            "NE (Δe>0,Δn>0 equal) must be 45°, got {}",
            pa_ne.0
        );
        // 南東: Δe>0, Δn<0 等量 → 135°。x=0.55,ξ=0.25→Δe=+0.30; y=0.10,η=0.40→Δn=−0.30。
        let pa_se = contact_position_angle(&elems_xy(0.55, 0.10), &obs_xieta(0.25, 0.40), false);
        assert!(
            angles_close(pa_se.0, 135.0, TOL_DEG),
            "SE (Δe>0,Δn<0 equal) must be 135°, got {}",
            pa_se.0
        );
        // 南西: Δe<0, Δn<0 等量 → 225°。x=0.10,ξ=0.40→Δe=−0.30; y=0.05,η=0.35→Δn=−0.30。
        let pa_sw = contact_position_angle(&elems_xy(0.10, 0.05), &obs_xieta(0.40, 0.35), false);
        assert!(
            angles_close(pa_sw.0, 225.0, TOL_DEG),
            "SW (Δe<0,Δn<0 equal) must be 225°, got {}",
            pa_sw.0
        );
        // 北西: Δe<0, Δn>0 等量 → 315°。x=0.05,ξ=0.35→Δe=−0.30; y=0.60,η=0.30→Δn=+0.30。
        let pa_nw = contact_position_angle(&elems_xy(0.05, 0.60), &obs_xieta(0.35, 0.30), false);
        assert!(
            angles_close(pa_nw.0, 315.0, TOL_DEG),
            "NW (Δe<0,Δn>0 equal) must be 315°, got {}",
            pa_nw.0
        );
    }

    /// 非等量の一般方向（手計算の atan2 ではなく、既知の三角値で角度を確定）。
    /// Δe=+1.0, Δn=+√3 → 北から東へ tan⁻¹(1/√3)=30° → PA=30°。
    /// x,y,ξ,η は相異な非対称値（x=1.2,ξ=0.2→Δe=1.0; y=0.5+√3,η=0.5→Δn=√3）。
    #[test]
    fn general_direction_known_trig_angle_30_degrees() {
        let sqrt3 = 3.0_f64.sqrt();
        let e = elems_xy(1.2, 0.5 + sqrt3);
        let o = obs_xieta(0.2, 0.5);
        let pa = contact_position_angle(&e, &o, false);
        assert!(
            angles_close(pa.0, 30.0, TOL_DEG),
            "Δe=1, Δn=√3 must be PA=30° (30° east of north), got {}",
            pa.0
        );
    }

    // ============================================================
    // σ（皆既内接）: umbral_interior が PA を 180° 反転させる
    // ============================================================

    /// 同一 (x,y,ξ,η) で umbral_interior=false/true を呼ぶと PA がちょうど 180° 反転（mod 360）。
    /// σ=+1（外接・金環内接）に対し σ=−1（皆既内接）は接触点が月中心の反対側 → 方向反転。
    /// 複数の方向で確認（取り違え・部分反転を撃破）。
    #[test]
    fn umbral_interior_flips_pa_by_180() {
        let cases = [
            (elems_xy(0.30, 0.50), obs_xieta(0.30, 0.20)), // 北
            (elems_xy(0.40, 0.25), obs_xieta(0.10, 0.25)), // 東
            (elems_xy(0.55, 0.10), obs_xieta(0.25, 0.40)), // 南東
            (elems_xy(0.05, 0.60), obs_xieta(0.35, 0.30)), // 北西
        ];
        for (e, o) in cases {
            let outer = contact_position_angle(&e, &o, false);
            let inner = contact_position_angle(&e, &o, true);
            assert!(
                angles_close(inner.0, outer.0 + 180.0, 1e-7),
                "umbral_interior must flip PA by 180° (mod 360): outer={} inner={}",
                outer.0,
                inner.0
            );
        }
    }

    // ============================================================
    // スケール不変・値域・取り違え/無関係フィールド
    // ============================================================

    /// ずれベクトル (Δe, Δn) を正の係数倍しても PA 不変（atan2 は方向のみ依存）。
    /// 基準: x=0.50,ξ=0.20,y=0.70,η=0.30 → (Δe,Δn)=(0.30,0.40)。これを 3 倍したずれを
    /// 別の非対称 (x,y,ξ,η) で再現（ξ,η を変えてオフセットを散らす）し PA 一致を縛る。
    #[test]
    fn pa_is_scale_invariant_in_offset() {
        let base = contact_position_angle(&elems_xy(0.50, 0.70), &obs_xieta(0.20, 0.30), false);
        // (Δe,Δn)=(0.30,0.40)×3=(0.90,1.20)。ξ=0.05,η=−0.10 を起点に x,y を作る。
        let xi2 = 0.05;
        let eta2 = -0.10;
        let scaled = contact_position_angle(
            &elems_xy(xi2 + 0.90, eta2 + 1.20),
            &obs_xieta(xi2, eta2),
            false,
        );
        assert!(
            angles_close(base.0, scaled.0, 1e-9),
            "PA must be invariant under positive scaling of the offset: base={} scaled={}",
            base.0,
            scaled.0
        );
    }

    /// 値域: 種々の (x,y,ξ,η) で PA ∈ [0,360)。負値・360 以上を出さない。
    #[test]
    fn pa_is_in_zero_to_360_range() {
        let cases = [
            (elems_xy(0.30, 0.50), obs_xieta(0.30, 0.20), false),
            (elems_xy(0.40, 0.25), obs_xieta(0.10, 0.25), true),
            (elems_xy(0.05, 0.33), obs_xieta(0.45, 0.33), false),
            (elems_xy(-0.7, -0.2), obs_xieta(0.3, 0.9), true),
            (elems_xy(0.12, -0.88), obs_xieta(-0.44, 0.05), false),
            (elems_xy(0.0, 0.0), obs_xieta(-0.3, -0.6), true),
        ];
        for (e, o, interior) in cases {
            let pa = contact_position_angle(&e, &o, interior);
            assert!(
                (0.0..360.0).contains(&pa.0),
                "PA must be in [0,360), got {} (interior={interior})",
                pa.0
            );
        }
    }

    /// 無関係フィールドの不変: declination/mu/l1/l2/tan_f1/tan_f2/time_tt（要素側）と
    /// zeta/xi_rate/eta_rate（観測者側）を大きく振っても、(x,y,ξ,η) が同じなら PA は不変。
    /// 取り違え（例: η の代わりに zeta を、y の代わりに declination を使う）変異を撃破する。
    #[test]
    fn unrelated_fields_do_not_affect_pa() {
        let e1 = InstantaneousBesselianElements {
            x: 0.40,
            y: 0.65,
            declination: Radians::new(0.10),
            mu: Radians::new(0.5),
            l1: 0.50,
            l2: -0.005,
            tan_f1: 0.004,
            tan_f2: 0.0039,
            time_tt: TtInstant::from_jd2(JulianDate2::new(2_451_545.0, 0.0)),
        };
        let e2 = InstantaneousBesselianElements {
            x: 0.40, // 同じ
            y: 0.65, // 同じ
            declination: Radians::new(-1.3),
            mu: Radians::new(5.9),
            l1: 0.61,
            l2: -0.02,
            tan_f1: 0.0048,
            tan_f2: 0.0047,
            time_tt: TtInstant::from_jd2(JulianDate2::new(2_460_000.0, 0.5)),
        };
        let o1 = ObserverFundamental {
            xi: 0.15,
            eta: 0.20,
            zeta: 0.90,
            xi_rate: 1.0e-6,
            eta_rate: 2.0e-6,
        };
        let o2 = ObserverFundamental {
            xi: 0.15,  // 同じ
            eta: 0.20, // 同じ
            zeta: -0.40,
            xi_rate: -9.9e-5,
            eta_rate: 7.7e-5,
        };
        for &interior in &[false, true] {
            let a = contact_position_angle(&e1, &o1, interior);
            let b = contact_position_angle(&e2, &o2, interior);
            assert!(
                angles_close(a.0, b.0, 1e-12),
                "PA must depend only on (x,y,ξ,η); unrelated fields changed result: {} vs {} (interior={interior})",
                a.0,
                b.0
            );
        }
    }

    /// フィールド取り違え撃破の集中ケース: 4 値すべて相異・象限が一意に決まる配置。
    /// (Δe,Δn)=(+0.20,−0.40)（南東寄り・第二象限の東側＝90〜180）。x=0.30,ξ=0.10→Δe=0.20;
    /// y=0.05,η=0.45→Δn=−0.40。手計算 PA = 90 + tan⁻¹(0.40/0.20)= 90+63.4349…=153.4349…°。
    /// tan⁻¹(2) を独立に評価（atan2 写経ではなく、東/南成分比の既知逆正接）。
    #[test]
    fn asymmetric_quadrant_pins_field_mapping() {
        let e = elems_xy(0.30, 0.05);
        let o = obs_xieta(0.10, 0.45);
        let pa = contact_position_angle(&e, &o, false);
        // 南東（東>0, 北<0）→ 90°〜180°。角度 = 90° + atan(|Δn|/Δe) を独立評価（東から南へ）。
        // ここで atan(2.0) は east 基準からの傾き。北0東回りでは PA = atan2(Δe,Δn) と一致するが、
        // テスト側は「東(90°)から南へ atan(|Δn|/|Δe|) 進む」という方向合成で独立に角度を出す。
        let want = 90.0 + (0.40_f64 / 0.20_f64).atan().to_degrees();
        assert!(
            angles_close(pa.0, want, TOL_DEG),
            "SE asymmetric offset must be {want}° (east then toward south), got {}",
            pa.0
        );
    }

    // ============================================================
    // 物理整合（第二義・緩い帯）: 2017-08-21 実日食の C1/C4
    // ============================================================
    //
    // 観測者・定数・窓・config は他モジュール test ヘルパを参照できないため自前で再掲する
    // （local_maximum.rs / local_contacts.rs の central_observer / tt_2017_max / make_source /
    //  config_tight / window_around 相当）。重い/不安定なら本テストは省略可だが、半球判定の
    //  緩い帯（厳密秒値・厳密角を避ける）なので決定的に通るはず。

    use crate::conjunction::RootConfig;
    use crate::local_contacts::solve_local_contacts;
    use crate::projection::project_observer_to_fundamental;
    use crate::source::{BesselianSource, DirectBesselianSource};
    use umbra_core::constants::{EARTH_EQUATORIAL_RADIUS_M, SOLAR_RADIUS_KM};
    use umbra_core::ellipsoid::{observer_geocentric, Ellipsoid};
    use umbra_core::{EspenakMeeusDeltaT, TimeInterval};

    const WGS84: Ellipsoid = Ellipsoid::WGS84;
    const R_SUN: f64 = SOLAR_RADIUS_KM;
    const R_MOON: f64 = 0.2725076 * (EARTH_EQUATORIAL_RADIUS_M / 1000.0);

    fn tt2(jd1: f64, jd2: f64) -> TtInstant {
        TtInstant::from_jd2(JulianDate2::new(jd1, jd2))
    }
    fn tt_jd(jd: f64) -> TtInstant {
        TtInstant::from_jd2(JulianDate2::from_jd(jd))
    }
    fn jd_of(t: TtInstant) -> f64 {
        t.jd2().jd()
    }
    /// 2017-08-21 皆既日食 最大食付近の TT エポック（local_maximum.rs と同一）。
    fn tt_2017_max() -> TtInstant {
        tt2(2_457_986.5, 7.685_322_222_222_222e-1)
    }
    fn source_interval() -> TimeInterval<TtInstant> {
        TimeInterval {
            start: tt2(2_457_985.0, 0.0),
            end: tt2(2_457_989.0, 0.0),
        }
    }
    fn config_tight() -> RootConfig {
        RootConfig {
            x_tolerance_days: 1e-9,
            max_iterations: 200,
        }
    }
    fn window_around(center: TtInstant, half_hours: f64) -> TimeInterval<TtInstant> {
        let c = jd_of(center);
        let half = half_hours / 24.0;
        TimeInterval {
            start: tt_jd(c - half),
            end: tt_jd(c + half),
        }
    }

    /// 物理整合（緩い帯）: 中心食地点（37.5°N, 西経89.2°）の C1 は太陽**西縁**側＝PA 西半分
    /// (180,360)、C4 は**東縁**側＝PA 東半分(0,180)。月は太陽の西から追いつき東へ抜けるので
    /// C1 接触点は西寄り・C4 接触点は東寄りになる。厳密角・秒値は使わず半球で縛る。
    /// C1/C4 は外接（umbral_interior=false, σ=+1）。
    #[test]
    fn physical_c1_is_western_half_c4_is_eastern_half_2017() {
        let dt = EspenakMeeusDeltaT;
        let src = DirectBesselianSource::new(R_SUN, R_MOON, &dt, source_interval());
        let lat = 37.5_f64.to_radians();
        let east_lon = Radians::new((-89.2_f64).to_radians());
        let obs = observer_geocentric(&WGS84, lat, 200.0);
        let w = window_around(tt_2017_max(), 3.0);

        let contacts = solve_local_contacts(&src, &obs, east_lon, w, config_tight())
            .expect("central site should yield a contact set");
        let c1 = contacts.c1.expect("central site must have C1").time_tt;
        let c4 = contacts.c4.expect("central site must have C4").time_tt;

        // C1 時点の要素・射影 → PA。
        let e1 = src.at(c1).expect("source.at(C1) should succeed");
        let p1 = project_observer_to_fundamental(&obs, east_lon, &e1);
        let pa_c1 = contact_position_angle(&e1, &p1, false).0;

        let e4 = src.at(c4).expect("source.at(C4) should succeed");
        let p4 = project_observer_to_fundamental(&obs, east_lon, &e4);
        let pa_c4 = contact_position_angle(&e4, &p4, false).0;

        assert!(
            (180.0..360.0).contains(&pa_c1),
            "C1 (sun's western limb) PA should be in western half (180,360), got {pa_c1}"
        );
        assert!(
            (0.0..180.0).contains(&pa_c4),
            "C4 (sun's eastern limb) PA should be in eastern half (0,180), got {pa_c4}"
        );
    }
}
