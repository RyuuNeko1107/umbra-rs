//! 日食候補フィルタ（ISSUE-018）。
//!
//! 合（ISSUE-017 [`crate::conjunction::Conjunction`]）が日食を起こしうるかを、安価な幾何量
//! （合時刻の月-太陽角距離・視半径・月地平視差）で早期判定する層。Meeus Ch.54 の食限を保守側
//! （偽陰性ゼロ）に評価し、明らかに食帯が地球を外す朔だけを棄却する。後段の高コストなベッセル
//! 要素計算（ISSUE-019+）に渡す候補を絞り込むのが目的。
//!
//! **契約の核**: `possible = (separation < s_sun + s_moon + π_moon + margin)`。実際の日食を 1 件も
//! 落とさない（偽陰性不可）。グレーゾーンは `possible=true` 側に倒す（偽陽性可）。`min_separation`/
//! `sum_semidiameters`/`approx_gamma` は情報用の付随量で、いずれもテスト側の独立式で検証する。
//!
//! 判定式（Meeus Ch.54 食限）: `possible = (separation < s_sun + s_moon + π_moon + margin)`。
//! `s_sun=asin(R_sun/d_sun)`, `s_moon=asin(k·Re/d_moon)`, `π_moon=asin(Re/d_moon)`（月地平視差, 偽陰性
//! 回避に必須）。margin は偽陰性ゼロのマージン（D6: 合の角距離が最大食付近の真の最小値を上回る分＝
//! 月相対角速度×合↔最大食ずれ ＋ 概算暦誤差）を保守側（広め）に取る。
//!
//! 注: 仕様の `eph: &impl Ephemeris` / `config: &EngineConfig` は省略。apparent が VSOP/ELP 直結
//! （ISSUE-037 と同じ繰延）、`EngineConfig` 未実装のため k は IauMean 既定（フィルタは部分食=半影基準ゆえ
//! k 差は判定境界に影響しない）。暦・config の差し替えは ISSUE-043 で統合する。本フィルタは無謬（Result 無し）。

// assess_eclipse_possibility は ISSUE-023（種別確定）/ search が消費する。結線され次第この許容は外す。
#![allow(dead_code)]

use umbra_core::constants::{EARTH_EQUATORIAL_RADIUS_M, SOLAR_RADIUS_KM};
use umbra_core::Radians;
use umbra_ephemeris::apparent::{moon_apparent_cirs, sun_apparent_cirs};

use crate::conjunction::Conjunction;

/// 月/太陽サイズ比 k（IauMean, conventions §9 既定）。
const MOON_SUN_SIZE_RATIO_K: f64 = 0.272_507_6;
/// 偽陰性ゼロのマージン \[rad\]（D6 の (2)+(3) 項; 月地平視差 (1) は別途加算）。合の角距離は最大食付近の
/// 真の最小値より大きくなりうる（月相対角速度 ~0.5°/h × 合↔最大食ずれ ~0.2h ≲ 0.1°）＋概算暦誤差 ~0.05°
/// を保守側に広げて 0.5°≈0.0087 rad。偽陰性ゼロ側＝広めに固定（conventions §11, 誤差を隠さない）。
/// D6 マージン実余裕統計（`filter_margins`）が `ECLIPSE_FILTER_SAFETY_MARGIN_RAD` として公開する。
pub(crate) const SAFETY_MARGIN_RAD: f64 = 0.008_7;

/// 棄却/採用の理由（デバッグ・テスト用）。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum PossibilityReason {
    /// 日食が起こりうる（後段のベッセル/全球判定で確定, ISSUE-021/023）。
    PossibleEclipse,
    /// 合の角距離 > 視半径和 + 月視差 + 保守マージン → 食帯が地球を外す。
    SeparationTooLarge,
}

/// 日食可能性の判定結果（判定根拠を保持し偽陽性/偽陰性のデバッグを可能にする）。
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct EclipsePossibility {
    /// 日食が起こりうるか（偽陽性可・偽陰性不可）。
    pub possible: bool,
    /// 元の合（ISSUE-017）。
    pub conjunction: Conjunction,
    /// 合付近の最小角距離（概算, = `conjunction.separation` を保守的に代用）。
    pub min_separation: Radians,
    /// 太陽視半径 + 月視半径（概算, rad）。
    pub sum_semidiameters: Radians,
    /// 食限（マージン抜き）= s_sun + s_moon + π_moon \[rad\]（D6 マージン実余裕統計の基準・`filter_margins`）。
    pub bare_limit: Radians,
    /// 影軸-地球中心 概算距離（Re 無次元, 情報用）。
    pub approx_gamma: f64,
    /// 棄却/採用の理由。
    pub reason: PossibilityReason,
}

/// 合が日食を起こしうるか安価な幾何量で早期判定する（偽陰性ゼロ。明らかに食帯が地球を外す朔だけ棄却）。
///
/// `possible = (conjunction.separation < s_sun + s_moon + π_moon + margin)`（Meeus Ch.54 食限, 保守側）。
/// 視半径・視差は合時刻の地心距離から、`approx_gamma` は地心→Sun-Moon 線の垂線距離 `|S×M|/(|M−S|·Re)` で
/// 概算する（情報用。フル gamma より必ず甘い＝採用寄り）。
pub(crate) fn assess_eclipse_possibility(conjunction: &Conjunction) -> EclipsePossibility {
    let t = conjunction.time_tt;
    let re_km = EARTH_EQUATORIAL_RADIUS_M / 1000.0;
    let sun = sun_apparent_cirs(t);
    let moon = moon_apparent_cirs(t);
    let d_sun = sun.norm();
    let d_moon = moon.norm();

    // 視半径 s = asin(R/d)、月地平視差 π_moon = asin(Re/d_moon)。asin 引数はクランプ（数値安全）。
    let s_sun = (SOLAR_RADIUS_KM / d_sun).clamp(-1.0, 1.0).asin();
    let s_moon = (MOON_SUN_SIZE_RATIO_K * re_km / d_moon)
        .clamp(-1.0, 1.0)
        .asin();
    let parallax_moon = (re_km / d_moon).clamp(-1.0, 1.0).asin();
    let sum_semidiameters = s_sun + s_moon;

    // 影軸-地球中心 概算距離（Re）: 地心→Sun-Moon 直線の垂線距離 = |S×M|/(|M−S|·Re)。
    let approx_gamma = sun.cross(moon).norm() / ((moon - sun).norm() * re_km);

    // Meeus Ch.54 食限（保守側）。separation がこの限界未満なら地球上のどこかで（部分）食。
    // bare_limit = マージン抜き食限（D6 マージン実余裕統計の基準）。limit = bare_limit + 保守マージン。
    let bare_limit = sum_semidiameters + parallax_moon;
    let limit = bare_limit + SAFETY_MARGIN_RAD;
    let possible = conjunction.separation.0 < limit;
    let reason = if possible {
        PossibilityReason::PossibleEclipse
    } else {
        PossibilityReason::SeparationTooLarge
    };

    EclipsePossibility {
        possible,
        conjunction: *conjunction,
        min_separation: conjunction.separation,
        sum_semidiameters: Radians(sum_semidiameters),
        bare_limit: Radians(bare_limit),
        approx_gamma,
        reason,
    }
}

#[cfg(test)]
mod tests {
    // 実装本体（同モジュール直下）が定義する型・関数。impl 担当の `use` 構成に依存しないよう、
    // テスト側で必要なシンボルを `super::*` で取り込む（conjunction.rs / candidates.rs と同手順）。
    use super::*;

    use crate::candidates::{new_moon_candidates, NewMoonCandidate};
    use crate::conjunction::{solve_conjunction, Conjunction, ConjunctionKind, RootConfig};

    use umbra_core::constants::{EARTH_EQUATORIAL_RADIUS_M, SOLAR_RADIUS_KM};
    use umbra_core::{JulianDate2, TimeRange, TtInstant, UtcInstant, Vector3};
    use umbra_ephemeris::apparent::{moon_apparent_cirs, sun_apparent_cirs};

    // ============================================================
    // 共通定数・ヘルパ
    // ============================================================

    /// 月/太陽サイズ比 k（IAU, conventions §9）。月視半径 s_moon = asin(k·Re / d_moon)。
    const MOON_SUN_SIZE_RATIO_K: f64 = 0.272_507_6;

    /// 地球赤道半径 \[km\]（視差・gamma 無次元化の基準 Re）。constants は m なので /1000。
    const EARTH_EQUATORIAL_RADIUS_KM: f64 = EARTH_EQUATORIAL_RADIUS_M / 1000.0;

    /// 許容つきスカラ比較（clippy::float_cmp 回避）。
    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    /// 単一 JD（UTC スケール）から UtcInstant。
    fn utc_from_jd(jd: f64) -> UtcInstant {
        UtcInstant::from_jd2(JulianDate2::from_jd(jd))
    }

    /// UTC の `[start_jd, end_jd]` 範囲。
    fn utc_range(start_jd: f64, end_jd: f64) -> TimeRange<UtcInstant> {
        TimeRange {
            start: utc_from_jd(start_jd),
            end: utc_from_jd(end_jd),
        }
    }

    /// TtInstant を JD（TT）で取り出す。
    fn jd_of(t: TtInstant) -> f64 {
        t.jd2().jd()
    }

    /// 実 ephemeris 用の標準 RootConfig（仕様指定値: x_tolerance_days=1e-7, max_iterations=100）。
    fn config_tight() -> RootConfig {
        RootConfig {
            x_tolerance_days: 1e-7,
            max_iterations: 100,
        }
    }

    // ---- 独立オラクル（テスト側で同式を再実装。実装の内部関数・閾値には依存しない）----

    /// 太陽の地心距離 \[km\]（合時刻の見かけ CIRS 位置のノルム。CIRS は回転ゆえ距離保存）。
    fn d_sun_km(t: TtInstant) -> f64 {
        sun_apparent_cirs(t).norm()
    }

    /// 月の地心距離 \[km\]。
    fn d_moon_km(t: TtInstant) -> f64 {
        moon_apparent_cirs(t).norm()
    }

    /// 太陽視半径 s_sun = asin(R_sun / d_sun) \[rad\]。
    fn sun_semidiameter_rad(t: TtInstant) -> f64 {
        (SOLAR_RADIUS_KM / d_sun_km(t)).clamp(-1.0, 1.0).asin()
    }

    /// 月視半径 s_moon = asin(k·Re / d_moon) \[rad\]。
    fn moon_semidiameter_rad(t: TtInstant) -> f64 {
        (MOON_SUN_SIZE_RATIO_K * EARTH_EQUATORIAL_RADIUS_KM / d_moon_km(t))
            .clamp(-1.0, 1.0)
            .asin()
    }

    /// 月の地平視差 π_moon = asin(Re / d_moon) \[rad\]（≈0.95°, 偽陰性回避に必須）。
    fn moon_parallax_rad(t: TtInstant) -> f64 {
        (EARTH_EQUATORIAL_RADIUS_KM / d_moon_km(t))
            .clamp(-1.0, 1.0)
            .asin()
    }

    /// 視半径和 s_sun + s_moon \[rad\]（独立オラクル, 契約2）。
    fn sum_semidiameters_oracle(t: TtInstant) -> f64 {
        sun_semidiameter_rad(t) + moon_semidiameter_rad(t)
    }

    /// 食限（保守マージン抜き）s_sun + s_moon + π_moon \[rad\]。判定式の下限境界。
    /// 実装はこれに保守マージンを足すので、`separation < limit_no_margin` なら必ず possible=true。
    fn eclipse_limit_no_margin_rad(t: TtInstant) -> f64 {
        sum_semidiameters_oracle(t) + moon_parallax_rad(t)
    }

    /// 影軸-地球中心の概算距離 approx_gamma（Re 無次元, 契約3）。
    /// = |S × M| / (|M − S| · Re)。S=太陽, M=月 の合時刻地心位置 \[km\]。CIRS/GCRS は回転差のみで
    /// |S×M|/|M−S| を保つため見かけ CIRS 位置で計算してよい。
    fn approx_gamma_oracle(t: TtInstant) -> f64 {
        let s: Vector3 = sun_apparent_cirs(t);
        let m: Vector3 = moon_apparent_cirs(t);
        let cross = s.cross(m).norm();
        let baseline = (m - s).norm();
        cross / (baseline * EARTH_EQUATORIAL_RADIUS_KM)
    }

    // ---- 合 fixture の構築（候補生成 → 対象日に最も近い候補 → solve_conjunction）----

    /// `target_utc_jd`（UTC-JD）の朔に対応する合を解いて返す。±12 日窓で朔候補を生成し、approx_tt が
    /// 対象日に最も近い候補を選び、黄経合を `config_tight` で解く（仕様の合の作り方）。
    fn conjunction_near(target_utc_jd: f64) -> Conjunction {
        let candidates = new_moon_candidates(utc_range(target_utc_jd - 12.0, target_utc_jd + 12.0))
            .expect("post-1972 range");
        // approx_tt(TT-JD) が対象 UTC-JD に最も近い候補を選ぶ（ΔT≈69s≈8e-4日 は朔選別に無影響）。
        let candidate: NewMoonCandidate = candidates
            .into_iter()
            .min_by(|a, b| {
                (jd_of(a.approx_tt) - target_utc_jd)
                    .abs()
                    .total_cmp(&(jd_of(b.approx_tt) - target_utc_jd).abs())
            })
            .expect("non-empty candidates near target");
        solve_conjunction(
            &candidate,
            ConjunctionKind::EclipticLongitude,
            config_tight(),
        )
        .expect("conjunction must solve near target date")
    }

    // ---- 既知の実日食 / 非日食の朔（fixtures, 出典: NASA 5000 年日食カタログ）----

    /// 既知の皆既/金環日食の new moon（UTC-JD）。すべて possible=true が期待（偽陰性ゼロの砦）。
    /// 出典: NASA Five Millennium Catalog of Solar Eclipses (eclipse.gsfc.nasa.gov)。
    /// - 2017-08-21 皆既 / 2019-07-02 皆既 / 2020-06-21 金環 / 2024-04-08 皆既。
    const KNOWN_ECLIPSE_NEW_MOONS_UTC_JD: [f64; 4] = [
        2_457_987.0, // 2017-08-21 皆既日食
        2_458_667.0, // 2019-07-02 皆既日食
        2_459_021.0, // 2020-06-21 金環日食
        2_460_408.0, // 2024-04-08 皆既日食
    ];

    /// 既知の **grazing 部分日食**（中心線を持たず半影のみ地球に触れる, gamma が食限近傍）の new moon。
    /// すべて possible=true が期待。**偽陰性ゼロの最難ケース**: 合の角距離が食限（s_sun+s_moon+π_moon, マージン
    /// 抜き）を上回りうるため、保守マージンが効いて初めて拾える＝マージンの存在を縛る（中心食では効かない）。
    /// 出典: NASA Five Millennium Catalog（部分日食）。
    /// - 2011-07-01 部分（gamma≈1.49, 極浅 grazing）/ 2018-07-13 部分（gamma≈1.35, 深い部分）。
    const KNOWN_PARTIAL_ECLIPSE_NEW_MOONS_UTC_JD: [f64; 2] = [
        2_455_743.9, // 2011-07-01 部分日食（gamma≈1.49）
        2_458_312.6, // 2018-07-13 部分日食（gamma≈1.35）
    ];

    /// 日食でない朔（possible=false 期待）。上記日食の「次の朔」＝半月分交点から離れる朔。NASA
    /// カタログにこれらの日付の日食は無い。月黄緯が大きく食帯が地球を外す（SeparationTooLarge）。
    /// - 2017-09-20 頃（2017-08-21 の次の朔）/ 2024-05-08 頃（2024-04-08 の次の朔）。
    const KNOWN_NON_ECLIPSE_NEW_MOONS_UTC_JD: [f64; 2] = [
        2_458_016.0, // 2017-09-20 頃の朔（日食なし）
        2_460_438.0, // 2024-05-08 頃の朔（日食なし）
    ];

    // ============================================================
    // 契約5（最重要）: 偽陰性ゼロ — 既知の実日食は必ず possible=true
    // ============================================================

    /// 既知の皆既/金環日食 4 件の合がすべて possible=true。1 件でも false なら fail。
    /// 実日食を 1 つも落とさないという ISSUE-018 の最重要契約の直接検証。
    #[test]
    fn known_real_eclipses_are_all_possible() {
        for &utc_jd in &KNOWN_ECLIPSE_NEW_MOONS_UTC_JD {
            let conj = conjunction_near(utc_jd);
            let assessment = assess_eclipse_possibility(&conj);
            assert!(
                assessment.possible,
                "known real eclipse at UTC-JD {utc_jd} (sep={} rad) classified as NOT possible \
                 (false negative — forbidden)",
                conj.separation.0
            );
        }
    }

    /// 偽陰性ゼロ（最難ケース）: grazing 部分日食もすべて possible=true。これらは合の角距離が食限
    /// （マージン抜き）を上回りうるため、**保守マージンが効いて初めて拾える**。マージンを縮める/符号反転
    /// する退行（偽陰性）を直接捕捉する（中心食では separation が小さくマージン無しでも拾えるため捕捉不能）。
    #[test]
    fn known_grazing_partial_eclipses_are_all_possible() {
        for &utc_jd in &KNOWN_PARTIAL_ECLIPSE_NEW_MOONS_UTC_JD {
            let conj = conjunction_near(utc_jd);
            let assessment = assess_eclipse_possibility(&conj);
            assert!(
                assessment.possible,
                "known grazing partial eclipse at UTC-JD {utc_jd} (sep={} rad, limit_no_margin={}) \
                 classified as NOT possible (false negative — margin insufficient/inverted?)",
                conj.separation.0,
                eclipse_limit_no_margin_rad(conj.time_tt)
            );
        }
    }

    /// 偽陰性ゼロの内訳: 各実日食で separation が食限（マージン抜き）を下回ること自体を独立式で確認。
    /// 実装の閾値定数に依存せず、separation < s_sun+s_moon+π_moon という物理下限が成り立つことを示す
    /// （実装はこれに保守マージンを足すだけなので、ここが成り立てば possible=true は必然）。
    #[test]
    fn real_eclipse_separation_is_below_geometric_limit() {
        for &utc_jd in &KNOWN_ECLIPSE_NEW_MOONS_UTC_JD {
            let conj = conjunction_near(utc_jd);
            let limit = eclipse_limit_no_margin_rad(conj.time_tt);
            assert!(
                conj.separation.0 < limit,
                "eclipse at UTC-JD {utc_jd}: separation {} rad not below geometric eclipse limit \
                 {limit} rad (s_sun+s_moon+π_moon)",
                conj.separation.0
            );
        }
    }

    // ============================================================
    // 契約6: 明らかな非日食の棄却 — possible=false かつ SeparationTooLarge
    // ============================================================

    /// 日食でない朔（交点から遠く月黄緯が大きい）の合は possible=false（SeparationTooLarge）。
    /// 棄却が 1 件も起きない（フィルタが素通し）退行を捕捉する。境界が緩すぎてどれも false に
    /// ならない場合はテスト失敗とし、より交点から遠い朔へ選び直す方針（仕様注記）。
    #[test]
    fn clear_non_eclipse_new_moons_are_rejected() {
        let mut rejected_count = 0usize;
        for &utc_jd in &KNOWN_NON_ECLIPSE_NEW_MOONS_UTC_JD {
            let conj = conjunction_near(utc_jd);
            let assessment = assess_eclipse_possibility(&conj);
            if !assessment.possible {
                rejected_count += 1;
                // 棄却理由は SeparationTooLarge でなければならない（契約4の possible=false 側）。
                assert_eq!(
                    assessment.reason,
                    PossibilityReason::SeparationTooLarge,
                    "non-eclipse at UTC-JD {utc_jd}: possible=false but reason != SeparationTooLarge"
                );
            }
        }
        // フィルタが少なくとも 1 件の明確な非日食を棄却すること（素通しフィルタの退行検出）。
        assert!(
            rejected_count >= 1,
            "filter rejected none of the clear non-eclipse new moons (filter passes everything?)"
        );
    }

    /// 非日食棄却の独立根拠: 棄却された朔では separation が食限（保守的に π_moon と現実的なマージン
    /// 上限を足してもなお）を上回ること。separation が十分大きい（食帯が地球を外す）→ false という
    /// 物理的妥当性をテスト側の式で確認する。マージンの厳密値は実装非依存なので、判定式の下限
    /// （s_sun+s_moon+π_moon）を「明確に」超える朔が少なくとも 1 件あることを要求する。
    #[test]
    fn at_least_one_new_moon_has_separation_clearly_above_limit() {
        // 非日食朔の中に、食限（マージン抜き）を「明確に」超える separation を持つものがあること。
        // 「明確に」= 食限の 1.5 倍以上（保守マージンが食限の数割でも possible=false に倒れる余地）。
        let clearly_above = KNOWN_NON_ECLIPSE_NEW_MOONS_UTC_JD.iter().any(|&utc_jd| {
            let conj = conjunction_near(utc_jd);
            let limit = eclipse_limit_no_margin_rad(conj.time_tt);
            conj.separation.0 > 1.5 * limit
        });
        assert!(
            clearly_above,
            "expected at least one clear non-eclipse new moon with separation > 1.5× eclipse limit"
        );
    }

    // ============================================================
    // 契約1+2: 判定式・付随量（独立オラクルと一致）
    // ============================================================

    /// 契約2: min_separation == conjunction.separation（保守的代用, 厳密一致）。
    #[test]
    fn min_separation_equals_conjunction_separation() {
        let conj = conjunction_near(KNOWN_ECLIPSE_NEW_MOONS_UTC_JD[0]); // 2017-08-21
        let assessment = assess_eclipse_possibility(&conj);
        assert!(
            close(assessment.min_separation.0, conj.separation.0, 1e-12),
            "min_separation {} != conjunction.separation {}",
            assessment.min_separation.0,
            conj.separation.0
        );
    }

    /// 契約2: sum_semidiameters == s_sun + s_moon（合時刻の距離から, 独立式と一致, tol 1e-6 rad）。
    #[test]
    fn sum_semidiameters_matches_independent_formula() {
        let conj = conjunction_near(KNOWN_ECLIPSE_NEW_MOONS_UTC_JD[0]); // 2017-08-21
        let assessment = assess_eclipse_possibility(&conj);
        let oracle = sum_semidiameters_oracle(conj.time_tt);
        assert!(
            close(assessment.sum_semidiameters.0, oracle, 1e-6),
            "sum_semidiameters {} != independent s_sun+s_moon {oracle}",
            assessment.sum_semidiameters.0
        );
    }

    /// 契約2（桁感）: sum_semidiameters が ~0.0087 rad（≈0.5°, 太陽視直径オーダー）にあること。
    /// 単位取り違え（度/ラジアン）や Re の m/km 取り違えを殺す。
    #[test]
    fn sum_semidiameters_is_about_half_degree() {
        let conj = conjunction_near(KNOWN_ECLIPSE_NEW_MOONS_UTC_JD[0]); // 2017-08-21
        let assessment = assess_eclipse_possibility(&conj);
        // 太陽+月の視半径和は ~0.0087 rad（0.5°）。幅を持たせて 0.006–0.012 rad で縛る。
        assert!(
            (0.006..=0.012).contains(&assessment.sum_semidiameters.0),
            "sum_semidiameters {} rad outside ~0.5° band (0.006–0.012 rad)",
            assessment.sum_semidiameters.0
        );
    }

    /// 契約3: approx_gamma が独立式 |S×M|/(|M−S|·Re) と一致（tol 1e-4, Re 無次元）。
    #[test]
    fn approx_gamma_matches_independent_formula() {
        let conj = conjunction_near(KNOWN_ECLIPSE_NEW_MOONS_UTC_JD[0]); // 2017-08-21
        let assessment = assess_eclipse_possibility(&conj);
        let oracle = approx_gamma_oracle(conj.time_tt);
        assert!(
            close(assessment.approx_gamma, oracle, 1e-4),
            "approx_gamma {} != independent |S×M|/(|M−S|·Re) {oracle}",
            assessment.approx_gamma
        );
    }

    /// 契約3（既知値）: 2017-08-21 で approx_gamma が NASA gamma=0.4367 近傍（0.40–0.47）。
    /// gamma の桁・符号スケール（影軸が地球中心からどれだけ外れるか）を実イベントで固定する。
    #[test]
    fn approx_gamma_for_2017_is_near_nasa_value() {
        let conj = conjunction_near(KNOWN_ECLIPSE_NEW_MOONS_UTC_JD[0]); // 2017-08-21, NASA gamma=0.4367
        let assessment = assess_eclipse_possibility(&conj);
        assert!(
            (0.40..=0.47).contains(&assessment.approx_gamma),
            "approx_gamma {} for 2017-08-21 not near NASA gamma 0.4367 (0.40–0.47)",
            assessment.approx_gamma
        );
    }

    // ============================================================
    // 契約4: reason 整合（possible ⟺ PossibleEclipse）
    // ============================================================

    /// possible ⟺ matches!(reason, PossibleEclipse)。possible=false なら必ず SeparationTooLarge。
    /// 実日食（possible=true）と非日食（possible=false 期待）の両系列で reason 整合を縛る。
    #[test]
    fn reason_is_consistent_with_possible_flag() {
        let all_utc_jds = KNOWN_ECLIPSE_NEW_MOONS_UTC_JD
            .iter()
            .chain(KNOWN_NON_ECLIPSE_NEW_MOONS_UTC_JD.iter());
        for &utc_jd in all_utc_jds {
            let conj = conjunction_near(utc_jd);
            let a = assess_eclipse_possibility(&conj);
            // possible ⟺ PossibleEclipse。
            assert_eq!(
                a.possible,
                matches!(a.reason, PossibilityReason::PossibleEclipse),
                "UTC-JD {utc_jd}: possible={} but reason={:?} (inconsistent)",
                a.possible,
                a.reason
            );
            // possible=false の唯一の理由は SeparationTooLarge。
            if !a.possible {
                assert_eq!(
                    a.reason,
                    PossibilityReason::SeparationTooLarge,
                    "UTC-JD {utc_jd}: possible=false but reason != SeparationTooLarge"
                );
            }
        }
    }

    /// 契約4: 実日食では reason==PossibleEclipse（possible=true 側の reason 取り違え検出）。
    #[test]
    fn real_eclipses_have_possible_eclipse_reason() {
        for &utc_jd in &KNOWN_ECLIPSE_NEW_MOONS_UTC_JD {
            let conj = conjunction_near(utc_jd);
            let a = assess_eclipse_possibility(&conj);
            assert_eq!(
                a.reason,
                PossibilityReason::PossibleEclipse,
                "real eclipse at UTC-JD {utc_jd}: reason != PossibleEclipse"
            );
        }
    }

    // ============================================================
    // 契約7: 保守性プロパティ — separation が小さいほど possible（単調性）
    // ============================================================

    /// 保守性: separation が食限（マージン抜き）を下回るイベントはすべて possible=true。実日食は
    /// すべてこれを満たす（契約5）。実装の閾値はこの下限以上（マージン≥0）でなければ偽陰性が出るため、
    /// 「separation < 食限 ⟹ possible」を全実日食で要求する＝判定式の保守側の下限を縛る。
    #[test]
    fn separation_below_limit_implies_possible() {
        for &utc_jd in &KNOWN_ECLIPSE_NEW_MOONS_UTC_JD {
            let conj = conjunction_near(utc_jd);
            let limit = eclipse_limit_no_margin_rad(conj.time_tt);
            let a = assess_eclipse_possibility(&conj);
            // 前提（実日食は食限内）を確認した上で、possible=true を要求する。
            assert!(
                conj.separation.0 < limit,
                "precondition: eclipse at {utc_jd} separation not below limit"
            );
            assert!(
                a.possible,
                "conservativeness: separation {} < limit {limit} at {utc_jd} but possible=false",
                conj.separation.0
            );
        }
    }

    /// 保守性（単調性, メタ）: より小さい separation を持つイベントが possible でないことは起き得ない。
    /// 実日食（最小級の separation）と非日食（より大きい separation）を比べ、実日食側が possible で
    /// あること（separation 順序と possible の整合）を一対で確認する。判定式が separation について
    /// 単調（小さいほど possible 寄り）であることの観測的検証。
    #[test]
    fn smaller_separation_event_is_never_less_possible() {
        let eclipse = conjunction_near(KNOWN_ECLIPSE_NEW_MOONS_UTC_JD[0]); // 2017-08-21（小 separation）
        let non_eclipse = conjunction_near(KNOWN_NON_ECLIPSE_NEW_MOONS_UTC_JD[0]); // 大 separation
                                                                                   // 前提: 実日食の separation は非日食より小さい（イベント選別が正しい）。
        assert!(
            eclipse.separation.0 < non_eclipse.separation.0,
            "precondition: eclipse separation {} not smaller than non-eclipse {}",
            eclipse.separation.0,
            non_eclipse.separation.0
        );
        let a_eclipse = assess_eclipse_possibility(&eclipse);
        // 小さい separation 側（実日食）は必ず possible（単調性）。
        assert!(
            a_eclipse.possible,
            "monotonicity: smaller-separation eclipse must be possible"
        );
    }

    // ============================================================
    // 健全性: 全フィールド有限・符号
    // ============================================================

    /// 全フィールドが有限、approx_gamma≥0、sum_semidiameters>0、min_separation≥0。
    /// reason と conjunction フィールドも整合（conjunction はそのまま埋め込まれる）。
    #[test]
    fn all_fields_are_finite_and_well_signed() {
        let all_utc_jds = KNOWN_ECLIPSE_NEW_MOONS_UTC_JD
            .iter()
            .chain(KNOWN_NON_ECLIPSE_NEW_MOONS_UTC_JD.iter());
        for &utc_jd in all_utc_jds {
            let conj = conjunction_near(utc_jd);
            let a = assess_eclipse_possibility(&conj);
            assert!(
                a.min_separation.0.is_finite(),
                "UTC-JD {utc_jd}: min_separation non-finite"
            );
            assert!(
                a.sum_semidiameters.0.is_finite(),
                "UTC-JD {utc_jd}: sum_semidiameters non-finite"
            );
            assert!(
                a.approx_gamma.is_finite(),
                "UTC-JD {utc_jd}: approx_gamma non-finite"
            );
            // approx_gamma は距離（Re 無次元）なので非負。
            assert!(
                a.approx_gamma >= 0.0,
                "UTC-JD {utc_jd}: approx_gamma {} negative",
                a.approx_gamma
            );
            // 視半径和は正（太陽・月とも有限距離で視半径>0）。
            assert!(
                a.sum_semidiameters.0 > 0.0,
                "UTC-JD {utc_jd}: sum_semidiameters {} not positive",
                a.sum_semidiameters.0
            );
            // min_separation は角距離なので非負。
            assert!(
                a.min_separation.0 >= 0.0,
                "UTC-JD {utc_jd}: min_separation {} negative",
                a.min_separation.0
            );
        }
    }

    /// 健全性: 埋め込まれた conjunction が入力と一致（separation/kind/time_tt をそのまま保持）。
    /// フィルタが合の情報を破壊・差し替えしないこと。
    #[test]
    fn embedded_conjunction_matches_input() {
        let conj = conjunction_near(KNOWN_ECLIPSE_NEW_MOONS_UTC_JD[0]); // 2017-08-21
        let a = assess_eclipse_possibility(&conj);
        assert_eq!(
            a.conjunction, conj,
            "embedded conjunction differs from the input conjunction"
        );
    }
}
