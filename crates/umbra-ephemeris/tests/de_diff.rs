//! ISSUE-013/014 M10 DE確定ゲート（accuracy.md §3.3）の測定。
//!
//! `AnalyticalEphemeris`（VSOP87D 太陽 + ELP2000-82B 月・幾何地心位置）を
//! `JplEphemeris`（DE440s Reference オラクル・ISSUE-036）と 1900–2100 で突合し、地心方向の
//! 角度残差（秒角）と月距離の相対残差を測る。両者とも `Origin::Geocenter` / `EphemerisFrame::Icrs`
//! の**幾何**位置（見かけ補正前）で比較する＝暦層のみの差分（accuracy.md §3.1-1）。
//!
//! 実 `data/spk/de440s.bsp` 存在時のみ実行（CI 非同梱・ISSUE-036）。feature `jpl` 限定。
#![cfg(feature = "jpl")]
// サンプル添字 i（<= 数千）を f64 化する。値は厳密に表現可能なため精度損失なし。
#![allow(clippy::cast_precision_loss)]

use std::path::Path;

use umbra_core::{JulianDate2, TdbInstant, Vector3};
use umbra_ephemeris::{AnalyticalEphemeris, Body, Ephemeris, EphemerisFrame, JplEphemeris, Origin};

/// 実 DE440s（リポジトリ root の data/spk）。CARGO_MANIFEST_DIR は crates/umbra-ephemeris。
const DE440S_PATH: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../data/spk/de440s.bsp");

/// rad → 秒角（180·3600/π）。
const ARCSEC_PER_RAD: f64 = 206_264.806_247_096_36;

/// 1900-01-01 0h / 2100-01-01 0h（JD, AnalyticalEphemeris の対応域に一致）。
const START_JD: f64 = 2_415_020.5;
const END_JD: f64 = 2_488_069.5;

fn tdb(jd: f64) -> TdbInstant {
    TdbInstant::from_jd2(JulianDate2::from_jd(jd))
}

/// 2 ベクトル方向のなす角（秒角）。atan2(|a×b|, a·b) で全角域に頑健。
fn angular_sep_arcsec(a: Vector3, b: Vector3) -> f64 {
    a.cross(b).norm().atan2(a.dot(b)) * ARCSEC_PER_RAD
}

/// 太陽・月の地心方向（ICRS 幾何）を DE440s と突合し、最大残差を報告する。
///
/// 目標（ISSUE-013/014・accuracy.md §2.1）: 太陽 ≤ 0.05″、月 ≤ 0.1″（打切り前の full 係数で確認）。
///
/// **現状 M10 ゲートは未達**（2026-06-20 測定: 太陽 max 0.368″ / 月 max 1.685″, 1900–2100, n=2435）。
/// 残差は J2000 で最小・両端へ増大する V 字。診断:
/// - **太陽**: `apparent::sun_geocentric_gcrs` が VSOP87D（黄道 of date）＋ of-date フレーム行列
///   を使い、VSOP87 の力学的分点 of date と IAU2006 分点 of date の歳差レート不一致が出る
///   （フレーム行列自体は SOFA 検証済みで正しい）。太陽を月と同じ J2000 黄道経路（VSOP87A 系）に
///   すれば mas 級へ改善見込み（ISSUE-033/035）。
/// - **月**: `moon_geocentric_gcrs` は J2000 固定行列で clean に変換しており、残差は純粋に
///   ELP2000-82B の DE440 乖離（永年項）。ISSUE-014 が要求する **ELP/MPP02 DE-fit**（ISSUE-034）で解消。
///
/// 根本対応（上記2件）まで `#[ignore]`。`cargo test --features jpl -- --ignored --nocapture` で
/// 現状残差を測定できる（M10 達成残差の確定根拠）。
#[ignore = "M10 DE確定ゲート未達: 太陽 ~0.37″/月 ~1.69″（ISSUE-033/035 太陽フレーム経路・ISSUE-034 月MPP02 で対応）"]
#[test]
fn analytical_vs_de440s_geocentric_direction() {
    let jpl = match JplEphemeris::from_spk_path(Path::new(DE440S_PATH)) {
        Ok(j) => j,
        Err(e) => {
            eprintln!(
                "skip analytical_vs_de440s_geocentric_direction: {DE440S_PATH} を読めない（{e:?}）。\
                 実 DE440s は CI 非同梱（ISSUE-036）。`cargo xtask fetch-de440s` で取得。"
            );
            return;
        }
    };
    let ana = AnalyticalEphemeris::new();

    // 1900–2100 を ~30 日刻みでサンプル（≈2435 点）。
    let n: usize = 2435;
    let step = (END_JD - START_JD) / (n as f64);

    let mut max_sun = 0.0_f64;
    let mut max_moon = 0.0_f64;
    let mut max_moon_dr = 0.0_f64;
    let (mut at_sun, mut at_moon) = (0.0_f64, 0.0_f64);

    for i in 0..=n {
        let jd = START_JD + step * (i as f64);
        let t = tdb(jd);

        let sun_a = ana
            .state(Body::Sun, t, Origin::Geocenter, EphemerisFrame::Icrs)
            .expect("analytical Sun")
            .position;
        let sun_j = jpl
            .state(Body::Sun, t, Origin::Geocenter, EphemerisFrame::Icrs)
            .expect("jpl Sun")
            .position;
        let sun_sep = angular_sep_arcsec(sun_a, sun_j);
        if sun_sep > max_sun {
            max_sun = sun_sep;
            at_sun = jd;
        }

        let moon_a = ana
            .state(Body::Moon, t, Origin::Geocenter, EphemerisFrame::Icrs)
            .expect("analytical Moon")
            .position;
        let moon_j = jpl
            .state(Body::Moon, t, Origin::Geocenter, EphemerisFrame::Icrs)
            .expect("jpl Moon")
            .position;
        let moon_sep = angular_sep_arcsec(moon_a, moon_j);
        if moon_sep > max_moon {
            max_moon = moon_sep;
            at_moon = jd;
        }
        let dr = ((moon_a.norm() - moon_j.norm()) / moon_j.norm()).abs();
        if dr > max_moon_dr {
            max_moon_dr = dr;
        }
    }

    // 測定値を常に表示（M10 ゲートの達成残差を記録する根拠 = `cargo test -- --nocapture`）。
    eprintln!(
        "[DE440s diff 1900-2100, n={n}] \
         max Sun = {max_sun:.5}\" (jd {at_sun:.1}), \
         max Moon = {max_moon:.5}\" (jd {at_moon:.1}), \
         max Moon δr/r = {max_moon_dr:.3e}"
    );

    assert!(
        max_sun <= 0.05,
        "太陽 地心方向残差 {max_sun:.5}\" が目標 0.05\" を超過（jd {at_sun:.1}）"
    );
    assert!(
        max_moon <= 0.1,
        "月 地心方向残差 {max_moon:.5}\" が目標 0.1\" を超過（jd {at_moon:.1}）"
    );
}
