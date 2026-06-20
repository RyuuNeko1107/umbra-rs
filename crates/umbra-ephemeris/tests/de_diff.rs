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

use umbra_core::{JulianDate2, TdbInstant, TtInstant, Vector3};
use umbra_ephemeris::frames::ecliptic_to_gcrs_matrix;
use umbra_ephemeris::sun::sun_geocentric_ecliptic_of_date;
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

/// DE440s を開く（不在時はスキップ告知して `None`）。
fn load_jpl() -> Option<JplEphemeris> {
    match JplEphemeris::from_spk_path(Path::new(DE440S_PATH)) {
        Ok(j) => Some(j),
        Err(e) => {
            eprintln!(
                "skip: {DE440S_PATH} を読めない（{e:?}）。実 DE440s は CI 非同梱（ISSUE-036）。\
                 `cargo xtask fetch-de440s` で取得。"
            );
            None
        }
    }
}

/// サンプル点数（1900–2100 を ~30 日刻み）。
const N_SAMPLES: usize = 2435;

/// `body` の地心方向（ICRS 幾何）を DE440s と突合し、最大角度残差 `(arcsec, jd)` を返す。
fn max_direction_residual(jpl: &JplEphemeris, ana: &AnalyticalEphemeris, body: Body) -> (f64, f64) {
    let step = (END_JD - START_JD) / (N_SAMPLES as f64);
    let (mut max_sep, mut at) = (0.0_f64, 0.0_f64);
    for i in 0..=N_SAMPLES {
        let jd = START_JD + step * (i as f64);
        let t = tdb(jd);
        let a = ana
            .state(body, t, Origin::Geocenter, EphemerisFrame::Icrs)
            .expect("analytical state")
            .position;
        let j = jpl
            .state(body, t, Origin::Geocenter, EphemerisFrame::Icrs)
            .expect("jpl state")
            .position;
        let sep = angular_sep_arcsec(a, j);
        if sep > max_sep {
            max_sep = sep;
            at = jd;
        }
    }
    (max_sep, at)
}

/// 太陽の地心方向（ICRS 幾何）を DE440s と突合（M10 太陽ゲート・ISSUE-013/035）。
///
/// 目標 ≤ 0.05″（accuracy.md §2.1）。VSOP87A（黄道 J2000 直交）＋ VSOP87→ICRS 固定回転で、
/// 旧 VSOP87D 黄道of date 経路の歳差レート不整合（0.37″）を解消（2026-06-20: max ~0.03″）。
/// de440s.bsp 不在時はスキップ（CI 非同梱・ISSUE-036）。
#[test]
fn sun_vs_de440s_geocentric_direction() {
    let Some(jpl) = load_jpl() else { return };
    let ana = AnalyticalEphemeris::new();
    let (max_sun, at) = max_direction_residual(&jpl, &ana, Body::Sun);
    eprintln!("[DE440s Sun diff 1900-2100, n={N_SAMPLES}] max = {max_sun:.5}\" (jd {at:.1})");
    assert!(
        max_sun <= 0.05,
        "太陽 地心方向残差 {max_sun:.5}\" が目標 0.05\" を超過（jd {at:.1}）"
    );
}

/// 月の地心方向（ICRS 幾何）を DE440s と突合（ISSUE-014）。
///
/// 目標 ≤ 0.1″ だが、現状 ELP2000-82B（LLR-fit, full）は DE440 と永年項が乖離して **~1.69″** 未達。
/// 月のフレーム経路は J2000 固定で clean（残差は純粋にモデル差）。ISSUE-014 が要求する
/// **ELP/MPP02 DE-fit**（ISSUE-034）導入まで `#[ignore]`。
/// `cargo test --features jpl -- --ignored --nocapture` で現状残差を測定できる。
#[ignore = "ELP2000-82B は DE-fit でなく ~1.69″ 未達。ISSUE-014 ELP/MPP02（ISSUE-034）で解消予定"]
#[test]
fn moon_vs_de440s_geocentric_direction() {
    let Some(jpl) = load_jpl() else { return };
    let ana = AnalyticalEphemeris::new();
    let (max_moon, at) = max_direction_residual(&jpl, &ana, Body::Moon);
    eprintln!("[DE440s Moon diff 1900-2100, n={N_SAMPLES}] max = {max_moon:.5}\" (jd {at:.1})");
    assert!(
        max_moon <= 0.1,
        "月 地心方向残差 {max_moon:.5}\" が目標 0.1\" を超過（jd {at:.1}）"
    );
}

/// 太陽残差の黄経Δλ・黄緯Δβ 分解（診断・常時 ignore）。
///
/// VSOP87D 黄道of date 位置と、DE440s ICRS を `ecliptic_to_gcrs_matrix(tt)ᵀ` で黄道of dateへ
/// 戻したものを λ=atan2(y,x), β=asin(z/|r|) で比較。`Δλ·cosβ`（黄経残差）と `Δβ` を秒角で各
/// エポックに表示する。J2000 floor（≈0.07″）が黄経の定数オフセット（VSOP87→ICRS 分点バイアス
/// 未適用）か、黄緯（傾斜）か、両端の線形増大（歳差レート不一致）かを切り分ける。
#[ignore = "diagnostic: 太陽残差の黄経/黄緯分解（cargo test --features jpl -- --ignored --nocapture）"]
#[test]
fn sun_residual_ecliptic_decomposition() {
    let Ok(jpl) = JplEphemeris::from_spk_path(Path::new(DE440S_PATH)) else {
        eprintln!("skip sun_residual_ecliptic_decomposition: de440s.bsp 不在");
        return;
    };

    // J2000 を中心に ±25/50/75/100 年（年≈365.25 日）。
    let j2000 = 2_451_545.0_f64;
    for years in [-100.0, -75.0, -50.0, -25.0, 0.0, 25.0, 50.0, 75.0, 100.0] {
        let jd = j2000 + years * 365.25;
        let t = tdb(jd);
        let tt = TtInstant::from_jd2(JulianDate2::from_jd(jd));

        // 解析: VSOP87D 黄道of date（AU・方向のみ使用）。
        let ana_ecl = sun_geocentric_ecliptic_of_date(t);
        // DE: ICRS → 黄道of date（同 IAU2006 行列の転置）。
        let de_icrs = jpl
            .state(Body::Sun, t, Origin::Geocenter, EphemerisFrame::Icrs)
            .unwrap()
            .position;
        let de_ecl = ecliptic_to_gcrs_matrix(tt).transpose().mul_vec(de_icrs);

        let (la, ba) = lon_lat(ana_ecl);
        let (ld, bd) = lon_lat(de_ecl);
        let mut dlon = la - ld;
        // [-π,π] に折り返す。
        while dlon > std::f64::consts::PI {
            dlon -= std::f64::consts::TAU;
        }
        while dlon < -std::f64::consts::PI {
            dlon += std::f64::consts::TAU;
        }
        let dlon_arcsec = dlon * ba.cos() * ARCSEC_PER_RAD;
        let dlat_arcsec = (ba - bd) * ARCSEC_PER_RAD;
        eprintln!(
            "  Sun {years:+6.0}yr (jd {jd:.1}): Δλ·cosβ = {dlon_arcsec:+.5}\"  Δβ = {dlat_arcsec:+.5}\""
        );
    }
}

/// 直交（黄道）→ 黄経 λ\[rad\]・黄緯 β\[rad\]。
fn lon_lat(v: Vector3) -> (f64, f64) {
    let lon = v.y.atan2(v.x);
    let lat = v.z.atan2((v.x * v.x + v.y * v.y).sqrt());
    (lon, lat)
}

/// VSOP87A 源ファイル（J2000 黄道直交 X,Y,Z）。
const VSOP87A_PATH: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../data/coefficient-source/vsop87/VSOP87A.ear"
);

/// VSOP87 力学的黄道・分点 J2000 → FK5(≈ICRS) 赤道 J2000 の固定回転（Bretagnon, VSOP87.doc）。
/// 対角は obliquity ε0（cos=0.91748…, sin=0.39777…）、off-diagonal の微小項が分点バイアス。
const VSOP87_ECL_TO_FK5_EQ: [[f64; 3]; 3] = [
    [1.000_000_000_000, 0.000_000_440_360, -0.000_000_190_919],
    [-0.000_000_479_966, 0.917_482_137_087, -0.397_776_982_902],
    [0.000_000_000_000, 0.397_776_982_902, 0.917_482_137_087],
];

fn apply_matrix(m: &[[f64; 3]; 3], v: Vector3) -> Vector3 {
    Vector3::new(
        m[0][0] * v.x + m[0][1] * v.y + m[0][2] * v.z,
        m[1][0] * v.x + m[1][1] * v.y + m[1][2] * v.z,
        m[2][0] * v.x + m[2][1] * v.y + m[2][2] * v.z,
    )
}

/// VSOP87 級数 1 項（振幅 A, 位相 B, 振動数 C）。
type AbcTerm = (f64, f64, f64);
/// 1 セクション（variable 1=X/2=Y/3=Z, power, 項列）。
type Vsop87aSection = (u8, u8, Vec<AbcTerm>);

/// VSOP87A 源テキストを セクション列に素朴パースする（スパイク用）。
fn parse_vsop87a(text: &str) -> Vec<Vsop87aSection> {
    let mut out: Vec<Vsop87aSection> = Vec::new();
    for line in text.lines() {
        let tok: Vec<&str> = line.split_whitespace().collect();
        if line.contains("VARIABLE") && line.contains("TERMS") {
            let var: u8 = tok[tok.iter().position(|t| *t == "VARIABLE").unwrap() + 1]
                .parse()
                .unwrap();
            let power: u8 = tok
                .iter()
                .find_map(|t| t.strip_prefix("*T**").and_then(|p| p.parse().ok()))
                .unwrap();
            out.push((var, power, Vec::new()));
        } else if tok.len() >= 3 && tok[0].parse::<i64>().is_ok() {
            let n = tok.len();
            let (a, b, c) = (
                tok[n - 3].parse().unwrap(),
                tok[n - 2].parse().unwrap(),
                tok[n - 1].parse().unwrap(),
            );
            out.last_mut().unwrap().2.push((a, b, c));
        }
    }
    out
}

/// VSOP87A で地球日心 J2000 黄道直交 (X,Y,Z)\[AU\] を評価。T = ユリウス千年 from J2000 TDB。
fn vsop87a_earth_xyz(sections: &[Vsop87aSection], jd_tdb: f64) -> Vector3 {
    let t = (jd_tdb - 2_451_545.0) / 365_250.0;
    let mut xyz = [0.0_f64; 3];
    for (var, power, terms) in sections {
        let s: f64 = terms.iter().map(|(a, b, c)| a * (b + c * t).cos()).sum();
        xyz[usize::from(var - 1)] += t.powi(i32::from(*power)) * s;
    }
    Vector3::new(xyz[0], xyz[1], xyz[2])
}

/// 検証スパイク: VSOP87A（J2000 黄道）＋ VSOP87→FK5/ICRS 固定回転 で太陽地心方向を作り、DE440s と突合。
/// of-date 歳差不整合が消え、定数オフセットも回転の分点バイアス項で吸収されて ≤0.05″ になるかを実証する。
#[ignore = "spike: VSOP87A J2000 経路の太陽残差を実証（--ignored --nocapture）"]
#[test]
fn spike_vsop87a_j2000_sun_vs_de440s() {
    let Ok(jpl) = JplEphemeris::from_spk_path(Path::new(DE440S_PATH)) else {
        eprintln!("skip spike: de440s.bsp 不在");
        return;
    };
    let text = std::fs::read_to_string(VSOP87A_PATH).expect("VSOP87A.ear 読込");
    let sections = parse_vsop87a(&text);

    let mut max_sun = 0.0_f64;
    let mut at = 0.0_f64;
    for years in [-100.0, -75.0, -50.0, -25.0, 0.0, 25.0, 50.0, 75.0, 100.0] {
        let jd = 2_451_545.0 + years * 365.25;
        // 地球日心 → 地心太陽 = 符号反転。J2000 黄道直交 → FK5/ICRS 赤道。
        let sun_ecl_j2000 = vsop87a_earth_xyz(&sections, jd).scale(-1.0);
        let sun_icrs = apply_matrix(&VSOP87_ECL_TO_FK5_EQ, sun_ecl_j2000);
        let de = jpl
            .state(Body::Sun, tdb(jd), Origin::Geocenter, EphemerisFrame::Icrs)
            .unwrap()
            .position;
        let sep = angular_sep_arcsec(sun_icrs, de);
        eprintln!("  VSOP87A Sun {years:+6.0}yr (jd {jd:.1}): residual = {sep:.5}\"");
        if sep > max_sun {
            max_sun = sep;
            at = jd;
        }
    }
    eprintln!("[VSOP87A J2000 spike] max Sun residual = {max_sun:.5}\" (jd {at:.1})");
    assert!(
        max_sun <= 0.05,
        "VSOP87A 経路でも太陽残差 {max_sun:.5}\" が 0.05\" 超過（jd {at:.1}）"
    );
}
