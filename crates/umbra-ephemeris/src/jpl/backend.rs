//! `JplEphemeris` — JPL DE（DE440s SPK）Reference バックエンド（ISSUE-036 S3・feature `jpl`）。
//!
//! S1（[`crate::jpl::daf`] 構造解析）と S2（[`crate::jpl::eval`] type2 評価）を配線し、
//! [`Ephemeris`] trait を実装する。SPK native は ICRF(≈ICRS)/SSB 原点。`EphemerisFrame::Icrs`
//! を native とし、`Origin::Geocenter` は body 差で算出する（地心 Sun = 10−(3+399)、
//! 地心 Moon = 301−399。data/spk/PROVENANCE.md）。`EphemerisFrame::EclipticOfDate` は
//! ISSUE-035 の変換経由（Reference は ICRS 直接利用が主）のため本バックエンドでは未提供
//! （[`EphemerisError::DataUnavailable`]）。
//!
//! DE データは crate 非同梱（利用者が `.bsp` を任意取得し `from_spk_path` で読む。
//! data-sources §2.3）。時刻は TDB 基準（SPK の ET=TDB と直接対応, conventions §6）。

use std::path::Path;

use umbra_core::constants::J2000_JD;
use umbra_core::{JulianDate2, TdbInstant, TimeRange, Vector3};

use crate::ephemeris::{
    Body, Ephemeris, EphemerisError, EphemerisFrame, EphemerisMetadata, Origin, StateVector,
};
use crate::jpl::daf::{parse_spk_segments, SpkSegment};
use crate::jpl::eval::eval_type2;

/// 1 日の秒数（ET 秒 ↔ 日）。J2000 基準は [`J2000_JD`]（umbra_core の単一権威定数）。
const SECONDS_PER_DAY: f64 = 86_400.0;

// 参照する NAIF body ID（data/spk/PROVENANCE.md）。
const ID_SSB: i32 = 0;
const ID_EMB: i32 = 3;
const ID_SUN: i32 = 10;
const ID_MOON: i32 = 301;
const ID_EARTH: i32 = 399;

/// SSB 基準位置を組むのに必要なセグメント (target, center)。Earth/Moon は EMB 経由で合成する。
/// 各 (target, center) は単一 ET スパン前提（DE440s 専用。複数スパン分割 SPK は非対応）。
const REQUIRED_SEGMENTS: [(i32, i32); 4] = [
    (ID_SUN, ID_SSB),
    (ID_EMB, ID_SSB),
    (ID_EARTH, ID_EMB),
    (ID_MOON, ID_EMB),
];

/// ET 秒（J2000 基準, TDB）から `TdbInstant` を作る。
fn tdb_from_et(et: f64) -> TdbInstant {
    TdbInstant::from_jd2(JulianDate2::new(J2000_JD, et / SECONDS_PER_DAY))
}

/// `TdbInstant` を ET 秒（J2000 基準, TDB）へ。2 部 JD の差分で精度を保つ。
fn et_from_tdb(time: TdbInstant) -> f64 {
    time.jd2().days_since(JulianDate2::new(J2000_JD, 0.0)) * SECONDS_PER_DAY
}

/// JPL DE（DE440s SPK）Reference バックエンド。
///
/// `from_spk_path` で SPK を読み、Sun/Moon/Earth/EMB のセグメントを保持する。`state` は
/// ICRS native／SSB・Geocenter 原点に対応（差分テストの第一義オラクル, accuracy.md §3.1）。
#[derive(Debug)]
pub struct JplEphemeris {
    /// SPK バイト列（DAF ワードアドレスで参照）。
    bytes: Vec<u8>,
    /// 解析済みセグメント記述子（target/center/addr…）。
    segments: Vec<SpkSegment>,
    /// 必要セグメントの被覆 ET 区間の交差（開始秒, J2000 基準）。
    coverage_start_et: f64,
    /// 同上（終了秒）。
    coverage_end_et: f64,
}

impl JplEphemeris {
    /// SPK ファイル（`.bsp`）を読み込み、必要セグメント（Sun/Moon/Earth/EMB）を解析する。
    ///
    /// IO 失敗・形式不正は [`EphemerisError::MalformedSpk`]、必要セグメント欠落は
    /// [`EphemerisError::DataUnavailable`]。被覆区間は必要セグメント被覆の交差。
    pub fn from_spk_path(path: &Path) -> Result<Self, EphemerisError> {
        let bytes = std::fs::read(path).map_err(|e| {
            EphemerisError::MalformedSpk(format!("failed to read SPK {}: {e}", path.display()))
        })?;
        let segments = parse_spk_segments(&bytes)?;

        // 必要セグメントの存在を確認し、被覆 ET 区間の交差を取る（全 body が揃う範囲）。
        let mut coverage_start_et = f64::NEG_INFINITY;
        let mut coverage_end_et = f64::INFINITY;
        for (target, center) in REQUIRED_SEGMENTS {
            let seg = segments
                .iter()
                .find(|s| s.target == target && s.center == center)
                .ok_or(EphemerisError::DataUnavailable)?;
            coverage_start_et = coverage_start_et.max(seg.start_et);
            coverage_end_et = coverage_end_et.min(seg.end_et);
        }

        Ok(Self {
            bytes,
            segments,
            coverage_start_et,
            coverage_end_et,
        })
    }

    /// (target, center) セグメントを `et` で評価し、位置/速度（km, km/s）を返す。
    /// セグメント欠落は `DataUnavailable`、被覆外は `OutOfSupportedRange`（eval_type2 由来）。
    fn eval(
        &self,
        target: i32,
        center: i32,
        et: f64,
    ) -> Result<(Vector3, Vector3), EphemerisError> {
        let seg = self
            .segments
            .iter()
            .find(|s| s.target == target && s.center == center)
            .ok_or(EphemerisError::DataUnavailable)?;
        let s = eval_type2(&self.bytes, seg, et)?;
        Ok((
            Vector3::new(s.position[0], s.position[1], s.position[2]),
            Vector3::new(s.velocity[0], s.velocity[1], s.velocity[2]),
        ))
    }

    /// `body` の SSB 基準 位置/速度（km, km/s）。Earth/Moon は EMB 経由で合成する。
    fn ssb_state(&self, body: Body, et: f64) -> Result<(Vector3, Vector3), EphemerisError> {
        match body {
            Body::Sun => self.eval(ID_SUN, ID_SSB, et),
            Body::EarthMoonBarycenter => self.eval(ID_EMB, ID_SSB, et),
            Body::Earth => {
                let (emb_p, emb_v) = self.eval(ID_EMB, ID_SSB, et)?;
                let (rel_p, rel_v) = self.eval(ID_EARTH, ID_EMB, et)?;
                Ok((emb_p + rel_p, emb_v + rel_v))
            }
            Body::Moon => {
                let (emb_p, emb_v) = self.eval(ID_EMB, ID_SSB, et)?;
                let (rel_p, rel_v) = self.eval(ID_MOON, ID_EMB, et)?;
                Ok((emb_p + rel_p, emb_v + rel_v))
            }
        }
    }
}

impl Ephemeris for JplEphemeris {
    fn state(
        &self,
        body: Body,
        time: TdbInstant,
        origin: Origin,
        frame: EphemerisFrame,
    ) -> Result<StateVector, EphemerisError> {
        // EclipticOfDate は本バックエンド未提供（ISSUE-035 変換経由・Reference は ICRS 直接利用が主）。
        if matches!(frame, EphemerisFrame::EclipticOfDate) {
            return Err(EphemerisError::DataUnavailable);
        }
        // 地心の地球は原点。時刻に依らず厳密 0（SSB 差の丸め残差を出さない）。
        if matches!(body, Body::Earth) && matches!(origin, Origin::Geocenter) {
            return Ok(StateVector {
                position: Vector3::ZERO,
                velocity: Some(Vector3::ZERO),
            });
        }

        let et = et_from_tdb(time);
        let (pos, vel) = self.ssb_state(body, et)?;
        let (position, velocity) = match origin {
            Origin::SolarSystemBarycenter => (pos, vel),
            Origin::Geocenter => {
                let (earth_p, earth_v) = self.ssb_state(Body::Earth, et)?;
                (pos - earth_p, vel - earth_v)
            }
        };
        Ok(StateVector {
            position,
            velocity: Some(velocity),
        })
    }

    fn supported_range(&self) -> TimeRange<TdbInstant> {
        TimeRange {
            start: tdb_from_et(self.coverage_start_et),
            end: tdb_from_et(self.coverage_end_et),
        }
    }

    fn metadata(&self) -> EphemerisMetadata {
        EphemerisMetadata {
            model: "JPL DE440s".to_string(),
            version: "DE440 短期版（de440s.bsp, ET 1849–2150）".to_string(),
            source: "NASA/JPL/NAIF; Park, Folkner, Williams, Boggs (2021), AJ 161:105, \
                     DOI:10.3847/1538-3881/abd414"
                .to_string(),
            license: "JPL/Caltech・非同梱・任意DL（data-sources §2.3/§6）".to_string(),
            supported: self.supported_range(),
            // Reference オラクル: Chebyshev 評価の数値誤差のみ（暦由来の許容誤差ではない）。
            max_residual_arcsec: 0.0,
        }
    }
}

#[cfg(test)]
mod tests {
    // 位置/速度と SPICE 基準を並列添字（成分 i）で比較するため許可する。
    #![allow(clippy::needless_range_loop)]

    use super::*;
    use umbra_core::{JulianDate2, TdbInstant, Vector3};

    // ==================================================================
    // 確定仕様（挙動契約）の要点（テスト設計の根拠）:
    //   時刻系: SPK の ET=TDB 秒（J2000=JD 2451545.0 基準）。ある ET 秒 et の
    //     TdbInstant は from_jd2(JulianDate2::new(2451545.0, et/86400.0))。
    //   座標/原点（NAIF SPK native = ICRF≈ICRS / SSB 原点。PROVENANCE.md）:
    //     NAIF ID: Sun=10, Moon=301, Earth=399, EMB=3, SSB=0。
    //     Origin::SolarSystemBarycenter → body の SSB 基準位置。
    //     Origin::Geocenter → body の SSB 基準 − Earth_SSB
    //       （地心 Sun=10−(3+399)、地心 Moon=301−399）。
    //     Body::Earth + Origin::Geocenter → 原点（位置 ZERO・速度 Some(ZERO)）。
    //   フレーム: Icrs は SPK native そのまま。EclipticOfDate は DataUnavailable。
    //   速度: JPL は Chebyshev 微分で提供するため velocity は Some(_)。
    //   範囲: 被覆外時刻は OutOfSupportedRange。supported_range() は DE440s 被覆
    //     （おおよそ 1849-12-26〜2150-01-21）。
    // ==================================================================

    /// J2000 基準 JD（SPK の ET=TDB 秒の基準）。
    const J2000_JD: f64 = 2451545.0;
    /// 1 日の秒数。
    const SECONDS_PER_DAY: f64 = 86400.0;

    /// ET 秒（J2000 基準, TDB）から `TdbInstant` を作る（確定仕様の構築式）。
    fn tdb_from_et(et: f64) -> TdbInstant {
        TdbInstant::from_jd2(JulianDate2::new(J2000_JD, et / SECONDS_PER_DAY))
    }

    /// 実 DE440s（data/spk/de440s.bsp）の絶対パス。CARGO_MANIFEST_DIR は
    /// crates/umbra-ephemeris。リポジトリ root の data/spk に置く。
    const DE440S_PATH: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../data/spk/de440s.bsp");

    /// 実 DE440s から `JplEphemeris` を構築する。不在なら告知して None（テストは早期 return）。
    fn load_de440s() -> Option<JplEphemeris> {
        let path = std::path::Path::new(DE440S_PATH);
        if !path.exists() {
            eprintln!(
                "skip: {DE440S_PATH} が存在しない。実 DE440s は CI 非同梱（ISSUE-036）。\
                 `cargo xtask fetch-de440s` で取得。"
            );
            return None;
        }
        match JplEphemeris::from_spk_path(path) {
            Ok(eph) => Some(eph),
            Err(e) => panic!("実 DE440s の解析に失敗: {e:?}"),
        }
    }

    /// 位置成分（km）の絶対許容。数 m 級 = 1e-2 km。
    const POS_TOL_KM: f64 = 1.0e-2;
    /// 速度成分（km/s）の絶対許容。mm/s 級 = 1e-6 km/s。
    const VEL_TOL_KMS: f64 = 1.0e-6;

    /// `StateVector` を SPICE state\[6\]（位置 km・速度 km/s）と成分照合する。
    /// 速度は Some であることを要求する（JPL は Chebyshev 微分で速度を提供）。
    fn assert_state_matches_spice(actual: &StateVector, want: [f64; 6], label: &str) {
        let pos = actual.position;
        let pa = [pos.x, pos.y, pos.z];
        for i in 0..3 {
            assert!(
                (pa[i] - want[i]).abs() < POS_TOL_KM,
                "{label} position[{i}]: got {} want {} (Δ={:.3e} km)",
                pa[i],
                want[i],
                (pa[i] - want[i]).abs()
            );
        }
        let vel = actual
            .velocity
            .unwrap_or_else(|| panic!("{label}: velocity は Some であるべき（JPL は速度を提供）"));
        let va = [vel.x, vel.y, vel.z];
        for i in 0..3 {
            assert!(
                (va[i] - want[3 + i]).abs() < VEL_TOL_KMS,
                "{label} velocity[{i}]: got {} want {} (Δ={:.3e} km/s)",
                va[i],
                want[3 + i],
                (va[i] - want[3 + i]).abs()
            );
        }
    }

    // ==================================================================
    // A. ユニット（実ファイル不要）— from_spk_path の異常系
    // ==================================================================

    /// 観点: 存在しないパスは Err（IO 失敗。仕様上は MalformedSpk だが、
    /// DataUnavailable も両許容で受ける＝重要なのは Ok にならないこと）。
    #[test]
    fn from_spk_path_nonexistent_returns_err() {
        let path = std::path::Path::new("/no/such/de440s/__definitely_missing__.bsp");
        let err = JplEphemeris::from_spk_path(path).expect_err("存在しないパスは Err であるべき");
        assert!(
            matches!(
                err,
                EphemerisError::MalformedSpk(_) | EphemerisError::DataUnavailable
            ),
            "存在しないパスは MalformedSpk か DataUnavailable: got {err:?}"
        );
    }

    /// 観点: 不正バイトの一時ファイルは Err(MalformedSpk(_))。
    /// temp_dir に DAF マジックを持たないゴミバイトを書き、解析失敗を確認する。
    /// テスト後にファイルを削除（成否に関わらず）。
    #[test]
    fn from_spk_path_garbage_bytes_is_malformed() {
        let mut path = std::env::temp_dir();
        // 衝突回避のためプロセス固有名にする。
        path.push(format!(
            "umbra_jpl_backend_garbage_{}.bsp",
            std::process::id()
        ));
        std::fs::write(
            &path,
            b"this is not a DAF/SPK kernel: \x00\x01\x02\xff garbage bytes",
        )
        .expect("一時ファイルの書き込みに失敗");

        let result = JplEphemeris::from_spk_path(&path);
        // 後始末（assert 前に削除し、失敗時もファイルを残さない）。
        let _ = std::fs::remove_file(&path);

        let err = result.expect_err("不正バイトは Err であるべき");
        assert!(
            matches!(err, EphemerisError::MalformedSpk(_)),
            "不正バイトは MalformedSpk: got {err:?}"
        );
    }

    // ==================================================================
    // B. 実 DE440s 受入（ゲート・SPICE オラクル）
    //   data/spk/de440s.bsp 存在時のみ実行（CI 非同梱・ISSUE-036）。
    //
    // オラクル基準値の生成（出典・逐語転記）:
    //   spiceypy 8.1.2 / SPICE toolkit CSPICE_N0067
    //   （Docker python:3.12-slim, `pip install spiceypy`、本リポジトリ
    //    data/spk/de440s.bsp を /spk にマウント）
    //   import spiceypy as sp; sp.furnsh('/spk/de440s.bsp')
    //   sp.spkgeo(target, et, 'J2000', obs) -> (state[0:3]=位置km, state[3:6]=速度km/s, lt)
    //   spkgeo は target を obs 基準・J2000(=ICRF) フレームで返す。
    //   検証 ET: 0.0 / 750000000.0 / -1000000000.0（いずれも DE440s 被覆内）。
    //
    //   Sun(10) wrt SSB(0):  sp.spkgeo(10, et, 'J2000', 0)
    //     et=0.0:
    //       [-1067706.8053809535, -396036.18479594623, -138065.18428688092,
    //         0.009312571926520472, -0.01170150612817771, -0.005251266205200356]
    //     et=750000000.0:
    //       [-1249290.2370180925, -324838.0277918544, -106018.13729732925,
    //         0.007071092471312386, -0.012387616579323522, -0.005423569393921095]
    //     et=-1000000000.0:
    //       [571997.7425724803, -213216.87848410784, -98740.71783691116,
    //         0.005839325981523217, 0.007990608420087036, 0.0033003327192824682]
    //   EMB(3) wrt SSB(0):  sp.spkgeo(3, et, 'J2000', 0)
    //     et=0.0:
    //       [-27570283.695094064, 132358140.38814957, 57417728.59655908,
    //        -29.7771282160576, -5.037847169534867, -2.1843063524408897]
    //     et=750000000.0:
    //       [143698463.23331022, 33293616.59558217, 14466107.830352716,
    //        -7.778829835256609, 26.38246248432241, 11.436475412009152]
    //     et=-1000000000.0:
    //       [-123101794.51170221, -78891925.12273566, -34217509.118687585,
    //        16.493834072453577, -22.551978309772274, -9.779473646391418]
    //   Sun(10) wrt Earth(399):  sp.spkgeo(10, et, 'J2000', 399)
    //     et=0.0:
    //       [26499033.677425094, -132757417.33833946, -57556718.47053819,
    //        29.79426007042197, 5.018052308786111, 2.1753938028266693]
    //     et=750000000.0:
    //       [-144950236.73972416, -33614783.430263974, -14570046.596782928,
    //        7.775578180222956, -26.400048664969695, -11.444240960154378]
    //     et=-1000000000.0:
    //       [123678631.84380072, 78678711.61371252, 34118603.96893796,
    //        -16.487378304914515, 22.570520812042336, 9.788464963040656]
    //   Moon(301) wrt Earth(399):  sp.spkgeo(301, et, 'J2000', 399)
    //     et=0.0:
    //       [-291608.38463343546, -266716.83339423337, -76102.48709990202,
    //         0.6435313877190327, -0.6660876840916304, -0.30132570498227307]
    //     et=750000000.0:
    //       [-204374.48232724814, 302141.2790077692, 171133.4039068906,
    //        -0.8495679852656511, -0.4278447767288556, -0.19274618199514645]
    //     et=-1000000000.0:
    //       [398300.9679511811, 277.30855189943094, -13532.839849738413,
    //         0.05073349046083185, 0.8684268596657754, 0.46837121117475516]
    //   Earth(399) wrt SSB(0):  sp.spkgeo(399, et, 'J2000', 0)
    //     et=0.0:
    //       [-27566740.482806046, 132361381.15354352, 57418653.286251314,
    //        -29.78494749849545, -5.029753814914289, -2.1806450690318697]
    //     et=750000000.0:
    //       [143700946.50270608, 33289945.40247212, 14464028.4594856,
    //        -7.768507087751644, 26.38766104839037, 11.438817390760457]
    //     et=-1000000000.0:
    //       [-123106634.10122824, -78891928.49219663, -34217344.68677487,
    //        16.493217630896037, -22.56253020362225, -9.785164630321374]
    // ==================================================================

    /// 観点: state(Sun, SSB, Icrs) が spkgeo(10, et, 'J2000', 0) と一致する。
    #[test]
    fn de440s_sun_ssb_matches_spice() {
        let Some(eph) = load_de440s() else { return };
        let cases: [(f64, [f64; 6]); 3] = [
            (
                0.0,
                [
                    -1067706.8053809535,
                    -396036.18479594623,
                    -138065.18428688092,
                    0.009312571926520472,
                    -0.01170150612817771,
                    -0.005251266205200356,
                ],
            ),
            (
                750000000.0,
                [
                    -1249290.2370180925,
                    -324838.0277918544,
                    -106018.13729732925,
                    0.007071092471312386,
                    -0.012387616579323522,
                    -0.005423569393921095,
                ],
            ),
            (
                -1000000000.0,
                [
                    571997.7425724803,
                    -213216.87848410784,
                    -98740.71783691116,
                    0.005839325981523217,
                    0.007990608420087036,
                    0.0033003327192824682,
                ],
            ),
        ];
        for (et, want) in cases {
            let st = eph
                .state(
                    Body::Sun,
                    tdb_from_et(et),
                    Origin::SolarSystemBarycenter,
                    EphemerisFrame::Icrs,
                )
                .unwrap_or_else(|e| panic!("Sun/SSB et={et} で失敗: {e:?}"));
            assert_state_matches_spice(&st, want, &format!("Sun/SSB et={et}"));
        }
    }

    /// 観点: state(EarthMoonBarycenter, SSB, Icrs) が spkgeo(3, et, 'J2000', 0) と一致する。
    #[test]
    fn de440s_emb_ssb_matches_spice() {
        let Some(eph) = load_de440s() else { return };
        let cases: [(f64, [f64; 6]); 3] = [
            (
                0.0,
                [
                    -27570283.695094064,
                    132358140.38814957,
                    57417728.59655908,
                    -29.7771282160576,
                    -5.037847169534867,
                    -2.1843063524408897,
                ],
            ),
            (
                750000000.0,
                [
                    143698463.23331022,
                    33293616.59558217,
                    14466107.830352716,
                    -7.778829835256609,
                    26.38246248432241,
                    11.436475412009152,
                ],
            ),
            (
                -1000000000.0,
                [
                    -123101794.51170221,
                    -78891925.12273566,
                    -34217509.118687585,
                    16.493834072453577,
                    -22.551978309772274,
                    -9.779473646391418,
                ],
            ),
        ];
        for (et, want) in cases {
            let st = eph
                .state(
                    Body::EarthMoonBarycenter,
                    tdb_from_et(et),
                    Origin::SolarSystemBarycenter,
                    EphemerisFrame::Icrs,
                )
                .unwrap_or_else(|e| panic!("EMB/SSB et={et} で失敗: {e:?}"));
            assert_state_matches_spice(&st, want, &format!("EMB/SSB et={et}"));
        }
    }

    /// 観点: state(Earth, SSB, Icrs) が spkgeo(399, et, 'J2000', 0) と一致する。
    /// Earth_SSB = EMB_SSB + eval(399,3)。
    #[test]
    fn de440s_earth_ssb_matches_spice() {
        let Some(eph) = load_de440s() else { return };
        let cases: [(f64, [f64; 6]); 3] = [
            (
                0.0,
                [
                    -27566740.482806046,
                    132361381.15354352,
                    57418653.286251314,
                    -29.78494749849545,
                    -5.029753814914289,
                    -2.1806450690318697,
                ],
            ),
            (
                750000000.0,
                [
                    143700946.50270608,
                    33289945.40247212,
                    14464028.4594856,
                    -7.768507087751644,
                    26.38766104839037,
                    11.438817390760457,
                ],
            ),
            (
                -1000000000.0,
                [
                    -123106634.10122824,
                    -78891928.49219663,
                    -34217344.68677487,
                    16.493217630896037,
                    -22.56253020362225,
                    -9.785164630321374,
                ],
            ),
        ];
        for (et, want) in cases {
            let st = eph
                .state(
                    Body::Earth,
                    tdb_from_et(et),
                    Origin::SolarSystemBarycenter,
                    EphemerisFrame::Icrs,
                )
                .unwrap_or_else(|e| panic!("Earth/SSB et={et} で失敗: {e:?}"));
            assert_state_matches_spice(&st, want, &format!("Earth/SSB et={et}"));
        }
    }

    /// 観点: state(Sun, Geocenter, Icrs) が spkgeo(10, et, 'J2000', 399) と一致する。
    /// 地心 Sun = 10 − (3 + 399)。
    #[test]
    fn de440s_sun_geocenter_matches_spice() {
        let Some(eph) = load_de440s() else { return };
        let cases: [(f64, [f64; 6]); 3] = [
            (
                0.0,
                [
                    26499033.677425094,
                    -132757417.33833946,
                    -57556718.47053819,
                    29.79426007042197,
                    5.018052308786111,
                    2.1753938028266693,
                ],
            ),
            (
                750000000.0,
                [
                    -144950236.73972416,
                    -33614783.430263974,
                    -14570046.596782928,
                    7.775578180222956,
                    -26.400048664969695,
                    -11.444240960154378,
                ],
            ),
            (
                -1000000000.0,
                [
                    123678631.84380072,
                    78678711.61371252,
                    34118603.96893796,
                    -16.487378304914515,
                    22.570520812042336,
                    9.788464963040656,
                ],
            ),
        ];
        for (et, want) in cases {
            let st = eph
                .state(
                    Body::Sun,
                    tdb_from_et(et),
                    Origin::Geocenter,
                    EphemerisFrame::Icrs,
                )
                .unwrap_or_else(|e| panic!("Sun/Geo et={et} で失敗: {e:?}"));
            assert_state_matches_spice(&st, want, &format!("Sun/Geo et={et}"));
        }
    }

    /// 観点: state(Moon, Geocenter, Icrs) が spkgeo(301, et, 'J2000', 399) と一致する。
    /// 地心 Moon = 301 − 399。
    #[test]
    fn de440s_moon_geocenter_matches_spice() {
        let Some(eph) = load_de440s() else { return };
        let cases: [(f64, [f64; 6]); 3] = [
            (
                0.0,
                [
                    -291608.38463343546,
                    -266716.83339423337,
                    -76102.48709990202,
                    0.6435313877190327,
                    -0.6660876840916304,
                    -0.30132570498227307,
                ],
            ),
            (
                750000000.0,
                [
                    -204374.48232724814,
                    302141.2790077692,
                    171133.4039068906,
                    -0.8495679852656511,
                    -0.4278447767288556,
                    -0.19274618199514645,
                ],
            ),
            (
                -1000000000.0,
                [
                    398300.9679511811,
                    277.30855189943094,
                    -13532.839849738413,
                    0.05073349046083185,
                    0.8684268596657754,
                    0.46837121117475516,
                ],
            ),
        ];
        for (et, want) in cases {
            let st = eph
                .state(
                    Body::Moon,
                    tdb_from_et(et),
                    Origin::Geocenter,
                    EphemerisFrame::Icrs,
                )
                .unwrap_or_else(|e| panic!("Moon/Geo et={et} で失敗: {e:?}"));
            assert_state_matches_spice(&st, want, &format!("Moon/Geo et={et}"));
        }
    }

    /// 観点: state(Moon, SSB, Icrs) が spkgeo(301, et, 'J2000', 0)（Moon wrt SSB）と一致する。
    /// Moon_SSB = EMB_SSB(3) + eval(301, 3)。SSB 合成パスを Geocenter 経由でなく直接突合する
    /// （既存 de440s_moon_geocenter_matches_spice は Geocenter 基準のみで SSB 合成を未検証）。
    ///
    /// オラクル基準値の生成（spiceypy 8.1.2 / SPICE toolkit CSPICE_N0067。
    ///   Docker python:3.12-slim, `pip install spiceypy`、本リポジトリ data/spk/de440s.bsp を
    ///   /spk にマウント）:
    ///   import spiceypy as sp; sp.furnsh('/spk/de440s.bsp')
    ///   [list(sp.spkgeo(301, et, 'J2000', 0)[0]) for et in (0.0, 750000000.0, -1000000000.0)]
    ///   （spkgeo(301, et, 'J2000', 0) は Moon を SSB(0) 基準・J2000(=ICRF) で返す。
    ///    state[0:3]=位置km, state[3:6]=速度km/s。逐語転記。）
    #[test]
    fn de440s_moon_ssb_matches_spice() {
        let Some(eph) = load_de440s() else { return };
        let cases: [(f64, [f64; 6]); 3] = [
            (
                0.0,
                [
                    -27858348.86743948,
                    132094664.32014929,
                    57342550.79915141,
                    -29.141416110776415,
                    -5.69584149900592,
                    -2.481970774014143,
                ],
            ),
            (
                750000000.0,
                [
                    143496572.02037886,
                    33592086.68147989,
                    14635161.86339249,
                    -8.618075073017295,
                    25.959816271661513,
                    11.246071208765311,
                ],
            ),
            (
                -1000000000.0,
                [
                    -122708333.13327706,
                    -78891651.18364473,
                    -34230877.52662461,
                    16.54395112135687,
                    -21.69410334395647,
                    -9.316793419146618,
                ],
            ),
        ];
        for (et, want) in cases {
            let st = eph
                .state(
                    Body::Moon,
                    tdb_from_et(et),
                    Origin::SolarSystemBarycenter,
                    EphemerisFrame::Icrs,
                )
                .unwrap_or_else(|e| panic!("Moon/SSB et={et} で失敗: {e:?}"));
            assert_state_matches_spice(&st, want, &format!("Moon/SSB et={et}"));
        }
    }

    /// 観点: Body::Earth + Origin::Geocenter は原点（位置 ZERO・速度 Some(ZERO)）。
    /// 厳密 0（許容ゼロ）。
    #[test]
    fn de440s_earth_geocenter_is_exact_zero() {
        let Some(eph) = load_de440s() else { return };
        let st = eph
            .state(
                Body::Earth,
                tdb_from_et(0.0),
                Origin::Geocenter,
                EphemerisFrame::Icrs,
            )
            .expect("Earth/Geocenter は成功すべき");
        assert_eq!(
            st.position,
            Vector3::ZERO,
            "Earth/Geocenter の位置は厳密 ZERO であるべき"
        );
        assert_eq!(
            st.velocity,
            Some(Vector3::ZERO),
            "Earth/Geocenter の速度は Some(ZERO) であるべき"
        );
    }

    /// 観点: velocity が Some であること（少なくとも 1 ケース・確定仕様）。
    #[test]
    fn de440s_velocity_is_some() {
        let Some(eph) = load_de440s() else { return };
        let st = eph
            .state(
                Body::Sun,
                tdb_from_et(0.0),
                Origin::SolarSystemBarycenter,
                EphemerisFrame::Icrs,
            )
            .expect("Sun/SSB は成功すべき");
        assert!(
            st.velocity.is_some(),
            "JPL は Chebyshev 微分で速度を提供するため velocity は Some"
        );
    }

    /// 観点: EphemerisFrame::EclipticOfDate は本バックエンド未提供 → DataUnavailable。
    /// 代表として Sun/Geocenter/EclipticOfDate。
    #[test]
    fn de440s_ecliptic_of_date_is_data_unavailable() {
        let Some(eph) = load_de440s() else { return };
        let err = eph
            .state(
                Body::Sun,
                tdb_from_et(0.0),
                Origin::Geocenter,
                EphemerisFrame::EclipticOfDate,
            )
            .expect_err("EclipticOfDate は Err であるべき");
        assert_eq!(
            err,
            EphemerisError::DataUnavailable,
            "EclipticOfDate は DataUnavailable: got {err:?}"
        );
    }

    /// 観点: 被覆外時刻は OutOfSupportedRange。
    /// supported_range の十分手前/後ろ（DE440s ~1849-12-26〜2150-01-21）を試す。
    /// et の基準は J2000=2000-01-01。1700 年相当 ≈ −9.46e9 s、2400 年相当 ≈ +1.26e10 s。
    #[test]
    fn de440s_out_of_range_is_out_of_supported_range() {
        let Some(eph) = load_de440s() else { return };
        // 被覆十分手前（西暦 ~1700）と十分後ろ（西暦 ~2400）。
        let far_before = -9.46e9_f64; // ≈ 300 年手前
        let far_after = 1.26e10_f64; // ≈ 400 年後ろ
        for et in [far_before, far_after] {
            let err = eph
                .state(
                    Body::Moon,
                    tdb_from_et(et),
                    Origin::Geocenter,
                    EphemerisFrame::Icrs,
                )
                .err()
                .unwrap_or_else(|| panic!("被覆外 et={et} は Err であるべき"));
            assert_eq!(
                err,
                EphemerisError::OutOfSupportedRange,
                "被覆外 et={et} は OutOfSupportedRange: got {err:?}"
            );
        }
    }

    /// 観点: supported_range() は start < end かつ 1900–2100 を内包する（おおよそ）。
    /// 1900-01-01 ≈ JD 2415020.5、2100-01-01 ≈ JD 2488069.5。
    #[test]
    fn de440s_supported_range_contains_1900_2100() {
        let Some(eph) = load_de440s() else { return };
        let range = eph.supported_range();
        let start_jd = range.start.jd2().jd();
        let end_jd = range.end.jd2().jd();
        assert!(
            start_jd < end_jd,
            "supported_range は start<end: start_jd={start_jd} end_jd={end_jd}"
        );
        // 1900-01-01 / 2100-01-01 の JD（おおよそ。被覆を内包すれば十分）。
        let jd_1900 = 2_415_020.5_f64;
        let jd_2100 = 2_488_069.5_f64;
        assert!(
            start_jd <= jd_1900,
            "supported_range の start は 1900 以前: start_jd={start_jd}"
        );
        assert!(
            end_jd >= jd_2100,
            "supported_range の end は 2100 以後: end_jd={end_jd}"
        );
    }

    /// 観点: supported_range() の被覆区間を二側（両端）の窓で厳密にピン留めする。
    ///
    /// 片側包含（start≤1900・end≥2100）だけでは、`supported_range()` 内部の
    /// ET 秒→`TdbInstant`（JD）変換が壊れても（例: 除算 et/86400 を乗算 et*86400 や
    /// 剰余 et%86400 に変えても）極端に外れた値が片側境界を満たして通過しうる。
    /// そこで実値を両側の窓で縛り、変換破損を検出する。
    ///
    /// 窓の根拠: DE440s の被覆は ET 1849-12-26〜2150-01-21（data/spk/PROVENANCE.md）。
    ///   1849-12-26 ≈ JD 2396758.5、2150-01-21 ≈ JD 2506332.5。
    ///   start_jd は 1849〜1850 年相当、end_jd は 2149〜2151 年相当の窓に収める。
    ///   窓幅は「正しい値は通すが、別の年へずれた／極端に外れた値は弾く」程度に狭く取る。
    #[test]
    fn de440s_supported_range_matches_de440s_span() {
        let Some(eph) = load_de440s() else { return };
        let range = eph.supported_range();
        let start_jd = range.start.jd2().jd();
        let end_jd = range.end.jd2().jd();
        // start_jd は 1849-12-26（≈ JD 2396758.5）近傍。1849〜1850 年相当の窓で縛る。
        assert!(
            2_396_000.0 < start_jd && start_jd < 2_397_500.0,
            "supported_range の start は DE440s 被覆開始（ET 1849-12-26 ≈ JD 2396758.5）\
             近傍であるべき: got start_jd={start_jd}"
        );
        // end_jd は 2150-01-21（≈ JD 2506332.5）近傍。2149〜2151 年相当の窓で縛る。
        assert!(
            2_505_500.0 < end_jd && end_jd < 2_507_000.0,
            "supported_range の end は DE440s 被覆終了（ET 2150-01-21 ≈ JD 2506332.5）\
             近傍であるべき: got end_jd={end_jd}"
        );
    }

    /// 観点: metadata() の緩い検証。model に "DE440" を含む・license に "非同梱" を含む・
    /// source が空でない・supported が supported_range() と一致する。
    /// max_residual_arcsec は実装の申告に依存するため検証しない（オラクルのため 0/NaN）。
    #[test]
    fn de440s_metadata_has_expected_fields() {
        let Some(eph) = load_de440s() else { return };
        let meta = eph.metadata();
        assert!(
            meta.model.contains("DE440"),
            "metadata.model は DE 版を含む（例 'JPL DE440s'）: got {:?}",
            meta.model
        );
        assert!(
            meta.license.contains("非同梱"),
            "metadata.license は「非同梱」を含む: got {:?}",
            meta.license
        );
        assert!(
            !meta.source.is_empty(),
            "metadata.source は空でない（NAIF / Park et al. 2021 など）"
        );
        assert_eq!(
            meta.supported,
            eph.supported_range(),
            "metadata.supported は supported_range() と一致すべき"
        );
    }

    /// 観点（任意・内部整合）: state(Sun, Geocenter) ≈ state(Sun, SSB) − state(Earth, SSB)。
    /// 原点差の定義（地心 = SSB 基準 − Earth_SSB）が成分で整合すること。
    /// 大きな量（~1.5e8 km）同士の差なので相対許容を緩める（位置 ~1 km、速度 ~1e-4 km/s）。
    #[test]
    fn de440s_geocenter_equals_ssb_difference() {
        let Some(eph) = load_de440s() else { return };
        let et = 750000000.0;
        let sun_geo = eph
            .state(
                Body::Sun,
                tdb_from_et(et),
                Origin::Geocenter,
                EphemerisFrame::Icrs,
            )
            .expect("Sun/Geo");
        let sun_ssb = eph
            .state(
                Body::Sun,
                tdb_from_et(et),
                Origin::SolarSystemBarycenter,
                EphemerisFrame::Icrs,
            )
            .expect("Sun/SSB");
        let earth_ssb = eph
            .state(
                Body::Earth,
                tdb_from_et(et),
                Origin::SolarSystemBarycenter,
                EphemerisFrame::Icrs,
            )
            .expect("Earth/SSB");

        let want_pos = sun_ssb.position - earth_ssb.position;
        let gp = [sun_geo.position.x, sun_geo.position.y, sun_geo.position.z];
        let wp = [want_pos.x, want_pos.y, want_pos.z];
        for i in 0..3 {
            assert!(
                (gp[i] - wp[i]).abs() < 1.0,
                "geocenter 位置整合[{i}]: got {} want {} (Δ={:.3e} km)",
                gp[i],
                wp[i],
                (gp[i] - wp[i]).abs()
            );
        }
        // 速度も同様に SSB 差と整合（両者 Some を前提）。
        let want_vel =
            sun_ssb.velocity.expect("Sun/SSB vel") - earth_ssb.velocity.expect("Earth/SSB vel");
        let gv_opt = sun_geo.velocity.expect("Sun/Geo vel");
        let gv = [gv_opt.x, gv_opt.y, gv_opt.z];
        let wv = [want_vel.x, want_vel.y, want_vel.z];
        for i in 0..3 {
            assert!(
                (gv[i] - wv[i]).abs() < 1.0e-4,
                "geocenter 速度整合[{i}]: got {} want {} (Δ={:.3e} km/s)",
                gv[i],
                wv[i],
                (gv[i] - wv[i]).abs()
            );
        }
    }
}
