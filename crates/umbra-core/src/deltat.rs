//! ΔT (= TT − UT1) モデルと UT1 変換（`docs/issues/ISSUE-007`、`docs/algorithms/01-time-scales.md`）。
//!
//! [`EspenakMeeusDeltaT`] は Espenak & Meeus（NASA TP-2006-214141, `eclipse.gsfc.nasa.gov`）の
//! **区分多項式**。公開された数式（バルクデータではない）なのでライセンス問題なく組み込める。
//! 厳密な近代 ΔT は EOP（UT1−UTC 実測, 将来の ISSUE-007 完全版）から得るが、本モデルは
//! その予測値・外挿値と不確実性帯（accuracy.md §0）を与える。

use crate::calendar::jd2_to_gregorian;
use crate::constants::TT_MINUS_TAI_SECONDS;
use crate::eop::EarthOrientation;
use crate::error::TimeError;
use crate::julian::JulianDate2;
use crate::time::{tai_minus_utc, TtInstant, Ut1Instant, UtcInstant};

const SECONDS_PER_DAY: f64 = 86_400.0;

/// ΔT = TT − UT1 モデル。
pub trait DeltaTModel {
    /// 十進年に対する ΔT（秒）。
    fn delta_t_seconds(&self, decimal_year: f64) -> f64;
    /// ΔT の不確実性（秒, 1σ 目安。accuracy.md §0 の不確実性帯）。
    fn uncertainty_seconds(&self, decimal_year: f64) -> f64;
    /// モデル名（`CalculationMetadata.delta_t_model` のレシピ識別子, accuracy.md §0）。
    fn model_name(&self) -> &'static str;
}

/// Espenak & Meeus の ΔT 区分多項式（NASA TP-2006-214141）。
#[derive(Debug, Clone, Copy, Default)]
pub struct EspenakMeeusDeltaT;

/// 1820 を基準とした長期放物線 ΔT = −20 + 32·u²（u = (y−1820)/100）。範囲外の既定。
fn long_term(y: f64) -> f64 {
    let u = (y - 1820.0) / 100.0;
    -20.0 + 32.0 * u * u
}

impl DeltaTModel for EspenakMeeusDeltaT {
    fn delta_t_seconds(&self, y: f64) -> f64 {
        if y < 1900.0 {
            long_term(y)
        } else if y < 1920.0 {
            let t = y - 1900.0;
            -2.79 + 1.494119 * t - 0.0598939 * t * t + 0.0061966 * t * t * t
                - 0.000197 * t * t * t * t
        } else if y < 1941.0 {
            let t = y - 1920.0;
            21.20 + 0.84493 * t - 0.076100 * t * t + 0.0020936 * t * t * t
        } else if y < 1961.0 {
            let t = y - 1950.0;
            29.07 + 0.407 * t - t * t / 233.0 + t * t * t / 2547.0
        } else if y < 1986.0 {
            let t = y - 1975.0;
            45.45 + 1.067 * t - t * t / 260.0 - t * t * t / 718.0
        } else if y < 2005.0 {
            let t = y - 2000.0;
            63.86 + 0.3345 * t - 0.060374 * t * t
                + 0.0017275 * t * t * t
                + 0.000651814 * t * t * t * t
                + 0.00002373599 * t * t * t * t * t
        } else if y < 2050.0 {
            let t = y - 2000.0;
            62.92 + 0.32217 * t + 0.005589 * t * t
        } else if y < 2150.0 {
            -20.0 + 32.0 * ((y - 1820.0) / 100.0) * ((y - 1820.0) / 100.0) - 0.5628 * (2150.0 - y)
        } else {
            long_term(y)
        }
    }

    fn uncertainty_seconds(&self, y: f64) -> f64 {
        // 粗い目安（要確認）。実測 EOP 域では UT1−UTC から精密に得るべきで、本値は
        // モデル予測の帯。順序しきい値（不連続）にして各境界・式を検証可能にする。
        if y < 1900.0 {
            5.0 + 0.1 * (1900.0 - y) // 古い年代ほど増大
        } else if y < 1955.0 {
            2.0
        } else if y < 2006.0 {
            0.5 // 近代（おおむね実測 EOP 域）
        } else {
            1.0 + 0.5 * (y - 2006.0) // 2006 以降の外挿は年々増大
        }
    }

    fn model_name(&self) -> &'static str {
        "Espenak-Meeus"
    }
}

/// 十進年 `year + (month − 0.5)/12`（Espenak 慣習）を JD から求める。
pub fn decimal_year(jd: JulianDate2) -> f64 {
    let (year, month, ..) = jd2_to_gregorian(jd);
    f64::from(year) + (f64::from(month) - 0.5) / 12.0
}

/// TT → UT1（`UT1 = TT − ΔT`）。
pub fn tt_to_ut1<M: DeltaTModel>(tt: TtInstant, model: &M) -> Ut1Instant {
    let dt = model.delta_t_seconds(decimal_year(tt.jd2()));
    Ut1Instant::from_jd2(tt.jd2().add_days(-dt / SECONDS_PER_DAY))
}

/// UT1 → TT（`TT = UT1 + ΔT`）。ΔT は UT1 の年から評価する（TT 年との差は ΔT に無視可能）。
pub fn ut1_to_tt<M: DeltaTModel>(ut1: Ut1Instant, model: &M) -> TtInstant {
    let dt = model.delta_t_seconds(decimal_year(ut1.jd2()));
    TtInstant::from_jd2(ut1.jd2().add_days(dt / SECONDS_PER_DAY))
}

/// UTC → UT1（`UT1 = UTC + (UT1−UTC)`）。`EarthOrientation` から UT1−UTC（秒）を引く。
///
/// coverage 外は [`TimeError::MissingEarthOrientationData`]（EOP の Missing を透過）。
/// `TimeScales`（ISSUE-042）が未整備のため当面は自由関数として供給する（ISSUE-007 §公開IF）。
pub fn utc_to_ut1<E: EarthOrientation>(utc: UtcInstant, eo: &E) -> Result<Ut1Instant, TimeError> {
    let ut1_minus_utc = eo.ut1_minus_utc(utc)?;
    Ok(Ut1Instant::from_jd2(
        utc.jd2().add_days(ut1_minus_utc / SECONDS_PER_DAY),
    ))
}

/// EOP 実測域での ΔT 不確実性（秒, 1σ 目安。accuracy.md §0 の「過去/近傍 <0.1 s」）。
///
/// IERS EOP C04 の UT1−UTC は数 ms 精度（accuracy.md §2.3）、閏秒 (TAI−UTC) は厳密値であり、
/// 日次線形補間の誤差も sub-ms に収まるため、恒等式由来 ΔT の不確実性は数 ms オーダ。
/// 保守的に 5 ms を採る（<0.1 s を満たす）。
pub const EOP_DELTA_T_UNCERTAINTY_SECONDS: f64 = 0.005;

/// 暦上の瞬時（UTC）に対する ΔT = TT − UT1 とその不確実性帯の供給（accuracy.md §0）。
///
/// 低水準の区分多項式 [`DeltaTModel`]（十進年入力）と異なり、本 trait は EOP 実測の恒等式と
/// 長期外挿を合成する高水準モデルで、瞬時 UTC を入力に取る。将来 `CalculationMetadata` へ
/// ΔT モデルと不確実性帯を供給する源（ISSUE-007 §目的）。
pub trait DeltaTSource: Send + Sync {
    /// `utc` における ΔT = TT − UT1（秒）。
    fn delta_t_seconds(&self, utc: UtcInstant) -> f64;
    /// `utc` における ΔT の不確実性（秒, 1σ 目安。accuracy.md §0 の不確実性帯）。
    fn uncertainty_seconds(&self, utc: UtcInstant) -> f64;
}

/// EOP 実測由来の高精度 ΔT（恒等式）と Espenak–Meeus 外挿を合成する [`DeltaTSource`]。
///
/// - **EOP coverage 内 かつ 閏秒テーブル域（≥1972）**: 恒等式
///   `ΔT = (TAI−UTC) + 32.184 − (UT1−UTC)`（conventions §6 / ISSUE-007 §数式）で高精度に算出。
///   不確実性は [`EOP_DELTA_T_UNCERTAINTY_SECONDS`]（IERS 実測帯 <0.1 s）。
/// - **それ以外**（EOP 範囲外、または 1972 以前で閏秒テーブル未定義）: [`EspenakMeeusDeltaT`] の
///   区分多項式へ外挿し、不確実性も [`EspenakMeeusDeltaT::uncertainty_seconds`] に従う（将来ほど増大）。
///
/// 1962–1972 は EOP に UT1−UTC があっても閏秒が未定義のため恒等式を使わず外挿する
/// （ISSUE-007 §実装メモ・1972 以前の橋渡し）。
#[derive(Debug)]
pub struct CompositeDeltaT<E: EarthOrientation> {
    eop: E,
    fallback: EspenakMeeusDeltaT,
}

impl<E: EarthOrientation> CompositeDeltaT<E> {
    /// EOP ソースから合成器を構築する（外挿器は既定の [`EspenakMeeusDeltaT`]）。
    pub fn new(eop: E) -> Self {
        Self {
            eop,
            fallback: EspenakMeeusDeltaT,
        }
    }

    /// 恒等式が適用可能（EOP coverage 内 ∧ 閏秒域）なら `Some(ΔT)`、不可なら `None`（外挿へ委ねる）。
    /// 恒等式 `ΔT = (TAI−UTC) + 32.184 − (UT1−UTC)`。両ソースのいずれかが欠けると `None`。
    fn identity_delta_t(&self, utc: UtcInstant) -> Option<f64> {
        let tai_utc = tai_minus_utc(utc).ok()?;
        let ut1_utc = self.eop.ut1_minus_utc(utc).ok()?;
        Some(tai_utc + TT_MINUS_TAI_SECONDS - ut1_utc)
    }
}

impl<E: EarthOrientation> DeltaTSource for CompositeDeltaT<E> {
    fn delta_t_seconds(&self, utc: UtcInstant) -> f64 {
        self.identity_delta_t(utc)
            .unwrap_or_else(|| self.fallback.delta_t_seconds(decimal_year(utc.jd2())))
    }

    fn uncertainty_seconds(&self, utc: UtcInstant) -> f64 {
        if self.identity_delta_t(utc).is_some() {
            EOP_DELTA_T_UNCERTAINTY_SECONDS
        } else {
            self.fallback.uncertainty_seconds(decimal_year(utc.jd2()))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants::TT_MINUS_TAI_SECONDS;
    use crate::eop::{EopRecord, IersEopData};
    use crate::error::TimeError;
    use crate::metadata::DataSetMetadata;
    use crate::time::{tai_minus_utc, UtcInstant};

    const DT: EspenakMeeusDeltaT = EspenakMeeusDeltaT;

    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    // ---- EOP 合成スキャフォールド（eop.rs テストのパターンを流用）-----------
    // bundled データは別 crate にあり core からは使えないため、coverage を決定的に
    // 制御できる合成 IersEopData をテスト内で構築する。

    /// IERS EOP 14 C04 実測オラクル（独立計算・eop.rs と同じ verbatim 値）。
    /// 2020-01-01 / MJD 58849。
    const MJD_20200101: i32 = 58849;
    const UT1_20200101: f64 = -0.1771222;
    const XP_20200101: f64 = 0.076609;
    const YP_20200101: f64 = 0.282358;
    /// 2020-01-02 / MJD 58850。
    const MJD_20200102: i32 = 58850;
    const UT1_20200102: f64 = -0.1775806;
    const XP_20200102: f64 = 0.074635;
    const YP_20200102: f64 = 0.282666;
    /// 1962-01-01 / MJD 37665（閏秒域外 = 恒等式不可・外挿になる境界）。
    const MJD_1962: i32 = 37665;
    const UT1_1962: f64 = 0.0326338;
    const XP_1962: f64 = -0.012700;
    const YP_1962: f64 = 0.213000;

    /// provenance 完全な代表 metadata（全フィールド非空）。
    fn eop_metadata() -> DataSetMetadata {
        DataSetMetadata {
            name: "iers-eop-c04".to_string(),
            version: "EOP 14 C04".to_string(),
            source: "IERS Earth Orientation Center, datacenter.iers.org".to_string(),
            license: "public-domain".to_string(),
            valid_from: "1962-01-01".to_string(),
            valid_to: "2020-01-02".to_string(),
            checksum: "0000000000000000000000000000000000000000000000000000000000000000"
                .to_string(),
        }
    }

    /// 2 点 {58849, 58850}（2020-01-01/02）のみの EOP。coverage = 2020 のみ。
    fn eop_2020() -> IersEopData {
        IersEopData::from_records(
            vec![
                EopRecord::new(MJD_20200101, UT1_20200101, XP_20200101, YP_20200101),
                EopRecord::new(MJD_20200102, UT1_20200102, XP_20200102, YP_20200102),
            ],
            "EOP 14 C04".to_string(),
            eop_metadata(),
        )
        .expect("two adjacent ascending 2020 records build")
    }

    /// 3 点 {37665, 58849, 58850}。coverage に 1962 を含む（閏秒域外の外挿分岐検証用）。
    fn eop_1962_to_2020() -> IersEopData {
        IersEopData::from_records(
            vec![
                EopRecord::new(MJD_1962, UT1_1962, XP_1962, YP_1962),
                EopRecord::new(MJD_20200101, UT1_20200101, XP_20200101, YP_20200101),
                EopRecord::new(MJD_20200102, UT1_20200102, XP_20200102, YP_20200102),
            ],
            "EOP 14 C04".to_string(),
            eop_metadata(),
        )
        .expect("three ascending records (1962..2020) build")
    }

    fn eop_utc(y: i32, mo: u8, d: u8, h: u8, mi: u8, s: f64) -> UtcInstant {
        UtcInstant::from_gregorian(y, mo, d, h, mi, s).expect("valid calendar date")
    }

    #[test]
    fn piece_values_match_formulas_at_interior_points() {
        // 各区分の内部点（非ゼロ t）で手計算オラクルと一致 → 区分内の演算子取り違えを検出。
        assert!(
            close(DT.delta_t_seconds(1910.0), 10.3884, 1e-3),
            "{}",
            DT.delta_t_seconds(1910.0)
        );
        assert!(
            close(DT.delta_t_seconds(1930.0), 24.1329, 1e-3),
            "{}",
            DT.delta_t_seconds(1930.0)
        );
        assert!(
            close(DT.delta_t_seconds(1955.0), 31.0468, 1e-3),
            "{}",
            DT.delta_t_seconds(1955.0)
        );
        assert!(
            close(DT.delta_t_seconds(1980.0), 50.5148, 1e-3),
            "{}",
            DT.delta_t_seconds(1980.0)
        );
        assert!(
            close(DT.delta_t_seconds(1995.0), 60.7954, 1e-3),
            "{}",
            DT.delta_t_seconds(1995.0)
        );
        assert!(
            close(DT.delta_t_seconds(2030.0), 77.6151, 1e-3),
            "{}",
            DT.delta_t_seconds(2030.0)
        );
        assert!(
            close(DT.delta_t_seconds(2100.0), 202.74, 1e-2),
            "{}",
            DT.delta_t_seconds(2100.0)
        );
        assert!(
            close(DT.delta_t_seconds(2200.0), 442.08, 1e-2),
            "{}",
            DT.delta_t_seconds(2200.0)
        );
        assert!(close(DT.delta_t_seconds(1850.0), long_term(1850.0), 1e-9));
    }

    #[test]
    fn piece_starts_match_constants_at_boundaries() {
        // 各区分の開始年（境界）で、その区分の定数項に一致 → 境界比較の向き（< vs <=）を検出。
        assert!(close(DT.delta_t_seconds(1900.0), -2.79, 1e-9));
        assert!(close(DT.delta_t_seconds(1920.0), 21.20, 1e-9));
        assert!(close(DT.delta_t_seconds(1950.0), 29.07, 1e-9)); // 1941–1961 区分 t=0
        assert!(close(DT.delta_t_seconds(1975.0), 45.45, 1e-9)); // 1961–1986 区分 t=0
        assert!(close(DT.delta_t_seconds(2000.0), 63.86, 1e-9)); // 1986–2005 区分 t=0
    }

    #[test]
    fn delta_t_2000_matches_known_value() {
        // ΔT(2000.0) ≈ 63.8 s（観測既知）に近い。
        assert!(close(DT.delta_t_seconds(2000.0), 63.8, 0.1));
    }

    #[test]
    fn uncertainty_branches_boundaries_and_growth() {
        // 各分岐の内部値（式の演算子取り違えを検出）。
        assert!(close(DT.uncertainty_seconds(1980.0), 0.5, 1e-12));
        assert!(close(DT.uncertainty_seconds(1930.0), 2.0, 1e-12));
        assert!(close(
            DT.uncertainty_seconds(1800.0),
            5.0 + 0.1 * 100.0,
            1e-9
        )); // 15.0
        assert!(close(
            DT.uncertainty_seconds(2106.0),
            1.0 + 0.5 * 100.0,
            1e-9
        )); // 51.0
            // 境界（不連続なので比較の向き < / <= を検出）。
        assert!(close(DT.uncertainty_seconds(1900.0), 2.0, 1e-12)); // y<1900 false
        assert!(close(DT.uncertainty_seconds(1955.0), 0.5, 1e-12)); // y<1955 false
        assert!(close(DT.uncertainty_seconds(2006.0), 1.0, 1e-12)); // y<2006 false
                                                                    // 将来へ増大。
        assert!(DT.uncertainty_seconds(2100.0) > DT.uncertainty_seconds(2030.0));
    }

    #[test]
    fn decimal_year_uses_mid_month_convention() {
        let jd = crate::calendar::gregorian_to_jd2(2000, 7, 2, 0, 0, 0.0).unwrap();
        // 7月 → 2000 + (7-0.5)/12 = 2000.5417
        assert!(close(decimal_year(jd), 2000.0 + 6.5 / 12.0, 1e-9));
    }

    #[test]
    fn tt_to_ut1_subtracts_delta_t() {
        let tt =
            TtInstant::from_jd2(crate::calendar::gregorian_to_jd2(2010, 1, 1, 0, 0, 0.0).unwrap());
        let ut1 = tt_to_ut1(tt, &DT);
        let dt = DT.delta_t_seconds(decimal_year(tt.jd2()));
        let diff_s = ut1.jd2().days_since(tt.jd2()) * SECONDS_PER_DAY;
        assert!(close(diff_s, -dt, 1e-6), "diff = {diff_s}, dt = {dt}");
        assert!(dt > 60.0 && dt < 75.0, "dt(2010) = {dt}");
    }

    #[test]
    fn ut1_tt_round_trip() {
        let tt =
            TtInstant::from_jd2(crate::calendar::gregorian_to_jd2(2035, 9, 2, 1, 30, 0.0).unwrap());
        let back = ut1_to_tt(tt_to_ut1(tt, &DT), &DT);
        assert!(back.jd2().days_since(tt.jd2()).abs() * SECONDS_PER_DAY < 1e-3);
    }

    // ===================================================================
    // ISSUE-007 EOP part P3: utc_to_ut1 / DeltaTSource / CompositeDeltaT
    // ===================================================================

    // ---- utc_to_ut1 自由関数 ------------------------------------------

    /// 正常系: 戻り値 UT1 の JD は UTC の JD に (UT1−UTC)/86400 日を加えたもの。
    /// 2020-01-01 0h は厳密日ルックアップで UT1−UTC = -0.1771222 s。
    /// 加算（符号・/86400 の除数・対象は UTC.jd2）を exact に固定する。
    /// 変異: `+`→`-`（符号反転）、`/86400`→`*86400` や除数取り違え、UT1−UTC を引かない、を殺す。
    #[test]
    fn utc_to_ut1_adds_ut1_minus_utc_to_jd() {
        let eo = eop_2020();
        let utc = eop_utc(2020, 1, 1, 0, 0, 0.0);
        let ut1 = utc_to_ut1(utc, &eo).expect("2020-01-01 is within EOP coverage");
        // 期待 JD = UTC.jd + (UT1−UTC)/86400。UT1−UTC は負なので JD は僅かに減る。
        let want_jd = utc.jd2().jd() + UT1_20200101 / SECONDS_PER_DAY;
        assert!(
            (ut1.jd2().jd() - want_jd).abs() < 1e-12,
            "ut1 jd = {}, want {want_jd}",
            ut1.jd2().jd()
        );
        // 差分を秒で取り直し、UT1−UTC（実測オラクル）に厳密一致することも固定する。
        let diff_s = ut1.jd2().days_since(utc.jd2()) * SECONDS_PER_DAY;
        assert!(
            (diff_s - UT1_20200101).abs() < 1e-6,
            "diff = {diff_s} s, want {UT1_20200101} s"
        );
    }

    /// coverage 外（最後のレコードより後 2020-01-03）は Err(MissingEarthOrientationData)。
    /// EOP の Missing を素通しせず正しい variant に写すことを固定する
    /// （別 variant への取り違え・`?` 抜けで Ok を返す変異を殺す）。
    #[test]
    fn utc_to_ut1_outside_coverage_is_missing() {
        let eo = eop_2020();
        let outside = eop_utc(2020, 1, 3, 0, 0, 0.0); // MJD 58851 > last(58850)
        assert_eq!(
            utc_to_ut1(outside, &eo).unwrap_err(),
            TimeError::MissingEarthOrientationData
        );
    }

    // ---- CompositeDeltaT::delta_t_seconds 恒等式分岐 -------------------

    /// 恒等式分岐（EOP 域 ∧ ≥1972, 例 2020-01-01）:
    /// ΔT = (TAI−UTC) + 32.184 − (UT1−UTC) = 37.0 + 32.184 − (−0.1771222) = 69.3611222 s。
    /// 独立に計算した既知 ΔT(2020)≈69.36 s と一致する厳密値。
    /// 変異: 32.184 の取り違え、UT1−UTC の符号（`−`→`+`）、閏秒項の脱落、恒等式と外挿の取り違え
    /// （外挿は ~71.62 s で 2.26 s 異なるため判別可能）を殺す。
    #[test]
    fn composite_delta_t_uses_identity_in_eop_and_leap_domain() {
        let comp = CompositeDeltaT::new(eop_2020());
        let utc = eop_utc(2020, 1, 1, 0, 0, 0.0);
        let got = comp.delta_t_seconds(utc);
        // 恒等式オラクル（各項を verbatim 値から独立に組み立てる）。
        let want = 37.0 + 32.184 - UT1_20200101; // = 69.3611222
        assert!(
            (got - want).abs() < 1e-9,
            "identity ΔT(2020-01-01) = {got}, want {want}"
        );
        // 既知 ΔT(2020) ≈ 69.36 s（独立既知値）と一致。
        assert!(
            (got - 69.3611222).abs() < 1e-7,
            "ΔT(2020-01-01) = {got}, want ≈ 69.3611222"
        );
    }

    /// 恒等式の各構成項を、ライブラリの一次ソース（tai_minus_utc / ut1_minus_utc /
    /// TT_MINUS_TAI_SECONDS 定数）から組み立て直し、CompositeDeltaT がその恒等式どおりに
    /// 算出していることを固定する（恒等式の係数・符号・項構成の取り違えを殺す）。
    #[test]
    fn composite_delta_t_identity_matches_term_by_term() {
        use crate::eop::EarthOrientation;
        let eo = eop_2020();
        let comp = CompositeDeltaT::new(eop_2020());
        let utc = eop_utc(2020, 1, 1, 0, 0, 0.0);
        let tai_utc = tai_minus_utc(utc).expect("2020 is in leap-second table");
        let ut1_utc = eo.ut1_minus_utc(utc).expect("2020 is in EOP coverage");
        let want = tai_utc + TT_MINUS_TAI_SECONDS - ut1_utc;
        assert!(
            (comp.delta_t_seconds(utc) - want).abs() < 1e-12,
            "ΔT = {}, term-by-term want {want}",
            comp.delta_t_seconds(utc)
        );
    }

    // ---- CompositeDeltaT::delta_t_seconds 外挿分岐 --------------------

    /// 外挿分岐（EOP coverage 外, 例 2025-06）→ Espenak–Meeus と一致。
    /// EOP のみを持つ coverage（2020 のみ）の外を引くと、恒等式ではなく
    /// EspenakMeeusDeltaT.delta_t_seconds(decimal_year(jd)) に厳密一致する。
    /// 変異: 外挿器の引数（decimal_year を渡さない）、外挿分岐の脱落、を殺す。
    #[test]
    fn composite_delta_t_extrapolates_outside_eop_coverage() {
        let comp = CompositeDeltaT::new(eop_2020());
        let utc = eop_utc(2025, 6, 15, 0, 0, 0.0); // 2020 coverage の外
        let got = comp.delta_t_seconds(utc);
        let want = EspenakMeeusDeltaT.delta_t_seconds(decimal_year(utc.jd2()));
        assert!(
            (got - want).abs() < 1e-12,
            "extrapolated ΔT = {got}, want EM {want}"
        );
    }

    /// 外挿分岐（EOP 域内だが 1972 以前 = 閏秒テーブル未定義, 例 1962-01-01）→
    /// 恒等式ではなく Espenak–Meeus と一致。
    /// この境界が「恒等式には閏秒域条件（tai_minus_utc が Ok）も必要」であることを固定する
    /// 重要テスト: EOP coverage 内であることだけで恒等式を選ぶ変異（閏秒域条件の脱落）を殺す。
    /// 1962 で恒等式を誤って使うと閏秒が未定義（Err）なので、外挿になっていることを exact 比較で示す。
    #[test]
    fn composite_delta_t_extrapolates_when_before_leap_seconds() {
        let comp = CompositeDeltaT::new(eop_1962_to_2020());
        let utc = eop_utc(1962, 1, 1, 0, 0, 0.0); // EOP 域内だが 1972 以前。
                                                  // 前提: 1962 は閏秒テーブル域外（恒等式は使えない）。
        assert!(
            tai_minus_utc(utc).is_err(),
            "前提: 1962 は閏秒テーブル未定義であること"
        );
        let got = comp.delta_t_seconds(utc);
        let want = EspenakMeeusDeltaT.delta_t_seconds(decimal_year(utc.jd2()));
        assert!(
            (got - want).abs() < 1e-12,
            "1962 ΔT = {got}, want EM extrapolation {want}"
        );
    }

    // ---- CompositeDeltaT::uncertainty_seconds -------------------------

    /// 恒等式域（EOP 実測 ∧ ≥1972, 例 2020-01-01）の不確実性は固定定数
    /// EOP_DELTA_T_UNCERTAINTY_SECONDS (=0.005) に exact 一致し、かつ < 0.1（IERS 実測帯）。
    /// 変異: 定数値の取り違え、外挿器の不確実性を返す（実測域で 0.5 等になる）、を殺す。
    #[test]
    fn composite_uncertainty_in_eop_domain_is_fixed_constant() {
        let comp = CompositeDeltaT::new(eop_2020());
        let utc = eop_utc(2020, 1, 1, 0, 0, 0.0);
        let got = comp.uncertainty_seconds(utc);
        assert!(
            (got - EOP_DELTA_T_UNCERTAINTY_SECONDS).abs() < 1e-12,
            "EOP-domain uncertainty = {got}, want {EOP_DELTA_T_UNCERTAINTY_SECONDS}"
        );
        // 公開定数そのものが 0.005 であることも固定する。
        assert!(
            (EOP_DELTA_T_UNCERTAINTY_SECONDS - 0.005).abs() < 1e-12,
            "EOP_DELTA_T_UNCERTAINTY_SECONDS = {EOP_DELTA_T_UNCERTAINTY_SECONDS}, want 0.005"
        );
        // IERS 実測帯（< 0.1 s）。
        assert!(got < 0.1, "EOP-domain uncertainty {got} must be < 0.1 s");
    }

    /// 外挿域（EOP coverage 外）の不確実性は
    /// EspenakMeeusDeltaT.uncertainty_seconds(decimal_year(jd)) に一致。
    /// 変異: 外挿域でも固定定数 0.005 を返す、引数の取り違え、を殺す。
    #[test]
    fn composite_uncertainty_outside_eop_matches_espenak_meeus() {
        let comp = CompositeDeltaT::new(eop_2020());
        let utc = eop_utc(2030, 1, 1, 0, 0, 0.0); // 2020 coverage の外（外挿域）。
        let got = comp.uncertainty_seconds(utc);
        let want = EspenakMeeusDeltaT.uncertainty_seconds(decimal_year(utc.jd2()));
        assert!(
            (got - want).abs() < 1e-12,
            "extrapolated uncertainty = {got}, want EM {want}"
        );
        // 外挿域の値は固定定数 0.005 とは別物であること（分岐取り違えを更に固定）。
        assert!(
            (got - EOP_DELTA_T_UNCERTAINTY_SECONDS).abs() > 1e-6,
            "外挿域の不確実性 {got} は実測帯定数と区別できること"
        );
    }

    /// 単調性プロパティ: 外挿域では将来ほど不確実（accuracy.md §0）。
    /// 2100 相当 > 2030 相当（いずれも 2020 coverage の外 = 外挿域）。
    /// 外挿域で Espenak の年依存不確実性に従っていること（定数化する変異）を殺す。
    #[test]
    fn composite_uncertainty_grows_into_the_future_in_extrapolation() {
        let comp = CompositeDeltaT::new(eop_2020());
        let near = comp.uncertainty_seconds(eop_utc(2030, 1, 1, 0, 0, 0.0));
        let far = comp.uncertainty_seconds(eop_utc(2100, 1, 1, 0, 0, 0.0));
        assert!(
            far > near,
            "future uncertainty {far} (2100) must exceed nearer {near} (2030)"
        );
    }

    // ---- trait 経由呼び出し / Send + Sync -----------------------------

    /// CompositeDeltaT は DeltaTSource として（trait 越しに）呼べる。
    /// dyn DeltaTSource 経由で delta_t_seconds を呼び、恒等式オラクルと一致することを固定する
    /// （trait 実装の取り違え・未配線を殺す）。
    #[test]
    fn composite_callable_through_delta_t_source_trait() {
        let comp = CompositeDeltaT::new(eop_2020());
        let src: &dyn DeltaTSource = &comp;
        let utc = eop_utc(2020, 1, 1, 0, 0, 0.0);
        let got = src.delta_t_seconds(utc);
        assert!(
            (got - (37.0 + 32.184 - UT1_20200101)).abs() < 1e-9,
            "trait-object ΔT = {got}, want 69.3611222"
        );
        // 不確実性も trait 越しに供給される。
        assert!((src.uncertainty_seconds(utc) - EOP_DELTA_T_UNCERTAINTY_SECONDS).abs() < 1e-12);
    }

    /// ジェネリック越し（`fn call<S: DeltaTSource>`）でも呼べることを固定する。
    #[test]
    fn composite_callable_through_generic_bound() {
        fn delta_via<S: DeltaTSource>(src: &S, utc: UtcInstant) -> f64 {
            src.delta_t_seconds(utc)
        }
        let comp = CompositeDeltaT::new(eop_2020());
        let utc = eop_utc(2020, 1, 1, 0, 0, 0.0);
        let got = delta_via(&comp, utc);
        assert!((got - (37.0 + 32.184 - UT1_20200101)).abs() < 1e-9, "{got}");
    }

    /// `CompositeDeltaT<IersEopData>: Send + Sync`（DeltaTSource: Send + Sync 制約）の
    /// コンパイル時アサーション。trait 境界から Send/Sync を外す変異を殺す。
    #[test]
    fn composite_delta_t_is_send_sync() {
        fn _assert_send_sync<T: Send + Sync>() {}
        _assert_send_sync::<CompositeDeltaT<IersEopData>>();
    }
}
