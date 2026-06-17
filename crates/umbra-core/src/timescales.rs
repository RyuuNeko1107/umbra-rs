//! 時刻系データ束 [`TimeData`] と変換 facade [`TimeScales`]（`docs/issues/ISSUE-042` S2）。
//!
//! 閏秒テーブル（[`crate::LeapSecondTable`]）と IERS EOP（[`crate::IersEopData`]）を 1 つの
//! 純粋データ束 [`TimeData`] に束ね、[`TimeScales`] が UTC ↔ TT ↔ UT1 変換・極運動・ΔT 不確実性
//! 帯を Result 付きで供給する。core 純粋型（データは外部供給・同梱は別 crate）。
//!
//! ΔT は EOP coverage 内 ∧ 閏秒域では高精度（実測帯 <0.1 s）、範囲外は Espenak–Meeus 外挿
//! （[`crate::deltat`]）。閏秒は最終値を据え置く（≥1972 では `utc_to_tt` が常に Ok）一方、
//! UT1・極運動は EOP coverage 外で `MissingEarthOrientationData`（ISSUE-042 §46:「範囲外でも
//! TtInstant は返せる」＝将来食の TT は供給可、UTC 絶対は予測律速）。

use crate::angle::Radians;
use crate::constants::TT_MINUS_TAI_SECONDS;
use crate::deltat::{
    decimal_year, DeltaTModel, EspenakMeeusDeltaT, EOP_DELTA_T_UNCERTAINTY_SECONDS,
};
use crate::eop::{EarthOrientation, IersEopData};
use crate::error::TimeError;
use crate::metadata::DataSetMetadata;
use crate::time::{LeapSecondTable, TimeRange, TtInstant, Ut1Instant, UtcInstant};

const SECONDS_PER_DAY: f64 = 86_400.0;

/// 閏秒テーブル + IERS EOP を束ねた時刻系データ束（純粋型・データを持つ, 確定B3）。
///
/// ΔT 外挿器は固定の [`EspenakMeeusDeltaT`] なので別途保持しない（不確実性帯は導出）。
/// 同梱バイトでの構築 `bundled()` は umbra-ephemeris（`bundled-data` feature）側に置く（B3）。
#[derive(Clone, Debug)]
pub struct TimeData {
    leap: LeapSecondTable,
    eop: IersEopData,
}

impl TimeData {
    /// in-memory のデータ部品（閏秒テーブル・EOP）から束ねる。
    pub fn new(leap: LeapSecondTable, eop: IersEopData) -> Self {
        Self { leap, eop }
    }

    /// 全サービス（TT + UT1 + 極運動）が利用可能な範囲＝閏秒域と EOP coverage の積。
    ///
    /// 閏秒は上側に上限を持たない（最終値据え置き）ため上端は EOP coverage 終了が律速し、
    /// 下端は閏秒の最初の発効日と EOP coverage 開始の遅い方（max）。
    pub fn valid_range(&self) -> TimeRange<UtcInstant> {
        let leap_start = self.leap.earliest_utc();
        let eop_range = self.eop.coverage();
        let start = if leap_start.jd2().jd() >= eop_range.start.jd2().jd() {
            leap_start
        } else {
            eop_range.start
        };
        TimeRange {
            start,
            end: eop_range.end,
        }
    }

    /// 各データセットの provenance（閏秒・EOP の順）。
    pub fn metadata(&self) -> Vec<&DataSetMetadata> {
        vec![self.leap.metadata(), self.eop.metadata()]
    }

    /// 束ねた IERS EOP（`EclipseEngine` の `EarthOrientation` 供給に使う, ISSUE-043）。
    pub fn eop(&self) -> &IersEopData {
        &self.eop
    }

    /// 束ねた閏秒テーブル。
    pub fn leap(&self) -> &LeapSecondTable {
        &self.leap
    }
}

/// 時刻系変換 facade（006/007 を束ねる, 確定B3）。[`TimeData`] から構築し、変換は Result。
#[derive(Clone, Debug)]
pub struct TimeScales {
    data: TimeData,
}

impl TimeScales {
    /// [`TimeData`] から変換 facade を構築する。
    pub fn new(data: TimeData) -> Self {
        Self { data }
    }

    /// 構築に用いた [`TimeData`]（coverage/metadata 参照用）。
    pub fn data(&self) -> &TimeData {
        &self.data
    }

    /// UTC → TT。`TT = UTC + (TAI−UTC) + 32.184 s`。閏秒不足（1972 前）は
    /// [`TimeError::MissingLeapSecondData`]。最終閏秒以降は据え置きで未来も Ok。
    pub fn utc_to_tt(&self, t: UtcInstant) -> Result<TtInstant, TimeError> {
        let dat = self.data.leap.tai_minus_utc(t)?;
        Ok(TtInstant::from_jd2(
            t.jd2()
                .add_days((dat + TT_MINUS_TAI_SECONDS) / SECONDS_PER_DAY),
        ))
    }

    /// UTC → UT1。`UT1 = UTC + (UT1−UTC)`。EOP coverage 外は
    /// [`TimeError::MissingEarthOrientationData`]。
    pub fn utc_to_ut1(&self, t: UtcInstant) -> Result<Ut1Instant, TimeError> {
        crate::deltat::utc_to_ut1(t, &self.data.eop)
    }

    /// TT → UTC（[`Self::utc_to_tt`] の逆）。閏秒不足は [`TimeError::MissingLeapSecondData`]。
    ///
    /// `TT → TAI`（定数 −32.184 s）後、TAI を UTC とみなして ΔAT を引く（`time::tai_to_utc` と
    /// 同方式：閏秒挿入の前後 1 s 以内でのみ最大 1 s（閏秒 1 回分）ずれうるが報告用途では十分）。
    pub fn tt_to_utc(&self, t: TtInstant) -> Result<UtcInstant, TimeError> {
        let tai_jd = t.jd2().add_days(-TT_MINUS_TAI_SECONDS / SECONDS_PER_DAY);
        let dat = self.data.leap.tai_minus_utc(UtcInstant::from_jd2(tai_jd))?;
        Ok(UtcInstant::from_jd2(
            tai_jd.add_days(-dat / SECONDS_PER_DAY),
        ))
    }

    /// 極運動 (xp, yp)（CIRS→ITRS の極運動段, conventions §5）。EOP coverage 外は
    /// [`TimeError::MissingEarthOrientationData`]。
    pub fn polar_motion(&self, t: UtcInstant) -> Result<(Radians, Radians), TimeError> {
        self.data.eop.polar_motion(t)
    }

    /// ΔT の不確実性帯（秒, accuracy.md §0）。EOP coverage 内 ∧ 閏秒域は実測帯
    /// [`EOP_DELTA_T_UNCERTAINTY_SECONDS`]（<0.1 s）、それ以外は Espenak–Meeus 外挿の年依存帯
    /// （将来ほど増大）。`CalculationMetadata.delta_t_uncertainty_seconds` の源。
    pub fn delta_t_uncertainty_seconds(&self, t: UtcInstant) -> f64 {
        if self.data.leap.tai_minus_utc(t).is_ok() && self.data.eop.ut1_minus_utc(t).is_ok() {
            EOP_DELTA_T_UNCERTAINTY_SECONDS
        } else {
            EspenakMeeusDeltaT.uncertainty_seconds(decimal_year(t.jd2()))
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::deltat::{decimal_year, EspenakMeeusDeltaT};
    use crate::eop::{EopRecord, IersEopData};
    use crate::time::{tai_minus_utc, LeapSecondTable};
    use crate::timescales::{TimeData, TimeScales};
    use crate::{
        DataSetMetadata, DeltaTModel, Radians, TimeError, TtInstant, UtcInstant,
        EOP_DELTA_T_UNCERTAINTY_SECONDS,
    };

    // ====================================================================
    // 定数・補助（eop.rs / deltat.rs テストの verbatim オラクルを流用）
    // ====================================================================

    const SECONDS_PER_DAY: f64 = 86_400.0;
    /// arcsec → radian 係数（1" = π/648000 rad ≈ 4.84813681109536e-6）。独立計算。
    const ARCSEC_TO_RAD: f64 = std::f64::consts::PI / 648_000.0;

    /// IERS EOP 14 C04 実測オラクル（verbatim・eop.rs と同一値）。
    /// 1962-01-01 / MJD 37665（閏秒域外 = 恒等式不可・valid_range の閏秒 binding 検証用）。
    const MJD_1962: i32 = 37665;
    const UT1_1962: f64 = 0.0326338;
    const XP_1962: f64 = -0.012700;
    const YP_1962: f64 = 0.213000;
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

    /// provenance 完全な代表 EOP metadata（全フィールド非空）。
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

    /// 2 点 {58849, 58850}（2020-01-01/02）のみの EOP。coverage = 2020-01-01..2020-01-02。
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

    /// 3 点 {37665, 58849, 58850}。coverage 開始は 1962（閏秒起点 1972 より前）。
    /// valid_range の lower bound が「閏秒の最初の発効日(1972)」になることを試すための合成。
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

    /// 同梱閏秒（1972-01-01..2017-01-01 据え置き）＋ 2020 のみ EOP の代表 TimeData。
    fn time_data_2020() -> TimeData {
        TimeData::new(LeapSecondTable::bundled(), eop_2020())
    }

    /// 同梱閏秒 ＋ 1962 起点 EOP の TimeData（valid_range の閏秒 binding 検証用）。
    fn time_data_1962() -> TimeData {
        TimeData::new(LeapSecondTable::bundled(), eop_1962_to_2020())
    }

    /// 同梱閏秒 ＋ 2020 のみ EOP の代表 TimeScales。
    fn scales_2020() -> TimeScales {
        TimeScales::new(time_data_2020())
    }

    fn utc(y: i32, mo: u8, d: u8, h: u8, mi: u8, s: f64) -> UtcInstant {
        UtcInstant::from_gregorian(y, mo, d, h, mi, s).expect("valid calendar date")
    }

    // ====================================================================
    // TimeScales::utc_to_tt
    // ====================================================================

    /// utc_to_tt(2020-01-01) の TT−UTC = ΔAT(2020) + 32.184 = 37 + 32.184 = 69.184 s（exact）。
    /// 変異: 32.184 定数の取り違え、閏秒項の脱落、加算の符号（`+`→`-`）を殺す。
    #[test]
    fn utc_to_tt_offset_2020_is_69_184_exact() {
        let scales = scales_2020();
        let u = utc(2020, 1, 1, 0, 0, 0.0);
        let tt = scales.utc_to_tt(u).expect("2020 is in leap-second table");
        let diff_s = tt.jd2().days_since(u.jd2()) * SECONDS_PER_DAY;
        assert!(
            (diff_s - 69.184).abs() < 1e-6,
            "TT-UTC(2020) = {diff_s} s, want 69.184"
        );
    }

    /// 据え置き仕様の核心: utc_to_tt は **EOP coverage 外の未来でも Ok(TT) を返す**。
    /// 閏秒は最終値(37s, 2017以降)を据え置くため EOP に依存しない。2030・2100 とも Ok かつ
    /// TT−UTC = 69.184（同じ据え置き ΔAT=37）。
    /// 変異: utc_to_tt を EOP coverage に依存させて未来で Missing にする、据え置きを落とす、を殺す。
    #[test]
    fn utc_to_tt_returns_ok_beyond_eop_coverage() {
        let scales = scales_2020(); // EOP coverage は 2020 のみ。
        for (y, label) in [(2030, "2030"), (2100, "2100")] {
            let u = utc(y, 1, 1, 0, 0, 0.0);
            let tt = scales
                .utc_to_tt(u)
                .unwrap_or_else(|_| panic!("utc_to_tt({label}) must be Ok despite EOP coverage"));
            let diff_s = tt.jd2().days_since(u.jd2()) * SECONDS_PER_DAY;
            assert!(
                (diff_s - 69.184).abs() < 1e-6,
                "TT-UTC({label}) = {diff_s} s, want 69.184 (clamped ΔAT=37)"
            );
        }
    }

    /// utc_to_tt(1971-12-31) は閏秒テーブル最初の発効日(1972-01-01)より前 = MissingLeapSecondData。
    /// 変異: 下側ガードの脱落、別 variant（MissingEarthOrientationData 等）への取り違えを殺す。
    #[test]
    fn utc_to_tt_before_1972_is_missing_leap_second() {
        let scales = scales_2020();
        let pre = utc(1971, 12, 31, 0, 0, 0.0);
        assert_eq!(
            scales.utc_to_tt(pre).unwrap_err(),
            TimeError::MissingLeapSecondData
        );
    }

    // ====================================================================
    // TimeScales::utc_to_ut1
    // ====================================================================

    /// utc_to_ut1(2020-01-01 0h) は UTC に (UT1−UTC=-0.1771222)/86400 日を加算（exact）。
    /// 変異: 加算の符号（`+`→`-`）、/86400 の除数取り違え、UT1−UTC を引かない、対象 jd の取り違えを殺す。
    #[test]
    fn utc_to_ut1_adds_ut1_minus_utc_2020_exact() {
        let scales = scales_2020();
        let u = utc(2020, 1, 1, 0, 0, 0.0);
        let ut1 = scales.utc_to_ut1(u).expect("2020-01-01 is in EOP coverage");
        let want_jd = u.jd2().jd() + UT1_20200101 / SECONDS_PER_DAY;
        assert!(
            (ut1.jd2().jd() - want_jd).abs() < 1e-12,
            "ut1 jd = {}, want {want_jd}",
            ut1.jd2().jd()
        );
        // 差分を秒で取り直し、実測オラクル UT1−UTC に厳密一致することも固定する。
        let diff_s = ut1.jd2().days_since(u.jd2()) * SECONDS_PER_DAY;
        assert!(
            (diff_s - UT1_20200101).abs() < 1e-6,
            "diff = {diff_s} s, want {UT1_20200101} s"
        );
    }

    /// utc_to_ut1 は EOP coverage 外（2030, 据え置き対象外）では MissingEarthOrientationData。
    /// utc_to_tt が未来で Ok でも、ut1 は EOP に依存するため Missing になることを固定する
    /// （utc_to_tt と同じ据え置きロジックを ut1 にも誤適用する変異・EOP 不参照変異を殺す）。
    #[test]
    fn utc_to_ut1_outside_eop_coverage_is_missing() {
        let scales = scales_2020();
        let future = utc(2030, 1, 1, 0, 0, 0.0); // EOP coverage(2020) の外。
        assert_eq!(
            scales.utc_to_ut1(future).unwrap_err(),
            TimeError::MissingEarthOrientationData
        );
        // coverage より過去側（2020-01-02 の翌日）でも Missing（上側境界 `>` の固定）。
        let after = utc(2020, 1, 3, 0, 0, 0.0);
        assert_eq!(
            scales.utc_to_ut1(after).unwrap_err(),
            TimeError::MissingEarthOrientationData
        );
    }

    // ====================================================================
    // TimeScales::tt_to_utc
    // ====================================================================

    /// tt_to_utc(utc_to_tt(u)) ≈ u（ラウンドトリップ恒等, 1e-6 s 以内）。2020 で検証。
    /// 変異: tt_to_utc が utc_to_tt の正しい逆でない（符号・閏秒項取り違え）を殺す。
    #[test]
    fn tt_to_utc_round_trips_2020() {
        let scales = scales_2020();
        let u = utc(2020, 1, 1, 0, 0, 0.0);
        let tt = scales.utc_to_tt(u).expect("forward conversion is Ok");
        let back = scales.tt_to_utc(tt).expect("inverse conversion is Ok");
        let err_s = back.jd2().days_since(u.jd2()).abs() * SECONDS_PER_DAY;
        assert!(err_s < 1e-6, "round-trip error = {err_s} s, want < 1e-6");
    }

    /// tt_to_utc は閏秒域外（1971-12-31 相当の TT）では MissingLeapSecondData。
    /// TT を一旦 UTC とみなして閏秒を引く実装でも、1971 域は閏秒未定義で Missing になる。
    /// 変異: 下側ガードの脱落・別 variant への取り違えを殺す。
    #[test]
    fn tt_to_utc_before_1972_is_missing_leap_second() {
        let scales = scales_2020();
        // 1971-12-31 0h を TT スケールの JD として与える（閏秒テーブル域外）。
        let tt = TtInstant::from_jd2(utc(1971, 12, 31, 0, 0, 0.0).jd2());
        assert_eq!(
            scales.tt_to_utc(tt).unwrap_err(),
            TimeError::MissingLeapSecondData
        );
    }

    // ====================================================================
    // TimeScales::polar_motion
    // ====================================================================

    /// polar_motion(2020-01-01 0h) は (xp, yp) arcsec を rad 変換した値（exact, ~1e-12）。
    /// 変異: arcsec→rad 係数の脱落・取り違え、xp/yp の入れ替え、UT1 列を返す、を殺す。
    #[test]
    fn polar_motion_2020_in_radians_exact() {
        let scales = scales_2020();
        let (Radians(xp), Radians(yp)) = scales
            .polar_motion(utc(2020, 1, 1, 0, 0, 0.0))
            .expect("2020-01-01 is in EOP coverage");
        let want_xp = XP_20200101 * ARCSEC_TO_RAD;
        let want_yp = YP_20200101 * ARCSEC_TO_RAD;
        assert!((xp - want_xp).abs() < 1e-12, "xp = {xp}, want {want_xp}");
        assert!((yp - want_yp).abs() < 1e-12, "yp = {yp}, want {want_yp}");
    }

    /// polar_motion は EOP coverage 外（2030）では MissingEarthOrientationData。
    /// utc_to_tt の据え置きと違い、極運動は EOP 必須であることを固定する。
    #[test]
    fn polar_motion_outside_eop_coverage_is_missing() {
        let scales = scales_2020();
        let future = utc(2030, 1, 1, 0, 0, 0.0);
        assert_eq!(
            scales.polar_motion(future).unwrap_err(),
            TimeError::MissingEarthOrientationData
        );
    }

    // ====================================================================
    // TimeScales::delta_t_uncertainty_seconds
    // ====================================================================

    /// EOP coverage 内 かつ 閏秒域（2020-01-01）は固定定数 0.005 (=EOP_DELTA_T_UNCERTAINTY_SECONDS)
    /// に exact 一致し、かつ < 0.1（IERS 実測帯）。
    /// 変異: 定数値の取り違え、実測域で外挿器の不確実性(0.5 等)を返す、を殺す。
    #[test]
    fn delta_t_uncertainty_in_eop_domain_is_fixed_constant() {
        let scales = scales_2020();
        let got = scales.delta_t_uncertainty_seconds(utc(2020, 1, 1, 0, 0, 0.0));
        assert!(
            (got - EOP_DELTA_T_UNCERTAINTY_SECONDS).abs() < 1e-12,
            "EOP-domain uncertainty = {got}, want {EOP_DELTA_T_UNCERTAINTY_SECONDS}"
        );
        // 公開定数そのものが 0.005 であることも固定する。
        assert!(
            (EOP_DELTA_T_UNCERTAINTY_SECONDS - 0.005).abs() < 1e-12,
            "EOP_DELTA_T_UNCERTAINTY_SECONDS = {EOP_DELTA_T_UNCERTAINTY_SECONDS}, want 0.005"
        );
        assert!(got < 0.1, "EOP-domain uncertainty {got} must be < 0.1 s");
    }

    /// EOP coverage 外（2030）は EspenakMeeusDeltaT.uncertainty_seconds(decimal_year(t)) に一致し、
    /// かつ > 0.1（外挿域は実測帯より大きい）。
    /// 変異: 外挿域でも固定定数 0.005 を返す、decimal_year を渡さない引数取り違え、を殺す。
    #[test]
    fn delta_t_uncertainty_outside_eop_matches_espenak_meeus() {
        let scales = scales_2020();
        let t = utc(2030, 1, 1, 0, 0, 0.0); // EOP coverage(2020) の外。
        let got = scales.delta_t_uncertainty_seconds(t);
        let want = EspenakMeeusDeltaT.uncertainty_seconds(decimal_year(t.jd2()));
        assert!(
            (got - want).abs() < 1e-12,
            "extrapolated uncertainty = {got}, want EM {want}"
        );
        // 外挿域は実測帯定数 0.005 とは明確に区別され、> 0.1。
        assert!(got > 0.1, "extrapolated uncertainty {got} must be > 0.1 s");
        assert!(
            (got - EOP_DELTA_T_UNCERTAINTY_SECONDS).abs() > 1e-6,
            "外挿域の不確実性 {got} は実測帯定数と区別できること"
        );
    }

    /// 外挿域の単調性: 将来ほど不確実（2100 > 2030, いずれも EOP coverage 外）。
    /// 外挿域で年依存性を定数化する変異・大小取り違えを殺す。
    #[test]
    fn delta_t_uncertainty_grows_into_the_future() {
        let scales = scales_2020();
        let near = scales.delta_t_uncertainty_seconds(utc(2030, 1, 1, 0, 0, 0.0));
        let far = scales.delta_t_uncertainty_seconds(utc(2100, 1, 1, 0, 0, 0.0));
        assert!(
            far > near,
            "future uncertainty {far} (2100) must exceed nearer {near} (2030)"
        );
    }

    // ====================================================================
    // TimeData::valid_range
    // ====================================================================

    /// EOP が 2020 起点（閏秒起点 1972 より後）→ valid_range = [EOP開始, EOP終了]
    /// = [2020-01-01, 2020-01-02]。max(1972, 2020-01-01) = 2020-01-01 が下端。
    /// 変異: lower bound に閏秒起点を常に使う、start/end の取り違え、を殺す。
    #[test]
    fn valid_range_eop_start_after_leap_start() {
        let data = time_data_2020();
        let range = data.valid_range();
        let want_start = utc(2020, 1, 1, 0, 0, 0.0);
        let want_end = utc(2020, 1, 2, 0, 0, 0.0);
        assert!(
            (range.start.jd2().jd() - want_start.jd2().jd()).abs() < 1e-9,
            "valid_range.start jd = {}, want {} (2020-01-01)",
            range.start.jd2().jd(),
            want_start.jd2().jd()
        );
        assert!(
            (range.end.jd2().jd() - want_end.jd2().jd()).abs() < 1e-9,
            "valid_range.end jd = {}, want {} (2020-01-02)",
            range.end.jd2().jd(),
            want_end.jd2().jd()
        );
    }

    /// EOP が 1962 起点（閏秒起点 1972 より前）→ valid_range.start = 1972-01-01（閏秒 binding）。
    /// max(1972-01-01, 1962-01-01) = 1972-01-01。end は EOP 終了 2020-01-02 のまま。
    /// 変異: 下端に min を使う・EOP 開始をそのまま使う（閏秒下限を無視する）を殺す。
    #[test]
    fn valid_range_leap_start_is_binding_lower_bound() {
        let data = time_data_1962();
        let range = data.valid_range();
        let want_start = utc(1972, 1, 1, 0, 0, 0.0); // 閏秒の最初の発効日。
        let want_end = utc(2020, 1, 2, 0, 0, 0.0); // EOP 終了。
        assert!(
            (range.start.jd2().jd() - want_start.jd2().jd()).abs() < 1e-9,
            "valid_range.start jd = {}, want {} (1972-01-01, leap binding)",
            range.start.jd2().jd(),
            want_start.jd2().jd()
        );
        assert!(
            (range.end.jd2().jd() - want_end.jd2().jd()).abs() < 1e-9,
            "valid_range.end jd = {}, want {} (2020-01-02, EOP end)",
            range.end.jd2().jd(),
            want_end.jd2().jd()
        );
    }

    /// valid_range の整合: その下端は閏秒域・EOP coverage の両方が成立し、両 UT1/TT が引ける。
    /// 下端で utc_to_tt（閏秒）と utc_to_ut1（EOP）がともに Ok であることを確認し、
    /// valid_range が「全サービス利用可能範囲」であるという定義を固定する。
    #[test]
    fn valid_range_start_is_serviceable_by_all() {
        let data = time_data_1962();
        let scales = TimeScales::new(data);
        let start = scales.data().valid_range().start;
        assert!(
            scales.utc_to_tt(start).is_ok(),
            "valid_range.start must have leap-second data (utc_to_tt Ok)"
        );
        assert!(
            scales.utc_to_ut1(start).is_ok(),
            "valid_range.start must have EOP data (utc_to_ut1 Ok)"
        );
    }

    // ====================================================================
    // TimeData::metadata
    // ====================================================================

    /// metadata() は 2 件（閏秒, EOP の順）を返し、各 provenance が完全。
    /// 変異: 件数の取り違え、順序の入れ替え、provenance 不完全なメタデータの混入を殺す。
    #[test]
    fn metadata_returns_leap_then_eop_with_complete_provenance() {
        let data = time_data_2020();
        let md = data.metadata();
        assert_eq!(md.len(), 2, "metadata must contain exactly 2 datasets");
        // 順序: [0] = 閏秒, [1] = EOP。
        assert_eq!(
            md[0],
            LeapSecondTable::bundled().metadata(),
            "metadata[0] must be the leap-second dataset"
        );
        assert_eq!(
            md[1],
            &eop_metadata(),
            "metadata[1] must be the EOP dataset"
        );
        // 両方とも provenance 完全。
        assert!(
            md[0].has_complete_provenance(),
            "leap-second metadata must have complete provenance"
        );
        assert!(
            md[1].has_complete_provenance(),
            "EOP metadata must have complete provenance"
        );
    }

    // ====================================================================
    // TimeScales::data / 構築の透過性
    // ====================================================================

    /// TimeScales::data() は構築時に渡した TimeData をそのまま返す
    /// （valid_range と metadata が構築入力と一致することで透過性を固定）。
    /// 変異: data() が別の TimeData を返す・未配線を殺す。
    #[test]
    fn data_returns_constructed_time_data() {
        let data = time_data_2020();
        let want_range = data.valid_range();
        let scales = TimeScales::new(data);
        let got = scales.data();
        // valid_range が一致。
        assert!(
            (got.valid_range().start.jd2().jd() - want_range.start.jd2().jd()).abs() < 1e-9
                && (got.valid_range().end.jd2().jd() - want_range.end.jd2().jd()).abs() < 1e-9,
            "data() valid_range must match the constructed TimeData"
        );
        // metadata の件数・内容が一致。
        assert_eq!(got.metadata().len(), 2);
        assert_eq!(got.metadata()[0], LeapSecondTable::bundled().metadata());
        assert_eq!(got.metadata()[1], &eop_metadata());
    }

    // ====================================================================
    // 既存自由関数との整合（型化が観測挙動を変えないことの回帰固定）
    // ====================================================================

    /// utc_to_tt は閏秒域の代表点で free fn tai_minus_utc から組み立てた TT−UTC と一致する。
    /// 2020/2017/2000 の各 ΔAT で TT−UTC = ΔAT + 32.184 を独立に組み立て、facade と照合する
    /// （閏秒テーブル参照を facade が正しく配線していることを固定）。
    #[test]
    fn utc_to_tt_matches_leap_table_plus_32_184() {
        let scales = scales_2020();
        for (y, mo, d) in [(2020, 1, 1), (2017, 1, 1), (2000, 1, 1)] {
            let u = utc(y, mo, d, 0, 0, 0.0);
            let dat = tai_minus_utc(u).expect("in leap-second table");
            let want = dat + 32.184;
            let tt = scales.utc_to_tt(u).expect("Ok in leap domain");
            let got = tt.jd2().days_since(u.jd2()) * SECONDS_PER_DAY;
            assert!(
                (got - want).abs() < 1e-6,
                "TT-UTC({y}) = {got}, want {want}"
            );
        }
    }

    // ====================================================================
    // Send + Sync（コンパイル時アサーション）
    // ====================================================================

    /// `TimeData: Send + Sync`（純粋データ束）のコンパイル時アサーション。
    #[test]
    fn time_data_is_send_sync() {
        fn _assert_send_sync<T: Send + Sync>() {}
        _assert_send_sync::<TimeData>();
    }

    /// `TimeScales: Send + Sync`（スレッド間で共有可能な facade）のコンパイル時アサーション。
    #[test]
    fn time_scales_is_send_sync() {
        fn _assert_send_sync<T: Send + Sync>() {}
        _assert_send_sync::<TimeScales>();
    }
}
