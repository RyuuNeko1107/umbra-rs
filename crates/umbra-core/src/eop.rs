//! 地球姿勢（EOP）: [`EarthOrientation`] trait と IERS EOP C04 データ [`IersEopData`]
//! （`docs/issues/ISSUE-007` part P1）。
//!
//! UT1−UTC と極運動 (xp, yp) を日次 IERS EOP C04 テーブルから供給する。**本型はデータを
//! 持たない純粋型**（`from_records` で外部から受ける・確定B3）。同梱バイトを返す
//! `bundled()` は後続 part（umbra-ephemeris の `bundled-data` feature）で提供する。
//!
//! 単位・規約（ISSUE-007 §単位 / conventions §1）: 極運動は IERS が秒角（arcsec）配布のため
//! [`EopRecord`] は**arcsec を数値事実として保持**し、[`EarthOrientation::polar_motion`] の
//! 出力境界で [`crate::constants::ARCSEC_TO_RAD`] により [`Radians`] へ変換する。UT1−UTC は秒。
//! 補間は日次 MJD の線形補間。coverage 外は [`TimeError::MissingEarthOrientationData`]。
//!
//! 注（補間と閏秒）: C04 の UT1−UTC は閏秒境界で 1 s 跳ぶ。隣接日が同一閏秒区間内なら線形補間で
//! 妥当だが、閏秒を跨ぐ日対では UT1−TAI（連続量）での補間が本来適切。本 P1 は補間機構と
//! coverage/異常系の確立が目的で、機構は線形補間とする。
//! TODO(ISSUE-007 EOP part P2 / 同梱データ整備): 閏秒跨ぎ日対は UT1−TAI で補間し直す
//! （`MissingLeapSecondData` の閏秒テーブルと突合）。沈黙した誤補間値を残さない（accuracy.md §0）。

use crate::angle::Radians;
use crate::constants::ARCSEC_TO_RAD;
use crate::error::TimeError;
use crate::julian::JulianDate2;
use crate::metadata::DataSetMetadata;
use crate::time::{TimeRange, UtcInstant};

/// MJD = JD − 2400000.5（修正ユリウス日のオフセット）。
const MJD_JD_OFFSET: f64 = 2_400_000.5;

/// 地球回転に関わる姿勢量（UT1−UTC、極運動）の供給抽象（ISSUE-007 §公開IF）。
pub trait EarthOrientation: Send + Sync {
    /// `utc` における UT1−UTC（秒）。データ coverage 外は [`TimeError::MissingEarthOrientationData`]。
    fn ut1_minus_utc(&self, utc: UtcInstant) -> Result<f64, TimeError>;
    /// `utc` における極運動 (xp, yp)（ラジアン）。coverage 外は [`TimeError::MissingEarthOrientationData`]。
    fn polar_motion(&self, utc: UtcInstant) -> Result<(Radians, Radians), TimeError>;
}

/// IERS EOP C04 の 1 日次レコード（0h UTC）。極運動は arcsec で保持（数値事実・変換は出力境界）。
///
/// 前方互換のため `#[non_exhaustive]`（将来 LOD・各誤差列の追加余地, ISSUE-007 §数式）。外部 crate は
/// 構造体リテラルではなく [`EopRecord::new`] で構築する（後続 part の EOP パイプラインが消費）。
#[non_exhaustive]
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EopRecord {
    /// 修正ユリウス日（整数, 0h UTC）。
    pub mjd: i32,
    /// UT1 − UTC（秒）。
    pub ut1_minus_utc_s: f64,
    /// 極運動 x 成分（arcsec）。
    pub x_pole_arcsec: f64,
    /// 極運動 y 成分（arcsec）。
    pub y_pole_arcsec: f64,
}

impl EopRecord {
    /// 1 日次レコードを構築する（`#[non_exhaustive]` のため外部 crate はこの経路を使う）。
    pub fn new(mjd: i32, ut1_minus_utc_s: f64, x_pole_arcsec: f64, y_pole_arcsec: f64) -> Self {
        Self {
            mjd,
            ut1_minus_utc_s,
            x_pole_arcsec,
            y_pole_arcsec,
        }
    }
}

/// IERS EOP C04 データ（日次レコード＋系列版＋provenance/checksum）。純粋型（データは外部供給）。
#[derive(Clone, Debug, PartialEq)]
pub struct IersEopData {
    /// 厳密昇順・mjd 一意の日次レコード（非空）。
    records: Vec<EopRecord>,
    /// 系列版（例 `"EOP 14 C04"`）。
    series_version: String,
    /// 出所・完全性メタデータ。
    metadata: DataSetMetadata,
}

impl IersEopData {
    /// 日次レコードから構築する。`records` は**非空・mjd 厳密昇順（一意）**であること。
    ///
    /// 空・非昇順・重複 mjd は [`TimeError::InvalidEopData`]（不正な EOP テーブル）。
    pub fn from_records(
        records: Vec<EopRecord>,
        series_version: String,
        metadata: DataSetMetadata,
    ) -> Result<Self, TimeError> {
        if records.is_empty() {
            return Err(TimeError::InvalidEopData);
        }
        for pair in records.windows(2) {
            if pair[1].mjd <= pair[0].mjd {
                return Err(TimeError::InvalidEopData);
            }
        }
        Ok(Self {
            records,
            series_version,
            metadata,
        })
    }

    /// 系列版（`CalculationMetadata` へ伝播, ISSUE-007 §受け入れテスト）。
    pub fn series_version(&self) -> &str {
        &self.series_version
    }

    /// データの有効範囲 [最初の日 0h, 最後の日 0h]（UTC）。
    pub fn coverage(&self) -> TimeRange<UtcInstant> {
        let first = self.records.first().expect("non-empty by construction").mjd;
        let last = self.records.last().expect("non-empty by construction").mjd;
        TimeRange {
            start: mjd_to_utc(f64::from(first)),
            end: mjd_to_utc(f64::from(last)),
        }
    }

    /// 出所・完全性メタデータ。
    pub fn metadata(&self) -> &DataSetMetadata {
        &self.metadata
    }

    /// `utc` における (UT1−UTC, x_arcsec, y_arcsec) を日次線形補間で返す。範囲外は Missing。
    fn sample(&self, utc: UtcInstant) -> Result<(f64, f64, f64), TimeError> {
        let mjd = utc.jd2().jd() - MJD_JD_OFFSET;
        let first = f64::from(self.records.first().expect("non-empty").mjd);
        let last = f64::from(self.records.last().expect("non-empty").mjd);
        if mjd < first || mjd > last {
            return Err(TimeError::MissingEarthOrientationData);
        }
        // mjd_record ≤ mjd を満たすレコード数（昇順）。mjd ≥ first より ≥ 1。
        let count = self.records.partition_point(|r| f64::from(r.mjd) <= mjd);
        let lo = self.records[count - 1];
        // count == len ⟺ mjd == last（範囲チェックで mjd ≤ last 既知）。末尾日は上側レコードが無いため
        // ここで返す（`records[count]` の領域外参照を防ぐ）。厳密な非末尾日は下の補間で t=0 となり
        // lo の値に厳密一致するため、専用分岐は不要（冗長分岐を撤去）。
        if count == self.records.len() {
            return Ok((lo.ut1_minus_utc_s, lo.x_pole_arcsec, lo.y_pole_arcsec));
        }
        let hi = self.records[count];
        let span = f64::from(hi.mjd - lo.mjd);
        let t = (mjd - f64::from(lo.mjd)) / span;
        Ok((
            lo.ut1_minus_utc_s + t * (hi.ut1_minus_utc_s - lo.ut1_minus_utc_s),
            lo.x_pole_arcsec + t * (hi.x_pole_arcsec - lo.x_pole_arcsec),
            lo.y_pole_arcsec + t * (hi.y_pole_arcsec - lo.y_pole_arcsec),
        ))
    }
}

impl EarthOrientation for IersEopData {
    fn ut1_minus_utc(&self, utc: UtcInstant) -> Result<f64, TimeError> {
        let (ut1, _, _) = self.sample(utc)?;
        Ok(ut1)
    }

    fn polar_motion(&self, utc: UtcInstant) -> Result<(Radians, Radians), TimeError> {
        let (_, xp_arcsec, yp_arcsec) = self.sample(utc)?;
        Ok((
            Radians(xp_arcsec * ARCSEC_TO_RAD),
            Radians(yp_arcsec * ARCSEC_TO_RAD),
        ))
    }
}

/// 整数/小数 MJD（0h 起算）→ UTC 瞬時（JD = MJD + 2400000.5）。
fn mjd_to_utc(mjd: f64) -> UtcInstant {
    UtcInstant::from_jd2(JulianDate2::from_jd(mjd + MJD_JD_OFFSET))
}

#[cfg(test)]
mod tests {
    use crate::eop::{EarthOrientation, EopRecord, IersEopData};
    use crate::{DataSetMetadata, Radians, TimeError, UtcInstant};

    // ---- 定数・補助 -------------------------------------------------------

    /// arcsec → radian 係数（1" = π/648000 rad ≈ 4.84813681109536e-6）。
    const ARCSEC_TO_RAD: f64 = std::f64::consts::PI / 648_000.0;

    /// IERS EOP 14 C04 実測オラクル（verbatim・本文の表より）。
    /// 1962-01-01 / MJD 37665。
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

    fn rec_1962() -> EopRecord {
        EopRecord {
            mjd: MJD_1962,
            ut1_minus_utc_s: UT1_1962,
            x_pole_arcsec: XP_1962,
            y_pole_arcsec: YP_1962,
        }
    }
    fn rec_20200101() -> EopRecord {
        EopRecord {
            mjd: MJD_20200101,
            ut1_minus_utc_s: UT1_20200101,
            x_pole_arcsec: XP_20200101,
            y_pole_arcsec: YP_20200101,
        }
    }
    fn rec_20200102() -> EopRecord {
        EopRecord {
            mjd: MJD_20200102,
            ut1_minus_utc_s: UT1_20200102,
            x_pole_arcsec: XP_20200102,
            y_pole_arcsec: YP_20200102,
        }
    }

    /// provenance 完全な代表 metadata（全フィールド非空）。
    fn metadata() -> DataSetMetadata {
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

    /// 2 点 {58849, 58850} のみから成る補間用データセット
    /// （巨大ギャップを避け、隣接日での線形補間を検証する）。
    fn two_day_data() -> IersEopData {
        IersEopData::from_records(
            vec![rec_20200101(), rec_20200102()],
            "EOP 14 C04".to_string(),
            metadata(),
        )
        .expect("two adjacent ascending records build")
    }

    /// 3 点 {37665, 58849, 58850} の coverage / exact-lookup 用データセット。
    fn three_record_data() -> IersEopData {
        IersEopData::from_records(
            vec![rec_1962(), rec_20200101(), rec_20200102()],
            "EOP 14 C04".to_string(),
            metadata(),
        )
        .expect("three ascending records build")
    }

    fn utc(y: i32, mo: u8, d: u8, h: u8, mi: u8, s: f64) -> UtcInstant {
        UtcInstant::from_gregorian(y, mo, d, h, mi, s).expect("valid calendar date")
    }

    // ---- 1. コンストラクタ検証 -------------------------------------------

    /// 空の Vec は構築不可（Err）。
    #[test]
    fn from_records_rejects_empty() {
        let result = IersEopData::from_records(vec![], "EOP 14 C04".to_string(), metadata());
        assert!(result.is_err(), "empty records must be Err");
    }

    /// mjd 重複（非厳密増加）は構築不可（Err）。
    #[test]
    fn from_records_rejects_duplicate_mjd() {
        let result = IersEopData::from_records(
            vec![rec_20200101(), rec_20200101()],
            "EOP 14 C04".to_string(),
            metadata(),
        );
        assert!(result.is_err(), "duplicate mjd must be Err");
    }

    /// mjd 降順（昇順でない）は構築不可（Err）。
    #[test]
    fn from_records_rejects_descending_mjd() {
        let result = IersEopData::from_records(
            vec![rec_20200102(), rec_20200101()],
            "EOP 14 C04".to_string(),
            metadata(),
        );
        assert!(result.is_err(), "descending mjd must be Err");
    }

    /// 厳密昇順・非空なら構築できる（Ok）。
    #[test]
    fn from_records_accepts_ascending() {
        let result = IersEopData::from_records(
            vec![rec_20200101(), rec_20200102()],
            "EOP 14 C04".to_string(),
            metadata(),
        );
        assert!(result.is_ok(), "ascending unique records must be Ok");
    }

    // ---- 2. 厳密日ルックアップ -------------------------------------------

    /// レコードの MJD と一致する暦日では、その日の UT1−UTC をそのまま返す。
    #[test]
    fn exact_day_ut1_matches_record() {
        let data = three_record_data();
        // 2020-01-01 0h UTC = MJD 58849。
        let v = data.ut1_minus_utc(utc(2020, 1, 1, 0, 0, 0.0)).unwrap();
        assert!(
            (v - UT1_20200101).abs() < 1e-9,
            "exact-day UT1 = {v}, want {UT1_20200101}"
        );
    }

    /// 厳密日では極運動はその日の arcsec を rad へ変換した値を返す。
    #[test]
    fn exact_day_polar_motion_matches_record_in_radians() {
        let data = three_record_data();
        let (Radians(xp), Radians(yp)) = data.polar_motion(utc(2020, 1, 1, 0, 0, 0.0)).unwrap();
        let want_xp = XP_20200101 * ARCSEC_TO_RAD;
        let want_yp = YP_20200101 * ARCSEC_TO_RAD;
        assert!((xp - want_xp).abs() < 1e-12, "xp = {xp}, want {want_xp}");
        assert!((yp - want_yp).abs() < 1e-12, "yp = {yp}, want {want_yp}");
    }

    /// 1962-01-01（最初のレコード）でも厳密日ルックアップが成り立つ。
    #[test]
    fn exact_day_first_record_1962() {
        let data = three_record_data();
        let v = data.ut1_minus_utc(utc(1962, 1, 1, 0, 0, 0.0)).unwrap();
        assert!(
            (v - UT1_1962).abs() < 1e-9,
            "1962 UT1 = {v}, want {UT1_1962}"
        );
    }

    // ---- 3. 線形補間 -----------------------------------------------------

    /// 隣接 2 日の中点（MJD 58849.5 = 2020-01-01 12:00 UTC）で UT1−UTC が
    /// 両日の平均（中点では厳密に平均）になる。
    #[test]
    fn linear_interpolation_ut1_at_midpoint() {
        let data = two_day_data();
        let v = data.ut1_minus_utc(utc(2020, 1, 1, 12, 0, 0.0)).unwrap();
        let want = (UT1_20200101 + UT1_20200102) / 2.0; // -0.1773514
        assert!(
            (v - want).abs() < 1e-9,
            "midpoint UT1 = {v}, want mean {want}"
        );
    }

    /// 中点での極運動は両日 arcsec の平均を rad へ変換した値。
    #[test]
    fn linear_interpolation_polar_motion_at_midpoint() {
        let data = two_day_data();
        let (Radians(xp), Radians(yp)) = data.polar_motion(utc(2020, 1, 1, 12, 0, 0.0)).unwrap();
        let want_xp = ((XP_20200101 + XP_20200102) / 2.0) * ARCSEC_TO_RAD;
        let want_yp = ((YP_20200101 + YP_20200102) / 2.0) * ARCSEC_TO_RAD;
        assert!(
            (xp - want_xp).abs() < 1e-12,
            "midpoint xp = {xp}, want {want_xp}"
        );
        assert!(
            (yp - want_yp).abs() < 1e-12,
            "midpoint yp = {yp}, want {want_yp}"
        );
    }

    // ---- 4. coverage -----------------------------------------------------

    /// coverage は [最初の日 0h, 最後の日 0h]。
    /// {37665, 58849, 58850} → start=1962-01-01, end=2020-01-02。
    #[test]
    fn coverage_spans_first_to_last_record() {
        let data = three_record_data();
        let range = data.coverage();
        let want_start = utc(1962, 1, 1, 0, 0, 0.0);
        let want_end = utc(2020, 1, 2, 0, 0, 0.0);
        assert!(
            (range.start.jd2().jd() - want_start.jd2().jd()).abs() < 1e-9,
            "coverage start jd = {}, want {}",
            range.start.jd2().jd(),
            want_start.jd2().jd()
        );
        assert!(
            (range.end.jd2().jd() - want_end.jd2().jd()).abs() < 1e-9,
            "coverage end jd = {}, want {}",
            range.end.jd2().jd(),
            want_end.jd2().jd()
        );
    }

    // ---- 5. coverage 外 --------------------------------------------------

    /// 最初のレコードより前（1961-12-31）は両 API とも MissingEarthOrientationData。
    #[test]
    fn before_coverage_is_missing() {
        let data = three_record_data();
        let before = utc(1961, 12, 31, 0, 0, 0.0);
        assert_eq!(
            data.ut1_minus_utc(before).unwrap_err(),
            TimeError::MissingEarthOrientationData
        );
        assert_eq!(
            data.polar_motion(before).unwrap_err(),
            TimeError::MissingEarthOrientationData
        );
    }

    /// 最後のレコードより後（2020-01-03）は両 API とも MissingEarthOrientationData。
    #[test]
    fn after_coverage_is_missing() {
        let data = three_record_data();
        let after = utc(2020, 1, 3, 0, 0, 0.0);
        assert_eq!(
            data.ut1_minus_utc(after).unwrap_err(),
            TimeError::MissingEarthOrientationData
        );
        assert_eq!(
            data.polar_motion(after).unwrap_err(),
            TimeError::MissingEarthOrientationData
        );
    }

    // ---- 6. series_version / metadata 透過 -------------------------------

    /// series_version() は構築時に渡した文字列をそのまま返す。
    #[test]
    fn series_version_passthrough() {
        let data = two_day_data();
        assert_eq!(data.series_version(), "EOP 14 C04");
    }

    /// metadata() は構築時に渡した DataSetMetadata を返し、provenance が完全。
    #[test]
    fn metadata_passthrough_with_complete_provenance() {
        let data = two_day_data();
        let md = data.metadata();
        assert_eq!(md, &metadata());
        assert!(
            md.has_complete_provenance(),
            "metadata must have complete provenance"
        );
    }

    // ---- 7. Send + Sync --------------------------------------------------

    /// `IersEopData: Send + Sync`（trait 制約）のコンパイル時アサーション。
    #[test]
    fn iers_eop_data_is_send_sync() {
        fn _assert_send_sync<T: Send + Sync>() {}
        _assert_send_sync::<IersEopData>();
    }

    // ---- 8. 末尾日 inclusive ＋ 2 日ギャップ補間 -------------------------

    /// 最後のレコードと一致する暦日（2020-01-02 / MJD 58850）は **inclusive** で Ok を返し、
    /// その日のレコード値（UT1_20200102）に厳密一致する。極運動も Ok。
    /// （上側 coverage 境界 `mjd > last` が `>=` でないこと、および `count == len` 末尾分岐が
    /// 末尾レコードの値を返すことを固定する。）
    #[test]
    fn last_record_exact_returns_ok() {
        let data = two_day_data();
        let last_day = utc(2020, 1, 2, 0, 0, 0.0); // MJD 58850 = 最後のレコード日。
        let v = data
            .ut1_minus_utc(last_day)
            .expect("last coverage day must be inclusive (Ok)");
        assert!(
            (v - UT1_20200102).abs() < 1e-9,
            "last-day UT1 = {v}, want {UT1_20200102}"
        );
        assert!(
            data.polar_motion(last_day).is_ok(),
            "last coverage day polar_motion must be Ok"
        );
    }

    /// 2 日ギャップ（span = 2 日）の合成データで線形補間 **公式**（span と t）を固定する。
    /// レコード値は IERS オラクルではなく、補間機構を検算しやすくするための合成スキャフォールド
    /// （実アンカーの厳密日テストとは別目的：ここは FORMULA の検証）。
    ///
    /// records: MJD 58849 (2020-01-01) と MJD 58851 (2020-01-03)。間の MJD 58850 (2020-01-02) は
    /// **欠落**しているため、その日を引くと span = 2 日での補間となる。
    #[test]
    fn interpolation_over_two_day_gap() {
        let data = IersEopData::from_records(
            vec![
                EopRecord::new(58849, 10.0, 1.0, 2.0), // 2020-01-01
                EopRecord::new(58851, 20.0, 5.0, 6.0), // 2020-01-03（MJD 58850 は欠落）
            ],
            "EOP 14 C04".to_string(),
            metadata(),
        )
        .expect("two ascending records with a 2-day gap build");

        // 中点 MJD 58850 (2020-01-02), t = 0.5 → 10 と 20 の平均 = 15.0。
        // span を `/` に変異させると span=1.0 となり結果 20.0、t を `%`/`*` に変異させると
        // 20.0/30.0 となり、いずれも 15.0 と異なるため変異を殺す。
        let mid = data
            .ut1_minus_utc(utc(2020, 1, 2, 0, 0, 0.0))
            .expect("midpoint of 2-day gap is in coverage");
        assert!(
            (mid - 15.0).abs() < 1e-9,
            "2-day-gap midpoint UT1 = {mid}, want 15.0"
        );

        // 非中点 MJD 58849.5 (2020-01-01 12:00), t = 0.25 → 10 + 0.25×(20−10) = 12.5。
        // t≠0.5 で線形係数を固定する（中点では拾えない span/t 変異も含めて殺す）。
        let quarter = data
            .ut1_minus_utc(utc(2020, 1, 1, 12, 0, 0.0))
            .expect("quarter point of 2-day gap is in coverage");
        assert!(
            (quarter - 12.5).abs() < 1e-9,
            "2-day-gap quarter-point UT1 = {quarter}, want 12.5"
        );
    }
}
