//! 時刻系と変換 UTC ↔ TAI ↔ TT（`docs/issues/ISSUE-006`、`docs/algorithms/01-time-scales.md`）。
//!
//! 各時刻は対応する時刻系の [`JulianDate2`] として保持する。変換は一様時系 TAI を経由する:
//! `TAI = UTC + ΔAT(UTC)`（ΔAT=閏秒）、`TT = TAI + 32.184 s`（定数, conventions §4.1）。
//!
//! ΔAT は IERS の閏秒テーブル（1972– の公開事実データ）を組み込む。**1972 年より前は
//! 本テーブルでは未定義**で `TimeError::MissingLeapSecondData` を返す（その領域は ΔT 経由
//! = ISSUE-007 で扱う）。UT1/TDB は EOP/ΔT が必要なため本モジュールでは扱わない。

use crate::calendar::{gregorian_to_jd2, jd2_to_gregorian};
use crate::constants::TT_MINUS_TAI_SECONDS;
use crate::error::{DomainError, TimeError};
use crate::julian::JulianDate2;
use crate::metadata::DataSetMetadata;

const SECONDS_PER_DAY: f64 = 86_400.0;
/// MJD = JD − 2400000.5（修正ユリウス日のオフセット, `eop.rs` と同値）。
const MJD_JD_OFFSET: f64 = 2_400_000.5;

/// 同梱 IERS 閏秒（各 0h UTC で発効する暦日の MJD と、その日以降の TAI−UTC[s]）。1972–2017。
/// MJD は各発効日 0h UTC の修正ユリウス日（例 1972-01-01 = 41317, 2017-01-01 = 57754）。
const BUNDLED_LEAP_SECONDS: &[LeapSecondEntry] = &[
    LeapSecondEntry {
        mjd: 41317,
        tai_minus_utc_s: 10.0,
    }, // 1972-01-01
    LeapSecondEntry {
        mjd: 41499,
        tai_minus_utc_s: 11.0,
    }, // 1972-07-01
    LeapSecondEntry {
        mjd: 41683,
        tai_minus_utc_s: 12.0,
    }, // 1973-01-01
    LeapSecondEntry {
        mjd: 42048,
        tai_minus_utc_s: 13.0,
    }, // 1974-01-01
    LeapSecondEntry {
        mjd: 42413,
        tai_minus_utc_s: 14.0,
    }, // 1975-01-01
    LeapSecondEntry {
        mjd: 42778,
        tai_minus_utc_s: 15.0,
    }, // 1976-01-01
    LeapSecondEntry {
        mjd: 43144,
        tai_minus_utc_s: 16.0,
    }, // 1977-01-01
    LeapSecondEntry {
        mjd: 43509,
        tai_minus_utc_s: 17.0,
    }, // 1978-01-01
    LeapSecondEntry {
        mjd: 43874,
        tai_minus_utc_s: 18.0,
    }, // 1979-01-01
    LeapSecondEntry {
        mjd: 44239,
        tai_minus_utc_s: 19.0,
    }, // 1980-01-01
    LeapSecondEntry {
        mjd: 44786,
        tai_minus_utc_s: 20.0,
    }, // 1981-07-01
    LeapSecondEntry {
        mjd: 45151,
        tai_minus_utc_s: 21.0,
    }, // 1982-07-01
    LeapSecondEntry {
        mjd: 45516,
        tai_minus_utc_s: 22.0,
    }, // 1983-07-01
    LeapSecondEntry {
        mjd: 46247,
        tai_minus_utc_s: 23.0,
    }, // 1985-07-01
    LeapSecondEntry {
        mjd: 47161,
        tai_minus_utc_s: 24.0,
    }, // 1988-01-01
    LeapSecondEntry {
        mjd: 47892,
        tai_minus_utc_s: 25.0,
    }, // 1990-01-01
    LeapSecondEntry {
        mjd: 48257,
        tai_minus_utc_s: 26.0,
    }, // 1991-01-01
    LeapSecondEntry {
        mjd: 48804,
        tai_minus_utc_s: 27.0,
    }, // 1992-07-01
    LeapSecondEntry {
        mjd: 49169,
        tai_minus_utc_s: 28.0,
    }, // 1993-07-01
    LeapSecondEntry {
        mjd: 49534,
        tai_minus_utc_s: 29.0,
    }, // 1994-07-01
    LeapSecondEntry {
        mjd: 50083,
        tai_minus_utc_s: 30.0,
    }, // 1996-01-01
    LeapSecondEntry {
        mjd: 50630,
        tai_minus_utc_s: 31.0,
    }, // 1997-07-01
    LeapSecondEntry {
        mjd: 51179,
        tai_minus_utc_s: 32.0,
    }, // 1999-01-01
    LeapSecondEntry {
        mjd: 53736,
        tai_minus_utc_s: 33.0,
    }, // 2006-01-01
    LeapSecondEntry {
        mjd: 54832,
        tai_minus_utc_s: 34.0,
    }, // 2009-01-01
    LeapSecondEntry {
        mjd: 56109,
        tai_minus_utc_s: 35.0,
    }, // 2012-07-01
    LeapSecondEntry {
        mjd: 57204,
        tai_minus_utc_s: 36.0,
    }, // 2015-07-01
    LeapSecondEntry {
        mjd: 57754,
        tai_minus_utc_s: 37.0,
    }, // 2017-01-01
];

/// 1 件の閏秒エントリ: 発効する暦日（0h UTC）の MJD と、その日 0h 以降の TAI−UTC（整数秒）。
///
/// 前方互換のため `#[non_exhaustive]`（`eop.rs::EopRecord` と対称）。外部 crate は
/// 構造体リテラルではなく [`LeapSecondEntry::new`] で構築する。
#[non_exhaustive]
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LeapSecondEntry {
    /// 発効日 0h UTC の修正ユリウス日（整数）。
    pub mjd: i32,
    /// その日 0h 以降の TAI − UTC（秒, 整数値）。
    pub tai_minus_utc_s: f64,
}

impl LeapSecondEntry {
    /// 1 件の閏秒エントリを構築する（`#[non_exhaustive]` のため外部 crate はこの経路を使う）。
    pub fn new(mjd: i32, tai_minus_utc_s: f64) -> Self {
        Self {
            mjd,
            tai_minus_utc_s,
        }
    }
}

/// IERS 閏秒テーブル（1972–, ΔAT = TAI − UTC）。発効日 MJD 厳密昇順・非空。
///
/// `eop.rs::IersEopData` と対称の純粋データ型（データは外部供給か同梱定数から）。
/// 最初の発効日より前は [`TimeError::MissingLeapSecondData`]、最終エントリ以降は
/// 最後の値を据え置く（閏秒は約半年前に告知され将来予報では最終値が最良推定）。
#[derive(Clone, Debug, PartialEq)]
pub struct LeapSecondTable {
    /// 厳密昇順・mjd 一意の発効エントリ（非空）。
    entries: Vec<LeapSecondEntry>,
    /// 出所・完全性メタデータ。
    metadata: DataSetMetadata,
}

impl LeapSecondTable {
    /// 発効エントリから構築する。`entries` は**非空・mjd 厳密昇順（一意）**であること。
    ///
    /// 空・非昇順・重複 mjd は [`TimeError::InvalidLeapSecondData`]。
    pub fn from_entries(
        entries: Vec<LeapSecondEntry>,
        metadata: DataSetMetadata,
    ) -> Result<Self, TimeError> {
        if entries.is_empty() {
            return Err(TimeError::InvalidLeapSecondData);
        }
        for pair in entries.windows(2) {
            if pair[1].mjd <= pair[0].mjd {
                return Err(TimeError::InvalidLeapSecondData);
            }
        }
        Ok(Self { entries, metadata })
    }

    /// 同梱 IERS 閏秒（1972-01-01 .. 2017-01-01, TAI−UTC 10 .. 37 s）。
    ///
    /// 同梱データは公開事実（IERS Bulletin C / IANA `leap-seconds.list`）。checksum は
    /// packed-LE-f64 `[n, then per entry: mjd, tai_minus_utc_s]` の SHA-256（EOP と同形式）。
    pub fn bundled() -> Self {
        Self::from_entries(BUNDLED_LEAP_SECONDS.to_vec(), bundled_leap_metadata())
            .expect("bundled leap-second entries are non-empty and strictly ascending")
    }

    /// `utc` における TAI − UTC（秒）。最初の発効日より前は [`TimeError::MissingLeapSecondData`]。
    pub fn tai_minus_utc(&self, utc: UtcInstant) -> Result<f64, TimeError> {
        lookup_tai_minus_utc(&self.entries, utc)
    }

    /// テーブル最初の発効日（0h UTC）。これより前は `tai_minus_utc` が Missing を返す
    /// （[`TimeData::valid_range`](crate::timescales::TimeData::valid_range) の下限算出に使う）。
    pub fn earliest_utc(&self) -> UtcInstant {
        let mjd = self.entries.first().expect("non-empty by construction").mjd;
        UtcInstant::from_jd2(JulianDate2::from_jd(f64::from(mjd) + MJD_JD_OFFSET))
    }

    /// 出所・完全性メタデータ。
    pub fn metadata(&self) -> &DataSetMetadata {
        &self.metadata
    }
}

/// 同梱閏秒の provenance。`checksum` は packed-LE-f64 `[n, then per entry: mjd, tai_minus_utc_s]`
/// の SHA-256 を**オフライン算出**した値（EOP `eopc04_14.bin` と同形式）。本 S1 では閏秒は core の
/// 定数 [`BUNDLED_LEAP_SECONDS`] でありファイル生成物ではないため、xtask による再生成・回帰検証
/// 経路は未整備（EOP の `verify-generated` 相当は後続で閏秒を生成物化する場合に追加）。
fn bundled_leap_metadata() -> DataSetMetadata {
    DataSetMetadata {
        name: "iers-leap-seconds".to_string(),
        version: "IERS Bulletin C (1972-2017)".to_string(),
        source: "IERS / IANA leap-seconds.list".to_string(),
        license: "public-domain".to_string(),
        valid_from: "1972-01-01".to_string(),
        valid_to: "2017-01-01".to_string(),
        checksum: "13ddc355a79126910081e0a33d1856dde530ca795e891c606d16899d92cd06bf".to_string(),
    }
}

/// 発効日 0h UTC の MJD → 閾値 JD。
fn leap_threshold_jd(entry: &LeapSecondEntry) -> f64 {
    f64::from(entry.mjd) + MJD_JD_OFFSET
}

/// 昇順エントリ列に対する TAI−UTC ルックアップ（自由関数・[`LeapSecondTable`] 双方が共用）。
///
/// 最初の発効日より前は [`TimeError::MissingLeapSecondData`]。最終エントリ以降は最後の値を据え置く。
fn lookup_tai_minus_utc(entries: &[LeapSecondEntry], utc: UtcInstant) -> Result<f64, TimeError> {
    let jd = utc.0.jd();
    let first = entries.first().expect("leap-second entries are non-empty");
    if jd < leap_threshold_jd(first) {
        return Err(TimeError::MissingLeapSecondData);
    }
    let mut dat = first.tai_minus_utc_s;
    for entry in entries {
        if jd >= leap_threshold_jd(entry) {
            dat = entry.tai_minus_utc_s;
        } else {
            break;
        }
    }
    Ok(dat)
}

/// 協定世界時 UTC の瞬時。
#[derive(Clone, Copy, Debug, PartialEq, PartialOrd)]
pub struct UtcInstant(JulianDate2);
/// 国際原子時 TAI の瞬時。
#[derive(Clone, Copy, Debug, PartialEq, PartialOrd)]
pub struct TaiInstant(JulianDate2);
/// 地球時 TT の瞬時。
#[derive(Clone, Copy, Debug, PartialEq, PartialOrd)]
pub struct TtInstant(JulianDate2);

/// JD を **0.1 秒へ丸めた**暦本体文字列 `YYYY-MM-DDThh:mm:ss.s`（時刻系マーカ無し）を作る
/// （iso 表示用, ISSUE-031 S31b）。
///
/// `jd2_to_gregorian` は丸め境界で ±eps を返しうる（例: 16:00:00 を 15:59:59.9995 と返す）。
/// これを `{:.1}` に素通しすると `15:59:60.0` のような不正表記になる。そこで秒を **整数 1/10 秒**
/// （`tenths`, 0..=600）へ丸め、60.0 到達分を 分→時→日 へ繰り上げる。日跨ぎ（23:59:59.95＋）は
/// JD の整数日に +1 して年月日のみ再導出する（時刻成分は 00:00:00）。整数で組み立てるため
/// `{:.1}` の浮動小数丸め（0.05→"0.1" 化）も回避する。lossless 値は別途 `jd` フィールドが持つ。
#[cfg(feature = "serde")]
#[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
fn iso_body_tenths(jd: JulianDate2) -> String {
    let (mut y, mut mo, mut d, mut h, mut mi, s) = jd2_to_gregorian(jd);
    // 秒を 1/10 秒（整数）へ丸める。s ∈ [0,60) なので tenths ∈ 0..=600。
    let mut tenths = (s * 10.0).round() as i64;
    if tenths >= 600 {
        tenths -= 600;
        mi += 1;
    }
    if mi >= 60 {
        mi -= 60;
        h += 1;
    }
    if h >= 24 {
        h -= 24;
        // 日跨ぎ: 整数日を 1 進めて年月日のみ採用（時刻は繰り上げ済みの 00:00:00.x）。
        let (ny, nmo, nd, ..) = jd2_to_gregorian(JulianDate2::new(jd.part1 + 1.0, jd.part2));
        y = ny;
        mo = nmo;
        d = nd;
    }
    let s_whole = (tenths / 10) as u8;
    let s_tenth = (tenths % 10) as u8;
    format!("{y:04}-{mo:02}-{d:02}T{h:02}:{mi:02}:{s_whole:02}.{s_tenth}")
}

/// `UtcInstant` の JSON 表現（ISSUE-031 S31b・api-draft §0/A7）。自己記述かつ可逆な
/// `{ "iso": <暦形式・末尾 Z>, "jd": { "part1", "part2" } }`。iso は人間可読の表示チャネル
/// （秒は 0.1 秒へ丸め）、jd は lossless チャネル（2 要素 JD をそのまま）。Serialize のみ（§0）。
#[cfg(feature = "serde")]
impl serde::Serialize for UtcInstant {
    fn serialize<S: serde::Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        use serde::ser::SerializeStruct;
        // UTC はオフセット 0 を表す `Z` を付ける（RFC3339）。
        let iso = format!("{}Z", iso_body_tenths(self.jd2()));
        let mut st = serializer.serialize_struct("UtcInstant", 2)?;
        st.serialize_field("iso", &iso)?;
        st.serialize_field("jd", &self.jd2())?;
        st.end()
    }
}

/// `TtInstant` の JSON 表現（ISSUE-031 S31b）。`{ "iso": <暦形式・Z なし>, "jd": {..} }`。
/// TT は UTC ではないため iso に UTC マーカ `Z` を付けない（時刻系はフィールド名 `time_tt` が表す）。
#[cfg(feature = "serde")]
impl serde::Serialize for TtInstant {
    fn serialize<S: serde::Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        use serde::ser::SerializeStruct;
        let iso = iso_body_tenths(self.jd2());
        let mut st = serializer.serialize_struct("TtInstant", 2)?;
        st.serialize_field("iso", &iso)?;
        st.serialize_field("jd", &self.jd2())?;
        st.end()
    }
}

impl UtcInstant {
    /// UTC スケールの JD から構築。
    pub fn from_jd2(jd: JulianDate2) -> Self {
        UtcInstant(jd)
    }
    /// グレゴリオ暦（UTC）から構築。
    pub fn from_gregorian(
        year: i32,
        month: u8,
        day: u8,
        hour: u8,
        minute: u8,
        second: f64,
    ) -> Result<Self, DomainError> {
        Ok(UtcInstant(gregorian_to_jd2(
            year, month, day, hour, minute, second,
        )?))
    }
    /// UTC スケールの JD。
    pub fn jd2(self) -> JulianDate2 {
        self.0
    }
    /// グレゴリオ暦（UTC）へ。
    pub fn to_gregorian(self) -> (i32, u8, u8, u8, u8, f64) {
        jd2_to_gregorian(self.0)
    }
}

impl TaiInstant {
    /// TAI スケールの JD から構築。
    pub fn from_jd2(jd: JulianDate2) -> Self {
        TaiInstant(jd)
    }
    /// TAI スケールの JD。
    pub fn jd2(self) -> JulianDate2 {
        self.0
    }
}

impl TtInstant {
    /// TT スケールの JD から構築。
    pub fn from_jd2(jd: JulianDate2) -> Self {
        TtInstant(jd)
    }
    /// TT スケールの JD。
    pub fn jd2(self) -> JulianDate2 {
        self.0
    }
}

/// 世界時 UT1 の瞬時（地球回転。TT − ΔT。ΔT は `crate::deltat` 参照）。
#[derive(Clone, Copy, Debug, PartialEq, PartialOrd)]
pub struct Ut1Instant(JulianDate2);

impl Ut1Instant {
    /// UT1 スケールの JD から構築。
    pub fn from_jd2(jd: JulianDate2) -> Self {
        Ut1Instant(jd)
    }
    /// UT1 スケールの JD。
    pub fn jd2(self) -> JulianDate2 {
        self.0
    }
}

/// 太陽系力学時 TDB の瞬時（Reference 暦用。TT との差は周期項で最大 ~1.7 ms。
/// TT↔TDB 変換は精度が要る段階で別途実装する）。
#[derive(Clone, Copy, Debug, PartialEq, PartialOrd)]
pub struct TdbInstant(JulianDate2);

impl TdbInstant {
    /// TDB スケールの JD から構築。
    pub fn from_jd2(jd: JulianDate2) -> Self {
        TdbInstant(jd)
    }
    /// TDB スケールの JD。
    pub fn jd2(self) -> JulianDate2 {
        self.0
    }
}

/// 時刻範囲 `[start, end]`（任意の時刻型に対する区間）。
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct TimeRange<T> {
    /// 開始。
    pub start: T,
    /// 終了。
    pub end: T,
}

/// 時間区間（フィット区間など。`TimeRange` と別用途で使い分ける）。
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize))]
pub struct TimeInterval<T> {
    /// 開始。
    pub start: T,
    /// 終了。
    pub end: T,
}

/// ΔAT = TAI − UTC（秒）。1972 年より前は `MissingLeapSecondData`。
///
/// 同梱閏秒（[`BUNDLED_LEAP_SECONDS`]）に対する薄いラッパ。データ型を取る版は
/// [`LeapSecondTable::tai_minus_utc`]（両者は同一の [`lookup_tai_minus_utc`] を共用）。
/// 最終エントリ以降は最後の値を据え置く（閏秒は約半年前に告知され、2017 以降は増えていない）。
pub fn tai_minus_utc(utc: UtcInstant) -> Result<f64, TimeError> {
    lookup_tai_minus_utc(BUNDLED_LEAP_SECONDS, utc)
}

/// UTC → TAI。
pub fn utc_to_tai(utc: UtcInstant) -> Result<TaiInstant, TimeError> {
    let dat = tai_minus_utc(utc)?;
    Ok(TaiInstant(utc.0.add_days(dat / SECONDS_PER_DAY)))
}

/// TAI → TT（定数 +32.184 s）。
pub fn tai_to_tt(tai: TaiInstant) -> TtInstant {
    TtInstant(tai.0.add_days(TT_MINUS_TAI_SECONDS / SECONDS_PER_DAY))
}

/// TT → TAI（定数 −32.184 s）。
pub fn tt_to_tai(tt: TtInstant) -> TaiInstant {
    TaiInstant(tt.0.add_days(-TT_MINUS_TAI_SECONDS / SECONDS_PER_DAY))
}

/// UTC → TT。
pub fn utc_to_tt(utc: UtcInstant) -> Result<TtInstant, TimeError> {
    Ok(tai_to_tt(utc_to_tai(utc)?))
}

/// TAI → UTC。ΔAT は UTC 依存だが、tai を UTC とみなして ΔAT を引く単純法を用いる
/// （閏秒挿入の前後 1 s 以内でのみ最大 1 s ずれうる。出力用途では十分。conventions / 要確認）。
pub fn tai_to_utc(tai: TaiInstant) -> Result<UtcInstant, TimeError> {
    let dat = tai_minus_utc(UtcInstant(tai.0))?;
    Ok(UtcInstant(tai.0.add_days(-dat / SECONDS_PER_DAY)))
}

/// TT → UTC。
pub fn tt_to_utc(tt: TtInstant) -> Result<UtcInstant, TimeError> {
    tai_to_utc(tt_to_tai(tt))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::metadata::DataSetMetadata;

    fn utc(y: i32, mo: u8, d: u8) -> UtcInstant {
        UtcInstant::from_gregorian(y, mo, d, 0, 0, 0.0).unwrap()
    }

    // ===== ISSUE-042 S1: LeapSecondTable / LeapSecondEntry のデータ型化 =====

    /// LeapSecondTable / from_entries 用の provenance 完全なメタデータ
    /// （eop.rs の metadata() ヘルパに倣う。全フィールド非空）。
    fn leap_metadata() -> DataSetMetadata {
        DataSetMetadata {
            name: "iers-leap-seconds".to_string(),
            version: "IERS Bulletin C".to_string(),
            source: "IERS Earth Orientation Center, datacenter.iers.org".to_string(),
            license: "public-domain".to_string(),
            valid_from: "1972-01-01".to_string(),
            valid_to: "2017-01-01".to_string(),
            checksum: "0000000000000000000000000000000000000000000000000000000000000000"
                .to_string(),
        }
    }

    // ---- LeapSecondEntry::new -------------------------------------------

    /// `LeapSecondEntry::new(mjd, tai_minus_utc_s)` が各フィールドを所定の位置へ
    /// 設定する。引数の取り違え（mjd と秒の入れ替え）や片方無視の変異を殺す。
    #[test]
    fn leap_second_entry_new_sets_fields() {
        let e = LeapSecondEntry::new(41317, 10.0);
        assert_eq!(e.mjd, 41317, "mjd field must be the first argument");
        assert_eq!(
            e.tai_minus_utc_s, 10.0,
            "tai_minus_utc_s field must be the second argument"
        );
    }

    // ---- LeapSecondTable::from_entries（検証） ---------------------------

    /// 空 Vec は構築不可（Err(InvalidLeapSecondData)）。
    /// 「非空」検査の削除変異を殺す。
    #[test]
    fn from_entries_rejects_empty() {
        let result = LeapSecondTable::from_entries(vec![], leap_metadata());
        assert_eq!(
            result.unwrap_err(),
            TimeError::InvalidLeapSecondData,
            "empty entries must be Err(InvalidLeapSecondData)"
        );
    }

    /// mjd 降順（昇順でない）は構築不可（Err(InvalidLeapSecondData)）。
    /// 昇順検査の比較の向きや検査削除の変異を殺す。
    #[test]
    fn from_entries_rejects_descending_mjd() {
        let result = LeapSecondTable::from_entries(
            vec![
                LeapSecondEntry::new(41499, 11.0),
                LeapSecondEntry::new(41317, 10.0),
            ],
            leap_metadata(),
        );
        assert_eq!(
            result.unwrap_err(),
            TimeError::InvalidLeapSecondData,
            "descending mjd must be Err(InvalidLeapSecondData)"
        );
    }

    /// mjd 重複（厳密増加でない）は構築不可（Err(InvalidLeapSecondData)）。
    /// 「`<=`」を「`<`」に弱める（重複を許す）変異を殺す。
    #[test]
    fn from_entries_rejects_duplicate_mjd() {
        let result = LeapSecondTable::from_entries(
            vec![
                LeapSecondEntry::new(41317, 10.0),
                LeapSecondEntry::new(41317, 11.0),
            ],
            leap_metadata(),
        );
        assert_eq!(
            result.unwrap_err(),
            TimeError::InvalidLeapSecondData,
            "duplicate mjd must be Err(InvalidLeapSecondData)"
        );
    }

    /// 厳密昇順・非空なら構築できる（Ok）。
    /// 上記の異常系テストが「常に Err」変異で通ってしまうのを防ぐ正常系の対。
    #[test]
    fn from_entries_accepts_ascending_nonempty() {
        let result = LeapSecondTable::from_entries(
            vec![
                LeapSecondEntry::new(41317, 10.0),
                LeapSecondEntry::new(41499, 11.0),
            ],
            leap_metadata(),
        );
        assert!(
            result.is_ok(),
            "ascending unique non-empty entries must be Ok"
        );
    }

    // ---- LeapSecondTable::tai_minus_utc（既知値） -----------------------

    /// 同梱テーブルの既知 ΔAT 値（IERS の公開事実）を固定する。
    /// 2020=37, 2000-01-01=32, 1985-07-01=23, 1972-07-01=11, 1972-01-01=10。
    /// テーブル参照・ステップ選択の取り違え変異を殺す。
    #[test]
    fn bundled_tai_minus_utc_known_values() {
        let t = LeapSecondTable::bundled();
        assert_eq!(t.tai_minus_utc(utc(2020, 6, 1)).unwrap(), 37.0);
        assert_eq!(t.tai_minus_utc(utc(2000, 1, 1)).unwrap(), 32.0);
        assert_eq!(t.tai_minus_utc(utc(1985, 7, 1)).unwrap(), 23.0);
        assert_eq!(t.tai_minus_utc(utc(1972, 7, 1)).unwrap(), 11.0);
        assert_eq!(t.tai_minus_utc(utc(1972, 1, 1)).unwrap(), 10.0);
    }

    /// 発効境界 2017-01-01 0h で 36 → 37 へ跳ぶ。境界当日は新しい値（37）。
    /// 比較が `<` か `<=` か（境界当日を新旧どちらに含めるか）を固定する。
    /// `>=` を `>` に弱めると当日 0h が 36 のままになり、これを殺す。
    #[test]
    fn bundled_tai_minus_utc_steps_exactly_on_boundary() {
        let t = LeapSecondTable::bundled();
        assert_eq!(
            t.tai_minus_utc(utc(2016, 12, 31)).unwrap(),
            36.0,
            "day before effective date keeps the old value"
        );
        assert_eq!(
            t.tai_minus_utc(utc(2017, 1, 1)).unwrap(),
            37.0,
            "effective day 0h takes the new value"
        );
    }

    /// 最終エントリ（2017-01-01）以降は最後の値 37 を据え置く（将来予報用）。
    /// 2026・2100 のような遠い未来でも 37。範囲外を Missing にしてしまう変異や
    /// 末尾据え置きを落とす変異を殺す。
    #[test]
    fn bundled_tai_minus_utc_clamps_to_last_value() {
        let t = LeapSecondTable::bundled();
        assert_eq!(t.tai_minus_utc(utc(2026, 1, 1)).unwrap(), 37.0);
        assert_eq!(t.tai_minus_utc(utc(2100, 1, 1)).unwrap(), 37.0);
    }

    /// 最初の発効日（1972-01-01 0h）より前（1971-12-31）は Missing。
    /// 下側ガードの削除や比較の向きの変異を殺す。
    #[test]
    fn bundled_tai_minus_utc_before_1972_is_missing() {
        let t = LeapSecondTable::bundled();
        let pre = UtcInstant::from_gregorian(1971, 12, 31, 0, 0, 0.0).unwrap();
        assert_eq!(
            t.tai_minus_utc(pre).unwrap_err(),
            TimeError::MissingLeapSecondData
        );
    }

    // ---- bundled() の主要エントリ・件数 ---------------------------------

    /// 同梱エントリの代表点（MJD, ΔAT秒）を MJD レベルで固定する。
    /// from_gregorian で構築した既知日の MJD（= jd2().jd() − 2400000.5）が
    /// 仕様値（1972-01-01=41317, 1972-07-01=41499, 2017-01-01=57754）と一致することで、
    /// 既知値テストが正しい暦日を引いていることを保証する（オラクル MJD の独立確認）。
    #[test]
    #[allow(clippy::cast_possible_truncation)] // MJD は小整数域、round 後の i32 化は安全。
    fn known_mjd_anchors() {
        const MJD_JD_OFFSET: f64 = 2_400_000.5;
        let mjd = |y, mo, d| (utc(y, mo, d).jd2().jd() - MJD_JD_OFFSET).round() as i32;
        assert_eq!(mjd(1972, 1, 1), 41317);
        assert_eq!(mjd(1972, 7, 1), 41499);
        assert_eq!(mjd(2017, 1, 1), 57754);
    }

    // ---- bundled() の metadata ------------------------------------------

    /// bundled() の metadata は provenance 完全で、valid_from/valid_to が
    /// 同梱範囲（1972-01-01 .. 2017-01-01）を表す。
    /// metadata 欠落・年代取り違えの変異を殺す。
    #[test]
    fn bundled_metadata_has_complete_provenance_and_range() {
        let t = LeapSecondTable::bundled();
        let md = t.metadata();
        assert!(
            md.has_complete_provenance(),
            "bundled metadata must have complete provenance"
        );
        assert_eq!(
            md.valid_from, "1972-01-01",
            "valid_from must be the first effective date"
        );
        assert_eq!(
            md.valid_to, "2017-01-01",
            "valid_to must be the last effective date"
        );
    }

    // ---- 同値性（型化リファクタが挙動を変えないことの固定・最重要） ----

    /// 任意の代表 utc で free fn `tai_minus_utc` と `bundled().tai_minus_utc` が
    /// 一致する（Ok 値も Err バリアントも）。境界・据え置き・1972前の Missing を含む。
    /// 自由関数を別実装へ差し替える/据え置きやガードの分岐を一方だけ変える変異を殺す
    /// （型化が観測可能な挙動を変えないことの回帰固定）。
    #[test]
    fn bundled_matches_free_function_everywhere() {
        let table = LeapSecondTable::bundled();
        let samples = [
            UtcInstant::from_gregorian(1971, 12, 31, 0, 0, 0.0).unwrap(), // Missing
            utc(1972, 1, 1),                                              // 最初の発効日
            utc(1972, 7, 1),
            utc(1985, 7, 1),
            utc(2000, 1, 1),
            utc(2016, 12, 31), // 境界の前日
            utc(2017, 1, 1),   // 最終エントリ当日（跳び）
            utc(2020, 6, 1),
            utc(2026, 1, 1), // 末尾据え置き
            utc(2100, 1, 1),
        ];
        for u in samples {
            assert_eq!(
                table.tai_minus_utc(u),
                tai_minus_utc(u),
                "table and free function must agree at jd = {}",
                u.jd2().jd()
            );
        }
    }

    #[test]
    fn delta_at_known_values() {
        assert_eq!(tai_minus_utc(utc(2020, 6, 1)).unwrap(), 37.0);
        assert_eq!(tai_minus_utc(utc(2017, 1, 1)).unwrap(), 37.0);
        assert_eq!(tai_minus_utc(utc(2016, 12, 31)).unwrap(), 36.0);
        assert_eq!(tai_minus_utc(utc(2000, 1, 1)).unwrap(), 32.0);
        assert_eq!(tai_minus_utc(utc(1985, 7, 1)).unwrap(), 23.0);
        assert_eq!(tai_minus_utc(utc(1972, 1, 1)).unwrap(), 10.0);
    }

    #[test]
    fn delta_at_steps_exactly_on_boundary() {
        // 2017-01-01 0h で 36 → 37 へ跳ぶ。境界当日は新しい値。
        assert_eq!(tai_minus_utc(utc(2016, 12, 31)).unwrap(), 36.0);
        assert_eq!(tai_minus_utc(utc(2017, 1, 1)).unwrap(), 37.0);
    }

    #[test]
    fn before_1972_is_missing() {
        let pre = UtcInstant::from_gregorian(1971, 12, 31, 0, 0, 0.0).unwrap();
        assert_eq!(
            tai_minus_utc(pre).unwrap_err(),
            TimeError::MissingLeapSecondData
        );
    }

    #[test]
    fn utc_to_tt_offset_is_dat_plus_32_184() {
        // 2020: ΔAT=37 → TT−UTC = 69.184 s。
        let u = utc(2020, 1, 1);
        let tt = utc_to_tt(u).unwrap();
        let diff_s = tt.jd2().days_since(u.jd2()) * SECONDS_PER_DAY;
        assert!((diff_s - 69.184).abs() < 1e-6, "diff = {diff_s}");
    }

    #[test]
    fn j2000_utc_to_tt_offset_is_64_184() {
        // 2000-01-01 12:00 UTC: ΔAT=32 → TT−UTC = 64.184 s。
        let u = UtcInstant::from_gregorian(2000, 1, 1, 12, 0, 0.0).unwrap();
        let tt = utc_to_tt(u).unwrap();
        let diff_s = tt.jd2().days_since(u.jd2()) * SECONDS_PER_DAY;
        assert!((diff_s - 64.184).abs() < 1e-6, "diff = {diff_s}");
    }

    #[test]
    fn tai_tt_offset_is_exactly_32_184() {
        let tai = TaiInstant::from_jd2(JulianDate2::from_jd(2_460_000.0));
        let tt = tai_to_tt(tai);
        assert!((tt.jd2().days_since(tai.jd2()) * SECONDS_PER_DAY - 32.184).abs() < 1e-9);
    }

    #[test]
    fn utc_instant_gregorian_round_trip() {
        // UtcInstant::from_gregorian → to_gregorian の往復（委譲先 calendar とは別に本型を検証）。
        let u = UtcInstant::from_gregorian(2035, 9, 2, 1, 30, 15.5).unwrap();
        let (y, mo, d, h, mi, s) = u.to_gregorian();
        assert_eq!((y, mo, d, h, mi), (2035, 9, 2, 1, 30));
        assert!((s - 15.5).abs() < 1e-4, "s = {s}");
    }

    #[test]
    fn round_trip_utc_tt_utc() {
        let u = UtcInstant::from_gregorian(2035, 9, 2, 1, 30, 15.0).unwrap();
        let back = tt_to_utc(utc_to_tt(u).unwrap()).unwrap();
        assert!(back.jd2().days_since(u.jd2()).abs() * SECONDS_PER_DAY < 1e-6);
    }
}
