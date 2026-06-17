//! 同梱／外部パスからの時刻系データ束 [`umbra_core::TimeData`] 構築（ISSUE-042 S3）。
//!
//! 2 つの公開関数を提供する:
//! - [`bundled_time_data`]: 同梱 IERS EOP 14 C04（`eop.rs` の [`crate::eop::bundled_eop`]）と
//!   同梱閏秒（[`umbra_core::LeapSecondTable::bundled`]）を束ねた [`TimeData`]。`bundled-data`
//!   feature（既定 on）でのみ提供される。
//! - [`time_data_from_path`]: 外部ディレクトリの EOP データ（`<dir>/eopc04_14.bin` +
//!   `<dir>/metadata.txt`、`eop.rs` の packed 形式と byte-for-byte 同一）を読み、SHA-256
//!   checksum を metadata と照合した上で、同梱閏秒据え置きで [`TimeData`] を構築する。
//!
//! packed 形式は flat little-endian f64
//! `[n_records, then per record: mjd, ut1_minus_utc_s, x_pole_arcsec, y_pole_arcsec]`。
//! metadata.txt は `key = value` 行（name/version/source/license/valid_from/valid_to/checksum）。
//!
//! `time_data_from_path` の検証順序: ファイル読込（Io）→ metadata パース（MissingMetadataField）
//! → packed 構造デコード（MalformedPacked）→ SHA-256 照合（ChecksumMismatch）→
//! `IersEopData::from_records`（Eop）。構造デコードを checksum 照合より前に置くことで、
//! 不正な長さ／レコード数の bin が（checksum で短絡されず）構造検証経路で弾かれる。

use std::path::Path;

use umbra_core::{DataSetMetadata, EopRecord, IersEopData, LeapSecondTable, TimeData, TimeError};

/// `from_path` の外部 EOP データ読込エラー。
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum TimeDataError {
    /// EOP データファイル（`eopc04_14.bin` / `metadata.txt`）の読込に失敗した。
    #[error("failed to read EOP data file: {0}")]
    Io(#[from] std::io::Error),
    /// bin の SHA-256 が metadata の `checksum` と一致しない（完全性違反）。
    #[error("EOP checksum mismatch: metadata says {expected}, computed {actual}")]
    ChecksumMismatch {
        /// metadata.txt が宣言する checksum。
        expected: String,
        /// bin から計算した SHA-256。
        actual: String,
    },
    /// packed バイト列の構造が不正（長さが 8 の倍数でない・レコード数とバイト長の不整合）。
    #[error("malformed EOP packed data: {0}")]
    MalformedPacked(&'static str),
    /// metadata.txt に必須フィールドが欠落している。
    #[error("EOP metadata is missing field '{0}'")]
    MissingMetadataField(String),
    /// 供給された EOP テーブルが不正（`IersEopData::from_records` 検証失敗）。
    #[error(transparent)]
    Eop(#[from] TimeError),
}

/// 同梱 IERS EOP 14 C04 ＋ 同梱閏秒で [`TimeData`] を構築する（`bundled-data` feature, 既定 on）。
///
/// `= TimeData::new(LeapSecondTable::bundled(), bundled_eop())`。実行時ネットワークは行わない
/// （accuracy.md §5）。off 時は本関数がシンボルごと消え、[`time_data_from_path`] が必須となる。
#[cfg(feature = "bundled-data")]
pub fn bundled_time_data() -> TimeData {
    TimeData::new(LeapSecondTable::bundled(), crate::eop::bundled_eop())
}

/// 外部ディレクトリの EOP データ（`<dir>/eopc04_14.bin` + `<dir>/metadata.txt`、`eop.rs` の
/// packed 形式と byte-for-byte 同一）を読み、同梱閏秒据え置きで [`TimeData`] を構築する。
///
/// 閏秒は [`LeapSecondTable::bundled`] 据え置き（EOP のみ外部供給）。bin の SHA-256 を
/// metadata の `checksum` と照合する（完全性。accuracy.md §5）。エラーは [`TimeDataError`]。
pub fn time_data_from_path(dir: &Path) -> Result<TimeData, TimeDataError> {
    let bytes = std::fs::read(dir.join("eopc04_14.bin"))?;
    let text = std::fs::read_to_string(dir.join("metadata.txt"))?;

    let metadata = parse_metadata(&text)?;
    let records = decode_eop_packed(&bytes)?;

    let actual = sha256_hex(&bytes);
    if actual != metadata.checksum {
        return Err(TimeDataError::ChecksumMismatch {
            expected: metadata.checksum.clone(),
            actual,
        });
    }

    // 系列識別子は metadata の version（EOP 14 C04）。`eop.rs::bundled_eop` と同方針。
    let series_version = metadata.version.clone();
    let eop = IersEopData::from_records(records, series_version, metadata)?;
    Ok(TimeData::new(LeapSecondTable::bundled(), eop))
}

/// `metadata.txt`（`key = value` 行）→ [`DataSetMetadata`]。欠落フィールドは
/// [`TimeDataError::MissingMetadataField`]（`eop.rs::parse_metadata` の fallible 版）。
fn parse_metadata(text: &str) -> Result<DataSetMetadata, TimeDataError> {
    let field = |key: &str| -> Result<String, TimeDataError> {
        let prefix = format!("{key} = ");
        text.lines()
            .find_map(|line| line.strip_prefix(&prefix))
            .map(str::to_string)
            .ok_or_else(|| TimeDataError::MissingMetadataField(key.to_string()))
    };
    Ok(DataSetMetadata {
        name: field("name")?,
        version: field("version")?,
        source: field("source")?,
        license: field("license")?,
        valid_from: field("valid_from")?,
        valid_to: field("valid_to")?,
        checksum: field("checksum")?,
    })
}

/// packed（`[n, then 各レコード = mjd, ut1, x, y]`）を fallible にデコードする。
///
/// 長さが 8 の倍数でない／先頭レコード数とバイト長が不整合な場合は
/// [`TimeDataError::MalformedPacked`]（`eop.rs::decode_packed` の検証付き版・外部データ用）。
#[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
fn decode_eop_packed(bytes: &[u8]) -> Result<Vec<EopRecord>, TimeDataError> {
    if bytes.len() % 8 != 0 {
        return Err(TimeDataError::MalformedPacked(
            "byte length is not a multiple of 8",
        ));
    }
    let n_f64 = bytes.len() / 8;
    if n_f64 == 0 {
        return Err(TimeDataError::MalformedPacked("empty (no record count)"));
    }
    // 先頭 f64 = レコード数。NaN/Inf/負値は `as usize` の飽和変換（Rust 1.45+, パニックなし）で
    // 0 または usize::MAX となり、いずれも後続の `checked_mul` オーバーフロー or `n_f64 != needed`
    // で MalformedPacked として弾かれる（明示の有限性ガードは冗長ゆえ置かない）。
    let count = packed_f64(bytes, 0).round() as usize;
    // 必要な f64 数 = 1（count）+ count*4（各レコード）。オーバーフローも構造不正として弾く。
    let needed = count
        .checked_mul(4)
        .and_then(|r| r.checked_add(1))
        .ok_or(TimeDataError::MalformedPacked("record count overflows"))?;
    if n_f64 != needed {
        return Err(TimeDataError::MalformedPacked(
            "record count does not match byte length",
        ));
    }
    let mut records = Vec::with_capacity(count);
    for i in 0..count {
        let base = 1 + i * 4;
        records.push(EopRecord::new(
            packed_f64(bytes, base).round() as i32,
            packed_f64(bytes, base + 1),
            packed_f64(bytes, base + 2),
            packed_f64(bytes, base + 3),
        ));
    }
    Ok(records)
}

/// packed の `index` 番目の f64（little-endian）。呼び出し側で範囲は保証済み。
fn packed_f64(bytes: &[u8], index: usize) -> f64 {
    let start = index * 8;
    let octet: [u8; 8] = bytes[start..start + 8]
        .try_into()
        .expect("index within bounds (length validated as a multiple of 8)");
    f64::from_le_bytes(octet)
}

/// バイト列の SHA-256 を 16 進小文字で返す。`xtask::checksum::sha256_hex` および metadata.txt の
/// `checksum` 生成と **byte-for-byte 同一**であること（契約・変更時は両実装を同期）。crate 越しの
/// ビルド依存を避けるため再実装だが、出力が違うと from_path が正データを ChecksumMismatch で誤拒否する。
fn sha256_hex(bytes: &[u8]) -> String {
    use sha2::{Digest, Sha256};
    let digest = Sha256::digest(bytes);
    let mut hex = String::with_capacity(digest.len() * 2);
    for byte in digest {
        use std::fmt::Write as _;
        let _ = write!(hex, "{byte:02x}");
    }
    hex
}

// ====================================================================
// テストモジュール。
//
// bundled 同等性テストが `bundled-data` feature の `bundled_eop`/`bundled_time_data` に依存する
// ため、テストモジュール自体を feature ゲートする（既定 features では on）。`time_data_from_path`
// 単体は feature 非依存だが、簡便のため同 mod 内にまとめる（ISSUE-042 S3 設計指示）。
// ====================================================================
#[cfg(all(test, feature = "bundled-data"))]
mod tests {
    // 検証値は IERS EOP 14 C04 のコミット済みデータの逐語転記（provenance 保持）。
    // 余剰桁は最近接 f64 へ丸められ値は不変のため、過剰精度リントを許可する。
    #![allow(clippy::excessive_precision)]

    use std::path::{Path, PathBuf};

    use umbra_core::{LeapSecondTable, Radians, TimeScales, UtcInstant};

    // 実装予定の公開 IF（crate 内テストのため `crate::time_data` 経由で参照）。再エクスポートは
    // 実装工程で `lib.rs` に追加される。
    use crate::time_data::{bundled_time_data, time_data_from_path, TimeDataError};

    // ---- 定数・補助 -------------------------------------------------------

    /// arcsec → radian 係数（1" = π/648000 rad ≈ 4.84813681109536e-6）。独立計算。
    const ARCSEC_TO_RAD: f64 = std::f64::consts::PI / 648_000.0;
    const SECONDS_PER_DAY: f64 = 86_400.0;

    /// 同梱 EOP 2020-01-01 0h の実測オラクル（eop.rs と同一値）。
    const UT1_20200101: f64 = -0.1771222;
    const XP_20200101: f64 = 0.076609;
    const YP_20200101: f64 = 0.282358;
    /// 同梱 EOP metadata の SHA-256 checksum（コミット済み bin・eop.rs と同一）。
    const EOP_CHECKSUM: &str = "4d12cc1d3dcae39a8db58690a0d6e669def893903ba74c23995ced1e5e354ac3";

    /// グレゴリオ暦（UTC）の瞬時を構築する補助。
    fn utc(y: i32, mo: u8, d: u8, h: u8, mi: u8, s: f64) -> UtcInstant {
        UtcInstant::from_gregorian(y, mo, d, h, mi, s).expect("valid calendar date")
    }

    /// 同梱 EOP 生成元ディレクトリ `generated/eop`（クレートからの相対）。
    fn generated_eop_dir() -> PathBuf {
        PathBuf::from(format!(
            "{}/../../generated/eop",
            env!("CARGO_MANIFEST_DIR")
        ))
    }

    /// テスト専用の一意な一時ディレクトリを作り、与えた bin/metadata を書き込む。
    /// `bin` が `None` の場合は eopc04_14.bin を書かない（欠落させる）。
    /// `metadata` が `None` の場合は metadata.txt を書かない。
    fn make_fixture_dir(name: &str, bin: Option<&[u8]>, metadata: Option<&str>) -> PathBuf {
        let dir = std::env::temp_dir().join(format!("umbra_s3_{name}"));
        // 既存の残骸を掃除してから作り直す（前回失敗時の残りで誤判定しないため）。
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).expect("create temp fixture dir");
        if let Some(bytes) = bin {
            std::fs::write(dir.join("eopc04_14.bin"), bytes).expect("write bin");
        }
        if let Some(text) = metadata {
            std::fs::write(dir.join("metadata.txt"), text).expect("write metadata");
        }
        dir
    }

    /// little-endian f64 列を packed バイト列へ。
    fn pack_f64s(values: &[f64]) -> Vec<u8> {
        let mut out = Vec::with_capacity(values.len() * 8);
        for v in values {
            out.extend_from_slice(&v.to_le_bytes());
        }
        out
    }

    /// n=1 の最小有効 packed（mjd/ut1/xp/yp 1 レコード = 先頭 count を含め 5 f64 = 40 バイト）。
    /// 構造としては正しい（長さは 8 の倍数・n と整合）が、checksum は別途指定する。
    fn minimal_valid_packed() -> Vec<u8> {
        pack_f64s(&[1.0, 58849.0, UT1_20200101, XP_20200101, YP_20200101])
    }

    /// 7 フィールド完備の metadata.txt 文面（checksum を引数で差し替え可能）。
    fn metadata_text(checksum: &str) -> String {
        format!(
            "name = eop-iers-c04\n\
             version = EOP 14 C04\n\
             source = test\n\
             license = test\n\
             valid_from = 2020-01-01\n\
             valid_to = 2020-01-01\n\
             checksum = {checksum}\n"
        )
    }

    // ==================================================================
    // 1. bundled_time_data: valid_range = [1972-01-01, 2026-01-05]。
    //    下端 = max(閏秒起点 1972, EOP 起点 1962) = 1972。上端 = EOP 終了 2026-01-05。
    //    変異: 下端に EOP 起点(1962)/min を使う、上端に閏秒・別日を使う、start/end 取り違えを殺す。
    // ==================================================================
    #[test]
    fn bundled_time_data_valid_range_is_1972_to_2026() {
        let data = bundled_time_data();
        let range = data.valid_range();
        let want_start = utc(1972, 1, 1, 0, 0, 0.0);
        let want_end = utc(2026, 1, 5, 0, 0, 0.0);
        assert!(
            (range.start.jd2().jd() - want_start.jd2().jd()).abs() < 1e-9,
            "valid_range.start jd = {}, want {} (1972-01-01, leap binding)",
            range.start.jd2().jd(),
            want_start.jd2().jd()
        );
        assert!(
            (range.end.jd2().jd() - want_end.jd2().jd()).abs() < 1e-9,
            "valid_range.end jd = {}, want {} (2026-01-05, EOP end)",
            range.end.jd2().jd(),
            want_end.jd2().jd()
        );
    }

    // ==================================================================
    // 2. bundled_time_data: metadata は 2 件（閏秒, EOP の順）で各 provenance 完全。
    //    閏秒 metadata は LeapSecondTable::bundled().metadata() と一致、
    //    EOP metadata は version="EOP 14 C04"・checksum=既知値。
    //    変異: 件数取り違え、順序入替（EOP→閏秒）、別メタデータ混入、provenance 欠落を殺す。
    // ==================================================================
    #[test]
    fn bundled_time_data_metadata_is_leap_then_eop_complete() {
        let data = bundled_time_data();
        let md = data.metadata();
        assert_eq!(md.len(), 2, "metadata must contain exactly 2 datasets");
        // [0] = 閏秒。
        assert_eq!(
            md[0],
            LeapSecondTable::bundled().metadata(),
            "metadata[0] must be the bundled leap-second dataset"
        );
        // [1] = EOP（同梱 EOP の provenance）。
        assert_eq!(md[1].version, "EOP 14 C04", "metadata[1].version (EOP)");
        assert_eq!(md[1].checksum, EOP_CHECKSUM, "metadata[1].checksum (EOP)");
        assert!(
            md[0].has_complete_provenance(),
            "leap-second metadata must have complete provenance"
        );
        assert!(
            md[1].has_complete_provenance(),
            "EOP metadata must have complete provenance"
        );
    }

    // ==================================================================
    // 3. bundled_time_data → TimeScales: utc_to_tt(2020-01-01) の TT−UTC = 69.184 s（exact）。
    //    ΔAT(2020)=37 + 32.184 = 69.184。閏秒テーブルが正しく束ねられていることを固定。
    //    変異: 閏秒テーブルを束ねない・別テーブル混入、32.184 取り違えを殺す。
    // ==================================================================
    #[test]
    fn bundled_time_data_scales_utc_to_tt_2020_is_69_184() {
        let scales = TimeScales::new(bundled_time_data());
        let u = utc(2020, 1, 1, 0, 0, 0.0);
        let tt = scales
            .utc_to_tt(u)
            .expect("2020 is in bundled leap-second table");
        let diff_s = tt.jd2().days_since(u.jd2()) * SECONDS_PER_DAY;
        assert!(
            (diff_s - 69.184).abs() < 1e-6,
            "TT-UTC(2020) = {diff_s} s, want 69.184"
        );
    }

    // ==================================================================
    // 4. bundled_time_data → TimeScales: utc_to_ut1(2020-01-01 0h) は
    //    UTC + (UT1−UTC=-0.1771222)/86400 日（exact）。EOP が正しく束ねられていることを固定。
    //    変異: EOP を束ねない・別 EOP 混入、UT1−UTC 加算の符号/除数取り違えを殺す。
    // ==================================================================
    #[test]
    fn bundled_time_data_scales_utc_to_ut1_2020_adds_known_offset() {
        let scales = TimeScales::new(bundled_time_data());
        let u = utc(2020, 1, 1, 0, 0, 0.0);
        let ut1 = scales
            .utc_to_ut1(u)
            .expect("2020-01-01 is in bundled EOP coverage");
        let diff_s = ut1.jd2().days_since(u.jd2()) * SECONDS_PER_DAY;
        assert!(
            (diff_s - UT1_20200101).abs() < 1e-6,
            "UT1-UTC(2020-01-01) = {diff_s} s, want {UT1_20200101}"
        );
    }

    // ==================================================================
    // 5. bundled_time_data → TimeScales: polar_motion(2020-01-01 0h) = (xp, yp) arcsec→rad（exact）。
    //    変異: arcsec→rad 係数の脱落/取り違え、xp/yp 入替、EOP 不参照を殺す。
    // ==================================================================
    #[test]
    fn bundled_time_data_scales_polar_motion_2020_in_radians() {
        let scales = TimeScales::new(bundled_time_data());
        let (Radians(xp), Radians(yp)) = scales
            .polar_motion(utc(2020, 1, 1, 0, 0, 0.0))
            .expect("2020-01-01 is in bundled EOP coverage");
        let want_xp = XP_20200101 * ARCSEC_TO_RAD;
        let want_yp = YP_20200101 * ARCSEC_TO_RAD;
        assert!((xp - want_xp).abs() < 1e-12, "xp = {xp}, want {want_xp}");
        assert!((yp - want_yp).abs() < 1e-12, "yp = {yp}, want {want_yp}");
    }

    // ==================================================================
    // 6. time_data_from_path(generated/eop): Ok を返す（正常系の基本）。
    //    変異: 正常データを誤って Err にする、metadata/bin 読込の未配線を殺す。
    // ==================================================================
    #[test]
    fn time_data_from_path_reads_generated_eop_ok() {
        let dir = generated_eop_dir();
        let data =
            time_data_from_path(&dir).expect("generated/eop directory must build a valid TimeData");
        // 構築できれば metadata は 2 件（閏秒, EOP）。
        assert_eq!(data.metadata().len(), 2, "metadata must contain 2 datasets");
    }

    // ==================================================================
    // 7. bundled 同等性（ISSUE-042 §73）: time_data_from_path(generated/eop) は
    //    bundled_time_data() と同じ EOP 値・coverage・metadata checksum になる（同一データ源）。
    //    変異: from_path が別の bin/metadata を読む・列順を取り違える・閏秒を据え置かない、を殺す。
    // ==================================================================
    #[test]
    fn time_data_from_path_generated_eop_matches_bundled() {
        let from_path = time_data_from_path(&generated_eop_dir())
            .expect("generated/eop builds a valid TimeData");
        let bundled = bundled_time_data();

        // (a) valid_range（coverage 由来の上端・閏秒 binding の下端）が一致。
        let rp = from_path.valid_range();
        let rb = bundled.valid_range();
        assert!(
            (rp.start.jd2().jd() - rb.start.jd2().jd()).abs() < 1e-9,
            "valid_range.start mismatch: from_path {}, bundled {}",
            rp.start.jd2().jd(),
            rb.start.jd2().jd()
        );
        assert!(
            (rp.end.jd2().jd() - rb.end.jd2().jd()).abs() < 1e-9,
            "valid_range.end mismatch: from_path {}, bundled {}",
            rp.end.jd2().jd(),
            rb.end.jd2().jd()
        );

        // (b) EOP metadata checksum（[1]=EOP）が一致＝同一 bin。
        assert_eq!(
            from_path.metadata()[1].checksum,
            bundled.metadata()[1].checksum,
            "EOP metadata checksum must match the bundled data"
        );
        assert_eq!(
            from_path.metadata()[1].checksum,
            EOP_CHECKSUM,
            "EOP metadata checksum must equal the committed value"
        );

        // (c) 2020-01-01 0h の UT1−UTC が両者で実測オラクルに一致（同一 EOP 値）。
        let p_scales = TimeScales::new(from_path);
        let b_scales = TimeScales::new(bundled);
        let u = utc(2020, 1, 1, 0, 0, 0.0);
        let p = p_scales.utc_to_ut1(u).expect("from_path covers 2020-01-01");
        let b = b_scales.utc_to_ut1(u).expect("bundled covers 2020-01-01");
        let p_s = p.jd2().days_since(u.jd2()) * SECONDS_PER_DAY;
        let b_s = b.jd2().days_since(u.jd2()) * SECONDS_PER_DAY;
        assert!(
            (p_s - UT1_20200101).abs() < 1e-6 && (b_s - UT1_20200101).abs() < 1e-6,
            "UT1-UTC(2020) from_path={p_s}, bundled={b_s}, want {UT1_20200101}"
        );
    }

    // ==================================================================
    // 8. 異常: 存在しないパス → TimeDataError::Io。
    //    変異: 読込失敗を別 variant にマップする、Result を握り潰す、を殺す。
    // ==================================================================
    #[test]
    fn time_data_from_path_nonexistent_is_io_error() {
        let err = time_data_from_path(Path::new("/nonexistent/umbra/xyz_s3"))
            .expect_err("nonexistent path must error");
        assert!(
            matches!(err, TimeDataError::Io(_)),
            "expected Io, got {err:?}"
        );
    }

    // ==================================================================
    // 9. 異常: checksum 不一致 → TimeDataError::ChecksumMismatch。
    //    構造的に有効な bin（n=1 の 40 バイト）＋ metadata.checksum をわざと別値にする。
    //    変異: checksum 照合の脱落、別 variant への取り違え、不一致でも Ok を返す、を殺す。
    // ==================================================================
    #[test]
    fn time_data_from_path_checksum_mismatch() {
        let bin = minimal_valid_packed();
        // bin の真の SHA-256 とは一致しないことが明らかな偽 checksum。
        let bogus = "deadbeefdeadbeefdeadbeefdeadbeefdeadbeefdeadbeefdeadbeefdeadbeef0";
        let dir = make_fixture_dir("checksum_mismatch", Some(&bin), Some(&metadata_text(bogus)));
        let err = time_data_from_path(&dir).expect_err("checksum mismatch must error");
        assert!(
            matches!(err, TimeDataError::ChecksumMismatch { .. }),
            "expected ChecksumMismatch, got {err:?}"
        );
        let _ = std::fs::remove_dir_all(&dir);
    }

    // ==================================================================
    // 10. 異常: bin 長が 8 の倍数でない → Err（構造不正）。
    //     7 バイトの bin。checksum 一致を作るのが難しいため、variant は固定せず
    //     「Err であること」のみを確認する（設計指示: MalformedPacked でも ChecksumMismatch でも可）。
    //     変異: 長さ検証の脱落（パニックや誤った Ok）を殺す。
    // ==================================================================
    #[test]
    fn time_data_from_path_non_multiple_of_8_is_error() {
        let bin = vec![0u8; 7]; // 8 の倍数でない。
        let dir = make_fixture_dir(
            "malformed_len7",
            Some(&bin),
            Some(&metadata_text(EOP_CHECKSUM)),
        );
        let result = time_data_from_path(&dir);
        assert!(
            result.is_err(),
            "7-byte bin must be Err (structural), got {result:?}"
        );
        let _ = std::fs::remove_dir_all(&dir);
    }

    // ==================================================================
    // 11. 異常: n（先頭 count）とバイト数が不整合 → Err（構造不正）。
    //     count=100 と宣言しつつ実体は 1 レコード分しかない（40 バイト = 8 の倍数だが
    //     必要な 1+100*4 f64 に満たない）。variant は固定せず「Err であること」のみ確認。
    //     変異: count とバイト長の整合検証の脱落（範囲外参照パニックや誤った Ok）を殺す。
    // ==================================================================
    #[test]
    fn time_data_from_path_count_byte_mismatch_is_error() {
        // count=100 だが続くレコードは 1 件分しかない。
        let bin = pack_f64s(&[100.0, 58849.0, UT1_20200101, XP_20200101, YP_20200101]);
        let dir = make_fixture_dir(
            "malformed_count",
            Some(&bin),
            Some(&metadata_text(EOP_CHECKSUM)),
        );
        let result = time_data_from_path(&dir);
        assert!(
            result.is_err(),
            "count/byte mismatch must be Err (structural), got {result:?}"
        );
        let _ = std::fs::remove_dir_all(&dir);
    }

    // ==================================================================
    // 12. 異常: metadata.txt のフィールド欠落 → TimeDataError::MissingMetadataField。
    //     7 フィールドのうち checksum 行を削る。構造的に有効な bin を与え、欠落で弾かれることを固定。
    //     変異: フィールド欠落検出の脱落（panic や別 variant）、欠落でも Ok を返す、を殺す。
    // ==================================================================
    #[test]
    fn time_data_from_path_missing_metadata_field() {
        let bin = minimal_valid_packed();
        // checksum 行を削った metadata（残り 6 フィールド）。
        let md = "name = eop-iers-c04\n\
                  version = EOP 14 C04\n\
                  source = test\n\
                  license = test\n\
                  valid_from = 2020-01-01\n\
                  valid_to = 2020-01-01\n";
        let dir = make_fixture_dir("missing_field", Some(&bin), Some(md));
        let err = time_data_from_path(&dir).expect_err("missing field must error");
        assert!(
            matches!(err, TimeDataError::MissingMetadataField(_)),
            "expected MissingMetadataField, got {err:?}"
        );
        let _ = std::fs::remove_dir_all(&dir);
    }

    // ==================================================================
    // 13. 異常: bin は存在するが metadata.txt が無い → TimeDataError::Io。
    //     片方のファイル欠落も読込失敗として Io にマップされることを固定。
    //     変異: metadata 欠落を別 variant にマップする・無視するを殺す。
    // ==================================================================
    #[test]
    fn time_data_from_path_missing_metadata_file_is_io_error() {
        let bin = minimal_valid_packed();
        let dir = make_fixture_dir("missing_metadata_file", Some(&bin), None);
        let err = time_data_from_path(&dir).expect_err("missing metadata.txt must error");
        assert!(
            matches!(err, TimeDataError::Io(_)),
            "expected Io (missing metadata.txt), got {err:?}"
        );
        let _ = std::fs::remove_dir_all(&dir);
    }
}
