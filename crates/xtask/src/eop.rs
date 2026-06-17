//! IERS EOP 14 C04 地球姿勢データの取り込み（ISSUE-007-EOP, part P2a）。
//!
//! 一次配布の C04 テキスト（`data/coefficient-source/eop/`・PROVENANCE.md）をパースし、
//! flat little-endian f64 の packed 形式へ直列化する（消費側 = umbra-ephemeris `bundled-data`
//! と byte-for-byte 契約・part P2b）。取り込む列は **MJD（整数, 0h UTC）/ x（arcsec）/ y（arcsec）/
//! UT1−UTC（秒）**（LOD・dX/dY・各誤差は v0.1 未使用, PROVENANCE.md）。
//!
//! packed レイアウト: `[ n_records, then 各レコード = [mjd, ut1_minus_utc_s, x_pole_arcsec,
//! y_pole_arcsec] ]`（全て LE f64, mjd は整数を無損失格納）。

use crate::checksum::sha256_hex;
use crate::error::XtaskError;
use crate::packed::{pack_f64_le, unpack_f64_le};
use std::fs;
use umbra_core::{DataSetMetadata, EopRecord, JulianDate2};

/// 原データ配置（リポジトリルートから実行する前提, accuracy.md §5）。
const SOURCE_DIR: &str = "data/coefficient-source/eop";
/// 一次原データファイル名（PROVENANCE.md）。
const SOURCE_NAME: &str = "eopc04_14_IAU2000A_1962-now.txt";
/// 生成物の出力先。
const GENERATED_DIR: &str = "generated/eop";
/// packed ファイル名。
const PACKED_NAME: &str = "eopc04_14.bin";
/// MJD = JD − 2400000.5。
const MJD_JD_OFFSET: f64 = 2_400_000.5;
/// packed ヘッダ件数の健全性上限（破損データ防御。C04 日次は ~2.3 万件で十分余裕）。
const MAX_EOP_RECORDS: f64 = 10_000_000.0;

/// IERS EOP 14 C04 テキストを日次レコード列へパースする（ファイル順＝昇順 MJD を保持）。
///
/// データ行は「先頭 4 トークンが year/month/day/MJD の整数（暦として妥当）」で識別し、ヘッダ・
/// 空行・FORMAT 行・`####` 行は読み飛ばす。x=col5 / y=col6 / UT1−UTC=col7（arcsec / arcsec / 秒）。
/// データ行が 0 件、または数値列が壊れている場合は [`XtaskError::MalformedSource`]。
pub fn parse_eop_c04(text: &str) -> Result<Vec<EopRecord>, XtaskError> {
    let mut records = Vec::new();
    for line in text.lines() {
        let tokens: Vec<&str> = line.split_whitespace().collect();
        if tokens.len() < 7 {
            continue;
        }
        // 先頭 4 整数（year, month, day, MJD）でデータ行を識別。
        let (Ok(year), Ok(month), Ok(day), Ok(mjd)) = (
            tokens[0].parse::<i32>(),
            tokens[1].parse::<i32>(),
            tokens[2].parse::<i32>(),
            tokens[3].parse::<i32>(),
        ) else {
            continue;
        };
        if !(1800..=2200).contains(&year) || !(1..=12).contains(&month) || !(1..=31).contains(&day)
        {
            continue;
        }
        let x_pole_arcsec = parse_value(tokens[4], "x_pole", line)?;
        let y_pole_arcsec = parse_value(tokens[5], "y_pole", line)?;
        let ut1_minus_utc_s = parse_value(tokens[6], "UT1-UTC", line)?;
        records.push(EopRecord::new(
            mjd,
            ut1_minus_utc_s,
            x_pole_arcsec,
            y_pole_arcsec,
        ));
    }
    if records.is_empty() {
        return Err(XtaskError::MalformedSource(
            "no EOP C04 data rows found (expected 'year month day MJD x y UT1-UTC ...')"
                .to_string(),
        ));
    }
    Ok(records)
}

/// データ行の数値列をパースする（壊れていれば [`XtaskError::MalformedSource`]）。
fn parse_value(token: &str, column: &str, line: &str) -> Result<f64, XtaskError> {
    token.parse::<f64>().map_err(|_| {
        XtaskError::MalformedSource(format!("invalid {column} value '{token}' in line: {line}"))
    })
}

/// レコード列を packed（flat LE f64）へ直列化する。
pub fn pack_records(records: &[EopRecord]) -> Vec<u8> {
    let count = u32::try_from(records.len()).expect("EOP record count fits in u32");
    let mut values: Vec<f64> = Vec::with_capacity(1 + records.len() * 4);
    values.push(f64::from(count));
    for r in records {
        values.push(f64::from(r.mjd));
        values.push(r.ut1_minus_utc_s);
        values.push(r.x_pole_arcsec);
        values.push(r.y_pole_arcsec);
    }
    pack_f64_le(&values)
}

/// packed を レコード列へ復号する（[`pack_records`] の逆。消費側と同一契約）。
///
/// バイト長が 8 の倍数でない、ヘッダ件数が残り長と矛盾する場合は [`XtaskError::MalformedPacked`]。
#[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
pub fn unpack_records(bytes: &[u8]) -> Result<Vec<EopRecord>, XtaskError> {
    let values = unpack_f64_le(bytes)?;
    let header = *values.first().ok_or_else(|| {
        XtaskError::MalformedPacked("packed EOP is empty (missing record-count header)".to_string())
    })?;
    // ヘッダ件数は有限・非負・整数・上限内（破損データの巨大/負/小数 count を弾く, nutation 同方針）。
    if !header.is_finite() || header < 0.0 || header.fract() != 0.0 || header > MAX_EOP_RECORDS {
        return Err(XtaskError::MalformedPacked(format!(
            "invalid EOP record-count header {header}"
        )));
    }
    let count = header as usize; // 検証済み: 0..=MAX_EOP_RECORDS の整数。
    if values.len() != 1 + count * 4 {
        return Err(XtaskError::MalformedPacked(format!(
            "packed EOP length {} inconsistent with header count {count} (expected {})",
            values.len(),
            1 + count * 4
        )));
    }
    let mut records = Vec::with_capacity(count);
    for chunk in values[1..].chunks_exact(4) {
        let mjd = chunk[0].round() as i32;
        records.push(EopRecord::new(mjd, chunk[1], chunk[2], chunk[3]));
    }
    Ok(records)
}

/// 原データ文字列から packed バイト列と `DataSetMetadata` を決定的に構築する。
fn build_artifact(source: &str) -> Result<(Vec<u8>, DataSetMetadata), XtaskError> {
    let records = parse_eop_c04(source)?;
    let bytes = pack_records(&records);
    let checksum = sha256_hex(&bytes);
    let first = records.first().expect("non-empty after parse");
    let last = records.last().expect("non-empty after parse");
    let metadata = DataSetMetadata {
        name: "eop-iers-c04".to_string(),
        version: "EOP 14 C04".to_string(),
        source:
            "IERS Earth Orientation Centre, EOP (IERS) 14 C04 time series (datacenter.iers.org)"
                .to_string(),
        license: "IERS published scientific data (attribution); not a GPL derivative".to_string(),
        valid_from: mjd_to_iso(first.mjd),
        valid_to: mjd_to_iso(last.mjd),
        checksum,
    };
    Ok((bytes, metadata))
}

/// 整数 MJD（0h UTC）→ `YYYY-MM-DD`。
fn mjd_to_iso(mjd: i32) -> String {
    let (year, month, day, _, _, _) =
        umbra_core::jd2_to_gregorian(JulianDate2::from_jd(f64::from(mjd) + MJD_JD_OFFSET));
    format!("{year:04}-{month:02}-{day:02}")
}

/// `SOURCE_DIR` 配下のファイルを読む（I/O エラーは [`XtaskError::Io`]）。
fn read_source(file_name: &str) -> Result<String, XtaskError> {
    let path = format!("{SOURCE_DIR}/{file_name}");
    fs::read_to_string(&path).map_err(|source| XtaskError::Io { path, source })
}

/// `metadata.txt` の人間可読レンダリング（nutation と同契約）。
fn render_metadata(metadata: &DataSetMetadata) -> String {
    format!(
        "name = {}\nversion = {}\nsource = {}\nlicense = {}\nvalid_from = {}\nvalid_to = {}\nchecksum = {}\n",
        metadata.name,
        metadata.version,
        metadata.source,
        metadata.license,
        metadata.valid_from,
        metadata.valid_to,
        metadata.checksum,
    )
}

/// 一次原データを読み、packed・metadata・NOTICE を `GENERATED_DIR` へ決定的に書き出す。
/// `cargo xtask generate-coefficients --dataset eop-c04`（リポジトリルートで実行）。
pub fn generate_to_disk() -> Result<DataSetMetadata, XtaskError> {
    let source = read_source(SOURCE_NAME)?;
    let (bytes, metadata) = build_artifact(&source)?;

    fs::create_dir_all(GENERATED_DIR).map_err(|source| XtaskError::Io {
        path: GENERATED_DIR.to_string(),
        source,
    })?;
    let write = |name: &str, content: &[u8]| -> Result<(), XtaskError> {
        let path = format!("{GENERATED_DIR}/{name}");
        fs::write(&path, content).map_err(|source| XtaskError::Io { path, source })
    };
    write(PACKED_NAME, &bytes)?;
    write("metadata.txt", render_metadata(&metadata).as_bytes())?;
    write(
        "NOTICE.md",
        format!(
            "# Generated IERS EOP 14 C04 data\n\n\
             Generated by `cargo xtask generate-coefficients --dataset eop-c04` from\n\
             `{SOURCE_DIR}/` (see PROVENANCE.md). Do not edit by hand; regenerate and verify with\n\
             `cargo xtask verify-generated --dataset eop-c04`.\n\n\
             {}\n",
            render_metadata(&metadata)
        )
        .as_bytes(),
    )?;
    Ok(metadata)
}

/// コミット済み packed が一次原データから決定的に再生成できることを検証する。
/// 差分があれば [`XtaskError::ChecksumMismatch`]（CI fail）。
pub fn verify_against_disk() -> Result<(), XtaskError> {
    let source = read_source(SOURCE_NAME)?;
    let (regenerated, _) = build_artifact(&source)?;
    let committed_path = format!("{GENERATED_DIR}/{PACKED_NAME}");
    let committed = fs::read(&committed_path).map_err(|source| XtaskError::Io {
        path: committed_path,
        source,
    })?;
    crate::compare_checksum("eop-c04", &sha256_hex(&committed), &regenerated)
}

#[cfg(test)]
mod tests {
    use super::*;
    use umbra_core::EopRecord;

    // ------------------------------------------------------------------
    // オラクル: IERS EOP 14 C04 の最小有効抜粋（ヘッダ + 3 データ行）。
    // 出典: data/coefficient-source/eop/PROVENANCE.md（IERS Earth Orientation
    // Centre, EOP (IERS) 14 C04 time series）。
    // データ行列: year month day MJD x" y" UT1-UTC[s] LOD ...
    // 取り込む列: MJD=col4 / x=col5 / y=col6 / UT1-UTC=col7（LOD 以降は無視）。
    // ------------------------------------------------------------------
    const EXCERPT: &str = "\
               EARTH ORIENTATION PARAMETER (EOP) PRODUCT CENTER CENTER (PARIS OBSERVATORY)
                      INTERNATIONAL EARTH ROTATION AND REFERENCE SYSTEMS SERVICE
                                    EOP (IERS) 14 C04 TIME SERIES
             FORMAT(3(I4),I7,2(F11.6),2(F12.7),2(F11.6),2(F11.6),2(F11.7),2(F12.6))
##################################################################################

      Date      MJD      x          y        UT1-UTC
                         \"          \"           s
     (0h UTC)

1962   1   1  37665  -0.012700   0.213000   0.0326338   0.0017230   0.000000   0.000000   0.030000   0.030000  0.0020000  0.0014000    0.004774    0.002000
1962   1   2  37666  -0.015900   0.214100   0.0320547   0.0016690   0.000000   0.000000   0.030000   0.030000  0.0020000  0.0014000    0.004774    0.002000
2020   1   1  58849   0.076609   0.282358  -0.1771222   0.0004455   0.000348   0.000003   0.000057   0.000039  0.0000123  0.0000126    0.000044    0.000057
";

    /// コミット済み実原データ（cwd 非依存に manifest 相対 include_str! で読む）。
    const REAL_C04: &str =
        include_str!("../../../data/coefficient-source/eop/eopc04_14_IAU2000A_1962-now.txt");

    /// コミット済み生成物（cwd 非依存に manifest 相対 include_bytes! で読む）。
    const COMMITTED_BIN: &[u8] = include_bytes!("../../../generated/eop/eopc04_14.bin");

    /// 1 レコードの全フィールドを厳密一致で確認するヘルパ。
    /// 振幅は表の 10 進が f64 へ一意変換されるため厳密 `==`（リテラル比較）でよい。
    #[track_caller]
    fn assert_record(rec: &EopRecord, mjd: i32, ut1: f64, x: f64, y: f64) {
        assert_eq!(rec.mjd, mjd, "mjd mismatch");
        assert_eq!(rec.ut1_minus_utc_s, ut1, "ut1_minus_utc_s mismatch (col7)");
        assert_eq!(rec.x_pole_arcsec, x, "x_pole_arcsec mismatch (col5)");
        assert_eq!(rec.y_pole_arcsec, y, "y_pole_arcsec mismatch (col6)");
    }

    // ------------------------------------------------------------------
    // (1) parse: レコード件数 & 各フィールド値。
    //     ヘッダ・空行・FORMAT 行・"#####" 行などの非データ行が skip され、
    //     ちょうど 3 レコードになること。各レコードの 4 フィールドを厳密一致で確認。
    // ------------------------------------------------------------------
    #[test]
    fn parse_excerpt_record_count_and_values() {
        let recs = parse_eop_c04(EXCERPT).expect("minimal valid C04 excerpt parses");
        assert_eq!(
            recs.len(),
            3,
            "exactly 3 data rows; header/blank/non-data lines skipped"
        );
        assert_record(&recs[0], 37665, 0.0326338, -0.012700, 0.213000);
        assert_record(&recs[1], 37666, 0.0320547, -0.015900, 0.214100);
        assert_record(&recs[2], 58849, -0.1771222, 0.076609, 0.282358);
    }

    // ------------------------------------------------------------------
    // (2) parse: 昇順 mjd がファイル順どおりに保たれる。
    //     抜粋は 37665 < 37666 < 58849 の昇順 → そのままの順で返る。
    // ------------------------------------------------------------------
    #[test]
    fn parse_excerpt_preserves_ascending_mjd_order() {
        let recs = parse_eop_c04(EXCERPT).expect("excerpt parses");
        let mjds: Vec<i32> = recs.iter().map(|r| r.mjd).collect();
        assert_eq!(mjds, vec![37665, 37666, 58849], "file order preserved");
        assert!(
            mjds.windows(2).all(|w| w[0] < w[1]),
            "mjd strictly ascending"
        );
    }

    // ------------------------------------------------------------------
    // (3) parse: コミット済み実ファイル全体。
    //     (a) Ok / (b) 数千件以上 / (c) 先頭 mjd==37665, ut1==0.0326338 /
    //     (d) 全体で mjd が厳密昇順 / (e) 内部既知行 mjd==58849 の x/y/ut1 一致。
    // ------------------------------------------------------------------
    #[test]
    fn parse_real_committed_file() {
        let recs = parse_eop_c04(REAL_C04).expect("real committed C04 file parses (Ok)");

        // (b) 数千件以上（系列は 1962–現在の日次 → 2 万件超）。
        assert!(
            recs.len() > 20_000,
            "expected > 20000 daily records, got {}",
            recs.len()
        );

        // (c) 先頭レコード = 1962-01-01。
        assert_eq!(recs[0].mjd, 37665, "first record mjd == 37665");
        assert_eq!(
            recs[0].ut1_minus_utc_s, 0.0326338,
            "first record ut1 == 0.0326338"
        );

        // (d) 全体で mjd が厳密昇順（隙間検査は不要、厳密昇順のみ）。
        assert!(
            recs.windows(2).all(|w| w[0].mjd < w[1].mjd),
            "mjd strictly ascending across the whole file"
        );

        // (e) 内部既知行 mjd==58849（2020-01-01）の全数値事実。
        let interior = recs
            .iter()
            .find(|r| r.mjd == 58849)
            .expect("interior record mjd == 58849 exists");
        assert_eq!(interior.ut1_minus_utc_s, -0.1771222, "interior ut1");
        assert_eq!(interior.x_pole_arcsec, 0.076609, "interior x");
        assert_eq!(interior.y_pole_arcsec, 0.282358, "interior y");
    }

    // ------------------------------------------------------------------
    // (4) pack/unpack 往復 + 長さ契約。
    //     unpack_records(pack_records(recs)) == recs（EopRecord PartialEq）。
    //     長さ = (1 ヘッダ件数 + 4 f64/レコード) * 8 バイト。
    // ------------------------------------------------------------------
    #[test]
    fn pack_unpack_round_trips_and_has_exact_length() {
        let recs = parse_eop_c04(EXCERPT).expect("excerpt parses");
        let bytes = pack_records(&recs);
        let restored = unpack_records(&bytes).expect("packed records round-trip");
        assert_eq!(restored, recs, "round-tripped records must equal original");

        let expected_len = (1 + 4 * recs.len()) * 8;
        assert_eq!(
            bytes.len(),
            expected_len,
            "packed length = (1 header n + 4 f64/record) * 8 bytes"
        );
    }

    // ------------------------------------------------------------------
    // (5) pack ヘッダ: 先頭 8 バイト（LE f64）はレコード件数（f64）。
    // ------------------------------------------------------------------
    #[test]
    fn pack_header_first_f64_is_record_count() {
        let recs = parse_eop_c04(EXCERPT).expect("excerpt parses");
        let bytes = pack_records(&recs);
        assert!(bytes.len() >= 8, "packed has at least the header f64");
        let header: [u8; 8] = bytes[..8].try_into().expect("8 bytes for header f64");
        let n = f64::from_le_bytes(header);
        // i32→f64 は無損失。3 レコードちょうど。
        let count = f64::from(u32::try_from(recs.len()).expect("record count fits in u32"));
        assert_eq!(n, count, "header f64 decodes to record count");
        assert_eq!(n, 3.0, "excerpt has 3 records");
    }

    // ------------------------------------------------------------------
    // (6) unpack 異常: 8 の倍数でない長さ / ヘッダ件数と残バイト長の不整合 → Err。
    //     具体的な variant は固定しない（is_err / Malformed* で緩く確認）。
    // ------------------------------------------------------------------
    #[test]
    fn unpack_rejects_non_f64_aligned_bytes() {
        let bytes = vec![0u8; 13]; // 8 の倍数でない
        assert!(
            unpack_records(&bytes).is_err(),
            "non-8-multiple length must error"
        );
    }

    #[test]
    fn unpack_rejects_header_length_mismatch() {
        // ヘッダで件数 1000 を宣言するが、後続データが伴わない（長さは 8 の倍数）。
        // 1 f64（ヘッダのみ）= 8 バイト。件数 1000 → 期待長 (1 + 4*1000)*8 と矛盾。
        let bytes = 1000.0_f64.to_le_bytes().to_vec();
        assert_eq!(bytes.len() % 8, 0, "fixture is f64-aligned");
        assert!(
            unpack_records(&bytes).is_err(),
            "header count vs data length mismatch must error"
        );
    }

    // ------------------------------------------------------------------
    // (7) コミット済み生成物が、コミット済み原データからの再生成と byte 一致。
    //     manifest 相対 include_str!/include_bytes! によるコンパイル時比較で
    //     cwd 非依存に検証する（runtime path 依存の verify_against_disk は
    //     CLI 用; ここでは原データ→pack の再現性のみを byte-for-byte 確認）。
    //     nutation テストの include_str!/include_bytes! パターンを踏襲。
    // ------------------------------------------------------------------
    #[test]
    fn committed_bin_matches_regeneration_from_source() {
        let records = parse_eop_c04(REAL_C04).expect("source parses");
        let regenerated = pack_records(&records);
        assert_eq!(
            regenerated.as_slice(),
            COMMITTED_BIN,
            "committed generated/eop/eopc04_14.bin must byte-match regeneration from the committed source"
        );
    }

    // ------------------------------------------------------------------
    // (8) build_artifact: メタデータ各フィールド & チェックサム契約。
    //     name/version/source/license の固定/非空、valid_from/to の日付境界、
    //     checksum が packed バイト列の sha256（64 桁小文字 hex）であること、
    //     provenance 完全。
    // ------------------------------------------------------------------
    #[test]
    fn build_artifact_metadata_and_checksum() {
        let (bytes, metadata) = build_artifact(EXCERPT).expect("excerpt builds an artifact");

        assert_eq!(metadata.name, "eop-iers-c04", "dataset name");
        assert_eq!(metadata.version, "EOP 14 C04", "dataset version");
        assert!(!metadata.source.is_empty(), "source non-empty");
        assert!(!metadata.license.is_empty(), "license non-empty");

        // 先頭 mjd 37665 = 1962-01-01, 末尾 mjd 58849 = 2020-01-01。
        assert_eq!(
            metadata.valid_from, "1962-01-01",
            "valid_from = first record"
        );
        assert_eq!(metadata.valid_to, "2020-01-01", "valid_to = last record");

        // checksum は packed バイト列上の sha256。
        assert_eq!(
            metadata.checksum,
            crate::checksum::sha256_hex(&bytes),
            "checksum is sha256 over the packed bytes"
        );
        assert_eq!(metadata.checksum.len(), 64, "sha256 hex is 64 chars");
        assert!(
            metadata
                .checksum
                .chars()
                .all(|c| c.is_ascii_hexdigit() && !c.is_ascii_uppercase()),
            "checksum is lowercase hex: {}",
            metadata.checksum
        );

        assert!(
            metadata.has_complete_provenance(),
            "all provenance fields non-empty: {metadata:?}"
        );
    }

    // ------------------------------------------------------------------
    // (9) build_artifact: 同一入力に対し決定的（bytes・metadata が一致）。
    // ------------------------------------------------------------------
    #[test]
    fn build_artifact_is_deterministic() {
        let (bytes_a, meta_a) = build_artifact(EXCERPT).expect("first build");
        let (bytes_b, meta_b) = build_artifact(EXCERPT).expect("second build");
        assert_eq!(bytes_a, bytes_b, "packed bytes are deterministic");
        assert_eq!(meta_a, meta_b, "metadata is deterministic");
    }

    // ------------------------------------------------------------------
    // (10) mjd_to_iso: JD オフセット & 年/月/日整形を既知日付で固定。
    //      月日の取り違えやオフセットの off-by-one は必ず落ちる。
    // ------------------------------------------------------------------
    #[test]
    fn mjd_to_iso_known_dates() {
        assert_eq!(mjd_to_iso(37665), "1962-01-01", "mjd 37665");
        assert_eq!(mjd_to_iso(58849), "2020-01-01", "mjd 58849");
        assert_eq!(mjd_to_iso(58850), "2020-01-02", "mjd 58850 (next day)");
    }

    // ------------------------------------------------------------------
    // (11) render_metadata: 全フィールドが `key = value` で出力に含まれる。
    //      各フィールドに別個の非空値を与え、取り違え/欠落を検出する。
    // ------------------------------------------------------------------
    #[test]
    fn render_metadata_includes_all_fields() {
        let metadata = DataSetMetadata {
            name: "the-name".to_string(),
            version: "the-version".to_string(),
            source: "the-source".to_string(),
            license: "the-license".to_string(),
            valid_from: "1900-01-02".to_string(),
            valid_to: "2099-12-31".to_string(),
            checksum: "the-checksum".to_string(),
        };
        let rendered = render_metadata(&metadata);
        assert!(rendered.contains("name = the-name"), "{rendered}");
        assert!(rendered.contains("version = the-version"), "{rendered}");
        assert!(rendered.contains("source = the-source"), "{rendered}");
        assert!(rendered.contains("license = the-license"), "{rendered}");
        assert!(rendered.contains("valid_from = 1900-01-02"), "{rendered}");
        assert!(rendered.contains("valid_to = 2099-12-31"), "{rendered}");
        assert!(rendered.contains("checksum = the-checksum"), "{rendered}");
    }

    // ------------------------------------------------------------------
    // (12) unpack 異常: ヘッダ件数が非有限/負/小数/上限超過 → MalformedPacked。
    //      NaN/負/巨大はヘッダ検証で（長さ検証の前に）弾かれる。
    // ------------------------------------------------------------------
    #[test]
    fn unpack_rejects_bad_header_count() {
        // (a) 負（後続レコードなし。長さ検証の前にヘッダ検証で弾く）。
        let negative = crate::packed::pack_f64_le(&[-1.0]);
        assert!(
            matches!(
                unpack_records(&negative),
                Err(XtaskError::MalformedPacked(_))
            ),
            "negative header count must be MalformedPacked"
        );

        // (b) 小数（fract != 0）。
        let fractional = crate::packed::pack_f64_le(&[1.5]);
        assert!(
            matches!(
                unpack_records(&fractional),
                Err(XtaskError::MalformedPacked(_))
            ),
            "fractional header count must be MalformedPacked"
        );

        // (c) NaN（非有限）。
        let nan = crate::packed::pack_f64_le(&[f64::NAN]);
        assert!(
            matches!(unpack_records(&nan), Err(XtaskError::MalformedPacked(_))),
            "NaN header count must be MalformedPacked"
        );

        // (d) 上限超過（MAX_EOP_RECORDS = 1e7 超）。
        let huge = crate::packed::pack_f64_le(&[2.0e9]);
        assert!(
            matches!(unpack_records(&huge), Err(XtaskError::MalformedPacked(_))),
            "absurdly large header count must be MalformedPacked"
        );
    }

    // ------------------------------------------------------------------
    // (13) parse: データ行がちょうど 7 トークン（year month day MJD x y UT1-UTC,
    //      末尾の追加列なし）でも 1 レコードとして取り込まれる（境界 len >= 7）。
    // ------------------------------------------------------------------
    #[test]
    fn parse_seven_token_data_row() {
        let recs = parse_eop_c04("1962   1   1  37665  -0.012700   0.213000   0.0326338")
            .expect("a 7-token data row is a valid data line");
        assert_eq!(
            recs.len(),
            1,
            "exactly one record from the single 7-token row"
        );
        assert_record(&recs[0], 37665, 0.0326338, -0.012700, 0.213000);
    }

    // ------------------------------------------------------------------
    // (14) parse: 暦範囲外（month/day/year が範囲外）の行は skip され、
    //      有効行のみ取り込まれる。range-guard の各 `||` 項を固定する。
    // ------------------------------------------------------------------
    #[test]
    fn parse_skips_out_of_calendar_range_rows() {
        let text = "\
1962   1   1  37665  -0.012700   0.213000   0.0326338   0.0017230   0.0   0.0   0.03   0.03   0.002   0.0014   0.0047   0.002
1962  13   1  37666  -0.015900   0.214100   0.0320547   0.0016690   0.0   0.0   0.03   0.03   0.002   0.0014   0.0047   0.002
1962   1  99  37667  -0.019000   0.215200   0.0315526   0.0015820   0.0   0.0   0.03   0.03   0.002   0.0014   0.0047   0.002
3000   1   1  37668  -0.022000   0.216300   0.0311435   0.0014960   0.0   0.0   0.03   0.03   0.002   0.0014   0.0047   0.002
2020   1   1  58849   0.076609   0.282358  -0.1771222   0.0004455   0.0   0.0   0.0001   0.0001   0.00001   0.00001   0.0001   0.0001
";
        let recs = parse_eop_c04(text).expect("two valid rows parse despite three bad rows");
        assert_eq!(
            recs.len(),
            2,
            "only the two in-calendar-range rows; month 13 / day 99 / year 3000 skipped"
        );
        assert_eq!(recs[0].mjd, 37665, "first kept row");
        assert_eq!(recs[1].mjd, 58849, "second kept row");
    }

    // ------------------------------------------------------------------
    // (15) unpack: ヘッダ件数 0 の有効バッファ（レコードなし）→ Ok(empty)。
    //      header == 0.0 は許容される（`< 0.0` 境界）。
    // ------------------------------------------------------------------
    #[test]
    fn unpack_accepts_zero_record_buffer() {
        let bytes = crate::packed::pack_f64_le(&[0.0]);
        let v = unpack_records(&bytes).expect("zero-record packed buffer is valid");
        assert!(v.is_empty(), "header count 0 yields no records");
    }

    // ------------------------------------------------------------------
    // (16) unpack: 長さ整合だが小数ヘッダ（1.5 + 4 レコード値）→ MalformedPacked。
    //      長さ検証を通過するため fract 検証のみが棄却できる（`||` の fract 項を固定）。
    // ------------------------------------------------------------------
    #[test]
    fn unpack_rejects_fractional_header_with_consistent_length() {
        let bytes = crate::packed::pack_f64_le(&[1.5, 58849.0, -0.1771222, 0.076609, 0.282358]);
        assert!(
            matches!(unpack_records(&bytes), Err(XtaskError::MalformedPacked(_))),
            "fractional header with length-consistent buffer must be MalformedPacked"
        );
    }
}
