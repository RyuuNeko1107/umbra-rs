//! 同梱 IERS EOP 14 C04 データの埋め込みコンシューマ（ISSUE-007 EOP part P2b）。
//!
//! xtask（part P2a）が生成・コミットした packed バイナリ
//! `generated/eop/eopc04_14.bin`（flat little-endian f64:
//! `[n_records, then per record: mjd, ut1_minus_utc_s, x_pole_arcsec, y_pole_arcsec]`）を
//! `include_bytes!` で取り込み、`generated/eop/metadata.txt` を `include_str!` で取り込んで
//! [`umbra_core::IersEopData`] を構築する。`nutation.rs` の packed コンシューマ慣習に倣う。
//!
//! 本コンシューマは `bundled-data` cargo feature（既定 on）でのみ提供される。off 時は
//! [`bundled_eop`] がシンボルごと消え、外部供給（`IersEopData::from_records`）が必須となる。

#[cfg(feature = "bundled-data")]
use umbra_core::{DataSetMetadata, EopRecord, IersEopData};

/// 同梱 packed バイナリ（P2a 生成物・byte-for-byte 契約）。
#[cfg(feature = "bundled-data")]
const PACKED: &[u8] = include_bytes!("../../../generated/eop/eopc04_14.bin");
/// 同梱メタデータ（provenance / checksum）。
#[cfg(feature = "bundled-data")]
const METADATA_TXT: &str = include_str!("../../../generated/eop/metadata.txt");

/// 同梱 IERS EOP 14 C04 を [`IersEopData`] として構築する（`bundled-data` feature）。
///
/// packed バイト列・metadata はビルド時に埋め込まれ、実行時ネットワークは行わない
/// （accuracy.md §5）。データの健全性は xtask `verify-generated` ゲートで担保されるため、
/// 復号は信頼前提（不整合は CI 段階で検出される）。
#[cfg(feature = "bundled-data")]
pub fn bundled_eop() -> IersEopData {
    let records = decode_packed(PACKED);
    let metadata = parse_metadata(METADATA_TXT);
    // EOP 14 C04 では系列識別子（series_version）と provenance 版（metadata.version）は同一。
    // 派生データで両者が分かれる場合はここを別供給にする。
    let series_version = metadata.version.clone();
    IersEopData::from_records(records, series_version, metadata)
        .expect("bundled EOP data is well-formed (xtask verify-generated gate)")
}

/// packed の `index` 番目の f64（little-endian）。
#[cfg(feature = "bundled-data")]
fn packed_f64(bytes: &[u8], index: usize) -> f64 {
    let start = index * 8;
    let octet: [u8; 8] = bytes[start..start + 8]
        .try_into()
        .expect("packed EOP length is a multiple of 8 (verify-generated gate)");
    f64::from_le_bytes(octet)
}

/// packed（`[n, then 各レコード = mjd, ut1, x, y]`）を復号する（xtask `pack_records` の逆・同一契約）。
#[cfg(feature = "bundled-data")]
#[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
fn decode_packed(bytes: &[u8]) -> Vec<EopRecord> {
    let count = packed_f64(bytes, 0).round() as usize;
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
    records
}

/// `metadata.txt`（`key = value` 行）を [`DataSetMetadata`] へ復元する（xtask `render_metadata` の逆）。
#[cfg(feature = "bundled-data")]
fn parse_metadata(text: &str) -> DataSetMetadata {
    let field = |key: &str| -> String {
        let prefix = format!("{key} = ");
        text.lines()
            .find_map(|line| line.strip_prefix(&prefix))
            .unwrap_or_else(|| panic!("bundled EOP metadata.txt is missing field '{key}'"))
            .to_string()
    };
    DataSetMetadata {
        name: field("name"),
        version: field("version"),
        source: field("source"),
        license: field("license"),
        valid_from: field("valid_from"),
        valid_to: field("valid_to"),
        checksum: field("checksum"),
    }
}

// ====================================================================
// テストモジュール。
//
// `bundled_eop()` は `bundled-data` feature でのみ提供されるため、テストモジュール自体を
// feature ゲートする（既定 features では on のため通常実行で評価される）。
// ====================================================================
#[cfg(all(test, feature = "bundled-data"))]
mod tests {
    // 検証値は IERS EOP 14 C04 のコミット済みデータの逐語転記（provenance 保持）。
    // 余剰桁は最近接 f64 へ丸められ値は不変のため、過剰精度リントを許可する。
    #![allow(clippy::excessive_precision)]

    use umbra_core::EarthOrientation;
    use umbra_core::UtcInstant;

    /// 同梱 EOP コンシューマ（実装予定の公開 IF）。crate 内テストのため `crate::eop` 経由で参照する
    /// （クレート名 `umbra_ephemeris::` は自クレート内からは解決されないため）。再エクスポートは
    /// 実装工程で `lib.rs` に追加され、外部からは `umbra_ephemeris::bundled_eop` で参照可能になる。
    use crate::eop::bundled_eop;

    // ---- 定数・補助 -------------------------------------------------------

    /// arcsec → radian 係数（1" = π/648000 rad ≈ 4.84813681109536e-6）。
    const ARCSEC_TO_RAD: f64 = std::f64::consts::PI / 648_000.0;

    /// グレゴリオ暦（UTC）の瞬時を構築する補助。
    fn utc(y: i32, mo: u8, d: u8, h: u8, mi: u8, s: f64) -> UtcInstant {
        UtcInstant::from_gregorian(y, mo, d, h, mi, s).expect("valid calendar date")
    }

    // ==================================================================
    // 1. 構築 ＆ coverage 境界
    //    bundled_eop() が値を返し、coverage が [1962-01-01 0h, 2026-01-05 0h]。
    //    （end はコミット済みデータの版に依存し、データ更新で追従する。）
    // ==================================================================
    #[test]
    fn bundled_eop_builds_with_expected_coverage() {
        let data = bundled_eop();
        let range = data.coverage();
        // coverage start = 最初のレコード日 1962-01-01 0h UTC。
        let want_start = utc(1962, 1, 1, 0, 0, 0.0);
        // coverage end = 最後のレコード日 2026-01-05 0h UTC（コミット済みデータ依存）。
        let want_end = utc(2026, 1, 5, 0, 0, 0.0);
        assert!(
            (range.start.jd2().jd() - want_start.jd2().jd()).abs() < 1e-9,
            "coverage start jd = {}, want {} (1962-01-01 0h)",
            range.start.jd2().jd(),
            want_start.jd2().jd()
        );
        assert!(
            (range.end.jd2().jd() - want_end.jd2().jd()).abs() < 1e-9,
            "coverage end jd = {}, want {} (2026-01-05 0h, tracks committed data)",
            range.end.jd2().jd(),
            want_end.jd2().jd()
        );
    }

    // ==================================================================
    // 2. series_version
    //    コミット済みデータの系列版 == "EOP 14 C04"。
    // ==================================================================
    #[test]
    fn bundled_eop_series_version_is_eop14c04() {
        let data = bundled_eop();
        assert_eq!(data.series_version(), "EOP 14 C04");
    }

    // ==================================================================
    // 3. 既知 UT1−UTC（フォーマット契約・強）
    //    2020-01-01 0h UTC で UT1−UTC == -0.1771222 s。これは厳密日ルックアップであり、
    //    packed のレイアウト（フィールド順 mjd/ut1/xp/yp）が xtask と byte-for-byte
    //    一致しなければ失敗する強い契約テスト。
    // ==================================================================
    #[test]
    fn bundled_eop_known_ut1_minus_utc_2020() {
        let data = bundled_eop();
        let v = data
            .ut1_minus_utc(utc(2020, 1, 1, 0, 0, 0.0))
            .expect("2020-01-01 is within coverage");
        let want = -0.1771222;
        assert!(
            (v - want).abs() < 1e-9,
            "UT1-UTC(2020-01-01 0h) = {v}, want {want} (byte-for-byte format contract)"
        );
    }

    // ==================================================================
    // 4. 既知 極運動
    //    2020-01-01 0h UTC で (xp, yp) == (0.076609, 0.282358) arcsec → rad。
    //    各成分 arcsec × (π/648000) と ~1e-12 以内で一致。
    // ==================================================================
    #[test]
    fn bundled_eop_known_polar_motion_2020() {
        let data = bundled_eop();
        let (umbra_core::Radians(xp), umbra_core::Radians(yp)) = data
            .polar_motion(utc(2020, 1, 1, 0, 0, 0.0))
            .expect("2020-01-01 is within coverage");
        let want_xp = 0.076609 * ARCSEC_TO_RAD;
        let want_yp = 0.282358 * ARCSEC_TO_RAD;
        assert!(
            (xp - want_xp).abs() < 1e-12,
            "xp(2020-01-01 0h) = {xp}, want {want_xp}"
        );
        assert!(
            (yp - want_yp).abs() < 1e-12,
            "yp(2020-01-01 0h) = {yp}, want {want_yp}"
        );
    }

    // ==================================================================
    // 5. coverage 外
    //    最初のレコード（1962-01-01）より前（1961-12-31）は
    //    MissingEarthOrientationData。
    // ==================================================================
    #[test]
    fn bundled_eop_before_coverage_is_missing() {
        let data = bundled_eop();
        let before = utc(1961, 12, 31, 0, 0, 0.0);
        assert_eq!(
            data.ut1_minus_utc(before).unwrap_err(),
            umbra_core::TimeError::MissingEarthOrientationData
        );
    }

    // ==================================================================
    // 6. metadata provenance
    //    name == "eop-iers-c04", version == "EOP 14 C04",
    //    checksum == コミット済み bin の SHA-256, provenance 完全。
    // ==================================================================
    #[test]
    fn bundled_eop_metadata_provenance() {
        let data = bundled_eop();
        let md = data.metadata();
        assert_eq!(md.name, "eop-iers-c04", "metadata.name");
        assert_eq!(md.version, "EOP 14 C04", "metadata.version");
        assert_eq!(
            md.checksum, "4d12cc1d3dcae39a8db58690a0d6e669def893903ba74c23995ced1e5e354ac3",
            "metadata.checksum"
        );
        assert!(
            md.has_complete_provenance(),
            "bundled metadata must have complete provenance"
        );
    }

    // ==================================================================
    // 7. 内挿サニティ（弱・任意）
    //    2020-01-01 12:00 の UT1−UTC は 2020-01-01 と 2020-01-02 の値の間にあり有限。
    //    強い固定は test 3。ここは線形補間が破綻していない（範囲内・有限）ことのみ確認。
    // ==================================================================
    #[test]
    fn bundled_eop_interior_interpolation_is_between_neighbors() {
        let data = bundled_eop();
        let day1 = data
            .ut1_minus_utc(utc(2020, 1, 1, 0, 0, 0.0))
            .expect("2020-01-01 within coverage");
        let day2 = data
            .ut1_minus_utc(utc(2020, 1, 2, 0, 0, 0.0))
            .expect("2020-01-02 within coverage");
        let mid = data
            .ut1_minus_utc(utc(2020, 1, 1, 12, 0, 0.0))
            .expect("2020-01-01 12:00 within coverage");
        assert!(
            mid.is_finite(),
            "interior UT1-UTC must be finite, got {mid}"
        );
        let (lo, hi) = if day1 <= day2 {
            (day1, day2)
        } else {
            (day2, day1)
        };
        assert!(
            (lo..=hi).contains(&mid),
            "interior UT1-UTC = {mid} must be between neighbors [{lo}, {hi}]"
        );
    }
}
