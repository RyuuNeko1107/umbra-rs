//! DAF/SPK 構造解析（ISSUE-036 S1）。
//!
//! DAF（Double Array File）= SPK の物理形式。1024 バイト固定長レコード列。先頭がファイルレコード
//! （マジック `DAF/SPK `・ND/NI・FWARD/BWARD/FREE・バイトオーダ `LTL-IEEE`）、FWARD から続く
//! サマリレコード列が各セグメントの記述子（2 double + 6 int）を持つ。本スライスは記述子の抽出のみ
//! （Chebyshev 係数の読取り・評価は S2 以降）。
//!
//! 仕様: NAIF DAF Required Reading / SPK Required Reading。最小版は **LTL-IEEE のみ対応**。
//! 実 DE440s（ND=2, NI=6, 全 type2, ET 1849–2150）で検証する。

// S1（構造解析）は単体で完結し、非テストビルドからの呼び出しは S2（Chebyshev 評価）/
// S3（JplEphemeris）で配線される。それまでは crate 内利用者が無く dead_code 警告となるため
// モジュール単位で許可する（配線後に解除）。テストモジュールは全項目を使用する。
#![allow(dead_code)]

use crate::ephemeris::EphemerisError;

/// DAF レコード長（バイト）。
const RECORD_BYTES: usize = 1024;

/// SPK 固定: 各サマリの double 成分数（ND）。初期/終期 ET。
const SPK_ND: i32 = 2;
/// SPK 固定: 各サマリの integer 成分数（NI）。target/center/frame/type/start_addr/end_addr。
const SPK_NI: i32 = 6;
/// サマリレコード先頭の制御 double（NEXT/PREV/NSUM）の総バイト。
const SUMMARY_CTRL_BYTES: usize = 24;
/// 1 サマリのバイト数 = (ND + ceil(NI/2)) double = (2 + 3) * 8 = 40。
const SUMMARY_BYTES: usize = 40;
/// 1 レコードに収まるサマリの上限 = floor((1024 - 24) / 40) = 25。
const MAX_SUMMARIES_PER_RECORD: usize = (RECORD_BYTES - SUMMARY_CTRL_BYTES) / SUMMARY_BYTES;

/// SPK セグメント記述子（DAF サマリ 1 件）。アドレスは 1 始まりの DAF ワード（8 バイト）番地。
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct SpkSegment {
    /// セグメント有効区間の開始時刻（ET=TDB 秒, J2000 基準）。
    pub start_et: f64,
    /// セグメント有効区間の終了時刻（ET=TDB 秒）。
    pub end_et: f64,
    /// 対象天体の NAIF ID（10=Sun, 301=Moon, 399=Earth, 3=EMB, …）。
    pub target: i32,
    /// 中心天体の NAIF ID（0=SSB, 3=EMB, …）。
    pub center: i32,
    /// 参照フレーム ID（1=J2000/ICRF）。
    pub frame: i32,
    /// SPK データ型（2=Chebyshev 位置, 3=Chebyshev 位置+速度）。
    pub data_type: i32,
    /// データ先頭の DAF ワード番地（1 始まり）。
    pub start_addr: i32,
    /// データ末尾の DAF ワード番地（1 始まり, 含む）。
    pub end_addr: i32,
}

fn malformed(msg: impl Into<String>) -> EphemerisError {
    EphemerisError::MalformedSpk(msg.into())
}

/// `bytes[off..off+4]` をリトルエンディアン i32 として読む。範囲外は MalformedSpk。
fn read_i32(bytes: &[u8], off: usize) -> Result<i32, EphemerisError> {
    let slice = bytes
        .get(off..off + 4)
        .ok_or_else(|| malformed(format!("truncated i32 at byte {off}")))?;
    Ok(i32::from_le_bytes(
        slice.try_into().expect("slice is 4 bytes"),
    ))
}

/// `bytes[off..off+8]` をリトルエンディアン f64 として読む。範囲外は MalformedSpk。
fn read_f64(bytes: &[u8], off: usize) -> Result<f64, EphemerisError> {
    let slice = bytes
        .get(off..off + 8)
        .ok_or_else(|| malformed(format!("truncated f64 at byte {off}")))?;
    Ok(f64::from_le_bytes(
        slice.try_into().expect("slice is 8 bytes"),
    ))
}

/// 制御 double（NEXT/PREV/NSUM・レコード番号）を非負整数のレコード番号へ。
/// 非有限・負・非整数は MalformedSpk（DAF はこれらを整数値として格納する）。
/// 非負整数であることを検証済みのためキャストは安全（truncation/sign-loss なし）。
#[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
fn record_number(value: f64, what: &str) -> Result<usize, EphemerisError> {
    if !value.is_finite() || value < 0.0 || value.fract() != 0.0 {
        return Err(malformed(format!(
            "{what} is not a non-negative integer: {value}"
        )));
    }
    Ok(value as usize)
}

/// DAF/SPK バイト列からセグメント記述子列を解析する（S1: 構造のみ）。
///
/// 検証: 最小長・マジック `DAF/SPK `・バイトオーダ `LTL-IEEE`・ND=2/NI=6・各レコード境界。
/// 不正/未対応は [`EphemerisError::MalformedSpk`]。
pub(crate) fn parse_spk_segments(bytes: &[u8]) -> Result<Vec<SpkSegment>, EphemerisError> {
    // --- ファイルレコード（先頭 1024B）の検証 ---
    if bytes.len() < RECORD_BYTES {
        return Err(malformed(format!(
            "file shorter than one DAF record ({} < {RECORD_BYTES})",
            bytes.len()
        )));
    }
    if &bytes[0..8] != b"DAF/SPK " {
        return Err(malformed("bad DAF magic (expected \"DAF/SPK \")"));
    }
    // バイトオーダ: 最小版は LTL-IEEE のみ対応。
    if &bytes[88..96] != b"LTL-IEEE" {
        return Err(malformed(format!(
            "unsupported byte order (expected \"LTL-IEEE\"): {:?}",
            String::from_utf8_lossy(&bytes[88..96])
        )));
    }
    let nd = read_i32(bytes, 8)?;
    let ni = read_i32(bytes, 12)?;
    if nd != SPK_ND || ni != SPK_NI {
        return Err(malformed(format!(
            "unsupported summary shape ND={nd}/NI={ni} (expected {SPK_ND}/{SPK_NI})"
        )));
    }
    // FWARD は最初のサマリレコード番号（1 始まり）。有効 DAF は空 SPK でも 1 件のサマリ
    // レコードを持つため必ず >= 1。0/負は不正（NEXT=0 のチェーン終端とは別概念）。
    let fward = read_i32(bytes, 76)?;
    if fward < 1 {
        return Err(malformed(format!("FWARD must be >= 1 (got {fward})")));
    }

    // ファイルに完全に存在するレコード数（末尾の端数は切り捨て）。
    let total_records = bytes.len() / RECORD_BYTES;

    // --- サマリレコードのチェーンを FWARD から辿る ---
    let mut segments = Vec::new();
    let mut rec = record_number(f64::from(fward), "FWARD")?;
    let mut steps = 0usize;
    while rec != 0 {
        // 1 始まりレコード番号。レコードが完全に存在することを要求。
        if rec < 1 || rec > total_records {
            return Err(malformed(format!(
                "summary record {rec} out of range (file has {total_records} records)"
            )));
        }
        // チェーン長はレコード総数で頭打ち（循環参照の保護）。
        steps += 1;
        if steps > total_records {
            return Err(malformed(
                "summary record chain does not terminate (cycle?)",
            ));
        }

        let base = (rec - 1) * RECORD_BYTES;
        // 制御 double: NEXT(base)・PREV(base+8, 前方走査では未使用)・NSUM(base+16)。
        let next = read_f64(bytes, base)?;
        let nsum = record_number(read_f64(bytes, base + 16)?, "NSUM")?;
        if nsum > MAX_SUMMARIES_PER_RECORD {
            return Err(malformed(format!(
                "NSUM={nsum} exceeds per-record capacity {MAX_SUMMARIES_PER_RECORD}"
            )));
        }

        for k in 0..nsum {
            let s = base + SUMMARY_CTRL_BYTES + k * SUMMARY_BYTES;
            segments.push(SpkSegment {
                start_et: read_f64(bytes, s)?,
                end_et: read_f64(bytes, s + 8)?,
                target: read_i32(bytes, s + 16)?,
                center: read_i32(bytes, s + 20)?,
                frame: read_i32(bytes, s + 24)?,
                data_type: read_i32(bytes, s + 28)?,
                start_addr: read_i32(bytes, s + 32)?,
                end_addr: read_i32(bytes, s + 36)?,
            });
        }

        rec = record_number(next, "NEXT")?;
    }

    Ok(segments)
}

#[cfg(test)]
mod tests {
    use super::*;

    // ------------------------------------------------------------------
    // 仕様出典: NAIF "DAF Required Reading" / "SPK Required Reading"。
    //   ファイルレコード（先頭 1024B）: LOCIDW(8) ND(i32) NI(i32) LOCIFN(60)
    //     FWARD(i32) BWARD(i32) FREE(i32) LOCFMT(8) PRENUL(603) FTPSTR(28) PSTNUL(297)。
    //   サマリレコード（FWARD レコードから, byte offset=(rec-1)*1024）:
    //     先頭 3 double = NEXT, PREV, NSUM。続いて NSUM 個のサマリ（offset 24 から）。
    //     SS = ND + ceil(NI/2) = 2 + 3 = 5 double = 40 バイト。
    //     各サマリ: [start_et:f64][end_et:f64][i32 x6 = target,center,frame,
    //       data_type,start_addr,end_addr]。NEXT を 0.0 まで辿る。
    //   LTL-IEEE のとき全 i32/f64 はリトルエンディアン（本最小版は LTL-IEEE のみ対応）。
    // 合成フィクスチャは全て上記レイアウトをリトルエンディアンで明示的に組み立てる。
    // ------------------------------------------------------------------

    /// 1 サマリのサイズ（バイト）= 5 double。
    const SUMMARY_BYTES: usize = 40;
    /// サマリレコード先頭の制御 double（NEXT/PREV/NSUM）の総バイト。
    const SUMMARY_CTRL_BYTES: usize = 24;

    /// 合成サマリの全フィールドを保持する素朴な記述（テスト入力の組立用）。
    #[derive(Clone, Copy)]
    struct Sum {
        start_et: f64,
        end_et: f64,
        target: i32,
        center: i32,
        frame: i32,
        data_type: i32,
        start_addr: i32,
        end_addr: i32,
    }

    /// 既定値（各テストで必要フィールドだけ上書きする）。
    fn sample_sum() -> Sum {
        Sum {
            start_et: -1.0e9,
            end_et: 1.0e9,
            target: 10,
            center: 0,
            frame: 1,
            data_type: 2,
            start_addr: 1,
            end_addr: 1000,
        }
    }

    /// 指定オフセットに i32 をリトルエンディアンで書く。
    fn put_i32(buf: &mut [u8], off: usize, v: i32) {
        buf[off..off + 4].copy_from_slice(&v.to_le_bytes());
    }

    /// 指定オフセットに f64 をリトルエンディアンで書く。
    fn put_f64(buf: &mut [u8], off: usize, v: f64) {
        buf[off..off + 8].copy_from_slice(&v.to_le_bytes());
    }

    /// 指定オフセットに ASCII バイト列を書く（長さは呼び出し側責任）。
    fn put_ascii(buf: &mut [u8], off: usize, s: &[u8]) {
        buf[off..off + s.len()].copy_from_slice(s);
    }

    /// ファイルレコード（1024B）を組み立てる。既定は健全な DAF/SPK・LTL-IEEE・ND=2/NI=6。
    /// 異常系テストは戻り値を後から上書きして使う。
    fn build_file_record(fward: i32, bward: i32, free: i32) -> Vec<u8> {
        let mut buf = vec![0u8; RECORD_BYTES];
        put_ascii(&mut buf, 0, b"DAF/SPK "); // LOCIDW（末尾スペース込み 8 文字）
        put_i32(&mut buf, 8, 2); // ND
        put_i32(&mut buf, 12, 6); // NI
        put_ascii(&mut buf, 16, b"UMBRA-TEST"); // LOCIFN（残りはゼロ）
        put_i32(&mut buf, 76, fward); // FWARD
        put_i32(&mut buf, 80, bward); // BWARD
        put_i32(&mut buf, 84, free); // FREE
        put_ascii(&mut buf, 88, b"LTL-IEEE"); // LOCFMT
        buf
    }

    /// サマリレコード（1024B）を組み立てる。`next` は次レコード番号（無ければ 0.0）。
    /// `prev` は前レコード番号。`sums` を offset 24 から 40B 刻みで並べる。
    /// `sums.len()` は最大 25（1 レコード容量）の小さな件数のため f64 化に精度損失なし。
    #[allow(clippy::cast_precision_loss)]
    fn build_summary_record(next: f64, prev: f64, sums: &[Sum]) -> Vec<u8> {
        let mut buf = vec![0u8; RECORD_BYTES];
        put_f64(&mut buf, 0, next); // NEXT
        put_f64(&mut buf, 8, prev); // PREV
        put_f64(&mut buf, 16, sums.len() as f64); // NSUM
        for (k, s) in sums.iter().enumerate() {
            let base = SUMMARY_CTRL_BYTES + k * SUMMARY_BYTES;
            put_f64(&mut buf, base, s.start_et);
            put_f64(&mut buf, base + 8, s.end_et);
            put_i32(&mut buf, base + 16, s.target);
            put_i32(&mut buf, base + 20, s.center);
            put_i32(&mut buf, base + 24, s.frame);
            put_i32(&mut buf, base + 28, s.data_type);
            put_i32(&mut buf, base + 32, s.start_addr);
            put_i32(&mut buf, base + 36, s.end_addr);
        }
        buf
    }

    /// 期待値（合成 Sum）と実値（SpkSegment）が全フィールド一致することを検証する。
    /// 合成では f64 もビット一致のため `==` で良い。
    fn assert_seg_eq(actual: &SpkSegment, expected: &Sum) {
        assert_eq!(actual.start_et, expected.start_et, "start_et");
        assert_eq!(actual.end_et, expected.end_et, "end_et");
        assert_eq!(actual.target, expected.target, "target");
        assert_eq!(actual.center, expected.center, "center");
        assert_eq!(actual.frame, expected.frame, "frame");
        assert_eq!(actual.data_type, expected.data_type, "data_type");
        assert_eq!(actual.start_addr, expected.start_addr, "start_addr");
        assert_eq!(actual.end_addr, expected.end_addr, "end_addr");
    }

    // ==================================================================
    // 正常系（合成）
    // ==================================================================

    /// 単一サマリレコードに 1 件のセグメント。全フィールドが正しく解析されること。
    #[test]
    fn parses_single_segment() {
        let file_rec = build_file_record(2, 2, 0); // サマリは 2 レコード目（FWARD=2）
        let expected = Sum {
            start_et: -4.0e9,
            end_et: 4.0e9,
            target: 301,
            center: 3,
            frame: 1,
            data_type: 2,
            start_addr: 769,
            end_addr: 12345,
        };
        let summary_rec = build_summary_record(0.0, 0.0, &[expected]); // NEXT=0 でチェーン終端

        let mut bytes = file_rec;
        bytes.extend_from_slice(&summary_rec);

        let segs = parse_spk_segments(&bytes).expect("健全な DAF/SPK は解析成功すべき");
        assert_eq!(segs.len(), 1);
        assert_seg_eq(&segs[0], &expected);
    }

    /// 単一サマリレコードに複数（3 件）のセグメント。順序・件数・各値が保たれること。
    #[test]
    fn parses_multiple_segments_in_one_record() {
        let file_rec = build_file_record(2, 2, 0);
        let s0 = Sum {
            target: 10,
            start_et: -1.0e9,
            end_et: 1.0e9,
            start_addr: 1001,
            end_addr: 2000,
            ..sample_sum()
        };
        let s1 = Sum {
            target: 399,
            center: 3,
            start_et: -2.0e9,
            end_et: 2.0e9,
            start_addr: 2001,
            end_addr: 3000,
            ..sample_sum()
        };
        let s2 = Sum {
            target: 3,
            start_et: -3.0e9,
            end_et: 3.0e9,
            start_addr: 3001,
            end_addr: 4000,
            ..sample_sum()
        };
        let sums = [s0, s1, s2];
        let summary_rec = build_summary_record(0.0, 0.0, &sums);

        let mut bytes = file_rec;
        bytes.extend_from_slice(&summary_rec);

        let segs = parse_spk_segments(&bytes).expect("解析成功すべき");
        assert_eq!(segs.len(), 3);
        for (got, want) in segs.iter().zip(sums.iter()) {
            assert_seg_eq(got, want);
        }
    }

    /// 2 レコードチェーン（NEXT で 2 レコード目のサマリレコードに続く）を最後まで辿ること。
    #[test]
    fn follows_summary_record_chain() {
        // レイアウト: rec1=ファイル, rec2=サマリ(NEXT→3), rec3=サマリ(NEXT=0)。
        let file_rec = build_file_record(2, 3, 0); // FWARD=2, BWARD=3

        let a = Sum {
            target: 10,
            start_addr: 100,
            end_addr: 200,
            ..sample_sum()
        };
        let b = Sum {
            target: 301,
            start_addr: 200,
            end_addr: 300,
            ..sample_sum()
        };
        // rec2: NEXT=3（次レコードへ）, PREV=0。
        let rec2 = build_summary_record(3.0, 0.0, &[a]);
        // rec3: NEXT=0（終端）, PREV=2。
        let rec3 = build_summary_record(0.0, 2.0, &[b]);

        let mut bytes = file_rec;
        bytes.extend_from_slice(&rec2);
        bytes.extend_from_slice(&rec3);

        let segs = parse_spk_segments(&bytes).expect("チェーンを辿り解析成功すべき");
        assert_eq!(segs.len(), 2);
        assert_seg_eq(&segs[0], &a);
        assert_seg_eq(&segs[1], &b);
    }

    /// FWARD と BWARD が同一の単一サマリレコードでも NEXT=0 で正しく終端すること。
    /// （サマリ 0 件のレコード = NSUM=0 を許容しエラーにしないこと。）
    #[test]
    fn parses_empty_summary_record() {
        let file_rec = build_file_record(2, 2, 0);
        let summary_rec = build_summary_record(0.0, 0.0, &[]); // NSUM=0

        let mut bytes = file_rec;
        bytes.extend_from_slice(&summary_rec);

        let segs = parse_spk_segments(&bytes).expect("空サマリレコードは 0 件で成功すべき");
        assert!(segs.is_empty());
    }

    // ==================================================================
    // 異常系（MalformedSpk）。メッセージ内容には依存せず変種のみ確認。
    // ==================================================================

    /// 最小長（1024B）未満は MalformedSpk。
    #[test]
    fn rejects_too_short() {
        let bytes = vec![0u8; RECORD_BYTES - 1];
        let err = parse_spk_segments(&bytes).unwrap_err();
        assert!(matches!(err, EphemerisError::MalformedSpk(_)));
    }

    /// 空バイト列も MalformedSpk。
    #[test]
    fn rejects_empty() {
        let err = parse_spk_segments(&[]).unwrap_err();
        assert!(matches!(err, EphemerisError::MalformedSpk(_)));
    }

    /// マジック（LOCIDW）不一致は MalformedSpk。
    #[test]
    fn rejects_bad_magic() {
        let mut file_rec = build_file_record(2, 2, 0);
        put_ascii(&mut file_rec, 0, b"NAIF/DAF"); // 誤マジック
        let summary_rec = build_summary_record(0.0, 0.0, &[sample_sum()]);
        let mut bytes = file_rec;
        bytes.extend_from_slice(&summary_rec);

        let err = parse_spk_segments(&bytes).unwrap_err();
        assert!(matches!(err, EphemerisError::MalformedSpk(_)));
    }

    /// LOCFMT が "BIG-IEEE"（未対応バイトオーダ）は MalformedSpk。
    #[test]
    fn rejects_big_endian_format() {
        let mut file_rec = build_file_record(2, 2, 0);
        put_ascii(&mut file_rec, 88, b"BIG-IEEE"); // 未対応
        let summary_rec = build_summary_record(0.0, 0.0, &[sample_sum()]);
        let mut bytes = file_rec;
        bytes.extend_from_slice(&summary_rec);

        let err = parse_spk_segments(&bytes).unwrap_err();
        assert!(matches!(err, EphemerisError::MalformedSpk(_)));
    }

    /// LOCFMT が未知文字列は MalformedSpk。
    #[test]
    fn rejects_unknown_format() {
        let mut file_rec = build_file_record(2, 2, 0);
        put_ascii(&mut file_rec, 88, b"WUT-IEEE"); // 未知
        let summary_rec = build_summary_record(0.0, 0.0, &[sample_sum()]);
        let mut bytes = file_rec;
        bytes.extend_from_slice(&summary_rec);

        let err = parse_spk_segments(&bytes).unwrap_err();
        assert!(matches!(err, EphemerisError::MalformedSpk(_)));
    }

    /// ND≠2 は MalformedSpk（SPK は ND=2 固定）。
    #[test]
    fn rejects_wrong_nd() {
        let mut file_rec = build_file_record(2, 2, 0);
        put_i32(&mut file_rec, 8, 3); // ND=3（不正）
        let summary_rec = build_summary_record(0.0, 0.0, &[sample_sum()]);
        let mut bytes = file_rec;
        bytes.extend_from_slice(&summary_rec);

        let err = parse_spk_segments(&bytes).unwrap_err();
        assert!(matches!(err, EphemerisError::MalformedSpk(_)));
    }

    /// NI≠6 は MalformedSpk（SPK は NI=6 固定）。
    #[test]
    fn rejects_wrong_ni() {
        let mut file_rec = build_file_record(2, 2, 0);
        put_i32(&mut file_rec, 12, 5); // NI=5（不正）
        let summary_rec = build_summary_record(0.0, 0.0, &[sample_sum()]);
        let mut bytes = file_rec;
        bytes.extend_from_slice(&summary_rec);

        let err = parse_spk_segments(&bytes).unwrap_err();
        assert!(matches!(err, EphemerisError::MalformedSpk(_)));
    }

    /// FWARD がファイル範囲外（存在しないレコード番号）を指す場合は MalformedSpk。
    #[test]
    fn rejects_fward_out_of_range() {
        // FWARD=5 だがファイルは 2 レコード分しか無い。
        let file_rec = build_file_record(5, 5, 0);
        let summary_rec = build_summary_record(0.0, 0.0, &[sample_sum()]);
        let mut bytes = file_rec;
        bytes.extend_from_slice(&summary_rec);

        let err = parse_spk_segments(&bytes).unwrap_err();
        assert!(matches!(err, EphemerisError::MalformedSpk(_)));
    }

    /// FWARD=0（< 1, 1 始まり番地として不正）は MalformedSpk。
    #[test]
    fn rejects_fward_zero() {
        let file_rec = build_file_record(0, 0, 0);
        let summary_rec = build_summary_record(0.0, 0.0, &[sample_sum()]);
        let mut bytes = file_rec;
        bytes.extend_from_slice(&summary_rec);

        let err = parse_spk_segments(&bytes).unwrap_err();
        assert!(matches!(err, EphemerisError::MalformedSpk(_)));
    }

    /// NEXT が範囲外レコードを指す場合（チェーンが境界を越える）は MalformedSpk。
    #[test]
    fn rejects_chain_next_out_of_range() {
        let file_rec = build_file_record(2, 2, 0);
        // rec2: NEXT=9（存在しないレコード）。ファイルは 2 レコード分のみ。
        let summary_rec = build_summary_record(9.0, 0.0, &[sample_sum()]);
        let mut bytes = file_rec;
        bytes.extend_from_slice(&summary_rec);

        let err = parse_spk_segments(&bytes).unwrap_err();
        assert!(matches!(err, EphemerisError::MalformedSpk(_)));
    }

    /// サマリレコードがファイル境界を越える（FWARD レコードの末尾が総バイトを超過）は MalformedSpk。
    #[test]
    fn rejects_summary_record_crosses_boundary() {
        // FWARD=2 を宣言するが、バイト列はファイルレコード（1 レコード）分しか無い。
        // → 2 レコード目（byte 1024..2048）が存在せず境界を越える。
        let bytes = build_file_record(2, 2, 0); // 1024 バイトのみ
        let err = parse_spk_segments(&bytes).unwrap_err();
        assert!(matches!(err, EphemerisError::MalformedSpk(_)));
    }

    /// NSUM が 1 レコードの上限（25 件）を超える値で、サマリ領域がレコードを越える場合は MalformedSpk。
    /// 1 レコードあたり最大 floor((1024-24)/40)=25 件。
    #[test]
    fn rejects_nsum_exceeds_record_capacity() {
        let file_rec = build_file_record(2, 2, 0);
        let mut summary_rec = build_summary_record(0.0, 0.0, &[]);
        put_f64(&mut summary_rec, 16, 26.0); // NSUM=26（容量 25 超過）
        let mut bytes = file_rec;
        bytes.extend_from_slice(&summary_rec);

        let err = parse_spk_segments(&bytes).unwrap_err();
        assert!(matches!(err, EphemerisError::MalformedSpk(_)));
    }

    /// NSUM（offset 16）が負値は MalformedSpk（件数は非負整数で負は不正）。
    #[test]
    fn rejects_nsum_negative() {
        let file_rec = build_file_record(2, 2, 0);
        let mut summary_rec = build_summary_record(0.0, 0.0, &[]);
        put_f64(&mut summary_rec, 16, -1.0); // NSUM=-1（非負整数違反）
        let mut bytes = file_rec;
        bytes.extend_from_slice(&summary_rec);

        let err = parse_spk_segments(&bytes).unwrap_err();
        assert!(matches!(err, EphemerisError::MalformedSpk(_)));
    }

    /// NEXT（offset 0）が非有限（NaN/∞）または非整数は MalformedSpk
    /// （NEXT は次サマリレコード番号で整数・有限でなければならない）。
    #[test]
    fn rejects_next_non_finite_or_non_integer() {
        for next_value in [f64::NAN, f64::INFINITY, 1.5] {
            let file_rec = build_file_record(2, 2, 0);
            // NEXT に不正値を書く（build_summary_record 第1引数が NEXT）。
            let summary_rec = build_summary_record(next_value, 0.0, &[sample_sum()]);
            let mut bytes = file_rec;
            bytes.extend_from_slice(&summary_rec);

            let err = parse_spk_segments(&bytes).unwrap_err();
            assert!(
                matches!(err, EphemerisError::MalformedSpk(_)),
                "NEXT={next_value} は MalformedSpk のはず"
            );
        }
    }

    /// FWARD が負値は MalformedSpk（FWARD は 1 始まりレコード番号で負は不正）。
    #[test]
    fn rejects_fward_negative() {
        let file_rec = build_file_record(-1, -1, 0); // FWARD=-1（1 始まり違反）
        let summary_rec = build_summary_record(0.0, 0.0, &[sample_sum()]);
        let mut bytes = file_rec;
        bytes.extend_from_slice(&summary_rec);

        let err = parse_spk_segments(&bytes).unwrap_err();
        assert!(matches!(err, EphemerisError::MalformedSpk(_)));
    }

    // ==================================================================
    // 正常系（実ファイル・ゲート）: data/spk/de440s.bsp 存在時のみ実行。
    //   不在時は eprintln! で告知して早期 return（CI 非同梱・ISSUE-036）。
    //   パス: CARGO_MANIFEST_DIR (crates/umbra-ephemeris) から ../../data/spk/de440s.bsp。
    // 既知特性（受入オラクル）:
    //   DAF/SPK・LTL-IEEE・ND=2・NI=6、全セグメント data_type==2、frame==1、
    //   target に 10/301/399/3 を含む、各 start_et<end_et、ET 範囲 1849–2150
    //   （おおよそ J2000 ET 秒で -4.74e9 〜 +4.74e9）。
    // ==================================================================

    /// 実 DE440s を解析し、既知特性を満たすこと（ファイル存在時のみ）。
    #[test]
    fn parses_real_de440s() {
        let path = concat!(env!("CARGO_MANIFEST_DIR"), "/../../data/spk/de440s.bsp");
        let bytes = match std::fs::read(path) {
            Ok(b) => b,
            Err(e) => {
                eprintln!(
                    "skip parses_real_de440s: {path} を読めない（{e}）。\
                     実ファイルは CI 非同梱（ISSUE-036）。"
                );
                return;
            }
        };

        let segs = parse_spk_segments(&bytes).expect("実 DE440s は解析成功すべき");
        assert!(!segs.is_empty(), "セグメントが 1 件以上あるべき");

        // ET 範囲の許容上下限（おおよそ 1849–2150 を J2000 ET 秒で覆う緩い境界）。
        const ET_LO: f64 = -4.8e9;
        const ET_HI: f64 = 4.8e9;

        let mut has_sun = false; // 10
        let mut has_moon = false; // 301
        let mut has_earth = false; // 399
        let mut has_emb = false; // 3
        for s in &segs {
            assert_eq!(s.data_type, 2, "全セグメント data_type==2 のはず");
            assert_eq!(s.frame, 1, "全セグメント frame==1（J2000/ICRF）のはず");
            assert!(s.start_et < s.end_et, "start_et<end_et のはず: {s:?}");
            assert!(
                s.start_et >= ET_LO && s.end_et <= ET_HI,
                "ET 範囲内のはず: {s:?}"
            );
            match s.target {
                10 => has_sun = true,
                301 => has_moon = true,
                399 => has_earth = true,
                3 => has_emb = true,
                _ => {}
            }
        }
        assert!(has_sun, "target=10(Sun) を含むはず");
        assert!(has_moon, "target=301(Moon) を含むはず");
        assert!(has_earth, "target=399(Earth) を含むはず");
        assert!(has_emb, "target=3(EMB) を含むはず");
    }

    // ==================================================================
    // 境界・堅牢性（ミューテーション生存変異の検出用）
    // ==================================================================

    /// 観点: NEXT が自己参照（範囲内だが巡回）でも無限ループせず MalformedSpk を返す（巡回ガード）。
    #[test]
    fn rejects_cyclic_summary_chain() {
        let file_rec = build_file_record(2, 2, 0); // FWARD=2
        let rec2 = build_summary_record(2.0, 0.0, &[sample_sum()]); // NEXT=2（自己参照・範囲内）
        let mut bytes = file_rec;
        bytes.extend_from_slice(&rec2);

        let err = parse_spk_segments(&bytes).unwrap_err();
        assert!(matches!(err, EphemerisError::MalformedSpk(_)));
    }

    /// 観点: NSUM=26（容量超過）が隣接レコードに食い込み読取り境界に達しなくても容量チェックで拒否する。
    #[test]
    fn rejects_nsum_exceeds_capacity_into_adjacent_record() {
        let file_rec = build_file_record(2, 2, 0);
        let mut rec2 = build_summary_record(0.0, 0.0, &[]);
        put_f64(&mut rec2, 16, 26.0); // NSUM=26（容量 25 超過）
        let rec3 = vec![0u8; RECORD_BYTES]; // 隣接パディング（26 件分の読取りが総バイト内に収まる）
        let mut bytes = file_rec;
        bytes.extend_from_slice(&rec2);
        bytes.extend_from_slice(&rec3);

        let err = parse_spk_segments(&bytes).unwrap_err();
        assert!(matches!(err, EphemerisError::MalformedSpk(_)));
    }

    /// 観点: NSUM=25（容量ちょうど）の境界正常系。25 件が解析され先頭・末尾が一致すること。
    #[test]
    fn parses_full_capacity_record() {
        let file_rec = build_file_record(2, 2, 0);
        let sums: Vec<Sum> = (0..25)
            .map(|i| Sum {
                target: i,
                ..sample_sum()
            })
            .collect();
        let rec2 = build_summary_record(0.0, 0.0, &sums);
        let mut bytes = file_rec;
        bytes.extend_from_slice(&rec2);

        let segs = parse_spk_segments(&bytes).expect("容量ちょうど 25 件は成功すべき");
        assert_eq!(segs.len(), 25);
        assert_seg_eq(&segs[0], &sums[0]);
        assert_seg_eq(&segs[24], &sums[24]);
    }
}
