//! packed 形式の f64 直列化（ISSUE-046）。
//!
//! 係数テーブルは little-endian f64 の連続配列として packed 化する（生成時に単位を確定し、
//! 消費側 ISSUE-013/035 等と byte-for-byte の契約を固定する）。本モジュールは
//! その最小プリミティブ（f64 配列 ⇄ little-endian バイト列）を提供する。

use crate::error::XtaskError;

/// f64 スライスを little-endian バイト列へ直列化する（決定的）。
pub fn pack_f64_le(values: &[f64]) -> Vec<u8> {
    let mut bytes = Vec::with_capacity(values.len() * 8);
    for &value in values {
        bytes.extend_from_slice(&value.to_le_bytes());
    }
    bytes
}

/// little-endian バイト列を f64 ベクタへ復元する。
/// 長さが 8 の倍数でなければ [`XtaskError::MalformedPacked`]。
pub fn unpack_f64_le(bytes: &[u8]) -> Result<Vec<f64>, XtaskError> {
    if bytes.len() % 8 != 0 {
        return Err(XtaskError::MalformedPacked(format!(
            "byte length {} is not a multiple of 8 (f64 boundary)",
            bytes.len()
        )));
    }
    let values = bytes
        .chunks_exact(8)
        .map(|chunk| {
            let octet: [u8; 8] = chunk.try_into().expect("chunks_exact(8) yields 8 bytes");
            f64::from_le_bytes(octet)
        })
        .collect();
    Ok(values)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// 非対称・全非ゼロ・特殊値（subnormal/極大/極小）を含む代表配列。
    /// 非対称なので LE バイト順の取り違え（要素内反転）でも往復が崩れて検出できる。
    const REPRESENTATIVE: &[f64] = &[1.0, -2.5, 0.0, f64::MIN_POSITIVE, 1e300, -1e-300];

    /// pack→unpack の往復で元の配列に厳密一致（ビット等価）。
    #[test]
    fn round_trip_preserves_values_exactly() {
        let xs = REPRESENTATIVE;
        let restored = unpack_f64_le(&pack_f64_le(xs)).expect("valid length round-trips");
        assert_eq!(restored.len(), xs.len(), "element count preserved");
        for (i, (&got, &want)) in restored.iter().zip(xs.iter()).enumerate() {
            // bit パターンで厳密一致を要求（±0 や subnormal も区別）。
            assert_eq!(
                got.to_bits(),
                want.to_bits(),
                "element {i}: got {got} want {want}"
            );
        }
    }

    /// NaN を含む配列の往復はビットパターンで保存される（NaN は == では不一致になるため bits 比較）。
    #[test]
    fn round_trip_preserves_nan_bits() {
        let xs = [f64::NAN, 1.0, f64::NEG_INFINITY, f64::INFINITY];
        let restored = unpack_f64_le(&pack_f64_le(&xs)).expect("valid length round-trips");
        for (i, (&got, &want)) in restored.iter().zip(xs.iter()).enumerate() {
            assert_eq!(got.to_bits(), want.to_bits(), "element {i} NaN/inf bits");
        }
    }

    /// pack は決定的: 同一入力 → 同一バイト列、長さ = 8 × 要素数。
    #[test]
    fn pack_is_deterministic_and_8x_length() {
        let xs = REPRESENTATIVE;
        let a = pack_f64_le(xs);
        let b = pack_f64_le(xs);
        assert_eq!(a, b, "pack must be deterministic");
        assert_eq!(a.len(), xs.len() * 8, "packed length = 8 * element count");
    }

    /// 空配列の往復: 空バイト列 → 空 Vec。
    #[test]
    fn empty_round_trips() {
        let packed = pack_f64_le(&[]);
        assert!(packed.is_empty(), "empty input packs to empty bytes");
        let restored = unpack_f64_le(&packed).expect("empty unpacks");
        assert!(restored.is_empty(), "empty bytes unpack to empty vec");
    }

    /// 長さが 8 の倍数でないバイト列は MalformedPacked（f64 境界違反）。
    #[test]
    fn unpack_rejects_non_multiple_of_8() {
        for len in [1usize, 7, 9, 15] {
            let bytes = vec![0u8; len];
            let err = unpack_f64_le(&bytes).expect_err("non-8-multiple must error");
            assert!(
                matches!(err, XtaskError::MalformedPacked(_)),
                "len {len}: expected MalformedPacked, got {err:?}"
            );
        }
    }
}
