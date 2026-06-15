//! 生成物バイト列の決定的ハッシュ（ISSUE-046・`docs/data-sources.md` §5）。
//!
//! アルゴリズムは SHA-256 に固定し、`DataSetMetadata.checksum` に 16 進小文字で記録する。

use sha2::{Digest, Sha256};

/// バイト列の SHA-256 を 16 進小文字文字列で返す（決定的）。
pub fn sha256_hex(bytes: &[u8]) -> String {
    let digest = Sha256::digest(bytes);
    let mut hex = String::with_capacity(digest.len() * 2);
    for byte in digest {
        // 16 進小文字 2 桁ずつ。`DataSetMetadata.checksum` の表現契約。
        hex.push(char::from_digit((byte >> 4) as u32, 16).expect("nibble is 0..=15"));
        hex.push(char::from_digit((byte & 0x0f) as u32, 16).expect("nibble is 0..=15"));
    }
    hex
}

#[cfg(test)]
mod tests {
    use super::*;

    // ------------------------------------------------------------------
    // SHA-256 既知応答（Known-Answer Test）の一次/権威ある出典:
    //   - "abc" の digest は NIST FIPS 180-4 "Secure Hash Standard"
    //     Appendix B.1（SHA-256 のワークド・エグザンプル）に逐語で示される値。
    //     https://nvlpubs.nist.gov/nistpubs/FIPS/NIST.FIPS.180-4.pdf
    //   - 空文字列 "" の digest は同規格の定義から導かれ、広く公開された
    //     標準テストベクトル集（di-mgt.com.au/sha_testvectors.html、NESSIE 等）
    //     と一致する。
    // 転記値（16 進小文字・連結形）:
    //   SHA-256("")    = e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855
    //   SHA-256("abc") = ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad
    // ------------------------------------------------------------------
    const SHA256_EMPTY: &str = "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855";
    const SHA256_ABC: &str = "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad";

    /// FIPS 180-4 既知応答（empty / "abc"）に一致する。
    #[test]
    fn matches_fips180_4_known_answers() {
        assert_eq!(sha256_hex(b""), SHA256_EMPTY);
        assert_eq!(sha256_hex(b"abc"), SHA256_ABC);
    }

    /// 出力は 64 文字の 16 進小文字（[0-9a-f]）。checksum 欄の表現契約。
    #[test]
    fn output_is_64_lowercase_hex_chars() {
        let h = sha256_hex(b"umbra-rs dataset bytes");
        assert_eq!(h.len(), 64, "digest hex length");
        assert!(
            h.chars()
                .all(|c| c.is_ascii_digit() || ('a'..='f').contains(&c)),
            "digest must be lowercase hex: {h}"
        );
    }

    /// 決定性: 同一入力を 2 回ハッシュすると同一出力。
    #[test]
    fn deterministic_for_same_input() {
        let bytes = b"deterministic input vector";
        assert_eq!(sha256_hex(bytes), sha256_hex(bytes));
    }

    /// 入力 1 バイトの変化で出力が変わる（衝突回避・改変検出の最小性質）。
    #[test]
    fn single_byte_change_changes_digest() {
        let a = sha256_hex(&[0u8, 1, 2, 3]);
        let b = sha256_hex(&[0u8, 1, 2, 4]); // 末尾 1 バイトのみ差
        assert_ne!(a, b, "1-byte change must change the digest");
    }
}
