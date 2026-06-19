//! 同梱フィクスチャの決定的チェックサム（変更検知・固定回帰の固定性, architecture §11）。
//!
//! アルゴリズムは SHA-256 に固定し 16 進小文字で表現する（xtask `sha256_hex` と同契約）。
//! [`FIXTURES_CHECKSUM`] は記録済みの期待値で、フィクスチャを意図的に変更したときのみ更新する
//! （意図しないドリフトをテストで検出する）。

use sha2::{Digest, Sha256};

use crate::FIXTURE_FILES;

/// 同梱フィクスチャ（[`crate::FIXTURE_FILES`] を固定順に連結した生バイト列）の SHA-256 を
/// 16 進小文字（64 文字）で返す。
pub fn fixtures_checksum() -> String {
    let mut hasher = Sha256::new();
    for file in FIXTURE_FILES {
        // 改行コード非依存にする（CRLF/LF どちらでチェックアウトされても同一ハッシュ）。
        // `include_str!` は作業ツリーのバイト列をそのまま埋め込むため、core.autocrlf 等で
        // CRLF 化された環境でもドリフト検出が誤発火しないよう CR を除去してから取り込む。
        // フィクスチャの「内容」を見るのが目的で、行末様式の差は変更とみなさない。
        let normalized: Vec<u8> = file.bytes().filter(|&b| b != b'\r').collect();
        // 各ファイルに長さプレフィックス（u64 LE）を前置し、連結境界の曖昧性
        // （ファイル A 末尾とファイル B 先頭の再分割）を排して衝突を防ぐ。
        hasher.update((normalized.len() as u64).to_le_bytes());
        hasher.update(&normalized);
    }
    let digest = hasher.finalize();
    let mut hex = String::with_capacity(digest.len() * 2);
    for byte in digest {
        hex.push(char::from_digit((byte >> 4) as u32, 16).expect("nibble is 0..=15"));
        hex.push(char::from_digit((byte & 0x0f) as u32, 16).expect("nibble is 0..=15"));
    }
    hex
}

/// 記録済みの期待チェックサム（フィクスチャ変更時のみ更新する＝ドリフト検出器）。
pub const FIXTURES_CHECKSUM: &str =
    "f74cf7b562b2ace4f16063f6a5c535918a1d393364db20d7fbb8fe0ad642c72e";
