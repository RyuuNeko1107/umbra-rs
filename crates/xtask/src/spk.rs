//! `cargo xtask fetch-de440s` / `verify-de440s` — JPL DE440s SPK カーネルの取得と検証
//! （ISSUE-036・data-sources §2.3/§6）。
//!
//! DE データは **crate 非同梱**（ライセンス・巨大データ回避, data-sources §2.3/§6）。本コマンドは
//! 利用者の**明示手順**として NAIF 公開 SPK を `data/spk/` へ取得し、出典の SHA-256 と照合する
//! （実行時ネットワーク禁止の例外＝この明示 DL 手順のみに隔離, accuracy.md §5）。取得した `.bsp`
//! は git 管理外（`.gitignore`）で、JplEphemeris（feature `jpl`）の検証テストが読む。

use std::path::Path;
use std::process::Command;

use crate::checksum::sha256_hex;
use crate::error::XtaskError;

/// NAIF 公開 DE440s SPK（短期版・ET 1849–2150・v0.1 範囲 1900–2100 を被覆）。
const DE440S_URL: &str =
    "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de440s.bsp";
/// 取得先（workspace ルート相対。`cargo xtask` の cwd 前提・eop.rs と同方式）。
const DE440S_PATH: &str = "data/spk/de440s.bsp";
/// 出典バイト数（Last-Modified 2020-12-21・固定配布物。data/spk/PROVENANCE.md）。
const DE440S_SIZE: u64 = 32_726_016;
/// 出典 SHA-256（取得日 2026-06-20。data/spk/PROVENANCE.md・改竄/破損検出）。
const DE440S_SHA256: &str = "c1c7feeab882263fc493a9d5a5b2ddd71b54826cdf65d8d17a76126b260a49f2";

fn io_err(path: &str, source: std::io::Error) -> XtaskError {
    XtaskError::Io {
        path: path.to_string(),
        source,
    }
}

/// `data/spk/de440s.bsp` を NAIF から取得（curl）し、サイズと SHA-256 を出典値と照合する。
pub fn fetch_de440s() -> Result<(), XtaskError> {
    if let Some(parent) = Path::new(DE440S_PATH).parent() {
        std::fs::create_dir_all(parent).map_err(|e| io_err(DE440S_PATH, e))?;
    }
    // 明示 DL（curl・実行時ネットワーク禁止の唯一の例外手順, accuracy.md §5）。
    let status = Command::new("curl")
        .args(["-fsSL", "--max-time", "600", DE440S_URL, "-o", DE440S_PATH])
        .status()
        .map_err(|e| io_err(DE440S_PATH, e))?;
    if !status.success() {
        return Err(io_err(
            DE440S_PATH,
            std::io::Error::other(format!("curl failed ({status}) for {DE440S_URL}")),
        ));
    }
    verify_de440s()?;
    println!("fetched {DE440S_PATH} (verified sha256 {DE440S_SHA256})");
    Ok(())
}

/// 既存の `data/spk/de440s.bsp` のサイズと SHA-256 を出典値と照合する（DL 不要の整合検査）。
/// 不在/破損/改竄は [`XtaskError`]。JplEphemeris テスト前提の正当性ゲート。
pub fn verify_de440s() -> Result<(), XtaskError> {
    let bytes = std::fs::read(DE440S_PATH).map_err(|e| io_err(DE440S_PATH, e))?;
    if bytes.len() as u64 != DE440S_SIZE {
        return Err(XtaskError::ChecksumMismatch {
            dataset: "de440s.bsp (size bytes)".to_string(),
            stored: DE440S_SIZE.to_string(),
            regenerated: bytes.len().to_string(),
        });
    }
    let actual = sha256_hex(&bytes);
    if actual != DE440S_SHA256 {
        return Err(XtaskError::ChecksumMismatch {
            dataset: "de440s.bsp".to_string(),
            stored: DE440S_SHA256.to_string(),
            regenerated: actual,
        });
    }
    println!("verified {DE440S_PATH} (sha256 {DE440S_SHA256})");
    Ok(())
}
