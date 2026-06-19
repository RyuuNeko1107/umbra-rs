//! xtask のエラー型（ISSUE-046）。

use thiserror::Error;
use umbra_eclipse::EclipseError;

/// xtask サブコマンド実行中のエラー。
#[derive(Debug, Error)]
pub enum XtaskError {
    /// 未知のサブコマンド。
    #[error("unknown subcommand: {0}")]
    UnknownSubcommand(String),
    /// 未知のデータセット指定。
    #[error("unknown dataset: {0}")]
    UnknownDataset(String),
    /// `--dataset` の引数が欠落。
    #[error("missing value for {0}")]
    MissingArgument(String),
    /// 当該データセットの生成ロジックは別 Issue で実装予定（033/034/040）。
    #[error("not yet implemented: {0}")]
    NotImplemented(String),
    /// 再生成物の checksum がコミット済み generated と一致しない。
    #[error("checksum mismatch for {dataset}: stored {stored}, regenerated {regenerated}")]
    ChecksumMismatch {
        /// 対象データセット。
        dataset: String,
        /// 記録済み checksum。
        stored: String,
        /// 再生成物の checksum。
        regenerated: String,
    },
    /// packed バイト列の長さが f64 境界（8 の倍数）でない等の不整合。
    #[error("malformed packed data: {0}")]
    MalformedPacked(String),
    /// 一次原データ（IERS 章動表等）のパース失敗・項数不整合。
    #[error("malformed source data: {0}")]
    MalformedSource(String),
    /// 原データファイルの入出力エラー。
    #[error("io error reading {path}: {source}")]
    Io {
        /// 対象パス。
        path: String,
        /// 元の I/O エラー。
        #[source]
        source: std::io::Error,
    },
    /// フラグに不正な値（`--format`/`--accuracy` の未知値, ISSUE-030 S30f validate）。
    #[error("invalid value '{value}' for {flag}")]
    InvalidArgument {
        /// 対象フラグ（例 `--format`）。
        flag: String,
        /// 与えられた不正値。
        value: String,
    },
    /// 日食エンジン側エラー（validate のゴールデン照合・透過, ISSUE-030 S30f）。
    #[error(transparent)]
    Eclipse(#[from] EclipseError),
    /// レポート JSON 整形失敗（serde_json 由来・透過, ISSUE-030 S30f）。
    #[error(transparent)]
    Json(#[from] serde_json::Error),
}
