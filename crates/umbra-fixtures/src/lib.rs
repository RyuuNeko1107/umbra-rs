//! `umbra-fixtures` — 検証用フィクスチャと許容誤差プロファイル（テスト専用）。
//!
//! 外部オラクル（NASA 5MCSE / USNO Solar Eclipse Calculator）の**数値事実のみ転記**した
//! 固定回帰データ（ゴールデン20）を保持する（`docs/issues/ISSUE-029`）。実装出力との一致比較は
//! ISSUE-030（`ToleranceProfile`）の領分で、本 crate は**データ整備とローダ**に徹する。
//! 通常の依存先には含めない（検証専用・`docs/architecture.md` §1）。
//!
//! ## seed フェーズ
//! v0.1 最終目標は 20 日食（皆既5/金環5/部分3/ハイブリッド2/境界5・accuracy §3.4）。本スライスは
//! **インフラ＋実データ seed**（型・TOML スキーマ・ローダ・チェックサム・被覆メタテストを、実在オラクル
//! から完全転記した少数日食で確立）であり、20 件への拡張は後続。ローダ正本は [`golden_eclipses`]
//! （現在同梱されている全件を返す。20 件到達時に "golden twenty" となる）。
//!
//! ## ハードコード禁止の規律（conventions §11 / data-sources §0/§4）
//! 数値は外部オラクルの**転記**であり、実装側へコピーしない。各 [`OracleSource`] に出典・取得日・
//! k/ΔT 慣習・ライセンス注記を必須付与する。

mod checksum;
mod loader;
mod report;
mod types;

pub use checksum::{fixtures_checksum, FIXTURES_CHECKSUM};
pub use loader::golden_eclipses;
pub use report::{
    aggregate_global, aggregate_local, compare_global, compare_local, ErrorStats, GlobalErrors,
    GlobalReport, LocalErrors, LocalReport, ToleranceProfile,
};
pub use types::{GoldenContact, GoldenEclipse, GoldenLocation, LocationClass, OracleSource};

/// 同梱フィクスチャ（生 TOML）の埋め込み。`include_str!` でビルド時に取り込み、ローダとチェックサムが
/// **同一バイト列・同一順序**で参照する（チェックサムの決定性のため順序固定）。
pub(crate) const FIXTURE_FILES: [&str; 5] = [
    include_str!("../data/golden/2017-08-21-total.toml"),
    include_str!("../data/golden/2023-04-20-hybrid.toml"),
    include_str!("../data/golden/2023-10-14-annular.toml"),
    include_str!("../data/golden/2024-04-08-total.toml"),
    include_str!("../data/golden/2025-03-29-partial.toml"),
];
