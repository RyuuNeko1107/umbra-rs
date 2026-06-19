//! `umbra-fixtures` — 検証用フィクスチャと許容誤差プロファイル（テスト専用）。
//!
//! 外部オラクル（NASA 5MCSE / USNO Solar Eclipse Calculator）の**数値事実のみ転記**した
//! 固定回帰データ（ゴールデン20）を保持する（`docs/issues/ISSUE-029`）。実装出力との一致比較は
//! ISSUE-030（`ToleranceProfile`）の領分で、本 crate は**データ整備とローダ**に徹する。
//! 通常の依存先には含めない（検証専用・`docs/architecture.md` §1）。
//!
//! ## ゴールデン20（完備）
//! v0.1 目標の **20 日食**（皆既7/金環7/ハイブリッド2/部分4・各 5 地点・計 100 地点・accuracy §3.4）を
//! 実在オラクルから**完全転記**して同梱する。全球パラメータは NASA 5MCSE（besselian element ページ）、
//! 地点別局地状況は USNO Solar Eclipse Calculator API（1800–2050 対応）から数値事実のみ転記。
//! 年代は 1999–2026、地理は 6 大陸＋極域、食条件は中心線/限界近傍/部分域/日の出/日没/高標高を被覆する。
//! ローダ正本は [`golden_eclipses`]（同梱の全 20 件を返す）。
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
    aggregate_global, aggregate_local, compare_global, compare_local, render_json, render_text,
    report_against_golden, ErrorStats, GlobalErrors, GlobalReport, GoldenComputer, GoldenReport,
    LocalErrors, LocalReport, ToleranceProfile,
};
pub use types::{GoldenContact, GoldenEclipse, GoldenLocation, LocationClass, OracleSource};

/// 同梱フィクスチャ（生 TOML）の埋め込み。`include_str!` でビルド時に取り込み、ローダとチェックサムが
/// **同一バイト列・同一順序**で参照する（チェックサムの決定性のため順序固定＝event_key 昇順）。
pub(crate) const FIXTURE_FILES: [&str; 20] = [
    include_str!("../data/golden/1999-08-11-total.toml"),
    include_str!("../data/golden/2002-12-04-total.toml"),
    include_str!("../data/golden/2003-05-31-annular.toml"),
    include_str!("../data/golden/2005-10-03-annular.toml"),
    include_str!("../data/golden/2006-03-29-total.toml"),
    include_str!("../data/golden/2009-07-22-total.toml"),
    include_str!("../data/golden/2010-01-15-annular.toml"),
    include_str!("../data/golden/2012-05-20-annular.toml"),
    include_str!("../data/golden/2013-11-03-hybrid.toml"),
    include_str!("../data/golden/2014-10-23-partial.toml"),
    include_str!("../data/golden/2017-08-21-total.toml"),
    include_str!("../data/golden/2018-08-11-partial.toml"),
    include_str!("../data/golden/2020-06-21-annular.toml"),
    include_str!("../data/golden/2021-06-10-annular.toml"),
    include_str!("../data/golden/2022-10-25-partial.toml"),
    include_str!("../data/golden/2023-04-20-hybrid.toml"),
    include_str!("../data/golden/2023-10-14-annular.toml"),
    include_str!("../data/golden/2024-04-08-total.toml"),
    include_str!("../data/golden/2025-03-29-partial.toml"),
    include_str!("../data/golden/2026-08-12-total.toml"),
];
