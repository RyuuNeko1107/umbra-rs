//! JPL DE（DE440s SPK）バックエンド＝Reference オラクル（ISSUE-036・feature `jpl`）。
//!
//! 自前 DAF/SPK パーサ（最小版・太陽/月/地球/EMB）。差分テストの第一義オラクル（accuracy.md §3.1）。
//! DE データは crate 非同梱（利用者が `.bsp` を任意取得, data-sources §2.3）。
//!
//! - [`daf`]: DAF/SPK の構造解析（ファイルレコード＋セグメント記述子）。S1。
//! - [`eval`]: SPK type 2 Chebyshev 評価（セグメント native 位置/速度）。S2。
//! - 後続: `JplEphemeris`＋`Ephemeris` 実装（S3）。

pub(crate) mod daf;
pub(crate) mod eval;
