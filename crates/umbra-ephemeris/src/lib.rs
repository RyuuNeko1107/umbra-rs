//! `umbra-ephemeris` — 天体暦バックエンドと見かけ位置補正。
//!
//! 抽象は [`Ephemeris`] trait（ISSUE-012）。テスト用の人工配置 [`MockEphemeris`]（ISSUE-038）を
//! 提供する。解析暦（VSOP87D+ELP/MPP02）・JPL DE バックエンドは Milestone 2 以降。
//! 設計は `docs/architecture.md` §4、`docs/algorithms/03-ephemeris.md` / `02-frames.md`。

pub mod cio;
pub mod ephemeris;
pub mod frames;
pub mod mock;
pub mod nutation;

pub use ephemeris::{
    Body, Ephemeris, EphemerisError, EphemerisFrame, EphemerisMetadata, Origin, StateVector,
};
pub use mock::MockEphemeris;
