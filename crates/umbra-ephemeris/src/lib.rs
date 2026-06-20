//! `umbra-ephemeris` — 天体暦バックエンドと見かけ位置補正。
//!
//! 抽象は [`Ephemeris`] trait（ISSUE-012）。テスト用の人工配置 [`MockEphemeris`]（ISSUE-038）を
//! 提供する。解析暦（VSOP87D+ELP2000-82B = [`AnalyticalEphemeris`]）・JPL DE バックエンドは
//! Milestone 2 以降。
//! 設計は `docs/architecture.md` §4、`docs/algorithms/03-ephemeris.md` / `02-frames.md`。

pub mod analytical;
pub mod apparent;
pub mod cio;
pub mod eop;
pub mod ephemeris;
pub mod frames;
#[cfg(feature = "jpl")]
pub mod jpl;
pub mod mock;
pub mod moon;
pub mod nutation;
pub mod sun;
pub mod time_data;

pub use analytical::AnalyticalEphemeris;
pub use apparent::{apparent_cirs, AstrometryOptions};
#[cfg(feature = "bundled-data")]
pub use eop::bundled_eop;
pub use ephemeris::{
    Body, Ephemeris, EphemerisError, EphemerisFrame, EphemerisMetadata, Origin, StateVector,
};
#[cfg(feature = "jpl")]
pub use jpl::backend::JplEphemeris;
pub use mock::MockEphemeris;
#[cfg(feature = "bundled-data")]
pub use time_data::bundled_time_data;
pub use time_data::{time_data_from_path, TimeDataError};
