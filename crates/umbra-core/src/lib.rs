//! `umbra-core` — `umbra-rs` の基盤プリミティブ。
//!
//! 時刻・角度・距離・ベクトル・物理定数・数値解法といった、天文暦や日食固有の
//! 概念に依存しない純粋な土台のみを提供する（`docs/architecture.md` §1）。
//!
//! 規約は `docs/conventions.md`、数値方針は `docs/numerical-policy.md` を参照。

pub mod angle;
pub mod calendar;
pub mod constants;
pub mod deltat;
pub mod ellipsoid;
pub mod error;
pub mod julian;
pub mod matrix;
pub mod metadata;
pub mod solver;
pub mod time;
pub mod vector;

pub use angle::{Degrees, Radians};
pub use calendar::{gregorian_to_jd2, jd2_to_gregorian};
pub use deltat::{tt_to_ut1, ut1_to_tt, DeltaTModel, EspenakMeeusDeltaT};
pub use ellipsoid::{Ellipsoid, GeocentricObserver};
pub use error::{DomainError, SolverError, TimeError};
pub use julian::JulianDate2;
pub use matrix::Matrix3;
pub use metadata::DataSetMetadata;
pub use solver::{brent_root, minimize_golden};
pub use time::{
    TaiInstant, TdbInstant, TimeInterval, TimeRange, TtInstant, Ut1Instant, UtcInstant,
};
pub use vector::{UnitVector3, Vector3};
