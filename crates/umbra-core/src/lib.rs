//! `umbra-core` — `umbra-rs` の基盤プリミティブ。
//!
//! 時刻・角度・距離・ベクトル・物理定数・数値解法といった、天文暦や日食固有の
//! 概念に依存しない純粋な土台のみを提供する（`docs/architecture.md` §1）。
//!
//! 規約は `docs/conventions.md`、数値方針は `docs/numerical-policy.md` を参照。

pub mod angle;
pub mod constants;
pub mod ellipsoid;
pub mod error;
pub mod julian;
pub mod solver;
pub mod vector;

pub use angle::{Degrees, Radians};
pub use ellipsoid::{Ellipsoid, GeocentricObserver};
pub use error::{DomainError, SolverError};
pub use julian::JulianDate2;
pub use solver::{brent_root, minimize_golden};
pub use vector::Vector3;
