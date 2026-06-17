//! `umbra-eclipse` — 日食エンジン（影幾何・ベッセル要素・全球/局地条件）。
//!
//! Milestone 4 以降。現状は影円錐幾何（ISSUE-019）と最小エラー型。基本面・ベッセル要素・
//! 全球/局地は順次追加（`docs/algorithms/` §5–§9）。
//! 幾何ロジックは [`umbra_ephemeris::MockEphemeris`] 上で検証する（accuracy.md §3.1）。

pub mod bessel_poly;
pub mod besselian;
pub mod calc_metadata;
pub mod candidates;
pub mod config;
pub mod conjunction;
pub mod eclipse_filter;
pub mod error;
pub mod fundamental;
pub mod global;
pub mod horizontal;
pub mod local_contacts;
pub mod local_maximum;
pub mod magnitude;
pub mod path;
pub mod polynomial;
pub mod projection;
pub mod results;
pub mod shadow;
pub mod source;

pub use bessel_poly::{BesselFitError, BesselianPolynomial};
pub use besselian::{
    besselian_elements, besselian_elements_at, besselian_mu, BesselianElements,
    InstantaneousBesselianElements,
};
pub use calc_metadata::CalculationMetadata;
pub use config::{AccuracyProfile, EngineConfig, LunarRadiusModel, SolarRadiusModel};
pub use error::EclipseError;
pub use fundamental::{fundamental_plane_basis, FundamentalPlaneBasis};
pub use global::{classify, SolarEclipseKind};
pub use horizontal::{Horizontal, RefractionModel, Visibility};
pub use magnitude::{EclipseMagnitude, Obscuration};
pub use path::{EclipsePath, PathOptions, PathSample};
pub use polynomial::Polynomial;
pub use projection::{project_observer_to_fundamental, ObserverFundamental};
pub use results::{
    GlobalCircumstances, GlobalContact, GreatestEclipse, LocalCircumstances, LocalContact,
    LocalContactSet, SolarEclipse, VisibleSolarEclipse,
};
pub use shadow::{shadow_cone, ShadowCone};
pub use source::{BesselianSource, DirectBesselianSource, InstantaneousEvaluator};
