//! `umbra-geo` — 中心線・限界線・部分食域・GeoJSON 出力。
//!
//! Milestone 9 で実装（v0.1 では path/GeoJSON は未実装、`docs/issues/ISSUE-045`）。
//! 地理座標の幾何プリミティブは [`geometry`] を参照（`docs/api-draft.md` §4）。

pub mod geometry;

pub use geometry::{GeoLine, GeoPoint, GeoPolygon};
