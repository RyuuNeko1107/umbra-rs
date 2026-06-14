//! `umbra-ephemeris` — 天体暦バックエンドと見かけ位置補正。
//!
//! Milestone 2 以降で実装（VSOP87D / ELP/MPP02 / JPL DE）。設計は `docs/architecture.md` §4、
//! `docs/algorithms/03-ephemeris.md`、`docs/algorithms/02-frames.md` を参照。
#![allow(unused_imports)]

use umbra_core as _;
