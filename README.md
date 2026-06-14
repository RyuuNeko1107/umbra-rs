# umbra-rs

**An experimental, pure-Rust solar eclipse prediction engine — where the design, implementation, and verification are all carried out by an AI agent.**

*(日本語: [README.ja.md](README.ja.md))*

> ⚠️ **Experimental / work in progress.** This is a research experiment in AI-driven software engineering. It is **not** ready for use and must **not** be used for research, surveying, navigation, or any safety-critical purpose.

---

## What this is

`umbra-rs` aims to compute solar eclipses end-to-end in pure Rust — from Sun/Moon positions through shadow-cone geometry, Besselian elements, and global/local circumstances — under one consistent set of conventions.

The distinguishing feature of this repository is **the process, not just the product**: every stage is performed by an AI agent (Claude) following an explicit, auditable workflow:

1. **Design first** — a full set of design documents (conventions, error budget, algorithms, numerical policy, physical models, operations) is written and cross-checked *before* any code.
2. **Independent adversarial review** — separate AI reviewers audit the design from multiple angles (celestial mechanics, numerical analysis, software architecture, verifiability) *without* being shown the author's conclusions, and try to break it.
3. **Test-driven implementation** — code is built against tests, verified entirely inside Docker.

All of this lives in the repo so the experiment is reproducible and inspectable.

## Precision philosophy: *verifiable* precision

This project commits only to **verifiable precision** — agreement with a validation oracle, under stated model assumptions (mean lunar limb, point-source Sun, specified ΔT/ephemeris/radius models). Concretely, the Standard profile *targets* (model-internal, oracle-compared):

| Quantity | Target |
|---|---|
| Greatest-eclipse time (TT) | ±1–2 s |
| Local contact time (geometric) | ±2 s |
| Magnitude | ±0.0005 |
| Central line position | sub-km |

Explicitly **not guaranteed**: real-world observed contact times (limited by lunar-limb topography, ±several seconds), absolute UTC of future eclipses (limited by ΔT/UT1 prediction), and dates before continuous EOP measurements. Targets are design goals and are **not published as guarantees** until validated against JPL DE.

See [`docs/accuracy.md`](docs/accuracy.md) for the full error budget and the design rationale.

## Design documents

The design is the heart of this experiment. Start here:

- [`docs/architecture.md`](docs/architecture.md) — crate layout, type design, public boundaries
- [`docs/conventions.md`](docs/conventions.md) — units, frames, time scales, sign conventions, constants
- [`docs/accuracy.md`](docs/accuracy.md) — precision profiles, error budget, verification strategy
- [`docs/numerical-policy.md`](docs/numerical-policy.md) — summation, differentiation, root-finding, polynomial fitting
- [`docs/algorithms.md`](docs/algorithms.md) (+ `docs/algorithms/`) — the mathematical specification, step by step
- [`docs/physical-models.md`](docs/physical-models.md) — refraction, visibility, eclipse-type classification
- [`docs/operations.md`](docs/operations.md) — performance, reproducibility, feature/MSRV policy
- [`docs/data-sources.md`](docs/data-sources.md) — ephemeris/EOP sources and licensing
- [`docs/reviews/`](docs/reviews/) — independent design audits
- [`docs/issues/`](docs/issues/) — per-responsibility implementation tickets (001–047)

## Status

- ✅ Milestone 0 — design complete and independently audited
- 🚧 Milestone 1 — math & time foundations (`umbra-core`) in progress

## Workspace

| crate | role |
|---|---|
| `umbra-core` | time, angle, distance, vector, constants, numerics |
| `umbra-ephemeris` | Sun/Moon ephemeris and apparent-position corrections (WIP) |
| `umbra-eclipse` | shadow geometry, Besselian elements, global/local circumstances (WIP) |
| `umbra-geo` | central line / limits / GeoJSON (WIP) |
| `umbra-cli` | command-line interface (WIP) |
| `umbra-fixtures` | validation fixtures and tolerances (test-only, WIP) |

## Building & testing

All experiments and verification run **inside Docker** (no host toolchain needed beyond Docker):

```sh
docker compose -p umbra-rs run --rm rust cargo test  --workspace
docker compose -p umbra-rs run --rm rust cargo clippy --workspace --all-targets -- -D warnings
docker compose -p umbra-rs run --rm rust cargo fmt    --all --check
```

## License (provisional)

This is a research experiment; **no formal release (e.g. crates.io publish) is planned**, and all crates are marked `publish = false`.

The *code* is provisionally offered under either [Apache-2.0](LICENSE-APACHE) or [MIT](LICENSE-MIT) at your option. However, the **final licensing is not settled** — it depends on the terms of the third-party scientific data this project relies on (VSOP87 / ELP-MPP02 coefficients, IERS EOP, etc.), whose redistribution terms are still being clarified (see [`docs/data-sources.md`](docs/data-sources.md) §6). No GPL-licensed code or data is vendored in.
