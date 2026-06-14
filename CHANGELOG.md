# Changelog

All notable changes to this project are documented here.
The format is based on [Keep a Changelog](https://keepachangelog.com/), and this
project adheres to [Semantic Versioning](https://semver.org/).

## [Unreleased]

### Added
- Milestone 0: complete design documentation (`docs/`) — conventions, accuracy/error
  budget, architecture, public API draft, numerical policy, physical models,
  operations, full algorithm specification, and data-source/licensing notes.
- Independent multi-perspective design audits (`docs/reviews/`).
- Per-responsibility implementation tickets `docs/issues/ISSUE-001`–`047` with an index.
- Milestone 1 (in progress): `umbra-core` foundations — physical/convention constants,
  angle newtypes (`Radians`/`Degrees`), two-part Julian date (`JulianDate2`), `Vector3`.
- Docker-based verification setup (`Dockerfile`, `docker-compose.yml`).

### Notes
- Pre-1.0; public API and numeric outputs may change. Precision figures are design
  targets ("verifiable precision", model-internal) and are not published guarantees
  until validated against JPL DE (see `docs/accuracy.md`).
