# 検証用 Rust ツールチェーン（rustfmt / clippy を同梱）。
# 実験・検証は全て本イメージ上の Docker 内で実行する（docs 方針）。
FROM rust:1.96-bookworm

RUN rustup component add rustfmt clippy

WORKDIR /work
