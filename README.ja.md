# umbra-rs

**設計・実装・検証のすべてを AI エージェントが行う、実験的な純 Rust 製 日食予報エンジン。**

[![CI](https://github.com/RyuuNeko1107/umbra-rs/actions/workflows/ci.yml/badge.svg)](https://github.com/RyuuNeko1107/umbra-rs/actions/workflows/ci.yml)

*(English: [README.md](README.md))*

> ⚠️ **実験的・開発中。** これは「AI 駆動のソフトウェア開発」の研究実験です。**実用段階ではなく**、研究・測地・航法・その他の安全上重要な用途には**使用しないでください**。

---

## これは何か

`umbra-rs` は、太陽・月の位置から影円錐幾何・ベッセル要素・全球/局地条件までを、**一貫した規約のもとで純 Rust で一気通貫に**計算することを目指しています。

このリポジトリの特徴は**成果物だけでなく過程**にあります。各工程を AI エージェント（Claude）が、明示された監査可能なワークフローに沿って実行します:

1. **設計優先** — コードを書く*前*に、設計文書一式（規約・誤差バジェット・アルゴリズム・数値方針・物理モデル・運用）を作成し相互検証する。
2. **独立した反証レビュー** — 別個の AI レビュアーが、設計を複数観点（天体力学・数値解析・ソフトウェア設計・検証可能性）から、**作成者の結論を渡されずに**監査し、破綻を探す。
3. **テスト駆動の実装** — テストに対してコードを実装し、検証は**すべて Docker 内**で行う。

これらすべてをリポジトリに残し、実験を再現・検査可能にしています。

## 精度の方針: *検証可能な*精度

本プロジェクトが約束するのは **検証可能な精度** のみです。すなわち、規定したモデル仮定（平均月縁・点光源の太陽・指定した ΔT/天体暦/半径モデル）のもとで、検証オラクルと一致する範囲です。Standard プロファイルの**目標**（モデル内・オラクル比較）:

| 量 | 目標 |
|---|---|
| 最大食時刻（TT 基準） | ±1〜2 秒 |
| 局地接触時刻（幾何分） | ±2 秒 |
| 食分 | ±0.0005 |
| 中心線位置 | sub-km |

**保証しない**もの: 現実に観測される接触時刻（月縁地形により±数秒）、将来日食の UTC 絶対時刻（ΔT/UT1 予測律速）、連続的な EOP 実測が無い年代。目標はあくまで設計目標であり、JPL DE による検証が済むまで**保証値として公開しません**。

詳細な誤差バジェットと設計根拠は [`docs/accuracy.md`](docs/accuracy.md) を参照。

## 設計文書

設計こそがこの実験の本体です。まずはこちらから:

- [`docs/architecture.md`](docs/architecture.md) — クレート構成・型設計・公開境界
- [`docs/conventions.md`](docs/conventions.md) — 単位・座標系・時刻系・符号規約・定数
- [`docs/accuracy.md`](docs/accuracy.md) — 精度プロファイル・誤差バジェット・検証戦略
- [`docs/numerical-policy.md`](docs/numerical-policy.md) — 級数和・微分・求根・多項式フィット
- [`docs/algorithms.md`](docs/algorithms.md)（＋ `docs/algorithms/`） — 数式仕様（手順ごと）
- [`docs/physical-models.md`](docs/physical-models.md) — 大気差・可視性・種別判定
- [`docs/operations.md`](docs/operations.md) — 性能・再現性・feature/MSRV 方針
- [`docs/data-sources.md`](docs/data-sources.md) — 天体暦/EOP の出典とライセンス
- [`docs/reviews/`](docs/reviews/) — 独立設計監査
- [`docs/issues/`](docs/issues/) — 責務単位の実装チケット（001〜047）

## 状況

- ✅ Milestone 0 — 設計完了・独立監査済み
- 🚧 Milestone 1 — 数学・時刻基盤（`umbra-core`）実装中

## ワークスペース

| クレート | 役割 |
|---|---|
| `umbra-core` | 時刻・角度・距離・ベクトル・定数・数値解法 |
| `umbra-ephemeris` | 太陽/月の天体暦と見かけ位置補正（WIP） |
| `umbra-eclipse` | 影幾何・ベッセル要素・全球/局地条件（WIP） |
| `umbra-geo` | 中心線・限界線・GeoJSON（WIP） |
| `umbra-cli` | コマンドラインインターフェース（WIP） |
| `umbra-fixtures` | 検証フィクスチャと許容誤差（テスト専用・WIP） |

## ビルドと検証

実験・検証は**すべて Docker 内**で実行します（Docker 以外のホスト側ツールチェーンは不要）:

```sh
docker compose -p umbra-rs run --rm rust cargo test  --workspace
docker compose -p umbra-rs run --rm rust cargo clippy --workspace --all-targets -- -D warnings
docker compose -p umbra-rs run --rm rust cargo fmt    --all --check
```

## ライセンス（暫定）

これは研究実験であり、**正式なリリース（crates.io への公開等）は予定していません**。全クレートは `publish = false` です。

*コード*は暫定的に [Apache-2.0](LICENSE-APACHE) または [MIT](LICENSE-MIT) のいずれかを選択可能として提供します。ただし**最終的なライセンスは未確定**で、本プロジェクトが依拠する第三者の科学データ（VSOP87 / ELP-MPP02 係数、IERS EOP 等）の再配布条件に依存します（[`docs/data-sources.md`](docs/data-sources.md) §6 を参照、現在確認中）。GPL ライセンスのコード・データは一切取り込んでいません。
