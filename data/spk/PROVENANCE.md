# JPL DE SPK kernels — provenance（ISSUE-036 / data-sources §2.3, §4, §6）

このディレクトリは JPL 開発暦（DE）の SPK バイナリ（`.bsp`）を置く場所です。
**バイナリは crate に同梱せず**（ライセンス・巨大データ回避, data-sources §2.3/§6）、`.gitignore`
で git 管理外（`/data/spk/*.bsp`）。利用者が下記の明示手順で取得します（実行時ネットワーク禁止＝
取得は xtask の明示コマンドのみに隔離, accuracy.md §5）。

## 取得手順

```
cargo xtask fetch-de440s     # NAIF から取得し SHA-256 照合
cargo xtask verify-de440s    # 取得済みファイルの SHA-256 整合のみ検査（DL 不要）
```

（Docker 検証環境では `docker compose -p umbra-rs run --rm rust cargo xtask fetch-de440s`）

## de440s.bsp

| 項目 | 値 |
|---|---|
| ファイル | `de440s.bsp`（DE440 短期版 SPK） |
| 出典 URL | https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de440s.bsp |
| サイズ | 32,726,016 bytes |
| SHA-256 | `c1c7feeab882263fc493a9d5a5b2ddd71b54826cdf65d8d17a76126b260a49f2` |
| Last-Modified | 2020-12-21（NAIF 配布・固定） |
| 取得日 | 2026-06-20 |
| ET 範囲 | 1849-12-26 〜 2150-01-21（v0.1 の 1900–2100 を被覆） |
| 形式 | DAF/SPK・LTL-IEEE（little-endian）・全セグメント type 2（Chebyshev 位置） |

日食用に参照する body（NAIF ID）: Sun(10)/Moon(301)/Earth(399)/EMB(3)/SSB(0)。
地心太陽 = 10−(3+399)、地心月 = 301−399（SSB 基準ベクトルの差で Geocenter 原点へ）。

## 出典・ライセンス（引用必須・data-sources §2.3/§6）

- **DE440/441**: Park, R. S., Folkner, W. M., Williams, J. G., & Boggs, D. H. (2021).
  *The JPL Planetary and Lunar Ephemerides DE440 and DE441.* The Astronomical Journal, 161:105.
  DOI: 10.3847/1538-3881/abd414。
- 提供: NASA/JPL/Caltech（NAIF）。**米政府/Caltech-JPL の著作物**。OSS 同梱可否が未確定のため
  **本リポジトリには同梱しない**（利用者が任意取得）。`JplEphemeris::metadata` の license に
  「JPL/Caltech・非同梱・任意DL」を明記する。
- 用途: **Reference オラクル専用**（差分テストの第一義オラクル, accuracy.md §3.1）。本番（Standard）
  経路は AnalyticalEphemeris（VSOP87D/ELP）で、DE は組み込まない。
