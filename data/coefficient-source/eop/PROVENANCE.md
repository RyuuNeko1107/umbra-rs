# IERS EOP 14 C04 原データ provenance（ISSUE-007 EOP / `docs/data-sources.md` §0/§3.1/§5）

地球姿勢パラメータ（UT1−UTC, 極運動 xp/yp）の IERS EOP **14 C04**（ITRF2014 対応・IAU2000A）系列。
日食の UT1/極運動供給（`umbra_core::IersEopData`）と将来 UTC 精度の不確実性帯の源（accuracy.md §0/§2.3）。

## 取得元（一次配布）

IERS Earth Orientation Centre（Paris Observatory）。

| ファイル | 内容 | URL |
|---|---|---|
| `eopc04_14_IAU2000A_1962-now.txt` | EOP 14 C04 日次（1962–現在）: MJD, x, y, UT1−UTC, LOD, dX, dY ＋各誤差 | https://datacenter.iers.org/data/latestVersion/EOP_14_C04_IAU2000A_one_file_1962-now.txt |

- 取得日: 2026-06-17（UTC）。TLS 検証あり（システム CA）。系列上端は本ファイル時点で 2026-01-03（MJD 61043）。
- 取得は開発時の明示手順（本記録）による一回限り。**ビルド時/実行時のネットワーク取得は行わない**
  （`docs/accuracy.md` §5）。`xtask generate-coefficients --dataset eop-c04` は本ディレクトリのローカル
  ファイルのみを入力とする。
- 系列の上端（valid_to）は配布更新で延びる。更新は xtask 再生成＋checksum 更新で行い、`DataSetMetadata`
  に取得日・系列版・valid_to を固定（data-sources §3.1）。

## 元データ checksum（SHA-256）

```
eopc04_14_IAU2000A_1962-now.txt  9e26da8bc2c8490828f0c1e6ef9587b5d40716e8bac590cf3ab4149f855e5531
```

## ライセンス・帰属

IERS 公開データ（科学データ）。明示的 SPDX ライセンスは付かないが IERS Earth Orientation Centre が
公開配布しており、帰属表示で再配布可（`docs/data-sources.md` §0/§6）。GPL 派生物ではない。
一次出典: IERS Earth Orientation Centre, EOP (IERS) 14 C04 time series, https://hpiers.obspm.fr/eoppc/eop/eopc04/。
記述: C04.guide.pdf（https://hpiers.obspm.fr/eoppc/eop/eopc04/C04.guide.pdf）。

## 系列の構造（パース仕様の正本）

- ヘッダ（テキスト）の後に日次データ行が続く。配布の FORTRAN フォーマット宣言:
  `FORMAT(3(I4),I7,2(F11.6),2(F12.7),2(F11.6),2(F11.6),2(F11.7),2(F12.6))`。
- データ行の列（空白区切りで抽出可能）:
  `year  month  day  MJD  x["]  y["]  UT1-UTC[s]  LOD[s]  dX  dY  x_err  y_err  UT1-UTC_err  LOD_err  dX_err  dY_err`。
- 取り込む列（数値事実のみ）: **MJD（整数, 0h UTC）/ x（arcsec）/ y（arcsec）/ UT1−UTC（秒）**。LOD・dX/dY・各誤差は
  v0.1 では未使用（`EopRecord` は MJD・UT1−UTC・x・y のみ。将来 LOD/誤差は `#[non_exhaustive]` で追加余地）。
- 0h UTC の整数 MJD のみ（日内補間は `IersEopData` の線形補間）。閏秒跨ぎ日対の補間は UT1−TAI 化で扱う
  （`umbra_core::eop` の TODO・後続整備）。
