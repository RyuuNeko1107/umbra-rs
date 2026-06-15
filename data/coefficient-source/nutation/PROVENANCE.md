# 章動係数 原データ provenance（ISSUE-040 / `docs/data-sources.md` §0/§5）

IAU 2000_R06 章動（IAU2006 整合・eclipse フレーム連鎖 `iauPnm06a` → `iauNut06a` が消費する系列）。

## 取得元（一次配布）

IERS Conventions Centre, IERS Conventions (2010), Chapter 5, additional information。

| ファイル | 内容 | URL |
|---|---|---|
| `tab5.3a.txt` | 黄経章動 Δψ（IAU 2000_R06） | https://iers-conventions.obspm.fr/content/chapter5/additional_info/tab5.3a.txt |
| `tab5.3b.txt` | 黄道傾斜章動 Δε（IAU 2000_R06） | https://iers-conventions.obspm.fr/content/chapter5/additional_info/tab5.3b.txt |

- 取得日: 2026-06-15（UTC）。TLS 検証あり（システム CA）。
- 取得は開発時の明示手順（本記録）による一回限り。**ビルド時/実行時のネットワーク取得は行わない**（`docs/accuracy.md` §5）。`xtask generate-coefficients --dataset nutation-iau2000a` は本ディレクトリのローカルファイルのみを入力とする。

## 元データ checksum（SHA-256）

```
tab5.3a.txt  6da73bfe10873ac815520d00fffd67114d647a34afebc5946cfc275e73693f32
tab5.3b.txt  f0dff02c78809b629cc64e2a9fbeffaea5ae20f67e1a62a0ed966f8624807557
```

## ライセンス・帰属

IERS 公開データ（科学データ）。明示的 SPDX ライセンスは付かないが IERS Conventions Centre が
公開配布しており、帰属表示で再配布可（`docs/data-sources.md` §0/§6）。**SOFA C は参照のみ・非移植**。
一次出典: IERS Conventions (2010), G. Petit & B. Luzum (eds.), IERS Technical Note No. 36。
章動理論: Mathews, Herring & Buffett (2002), *JGR* 107(B4); IAU2006 調整 = Capitaine, Wallace & Chapront (2003)。

## 系列の構造（パース仕様の正本）

- 両表とも各データ行 = `i, <sin 係数>, <cos 係数>, l l' F D Ω L_Me L_Ve L_E L_Ma L_J L_Sa L_U L_Ne p_A`（14 整数乗数）。
  - `tab5.3a`（Δψ）: 見出し `A_i A"_i`。Δψ = Σ[A_i·sin(ARG) + A"_i·cos(ARG)] + t·Σ[A'_i·sin + A"'_i·cos]。
  - `tab5.3b`（Δε）: 見出し `B"_i B_i`。Δε = Σ[B_i·cos(ARG) + B"_i·sin(ARG)] + t·Σ[B'_i·cos + B"'_i·sin]。
  - **位置規約（両表共通）**: col2 = sin の係数、col3 = cos の係数（5.3b は名前順が `B"(sin), B(cos)` で col 位置と一致）。
- ブロック（`j = 0` 定数項 / `j = 1` ×t 項）と項数:
  - Δψ: j=0 → **1320**、j=1 → **38**。
  - Δε: j=0 → **1037**、j=1 → **19**。
- 振幅単位 = マイクロ秒角（µas）。打切り 0.1 µas。基本引数式は IERS Conventions 2003（評価は ISSUE-035）。
