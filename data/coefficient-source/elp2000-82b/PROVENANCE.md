# 月 ELP2000-82B 原データ provenance（ISSUE-034/014 / `docs/data-sources.md`）

月の地心黄道座標（黄経 V / 黄緯 U / 距離 r）の解析理論 **ELP2000-82B**（Chapront-Touzé & Chapront 1983/1988）。

> **モデル選定（ユーザー判断 2026-06-15）**: issue 当初は ELP/MPP02（~mas、LLR fit）を想定したが、その
> 原データの清浄な一次配布が見つからず（GPL 再実装 ytliu0 等は取り込み禁止）、IMCCE FTP に即入手可能な
> **ELP2000-82B を採用**。1900–2100・月バジェット 0.40″（`docs/accuracy.md` §2.1）には十分な見込み
> （要・評価器での実測検証）。R06 章動と同様の「実用モデル判断」。将来 ELP/MPP02 入手時に差し替え可。

## 取得元（一次配布）

IMCCE（Observatoire de Paris）, `ftp://ftp.imcce.fr/pub/ephem/moon/elp82b/`。Chapront-Touzé M. & Chapront J.,
*A&A* 124, 50 (1983) / *A&A* 190, 342 (1988)。取得日 2026-06-15（UTC）。全 36 ファイル `ELP1`–`ELP36`
（計 2,500,852 B）。**ビルド時/実行時ネットワーク取得は行わない**（`docs/accuracy.md` §5）。

## ファイル構成（系列マップ）

| ファイル | 内容 | 形式 |
|---|---|---|
| ELP1/2/3 | 主問題 経度(sin)/緯度(sin)/距離(cos) | 主問題形式 |
| ELP4-6 / 7-9 | 地球形状摂動 L/B/R（と ×T） | 摂動形式 |
| ELP10-15 | 惑星摂動 表1 L/B/R（と ×T） | 摂動形式 |
| ELP16-21 | 惑星摂動 表2 L/B/R（と ×T） | 摂動形式 |
| ELP22-27 | 潮汐 L/B/R（と ×T） | 摂動形式 |
| ELP28-30 | 月形状摂動 L/B/R | 摂動形式 |
| ELP31-33 | 相対論摂動 L/B/R | 摂動形式 |
| ELP34-36 | 惑星摂動(太陽離心率) L/B/R | 摂動形式 |

`...LONGITUDE/T` 等の `/T` 系列は係数を T（世紀）倍する項。主問題 2,645 項 + 摂動 35,227 項 ≈ **37,872 項**。

## 行フォーマット（パース仕様の正本）

- **主問題（ELP1-3）**: `i1 i2 i3 i4   A   B1 B2 B3 B4 B5 B6`（`4i3,2x,f13.5,6(2x,f10.2)`）。
  - `i1..i4` = Delaunay 引数 `D, l', l, F` の整数乗数。`coef(1)=A`=振幅、`coef(2..7)=B1..B6`=定数偏微分。
  - **B 列は評価器(ISSUE-014)が DE200/LE200 フィット補正に使用**（`elp82b_1`: `tgv=B1+dtasm·B5`,
    `A' = A + tgv·(delnp−am·delnu) + B2·delg + B3·dele + B4·delep`、距離は `A−2A·delnu/3` も）。
    coef(7)=B6 は未使用だが**全 7 実数を忠実保存**（パーサは評価知識を持たない）。
  - 項 = `A'·sin(arg)`（ELP1 経度/ELP2 緯度）/ `A'·cos(arg)`（ELP3 距離）、`arg = i1·D + i2·l' + i3·l + i4·F`。
  - 単位は経度/緯度=秒角、距離=km（評価時に確定。ISSUE-014）。
- **摂動（ELP4-36）**: `m1 .. mK   φ   A   period`。
  - `m1..mK` = 当該系列の引数（Delaunay + 惑星平均黄経 + 歳差等）の整数乗数（K は系列ごとに異なる:
    地球形状=5、惑星表1/2=11、潮汐=5 等。引数の正確な定義は ELP2000-82B 文書、ISSUE-014 で確定）。
  - `φ` = 位相（**度**）、`A` = 振幅、`period` = 周期（日, 参考）。
  - 項 = `A·sin(Σ_k m_k·arg_k + φ·π/180)`。`/T` 系列は全体に T を掛ける。

> 注: 本 Issue（034）は **36 ファイルのパースと packed 化**まで。引数の定義・単位・sin/cos・T 倍の適用
> （= 評価式）は ISSUE-014。パーサは各ファイルの「乗数群 + 係数（主=A / 摂動=φ,A）」を忠実に取り込む。

## ライセンス・帰属

ELP2000-82B 係数は IMCCE 由来の科学データ（明示的 SPDX ライセンスなし）。一次配布元から自前生成し
帰属表示で再配布（`docs/data-sources.md` §0/§6）。Chapront-Touzé & Chapront を引用。**GPL 再実装は不使用**。

## 元データ checksum（SHA-256）

```
ELP1  ae30cbffb83a7bd4582a83a32a322d08a48ba057a4df7bf9dd5df9f06b1688fa
ELP2  c91e5585b0a9e7bd091304b164ce89a6461acd0e439d47957c890aec1e031e08
ELP3  862a8e4c8e70ce8b28383be4c9f2e2c025a8f633d7b7a811eb3afdab4ed9f354
ELP4  f27ea439bf8f4fd35bed31c0a42de5414db07f9fcc3587237f6908891d43d773
ELP5  6803422481e4decae4a59f89d4f94c7b33af21d293d9bc807d565f51c29b9915
ELP6  2a6be4d33dfce4cf2351d295b4aade8f34cc476a97747b3eccb12763658e1fd1
ELP7  35491a0c73ff6bcb136d8f54db89d8df2fb741af43aed9707618e3a925df474d
ELP8  f3e7f4c851e7f9ac1a0556fe7e685e44612f3dca317ab605eb35f565251e020e
ELP9  574347346363df52c7602f56747b790e9cbe60152127d24451c6a1633bc79f0e
ELP10 dbd82ddc6064e4cc7b4f08fa27b2fcb48f82456ad36a850a0d3ddae098c3e2e6
ELP11 0ad7a914c9f98008a648881783c9dd4a14692ec14e2e1bfce4709c68d17bd659
ELP12 8ed7be0ab70f4ffae6b1f711cc4e915257ae5e269fbbfc5a7060f7e952728ba8
ELP13 643295b3894023b4b1bd6ee2b0ecf5d3ff23d703baccc6302caa05fa8b84f76c
ELP14 b59d8b9bbef282f2bead538d6906781257a7fb5b8699bb6adcacb070a76f1e89
ELP15 17ab0d521c178187a5de4847b6696fbcb7a55d77d776568a0c16f33fd3be342a
ELP16 2bef867d8aad4075bc2711559cf1bc42757501bb10c307ff121152bddd344a66
ELP17 6cf0746d034ac75ed60d4d16ed0de790fe5b7aad9df7b462091c38020e1b1bfc
ELP18 b1d93931f6016023c83354cd54a2614978de3b6bc5b7234537d461200f7f4753
ELP19 dd0b0bd5f5c354683f035ee8a09d82e9c138baaf27758ca311d07c76983bdd2e
ELP20 0f1d571879dc9b1a6f7698b403bef26ab151c3ec2ed42cbdec5b297cb0464a8e
ELP21 1546d0e8af01f759f9dfb7a3bf4334a3e579f59e5c22edab29d660dede2ed4b4
ELP22 44263bb254c2b6c0df963bddfd1ecfd01950668fd913f9e6ba6da6ff729f3e41
ELP23 38917cf2cbe0f9afcd271444a47066b098f9a8792836ec85dac359f4c11b9464
ELP24 ca67e5db5933887130709767c1eb3ac009e4bce596978bad7995c416ee71c59f
ELP25 fd0cb03d496cbf23bf7bebef40aa09019c6a50072fc7fd2573645f26a56ab635
ELP26 d5b2a33974b099448a35536987b2f08aff5b11d5801ffa55c87ae8a1d26e9f4f
ELP27 648379d85e1753cc37bc3852b63899d414a05d3aea1c8621dc44e9fdfff234f1
ELP28 0785b8e002887799bc303e8be1abd71c37ae9de671bffbdeb8f50998f892f18e
ELP29 816fd1a94b1cb4e2e5e6cec72971f27ad0f2bf987735d184586fdadd033c3d02
ELP30 cff5ab4c84a6a36855e5b2e2f47e1e0e1d605e789ff2954755dc64d188067a13
ELP31 c2fc53c2442c1b61404991f31c859cc5eb300eb66bb764e1f16c70fe8d199dc3
ELP32 7a07397b63d1ade0909c12be9024632de6ff27fd9e8e410e5f97f821d2390a60
ELP33 459ea9eff9a9d7b5c224245f5edb113060174b4991d3adf7d27079719bcf2339
ELP34 b83178e98bd33e8f26ef6662e03455761ffac7dae399ad1ccb4bf028c5f0e774
ELP35 692d1752a7ea28c7157dbf694750af6ea2cbf74c744c4b792a7e8bf1bb5ad7d7
ELP36 1f8eec292def5ceb4fff9a09ca678bd81e9cdc23ff7bdc802f2281c837357bda
```
