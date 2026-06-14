# ISSUE-034: ELP/MPP02 係数生成パイプライン（xtask、GPL 非混入）

- crate: xtask（生成物は umbra-ephemeris に組込み）
- 依存: ISSUE-012/014（消費側の packed 形式要件）, ISSUE-033（xtask 基盤・DataSetMetadata 共有）
- モード(tdd-workflow): **strict**（データ生成・永続形式・provenance/checksum・**ライセンス隔離**は再現性と法的整合に直結。data-sources §0/§2.2 で strict）

## 目的
ELP/MPP02（DE fit 版）の原データ（IMCCE 一次配布）を取り込み、打切り・packed 化・checksum/provenance 付与までを行う決定的生成パイプラインを `cargo xtask generate-coefficients`（ELP 部分）として実装する。**GPL 再実装/GPL 派生係数を一切混入させない**ことを工程と検査で保証する。

## 非目的
- 級数評価ロジック（ISSUE-014）。
- VSOP87 側（ISSUE-033）。
- 実行時/ビルド時ネットワーク（禁止。accuracy.md §5）。
- LLR 版系列の採用（DE fit 版を採用。data-sources §2.2）。

## 公開インターフェース
- `cargo xtask generate-coefficients --target elp-mpp02`（入力 `data/coefficient-source/elp-mpp02/`、出力 `generated/elp-mpp02/`）。
- `cargo xtask verify-generated --target elp-mpp02`（再生成 → checksum 突合、差分で非0終了）。
- 生成物に `DataSetMetadata`（name/version/source/license/valid_from/valid_to/checksum）。packed 形式 = 主問題級数 + 摂動級数（整数 Delaunay 係数テーブル + 実振幅 + 基本引数多項式係数 + 時間単位/黄道版メタ）。ISSUE-014 と byte 契約を固定。

## 数式・アルゴリズム（形式・出典）
- **ELP/MPP02**: Chapront & Francou (2002), *A&A* 387, 700。構造 = 主問題（Delaunay 引数 D,l,l′,F の整数結合に対する正弦/余弦項）＋ 摂動級数（惑星摂動・地球扁平・潮汐・相対論・固有摂動）。経度 V・緯度 U・距離 r。data-sources §2.2。
- 採用版 **MPP02 DE fit 版**（Reference=JPL DE と整合）。fit 系列識別子をメタに固定。
- 基本引数（D,l,l′,F,平均経度）の時間多項式は MPP02 採用版を原データから抽出。時間単位（世紀 vs 千年）・黄道版（J2000 vs of date）を確定しメタに記録（ISSUE-014 と一致必須。不一致は系統誤差）。
- 打切り: 振幅寄与順 or DE 差分実測（accuracy.md §3.3）で **残差 0.1″ 級**を切る最小項数。

## 単位 / 時刻系 / 座標系
- 係数単位: 振幅 = 角度（V,U は arcsec or rad、r は km。生成時に正規化し確定）、引数係数 = 整数、基本引数多項式 = rad とその時間べき係数。
- 座標 = 黄道（採用版に従う）、時刻基準 = J2000 TDB。単位/版を packed ヘッダに明記。

## アルゴリズム概要
1. 一次原データ（IMCCE ELP/MPP02 DE fit ファイル）をローカル `data/coefficient-source/elp-mpp02/` から読む（取得手順・取得日・元 checksum 記録）。
2. パース → 主問題/摂動/基本引数に分離・正規化。
3. 打切り（しきい値 or 実測項数。ISSUE-014 の DE 差分と往復）。
4. packed シリアライズ + provenance（Chapront & Francou 引用・取得日・URL）+ SHA-256 checksum。
5. `generated/elp-mpp02/` 出力 + NOTICE 併置。**ライセンス区分が明確な隔離モジュール**へ配置。

## 受け入れテスト
- **決定性**: 2回生成 byte-identical。`verify-generated` 一致0/改変非0。
- **ラウンドトリップ**: packed 読み戻しで原データ代表項（最大振幅項）と数値一致（欠落なし）。
- **打切りガード（accuracy.md §3.3 を nightly に）**: 生成 packed を ISSUE-014 評価器へ通し DE440 差分で **残差 0.1″ 級**確認。未達なら項数増やす手順を文書化。距離残差も別途ガード（食分に効く）。
- **GPL 非混入検査（必須・本 Issue の核）**:
  - 原データの provenance が IMCCE 一次配布であること（ytliu0/ElpMpp02 等 GPL 由来でないこと）をメタとレビューで確認。
  - 生成物 checksum が GPL 実装同梱係数と一致しないことの確認手順（参照係数との偶然一致は理論可だが、由来は一次データであることを記録）。
  - `cargo-deny` allow-list に GPL を許可しない設定で CI が通ること。
- **メタ完全性**: 全フィールド非空、license 欄に「科学データ・明示OSSなし・GPL 再実装非混入・出典明記」。

## 許容誤差
- 生成データ起因の月位置残差 = **0.1″ 級**（accuracy.md §2.4 月側打切り残差）。最小項数で達成（過剰打切り回避）。
- 数値根拠: 月 0.40″ 総配分（§2.1）の支配項。打切り 0.1″・補正/フレーム残りは ISSUE-015/035 と RSS。月は相対角速度の主因（1″≈2s）。

## 実装メモ（ライセンス必須記載 — data-sources §0/§2.2）
- **ELP/MPP02 一次データは明示 OSS ライセンスなしの科学データ**（IMCCE 由来）。→ **IMCCE 原データから自前生成**、係数を事実として扱い Chapront & Francou (2002) を引用。
- **Yuk Tung Liu（ytliu0/ElpMpp02）の C++ 実装・生成係数は GPL-3.0。本 crate（MIT/Apache-2.0 想定）に取り込み禁止**。参照（理解・検証）のみ可。コード/数値テーブルを移植しない。GPL 実装の出力を期待値としてテストに貼らない（基準は DE）。
- 生成物に **出典・取得日・checksum・NOTICE** を付し、**ライセンス区分が明確な隔離クレート/モジュール**へ（architecture §11, data-sources §5）。
- OSS 公開前チェック（data-sources §6）: ELP/MPP02 原データ再配布可否を IMCCE で一次確認 + GPL 非混入保証まで「確認済み」としない。
- ネットワーク自動化しない（ローカル配置原データを入力）。
