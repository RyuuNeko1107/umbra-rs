# ISSUE-033: VSOP87 係数生成パイプライン（xtask）

- crate: xtask（生成物は umbra-ephemeris に組込み）
- 依存: ISSUE-012/013（消費側の packed 形式要件）, umbra-core（checksum/メタ型は共有 or xtask 内）
- モード(tdd-workflow): **strict**（データ生成・永続形式・provenance/checksum は再現性とライセンスに直結する公開仕様。data-sources §0/§5 で strict 指定）

## 目的
VSOP87D 原データ（IMCCE / CDS VizieR 一次配布）を取り込み、打切り・packed 化・checksum/provenance 付与までを行う決定的な生成パイプラインを `cargo xtask generate-coefficients`（VSOP87部分）として実装する。生成物を crate 組込みバイナリにし、`cargo xtask verify-generated` で CI 差分検査できるようにする。

## 非目的
- 級数評価ロジック（ISSUE-013）。
- ELP/MPP02 側（ISSUE-034）。
- 実行時ネットワーク（禁止。取得は xtask の明示手順 or 手動配置。accuracy.md §5）。
- 惑星係数（v1.0 は地球のみ採用。他は取り込まないか打切りで除外）。

## 公開インターフェース
- `cargo xtask generate-coefficients --target vsop87`（入力: `data/coefficient-source/vsop87/`、出力: `generated/vsop87/`）。
- `cargo xtask verify-generated --target vsop87`（再生成 → checksum 突合、差分で非ゼロ終了）。
- 生成物に `DataSetMetadata { name, version, source, license, valid_from, valid_to, checksum }`（architecture §11）を同梱。packed 形式（f64 振幅/位相/振動数の連続配列 + 変数/べき/項数オフセット表）。消費側（ISSUE-013）と byte-for-byte の契約を固定。

## 数式・アルゴリズム（形式・出典）
- **VSOP87**: Bretagnon & Francou (1988), *A&A* 202, 309。原ファイル形式 = 各変数（L, B, R）×べき α(0..5) ごとの項リスト（各項 S, K, A, B, C … 実際の配布カラムは VSOP87 ヘッダ仕様に従う。A=振幅, B=位相, C=振動数）。data-sources §2.1。
- 採用版 **VSOP87D**（黄道・平均分点 of date・球面、heliocentric）。地球の L/B/R のみ抽出。
- 打切り: 振幅 |A| の寄与順 or DE 差分実測（accuracy.md §3.3）で **残差 0.05″ 級**を切る最小項数。打切り後の項数・推定残差をメタに記録。
- packed: little-endian f64、位相/振動数はラジアン・ラジアン/千年へ正規化（生成時に単位確定）。T 単位 = ユリウス千年（ISSUE-013 と一致）。

## 単位 / 時刻系 / 座標系
- 係数単位: A = 無次元(L,B) / AU(R)、B = rad、C = rad/千年。座標 = 黄道 of date。時刻基準 T = J2000 TDB からのユリウス千年。
- これらを生成物メタ（packed ヘッダ）に明記し、ISSUE-013 が読む。

## アルゴリズム概要
1. 一次原データ（VSOP87D 地球ファイル）を `data/coefficient-source/vsop87/` から読む（取得手順・取得日・元 checksum を記録）。
2. パース → 変数×べき×項に正規化。
3. 打切り適用（しきい値 or 実測項数。ISSUE-013 の DE 差分結果と往復）。
4. packed バイナリへシリアライズ + provenance（source URL/論文引用/取得日）+ checksum(SHA-256) 算出。
5. `generated/vsop87/` へ書き出し、NOTICE（引用 Bretagnon & Francou）を併置。

## 受け入れテスト
- **決定性**: 同一原データから2回生成 → byte-identical（checksum 一致）。`verify-generated` が一致時 0 / 改変時 非0。
- **ラウンドトリップ**: packed を読み戻し、原データの代表項（最大振幅項など）と数値一致（パース欠落なし）。
- **打切りガード（accuracy.md §3.3 を CI/nightly に落とす）**: 生成 packed を ISSUE-013 評価器へ通し、DE440 差分で **残差 0.05″ 級**を満たすことを確認（手動/nightly、巨大データのため。§3.1）。未達なら打切り項数を増やす手順を文書化。
- **版D・地球系列の明示チェック（B4(c)）**: ロードした系列が **VSOP87 版D（黄道・平均分点 of date・球面、heliocentric）かつ地球系列（VSOP87D.ear に相当, body=Earth）**であることを明示的に検証する。EMB（地球–月重心）主系列や他版（A/B/C/E）を取り違えると **6.4″/月オーダーの系統誤差**になる（地球と EMB の差は月軌道による月次振動）。
  - 原ファイルのヘッダ/ファイル名（版識別子・body 識別子）をパース時に検査し、版!=D または body!=Earth なら生成を**非ゼロ終了で失敗**させる。
  - **反転テストでは検出不可**: 太陽地心＝地球日心の符号反転（ISSUE-013）は系列の絶対値の取り違えを打ち消さないため、版/body 誤用を検出できない。よって**版D・地球系列であることはここ（生成時）で明示チェックする**。生成メタ（DataSetMetadata.name/version/source）に版識別子・body を記録し、ISSUE-013 がロード時に再検査できるようにする。
- **メタ完全性**: DataSetMetadata 全フィールド非空、license 欄に「科学データ・明示OSSなし・出典明記」記述。**version 欄に「VSOP87D」、name/source に地球系列であることを明記**。

## 許容誤差
- 生成データ起因の太陽位置残差 = **0.05″ 級**（accuracy.md §2.4 太陽側打切り残差）。これを満たす最小項数を採用（過剰打切り回避 §2.4）。
- 数値根拠: 太陽 0.20″ 総配分（§2.1）の内、打切り 0.05″・光行差等 0.10″・フレーム 0.05″ の RSS 内訳。

## 実装メモ（ライセンス必須記載 — data-sources §0）
- **VSOP87 係数は明示的 OSS ライセンスを持たない科学データ**（IMCCE / 旧 Bureau des Longitudes 由来）。→ **一次配布元の原データから自前生成し、係数を事実として扱う**。Bretagnon & Francou (1988) を引用。
- **GPL コード/GPL 派生係数を取り込まない**（本 crate は MIT/Apache-2.0 想定）。他者再実装は参照のみ。
- 生成物に **出典・取得日・checksum・NOTICE** を必ず付す。生成係数は**ライセンス区分が明確なモジュール**に隔離（architecture §11, data-sources §5）。
- OSS 公開前チェック（data-sources §6）: VSOP87 原データの再配布可否を IMCCE/CDS で一次確認するまで「ライセンス確認済み」としない。`cargo-deny` allow-list に載せる。
- 取得はネットワーク自動化せず、xtask は「ローカル配置済み原データ」を入力とする（実行時/ビルド時ネットワーク禁止、accuracy.md §5）。
