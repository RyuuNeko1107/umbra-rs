# ISSUE-040: IAU2000A 章動係数データ管理（nut00a 1365項・provenance + checksum）

- crate: xtask（生成物は umbra-ephemeris に組込み）
- 依存: ISSUE-035（章動係数の消費側 = `iauNut00a` 相当の Δψ/Δε 評価器・packed 形式要件）, ISSUE-033/034（同格の生成パイプライン形式・`DataSetMetadata`/checksum 規約を踏襲）, umbra-core（checksum/メタ型は共有 or xtask 内）
- milestone: M3（フレーム・章動層。035 と同時に必要。033/034 と同格のデータ生成）
- モード(tdd-workflow): **strict**（データ生成・永続形式・provenance/checksum は再現性に直結する公開仕様。1365項の手書き直書き禁止を担保する基盤。033/034 と同格 strict）

## 目的
IAU2000A 章動の係数（`nut00a`、Δψ/Δε の **1365 項級数**）を、ISSUE-033（VSOP87）/ISSUE-034（ELP/MPP02）と**同格の provenance + checksum 付き生成パイプライン**で管理する。**手書き配列直書きを禁止**（architecture §11, conventions §11）し、一次出典（IERS Conventions 2010 表）から決定的に生成・packed 化して crate へ組込む。
- `cargo xtask generate-coefficients --target nut00a`（IAU2000A 章動係数）。
- `cargo xtask verify-generated --target nut00a`（再生成 → checksum 突合、差分で非ゼロ終了）。
- 生成物に `DataSetMetadata`（name/version/source/license/valid_from/valid_to/checksum）を同梱し、ISSUE-035 の章動評価器（`iauNut00a` 相当）が byte-for-byte 契約で読む。
- 併せて内部粗スキャン用 IAU2000B（`nut00b`, 77 項・非公開・存置）も同パイプラインで管理可能にする（既定は 2000A。公開出力は常に 2000A。要確認: 2000B を同 target に含めるか別 target か）。

## 非目的
- 章動級数の評価ロジック（Δψ, Δε の計算。ISSUE-035）。本 Issue は係数の取り込み・packed 化・checksum まで。
- 歳差・frame bias・ERA・極運動（ISSUE-035）。
- VSOP87/ELP-MPP02 係数（ISSUE-033/034。本 Issue はその形式・規約を踏襲する別 target）。
- 実行時/ビルド時ネットワーク（禁止。取得は xtask の明示手順 or 手動配置。accuracy.md §5）。

## 公開インターフェース
ISSUE-033/034 と同形式（xtask サブコマンド + packed 形式 + メタ）:
- `cargo xtask generate-coefficients --target nut00a`（入力: `data/coefficient-source/nutation/`、出力: `generated/nutation/`）。
- `cargo xtask verify-generated --target nut00a`（再生成 → checksum 突合、差分で非ゼロ終了）。
- 生成物に `DataSetMetadata { name, version, source, license, valid_from, valid_to, checksum }`（architecture §11）。packed 形式（各項の lunisolar/planetary 引数の整数係数 + 振幅係数 A,A',B,B' 等を連続配列 + 項数/区分オフセット表）。消費側（ISSUE-035）と byte-for-byte の契約を固定。

## 数式・アルゴリズム（形式・出典）
- **IAU2000A 章動係数**: Mathews, Herring, Buffett (2002) / IAU2000A。配布表 = **IERS Conventions 2010, ch.5 の章動級数表**（lunisolar 678 項 + planetary 687 項 = 1365 項）。SOFA `iauNut00a` が同係数を内蔵（**SOFA C は参照のみ・移植/コピーしない**、conventions §11 / data-sources §0）。
- 各項の形式: 5 つの基本引数（Delaunay 引数 l, l', F, D, Ω）＋惑星引数の**整数倍係数**と、Δψ・Δε の **振幅係数（in-phase / out-of-phase, sin/cos）**。data-sources §3 相当（章動表）。
- 単位: 振幅は配布表の単位（0.1 µas または µas。取り込み境界でラジアン化方針を生成時に確定し packed ヘッダに明記。要確認: 配布表のスケール単位）。引数係数は整数。
- 打切り: IAU2000A は**全項採用**（accuracy.md §2.1 フレーム 0.05″ 配分に対し実力 ~1mas で余裕、打切りは内部粗スキャン（非公開）の 2000B 側で対応）。2000A は項を間引かない（精度最優先方針）。

## 単位 / 時刻系 / 座標系
- 係数単位: 振幅 = µas（or 0.1µas）→ packed 化時にラジアン基準へ正規化（生成時単位確定、packed ヘッダに明記）。引数係数 = 整数。時間引数（基本引数の多項式）= TT のユリウス世紀 from J2000（ISSUE-035 と一致）。
- 座標系: 章動は CIP/真黄道に関わる量。本 Issue は係数データのみで座標変換はしない（ISSUE-035）。
- これらを生成物メタ（packed ヘッダ）に明記し、ISSUE-035 が読む。

## アルゴリズム概要
1. 一次原データ（IERS Conventions 2010 章動表の機械可読版、または SOFA 配布の係数表を**出典としてのみ参照しデータは IERS 一次から**）を `data/coefficient-source/nutation/` から読む（取得手順・取得日・元 checksum 記録）。
2. パース → lunisolar/planetary 区分 × 引数整数係数 × 振幅係数（in/out-phase）に正規化。1365 項の欠落なし検査。
3. packed バイナリへシリアライズ + provenance（IERS Conventions 2010 ch.5 引用・取得日）+ checksum(SHA-256)。
4. `generated/nutation/` へ書き出し、NOTICE（引用）を併置。
5. （任意）`nut00b` 77 項も同手順（内部粗スキャン用・非公開・存置、ISSUE-035 `iauNut00b` 相当）。
- 数値安定性: 単位スケールの取り違え（µas vs 0.1µas）を生成時テストでガード。整数係数のオーバーフロー/型を固定。

## 受け入れテスト
accuracy.md テストレベル **L3 の基盤**（係数データ完全性）。基準値は IERS/SOFA 公開ベクトルから（実装コピー禁止、conventions §11）。
- **決定性**: 同一原データから2回生成 → byte-identical（checksum 一致）。`verify-generated` が一致時 0 / 改変時 非0（ISSUE-033/034 と同等）。
- **項数完全性**: lunisolar 678 + planetary 687 = **1365 項**を読み出せる（欠落・重複なし）。
- **ラウンドトリップ**: packed を ISSUE-035 の評価器へ通し、既知 TT の Δψ/Δε が SOFA `iauNut00a` 参照値と一致（**~1mas 級 → 実際は µas 級一致**、accuracy.md §2.1）。
- **単位整合**: 代表項（最大振幅 = 主章動項 18.6 年）の振幅が公開値と単位込みで一致（µas スケールの取り違えを検出）。
- **メタ完全性**: `DataSetMetadata` 全フィールド非空、license 欄に「IERS 公開・出典明記・SOFA は参照のみ非移植」記述。
- **2000B（任意）**: 77 項版を同手順で生成・検証できる。

## 許容誤差
- 生成データ起因のフレーム残差は **0.05″ 配分内**（accuracy.md §2.1）。IAU2000A 全項採用のため理論残差は ~1mas 級で大幅余裕。本 Issue は係数の**完全性・単位正確性**を担保（評価誤差は ISSUE-035）。
- SOFA 参照値との Δψ/Δε 一致は **µas 級**を目標（ISSUE-035 のラウンドトリップで実測）。
- 数値根拠: §2.1 RSS（月0.40/太陽0.20/光行差0.10/フレーム0.05/solver0.05）。

## 実装メモ（ライセンス・データ管理 — architecture §11 / conventions §11）
- **手書き配列直書き禁止**（architecture §11）: 1365 項を Rust 配列にハードコードしない。ISSUE-033/034 と**同格のデータ生成パイプライン**で管理（本レビュー Important 040）。
- **SOFA は参照のみ・移植しない**（SOFA は C・独自ライセンス）。係数は **IERS Conventions 2010 表を一次出典**として扱い、provenance + checksum を付す（data-sources §0/§5、conventions §11 magic number 禁止）。
- 生成係数は**ライセンス区分が明確なモジュール**に隔離し NOTICE・引用を同梱（architecture §11）。
- 取得はネットワーク自動化せず xtask は「ローカル配置済み原データ」を入力（実行時/ビルド時ネットワーク禁止、accuracy.md §5）。
- `cargo-deny` allow-list / OSS 公開前ライセンス確認（data-sources §6）の対象に含める。
- レビュー重点: 1365 項完全性、µas 単位スケール、決定性 checksum、消費側（035）byte-for-byte 契約、SOFA 非移植の明記。
