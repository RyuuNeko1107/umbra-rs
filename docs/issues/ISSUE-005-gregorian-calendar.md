# ISSUE-005: Gregorian calendar conversion（暦⇔JD・負の年規約）

- crate: umbra-core
- 依存: ISSUE-001, ISSUE-004（`JulianDate2`）
- モード(tdd-workflow): strict（公開仕様＝`UtcInstant::from_gregorian`/`to_gregorian` の基礎。負の年・暦の境界規約は誤りやすく永続的な前提となるため strict）

## 目的
グレゴリオ暦（年月日時分秒）⇔ `JulianDate2` の双方向変換を提供する。
- 負の年規約を一意に固定（天文学的年番号 = 0 年あり、ISO 8601 準拠）。
- 暦から JD、JD から暦（小数秒まで）への高精度復元。
- 本プロジェクト v1.0 範囲 1900–2100 を主対象とし、DE441 範囲程度まで破綻しないこと。

## 非目的
- ユリウス暦（1582 改暦前）の扱い。v1.0 は**プロレプティック・グレゴリオ暦**で統一する（要確認: NASA カタログ等が改暦境界でユリウス暦を使う点は、照合時の慣習差として accuracy.md に記録。本変換はグレゴリオ一本）。
- 閏秒（UTC の 60 秒目）。これは ISSUE-006 の UTC↔TAI 側で扱う。本 issue の `to/from_gregorian` は名目上の暦演算。
- タイムゾーン（常に時刻系内の値、ローカル時刻変換は対象外）。

## 公開インターフェース
api-draft §1.3（`UtcInstant::from_gregorian` 等）の下層となる純暦変換:
```rust
// umbra-core 内（時刻系非依存の暦⇔JD）
pub fn gregorian_to_jd(y: i32, mo: u8, d: u8, h: u8, mi: u8, s: f64)
    -> Result<JulianDate2, TimeError>;     // InvalidDate を返しうる
pub fn jd_to_gregorian(jd: JulianDate2) -> (i32, u8, u8, u8, u8, f64);
```
- これを `UtcInstant::from_gregorian` / `to_gregorian`（api-draft §1.3）が利用する。生 f64 JD を公開せず `JulianDate2` を介す（conventions §6）。

## 数式・アルゴリズムの出典
- **Meeus, "Astronomical Algorithms", 2nd ed., Chapter 7 "Julian Day"**:
  - 暦→JD: 式 (7.1)（B = 2 − A + A/4 補正を含むグレゴリオ分岐）。本プロジェクトはプロレプティック・グレゴリオなので**常にグレゴリオ分岐**を採用（Meeus のユリウス/グレゴリオ判定は使わない旨コメント）。
  - JD→暦: Chapter 7 の逆変換手順（Z, F, α, A, B, C, D, E から復元）。
- 妥当性確認用: **IAU SOFA `iauCal2jd` / `iauJd2cal`**（こちらは MJD 基準・グレゴリオ。範囲チェックあり）。SOFA と Meeus の年番号規約差（負年・0 年）を実装コメントで明記。
- 年規約: ISO 8601 / 天文学的年番号（1 BC = 0 年, 2 BC = −1 年）。Meeus も天文学的番号を用いる。

## 単位 / 時刻系 / 座標系
- 入力: 年(i32)・月(1–12)・日(1–31)・時(0–23)・分(0–59)・秒(f64, 0–60 未満。閏秒は本層では非対応＝60 は InvalidDate)。
- 出力: `JulianDate2`（日）/ 暦 6 タプル。
- 時刻系: なし（暦演算は時刻系非依存。UTC の意味は上位で付与）。
- 座標系: 無関係。

## アルゴリズム概要
1. `gregorian_to_jd`: 入力検証（月・日・時分秒の範囲、月ごとの日数・閏年）。Meeus 7.1 のグレゴリオ式で整数部 JD を計算し、時分秒を `(h*3600+mi*60+s)/86400` 日として `JulianDate2` の part2 に格納（part1 は 0.5 境界＝JD は正午起算に注意）。
2. `jd_to_gregorian`: Meeus 逆手順。小数秒は part2 から復元し、丸めで 60 秒に達しないようガード。
3. 閏年判定はグレゴリオ規則（4 で割れ、100 で割れず、400 で割れる）。定数化、magic number コメント。
- 数値安定性: 時分秒は part2 に分離して `JulianDate2`（ISSUE-004）の精度を活かす。JD の 0.5 起算（正午）と暦日の境界処理を厳密に。禁止: 単一 f64 JD への押し込み（conventions §6/§11）。

## 受け入れテスト
accuracy.md テストレベル **L2（時刻）**。
- 既知値（オラクル＝Meeus AA 2nd ed. Ch.7 の worked examples、および SOFA 出力。実装からコピーしない）:
  - 2000-01-01 12:00:00 → JD 2451545.0（J2000.0、Meeus 記載値）。
  - 1987-01-27 00:00 → JD 2446822.5（Meeus 例）。
  - 1957-10-04.81（Sputnik、Meeus 例）。
  - 往復: 1900–2100 の多数日付で `jd_to_gregorian(gregorian_to_jd(...))` が元に戻る。
- 負の年: 0 年（=1 BC）、−1 年（=2 BC）の往復。Meeus/SOFA の天文学的年番号と一致。
- 境界値: 月末・閏日（2000-02-29 有, 1900-02-29 無）、年跨ぎ 12-31→01-01、JD 0.5 境界（暦日の正午/真夜中）。
- 異常系: 2 月 30 日, 13 月, 秒=60（閏秒は本層非対応）→ `TimeError::InvalidDate`。

## 許容誤差
accuracy.md §2.1 の時刻バジェット観点（最終 µs オーダ目標）から:
- 往復誤差: 現代日付で **≤ 1µs（≈1.16e-11 日）**。根拠: `JulianDate2`（ISSUE-004）の精度上限と一致させ、暦変換が時刻精度の律速にならないこと。
- 既知 JD（J2000 等）との一致: 整数日部は**厳密一致**、小数秒部は ≤ 1µs。
- 負年・遠隔年代は丸め余裕を持たせるが、v1.0 範囲 1900–2100 では上記を厳守。

## 実装メモ
- プロレプティック・グレゴリオ統一の判断は「要確認」事項として記録: NASA 5千年カタログ（data-sources §4.1）は歴史日付でユリウス暦を使う場合があり、照合時に系統差が出る。照合は慣習差を明記した第二義チェックに留める（accuracy.md §3.1）。
- 閏秒（UTC の 23:59:60）は本層で `InvalidDate`。UTC としての 60 秒目は ISSUE-006 の `TimeScales` が扱う（暦変換と閏秒適用の責務分離）。
- レビュー重点: 負の年・0 年の往復、JD 0.5（正午起算）境界、グレゴリオ分岐固定（ユリウス暦に落ちない）こと。
