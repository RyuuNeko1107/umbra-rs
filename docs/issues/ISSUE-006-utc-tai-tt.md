# ISSUE-006: UTC / TAI / TT 変換（閏秒テーブル・境界）

- crate: umbra-core（時刻系型）／ 閏秒データ取り込みは architecture §11 の data/ + xtask 方針に従う
- 依存: ISSUE-001, ISSUE-004（`JulianDate2`）, ISSUE-005（暦変換）
- モード(tdd-workflow): strict（公開仕様・永続データ形式・閏秒境界の安全性。誤ると全 UTC 接触時刻が秒単位でずれるため strict）

## 目的
UTC ↔ TAI ↔ TT の相互変換と閏秒テーブルを実装する。
- 閏秒テーブル（TAI−UTC、1972– 積算。data-sources §3.2）を versioned + checksum で取り込み、期限切れ検知（valid_to）を持つ。
- UTC の閏秒境界（23:59:60）の正しい処理。
- `TtInstant` / `TaiInstant` / `UtcInstant` 型と相互変換（api-draft §1.3 / §3.2 `TimeScales`）。

## 非目的
- UT1 / ΔT / EOP（ISSUE-007）。本 issue は UTC↔TAI↔TT のみ（UT1 を除く決定論的部分）。
- TDB（ISSUE 範囲外。ephemeris 側で TT→TDB 近似）。ただし型 `TtInstant` は本 issue で確立。
- 1972 年以前の UTC（当時の周波数オフセット/レート閏秒）。v1.0 範囲は 1900–2100 だが、1972 以前は TAI−UTC が階段でないため**ΔT 経由（ISSUE-007）で扱う**方針とし、本 issue の閏秒テーブルは 1972– を対象（要確認: 1900–1972 の UTC 入力時の挙動を §実装メモで規定）。

## 公開インターフェース
api-draft §1.3 / §3.2 を転記・具体化:
```rust
#[derive(Clone, Copy, Debug, PartialEq, PartialOrd)] pub struct UtcInstant(/* JulianDate2 ベース */);
#[derive(Clone, Copy, Debug, PartialEq, PartialOrd)] pub struct TaiInstant(/* JulianDate2 */);
#[derive(Clone, Copy, Debug, PartialEq, PartialOrd)] pub struct TtInstant(/* JulianDate2 */);

pub struct LeapSecondTable { /* versioned + checksum, valid_to */ }
impl LeapSecondTable {
    pub fn bundled() -> Self;                       // data/leap-seconds 同梱
    pub fn tai_minus_utc(&self, utc: UtcInstant) -> Result<f64, TimeError>; // 秒(整数)
    pub fn valid_to(&self) -> UtcInstant;
}

// 変換（TimeScales の一部。api-draft §3.2）
//   utc_to_tai / tai_to_utc（閏秒テーブル）
//   tai_to_tt  / tt_to_tai （定数 32.184 s）
//   utc_to_tt  / tt_to_utc （合成）
impl UtcInstant { pub fn from_gregorian(..) -> Result<Self, TimeError>; pub fn to_gregorian(self) -> (..); }
```
- 生 f64 JD を渡さない（conventions §6 / §11）。閏秒未取得時は `TimeError::MissingLeapSecondData`。

## 数式・アルゴリズムの出典
- **TT − TAI = 32.184 s（定数・定義値）**: IAU 1991 決議、`IERS Conventions 2010 (TN 36) §1`、SOFA `iauTaitt`/`iauTttai` と同値。magic number ではなく定義定数として明記。
- **TAI − UTC（閏秒積算）**: IERS Bulletin C / IANA `leap-seconds.list`（data-sources §3.2）。SOFA `iauDat`（うるう秒テーブル）に対応する手順。
- **UTC↔TAI 境界処理**: SOFA `iauUtctai` / `iauTaiutc`。UTC は閏秒挿入で日が 86401 秒を持つことを考慮（SOFA は "quasi-JD" 規約）。
- 要確認: UTC を `JulianDate2` で持つ際の閏秒日表現（SOFA の `iauUtctai` は閏秒日に非線形性を持つ）。SOFA 手順をコメントに転記し採用。

## 単位 / 時刻系 / 座標系
- 入力/出力: 各 `*Instant`（内部 `JulianDate2`、日）。`tai_minus_utc` は秒（整数）。
- 時刻系: UTC / TAI / TT（conventions §6）。UT1 は本 issue 範囲外。
- 座標系: 無関係。

## アルゴリズム概要
1. 閏秒テーブルを `data/leap-seconds/` から versioned + checksum で読み込み（実行時ネットワーク禁止・accuracy.md §5）。各エントリ: 発効 UTC 日 → 積算 TAI−UTC（秒）。
2. `utc_to_tai`: 該当区間の TAI−UTC を引き、UTC に加算（閏秒日の 86401 秒境界を SOFA 手順で処理）。
3. `tai_to_utc`: 逆引き。閏秒の瞬間（23:59:60）への写像を一意化。
4. `tai_to_tt`: `+32.184 s`（定数）。`tt_to_tai`: `−32.184 s`。
5. `utc_to_tt` / `tt_to_utc`: 上記を合成。
- 数値安定性: 全加減算は `JulianDate2::add_seconds`（ISSUE-004）で µs を保つ。閏秒境界は分岐を明示しテストで固定。禁止: 閏秒の magic 値直書き（テーブル化）、生 f64 時刻。

## 受け入れテスト
accuracy.md テストレベル **L2（時刻）**。
- 既知値（オラクル＝IERS leap-seconds.list と SOFA `iauUtctai`/`iauTaitt`。実装からコピーしない）:
  - 2017-01-01 00:00:00 UTC で TAI−UTC = 37 s（2017 閏秒後）。
  - 任意の現代 UTC で `tt = utc + (TAI−UTC) + 32.184 s`。
  - 往復: `tt_to_utc(utc_to_tt(t)) ≈ t`（多数の日付、許容内）。
- 閏秒境界: 2016-12-31 23:59:60 UTC（直近の挿入閏秒）が表現でき、その前後 1 秒で TAI が連続増加すること。23:59:59 → 23:59:60 → 00:00:00 の TAI 差が各 1 s。
- 境界値: 閏秒発効日の 0 時直前/直後、テーブル先頭（1972）、`valid_to` 超過。
- 異常系: `valid_to` を超える未来 UTC → `MissingLeapSecondData`（または明示的な期限切れエラー）。テーブル未ロード時 → 同。秒=60 が閏秒日以外 → `InvalidDate`。

## 許容誤差
accuracy.md §2.3（絶対 UTC: ΔT は別、UT1−UTC <数 ms）と §2.1（最大食 ±1.5s）から:
- UTC↔TAI↔TT 変換は**決定論的**（閏秒は整数秒、32.184 は定義値）。よって変換誤差は丸めのみ: 現代時刻で **≤ 1µs**。根拠: ΔT/UT1（ISSUE-007）が秒オーダの不確実性を持つのに対し、この決定論部は µs まで詰められ、詰めるべき（誤差を ΔT 側に隠さない・conventions §11）。
- 閏秒整数値・32.184 s は**厳密一致**（オフバイワン無し）を合格条件とする。

## 実装メモ
- 1972 以前 UTC の扱いは「要確認」: 当時 UTC は階段閏秒でない。v1.0 では 1900–1972 の UTC 入力を ΔT 経由（ISSUE-007）で UT1/TT に橋渡しする方針をコメントで固定するか、`UnsupportedTimeRange` とするかをレビューで決める。
- `LeapSecondTable` は `DataSetMetadata`（name/version/source/license/valid_from/valid_to/checksum, architecture §11）を保持。期限切れ検知を必須に。
- レビュー重点: 閏秒境界（23:59:60）の双方向一意性、`valid_to` 超過時のエラー、32.184 と閏秒積算の合成順序。SOFA `iauUtctai` 規約との整合コメント。
