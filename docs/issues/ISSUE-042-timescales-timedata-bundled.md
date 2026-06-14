# ISSUE-042: TimeScales / TimeData 同梱 API（bundled + from_path・閏秒/EOP/ΔT 束ね）

- crate: **umbra-core（`TimeData`/`TimeScales`/`DataSetMetadata` 型定義＝純粋・no_std 互換・データを持たない, B3）。`TimeData::bundled()`（同梱バイト）は umbra-ephemeris に `bundled-data` feature ゲートで置き re-export（core にデータを置かない, B3）**。同梱データ実体は data/ + xtask（architecture §11, data-sources §3）
- 依存: ISSUE-006（UTC/TAI/TT・閏秒テーブル `LeapSecondTable`）, ISSUE-007（`DeltaTModel` / `EarthOrientation` / `IersEopData` / `EspenakMeeusDeltaT` / `utc_to_ut1`）, ISSUE-004（`JulianDate2`）, ISSUE-044（`MissingLeapSecondData`/`MissingEarthOrientationData` を含むエラー集約）
- milestone: M2（時刻層の束ね。006/007 完了後。エンジン組立 043 と standard_engine の前提）
- モード(tdd-workflow): **strict**（**確定A2**: `TimeData::bundled()` の同梱データ API 形と valid_to 超過時の挙動は公開仕様・永続データ契約。誤ると将来 UTC 精度の取り扱いを誤るため strict）

## 目的
**確定A2**（milestone0-review.md）＋**確定B3**に従い、閏秒・EOP・ΔT を束ねる `TimeData`（データ束ね型）と、それを使う `TimeScales`（変換 facade）の同梱 API を提供する。**2 型構成（B3）: `TimeData`=データ束ね、`TimeScales`=変換 facade（TimeData から構築）。変換は Result。両型は umbra-core の純粋型（no_std 互換）**。
- `TimeData::bundled()`（crate **埋込み**データ）＋ `TimeData::from_path()`（外部ファイル）の 2 系統（A2）。**`bundled()` は umbra-ephemeris に `bundled-data` feature ゲートで置き re-export（core にデータを置かない, B3）。off 時は `bundled()` を `#[cfg(feature = "bundled-data")]` で消す**。
- `TimeScales`（変換 facade, B3）: `TimeData` から構築し、`utc_to_tt` / `utc_to_ut1` / `tt_to_utc` を提供（いずれも **Result**, api-draft §3.2、ISSUE-006/007 を束ねる）。
- **valid_to 超過時**は `MissingLeapSecondData` / `MissingEarthOrientationData` を返し、**`CalculationMetadata` に記録**（A2: valid_to 超過は Missing*Data を返し metadata 記録）。
- ΔT は EOP 範囲内は EOP 由来高精度、範囲外は Espenak–Meeus 外挿（ISSUE-007）。不確実性帯を伝播。

## 非目的
- 閏秒/EOP/ΔT データ自体の取り込み実装（ISSUE-006/007。本 Issue はそれらを**束ねる API** と同梱データ供給形）。
- ERA・見かけ恒星時の構成（ISSUE-039。本 Issue は UT1 値までの時刻系変換）。
- データ生成パイプライン（xtask。本 Issue は生成済み同梱データを埋込む側、033/034/040 と同様の versioned + checksum を前提）。
- エラー型の定義そのもの（ISSUE-044。本 Issue は valid_to 超過時に適切な Missing*Data を返す利用側）。

## 公開インターフェース
api-draft §3.2 / §5 を転記・具体化（**確定A2**）:
```rust
/// 閏秒 + EOP + ΔT を束ねた時刻系データ束（A2）。型定義は umbra-core（純粋・データを持たない, B3）。
pub struct TimeData { /* LeapSecondTable + IersEopData + DeltaTModel + DataSetMetadata 群 */ }
impl TimeData {
    /// 埋込み同梱データ。**umbra-ephemeris で `bundled-data` feature ゲート＋re-export**（core にデータを置かない, B3）。
    /// off 時は本コンストラクタが API から消える（コンパイル時不在）。
    #[cfg(feature = "bundled-data")]
    pub fn bundled() -> Self;                                   // crate 埋込みデータ（A2, B3 ephemeris ゲート）
    pub fn from_path(dir: &std::path::Path) -> Result<Self, TimeError>; // 外部ファイル（A2, off 時はこちら必須）
    pub fn coverage(&self) -> TimeRange<UtcInstant>;            // 各データの有効範囲の積
    pub fn metadata(&self) -> &[DataSetMetadata];               // version/source/valid_to/checksum（DataSetMetadata は core 純粋型, B3）
}

/// 時刻系変換 facade（006/007 を束ねる, B3）。TimeData から構築。変換は Result。型定義は umbra-core（純粋）。
pub struct TimeScales { /* TimeData を保持 */ }
impl TimeScales {
    pub fn new(data: TimeData) -> Self;                         // api-draft §3.2 を A2 形へ
    pub fn utc_to_tt(&self, t: UtcInstant) -> Result<TtInstant, TimeError>;   // 閏秒必要
    pub fn utc_to_ut1(&self, t: UtcInstant) -> Result<Ut1Instant, TimeError>; // EOP 必要
    pub fn tt_to_utc(&self, t: TtInstant) -> Result<UtcInstant, TimeError>;
    pub fn delta_t_uncertainty_seconds(&self, t: UtcInstant) -> f64;          // metadata へ（accuracy.md §0）
}
```
- valid_to 超過: 該当変換は `TimeError::MissingLeapSecondData` / `MissingEarthOrientationData`（A2、api-draft §1.6）。範囲外でも `TtInstant`（閏秒のみで足りる近傍）は返せる箇所と UT1 必須箇所を区別。
- 要確認: api-draft §3.2 旧 `TimeScales::new(leap, eop, dt)` 3 引数形を `TimeData` 1 引数形へ集約する破壊的変更（A2 確定に伴う改訂。api-draft §6 未確定 TimeData 項を確定）。

## 数式・アルゴリズムの出典
- 時刻系恒等式（ISSUE-006/007・conventions §6）:
  - `TT = UTC + (TAI−UTC) + 32.184 s`（TT−TAI、ISSUE-041 定数）。
  - `UT1 = UTC + (UT1−UTC)`（EOP、ISSUE-007）。
  - `ΔT = TT − UT1 = (TAI−UTC) + 32.184 − (UT1−UTC)`（EOP 範囲内高精度、範囲外は Espenak–Meeus 外挿）。
- 出典: IERS Conventions 2010 ch.5（時刻系定義）、Espenak–Meeus ΔT 多項式（NASA、ISSUE-007/data-sources §3.3）。閏秒 = IERS Bulletin C / IANA `leap-seconds.list`（data-sources §3.2）。EOP C04（data-sources §3.1）。
- 1972 以前の UTC 非階段年代の橋渡しは ISSUE-007 の方針に従う（ΔT 経由）。要確認: TimeData 束ね時の年代境界での coverage 表現。

## 単位 / 時刻系 / 座標系
- 単位: 時刻は型（`UtcInstant`/`TaiInstant`/`TtInstant`/`Ut1Instant`）。内部 JD は `JulianDate2`（生 f64 禁止、conventions §6）。ΔT/UT1−UTC は秒。
- 時刻系: UTC 入力、TT/UT1 出力。閏秒は UTC↔TAI↔TT、EOP は UT1。conventions §6。
- 座標系: 該当なし（時刻層）。極運動値の供給は ISSUE-007（本 Issue は時刻系のみ束ねる。極運動は EOP からフレーム連鎖 035 が直接取得 ※要確認: 極運動も TimeData 経由で出すか）。

## アルゴリズム概要
1. `TimeData::bundled()`: 埋込みの閏秒/EOP C04/ΔT（versioned + checksum、033/034/040 同様）をロード（実行時ネットワーク禁止、accuracy.md §5）。
2. `TimeData::from_path()`: 外部ディレクトリから同形式を読む（更新運用向け）。checksum 検証、不正は `TimeError`。
3. `TimeScales::utc_to_tt`: 閏秒テーブルで TAI−UTC を取り `+32.184`。閏秒 valid_to 超過は `MissingLeapSecondData` 並びに metadata 記録。
4. `utc_to_ut1`: EOP の UT1−UTC で UT1。EOP coverage 外は `MissingEarthOrientationData` 並びに metadata 記録。
5. `delta_t_uncertainty_seconds`: ISSUE-007 の不確実性帯（近傍 <0.1s、将来数秒）を `CalculationMetadata.delta_t_uncertainty_seconds` へ（accuracy.md §0/§2.3）。
6. `coverage()` / `metadata()`: 各データの valid_from/valid_to/version/checksum を公開（A2: 超過記録の源）。
- 数値安定性: 補間端点外挿を弾く（ISSUE-007）。禁止: 生 f64 時刻、valid_to 無視の沈黙外挿、不確実性 0 の誤報。

## 受け入れテスト
accuracy.md テストレベル **L2（時刻）**。基準値は IERS/NASA 公開値（実装コピー禁止、conventions §11）。
- **bundled / from_path 同等性**: 同一データで `bundled()` と `from_path()` が同一変換結果（A2）。
- **変換正当性**: 近傍年（例 2020）で `utc_to_tt` / `utc_to_ut1` が IERS 由来値と整合（ISSUE-006/007 の許容内）。`tt_to_utc` ラウンドトリップ恒等。
- **valid_to 超過挙動（A2 必須）**: 閏秒 valid_to 超過 UTC で `utc_to_tt` が `MissingLeapSecondData`、EOP coverage 外 UTC で `utc_to_ut1` が `MissingEarthOrientationData`。**いずれも metadata に記録**されること。
- **不確実性帯伝播**: 将来年で `delta_t_uncertainty_seconds` が数秒オーダ、近傍で <0.1s（accuracy.md §0、ISSUE-007 と整合）。
- **メタ完全性**: `metadata()` の各 `DataSetMetadata` が version/source/valid_to/checksum を持つ。
- **データ未ロード/破損**: `from_path` の checksum 不一致・欠落で `TimeError`。

## 許容誤差
accuracy.md §0/§2.3 から（ISSUE-006/007 を束ねるため許容は両 Issue 準拠）:
- TT 変換: 閏秒は整数秒テーブル＋32.184（厳密）。UT1: UT1−UTC 補間 <数 ms（accuracy.md §2.3）。
- ΔT 不確実性帯: 近傍 <0.1s、将来（〜2100）数秒（**誤差でなく不確実性として出力必須**、accuracy.md §0）。
- valid_to 超過は**誤差ではなくエラー＋metadata 記録**で扱う（A2）。沈黙外挿で精度を偽らない。

## 実装メモ
- **確定A2 厳守**: `bundled()`（埋込）＋ `from_path()`、valid_to 超過は Missing*Data を返し metadata 記録。api-draft §6 の TimeData 未確定項を本 Issue で確定（旧 3 引数 `TimeScales::new` を `TimeData` 集約へ改訂、破壊的変更を明記）。
- **crate 配置（確定 B3）**: `TimeData`/`TimeScales`/`DataSetMetadata` **型定義は umbra-core（純粋・no_std 互換・データを持たない）**。`bundled()`（埋込みバイトを返すコンストラクタ）は **umbra-ephemeris に `bundled-data` feature ゲートで置き、上位 crate へ re-export**（core にデータを置かない原則）。**off 時は `#[cfg(feature = "bundled-data")]` で `bundled()` がシンボルごと消える**（from_path 必須）。api-draft §5 の例 `standard_engine(TimeData::bundled())` は `bundled-data` on（既定）前提。ISSUE-043 `standard_engine` が消費。
- 同梱データは versioned + checksum（033/034/040 と同格、architecture §11）。更新は xtask、実行時ネットワーク禁止（accuracy.md §5）。
- エラーは ISSUE-044 集約（thiserror, A6）と整合（`From<TimeError>`）。
- レビュー重点: A2 の bundled/from_path 形、valid_to 超過の Missing*Data + metadata、不確実性帯伝播、**B3 の crate 帰属（型は core 純粋・bundled() は ephemeris の bundled-data ゲート＋re-export・off 時 cfg で消える・DataSetMetadata は core 純粋型）・2 型構成（TimeData/TimeScales）・変換 Result**、旧 API 改訂の破壊的変更明記。
