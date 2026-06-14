# ISSUE-007: UT1 / ΔT abstraction（DeltaTModel / EarthOrientation・IERS EOP・不確実性帯）

- crate: umbra-core（trait・型）／ EOP/ΔT データ取り込みは data/ + xtask（architecture §11, data-sources §3）
- 依存: ISSUE-001, ISSUE-004（`JulianDate2`）, ISSUE-006（UTC/TAI/TT）, ISSUE-002（極運動角 `Radians`）
- モード(tdd-workflow): strict（公開仕様・永続データ形式・将来 UTC 精度を律速する不確実性帯の出力。誤ると約束精度を誤るため最重要 strict）

## 目的
地球回転に関わる時刻・姿勢量を抽象化し、EOP 実データと将来予測の不確実性帯を扱う。
- `DeltaTModel` trait: ΔT = TT − UT1 と**その不確実性帯**（accuracy.md §0/§2.3）。
- `EarthOrientation` trait: UT1−UTC、極運動 (xp, yp)。
- `IersEopData`: IERS EOP **C04** 実データ（履歴 1962–現在＋短期予測）の versioned + checksum 取り込み。
- `EspenakMeeusDeltaT`: 長期 ΔT 外挿（1972 以前・将来）。
- UT1 への変換（`utc_to_ut1` → `Ut1Instant`）。
- 将来・遠隔年代の**不確実性帯を `CalculationMetadata.delta_t_uncertainty_seconds` に供給**できる形（accuracy.md §0）。

## 非目的
- ERA・恒星時の計算（UT1 を消費する側 = umbra-ephemeris）。本 issue は UT1 値と極運動値の供給まで。
- ITRS 変換そのもの（極運動を使う側 = ephemeris フレーム連鎖）。極運動値の提供までが責務。
- IERS Bulletin A 短期予測の自前生成（取り込むのは C04 系列＋ Espenak–Meeus 外挿。Bulletin A は将来短期の選択肢として記録に留め、実装はレビューで判断）。

## 公開インターフェース
api-draft §2 末尾を転記・具体化:
```rust
pub trait DeltaTModel: Send + Sync {
    fn delta_t_seconds(&self, utc: UtcInstant) -> f64;        // TT − UT1
    fn uncertainty_seconds(&self, utc: UtcInstant) -> f64;     // accuracy.md §0 不確実性帯
}
pub trait EarthOrientation: Send + Sync {
    fn ut1_minus_utc(&self, utc: UtcInstant) -> Result<f64, TimeError>;            // 秒
    fn polar_motion(&self, utc: UtcInstant) -> Result<(Radians, Radians), TimeError>; // (xp, yp)
}

pub struct IersEopData { /* versioned + checksum, series 版を保持 */ }
impl IersEopData {
    pub fn bundled() -> Self;
    pub fn series_version(&self) -> &str;   // "EOP 14 C04" / "EOP 20 C04"
    pub fn coverage(&self) -> TimeRange<UtcInstant>;
}
impl EarthOrientation for IersEopData { /* 補間 */ }

pub struct EspenakMeeusDeltaT;            // 長期外挿
impl DeltaTModel for EspenakMeeusDeltaT { /* 多項式 + 不確実性 */ }

// 合成: EOP 範囲内は EOP 由来 ΔT（高精度）、範囲外は Espenak–Meeus 外挿（要確認: 合成器を別型で）
// utc_to_ut1 は TimeScales（api-draft §3.2）が EarthOrientation を使って実装
```

## 数式・アルゴリズムの出典
- **ΔT 定義**: ΔT = TT − UT1（conventions §6, IERS Conventions 2010）。UT1 = UTC + (UT1−UTC)、TT = UTC + (TAI−UTC) + 32.184（ISSUE-006）。よって ΔT = (TAI−UTC) + 32.184 − (UT1−UTC)。EOP 範囲内はこの恒等式で高精度に導出（accuracy.md §3.3 ΔT 履歴は EOP 由来）。
- **長期/外挿 ΔT**: **Espenak & Meeus 多項式**（NASA、"Five Millennium Canon of Solar Eclipses" / `eclipse.gsfc.nasa.gov` の ΔT polynomial expressions、区間別多項式 -1999〜+3000）。data-sources §3.3。各区間の係数を出典明記で取り込む。
- **EOP C04**: IERS Earth Orientation Center, `datacenter.iers.org`。系列 = **EOP 14 C04**（ITRF2014 対応・従来運用）と **EOP 20 C04**（ITRF2020 対応・新運用、移行中）。data-sources §3.1。**採用版を `IersEopData.series_version()` と `CalculationMetadata` に固定**。1 日間隔 → 補間（線形 or Lagrange。要確認: IERS 推奨の補間と日内潮汐補正の要否）。
- **不確実性帯**: accuracy.md §0/§2.3 = 「過去/近傍 <0.1 s（IERS 実測）、〜2100 方向 数秒（予測外挿）」。Espenak–Meeus は外挿の標準誤差式（NASA 記載の不確実性）を採用。要確認: 採用する不確実性式の一次出典。

## 単位 / 時刻系 / 座標系
- 入力: `UtcInstant`。
- 出力: ΔT・UT1−UTC は秒（f64）。極運動 (xp, yp) は `Radians`（IERS は秒角 arcsec 配布 → 取り込み境界で rad へ。conventions §1）。`utc_to_ut1` は `Ut1Instant`。
- 時刻系: UTC 入力、UT1 出力。ΔT は TT−UT1。
- 座標系: 極運動は ITRS 関連（TIRS→ITRS で使用、conventions §5）。本 issue は値供給のみ。

## アルゴリズム概要
1. `IersEopData`: C04 系列（採用版を固定）を versioned + checksum で読み込み（実行時ネットワーク禁止・accuracy.md §5）。日次テーブル（UT1−UTC, xp, yp）。
2. `ut1_minus_utc` / `polar_motion`: テーブル補間。範囲外（coverage 外）は `TimeError::MissingEarthOrientationData`。極運動は arcsec→`Radians` 変換（境界）。
3. `EspenakMeeusDeltaT::delta_t_seconds`: 該当年区間の多項式を Horner 評価（ISSUE-009 周辺の多項式評価方針）。区間境界の連続性に注意。
4. `uncertainty_seconds`: 近傍（EOP 実測内）は <0.1 s、外挿は年数に応じ増大する不確実性式（NASA）。EOP 範囲内/外で切替。
5. ΔT 合成器（推奨・別型）: EOP coverage 内は恒等式（手順の出典欄）で高精度、coverage 外は Espenak–Meeus。切替境界の不連続を許容内に。
6. `utc_to_ut1`: `UT1 = UTC + (UT1−UTC)`（`JulianDate2::add_seconds`）。
- 数値安定性: 補間は端点外挿を弾く。多項式は Horner。禁止: ΔT の magic 値直書き、生 f64 時刻、無条件の線形外挿で不確実性を 0 と誤報すること。

## 受け入れテスト
accuracy.md テストレベル **L2（時刻）**。
- 既知値（オラクル＝IERS EOP C04 配布値 / Espenak–Meeus 公表 ΔT。実装からコピーしない）:
  - 近傍年（例 2020）の ΔT が IERS 由来で **±0.1 s 以内**に再現（恒等式経由）。
  - Espenak–Meeus: 1900・2000・2050 の ΔT が NASA 公表多項式値と一致（区間別）。
  - UT1−UTC・極運動: C04 の特定日の値を補間で再現（補間許容内）。
- 不確実性帯: 近傍年で `uncertainty_seconds < 0.1`、2100 方向で**数秒オーダ**に増大すること（accuracy.md §0 の挙動を満たす）。単調性（将来ほど不確実）をプロパティ確認。
- 系列版: `series_version()` が "EOP 14 C04" / "EOP 20 C04" のいずれかを返し、`CalculationMetadata` に伝播すること。両版で coverage と値差が記録できる。
- 境界値: EOP coverage 端、Espenak–Meeus 区間境界（連続性）、1972 前後（EOP→外挿切替）。
- 異常系: coverage 外で `MissingEarthOrientationData`、データ未ロード時のエラー、`valid_to` 期限切れ検知。

## 許容誤差
accuracy.md §0 / §2.3 から直接引く:
- ΔT (TT−UT1): 過去/近傍は **<0.1 s**（IERS 実測由来）。将来（〜2100）は **数秒**（予測外挿）— これは誤差ではなく**不確実性帯として出力必須**。
- UT1−UTC: **<数 ms（EOP）**。δUT1 1 ms ≈ 地球回転 0.46″（赤道）≈ 観測者 ~14 m（accuracy.md §2.3）。補間誤差はこれ以下を目標。
- 極運動 → 地上: **~10 m**（accuracy.md §2.3）。極運動補間誤差はこれに収まること。
- 不確実性帯の出力自体は「欠落させない」が厳格基準（accuracy.md §0: README/非保証事項にも明記）。

## 実装メモ
- **2 系列の明示**: EOP 14 C04 / 20 C04 のどちらを既定にするか（移行中）はレビュー判断。採用版・取得日・checksum を `DataSetMetadata` と `CalculationMetadata` に固定（data-sources §3.1）。両版差は accuracy.md に記録。
- 将来 UTC は ΔT/UT1 予測律速（accuracy.md §0(b)）。接触時刻は TT 一級保持・UTC 併記（conventions §6）。本 issue の不確実性帯が `delta_t_uncertainty_seconds` の源。
- 1972 以前は UTC 階段でない（ISSUE-006 §実装メモ）。この年代は ΔT（Espenak–Meeus）経由で UT1/TT を扱う橋渡しをここで担う想定。境界仕様をレビューで確定。
- データはバージョン + checksum、実行時ネットワーク禁止、更新は xtask 隔離（accuracy.md §5 / architecture §11）。
- レビュー重点: 不確実性帯の増大挙動が accuracy.md §0 と整合するか、EOP↔外挿切替の連続性、極運動 arcsec→rad 変換、系列版の metadata 伝播。
