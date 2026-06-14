# 公開API草案 (api-draft)

`umbra-rs` の公開境界（pub 型・trait・関数）の草案。**型レベルで固めることが目的**で、コンパイルは通さない。
architecture.md / conventions.md / accuracy.md と整合。確定はレビューゲート（仕様→数式→実装）を通す。

> 状態: ドラフト（Milestone 0 最終成果物）。`/* ... */` は内部表現（非公開）。
> 安定化の単位は **crate**。SemVer 上の公開境界は各 crate の `pub` のみ。

---

## 0. 公開境界の原則

- 公開APIで生 `f64` を「単位付き量」として渡さない（薄い newtype）。conventions §1。
- 時刻は型で時刻系を持つ。フレームは型/名前で判別可能。conventions §5/§6。
- 「該当日食なし」はエラーにせず `Result<Option<_>, EclipseError>`。
- 結果型は必ず `CalculationMetadata` を伴う（accuracy.md §0）。
- 破壊的変更を避けるため、列挙型は `#[non_exhaustive]`、設定型は builder か `Default` + フィールド更新で前方互換を確保。
- `serde` は feature ゲート（`features = ["serde"]`）。コア計算は serde 非依存。

---

## 1. umbra-core — 基盤型

### 1.1 量（newtype）
```rust
#[derive(Clone, Copy, Debug, PartialEq, PartialOrd)] pub struct Radians(pub f64);
#[derive(Clone, Copy, Debug, PartialEq, PartialOrd)] pub struct Degrees(pub f64);
#[derive(Clone, Copy, Debug, PartialEq, PartialOrd)] pub struct Kilometers(pub f64);
#[derive(Clone, Copy, Debug, PartialEq, PartialOrd)] pub struct Meters(pub f64);
#[derive(Clone, Copy, Debug, PartialEq, PartialOrd)] pub struct AstronomicalUnits(pub f64);

impl Radians { pub fn to_degrees(self) -> Degrees; pub fn normalized_signed(self) -> Self; pub fn normalized_two_pi(self) -> Self; }
impl Degrees { pub fn to_radians(self) -> Radians; }
// 角度の正規化は用途別（signed [-π,π) / two_pi [0,2π)）。conventions §2。
```

### 1.2 緯度経度（観測者用の安全な構築）
```rust
#[derive(Clone, Copy, Debug)] pub struct GeodeticLatitude(/* Radians, [-π/2, π/2] */);
#[derive(Clone, Copy, Debug)] pub struct EastLongitude(/* Radians, [-π, π) */);

impl GeodeticLatitude { pub fn from_degrees(deg: f64) -> Result<Self, DomainError>; pub fn degrees(self) -> f64; }
impl EastLongitude   { pub fn from_degrees(deg: f64) -> Result<Self, DomainError>; pub fn degrees(self) -> f64;
                       /// 西経入力（負＝西）も受理し東経正へ正規化
                       pub fn from_signed_degrees(deg: f64) -> Result<Self, DomainError>; }
```

### 1.3 時刻
```rust
#[derive(Clone, Copy, Debug, PartialEq, PartialOrd)] pub struct JulianDate2 { pub part1: f64, pub part2: f64 }

#[derive(Clone, Copy, Debug, PartialEq, PartialOrd)] pub struct UtcInstant(/* ... */);
#[derive(Clone, Copy, Debug, PartialEq, PartialOrd)] pub struct TaiInstant(/* ... */);
#[derive(Clone, Copy, Debug, PartialEq, PartialOrd)] pub struct TtInstant(/* ... */);
#[derive(Clone, Copy, Debug, PartialEq, PartialOrd)] pub struct Ut1Instant(/* ... */);
#[derive(Clone, Copy, Debug, PartialEq, PartialOrd)] pub struct TdbInstant(/* ... */);

impl UtcInstant {
    pub fn from_gregorian(y: i32, mo: u8, d: u8, h: u8, mi: u8, s: f64) -> Result<Self, TimeError>;
    pub fn from_rfc3339(s: &str) -> Result<Self, TimeError>;
    pub fn to_gregorian(self) -> (i32, u8, u8, u8, u8, f64);
}
// UTC↔TAI↔TT は閏秒テーブル、UT1 は EOP（ΔT/UT1-UTC）が必要 → 変換は TimeScales 経由（§3.2、Result<_, TimeError>）。

#[derive(Clone, Copy, Debug)] pub struct TimeRange<T> { pub start: T, pub end: T }
#[derive(Clone, Copy, Debug)] pub struct TimeInterval<T> { pub start: T, pub end: T }
pub type UtcRange = TimeRange<UtcInstant>;
```

### 1.4 線形代数・座標（公開は最小限）
```rust
#[derive(Clone, Copy, Debug)] pub struct Vector3 { pub x: f64, pub y: f64, pub z: f64 }
#[derive(Clone, Copy, Debug)] pub struct UnitVector3(/* 正規化済み */);
// Matrix3 は内部利用中心。必要分のみ pub。

pub trait ReferenceFrame {}
pub enum Icrs {} pub enum Gcrs {} pub enum Cirs {} pub enum Tirs {} pub enum Itrs {} pub enum FundamentalPlane {}
impl ReferenceFrame for Icrs {} /* …他も */

#[derive(Clone, Copy, Debug)] pub struct Position<F: ReferenceFrame> { pub vector: Vector3, /* _frame: PhantomData<F> */ }
```

### 1.5 地球モデル・観測者
```rust
#[non_exhaustive] #[derive(Clone, Copy, Debug)] pub enum EarthModel { Wgs84 }   // 既定 Wgs84。conventions §4

#[derive(Clone, Copy, Debug)]
pub struct Observer {
    pub latitude: GeodeticLatitude,   // 測地緯度
    pub longitude: EastLongitude,     // 東経正
    pub elevation: Meters,            // 楕円体高
}
impl Observer { pub fn new(lat: GeodeticLatitude, lon: EastLongitude, elevation: Meters) -> Self; }
```

### 1.6 エラー（基盤）
```rust
#[non_exhaustive] #[derive(Debug)] pub enum DomainError { OutOfRange { what: &'static str } }
#[non_exhaustive] #[derive(Debug)] pub enum TimeError { InvalidDate, MissingLeapSecondData, MissingEarthOrientationData }
#[non_exhaustive] #[derive(Debug)] pub enum SolverError { RootNotBracketed, DidNotConverge, NumericalInstability }
// std::error::Error / Display を実装。
```

---

## 2. umbra-ephemeris — 天体暦と見かけ位置

```rust
#[non_exhaustive] #[derive(Clone, Copy, Debug, PartialEq, Eq)] pub enum Body { Sun, Earth, Moon, EarthMoonBarycenter }
#[non_exhaustive] #[derive(Clone, Copy, Debug)] pub enum Origin { SolarSystemBarycenter, Geocenter }
#[non_exhaustive] #[derive(Clone, Copy, Debug)] pub enum EphemerisFrame { Icrs, EclipticOfDate }

#[derive(Clone, Copy, Debug)] pub struct StateVector { pub position: Vector3, pub velocity: Option<Vector3> }

#[derive(Clone, Debug)] pub struct EphemerisMetadata {
    pub model: String,        // 例 "VSOP87D+ELP/MPP02"
    pub version: String,      // 採用打切り次数・達成残差を含む識別子
    pub source: String, pub license: String,
    pub supported: TimeRange<TdbInstant>,
    pub max_residual_arcsec: f64,   // accuracy.md §2.4 で実測した値
}

pub trait Ephemeris: Send + Sync {
    fn state(&self, body: Body, time: TdbInstant, origin: Origin, frame: EphemerisFrame)
        -> Result<StateVector, EphemerisError>;
    fn supported_range(&self) -> TimeRange<TdbInstant>;
    fn metadata(&self) -> EphemerisMetadata;
}

pub struct AnalyticalEphemeris { /* VSOP87D + ELP/MPP02 生成係数 */ }
impl AnalyticalEphemeris { pub fn new() -> Self; }            // 純Rust・本番標準
impl Default for AnalyticalEphemeris { fn default() -> Self; }

#[cfg(feature = "jpl")] pub struct JplEphemeris { /* SPK reader */ }
#[cfg(feature = "jpl")] impl JplEphemeris { pub fn from_spk_path(path: &std::path::Path) -> Result<Self, EphemerisError>; }

pub struct MockEphemeris { /* 人工配置 */ }
impl MockEphemeris {
    pub fn central_total() -> Self; pub fn clear_annular() -> Self;
    pub fn clear_partial() -> Self; pub fn shadow_misses_earth() -> Self;
}

// 見かけ地心位置の補正
#[derive(Clone, Copy, Debug)]
pub struct AstrometryOptions {
    pub light_time: bool, pub aberration: bool,
    pub precession_nutation: bool, pub relativistic_deflection: bool,
}
impl AstrometryOptions { pub fn standard() -> Self; pub fn fast() -> Self; }  // standard は前3つ true

#[non_exhaustive] #[derive(Debug)] pub enum EphemerisError {
    OutOfSupportedRange, DataUnavailable, Io(/* ... */),
}

// ΔT / UT1 / 閏秒
pub trait DeltaTModel: Send + Sync {
    fn delta_t_seconds(&self, utc: UtcInstant) -> f64;            // TT - UT1
    fn uncertainty_seconds(&self, utc: UtcInstant) -> f64;        // accuracy.md §0
}
pub trait EarthOrientation: Send + Sync {
    fn ut1_minus_utc(&self, utc: UtcInstant) -> Result<f64, TimeError>;
    fn polar_motion(&self, utc: UtcInstant) -> Result<(Radians, Radians), TimeError>;  // (xp, yp)
}
pub struct IersEopData { /* versioned + checksum */ }
pub struct EspenakMeeusDeltaT;   // 長期外挿
```

---

## 3. umbra-eclipse — 日食エンジン（中核公開API）

### 3.1 精度・設定
```rust
#[non_exhaustive] #[derive(Clone, Copy, Debug, PartialEq, Eq)] pub enum AccuracyProfile { Fast, Standard, Reference }

#[non_exhaustive] #[derive(Clone, Copy, Debug)] pub enum LunarRadiusModel {
    IauMean,            // k = 0.2725076（既定）
    EspenakUmbral,      // k = 0.272281
    EspenakPenumbral,   // k = 0.2725076
}
#[non_exhaustive] #[derive(Clone, Copy, Debug)] pub enum SolarRadiusModel { Iau2015 /* 696000 km */ }
#[non_exhaustive] #[derive(Clone, Copy, Debug)] pub enum RefractionModel { None, Standard }

#[derive(Clone, Copy, Debug)]
pub struct EngineConfig {
    pub accuracy: AccuracyProfile,
    pub earth_model: EarthModel,
    pub lunar_radius_model: LunarRadiusModel,
    pub solar_radius_model: SolarRadiusModel,
    pub refraction: RefractionModel,
    pub root_tolerance_seconds: f64,        // 目標の 1/10 以下
    pub path_sample_interval_seconds: f64,
}
impl EngineConfig { pub fn standard() -> Self; pub fn fast() -> Self; pub fn reference() -> Self; }
impl Default for EngineConfig { fn default() -> Self; }   // = standard()
```

### 3.2 エンジン
```rust
/// 時刻系・地球姿勢に必要なデータ（閏秒 + EOP + ΔT）を**束ねる**純粋データ型（B3）。
/// umbra-core 所属（純粋／no_std 互換、データ注入）。変換ロジックは持たず、
/// 変換 facade `TimeScales` の構築材料になる。旧 api-draft A2「TimeScales を TimeData に集約」は撤回。
pub struct TimeData { /* leap: LeapSecondTable, eop: IersEopData, dt: Box<dyn DeltaTModel> */ }
impl TimeData {
    /// 外部パスのデータで構築（更新は xtask 隔離・実行時ネットワーク禁止）。B3。
    pub fn from_path(dir: &std::path::Path) -> Result<Self, TimeError>;
    /// データの有効範囲。valid_to 超過の問い合わせは TimeScales 変換が Missing*Data を返す。
    pub fn valid_range(&self) -> TimeRange<UtcInstant>;
}
// 同梱データ（埋込 EOP/閏秒/ΔT）での構築 `TimeData::bundled()` は
// **`bundled-data` feature ゲート**で **umbra-ephemeris に置き re-export**（同梱バイトの帰属、B3）。
// umbra-core は no_std 互換の純粋データ注入を保つため、同梱バイトは core に置かない。

/// 時刻系変換の facade（B3）。`TimeData` から構築する。
/// 変換 (utc_to_tt / utc_to_ut1 等) はすべて `Result<_, TimeError>`。
/// umbra-core 所属（純粋／no_std 互換）。
pub struct TimeScales { /* TimeData から構築した変換器 */ }
impl TimeScales {
    /// `TimeData` から変換 facade を構築。
    pub fn new(data: TimeData) -> Self;

    /// 閏秒テーブル不足時は `TimeError::MissingLeapSecondData`。
    pub fn utc_to_tt(&self, t: UtcInstant) -> Result<TtInstant, TimeError>;
    /// valid_to 超過時は `TimeError::MissingEarthOrientationData`。
    pub fn utc_to_ut1(&self, t: UtcInstant) -> Result<Ut1Instant, TimeError>;
    /// 閏秒テーブル不足時は `TimeError::MissingLeapSecondData`。
    pub fn tt_to_utc(&self, t: TtInstant) -> Result<UtcInstant, TimeError>;
    /// 極運動 (xp, yp) も EOP 経由で供給する（CIRS→ITRS の極運動段、conventions §5）。
    /// valid_to 超過時は `TimeError::MissingEarthOrientationData`。
    pub fn polar_motion(&self, t: UtcInstant) -> Result<(Radians, Radians), TimeError>;
}
// 補足: 閏秒は MissingLeapSecondData、EOP/極運動は MissingEarthOrientationData を valid_to 超過時に返す。

pub struct EclipseEngine<E: Ephemeris, D: DeltaTModel, O: EarthOrientation> {
    /* ephemeris: E, delta_t: D, earth_orientation: O, time_scales: TimeScales, config */
    // B3: EclipseEngine は TimeScales（変換 facade）を保持する。
}

impl<E: Ephemeris, D: DeltaTModel, O: EarthOrientation> EclipseEngine<E, D, O> {
    pub fn new(ephemeris: E, delta_t: D, earth_orientation: O, config: EngineConfig) -> Self;

    /// 期間内の太陽食を列挙（偽陰性なし方針。plan Step1-3）
    pub fn search(&self, range: UtcRange) -> Result<Vec<SolarEclipse>, EclipseError>;

    /// 観測地点の局地条件。**既定 Standard は直接瞬時計算 `InstantaneousEvaluator`（037, fit 誤差ゼロ）**（B2）。
    /// 多項式(022)は「M2 で fit 残差<0.01″ を実測実証した後のバッチ局地の任意最適化（未達はフォールバック）」。
    pub fn local_circumstances(&self, eclipse: &SolarEclipse, observer: Observer)
        -> Result<LocalCircumstances, EclipseError>;

    /// 全球経路（中心線・限界線・部分食域・GeoJSON 用）。
    /// v0.1 未実装は panic でなく `Err(EclipseError::NotImplemented)` を返す（PATH、ISSUE-045）。
    pub fn path(&self, eclipse: &SolarEclipse, options: PathOptions)
        -> Result<EclipsePath, EclipseError>;

    /// 指定地点で次に「見える」日食
    pub fn next_visible_eclipse(&self, after: UtcInstant, observer: Observer)
        -> Result<Option<VisibleSolarEclipse>, EclipseError>;

    /// 中間値の検査（CLI inspect / 検証用）
    pub fn instantaneous_elements(&self, time: TtInstant)
        -> Result<InstantaneousBesselianElements, EclipseError>;
}
// 既定構成（A1）: ジェネリック(E,D,O)は維持。`Box<dyn>` は採用しない（性能・単型化優先）。
// 利便性は型エイリアス + ショートカット関数で確保する。
pub type StandardEngine = EclipseEngine<AnalyticalEphemeris, EspenakMeeusDeltaT, IersEopData>;
// B3: TimeData から TimeScales を構築して engine に保持。bundled() は umbra-ephemeris（bundled-data feature）由来。
pub fn standard_engine(time_data: TimeData) -> StandardEngine;   // 例: standard_engine(TimeData::bundled())
```

### 3.3 ベッセル要素と供給源（直接 / 多項式の両対応）
```rust
#[derive(Clone, Copy, Debug)]
pub struct InstantaneousBesselianElements {
    pub time_tt: TtInstant,
    pub x: f64, pub y: f64,                 // 単位 Re
    pub declination: Radians, pub hour_angle: Radians,
    pub l1: f64, pub l2: f64, pub tan_f1: f64, pub tan_f2: f64,
}

pub trait BesselianSource {
    fn at(&self, time: TtInstant) -> Result<InstantaneousBesselianElements, EclipseError>;
    fn fit_interval(&self) -> TimeInterval<TtInstant>;
}
// B2: Standard 局地の**既定供給源**＝直接瞬時計算 `InstantaneousEvaluator`（037, fit 誤差ゼロ）。
pub struct InstantaneousEvaluator</* &engine */> { /* ... */ }   // 既定（直接, 037, fit誤差ゼロ）
impl BesselianSource for InstantaneousEvaluator { /* ... */ }

#[derive(Clone, Debug)]
pub struct BesselianPolynomial {
    pub epoch_tt: TtInstant,
    pub x: Polynomial, pub y: Polynomial, pub d: Polynomial, pub mu: Polynomial,
    pub l1: Polynomial, pub l2: Polynomial, pub tan_f1: f64, pub tan_f2: f64,
    pub fit_interval: TimeInterval<TtInstant>,
    pub fit_error: BesselFitError,
}
// 多項式(022)は経路/GeoJSON/NASA エクスポート用＋「M2 で fit 残差<0.01″ を実測実証した後の
// バッチ局地の**任意最適化**（未達はフォールバック）」。局地の既定供給源ではない（B2）。
impl BesselianSource for BesselianPolynomial { /* ... */ }
#[derive(Clone, Copy, Debug)] pub struct BesselFitError { pub max_x: f64, pub max_y: f64, pub max_l1: f64, pub max_l2: f64 }
#[derive(Clone, Debug)] pub struct Polynomial { pub coefficients: Vec<f64> } // Horner 評価
```

### 3.4 結果型
```rust
#[non_exhaustive] #[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SolarEclipseKind { Partial, Annular, Total, Hybrid, NonCentralAnnular, NonCentralTotal }

#[derive(Clone, Debug)]
pub struct SolarEclipse {                 // search の各要素
    pub event_key: String,                // 安定キー（DB 用）。A4: 最大食 UTC 日付 + lunation 番号
    pub kind: SolarEclipseKind,
    pub global: GlobalCircumstances,
    pub bessel: BesselianPolynomial,      // 経路/エクスポート用
    pub metadata: CalculationMetadata,
}

#[derive(Clone, Debug)]
pub struct GlobalCircumstances {
    pub kind: SolarEclipseKind,
    pub partial_begin: Option<GlobalContact>, pub central_begin: Option<GlobalContact>,
    pub greatest: GreatestEclipse,
    pub central_end: Option<GlobalContact>,   pub partial_end: Option<GlobalContact>,
    pub gamma: f64,
}
#[derive(Clone, Copy, Debug)] pub struct GlobalContact { pub time_utc: UtcInstant, pub time_tt: TtInstant, pub position: GeoPoint }
#[derive(Clone, Copy, Debug)]
pub struct GreatestEclipse {
    pub time_utc: UtcInstant, pub time_tt: TtInstant,
    pub position: GeoPoint, pub magnitude: EclipseMagnitude, pub obscuration: Obscuration,
    pub path_width: Option<Kilometers>, pub central_duration: Option<f64 /* s */>,
    pub sun_altitude: Degrees,
}

/// 局地接触集合（A3）: 部分食地点では c1..c4 が None になりうる。
/// maximum は常に存在（非 Option）— どの地点でも「最大食」は定義される。
#[derive(Clone, Copy, Debug)]
pub struct LocalContactSet {
    pub c1: Option<LocalContact>, pub c2: Option<LocalContact>,
    pub maximum: LocalContact,                 // 非 Option
    pub c3: Option<LocalContact>, pub c4: Option<LocalContact>,
}

#[derive(Clone, Debug)]
pub struct LocalCircumstances {
    pub contacts: LocalContactSet,             // A3: c1..c4=Option, maximum=非Option
    pub magnitude: EclipseMagnitude, pub obscuration: Obscuration,
    pub maximum_altitude: Degrees,
    pub visibility: Visibility,
    pub metadata: CalculationMetadata,
}
#[derive(Clone, Copy, Debug)]
pub struct LocalContact {
    pub time_utc: UtcInstant, pub time_tt: TtInstant,    // accuracy.md §0 両方返す
    pub sun_altitude: Degrees, pub sun_azimuth: Degrees,
    pub position_angle: Degrees, pub visible: bool,
}
#[derive(Clone, Debug)] pub struct VisibleSolarEclipse { pub eclipse: SolarEclipse, pub local: LocalCircumstances }

#[non_exhaustive] #[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Visibility { NotVisible, BelowHorizon, SunriseEclipse, SunsetEclipse, PartialVisible, FullyVisible }

#[derive(Clone, Copy, Debug, PartialEq, PartialOrd)] pub struct EclipseMagnitude(pub f64);  // 皆既で 1 超を許容
#[derive(Clone, Copy, Debug, PartialEq, PartialOrd)] pub struct Obscuration(pub f64);        // 0..1
```

### 3.5 メタデータ・エラー
```rust
#[derive(Clone, Debug)]
pub struct CalculationMetadata {
    pub library_version: String,
    pub ephemeris_model: String, pub ephemeris_version: String,
    pub delta_t_model: String, pub delta_t_uncertainty_seconds: f64,   // 将来UTC律速
    pub earth_model: String, pub lunar_radius_model: String,
    pub accuracy_profile: AccuracyProfile,
    pub generated_at: UtcInstant,
}
impl CalculationMetadata { pub fn fingerprint(&self) -> String; }  // DB 差分再生成用ハッシュ（plan §22）

// A6: thiserror 採用（cargo-deny allow-list に追加）。
// variant 重複を一意化: 下位エラーは From 経由で **ラップ** し、直 variant は
// 日食固有の失敗に限定する（旧 RootNotBracketed/MissingLeapSecondData 等の二重定義を解消）。
#[non_exhaustive] #[derive(Debug, thiserror::Error)]
pub enum EclipseError {
    // 下位エラーは透過ラップ（#[from]）。時刻/EOP/閏秒系は Time、求根系は Solver に集約。
    #[error(transparent)] Time(#[from] TimeError),            // InvalidDate / MissingLeapSecondData / MissingEarthOrientationData
    #[error(transparent)] Ephemeris(#[from] EphemerisError),  // OutOfSupportedRange / DataUnavailable / Io
    #[error(transparent)] Solver(#[from] SolverError),        // RootNotBracketed / DidNotConverge / NumericalInstability
    #[error(transparent)] Domain(#[from] DomainError),        // 観測者/範囲などの定義域違反

    // 日食固有（下位に無い失敗のみ直 variant）
    #[error("degenerate eclipse geometry")] DegenerateGeometry,
    #[error("besselian fit exceeded tolerance")] BesselFitExceededTolerance,
    #[error("solver root not bracketed")] RootNotBracketed,  // 要確認: 下位 SolverError::RootNotBracketed と語義重複しないか（ISSUE-044）
    /// 未実装機能（PATH: path() 等。panic でなくこれを返す。UnsupportedTimeRange を流用しない）。
    #[error("not implemented")] NotImplemented,
}
// 方針(ERR): 「solver が収束しない」「EOP が無い」等は下位 variant をそのまま From ラップし、
// EclipseError 側に同義の直 variant を重複定義しない（透過ラップ版に一本化）。
// 直 variant は日食固有（DegenerateGeometry / BesselFitExceededTolerance / RootNotBracketed / NotImplemented 等、
// 下位と重複しないもの）のみ。未実装は NotImplemented（UnsupportedTimeRange を未実装の意味に流用しない）。
```

---

## 4. umbra-geo — 経路と GeoJSON

```rust
#[derive(Clone, Copy, Debug)] pub struct GeoPoint { pub lat: GeodeticLatitude, pub lon: EastLongitude }
#[derive(Clone, Debug)] pub struct GeoLine { pub points: Vec<GeoPoint> }
#[derive(Clone, Debug)] pub struct GeoPolygon { pub rings: Vec<Vec<GeoPoint>> }

#[derive(Clone, Debug)]
pub struct EclipsePath {
    pub center_line: Option<GeoLine>,
    pub northern_limit: Option<GeoLine>,
    pub southern_limit: Option<GeoLine>,
    pub partial_limit: Option<GeoPolygon>,
    pub greatest_point: GeoPoint,
    pub samples: Vec<PathSample>,    // 各点のプロパティ
}
#[derive(Clone, Copy, Debug)]
pub struct PathSample {
    pub time_utc: UtcInstant, pub center: GeoPoint,
    pub duration_seconds: f64, pub sun_altitude: Degrees,
    pub path_width: Kilometers, pub kind: SolarEclipseKind,
}
#[derive(Clone, Copy, Debug)]
pub struct PathOptions { pub sample_interval_seconds: f64, pub include_limits: bool, pub split_antimeridian: bool }

impl EclipsePath {
    #[cfg(feature = "geojson")] pub fn to_geojson(&self) -> String;  // 日付変更線分割・極域対応
}
```

---

## 5. 使用イメージ（doctest 候補）

```rust
use umbra_eclipse::*;
let engine = standard_engine(TimeData::bundled());           // 同梱 EOP/閏秒/ΔT
let eclipses = engine.search(UtcRange { start: utc(2026,1,1), end: utc(2040,1,1) })?;
let okayama = Observer::new(
    GeodeticLatitude::from_degrees(34.507)?,
    EastLongitude::from_degrees(133.508)?,
    Meters(10.0),
);
if let Some(v) = engine.next_visible_eclipse(utc(2026,1,1), okayama)? {
    println!("max(UTC)={:?} max(TT)={:?} mag={:?} alt={:?}",
        v.local.contacts.maximum.time_utc, v.local.contacts.maximum.time_tt,
        v.local.magnitude, v.local.maximum_altitude);
}
```

---

## 6. 設計判断（確定済み / 残課題）

Milestone 0 独立レビューの確定事項 A1–A7（reviews/milestone0-review.md）を反映。

### 6.1 確定済み

- **A1 ジェネリック維持**: `EclipseEngine<E, D, O>` のジェネリック3パラメータを保つ。`Box<dyn>` は採用しない（性能・単型化優先）。利便性は `StandardEngine` 型エイリアス + `standard_engine(TimeData)` ショートカットで確保（§3.2）。
- **A2（撤回・B3 で改訂）TimeData + TimeScales の 2 型**: 旧 A2「TimeScales を TimeData に集約」は **撤回**。確定 B3 に従い 2 型を維持する。
  - `TimeData`（データ束ね、umbra-core・純粋/no_std 互換・データ注入）: `from_path()`（外部）/ `valid_range()`。
  - `TimeScales`（変換 facade、umbra-core・`TimeData` から `TimeScales::new(data)` で構築）: `utc_to_tt` / `utc_to_ut1` / `tt_to_utc` / `polar_motion` をすべて **`Result<_, TimeError>`** で提供。
  - **`TimeData::bundled()`（同梱バイト）は `bundled-data` feature ゲートで umbra-ephemeris に置き re-export**（同梱バイトの帰属。core は no_std 純粋を保つ）。
  - `EclipseEngine` は `TimeScales` を保持。`standard_engine(TimeData::bundled())` で構築。
  - valid_to 超過は `Missing*Data`（閏秒=MissingLeapSecondData / EOP=MissingEarthOrientationData）を返し metadata に記録。**極運動 (xp, yp) も EOP 経由**で供給（§3.2）。
- **A3 接触集合**: `LocalContactSet { c1..c4: Option, maximum: 非Option }`。`LocalCircumstances.contacts` に反映（§3.4）。部分食地点で c2/c3 が None になる前提を型で表現。
- **A4 キー生成規則**:
  - `event_key` = **最大食 UTC 日付 + lunation 番号**（§3.4 SolarEclipse）。
  - `location_key` = **指定値があればそれ**、無ければ**緯度経度を丸めてハッシュ**した安定キー（局地結果の DB キー、plan §22 と整合）。
- **A5 NonCentral 公開**: `SolarEclipseKind` の NonCentral 系は enum 公開を維持（`#[non_exhaustive]`）。ただし **v0.1 CLI は主に Partial / 中心食**を対象（中心食詳細は M8）。
- **A6 thiserror 採用（透過ラップ版に一本化、ERR）**: エラーは `thiserror` で実装（cargo-deny allow-list 追加）。`EclipseError` は **透過ラップ版に一本化**—下位エラー（Time/Ephemeris/Solver/Domain）は `#[from]`（`#[error(transparent)]`）で **ラップ**し、同義の二重 variant を作らない。直 variant は**日食固有**で下位と重複しないもの（`DegenerateGeometry` / `BesselFitExceededTolerance` / `RootNotBracketed` / **`NotImplemented`**）に限定（§3.5）。未実装（path() 等）は `Err(EclipseError::NotImplemented)`（panic 不使用、`UnsupportedTimeRange` を未実装の意味に流用しない、PATH）。
- **A7 serde 安定表現**: 列挙型は `#[serde(tag = "type")]`、**数値の単位はフィールド名で明示**（例 `central_duration_seconds`、`path_width_km`）。`serde` は feature gate（コア計算は serde 非依存、§0）。

### 6.2 残課題（要確認）

- [ ] `location_key` の丸め桁・ハッシュ関数の具体（DB 衝突率と再現性のトレードオフ、plan §22 と詰める）。
- [ ] serde スキーマの版管理（タグ名・フィールド名の SemVer 上の扱い）。一次的なフィールド名規約は A7 で確定済みだが、互換ポリシは**要確認**。
- [ ] `relativistic_deflection` を Standard で既定 ON にするかの上限テスト（accuracy M1）— API 形は確定、既定値は実測後に確定。

---

## 関連
architecture.md（レイヤー/抽象） / conventions.md（単位・符号） / accuracy.md（精度・誤差バジェット） / data-sources.md（出典・ライセンス）。
