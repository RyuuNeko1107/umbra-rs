# アーキテクチャ (architecture)

`umbra-rs` は、時刻系・天体暦・見かけ位置・影円錐・ベッセル要素・全球/局地予報・経路生成・精度検証を、
**一つの整合した規約（conventions.md）で接続する純Rustの日食予報エンジン**である。

> 状態: ドラフト（Milestone 0）。型シグネチャは草案で、実装時にレビューゲートを通す。
> 関連: conventions.md（規約） / accuracy.md（精度・誤差バジェット）。

---

## 1. クレート構成

初期から Cargo workspace。細分化しすぎない。

```
umbra-rs/
├─ crates/
│  ├─ umbra-core/        時刻・角度・距離・ベクトル/行列・座標系・数値解法・定数・エラー
│  ├─ umbra-ephemeris/   太陽/月位置・歳差章動・光行差・地球回転・Ephemeris trait・係数
│  ├─ umbra-eclipse/     新月探索・候補判定・影幾何・ベッセル要素・全球/局地・接触・食分
│  ├─ umbra-geo/         中心線・限界線・部分食域・GeoJSON・経路サンプリング
│  ├─ umbra-cli/         search/local/path/bessel/inspect/validate/bench
│  └─ umbra-fixtures/    NASA/USNO/JPL 比較値・既知ベッセル・テスト地点・許容誤差（通常依存に含めない）
├─ data/                 coefficients / leap-seconds / earth-orientation / test-vectors
├─ docs/                 architecture / conventions / accuracy / algorithms / data-sources / validation
├─ tests/                regression / property / differential
└─ examples/             search / local / path / bessel
```

- `umbra-core` に天文暦・日食固有処理を入れない（純粋基盤）。
- `umbra-fixtures` は通常の依存先に含めない（検証専用）。
- 係数生成・データ更新は `xtask` に隔離。**ライブラリ本体は実行時にネットワークしない**（accuracy.md §5）。

---

## 2. 型設計の原則

- 公開APIで生 `f64` を「単位付き量」として渡さない。最低限の newtype:
  `Radians / Degrees / Kilometers / Meters / AstronomicalUnits / JulianDate`。
  ただし内部の性能・可読性を落とさぬよう **newtype は薄く保つ**。
- 時刻は型で時刻系を持つ: `UtcInstant / TaiInstant / TtInstant / Ut1Instant / TdbInstant`。
  `fn solar_position(jd: f64)` は禁止、`fn solar_position(time: TtInstant)` とする。
- JD は2要素: `JulianDate2 { part1, part2 }`。
- フレームは型で区別: `Position<F: ReferenceFrame>`（`Icrs/Gcrs/Cirs/Tirs/Itrs/FundamentalPlane`）。
  重すぎる箇所は内部限定で可だが、変数名・型名で座標系が判別できることは必須。
- 観測者:
  ```rust
  pub struct Observer {
      pub latitude: GeodeticLatitude,   // 測地緯度
      pub longitude: EastLongitude,     // 東経正 [-180,180)
      pub elevation: Meters,            // 楕円体高
  }
  ```

---

## 3. レイヤーとデータフロー

```
期間(UtcRange)
  → [umbra-ephemeris] 新月概算 → 朔の精解(合solver)
  → [umbra-eclipse]   日食可能性 早期棄却（偽陽性可・偽陰性不可）
  → [umbra-ephemeris] 太陽/月 見かけ地心状態（light-time→歳差章動→光行差）
  → [umbra-eclipse]   月影円錐 → ベッセル基本面 → 瞬時ベッセル要素
  → [umbra-eclipse]   全球条件（P1/U1/最大/U4/P4・種別・gamma・帯幅）
  → [umbra-eclipse]   局地条件（C1/C2/最大/C3/C4・食分・食面積・高度方位・可視性）
  → [umbra-geo]       中心線/限界線/部分食域 → GeoJSON
```

各段の数式・解法の詳細は algorithms.md（別途）。本書は境界とインターフェースを定める。

---

## 4. 天体暦バックエンド

太陽と月を別 trait にせず統合する（モデル整合性のため）。

```rust
pub trait Ephemeris: Send + Sync {
    fn state(&self, body: Body, time: TdbInstant, origin: Origin, frame: EphemerisFrame)
        -> Result<StateVector, EphemerisError>;
    fn supported_range(&self) -> TimeRange;
    fn metadata(&self) -> EphemerisMetadata;
}

pub enum Body { Sun, Earth, Moon, EarthMoonBarycenter }

pub struct StateVector { pub position: Vector3, pub velocity: Option<Vector3> }
```

- 速度を直接持たない解析モデル向けに、速度は **解析微分 / 対称差分 / 補間微分** のいずれかをバックエンドが提供。
- 実装:
  - `AnalyticalEphemeris` … **本番標準**。VSOP87D（太陽）+ ELP/MPP02（月）フル系列。係数を crate 内に保持、ネットワーク不要、再現可能、純Rust、対応期間限定。
  - `JplEphemeris` … feature `jpl`。外部 DE データを読む。Reference / 差分テスト用。本体に巨大データを同梱しない。
  - `MockEphemeris` … 単体テスト用の人工配置（accuracy.md §3.1）。
- 採用項数・達成残差は `EphemerisMetadata`（ephemeris_version）に記録（accuracy.md §2.4）。

---

## 5. 見かけ位置の補正パイプライン

```rust
pub struct AstrometryOptions {
    pub light_time: bool,
    pub aberration: bool,
    pub precession_nutation: bool,
    pub relativistic_deflection: bool,
}
```

- **Standard プロファイルは light_time / aberration / precession_nutation を必須 ON で固定**（accuracy.md §1）。
- 歳差章動は **IAU2006 歳差 + IAU2000A 章動**（Fast のみ IAU2000B 許容）。
- フレーム連鎖（conventions §5）: `GCRS → CIRS → TIRS → ITRS`（ERA は UT1、極運動は EOP）。
- 相対論偏向は初期は省略可。**省略は metadata に残す**。

---

## 6. 日食幾何とベッセル要素

```rust
pub struct ShadowCone {
    pub axis_origin: Vector3, pub axis_direction: UnitVector3,
    pub umbra_apex: Vector3,  pub penumbra_apex: Vector3,
    pub umbra_half_angle: Radians, pub penumbra_half_angle: Radians,
}

pub struct InstantaneousBesselianElements {
    pub time_tt: TtInstant,
    pub x: f64, pub y: f64,           // 単位: 地球赤道半径 Re
    pub declination: Radians, pub hour_angle: Radians,
    pub l1: f64, pub l2: f64,         // 単位: Re
    pub tan_f1: f64, pub tan_f2: f64,
}
```

- 金環食は本影頂点の地球側/反地球側を明示判定。極端配置でも基本面基底が壊れないようにする。
- x, y, l1, l2 の単位は Re（conventions §1）。NASA 表記との対応は本書/algorithms に文書化。

### 6.1 BesselianSource 抽象（直接 vs 多項式の両対応）

精度最優先のため、局地 solver は瞬時要素の**供給源にジェネリック**にする。

```rust
pub trait BesselianSource {
    /// 任意 TT における瞬時ベッセル要素を返す
    fn at(&self, time: TtInstant) -> Result<InstantaneousBesselianElements, EclipseError>;
    fn fit_interval(&self) -> TimeInterval<TtInstant>;
}

/// 各時刻で暦を直接再評価（fit 誤差ゼロ）。**Standard 局地の既定供給源**（B2）。多項式ゲート未達時のフォールバックも兼ねる。
pub struct InstantaneousEvaluator<E: Ephemeris> { /* engine refs */ }

/// 最大食付近を多項式近似。経路/GeoJSON/NASA エクスポート用・高速。局地転用は M2 で残差 <0.01″ 実証後の任意最適化（B2）。
pub struct PolynomialModel { /* coefficients + residuals */ }
```

```rust
pub struct BesselianPolynomial {
    pub epoch_tt: TtInstant,
    pub x: Polynomial, pub y: Polynomial, pub d: Polynomial, pub mu: Polynomial,
    pub l1: Polynomial, pub l2: Polynomial, pub tan_f1: f64, pub tan_f2: f64,
    pub fit_interval: TimeInterval<TtInstant>,
    pub fit_error: BesselFitError,   // 残差を必ず保持
}
pub struct BesselFitError { pub max_x: f64, pub max_y: f64, pub max_l1: f64, pub max_l2: f64 }
```

- 既定: **Standard 局地計算 = `InstantaneousEvaluator`（直接, fit 誤差ゼロ, 037, B2）**。精度最優先のため各地点で暦を直接再評価する。
- `PolynomialModel`（多項式, 022）は **経路/GeoJSON/NASA エクスポート用が本務**。**局地への転用は「M2 で fit 残差 <0.01″ を実測実証後のバッチ局地（多地点）の任意最適化」に限定**し、未達は直接（037）へフォールバック（operations.md §D1）。実証前は局地で多項式を既定にしない。
- accuracy.md L7 サブテストで「直接 vs 多項式」残差を実測し、局地転用の採否（ゲート閾値の妥当性）を裏取りする。
- 多項式次数は固定せず、NASA 形式に合わせた低次から開始。fit_error は許容内をガード。

---

## 7. 全球・局地・可視性

```rust
pub struct GlobalCircumstances {
    pub kind: SolarEclipseKind,
    pub partial_begin: Option<GlobalContact>, pub central_begin: Option<GlobalContact>,
    pub greatest: GreatestEclipse,
    pub central_end: Option<GlobalContact>,   pub partial_end: Option<GlobalContact>,
}

pub struct LocalContact {
    pub time_utc: UtcInstant, pub time_tt: TtInstant,   // accuracy.md §0: 両方返す
    pub sun_altitude: Degrees, pub sun_azimuth: Degrees,
    pub position_angle: Degrees, pub visible: bool,
}

pub enum Visibility {
    NotVisible, BelowHorizon, SunriseEclipse, SunsetEclipse, PartialVisible, FullyVisible,
}

pub struct EclipseMagnitude(pub f64);   // 部分0..1、皆既で1超を許容
pub struct Obscuration(pub f64);        // 0..1、境界条件を明示処理
```

- 接触時刻は UTC と TT を併記（将来 UTC は ΔT/UT1 律速、accuracy.md §0/§2.3）。
- 大気差は `RefractionModel { None, Standard }` で分離。通知用途は補正前後の高度を両方返せる。

---

## 8. エンジン API

```rust
pub struct EclipseEngine<E, D, O> {
    ephemeris: E, delta_t: D, earth_orientation: O, config: EngineConfig,
}

pub struct EngineConfig {
    pub accuracy: AccuracyProfile,
    pub earth_model: EarthModel,                 // 既定 WGS84
    pub lunar_radius_model: LunarRadiusModel,    // 既定 IauMean(k=0.2725076)
    pub solar_radius_model: SolarRadiusModel,    // 既定 IAU2015 696000km
    pub root_tolerance_seconds: f64,             // 目標の1/10以下
    pub path_sample_interval_seconds: f64,
}

impl<E, D, O> EclipseEngine<E, D, O> {
    pub fn search(&self, range: UtcRange) -> Result<Vec<SolarEclipse>, EclipseError>;
    pub fn local_circumstances(&self, eclipse: &SolarEclipse, observer: Observer)
        -> Result<LocalCircumstances, EclipseError>;
    pub fn path(&self, eclipse: &SolarEclipse, options: PathOptions)
        -> Result<EclipsePath, EclipseError>;
    pub fn next_visible_eclipse(&self, after: UtcInstant, observer: Observer)
        -> Result<Option<VisibleSolarEclipse>, EclipseError>;
}
```

「該当日食なし」はエラーにせず `Result<Option<_>, EclipseError>`。

---

## 9. 計算メタデータ（結果に必ず付与）

```rust
pub struct CalculationMetadata {
    pub library_version: String,
    pub ephemeris_model: String, pub ephemeris_version: String,
    pub delta_t_model: String,
    pub delta_t_uncertainty_seconds: f64,   // accuracy.md §0/§2.3: 将来UTCの不確実性帯
    pub earth_model: String,
    pub lunar_radius_model: String,
    pub accuracy_profile: AccuracyProfile,
    pub generated_at: UtcInstant,
}
```

- DB 保存時も保持（fingerprint 化）。モデル更新後に古い結果と混ざらないようにする（plan §22）。

---

## 10. エラー設計

```rust
pub enum EclipseError {
    InvalidDate, InvalidObserver, UnsupportedTimeRange,
    EphemerisUnavailable, MissingLeapSecondData, MissingEarthOrientationData,
    RootNotBracketed, SolverDidNotConverge, DegenerateGeometry,
    BesselFitExceededTolerance, NumericalInstability,
}
```

---

## 11. データ管理

- コードと係数を分離。大量係数を手書き Rust 配列に直書きしない。
  `data/coefficient-source/ → (生成ツール) generated/ → crate 組込み`。
- 再現性: `cargo xtask generate-coefficients` / `cargo xtask verify-generated`。CI で生成済みとの差分確認。
- 各データセットに `DataSetMetadata`（name/version/source/license/valid_from/valid_to/checksum）。
- EOP/ΔT データも同様に versioned + checksum（accuracy.md §5）。

---

## 12. 数値解法（algorithms.md に詳細）

- 求根: **Brent 法**（導関数不要・二分法の堅牢性）。Newton 単独は禁止（conventions §11）。
- 最小化（最大食）: Brent 1次元最小化 or 黄金分割。初期は黄金分割可。
- 数値微分: 中心差分 / 5点差分。差分幅は固定せずテストで決める。
- 多項式評価: Horner 法。
- 合 solver: 粗走査 → 符号変化/極小区間検出 → Brent。角度は ±π 折返しを除いた連続関数へ。

---

## 13. レビューゲート（plan §25）

各層を通過しない限り先へ進まない: 仕様→数式→実装→単体テスト→基準値比較→誤差評価。
誤差は必ず層ごとに分解（accuracy.md §4）。天体位置が不正確なまま日食側で誤差を打ち消さない。
