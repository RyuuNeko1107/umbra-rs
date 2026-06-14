# ISSUE-043: EclipseEngine 組立・オーケストレーション（search/local/next_visible/instantaneous の結線）

- crate: umbra-eclipse
- 依存: ISSUE-012（`Ephemeris`）, ISSUE-007（`DeltaTModel`/`EarthOrientation`）, ISSUE-042（`TimeData`/`TimeScales`・standard_engine 前提）, ISSUE-016/017（新月候補・合 solver）, ISSUE-018（候補棄却フィルタ）, ISSUE-015（見かけ地心位置）, ISSUE-019/020（影円錐・基本面基底）, ISSUE-021/037（瞬時ベッセル要素・供給源）, ISSUE-022（ベッセル多項式 = SolarEclipse.bessel）, ISSUE-023（全球分類・gamma）, ISSUE-024/025/026/027/028（局地射影・接触・最大食・食分・高度方位可視性）, ISSUE-013/014（解析暦 = standard_engine の E）, ISSUE-044（`EclipseError` 集約）, ISSUE-041（定数）
- milestone: M9（全球・局地が出揃った後の最終結線・オーケストレーション。中核公開 API）
- モード(tdd-workflow): **strict**（**確定A1**: ジェネリック(E,D,O) 維持・StandardEngine エイリアス・公開 API 結線は SemVer 境界の中核。api-draft §3.2 を確定するため strict）

## 目的
**確定A1**（milestone0-review.md）に従い、`EclipseEngine<E,D,O>` を組み立て、`search` / `local_circumstances` / `next_visible_eclipse` / `instantaneous_elements` のパイプラインを結線する（architecture §3/§8, api-draft §3.2）。
- パイプライン結線: **候補生成（016）→ 合 solver（017）→ 早期棄却（018）→ 見かけ暦（015）→ 影円錐/基本面（019/020）→ 瞬時ベッセル要素（021/037）→ ベッセル多項式 fit（022）→ 全球分類/gamma（023）→ 局地（024-028）**。
- **確定A1**: ジェネリック `EclipseEngine<E: Ephemeris, D: DeltaTModel, O: EarthOrientation>` を維持（`Box<dyn>` 不使用）。`StandardEngine` 型エイリアス＋`standard_engine(TimeData)` ショートカットで利便確保。
- **v0.1 で `path()` は `Err(EclipseError::NotImplemented)`**（PATH/確定・レビュー minor 045: 経路 = umbra-geo はスコープ外、結果型 bessel 多項式は持つが path 本体は未実装。panic/`unimplemented!` でなく `Err(NotImplemented)`。`UnsupportedTimeRange` を未実装の意味に流用しない）。

## 非目的
- 各段の数式実装（016-028/015/019-022。本 Issue は**結線とオーケストレーション**のみ、計算は委譲）。
- 経路生成 `path()` の中身（ISSUE-045/umbra-geo。v0.1 は `Err(EclipseError::NotImplemented)`、結果型 `BesselianPolynomial` までは保持）。
- 時刻系データの取り込み（ISSUE-042。本 Issue は `TimeScales` を保持・利用）。
- CLI（ISSUE-031/032。エンジンを呼ぶ側）。
- エラー型定義（ISSUE-044。本 Issue は `?` で集約エラーへ変換する利用側）。

## 公開インターフェース
api-draft §3.2 を転記・**確定A1** で具体化:
```rust
pub struct EclipseEngine<E: Ephemeris, D: DeltaTModel, O: EarthOrientation> {
    /* ephemeris: E, delta_t: D, earth_orientation: O, time_scales: TimeScales, config: EngineConfig */
    // B3: engine は TimeScales（変換 facade）を保持する。
}
impl<E: Ephemeris, D: DeltaTModel, O: EarthOrientation> EclipseEngine<E, D, O> {
    /// B3: 引数の `TimeData` から `TimeScales::new(time)` を構築して保持する。
    pub fn new(ephemeris: E, delta_t: D, earth_orientation: O, time: TimeData, config: EngineConfig) -> Self;

    pub fn search(&self, range: UtcRange) -> Result<Vec<SolarEclipse>, EclipseError>;
    pub fn local_circumstances(&self, eclipse: &SolarEclipse, observer: Observer)
        -> Result<LocalCircumstances, EclipseError>;
    pub fn next_visible_eclipse(&self, after: UtcInstant, observer: Observer)
        -> Result<Option<VisibleSolarEclipse>, EclipseError>;
    pub fn instantaneous_elements(&self, time: TtInstant)
        -> Result<InstantaneousBesselianElements, EclipseError>;
    /// v0.1 未実装（PATH/レビュー minor 045: 経路は umbra-geo・スコープ外）
    pub fn path(&self, eclipse: &SolarEclipse, options: PathOptions)
        -> Result<EclipsePath, EclipseError>;   // v0.1: Err(EclipseError::NotImplemented)
}

// 確定A1: 型エイリアス + ショートカット（dyn は使わない）
pub type StandardEngine = EclipseEngine<AnalyticalEphemeris, EspenakMeeusDeltaT, IersEopData>;
// B3: TimeData から TimeScales を構築して保持。bundled() は umbra-ephemeris（bundled-data feature）由来。
pub fn standard_engine(time: TimeData) -> StandardEngine;   // 既定 EngineConfig::standard()、例: standard_engine(TimeData::bundled())
```
- 「該当日食なし」はエラーにせず `Result<Option<_>, _>`（`next_visible_eclipse`）。結果型は `CalculationMetadata` を伴う（accuracy.md §0）。
- B3 確定: `new` の引数に `TimeData` を取り、内部で `TimeScales::new(time)` を構築して保持する（旧 A2「TimeScales を TimeData に集約」は撤回。TimeData=データ束ね／TimeScales=変換 facade の 2 型、変換は `Result<_, TimeError>`、api-draft §3.2）。

## 数式・アルゴリズムの出典
- 本 Issue は**結線**で固有数式なし。各段の出典は委譲先 Issue（016-028 等）。
- パイプライン定義: architecture §3（データフロー）/ §8（エンジン API）。偽陽性可・偽陰性不可の早期棄却方針（architecture §3, ISSUE-018, 確定 D6 マージン）。
- ベッセル多項式 fit を `SolarEclipse.bessel` に保持（ISSUE-022、経路/エクスポート用）。局地は既定で直接瞬時計算 `InstantaneousEvaluator`（ISSUE-037、fit 誤差ゼロ、architecture §6.1）。
- 最大食は dm/dt=0 求根（確定 D2、ISSUE-026）、m²=u²+v² 統一（D2）。本 Issue は solver を正しい順序で呼ぶ。

## 単位 / 時刻系 / 座標系
- 単位/座標は各段に委譲（ベッセル要素 Re、角度 rad 等。conventions §1/§5）。
- 時刻系: 入力 `UtcRange`/`UtcInstant`、内部は `TimeScales`（ISSUE-042）で TT/UT1 へ。接触は UTC と TT 両方返す（conventions §6, accuracy.md §0）。`instantaneous_elements` は `TtInstant` 入力。
- 座標系: パイプライン全体（GCRS→…→FundamentalPlane→局地）を結線（conventions §5）。

## アルゴリズム概要
1. `new`/`standard_engine`: E,D,O と `TimeData`→`TimeScales`、`EngineConfig` を保持（A1 ジェネリック維持）。
2. `search(range)`:
   候補生成（016）→ 各候補で合 solver（017）→ 早期棄却（018、偽陰性不可）→ 残候補で見かけ暦（015）→ 影円錐/基本面（019/020）→ 瞬時要素（021/037）→ ベッセル多項式 fit（022）→ 全球分類/gamma/接触（023）→ `SolarEclipse`（event_key=最大食 UTC 日付+lunation, A4）を構成・metadata 付与。
3. `local_circumstances(eclipse, observer)`:
   `InstantaneousEvaluator`（037）を供給源に、観測者射影（024）→ C1-C4 接触（025）→ 最大食（026, dm/dt=0）→ 食分/食面積（027）→ 高度方位/可視性（028）→ `LocalCircumstances`（接触は UTC/TT 両方、A3 LocalContactSet）。
4. `next_visible_eclipse(after, observer)`:
   `after` 以降を search で走査し、各 eclipse の `local_circumstances` を評価、可視（Visibility が見える種別）になる最初を `Some`、無ければ `None`。
5. `instantaneous_elements(time_tt)`: 021/037 を 1 時刻評価（CLI inspect/検証用）。
6. `path`: **v0.1 `Err(EclipseError::NotImplemented)`**（PATH/045/umbra-geo）。doc に v0.1 非対応・将来 umbra-geo 実装を明記。
- 数値安定性: solver は Brent（無条件 Newton 禁止、conventions §11, ISSUE-008）。角度連続化は各 solver（conventions §2）。誤差は層分解で保持（accuracy.md §4、誤差を日食側で打ち消さない conventions §11）。

## 受け入れテスト
accuracy.md テストレベル **L5（全球）＋ L6（局地）の結線・統合**。基準は fixtures/Mock（実装コピー禁止、conventions §11）。
- **MockEphemeris 統合（ISSUE-038、DE 無し CI）**: `central_total`/`clear_annular`/`clear_partial`/`shadow_misses_earth` をエンジンに通し、`search` の種別・gamma・有無が設計オラクルと一致。`shadow_misses_earth` は search で当該日食なし（または非該当）。
- **A1 構成**: `StandardEngine` 型エイリアス・`standard_engine(TimeData::bundled())` がコンパイル・動作（dyn 不使用、ジェネリック単相化。`bundled()` は umbra-ephemeris の bundled-data feature 由来、B3）。
- **next_visible_eclipse**: 該当なしで `Ok(None)`（エラーにしない、api-draft §0）。api-draft §5 例（岡山）でスモーク。
- **path v0.1**: `path()` 呼び出しが `Err(EclipseError::NotImplemented)` を返すこと（panic でなく Result）を `assert!(matches!(.., Err(EclipseError::NotImplemented)))` で固定し、v0.1 スコープ外を明示（PATH/045）。
- **接触 UTC/TT 併記**: `LocalContact`/`GlobalContact` が time_utc と time_tt 両方を持つ（accuracy.md §0）。
- **ゴールデン（ISSUE-029）連携**: ゴールデン 20 / 1900–2100 一括をエンジン経由で回し誤差統計（ISSUE-030）を生成できる結線。
- **層分解**: DE 差分（036）と解析暦で同一パイプラインを通し誤差を層帰属できる（accuracy.md §3.1/§4、nightly）。

## 許容誤差
- 本 Issue は結線で固有許容なし。**最終目標を委譲先で満たす**: 最大食 ±1.5s（幾何, accuracy.md §2.1）、局地接触 ±2s（幾何分、§1）、食分 ±0.0005（§2.2）。
- solver 収束は `root_tolerance_seconds` を目標の 1/10 以下（EngineConfig、accuracy.md §2.1 solver 0.05″ 配分）。
- 接触 UTC は ΔT/UT1 予測律速（accuracy.md §0(b)、不確実性帯を metadata 出力）。
- 誤差は層分解で保持（accuracy.md §4）。日食側で打ち消さない（conventions §11）。

## 実装メモ
- **確定A1 厳守**: ジェネリック(E,D,O) 維持・`Box<dyn>` 不使用、`StandardEngine` 型エイリアス + `standard_engine()`。
- **v0.1 path() = `Err(EclipseError::NotImplemented)`**（PATH/レビュー minor 045）。doc・テスト（matches! で Err 固定）で明示。panic/`unimplemented!` 不使用、`UnsupportedTimeRange` を未実装に流用しない。経路は umbra-geo（M11 想定）。
- event_key=最大食 UTC 日付+lunation 番号、location_key=指定 or 緯度経度丸めハッシュ（確定 A4、ISSUE-029/DB plan §22 と整合）。
- 局地既定は直接瞬時計算（037、fit 誤差ゼロ）、`SolarEclipse.bessel` 多項式（022）は経路/エクスポート用（architecture §6.1）。
- 早期棄却は偽陰性不可（018、確定 D6 マージン）。最大食は dm/dt=0・m² 統一（確定 D2）。μ は CIO 統一（確定 D4、039）。
- エラーは ISSUE-044 集約（thiserror、透過ラップ版 ERR）へ `?` 変換（`#[from]` で TimeError/EphemerisError/SolverError/DomainError をラップ）。
- レビュー重点: A1 ジェネリック維持、path = `Err(NotImplemented)`、結線順序（候補→合→棄却→暦→影→ベッセル→分類→局地）、UTC/TT 併記、誤差層分解、Mock CI 経路、TimeData 引数→TimeScales 保持（B3 整合）。
