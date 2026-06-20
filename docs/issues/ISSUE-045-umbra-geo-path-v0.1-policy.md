# ISSUE-045: umbra-geo / path の v0.1 方針（型・境界のみ定義／本実装は Milestone 9）

- crate: umbra-geo
- 依存: ISSUE-022（`BesselianPolynomial`＝経路計算の供給源・公開型）, ISSUE-023（`GlobalCircumstances`・種別・最大食地点）, ISSUE-021（瞬時ベッセル要素 x,y,d,μ,l1,l2,tan f）, ISSUE-002（角度・緯度経度 newtype）, ISSUE-010（WGS84・測地/地心緯度）, ISSUE-044（`EclipseError::NotImplemented` 等のエラー集約）, ISSUE-001（規約）
- milestone: M9 経路（v0.1 では型・境界のみ。本実装は Milestone 9）
- モード(tdd-workflow): standard（v0.1 では公開型シグネチャと未実装スタブの契約のみを固定する。`SolarEclipse.bessel`（`BesselianPolynomial`＝ISSUE-022）は v0.1 で必須だが、`path()`/中心線/限界線/GeoJSON の数式本体は Milestone 9 のため、型境界の前方互換性確保が要点で standard。本実装時に strict へ昇格）

## M9 実装状況
- **M9.1 中心線トラック**（2026-06-21・strict）: `EclipseEngine::path()` を実装（旧 `Err(NotImplemented)` スタブから昇格）。中心食（全球 U1/U4 接触が両方 Some）で `center_line` を生成＝`[U1,U4]` を `PathOptions::sample_interval_seconds` 刻みでサンプルし、各時刻のベッセル要素（`BesselianSource::at`）から影軸地表貫通点（`axis_intercept::shadow_axis_surface_point`・WGS84）を結んだ `GeoLine`。軸が地球を外す端（`RootNotBracketed`）はスキップ。非中心は `center_line=None`。`greatest_point` は `global.greatest.position` passthrough。**北/南限界線・部分食域・`samples`（帯幅/継続）・GeoJSON は未実装**（後続スライス）。5 テスト（FAST 4＋SLOW 1: 実 2017-08-21 皆既で中心線が太平洋〜大西洋を延び北米を横断・最大食点近傍を通ることを実証）。mutation 12 中 8 caught・2 unviable・2 timeout（ループ終端変異＝ハング検出）・生存0。
- **残（後続スライス）**: (2) 帯幅・中心食継続（`samples`・算法 §8.11/8.12＝**要一次資料確認**・最大食点の `GreatestEclipse.path_width`/`central_duration` も現状 None）、(3) 北/南限界線・部分食域（錐縁追跡）、(4) `EclipsePath::to_geojson()`（umbra-geo・日付変更線分割）。

## 目的
`umbra-geo` の経路 API（中心線・限界線・部分食域・GeoJSON）の **公開型と境界のみを v0.1 で確定**し、**本実装を Milestone 9 へ明示的に後回し**する方針を文書化する（レビュー minor 確定事項 / milestone0-review §Minor「045 umbra-geo/path はv0.1スコープ外だが結果型が bessel多項式(022)必須 → v0.1は path未実装方針を明文化」）。
- v0.1 完成条件（search・種別・最大食時刻・C1/最大/C4・食分食面積・50地点誤差レポート）に **`path()` は含まれない**。一方、`SolarEclipse.bessel: BesselianPolynomial`（api-draft §3.4）は v0.1 でも必須フィールドであり、ISSUE-022 が供給する。
- 本 Issue では `umbra-geo` の公開型（`GeoPoint`/`GeoLine`/`GeoPolygon`/`EclipsePath`/`PathSample`/`PathOptions`、api-draft §4）を **型として定義**し、`EclipseEngine::path()` は **v0.1 では未実装スタブ＝`Err(EclipseError::NotImplemented)`**（PATH 確定）とする。`EclipsePath::to_geojson()` も v0.1 未実装（Milestone 9）。
- 型と境界（フレーム規約・単位・日付変更線/極域の扱い方針）だけを固定し、中心線/限界線/GeoJSON の**数式本体は Milestone 9** であることを明文化する。

## 非目的
- 中心線・北限/南限・部分食域・GeoJSON の**実計算**（Milestone 9。本 Issue はスタブと型のみ）。
- 経路サンプリングの数式・日付変更線分割・極域特異点処理の実装（Milestone 9、accuracy.md / algorithms.md で別途）。
- `BesselianPolynomial`（ISSUE-022）・全球分類（ISSUE-023）の実装。本 Issue はそれらを**消費する境界**を置くのみ。
- v0.1 CLI の `path` サブコマンド本体（umbra-cli。スタブ呼出しで「未実装」を明示する整形のみ許容）。

## 公開インターフェース
api-draft §4 をそのまま型として確定（実装は Milestone 9）。
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
    pub samples: Vec<PathSample>,
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
    /// v0.1 未実装。Milestone 9 で実装。
    #[cfg(feature = "geojson")] pub fn to_geojson(&self) -> String;   // v0.1: 未実装（戻り型が String のため呼出経路に乗せない。CLI は「未実装」整形表示）
}
```
- `EclipseEngine::path(&self, eclipse: &SolarEclipse, options: PathOptions) -> Result<EclipsePath, EclipseError>`（api-draft §3.2）は v0.1 では**未実装スタブ**。
- **v0.1 スタブの戻り方（統一規則・確定 PATH）**: 「対応年代外」ではなく「機能未提供」を表すため、**`Err(EclipseError::NotImplemented)` を返す**（panic/`unimplemented!` は採用しない）。`UnsupportedTimeRange` は「対応年代外」専用語義に保ち、**未実装の意味に流用しない**。`NotImplemented` variant は ISSUE-044 で追加。CLI など実行経路に乗っても `Err(NotImplemented)` を「Milestone 9 で対応予定」と整形表示し、空 `EclipsePath` を成功として返さない。

## 数式・アルゴリズムの出典
- 本 Issue は**型・境界の確定のみ**で数式を持たない（数式本体は Milestone 9）。
- 参照（Milestone 9 で使う出典の予約・本 Issue では実装しない）:
  - 中心線・限界線・部分食域: ベッセル要素からの地上投影（Explanatory Supplement to the Astronomical Almanac、Espenak/NASA の経路生成手順）。**要確認**（一次資料の式番号は Milestone 9 で確定）。
  - GeoJSON: RFC 7946（日付変更線をまたぐ線分の分割規約・§3.1.9）。**要確認**。

## 単位 / 時刻系 / 座標系
- 角度: 公開は度（`GeodeticLatitude`/`EastLongitude`、conventions §3）。経度は**東経正** `[-180°,180°)`。
- 時刻: `PathSample.time_utc` は UTC（accuracy.md §0。TT 併記が必要なら Milestone 9 で `PathSample` を拡張、`#[non_exhaustive]` 検討）。
- 座標系: 地上点は ITRS→測地座標（WGS84、conventions §4/§5）。フレーム連鎖は ISSUE-035（GCRS→CIRS→TIRS→ITRS）に従う。
- 距離: 食帯幅は km（`Kilometers`、conventions §1）。継続時間は秒。

## アルゴリズム概要
v0.1（本 Issue のスコープ）:
1. api-draft §4 の公開型を `umbra-geo` に定義（フィールド・単位・フレーム規約を確定）。
2. `EclipseEngine::path()` を**未実装スタブ＝`Err(EclipseError::NotImplemented)`**として置く（前項「戻り方」規則、PATH）。`EclipsePath::to_geojson()` も v0.1 未実装。
3. `SolarEclipse.bessel`（`BesselianPolynomial`、ISSUE-022）は v0.1 で必須のため、**型として参照可能**にする（umbra-eclipse 側で生成。本 Issue は経路側で消費する境界の型整合のみ）。
4. ドキュメント（本 Issue・README・api-draft §4 注記）に「path は Milestone 9」と明記。

Milestone 9（本 Issue の非目的・予約）: ベッセル多項式から中心線/限界線/部分食域をサンプリングし、`EclipsePath` を構築。日付変更線分割・極域特異点処理・GeoJSON 出力。

## 受け入れテスト
v0.1（本 Issue）:
- 型整合: `EclipsePath`/`PathSample`/`PathOptions`/`GeoPoint`/`GeoLine`/`GeoPolygon` が api-draft §4 のフィールド・単位で定義され、`SolarEclipse.bessel: BesselianPolynomial` を含む `SolarEclipse` がコンパイル可能（型レベル検証）。
- スタブ契約: `path()` 呼出しが「未実装」を表す＝`Err(EclipseError::NotImplemented)` を返す（panic でなく Result、`UnsupportedTimeRange` を流用しない）。`assert!(matches!(.., Err(EclipseError::NotImplemented)))` で固定。**v0.1 の通常経路（search/local/next_visible）から path が呼ばれないこと**もテストで保証。
- CLI 整合: `umbra path`（あれば）は「Milestone 9 で対応予定」を表示し、誤った経路（空 `EclipsePath` を成功として返す等）を作らない。
- 前方互換: 列挙・設定型は `#[non_exhaustive]`/`Default` で Milestone 9 拡張時に破壊的変更を避けられる（api-draft §0）。
- 二段オラクルゲート（ISSUE-047 連動）: **本 Issue は v0.1 でスタブのため数値ゲート対象外**。Milestone 9 実装時に「M2 暫定ゲート（Mock+SOFA+NASA 経路値）」と「M10 最終ゲート（DE 差分）」を付す（ISSUE-047 の二段方針を継承）。

## 許容誤差
- v0.1（本 Issue）: 数値計算を行わないため**許容誤差なし**。
- Milestone 9（予約・本 Issue では保証しない）: 中心線位置 sub-km（≲0.5 km、幾何分、accuracy.md §1 Standard）。fit 残差は `BesselianPolynomial.fit_error`（ISSUE-022）でガード。

## 実装メモ
- 本 Issue は **milestone0-review §Minor の確定事項**の反映: 「v0.1 は path 未実装方針を明文化」。型と境界だけ定義し、本実装は Milestone 9 へ後回しを明示する。
- `SolarEclipse.bessel`（`BesselianPolynomial`）は **v0.1 で必須**（api-draft §3.4）。これは ISSUE-022 が供給し、本 Issue（umbra-geo）はそれを消費する経路側の型のみを持つ。両者の責務を混同しない。
- スタブの戻り方は ISSUE-044（`EclipseError` 集約）と整合し、**`Err(EclipseError::NotImplemented)` に統一**（確定 PATH）。`UnsupportedTimeRange` は「対応年代外」専用語義に保ち未実装に流用しない。`unimplemented!`（panic）は撤回。
- `umbra-geo` は v0.1 では実質スケルトン。ただし公開型は SemVer 境界なので、Milestone 9 で破壊しないよう `#[non_exhaustive]` とフィールド追加余地を意識する。
- レビュー重点: 「v0.1 で path を呼ばせない」保証、型の前方互換、`bessel` 必須と path 未実装の責務分離、スタブ語義の一貫性（**`Err(EclipseError::NotImplemented)` に統一**、panic/`unimplemented!` と `UnsupportedTimeRange` 流用を排除、PATH）。
