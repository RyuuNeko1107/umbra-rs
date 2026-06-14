# ISSUE-029: Golden fixtures（ゴールデン20・固定回帰・umbra-fixtures crate）

- crate: umbra-fixtures
- 依存: ISSUE-024〜028（局地条件の出力型一式）, ISSUE-021〜023（ベッセル要素・多項式・全球分類）, ISSUE-016〜019（新月・合・候補・影円錐＝search 経路）, ISSUE-001（規約）
- モード(tdd-workflow): strict（検証基盤。基準値の出典・k/ΔT 慣習の取り違えが全テストの信頼性を毀損する。オラクルのハードコード禁止＝数値事実のみ転記の規律が必須なため strict）

## 目的
固定回帰用の **ゴールデン20** データセットを `umbra-fixtures` crate に整備する（accuracy.md §3.4, architecture §1）。
- 構成: **皆既5 / 金環5 / 部分3 / ハイブリッド2 / 境界・日の出日没・極域5**（計20日食）、各日食に**地点5〜10**（accuracy.md §3.4）。
- 各日食・各地点の基準値（種別・最大食時刻・gamma・食分・C1〜C4・高度方位）を **NASA 5千年カタログ / USNO（data-sources §4）から数値事実のみ転記**（出典・取得日・系列バージョン明記、ライセンス配慮・data-sources §0/§4）。
- 固定回帰: 実装出力とゴールデン値を `ToleranceProfile`（ISSUE-030）で比較。`umbra-fixtures` は通常依存に含めない（検証専用・architecture §1）。

## 非目的
- 誤差統計レポートの生成（ISSUE-030）。本 issue はデータ整備とローダ。
- JPL DE 差分（第一義オラクル・accuracy.md §3.1）。本 issue は第二義（NASA/USNO 整合）の固定フィクスチャ。差分テストは別 issue（L7）。
- 1900〜2100 全日食一括比較（accuracy.md §3.4 後段）。本 issue は厳選 20 の固定回帰。
- オラクル値の再計算・補間（**数値事実のみ転記**。conventions §11 ハードコード禁止の精神＝実装側へコピー禁止だが、フィクスチャ側は外部オラクルの**転記**として出典付きで保持）。

## 公開インターフェース
architecture §1（umbra-fixtures）、api-draft 各結果型に整合:
```rust
/// 出典付きのゴールデン日食フィクスチャ。
#[derive(Clone, Debug)]
pub struct GoldenEclipse {
    pub event_key: String,
    pub kind_expected: SolarEclipseKind,
    pub greatest_time_tt: TtInstant, pub greatest_time_utc: UtcInstant,
    pub gamma: f64, pub magnitude: f64,
    pub locations: Vec<GoldenLocation>,
    pub source: OracleSource,           // 出典・取得日・k/ΔT 慣習
}
#[derive(Clone, Debug)]
pub struct GoldenLocation {
    pub name: String, pub observer: Observer,
    pub c1: Option<GoldenContact>, pub c2: Option<GoldenContact>,
    pub maximum: GoldenContact,
    pub c3: Option<GoldenContact>, pub c4: Option<GoldenContact>,
    pub magnitude: f64, pub obscuration: f64,
    pub max_altitude: f64, pub max_azimuth: f64,
    pub visibility_expected: Visibility,
}
#[derive(Clone, Copy, Debug)] pub struct GoldenContact { pub time_utc: UtcInstant /* +TT if oracle gives */, pub altitude: f64 }
#[derive(Clone, Debug)]
pub struct OracleSource {
    pub name: String,          // 例 "NASA Five Millennium Catalog of Solar Eclipses"
    pub url: String, pub retrieved: String,   // 取得日
    pub delta_t_convention: String,           // 慣習（accuracy.md §3.1）
    pub k_convention: String,                  // Espenak 2値 等（conventions §9）
    pub license_note: String,                  // data-sources §4 ライセンス配慮
}

pub fn golden_twenty() -> Vec<GoldenEclipse>;   // 固定回帰の正本
```
- `serde` feature でフィクスチャをデータファイル（TOML/JSON）として保持し、ローダで読む（手書き Rust 配列に直書きしない・architecture §11）。

## 数式・アルゴリズムの出典
- 本 issue は計算ではなく**検証データ整備**。基準値の出典:
  - **NASA 5千年日食カタログ（Espenak & Meeus, NASA/TP-2006-214141, eclipse.gsfc.nasa.gov）**（data-sources §4.1）。種別・最大食時刻・gamma・食分・（地点別の場合）接触時刻。
  - **USNO / 各国機関の地点別予報値**（data-sources §4.2）。局地接触・高度。
- **慣習の明記**（accuracy.md §3.1, conventions §9）: NASA は固有 ΔT・**Espenak 2値 k**（`EspenakUmbral`/`EspenakPenumbral`、conventions §9）。フィクスチャ比較時は本 crate も同慣習へ切替。慣習を `OracleSource` に記録し、整合チェック（絶対基準にしない・accuracy.md §3.1）。
- ライセンス（data-sources §0/§4）: **数値事実のみ転記**、サイト掲載物の体裁・著作性ある表現は転記しない。出典・取得日・系列バージョンを併記。OSS 公開前ライセンス確認チェックリスト（data-sources §6）に従う。

## 単位 / 時刻系 / 座標系
- 時刻: 接触・最大は UTC+TT（accuracy.md §0、オラクルが TT を持たない場合は UTC のみ転記し、その旨注記）。
- 角度: 度（公開・転記しやすさ。conventions §1 では内部ラジアンだがフィクスチャは度で保持し境界変換）。方位北0東回り（conventions §7）。
- 座標: 観測者は測地緯度・東経・楕円体高（conventions §3/§4）。西経表記のオラクルは東経正へ変換し記録（conventions §3）。
- gamma/食分は無次元。

## アルゴリズム概要（データ整備手順）
1. 20 日食を選定: 皆既5/金環5/部分3/ハイブリッド2/境界・日の出日没・極域5（accuracy.md §3.4）。多様な年代（1900〜2100、accuracy.md §6）・gamma・食種を含める。
2. 各日食に地点5〜10: 中心線上/付近/北南限/部分食域/限界/可視域外/日の出中/日没中/標高差（accuracy.md L6 の地点分類を被覆）。
3. NASA/USNO から数値事実を転記、`OracleSource`（出典・取得日・系列・慣習・ライセンス注記）を付与。
4. データファイル（TOML/JSON）に保存、ローダ `golden_twenty()` を実装。
5. checksum/version で固定（architecture §11、変更検知）。
- 注意: 西経入力の東経正変換、k/ΔT 慣習の記録、TT 有無の明示。偽の精度を作らない（オラクルが秒単位なら秒で保持）。

## 受け入れテスト
accuracy.md §3.4（ゴールデン20）、テストレベル **L6（局地）＋回帰**:
- ローダ: `golden_twenty()` が 20 件・各 5〜10 地点・各 `OracleSource` 充足を返す（構造検証）。
- 被覆検証: 食種（皆既/金環/部分/ハイブリッド/境界）と地点分類（中心線上〜標高差〜可視域外〜日の出日没）が漏れなく含まれる（メタテスト）。
- 出典完全性: 各 `OracleSource` に name/url/retrieved/慣習/ライセンス注記が非空（ハードコード防止・data-sources §4）。
- 慣習整合: k_convention が Espenak 2値 or 明示。比較は ISSUE-030 が慣習を揃えて実施。
- 値妥当性（サニティ）: gamma ∈ [-1.6,1.6]、食分 >0、皆既で食分 ≥1、接触順序 c1<max<c4。**実装出力との一致は ISSUE-030 の統計で評価**（本 issue はフィクスチャの内部整合のみ）。
- 固定性: checksum 一致（変更検知）。

## 許容誤差
本 issue は基準データの**整備**であり計算許容は持たない。比較許容は ISSUE-030 の `ToleranceProfile`（accuracy.md §2: 接触±2s/食分±0.0005、UTC は ΔT 律速 §0/§2.3）が定義。
- フィクスチャ自身の品質基準: オラクルの**有効桁を超える精度を捏造しない**（秒単位オラクルは秒で保持）。出典の系列・取得日でトレーサブル。
- 慣習差（k/ΔT、conventions §9, accuracy.md §3.1）は系統差として ISSUE-030 で分離報告（絶対基準にしない）。

## 実装メモ
- `umbra-fixtures` は通常依存に含めない（architecture §1）。dev-dependency / 検証 crate。
- データは手書き Rust 配列直書き禁止（architecture §11）。TOML/JSON + ローダ + checksum。
- 数値事実のみ転記・出典/取得日/慣習/ライセンス注記必須（data-sources §0/§4/§6）。GPL 等の取り込み不可データに触れない（data-sources §0）。
- オラクルが UTC のみ・TT なしの場合、TT は導出せず空にし注記（accuracy.md §0 の UTC/TT 分離を捏造しない）。
- 西経表記オラクルは東経正へ変換し元値も注記（conventions §3）。
- レビュー重点: 食種/地点分類の被覆、出典完全性、慣習記録、有効桁の誠実さ、ライセンス配慮、checksum 固定。
