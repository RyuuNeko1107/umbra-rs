# ISSUE-030: Error statistics report（誤差統計・層分解・ToleranceProfile）

- crate: umbra-fixtures
- 依存: ISSUE-029（ゴールデン20）, ISSUE-024〜028（局地出力）, ISSUE-036（JPL DE 差分・第一義オラクル）, ISSUE-001
- モード(tdd-workflow): strict（検証基盤。pass/fail だけでなく誤差統計と層分解（accuracy.md §4）を出す設計。誤差を隠さない＝レビューゲートの正本のため strict）

## 実装状況（2026-06-20）
- **S30a プリミティブ**: `ErrorStats`（max/mean/p95・R-7）/ `ToleranceProfile`（standard/reference）— 実装済み。
- **S30b/c/d ゴールデン比較**: `compare_global`/`compare_local`/`aggregate_*`/`GoldenComputer`/`report_against_golden`/`render_text`/`render_json` — 実装済み。
- **S30f CLI 統合**: `xtask validate`（実エンジンで同梱ゴールデン照合・text/json）— 実装済み。
- **DE 差分・誤差層分解の純コア**（strict 全工程）: `LayeredError`/`DifferentialReport`/`report_differential`/`render_differential_{text,json}`（umbra-fixtures）— 実装済み（mutation 31 中 29 caught・2 unviable・生存0）。**層分解は 2 バケット（暦層 ephemeris / 幾何・数値層 geometry）へ集約**＝6 物理層スケッチ（§公開IF）からの確定逸脱を accuracy.md §4.1 に記録。
- **DE 実エンジン結線ハーネス**（strict 全工程）: `xtask differential` サブコマンド＝`JplGoldenComputer`（DE440s Reference 版 `GoldenComputer`・feature `jpl`）＋`differential_report`（注入2エンジンで層分解レンダ）。FAST 5＋SLOW 1（実 DE440s×解析暦で 2017-08-21-total を end-to-end 層分解・compared==1・105s）。mutation（differential_report）2 中 2 caught。
- **層別誤差統計 `ErrorReport`**（strict 全工程）: `report_stratified`（umbra-fixtures）＝`by_metric`（全 metric 全体統計・固定 7 件）＋`by_era`/`by_kind`/`by_location_class`（**層別 metric = 局地最大食接触時刻誤差 秒**）＋`pass_fail`（ToleranceProfile 判定）＋`render_stratified_{text,json}`。年代は 50 年バケット `[start,start+50)`（`div_euclid`・半開区間）、食種/地点条件は出現分のみ・Debug 文字列昇順（BTreeMap で順序保証）。JSON は enum キーを Debug 文字列へ写像（`serialize_with`・SolarEclipseKind に serde 非依存）。19 テスト / mutation 15 中 13 caught・2 unviable・生存0。
- **全食スイープ 自己カタログ集計＋完備性突合 純コア**（strict 全工程）: `summarize_sweep`（umbra-fixtures）＝`total`＋`by_kind`（raw 種別件数・非中心は畳まず・Debug 昇順）＋`gamma_abs`/`magnitude`（`RangeStats` min/max/mean）＋`completeness`（NASA 4 区分＝非中心は中心へ畳む検出 vs 期待・`all_match`）＋`render_sweep_{text,json}`。19 テスト / mutation 28 中 26 caught・2 unviable・生存0（`RangeStats` の min/max は `f64::min`/`max` 化で `<`/`>` の等価変異を構造排除）。
- **`xtask sweep` 実走ランナー**（strict 全工程）: `parse_range`（`--from/--to` 年・既定 1900-2100）＋`parse_expected_counts`（`--expected-*`・既定0）＋`sweep_report`（summarize_sweep＋render）＋`run_sweep`（解析暦 `standard_engine` で範囲 search→集計→印字）。完備性 expected は CLI フラグ供給＝オラクル件数は利用者が出典付き指定（ハードコードしない・参考 NASA 5MCSE 2001-2100=全224件）・既定カタログのみ。FAST 15＋SLOW 1（実エンジンで 2024-04-08 皆既を狭窓 search→集計・total≥1/Total 検出・98s）。mutation（parse/sweep_report）13 中 10 caught・3 unviable・生存0。

- **D6 偽陰性ゼロ・マージン実余裕統計**（strict 全工程）: `umbra_eclipse::scan_filter_margins`／`aggregate_filter_margins`／`FilterMarginStats`／`MarginSample`／`ECLIPSE_FILTER_SAFETY_MARGIN_RAD`。candidate→合→フィルタ（engine.search 前段と同一経路・Besselian/分類省く軽量）を範囲全体で回し、採用候補のマージン実消費 `max(0, separation − bare_limit)` と実余裕 `margin − 実消費` を統計化。`eclipse_filter` の `SAFETY_MARGIN_RAD` を `pub(crate)` 化＋`EclipsePossibility.bare_limit` 追加。16 テスト（FAST 14＋SLOW 2）/ mutation 8 中 6 caught・2 unviable・生存0。**実測（1972-2100・1583朔）: 採用 414、実余裕最小 0.000032 rad（マージン 0.0087 の 0.4%）＝マージンほぼ枯渇寸前で保持、縮小不可**（accuracy.md §3.4）。

## 完了状況（2026-06-21・完全クローズ）
v0.1 検証層 ISSUE-030 は **完了**（プリミティブ／ゴールデン比較／CLI validate／DE 差分・層分解 純コア＋xtask differential／層別 ErrorReport／全食スイープ 純コア＋xtask sweep／**D6 マージン実余裕統計**）。当初 §公開IF スケッチからの確定逸脱（6物理層→暦層/幾何層の 2 バケット・accuracy.md §4.1）は記録済み。残課題なし。

## 目的
ゴールデン20（ISSUE-029）および 1900〜2100 全日食比較に対し、**pass/fail だけでなく誤差統計を生成**するレポータを実装する（accuracy.md §3.4/§4）。
- 統計項目: **最大 / 平均絶対 / 95%ile**、さらに**年代別 / 食種別 / 地点条件別**に層別（accuracy.md §3.4）。
- **誤差の層分解**（accuracy.md §4, conventions §11）: 時刻変換 / 太陽位置 / 月位置 / 影幾何 / 多項式近似 / 接触 solver の各層へ帰属。JPL DE 差分（ISSUE-036・第一義）で暦誤差を切り出し、残りを幾何/数値層へ。
- 比較許容は `ToleranceProfile`（plan §18、モデル別管理）で定義。

## 非目的
- フィクスチャ整備（ISSUE-029）。本 issue は比較・統計・レポート。
- 暦・幾何・solver の実装そのもの（ISSUE-013〜028）。本 issue はそれらの出力を評価。
- CI 常時実行の強制（JPL DE 差分は巨大データ・nightly/手動・accuracy.md §3.1）。

## 公開インターフェース
accuracy.md §3.4/§4、architecture §1 に整合:
```rust
#[derive(Clone, Debug)]
pub struct ToleranceProfile {            // plan §18・モデル別
    pub contact_seconds: f64,            // 接触 ±2s（accuracy.md §2）
    pub maximum_seconds: f64,            // 最大食 ±1〜2s
    pub magnitude: f64,                  // 食分 ±0.0005
    pub obscuration: f64,
    pub altitude_degrees: f64,
    pub note_utc_is_delta_t_limited: bool, // §0(b): UTC 絶対は ΔT 律速
}
impl ToleranceProfile { pub fn standard() -> Self; pub fn reference() -> Self; }

/// 1 項目の誤差統計。
#[derive(Clone, Debug)]
pub struct ErrorStats { pub n: usize, pub max_abs: f64, pub mean_abs: f64, pub p95: f64, pub units: &'static str }

/// 誤差の層分解（accuracy.md §4）。
#[derive(Clone, Debug)]
pub struct LayeredError {
    pub time_conversion: ErrorStats,
    pub sun_position: ErrorStats,
    pub moon_position: ErrorStats,
    pub shadow_geometry: ErrorStats,
    pub polynomial_fit: ErrorStats,
    pub contact_solver: ErrorStats,
}

/// 層別レポート（年代別/食種別/地点条件別）。
#[derive(Clone, Debug)]
pub struct ErrorReport {
    pub by_metric: Vec<(String, ErrorStats)>,        // 接触/最大/食分/食面積/高度/方位
    pub by_era: Vec<(String, ErrorStats)>,           // 例 1900-1950, ...
    pub by_kind: Vec<(SolarEclipseKind, ErrorStats)>,
    pub by_location_class: Vec<(LocationClass, ErrorStats)>, // 中心線上/北南限/部分食域/...
    pub layered: LayeredError,                        // accuracy.md §4
    pub pass_fail: Vec<(String, bool)>,              // ToleranceProfile 判定
}

pub enum LocationClass { OnCenterLine, NearCenterLine, NorthLimit, SouthLimit, PartialZone, GrazingLimit, OutOfVisibility, Sunrise, Sunset, HighElevation }

pub fn report_against_golden(golden: &[GoldenEclipse], profile: &ToleranceProfile) -> ErrorReport;
pub fn report_differential_jpl(/* engine refs */ profile: &ToleranceProfile) -> ErrorReport; // ISSUE-036
```

## 数式・アルゴリズムの出典
- 統計: 最大絶対 `max|e|`、平均絶対 `mean|e|`、95 パーセンタイル（線形補間 percentile）。標準的記述統計（出典＝統計の定義、特定文献不要・要確認: percentile 補間規約を 1 つに固定し doc 明記）。
- **層分解の方法論**（accuracy.md §4）: 同一のベッセル/接触パイプラインに解析暦と JPL DE を通し差分（accuracy.md §3.1 第一義）。
  - 暦誤差 = 解析暦 vs DE の太陽/月位置差（→ sun_position/moon_position 層）。
  - 幾何/数値誤差 = DE 入力でのパイプライン出力 vs オラクル（→ shadow_geometry/polynomial_fit/contact_solver 層）。
  - 時刻変換誤差 = ΔT/UT1/閏秒モデル差（→ time_conversion 層、accuracy.md §2.3）。
- **慣習差の分離**（conventions §9, accuracy.md §3.1）: k（Espenak 2値）・ΔT 慣習を揃えた上で比較し、揃わない分は系統差として別掲（絶対基準にしない）。

## 単位 / 時刻系 / 座標系
- 誤差単位: 接触/最大=秒、食分/食面積=無次元、高度/方位=度（`ErrorStats.units` に明示）。
- 時刻: **TT 基準の幾何誤差**（計算律速・accuracy.md §0(a)）と **UTC 絶対誤差**（ΔT/UT1 律速・§0(b)/§2.3）を**分けて報告**。UTC 列には不確実性帯を併記（`delta_t_uncertainty_seconds`）。
- 座標: 中心線位置誤差は km（sub-km 目標・accuracy.md §1）。

## アルゴリズム概要
1. ゴールデン20（ISSUE-029）/ 全日食について実装出力とオラクルを比較、metric 別誤差列を収集。
2. 層分解: JPL DE 差分（ISSUE-036）で暦層を切り出し、残差を幾何/数値層へ帰属（accuracy.md §4）。
3. 各 metric・各層を 年代別/食種別/地点条件別に集計（max/mean_abs/p95）。
4. `ToleranceProfile` で pass/fail 判定（接触±2s/食分±0.0005、UTC は ΔT 律速注記・accuracy.md §2/§0）。
5. レポート出力（人間可読＋機械可読 JSON。CLI validate から呼べる）。
- 注意: TT 基準と UTC 絶対を混ぜない（accuracy.md §0）。慣習差（k/ΔT）を系統差として分離。**誤差を層へ分解し、日食側で暦誤差を打ち消して隠さない**（conventions §11, accuracy.md §4, レビューゲート §25/§13）。

## 受け入れテスト
accuracy.md §3.4/§4、テストレベル **L7（差分）＋回帰**:
- 統計の正しさ（L1 メタテスト）: 既知誤差列（手計算）に対し max/mean_abs/p95 が一致。percentile 補間規約の境界（n 小・同値）で安定。
- 層分解: MockEphemeris（accuracy.md §3.1）で暦誤差ゼロの人工配置 → 暦層 ≈0、幾何/solver 層のみに誤差が出ることを確認（層帰属の正当性）。
- ToleranceProfile 判定: ゴールデン20 で各 metric の pass/fail とともに**統計が必ず出る**（pass でも数値を出す・accuracy.md §3.4）。
- 層別: 年代別/食種別/地点条件別の集計が漏れなく出る（被覆メタテスト）。
- UTC/TT 分離: UTC 列に不確実性帯が付き、TT 列は計算律速として別掲（accuracy.md §0）。
- 慣習差: k/ΔT 慣習を揃えた場合と揃えない場合で系統差が分離表示される。

## 許容誤差
本 issue は**統計レポータ**であり自身の計算許容は持たない。判定に使う許容は `ToleranceProfile`＝accuracy.md §2:
- 接触 **±2 s**（TT 基準・幾何分）、最大食 **±1〜2 s**、食分 **±0.0005**、中心線 sub-km。
- **UTC 絶対は ΔT/UT1 予測律速**（accuracy.md §0(b)/§2.3）。UTC 判定には不確実性帯を考慮し、将来日食では帯を超える誤差を「予測律速」として分類（計算誤差と混同しない）。
- 統計自体の数値精度: f64 丸めのみ（≤1e-12 相対）。
- 根拠: pass/fail を通すために許容を拡大しない（conventions §11）。許容超過は層分解で原因層を特定（accuracy.md §4）。

## 実装メモ
- レポートは人間可読（表）＋ JSON（CI/履歴比較）。CLI `validate`（umbra-cli）から呼べる形に。
- JPL DE 差分（ISSUE-036）は feature `jpl`・nightly/手動（巨大データ・accuracy.md §3.1）。DE 無し時はゴールデン20（第二義）のみで統計を出し、層分解は縮退（暦層を「未測定」と明示）。
- UTC と TT を必ず分離（accuracy.md §0）。将来日食の UTC 誤差は不確実性帯付きで「予測律速」ラベル。
- 慣習差（k/ΔT、conventions §9）は系統差として別掲・絶対基準にしない（accuracy.md §3.1）。
- レビュー重点: 層分解の帰属正当性、TT/UTC 分離、percentile 規約固定、pass でも統計を出す、誤差を隠さない（レビューゲート §13/§25）。
