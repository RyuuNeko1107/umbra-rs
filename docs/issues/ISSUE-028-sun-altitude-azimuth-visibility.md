# ISSUE-028: Sun altitude/azimuth and visibility（方位北0東回り・幾何学的高度・大気差分離・可視性判定）

- crate: umbra-eclipse
- 依存: ISSUE-015（太陽見かけ地心位置）, ISSUE-011（観測者→ITRS）, ISSUE-039（時角・恒星時供給: ERA(iauEra00)経由 CIO ベース。分点 GST 禁止・D4）, ISSUE-007（UT1/ΔT・UTC↔TT。恒星時は ISSUE-039 へ移管）, ISSUE-024（赤緯・時角・恒星時を共有）, ISSUE-001（規約）
- モード(tdd-workflow): strict（方位北0東回り・幾何学的高度・大気差分離は conventions §7 の固定規約。可視性判定は通知用途の正否を決めるため strict）

## 目的
各接触/最大時点での **太陽の地平座標（方位角・高度）** を計算し、**可視性 Visibility** を判定する（conventions §7, architecture §7）。
- 方位角: **北=0°、東回り**（北→東→南→西。conventions §7）。
- 高度: 既定**幾何学的高度**。大気差補正は `RefractionModel { None, Standard }` で分離し、補正前後を両方返せる（conventions §7, api-draft §3.1）。
- 位置角 (position angle): 天の北=0、東回り（接触点の向き。conventions §7）。
- 可視性: `BelowHorizon / SunriseEclipse / SunsetEclipse / PartialVisible / FullyVisible / NotVisible` を網羅（api-draft §3.4）。

## 非目的
- 接触・最大・食分の求解（ISSUE-025/026/027）。本 issue はそれら時点の alt/az/PA/visible を埋める。
- 厳密な地平線付近大気差・月縁地形（accuracy.md §6 非保証）。`Standard` は標準大気差式まで。
- 太陽位置の暦計算そのもの（ISSUE-015）。本 issue はその出力（見かけ地心赤経・赤緯）を地平座標へ。

## 公開インターフェース
conventions §7、api-draft §3.1/§3.4 に整合:
```rust
#[non_exhaustive] #[derive(Clone, Copy, Debug)] pub enum RefractionModel { None, Standard }

#[derive(Clone, Copy, Debug)]
pub struct Horizontal {
    pub altitude_geometric: Degrees,    // 幾何学的高度（既定）
    pub altitude_apparent: Degrees,     // 大気差補正後（RefractionModel::Standard 時）
    pub azimuth: Degrees,               // 北0東回り（conventions §7）
}

/// 太陽の地平座標。観測者・時刻・大気差モデルから。
pub(crate) fn sun_horizontal(
    observer: Observer, time: TtInstant, refraction: RefractionModel, time_scales: &TimeScales,
) -> Result<Horizontal, EclipseError>;

#[non_exhaustive] #[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Visibility { NotVisible, BelowHorizon, SunriseEclipse, SunsetEclipse, PartialVisible, FullyVisible }

/// 接触集合＋各接触の高度から可視性を分類。
pub(crate) fn classify_visibility(contacts: &LocalContactSet, maximum: &LocalContact) -> Visibility;
```
- `Degrees`（api-draft §1.1）。生 f64 で角度を渡さない（conventions §1）。
- `position_angle` は `LocalContact`（api-draft §3.4）に格納。

## 数式・アルゴリズムの出典
- **赤道座標 → 地平座標**（球面天文の標準式）: 局地時角 `H = μ − λ`（**D4: μ・恒星時は ERA(iauEra00)経由 CIO ベースで構成・ISSUE-039 供給、分点 GST 禁止**。α=見かけ赤経は CIRS/CIO 基準, ISSUE-015）。**時角符号 H = μ − λ（東経正）と方位（北0東回り）の符号はどちらか一方を正本化し、もう一方で交差検証する**（ISSUE-024 の時角符号と整合, conventions §3/§7）。本実装は CIO 系時角を正本とし、旧記述の地方恒星時 LST（分点系）直接構成は採らない。
  - 高度: `sin a = sin φ·sin δ + cos φ·cos δ·cos H`
  - 方位（北0東回り）: `tan A = sin H / (cos H·sin φ − tan δ·cos φ)` を **`atan2` で象限確定**、北 0・東回りへ規約変換（conventions §7。要確認: 南基準→北基準の 180° オフセット符号を実装で固定）。
  - 出典: **Explanatory Supplement §7（座標系）/ Meeus AA 2nd ed. Ch.13「Transformation of Coordinates」式 13.5/13.6**。Meeus は方位を**南基準**で定義するため**北0東回りへ変換**（conventions §7）。
- **SOFA 関数名**（参照・コード移植しない）: 観測者地平座標は `iauHd2ae`（hour angle/dec → azimuth/elevation）、`iauAe2hd` 逆変換。**時角・恒星時は D4 で ERA(iauEra00)経由 CIO ベースに統一（分点 GST 禁止）。`iauGst06a`（分点 GAST）は使わない**（供給源 ISSUE-039）。要確認: SOFA 流の完全な observed place（`iauAtio13` 等）は大気差込みで、本 issue は幾何高度＋分離大気差のため alt/az 変換のみ採用。
- **大気差（RefractionModel::Standard）**: Bennett (1982) or Sæmundsson の標準式（Meeus AA Ch.16「Atmospheric Refraction」式 16.3/16.4）。出典を実装コメントに式番号で。地平線付近の厳密大気差は非保証（accuracy.md §6）。
- **位置角 (PA)**: 接触点の天の北からの角（東回り）。Explanatory Supplement §11 / Meeus Ch.54 の接触点方向式。

## 単位 / 時刻系 / 座標系
- 入力: `Observer`（測地緯度・東経・楕円体高）、太陽見かけ地心 α,δ（ISSUE-015、`Radians`）、TT。
- 恒星時/時角: UT1 由来（**D4: ERA(iauEra00)経由 CIO ベース、ISSUE-039 供給。分点 GST 禁止**）。時角 `H = μ − λ`（東経・CIO 系見かけ時角から, conventions §3/§6）。
- 出力: `Degrees`（高度・方位・位置角）。方位北0東回り、高度幾何学的（補正前後両方）、PA 天の北0東回り（conventions §7）。
- 座標系: 地平座標（観測者局所）。赤道→地平変換は CIRS/of-date 赤経赤緯を使用（ISSUE-015 の見かけ位置）。

## アルゴリズム概要
1. ISSUE-015 から太陽見かけ地心 α,δ を取得（光行時間・光行差・歳差章動込み、Standard）。
2. **（D4）** ISSUE-039 から見かけ時角（ERA(iauEra00)経由 CIO ベース。分点 GST 禁止）を取得し、`H = μ − λ`（東経正・conventions §3。μ は CIO 系の見かけ時角。符号正本は時角符号 H=μ−λ に固定し方位側で交差検証）。
3. 高度 `a_geom`、方位 `A`（atan2）→ 北0東回りへ規約変換（conventions §7）。
4. `RefractionModel::Standard` なら大気差 ΔR を加え `a_apparent`（補正前後両方返す。conventions §7）。
5. 位置角 PA を計算（接触点方向）。
6. 可視性分類: 各接触/最大の高度から `BelowHorizon`（最大時も地平下）/ `SunriseEclipse`（C1〜最大が地平下で最大以降可視 or 日の出中に進行）/ `SunsetEclipse` / `PartialVisible`（一部接触が地平下）/ `FullyVisible`（C1〜C4 すべて地平上）/ `NotVisible`（その地点で食域外）。
- 数値安定性: 高度は φ,δ 全域で安定。方位は `atan2`（極・天頂で象限保持）。時角は `[-π,π)` signed 正規化（conventions §2）。禁止: 西経正持ち込み（conventions §3）、南基準方位の混入、幾何/見かけ高度の混同。
- 日の出/日没境界: 高度 0（幾何）or −0.83°（大気差込み太陽縁）の閾値で SunriseEclipse/SunsetEclipse を判定。閾値の選択（幾何 0 か縁＋大気差か）を doc で固定（conventions §7 幾何学的高度既定）。
- 部分食地点で c2/c3=None でも可視性は C1/最大/C4 の高度で判定（ISSUE-025 と整合）。

## 受け入れテスト
accuracy.md テストレベル **L3（天体位置）依存＋L6（局地）**。基準値は実装へコピー禁止。
- 地平座標（オラクル＝独立計算 / SOFA `iauHd2ae` 参照値 / Meeus Ch.13 例題、出典明記）:
  - 既知地点・既知時刻の太陽 alt/az を Meeus 例題と照合（北0東回りへ変換後）。
  - 方位規約: 真南正中で A=180°、真東で A=90°、真西で A=270°（北0東回り確認）。
- **時角符号×方位符号 交差検証（D4・必須）**: 時角符号 `H = μ − λ`（東経正・正本）と方位（北0東回り）の符号を、一方を正本としてもう一方で交差検証する受入テスト。例: 正中（H=0）で A=180°、東の地平（日の出側, H<0）と西の地平（日没側, H>0）で A・高度の符号が整合すること、ISSUE-024 の時角符号と一致すること。CIO 系時角（ISSUE-039 供給, 分点 GST 不使用）で構成した値で固定。
- 大気差（Meeus Ch.16 例題）: 高度 0/5/45° で ΔR が標準式値と一致。補正前後の両方が返る。
- 可視性（accuracy.md L6 地点分類）: 中心線上（FullyVisible）/ 日の出中（SunriseEclipse）/ 日没中（SunsetEclipse）/ 一部地平下（PartialVisible）/ 最大時も地平下（BelowHorizon）/ 食域外（NotVisible）の 6 値すべてを網羅するフィクスチャ。
- オラクル（第二義・整合・accuracy.md §3.1）: NASA/USNO の地点別最大食高度（data-sources §4）。fixtures 転記・出典/取得日明記（ISSUE-029）。
- プロパティ（L8）: 高度 ≤ 90°、方位 `[0,360)`、λ→λ+2π で不変。

## 許容誤差
conventions §7（方位北0東回り・幾何高度・大気差分離）、accuracy.md §2 から:
- **高度・方位（幾何学的）: ≤ 0.01°（目標）**。太陽位置（ISSUE-015）と恒星時（ISSUE-007）の精度に従属。最大食高度は `GreatestEclipse.sun_altitude`/`maximum_altitude` の表示用途で 0.01° 級で十分。
- **大気差補正後**: 標準式の妥当域（高度 >5°）で式値一致。地平線付近は非保証（accuracy.md §6）。
- 可視性分類: 高度閾値（0 or −0.83°）の境界はフィクスチャで決定的に一致（pass/fail）。
- 根拠: alt/az は表示・可視性判定用途で接触時刻・食分の主精度（accuracy.md §2）には直接律速しない。ただし可視性の誤分類は通知の正否に直結するため境界を厳密に固定。

## 実装メモ
- 方位は SOFA/天文標準が南基準のことがある。**北0東回りへ必ず変換**（conventions §7）。変換オフセットの符号をテストで固定。
- 高度は幾何学的が既定。`RefractionModel::Standard` で補正前後両方返す（通知用途・conventions §7）。
- 恒星時・赤緯・時角は ISSUE-024（基本面射影）と共有し重複評価を避ける。
- 位置角（PA）は天の北0東回り（conventions §7）。接触点方向の符号を実装コメントに固定。
- Visibility 6 値の判定木を明文化（最大時の高度 → 接触ごとの高度 → 食域内外）。日の出/日没閾値を 1 箇所で定義。
- レビュー重点: 北0東回り変換、幾何/見かけ高度の分離、可視性 6 値網羅、日の出日没境界閾値、SOFA コード非移植（参照のみ）。
