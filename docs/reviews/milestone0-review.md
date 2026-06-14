# Milestone 0 独立レビュー記録（2026-06-14）

37 Issue + docs5本に対し、起票者の自己評価を渡さず独立2観点でレビューした結果の要約。
このレビューに基づく確定判断は本ファイル末尾「確定事項」を参照（docs本体・各Issueへ反映済み）。

## 観点1: 完全性・ロードマップ

### Critical
- **043 EclipseEngine組立** が無い（search/local/next_visible の結線本体）。
- **038 MockEphemeris実装** が無い（多数Issueがテスト前提・DE無しの幾何検証足場）。
- **039 GAST/恒星時供給** が宙吊り（μ=GAST−α の供給元不在）。
- **第一義オラクル=JPL DE(036)がM10なのにM2暦確定がDE依存** → 047 暫定オラクル戦略＋036のM2前倒し（最小SPK reader）で二段ゲート化。

### Important
- 040 章動係数(IAU2000A)データ管理、041 物理定数集約、042 TimeData/TimeScales、044 EclipseError集約。
- 013/014 の依存行に 036（テストオラクル）が未記載。
- 全Issueに milestone 行が無くplanとの突合が機械化できない → ISSUE-INDEX を新設。

### Minor
- 045 umbra-geo/path はv0.1スコープ外だが結果型が bessel多項式(022)必須 → v0.1は path未実装方針を明文化。
- 046 xtask骨子・CI実体・cargo-deny/audit運用。

## 観点2: 天文・数値精度

### Critical
- **C1 局地 topocentric視差がバジェットに無い**: 月距離誤差→地平視差(≈3422″)→局地接触時刻 の経路が§2.1に未配分。

### Important
- **I1 μがUT1依存** → 「TT基準の局地幾何精度」にδUT1が混入（§0(a)と形式矛盾）。局地バジェットへ δUT1→μ 項を追加。
- **I2 最大食=距離最小化の√eps律速＋皆既帯の平底** → dm/dt=0 求根へ。m最小化(26)とm²(25)の混在を m² に統一。
- **I3 光行差の適用順序** がSOFA(iauAtciq)と逆（章動後にaberration）。順序確定が必要。
- **I4 NASAベッセルμの単位(度/時)** 対応が未定義 → 系統差を幾何誤差に誤帰属する恐れ。
- **I5 GAST(分点)とERA(CIO)の二重定義** → ERA経由CIOへ統一。
- **I6 偽陰性ゼロのマージン導出式** が未定量。

### Minor
- M1 相対論偏向の省略上限を受入テスト化、M2 章動係数データ管理の独立化、M3 速度のTDB/TT、M4 α_zのCIO一貫、M5 極域経路の特異点fixture。

### 精度バジェットの結論
地心側（全球gamma・最大食時刻TT）は ±1〜2s で閉じる見込みは妥当（RSS≈1.0s、目標1.5sに余裕）。ただし**局地接触・局地最大食は (a)月距離→視差、(b)δUT1→μ の2経路が未配分で、現状docのままでは過小評価**。別建てバジェットでRSS再合成すれば閉じる可能性は高い。一次資料確認を要するのは I3/I4/I5/I2。実測（M2 DE差分）前に保証値化しない方針は正しい。

---

## 確定事項（本レビューを受けた決定・docs/Issueへ反映）

### 精度決定
- D1 局地 topocentric バジェットを別建て（月距離→視差、δUT1→μ を項に追加しRSS再合成）。全球gamma・最大食時刻TTは純TT、局地接触時刻はδUT1混入を含むと脚注。
- D2 最大食時刻は dm/dt=0 の求根（Brent）を正式手法、距離最小化は粗ブラケット用。最小化対象は m²(=u²+v²) に統一（中心線尖点回避）。皆既帯平底fixtureを追加。
- D3 光行差適用順序を SOFA iauAtciq に固定（GCRS内で deflection→aberration→その後 bias+歳差+章動でCIRS）。
- D4 μ・時角は ERA(iauEra00)経由のCIOベースで統一。分点GST禁止（Standard）。039 が見かけ時角を供給。α_z もCIRS基準。
- D5 NASAベッセルμ単位(度・hour, μ'≈15°/hour)対応表を 021/022/029 に必須記載。
- D6 偽陰性マージン = 月地平視差 + 月最大黄緯速度×合↔最大食ずれ + 概算暦誤差上限。029で実余裕を統計出力。

### API決定（api-draft §6 を確定・改訂可）
- A1 EclipseEngine はジェネリック(E,D,O)維持。`StandardEngine` 型エイリアス＋`standard_engine()` で利便確保。dynは使わない。
- A2 TimeData::bundled()（埋込）+ from_path()。valid_to超過は Missing*Data を返しmetadataに記録。
- A3 接触は LocalContactSet（c1..c4=Option、maximum非Option）。
- A4 event_key=最大食UTC日付+lunation番号、location_key=指定 or 緯度経度丸めハッシュ。
- A5 NonCentral種別はenum公開(non_exhaustive)。v0.1 CLIは主にPartial/中心。
- A6 thiserror採用（deny allow-list追加）。
- A7 serde: 列挙型 `#[serde(tag="type")]`、単位はフィールド名で明示、feature gate。
