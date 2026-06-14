# ISSUE-039: 見かけ時角 μ 素材の供給（ERA/CIO ベース・GAST 相当）

- crate: umbra-ephemeris
- 依存: ISSUE-007（UT1/ΔT・`utc_to_ut1` → `Ut1Instant`・EOP）, ISSUE-035（CIO ベース歳差章動・`era_rotation`・CIP X,Y,s）, umbra-core（`Ut1Instant` / `TtInstant` / `Radians` / `JulianDate2`）
- milestone: M3（フレーム・地球回転層。035 のフレーム連鎖と同時期。021/024/028 が依存する見かけ時角の供給元）
- モード(tdd-workflow): **strict**（**確定D4**: μ・時角の供給規約は全幾何の基準。分点 GST 禁止・CIO ベース統一を型と関数で固定する公開仕様。誤ると I5 二重定義に逆戻りするため strict）

## 目的
**確定D4**（milestone0-review.md）に従い、ベッセル要素の時角 μ（= 見かけグリニッジ恒星時相当 − 影軸赤経）を構成するための**見かけ恒星時/見かけ時角の素材を ERA(iauEra00) 経由の CIO ベースで供給**する。分点 GST（iauGst06a 等）は**禁止**（I5 二重定義の解消、conventions §2/§5、architecture §5）。
- ERA（Earth Rotation Angle, UT1）と CIO ベースの章動・歳差（ISSUE-035 の X,Y,s）から、**CIO ベースの見かけ恒星時量 μ_素材**を提供。
- ISSUE-021（瞬時ベッセル要素 μ）・ISSUE-024（局地射影の局地時角 H=μ−λ）・ISSUE-028（高度方位の時角）がこれを消費。
- α_z（影軸赤経）も CIRS 基準で扱い、μ = 見かけ恒星時相当 − α_z を CIO 一貫で構成できる素材を返す（D4: α_z も CIRS 基準）。

## 非目的
- 影軸赤経 α_z 自体の算出（ISSUE-020 基本面基底）。本 Issue は μ 構成に必要な恒星時/時角素材の供給まで。
- 瞬時ベッセル要素の組立（ISSUE-021）。本 Issue は μ の「恒星時側」素材を提供する層。
- フレーム回転行列そのもの（ISSUE-035。本 Issue は ERA・CIO 量を恒星時/時角という**角度スカラ**として供給する責務）。
- UT1/EOP データ取り込み（ISSUE-007）。本 Issue は UT1 値を受けて角度を構成する。
- 分点ベース GST の提供（**禁止**。D4/I5）。

## 公開インターフェース（※署名はレビュー確定）
api-draft §3.2 周辺・architecture §5 と整合。CIO ベース見かけ時角の供給:
```rust
/// Earth Rotation Angle（CIO ベース・UT1）。SOFA iauEra00 相当。
pub fn earth_rotation_angle(time_ut1: Ut1Instant) -> Radians;

/// CIO ベースの「見かけグリニッジ恒星時相当」量。
/// ERA(UT1) に CIO locator s と CIP/章動寄与を合成した、μ 構成用の恒星時素材。
/// 分点 GST ではない（D4: CIO 統一）。time_tt は CIP/s の評価に使用。
pub fn cio_apparent_sidereal_angle(time_ut1: Ut1Instant, time_tt: TtInstant) -> Radians;

/// 見かけ時角 μ = cio_apparent_sidereal_angle − α_z（CIRS 基準の影軸赤経）。
/// 呼び出し側（ISSUE-021）が α_z を渡して構成するヘルパ。[0,2π) 正規化。
pub fn ephemeris_hour_angle(
    sidereal: Radians,        // cio_apparent_sidereal_angle の出力
    right_ascension_cirs: Radians, // α_z（CIRS 基準, ISSUE-020）
) -> Radians;
```
- 戻り値は `Radians`。恒星時/時角は `[0, 2π)` 正規化（conventions §2、赤経基準の用途）。生 f64 で角度を渡さない。
- UT1 は ISSUE-007 の `utc_to_ut1` 由来。TT は CIP/s の TT 引数（ISSUE-035）。

## 数式・アルゴリズムの出典（SOFA 関数名まで特定）
- **ERA**: IAU2000 定義。SOFA `iauEra00(UT1)`。`ERA = 2π(0.7790572732640 + 1.00273781191135448·(JD_UT1 − 2451545.0))`（IERS Conventions 2010 ch.5, eq. 5.15 相当）。ISSUE-035 の `era_rotation` と同一定義。
- **CIO ベースの恒星時/見かけ角**: IERS Conventions 2010 ch.5（CIO based transformation）。CIP 座標 X,Y と CIO locator s（SOFA `iauXys06a` / `iauS06`、ISSUE-035 由来）を用いる CIO ベース系。**分点 GST（iauGst06a / iauGst06）は採用しない**（D4: CIO 統一、I5 解消）。
  - 要確認: 「見かけグリニッジ恒星時相当」を CIO 系で μ に渡す際の正確な合成式（ERA と α_z を CIRS で揃え、s/CIP 寄与をどちらの項に含めるか）の一次式。SOFA iauC2i06a/iauEra00 の構成から導出し、実装コメントに式番号転記。NASA ベッセル μ（分点 GST 由来・度/hour, D5）との**系統差を accuracy.md に記録**（誤差に誤帰属しない、I4/I5）。
- **μ 定義**: μ = 見かけ恒星時 − α（影軸赤経）。conventions §2（時角/赤経の two_pi 正規化）。NASA/Espenak の μ（基本面経度基準, ISSUE-021 §符号規約）との対応を D5 対応表に従い記録。

## 単位 / 時刻系 / 座標系
- 単位: 角度 `Radians`、`[0, 2π)` 正規化（conventions §2）。
- 時刻系: ERA = **UT1**（ISSUE-007）、CIP/s/章動寄与 = **TT**（ISSUE-035）。conventions §6。μ は UT1（恒星時）と TT 幾何（α_z）の差（ISSUE-021 の μ は UT1 混入を明記）。
- 座標系: **CIO ベース（CIRS）で統一**（conventions §5）。分点ベースは持ち込まない（D4/I5）。α_z は CIRS 基準（D4）。

## アルゴリズム概要
1. `earth_rotation_angle`: UT1 から `iauEra00` 式で ERA を `[0,2π)` 正規化して返す。
2. `cio_apparent_sidereal_angle`: ERA に CIP X,Y・CIO locator s（ISSUE-035）の寄与を CIO ベースで合成し「見かけ恒星時相当」を構成（分点経路を経ない）。
3. `ephemeris_hour_angle`: 恒星時 − α_z（CIRS）を `[0,2π)` 正規化して μ を返す。
4. ISSUE-021 が μ を瞬時要素へ、ISSUE-024/028 が局地時角 H=μ−λ_east（東経正）へ利用。
- 数値安定性: 角度正規化は用途別（時角/赤経=two_pi、conventions §2）。連続性が要る求解側（合・接触）は呼び出し側で連続化（conventions §2, ISSUE-024）。禁止: 分点 GST 混在（D4/I5）、生 f64 時刻、出典なき係数。

## 受け入れテスト
accuracy.md テストレベル **L2/L3（時刻・地球回転）**。基準値は SOFA/IERS 公開ベクトルから（実装コピー禁止、conventions §11）。
- **ERA**: 既知 UT1（例 J2000・特定日）で `earth_rotation_angle` が SOFA `iauEra00` 参照値と一致（≲ µas 級）。
- **CIO 統一の検証（D4 必須）**: `cio_apparent_sidereal_angle` を分点 GST（iauGst06a）由来量と比較し、**差が CIO−分点の既知量（s, 赤経起点差）に一致**することを確認（=分点経路を内部で使っていないことの証明）。差は誤差でなく系統差として accuracy.md に記録（I4/I5）。
- **μ 構成**: 既知 α_z で `ephemeris_hour_angle` が手計算 μ と一致。`[0,2π)` 正規化、μ→μ+2π 不変（プロパティ, L8）。
- **下流整合**: ISSUE-021 の瞬時 μ、ISSUE-024 の H=μ−λ が本供給を使って NASA ベッセル μ（D5 単位 度/hour 対応）と慣習を揃えた上で整合（系統差記録）。
- **時刻系分離**: ERA は UT1、CIP/s は TT で評価していること（誤って TT を ERA に渡すと既知量ずれる回帰）。

## 許容誤差
- accuracy.md §2.1「歳差章動 + フレーム = 0.05″」内（IAU2006/2000A 実力 ~1mas、余裕）。ERA/CIO 量の数値再現は SOFA 参照に対し **µas 級**を目標。
- δUT1 → μ 経路: accuracy.md §2.3「δUT1 1 ms ≈ 0.46″（赤道）」。本 Issue は UT1 値を正しく角度化する責務で、δUT1 不確実性自体は ISSUE-007 の不確実性帯（局地バジェット D1/I1 へ寄与）。
- NASA μ（分点・度/hour, D5）との系統差は**誤差化せず記録**（accuracy.md §0, I4/I5）。

## 実装メモ
- **確定D4 厳守**: ERA(iauEra00) 経由の CIO ベースで μ 素材を供給。**分点 GST(iauGst06a) 禁止**。α_z も CIRS 基準（D4）。ISSUE-035 の CIO 統一方針と完全一致させる（I5 二重定義の最終解消）。
- 純Rust で SOFA を移植せず**式から実装**（ISSUE-035 同様）。ERA 係数・CIO 合成式は出典（IERS Conventions 2010 ch.5）コメント必須、magic number 禁止（conventions §11）。係数は ISSUE-041 定数 or 035 共有。
- **D5 対応表**: NASA ベッセル μ は分点基準・単位 度/hour（μ'≈15°/hour）。本供給の CIO ラジアン量との対応・系統差を ISSUE-021/022/029 と共有する対応表に記載（系統差を幾何誤差へ誤帰属しない、I4）。
- μ の UT1 混入（TT 要素に UT1 が入る）を doc に明記（ISSUE-021 §92・accuracy.md §0 と整合）。
- レビュー重点: CIO 統一（分点不使用の証明テスト）、ERA=UT1/CIP=TT の時刻系分離、two_pi 正規化、NASA μ 系統差の記録、出典コメント。
