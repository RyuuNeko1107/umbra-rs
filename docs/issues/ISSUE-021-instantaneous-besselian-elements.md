# ISSUE-021: Instantaneous Besselian elements（x,y,d,μ,l1,l2,tan f1,tan f2・単位 Re・NASA 表記対応）

- crate: umbra-eclipse
- 依存: ISSUE-019（影円錐・半角・頂点）, ISSUE-020（基本面基底・d・α_z）, ISSUE-015（見かけ地心位置）, ISSUE-039（μ・見かけ時角の供給源: ERA(iauEra00)経由 CIO ベース。分点 GST 禁止・D4）, ISSUE-007（UT1/ΔT・UTC↔TT。恒星時供給は ISSUE-039 へ移管）, ISSUE-010（WGS84 Re）, umbra-core（TtInstant, Radians）
- モード(tdd-workflow): strict（瞬時ベッセル要素の**定義と1時刻計算**。x,y,d,μ,l1,l2,tan f1,tan f2 の符号規約・Re 無次元化・NASA 表記対応がライブラリ全体の幾何の基準。誤差を隠さないための符号文書化が公開仕様級。strict）

## 目的
**1つの TT 時刻**における瞬時ベッセル要素 `InstantaneousBesselianElements`（x, y, d, μ, l1, l2, tan f1, tan f2）を、影円錐（ISSUE-019）と基本面基底（ISSUE-020）から算出する（architecture §6, api-draft §3.3）。
- x, y, l1, l2 を**地球赤道半径 Re で無次元化**（conventions §1, §4）。
- **符号規約と NASA/Espenak 表記との対応を必ず文書化**（誤差を隠さないため, 品質基準）。
- **責務**: 本 Issue は「瞬時要素の**定義**と**1時刻の値計算**」。任意 TT で連続供給する層（BesselianSource）は **ISSUE-037**。多項式近似は **ISSUE-022**。

## 非目的
- 任意 TT での供給インターフェース `BesselianSource` の実装（= ISSUE-037。本 Issue は ISSUE-037 が各時刻で呼ぶ「値計算カーネル」）。
- 多項式 fit（= ISSUE-022）。
- 影円錐・基底の構成（ISSUE-019/020 を利用するクライアント）。
- 全球種別判定・gamma（ISSUE-023。本 Issue は瞬時 x,y で gamma の素は与えるが分類はしない）。
- 局地接触・観測者投影（別 Issue）。

## 公開インターフェース
api-draft §3.3 / architecture §6 に準拠（公開型）。

```rust
#[derive(Clone, Copy, Debug)]
pub struct InstantaneousBesselianElements {
    pub time_tt: TtInstant,
    pub x: f64, pub y: f64,                 // 影軸の基本面交点座標。単位 Re（無次元）
    pub declination: Radians,               // d: 影軸の赤緯
    pub hour_angle: Radians,                // μ: エフェメリス時角。D4: ERA(iauEra00)経由 CIO ベースで構成（μ = 見かけ時角 − α_影軸。分点 GAST 禁止）。供給源 ISSUE-039
    pub l1: f64, pub l2: f64,               // 基本面での半影/本影 円錐半径。単位 Re
    pub tan_f1: f64, pub tan_f2: f64,       // 半影/本影 円錐半角の tan（無次元）
}

/// 影円錐＋基底から 1 時刻の要素を組む値計算カーネル（ISSUE-037 が各 TT で呼ぶ）。
pub(crate) fn besselian_elements_at(
    cone: &ShadowCone,                  // ISSUE-019
    basis: &FundamentalPlaneBasis,      // ISSUE-020（d, α_z 含む）
    moon_geocentric_km: Vector3,        // 影軸の基本面交点算出に使用
    earth_rotation_angle: Radians,      // D4: μ 構成。ERA(iauEra00)経由 CIO ベース（ISSUE-039 供給, UT1 依存）。分点 GST 禁止
    earth_radius_km: f64,               // Re = WGS84 a（conventions §4）
) -> InstantaneousBesselianElements;
```

- 単位は Re 無次元（x,y,l1,l2）。tan f1/tan f2 は無次元。d, μ はラジアン。
- `time_tt` を保持（TT 基準, conventions §6）。μ の恒星時は UT1 由来だが要素自体は TT 時刻でラベル。

## 数式・アルゴリズムの出典
- **第一義: Explanatory Supplement to the Astronomical Almanac (3rd ed.), Ch.11「Eclipses」の Besselian elements 節**（x, y, d, μ, l1, l2, f1, f2 の正式定義）。
- **補助: Meeus, *Astronomical Algorithms* (2nd ed.), Ch.54「Eclipses」**（同要素の実用式・係数）。
- **NASA 定義: Espenak の Besselian Elements 解説**（GSFC eclipse site / NASA TP-2006-214141, data-sources §4.1）。NASA 表記と本実装の対応を §符号規約で明示。
- **定義式（Explanatory Supplement Ch.11 / Meeus Ch.54, 符号・Re 無次元化を特定）**:
  - 影軸単位ベクトル方向の赤経 α、赤緯 d（ISSUE-020）。
  - **x, y** = 月中心（または影軸）の地心位置ベクトルを基本面 (x=東, y=北) へ射影し **Re で割って無次元化**。x = 東向き成分/Re, y = 北向き成分/Re（Meeus Ch.54 の x,y 定義, 符号＝東正/北正）。
  - **d** = 影軸赤緯（ISSUE-020）。
  - **μ**（D4 確定）= **見かけ時角 − α（影軸赤経）**。見かけ時角は **ERA(iauEra00)経由の CIO ベースで構成**（分点 GST 禁止・Standard）。α は CIRS 基準（ISSUE-020 の α_z, CIO 一貫）。供給元は **ISSUE-039**。NASA の μ（基本面の経度基準）に一致するよう構成する。**μ は UT1（ERA）依存と α（TT 幾何）の差で、慣習を明記**。旧記述の分点 GAST 直接構成は採らない。
  - **tan f1** = 半影錐半角の tan（ISSUE-019 の penumbra_half_angle）。**tan f2** = 本影錐半角の tan（umbra_half_angle）。Espenak は tan f1/tan f2 表記。符号: f2 は本影（収束）で正の小さい値、f1 は半影（発散）で正の大きい値。
  - **l1, l2** = 基本面における半影/本影 円錐半径（Re 単位）。基本面（z=0 平面）での錐の切口半径。`l1 = z_apex_dist1·tan f1 + ...`（Explanatory Supplement Ch.11 の `l1, l2` 定義式に従い、基本面での半影外縁半径 l1・本影縁半径 l2 を Re 無次元で算出）。**l2 符号の正本（B1）: l2 < 0 = 皆既（Total, 本影が地表に到達）、l2 > 0 = 金環（Annular, 本影頂点が基本面の地球側 = 反本影が地表）、経路上で符号反転 = ハイブリッド**。algorithms.md §0 を正本とし必ず実装し文書化（NASA 表記との符号対応は §NASA 対応表で系統差として明記）。

## NASA 対応表（D5・必須記載）
NASA/Espenak ベッセル要素と本実装の単位・規約対応を必ず本節と実装コメントに記載する（系統差を幾何誤差に誤帰属しないため, accuracy.md §0）。

| 量 | NASA/Espenak 表記・単位 | 本実装の単位 | 変換・規約 |
|---|---|---|---|
| x, y | 度ではなく Re 無次元（基本面座標） | Re 無次元 | 同一（東正/北正） |
| d | 度（影軸赤緯） | Radians | 度→rad |
| **μ** | **度**（基本面の経度基準、hour ではなく度で公開されることが多い ※要確認） | Radians | 度→rad。D4: ERA(iauEra00)経由 CIO ベースで構成（分点 GST 禁止） |
| **μ'（μ の時間変化率）** | **≈15°/hour**（地球自転 ≈15.041°/hour + 影軸赤経変化分。NASA は度・hour 系 ※最終単位は一次資料で要確認） | Radians/SI秒（内部微分） | hour→SI秒。μ'≈15°/hour = 15°/3600s を基準換算 |
| l1, l2 | Re 無次元（基本面の半影/本影半径） | Re 無次元 | **正本(B1): l2<0=皆既 / l2>0=金環 / 経路上符号反転=ハイブリッド**（algorithms.md §0）。NASA 表記と符号が逆の場合は系統差として記録（誤帰属しない） |
| tan f1, tan f2 | 無次元 | 無次元 | 同一 |
| ΔT 規約 | NASA カタログ採用の ΔT モデル（Espenak/Meeus）に依存 | conventions §9 / ISSUE-007 | 比較時に NASA の ΔT 規約へ揃える（系統差を記録） |
| GAST/時角 供給源 | NASA は分点系の場合あり | **本実装は CIO(ERA)基準・ISSUE-039 供給** | 分点 GST 禁止（D4）。NASA 値比較時は CIO↔分点の系統差を記録 |

- **μ の NASA 単位（度 vs hour）と μ'≈15°/hour の正確な定義は一次資料（NASA TP-2006-214141 / GSFC eclipse site）で最終確認（要確認）**。系統差は誤差として隠さず accuracy.md に記録。

## 単位 / 時刻系 / 座標系
- **時刻系: TT 基準**（conventions §6, accuracy.md §0(a)）。μ の恒星時（ERA）部分のみ UT1 由来（D4: 分点 GAST ではなく CIO ベース）で、要素は TT 時刻でラベル。
- **座標系: FundamentalPlane フレーム**（conventions §5）。x=東/y=北/z=影軸（ISSUE-020 基底）。
- **単位: x, y, l1, l2 は Re（WGS84 a, conventions §1/§4）で無次元化**。tan f1/f2 無次元。d, μ ラジアン。
- Re の値・NASA 慣習（NASA は赤道半径基準 6378.137km）との差異は conventions §4/accuracy.md に記録。

## アルゴリズム概要
1. ISSUE-019 で影円錐（半角・頂点）、ISSUE-020 で基本面基底（d, α_z）を取得。
2. 影軸（または月中心）の地心位置を基本面 (x_axis, y_axis) へ射影 → x, y（km）→ **Re で割り無次元化**。
3. d = 基底の declination。**μ = 見かけ時角 − α_z（D4: ERA(iauEra00)経由 CIO ベース, ISSUE-039 供給。分点 GST 禁止）** を `[0,2π)` 正規化（μ は赤経基準で two_pi 規約, conventions §2）。
4. tan f1 = tan(penumbra_half_angle), tan f2 = tan(umbra_half_angle)（ISSUE-019）。
5. l1, l2 = 基本面での半影/本影円錐半径（Explanatory Supplement Ch.11 式, Re 無次元）。**符号正本(B1): 皆既で l2<0、金環で l2>0（本影頂点が基本面の地球側）、経路上符号反転でハイブリッド**。
6. `InstantaneousBesselianElements` を構成して返す。
- 数値安定性: asin/acos クランプ（accuracy.md §2.2）。μ の正規化規約を固定。l2 符号の境界（皆既↔金環）で連続に符号反転すること（ISSUE-023 の種別境界の素）。

## 受け入れテスト
accuracy.md テストレベル **L4（ベッセル要素）**。**NASA 公開ベッセル値との比較＋瞬時値比較の両方**（品質基準）。
- **NASA ベッセル値比較（第二義, data-sources §4.1）**: NASA 5千年カタログの既知日食（皆既/金環/部分/ハイブリッド各）について、最大食付近の x,y,d,μ,l1,l2,tan f1,tan f2 を NASA 公開値と比較。**k 慣習を Espenak（EspenakUmbral/Penumbral, conventions §9）に揃え、ΔT も合わせる**（系統差を accuracy.md に記録, accuracy.md §3.1）。基準は fixtures（conventions §11, 実装非コピー）。
- **DE 差分（第一義, accuracy.md §3.1）**: 解析暦 vs DE440 で同一ベッセル計算を通し x,y,l1,l2 の差を層分解（暦残差 vs 幾何残差, §4）。
- **MockEphemeris 人工ケース（accuracy.md §3.1）**: 完全中心（x=y=0）、明確な金環（l2>0）、明確な皆既（l2<0）、影が地球を外す（|gamma|>1+l1）。各で要素の符号・値を解析オラクルと照合。
- **符号規約テスト（必須・品質基準, B1）**: x=東正、y=北正、皆既で l2<0、金環で l2>0 を境界をまたいで連続に確認。NASA 表記対応表（実装コメント/本 Issue）と一致。
- 単位テスト: x,y,l1,l2 が Re 無次元（km 値/Re と一致）。

## 許容誤差
- accuracy.md §2.1 幾何バジェット（影幾何誤差, §4 層分解）。最終 ±1.5s（最大食）・食分 ±0.0005 へ寄与。
- x,y の誤差は最大食時刻・gamma に直結（感度 0.5″/s, 1″≈2s）。l1,l2 の誤差は食分（0.001食分≈1.9″, accuracy.md §2.2）。
- **直接計算（本 Issue/ISSUE-037）の fit 誤差はゼロ**（暦再評価）。多項式（ISSUE-022）との残差は L7 サブテスト（accuracy.md §3.2）で比較。
- NASA 値との一致目標: 慣習を揃えた上で x,y は <0.001 Re 級（要実測, Milestone 2）。系統差（k/ΔT 慣習差）は誤差として隠さず記録（accuracy.md §0）。

## 実装メモ
- **責務分担（品質基準）**: 本 Issue = 瞬時要素の**定義＋1時刻カーネル `besselian_elements_at`**。ISSUE-037 = 任意 TT で暦を再評価して供給する `BesselianSource`（`InstantaneousEvaluator`）。ISSUE-022 = 多項式。Standard 局地の既定は ISSUE-037（直接, fit 誤差ゼロ）、経路は ISSUE-022。
- **NASA 表記対応表を実装コメントと本 Issue に必ず残す**（x,y,d,μ,l1,l2,tan f1,tan f2 の各定義・符号・単位・Re 無次元化・μ の恒星時規約）。誤差を隠さない（accuracy.md §0, conventions §10）。
- l2 の符号（**正本 B1: 皆既で負・金環で正**）は ISSUE-023 種別判定の境界条件。連続な符号反転を保証。
- **μ の恒星時部分は UT1 由来（D4: ERA(iauEra00)経由 CIO ベース, 供給源 ISSUE-039）。TT 要素に UT1 が混じる点を明記**（accuracy.md §0 絶対 UTC 律速と整合）。
- **I1（実装メモ・必須）: μ が UT1 依存ゆえ、局地接触時刻に δUT1 が混入する**。μ は UT1（ERA）で構成されるため、全球 gamma・最大食時刻 TT が純 TT で閉じるのに対し、μ を介する局地時角→局地接触時刻には δUT1 が直接効く（D1 局地 topocentric バジェットの δUT1→μ 項）。この経路を実装メモと metadata に明記し、局地接触時刻は δUT1 混入を含む脚注を付す（全球量は純 TT, 局地量は δUT1 混入）。
- レビュー重点: 全要素の符号規約・NASA 対応、Re 無次元化、**l2 符号正本（B1: 皆既で負・金環で正）**、μ の恒星時/赤経慣習、acos/asin クランプ。
