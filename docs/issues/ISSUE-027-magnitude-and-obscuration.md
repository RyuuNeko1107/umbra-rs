# ISSUE-027: Magnitude and obscuration（食分＝食い込み量・食面積＝円交差面積）

- crate: umbra-eclipse
- 依存: ISSUE-024（基本面射影で m, L1, L2 を供給）, ISSUE-021（瞬時要素 l1,l2,tan f）, ISSUE-001（規約）, ISSUE-026（最大食時点で本式を評価）
- モード(tdd-workflow): strict（食分 ±0.0005・食面積の境界条件と acos クランプが精度と健全性を律速。境界処理の取りこぼしが NaN/誤値になるため strict）

## 目的
観測地点の **食分 (magnitude＝太陽直径に対する食い込み量)** と **食面積 (obscuration＝太陽面積に対する月が覆う割合＝2 円の重なり面積比)** を計算する（conventions §8/§9, architecture §7）。
- 食分: `magnitude = (L1 − m)/(L1 + L2)`（基本面上の外接縁から内接縁への食い込み割合）。皆既で 1 超を許容、金環で ≈1。
- 食面積: 太陽円と月円の**重なり面積 / 太陽円面積**。2 円交差面積の幾何式（acos 引数を `[-1,1]` にクランプ、境界＝離隔/内包/外接/内接/同半径を明示処理。accuracy.md §2.2）。

## 非目的
- 最大食時刻の求解（ISSUE-026）。本 issue は与えられた瞬時要素（m, L1, L2 等）から食分・食面積を返す純関数。
- 接触判定（ISSUE-025）。
- 太陽視半径・月視半径の生成（瞬時要素 L1/L2 = 影半径、太陽/月見かけ半径は ISSUE-015/018 由来）。

## 公開インターフェース
api-draft §3.4 `EclipseMagnitude` / `Obscuration` に整合:
```rust
#[derive(Clone, Copy, Debug, PartialEq, PartialOrd)] pub struct EclipseMagnitude(pub f64); // 皆既で1超
#[derive(Clone, Copy, Debug, PartialEq, PartialOrd)] pub struct Obscuration(pub f64);       // 0..1

/// 基本面の中心間距離 m と影半径 L1/L2 から食分。
/// 太陽見かけ半径 r_sun・月見かけ半径 r_moon（同単位）から食面積。
pub(crate) fn eclipse_magnitude(m: f64, l1: f64, l2: f64) -> EclipseMagnitude;

pub(crate) fn eclipse_obscuration(
    separation: f64,  // 太陽-月 見かけ中心間距離（同単位）
    r_sun: f64, r_moon: f64,
) -> Obscuration;
```
- 単位は呼出側で統一（基本面 Re or 角度。食面積は太陽/月見かけ半径と中心離隔を同単位で）。生 f64 の単位はコメントで固定（conventions §1 の精神／本関数は内部 pub(crate)）。

## 数式・アルゴリズムの出典
- **食分**: Explanatory Supplement §11 / Meeus Ch.54 — `magnitude = (L1 − m)/(L1 + L2)`（部分食 0..1、皆既で >1、金環で m≈L2 付近）。出典式を実装コメントに章・式番号転記。
- **食面積（2 円の重なり面積）**: 円-円交差面積（lens area）の標準幾何式。半径 R, r、中心間距離 d:
  - `A = R²·acos((d²+R²−r²)/(2dR)) + r²·acos((d²+r²−R²)/(2dr)) − ½·√((−d+r+R)(d+r−R)(d−r+R)(d+r+R))`
  - obscuration = `A_overlap / (π·R_sun²)`。出典: 標準幾何（Weisstein, MathWorld "Circle-Circle Intersection"。要確認: 章・式番号ではなく標準公式のため出典は公式名で明記）。日食文脈の食面積定義は Explanatory Supplement §11 / Meeus Ch.54。
- **acos 引数クランプ**: 丸め誤差で引数が ±1 をわずかに超える → `[-1,1]` にクランプ（accuracy.md §2.2）。
- **境界条件**（accuracy.md §2.2 明示処理）:
  - 離隔（`d ≥ R+r`）: overlap=0 → obscuration=0。
  - 内包（`d ≤ |R−r|`）: overlap = π·min(R,r)² → 小円が大円に完全内包。月が太陽より大（皆既近傍）は obscuration=1、太陽が大（金環）は `r_moon²/R_sun²`。
  - 外接（`d = R+r`）/ 内接（`d = |R−r|`）/ 同半径（`R=r`）の縮退で 0 除算（2dR, 2dr, d=0）回避。

## 単位 / 時刻系 / 座標系
- 入力: 食分は基本面 Re（m, L1, L2）。食面積は太陽/月見かけ半径と中心離隔（角度 `Radians` or Re、同単位で固定）。
- 時刻系: なし（瞬時量。呼出側が時刻を管理）。
- 座標系: FundamentalPlane（食分）/ 視半径平面（食面積）。
- 月半径係数 k（conventions §9）の選択は L1/L2（瞬時要素）側で反映済み。本 issue は供給値を使う（皆既/金環境界の系統差は metadata、conventions §9）。

## アルゴリズム概要
1. 食分: `(L1 − m)/(L1 + L2)`。`m ≥ L1`（離隔）なら 0 にクランプ（食なし）。皆既/金環で 1 跨ぎを許容（EclipseMagnitude は 1 超可）。
2. 食面積: d, R(=r_sun), r(=r_moon) で場合分け（離隔/内包/部分重なり）。部分重なりのみ lens 公式、acos 引数を `[-1,1]` クランプ。
3. obscuration = overlap / (π R_sun²)、結果を `[0,1]` にクランプ（丸め対策）。
- 数値安定性: `d→0`（同心、皆既/金環の最大食近傍）で `2dR, 2dr` の 0 除算を内包分岐で回避。`R≈r`（ハイブリッド境界）で判別式 `√(...)` の負値（丸め）を 0 クランプ。acos クランプ必須（accuracy.md §2.2）。禁止: クランプ無し acos（NaN）、境界の場合分け漏れ。
- 部分食地点で c2/c3=None でも食分・食面積は定義される（最大時点 m で評価。ISSUE-025/026 と整合）。

## 受け入れテスト
accuracy.md テストレベル **L1（純数学・円交差）＋L6（局地）**。基準値は実装へコピー禁止。
- 純幾何（L1、オラクル＝解析解）:
  - `d=0, R=r` → overlap=πR²（完全重なり、obscuration=1）。
  - `d ≥ R+r` → 0。
  - `d=|R−r|`（内接）→ overlap=π·min² 。
  - `R=r, d=R`（既知の lens 面積の解析値）と一致。
  - acos 引数が丸めで 1.0000001 になるケース → クランプで NaN を出さない。
- 食分（L6）: 中心線上 → magnitude ≥1（皆既）/ ≈1（金環）。部分食域 → 0<mag<1。離隔 → 0。
- オラクル（第二義・整合・accuracy.md §3.1）: NASA カタログ / USNO の地点別食分・食面積（data-sources §4）。k・ΔT 慣習を揃える（conventions §9）。fixtures 転記・出典/取得日明記（ISSUE-029）。
- 境界（accuracy.md §2.2）: 離隔/内包/外接/内接/同半径の 5 縮退すべてで有限値・場合分け正当。
- プロパティ（L8）: obscuration は d 単調減少で単調増加、`[0,1]` に収まる。magnitude と obscuration の整合（部分食で両者 0..1）。

## 許容誤差
accuracy.md §2.2「食分 ±0.0005／0.001食分≈1.9″」、§2.1（位置精度）から:
- **食分: ±0.0005**（相対位置 ≲1″、§2.2）。L1/L2/m の精度（ISSUE-021/024）に従属。
- **食面積（obscuration）: ±0.0005 相当**（食分同等の相対位置精度から。要確認: accuracy.md は食面積の独立許容を明記せず、§2.2 食分基準に準拠と解釈）。
- 純幾何（円交差）の数値一致: ≤ 1e-10（解析解に対し f64 丸めのみ）。
- 根拠: 食分・食面積は計算律速（accuracy.md §0(a)）。acos クランプ・境界処理で NaN/誤値を排除し、相対位置精度をそのまま反映。許容を通すためだけに拡大しない（conventions §11）。

## 実装メモ
- acos クランプは `x.clamp(-1.0, 1.0)`（accuracy.md §2.2 必須）。判別式の負値も 0 クランプ。
- 月半径 k 選択（conventions §9）で皆既/金環境界が動く。L1/L2 側で反映され本 issue は供給値依存。照合時は Espenak 慣習へ（conventions §9, accuracy.md §2.2）。系統差は metadata。
- 食分（基本面 Re）と食面積（視半径平面）の単位系の違いを doc で明示。混在禁止（conventions §1）。
- EclipseMagnitude は 1 超を許容（皆既）。Obscuration は厳密 `[0,1]`。
- レビュー重点: 5 縮退境界の網羅、acos/判別式クランプ、0 除算回避、k 慣習の反映箇所、単位系の分離。
