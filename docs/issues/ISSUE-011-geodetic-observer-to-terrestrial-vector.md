# ISSUE-011: Geodetic observer to terrestrial vector（観測者→ITRS/ECEF）

- crate: umbra-core
- 依存: ISSUE-001, ISSUE-002（角度）, ISSUE-003（`Vector3`）, ISSUE-010（WGS84・地心緯度）
- モード(tdd-workflow): strict（公開仕様＝観測者位置の地球固定ベクトル化。局地接触・高度方位・中心線の基礎で、誤差が直接位置誤差になるため strict）

## 目的
観測者（測地緯度・東経・楕円体高）を地球固定座標 ITRS/ECEF の `Position<Itrs>`（直交ベクトル）へ変換する（conventions §3/§4/§5）。
- 入力は `Observer { latitude: GeodeticLatitude, longitude: EastLongitude, elevation: Meters }`（api-draft §1.5）。
- 出力は ITRS の直交ベクトル（単位: km、conventions §1 幾何計算は km）。
- ベッセル無次元化（Re 単位）への橋渡しも提供可能にする。

## 非目的
- ITRS→TIRS→CIRS→GCRS のフレーム連鎖（極運動・ERA・歳差章動。conventions §5）。これは umbra-ephemeris の責務。本 issue は ITRS（地球固定）ベクトルの生成まで。
- 局地高度・方位の計算（umbra-eclipse）。
- 地心緯度変換そのもの（ISSUE-010）。本 issue はそれを利用してベクトル化。

## 公開インターフェース
api-draft §1.4/§1.5、conventions §5 に整合:
```rust
/// 観測者 → ITRS（地球固定・ECEF）直交ベクトル。単位 km。
pub fn observer_to_itrs(model: EarthModel, observer: Observer) -> Position<Itrs>;

/// 必要なら Re 無次元化版（ベッセル面投影用）。単位 Re。
pub fn observer_to_itrs_re(model: EarthModel, observer: Observer) -> Position<Itrs>;
```
- `Position<Itrs>`（api-draft §1.4: `Position<F: ReferenceFrame>`、`Itrs` フレーム型）。
- 東経正 `EastLongitude`、測地緯度 `GeodeticLatitude`、楕円体高 `Meters`（conventions §3/§4）。生 f64 で渡さない。

## 数式・アルゴリズムの出典
- **測地座標 → ECEF 直交座標**: 標準楕円体公式
  `X = (N + h)·cos φ·cos λ`,
  `Y = (N + h)·cos φ·sin λ`,
  `Z = (N(1−e²) + h)·sin φ`,
  ここで `N = a / √(1 − e²·sin²φ)`（卯酉線曲率半径）。出典: **IERS Conventions 2010 (TN 36) §4**、および測地学標準（Torge "Geodesy" / Hofmann-Wellenhof "GNSS")。
- 等価表現として **Meeus AA 2nd ed. Ch.11 の `ρsinφ'`/`ρcosφ'`**（ISSUE-010）から `Z = a·ρsinφ'`, `√(X²+Y²) = a·ρcosφ'` も可（実装はどちらか一方を正本とし、もう一方を検証に）。
- 東経正・λ は東経（conventions §3）。`Z` 軸 = 地球自転軸（ITRS、conventions §5）。
- 要確認: ITRS の x 軸（経度 0=グリニッジ）方向定義を IERS と一致させる（極運動適用前の TIRS との差は ISSUE-007/ephemeris で扱う）。

## 単位 / 時刻系 / 座標系
- 入力: 測地緯度 `GeodeticLatitude`、東経 `EastLongitude`、楕円体高 `Meters`。
- 出力: `Position<Itrs>`、成分単位 km（既定）/ Re（無次元化版）。
- 時刻系: なし（地球固定。時刻依存の回転は上位フレーム連鎖）。
- 座標系: ITRS（≈WGS84 軸、右手系、Z=自転軸、X=グリニッジ子午線方向。conventions §5）。

## アルゴリズム概要
1. `EarthModel`（ISSUE-010）から a, e² を取得。
2. `N = a / √(1 − e²·sin²φ)` を計算（φ=測地緯度）。
3. `X, Y, Z` を上式で計算（h=楕円体高、λ=東経）。単位は m → km へ境界変換（conventions §1）。
4. `Position<Itrs>` を構築（フレーム型で ITRS を明示）。
5. Re 版は m を Re（=WGS84 a）で無次元化（conventions §4）。
- 数値安定性: `√(1−e²sin²φ)` は φ 全域で安定。極でも `cos φ=0` を `atan2`/直接式で安全に。単位変換は 1 箇所に集約。禁止: 西経正の内部持ち込み（conventions §3）、地心/測地取り違え、m と km の混在。

## 受け入れテスト
accuracy.md テストレベル **L1（純数学）→ L6（局地条件）の前提**。
- 既知値（オラクル＝独立計算 / IERS 例 / 既知測地点の ECEF。実装からコピーしない）:
  - 赤道・λ=0・h=0 → `(a, 0, 0)` km（X=a, Y=Z=0）。
  - 北極・h=0 → `(0, 0, b)` km（b=極半径）。
  - λ=90°E・赤道 → `(0, a, 0)`（東経正の符号確認）。
  - 既知地点（例 api-draft 例の岡山 34.507°N, 133.508°E, 10 m）の ECEF を独立公式で照合。
  - ノルム `√(X²+Y²+Z²)` ≈ 地心距離（ISSUE-010 の ρ·a と一致）。
- 西経入力: `EastLongitude::from_signed_degrees(-120)` が東経 240° 相当のベクトルになる（conventions §3 の西経吸収）。
- 境界値: 極（±90°）、日付変更線（±180°）、高標高/負標高。
- プロパティ（L8）: 経度を λ→λ+2π でベクトル不変、緯度符号反転で Z 符号反転。

## 許容誤差
accuracy.md §2.3「観測者/楕円体/標高 sub-m」「WGS84 で十分」から:
- ECEF 位置: **sub-m（≤ 1 m、目標 ≪0.1 m）**。根拠: §2.3 が観測者位置 sub-m。δUT1 1 ms ≈ ~14 m に対し観測者起点誤差はそれより十分小さく。中心線 sub-km（§1 Standard）に効かないこと。
- 既知点 ECEF 一致: ≤ 0.01 m 目標（純幾何変換の丸めのみ）。
- Re 無次元化: 相対 ≤ 数 ULP。
- 根拠: 観測者位置誤差は接触時刻・高度方位に直結。§2.3 の sub-m を満たし、誤差を 2 桁余裕で確保。

## 実装メモ
- ITRS（極運動適用後・地球固定）と TIRS（極運動前）の区別: 本 issue は WGS84 軸の ITRS ベクトルを生成。極運動補正（xp,yp、ISSUE-007）の適用は ephemeris のフレーム連鎖側（conventions §5）で、ここでは行わない旨を doc 明記。
- 単位は km 既定（conventions §1 幾何計算）。Re 版はベッセル投影用（architecture §6、x,y,l1,l2 が Re 単位）。混在させない。
- 東経正・楕円体高（conventions §3/§4）。西経・正高入力は `EastLongitude::from_signed_degrees` / 別変換で境界吸収。
- レビュー重点: 東経正の符号、極/日付変更線境界、m↔km↔Re 単位境界の一元化、ITRS と TIRS の責務分離、地心/測地の取り違え無し。
