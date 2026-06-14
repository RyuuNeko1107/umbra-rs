# ISSUE-024: Local observer projection（観測者を ITRS→基本面座標へ射影・各時刻）

- crate: umbra-eclipse
- 依存: ISSUE-001（規約）, ISSUE-003（`Vector3`）, ISSUE-011（観測者→ITRS/ECEF ベクトル）, ISSUE-021（瞬時ベッセル要素 x,y,d,μ,l1,l2,tan f）, ISSUE-020（基本面基底）, ISSUE-039（時角・恒星時供給: ERA(iauEra00)経由 CIO ベース。分点 GST 禁止・D4）, ISSUE-007（UT1/ΔT・UTC↔TT。恒星時は ISSUE-039 へ移管）
- モード(tdd-workflow): strict（局地接触・食分・高度方位すべての基礎。射影誤差が ξ,η,ζ の位置誤差→接触時刻・食分誤差へ直結するため strict）

## 目的
観測者（ITRS ベクトル）を、各時刻の月影軸を z 軸とする **ベッセル基本面 (FundamentalPlane)** の座標 `(ξ, η, ζ)` へ射影する（conventions §5, architecture §6）。接触・最大食・食分 solver が時刻関数として呼ぶ最下層プリミティブ。
- 入力: `Observer`（ITRS ベクトル化は ISSUE-011）と、ある TT における瞬時ベッセル要素（赤緯 d・時角 μ。ISSUE-021）。
- 出力: 基本面座標 `(ξ, η, ζ)`（単位 Re、conventions §1/§5）と、必要な時間微分 `(ξ', η')`（接触 solver の感度・継続時間に使用）。
- 「各時刻」= 時刻 t を与えると観測者の基本面座標を返す純関数として実装（solver が繰り返し評価）。

## 非目的
- 接触のゼロ点求解（ISSUE-025）、最大食最小化（ISSUE-026）、食分（ISSUE-027）。本 issue は座標射影のみ。
- 高度・方位（ISSUE-028）。ただし射影に使う赤緯・時角・恒星時は共有する（重複計算を避ける設計メモを残す）。
- 瞬時ベッセル要素そのものの生成（ISSUE-021 / 供給源 ISSUE-037 `BesselianSource`・architecture §6.1）。本 issue は供給された d, μ を使う。
- 観測者→ITRS 変換（ISSUE-011）。本 issue はその出力（Re 無次元化版）を入力とする。

## 公開インターフェース
architecture §6/§6.1、api-draft §1.4/§3.3 に整合。基本面射影は局地計算の内部プリミティブ（pub(crate) 中心、検証のため pub も検討）:
```rust
/// 観測者の基本面座標（月影軸 z のベッセル系）。単位 Re。
#[derive(Clone, Copy, Debug)]
pub struct ObserverFundamental {
    pub xi: f64, pub eta: f64, pub zeta: f64,   // 単位 Re（conventions §1/§5）
    pub xi_rate: f64, pub eta_rate: f64,        // d/dt（D1/conventions §1: 時間単位は SI秒基準で固定 = Re/SI秒）
}

/// ある TT における観測者の基本面座標へ射影。
/// elements は ISSUE-021 の瞬時ベッセル要素（d, μ を含む）。
pub fn project_observer_to_fundamental(
    observer_itrs_re: Position<Itrs>,   // ISSUE-011 の Re 無次元化版
    elements: &InstantaneousBesselianElements,
) -> ObserverFundamental;
```
- `InstantaneousBesselianElements`（api-draft §3.3）から `declination`(d), `hour_angle`(μ) を使用。
- `Position<Itrs>`（api-draft §1.4）。生 f64 を「単位付き量」として渡さない（conventions §1）。
- ζ は観測者が基本面より月側か地球内側かの符号判定に使う（接触の幾何条件・高度計算）。

## 数式・アルゴリズムの出典
- **Explanatory Supplement to the Astronomical Almanac (Seidelmann ed., 3rd ed.), §11 "Eclipses"（旧 §8）** の局地予報式。観測者地心座標 `(ρ sinφ', ρ cosφ')` と恒星時から基本面座標を作る標準形:
  - `ξ = ρ cosφ' · sin(μ − λ)`（要確認: 文献により H = μ + λ_east の符号定義が異なる。東経正・conventions §3 に合わせて固定し実装コメントに式番号転記）
  - `η = ρ sinφ' · cos d − ρ cosφ' · sin d · cos(μ − λ)`
  - `ζ = ρ sinφ' · sin d + ρ cosφ' · cos d · cos(μ − λ)`
  ここで d=月影軸赤緯、μ=エフェメリス時角（基本面の x 軸方向と観測点の関係。**D4: ERA(iauEra00)経由 CIO ベース・ISSUE-039 供給、分点 GST 禁止**）、λ=観測者東経、`ρ sinφ'/ρ cosφ'`=地心緯度成分（ISSUE-010/011）。
- **地球扁平（oblateness）の効き**: ξ, η, ζ には地心緯度成分 `ρ sinφ', ρ cosφ'` を通して **WGS84 楕円体の扁平**が直接効く（測地緯度 φ → 地心緯度 φ' の差・ρ の緯度依存）。地心/測地の取り違えや球近似は ξ,η,ζ → 接触時刻・食分の系統誤差になるため、`ρ sinφ', ρ cosφ'` は ISSUE-010/011 の WGS84 値を用い、扁平を含む点を実装コメントに明記する。
- **Meeus, Astronomical Algorithms 2nd ed., Ch.54 "Solar Eclipses"（式 54.1〜）** の `u, v`（=ξ, η）構成式。本式は Meeus の局地予報手順と等価。出典は章・式番号まで実装コメントに記す。
- 時間微分 `(ξ', η')`: ベッセル要素の時間変化（x', y', d', μ'）と地球自転から解析微分。出典: Explanatory Supplement の継続時間導出（ξ'≈ d/dt(ξ−x) に使う成分）。**時間単位は SI秒基準で固定（Re/SI秒、conventions §1 整合）**。μ の時間変化率は ≈ 地球自転角速度 + 影軸赤経変化（D4: μ は ERA(iauEra00)経由 CIO ベース・ISSUE-039 供給）。符号は要確認だが単位は SI秒で固定する。

## 単位 / 時刻系 / 座標系
- 入力: 観測者 ITRS（Re 無次元化、ISSUE-011）、瞬時ベッセル要素（d=`Radians`, μ=`Radians`、x,y,l1,l2=Re。架構 §6）。
- 出力: `(ξ, η, ζ)` 単位 Re（conventions §1/§5、ベッセル無次元化）。**微分の時間単位は SI秒基準で固定（Re/SI秒、conventions §1 整合）**。Meeus 慣習（分）で式を引く場合も内部は SI秒に換算して保持する。
- 時刻系: TT（瞬時要素は TT 基準。conventions §6）。恒星時・μ は UT1 由来（ISSUE-007）。時刻系の橋渡しは ISSUE-021/037 側で済み、本 issue は要素を受け取る。
- 座標系: FundamentalPlane（月影軸 z、conventions §5）。x 軸=春分点方向の影軸赤経基準（Explanatory Supplement 定義）。

## アルゴリズム概要
1. 観測者 ITRS（Re）から地心緯度成分 `ρ sinφ', ρ cosφ'` を取得（ISSUE-010/011。標高込み）。
2. 局地時角 `H = μ − λ_east`（東経正・conventions §3。符号は文献に合わせ固定）。
3. 上式で `ξ, η, ζ` を計算。
4. 微分 `ξ', η'` を解析式 or 与えられた要素微分から計算（接触・継続時間で使用）。
5. `ObserverFundamental` を返す。
- 数値安定性: 三角関数は φ', d 全域で安定。極・日付変更線でも `sin/cos(μ−λ)` は連続。μ−λ は `[-π,π)` の signed 正規化（conventions §2、折返しで微分が壊れないよう連続化）。禁止: 西経正の内部持ち込み（conventions §3）、地心/測地取り違え、km/Re 混在。
- 部分食地点で c2/c3=None になる扱いは ISSUE-025/027 の責務。本 issue は ζ 符号と (ξ,η) を返すのみ。

## 受け入れテスト
accuracy.md テストレベル **L6（局地条件）** の最下層。基準値は実装からコピーしない。
- MockEphemeris（accuracy.md §3.1）の人工配置で幾何検証:
  - 完全中心配置: 影軸直下の観測者で `ξ=η≈0`（中心線上）。
  - 影軸から既知の横ずれ点で `√(ξ²+η²)` が独立計算と一致。
- 既知地点（api-draft 例 岡山 34.507°N,133.508°E,10m）で、独立に手計算した `(ξ,η,ζ)` と照合（オラクル＝Explanatory Supplement 式の独立実装 or 公開ワークシート、出典・取得日明記）。
- 地点分類（accuracy.md L6）: 中心線上 / 付近 / 北南限 / 部分食域 / 可視域外 / 標高差（h=0 と h=4000m で ζ が標高分変化）。
- 西経入力: `EastLongitude::from_signed_degrees(-100)` が東経 260° 相当と一致（conventions §3 西経吸収）。
- プロパティ（L8）: λ→λ+2π で不変。μ→μ+2π で不変。
- 微分テスト: `ξ', η'` を中心差分（数値微分）と比較し一致（ISSUE-009 系の差分幅で）。
- **D1 ζ 誤差・月距離→視差 感度テスト（局地 topocentric バジェット）**: ζ（観測者の基本面前後位置, Re）の誤差、および月距離 r の誤差 δr/r → 地平視差 **π ≈ 3422″·(δr/r)** が、局地接触時刻に与える感度を実測する。ζ を微小に振り、また視差相当のずれを (ξ,η,ζ) に注入して局地接触時刻（ISSUE-025 連携 or 本層の m=√(ξ²+η²) 経路）の変化を測り、局地バジェット（別建て RSS, accuracy.md D1）の許容内に収まることを確認。月距離→視差→局地時刻の経路が地心量と別に効く点を fixture 化。

## 許容誤差
accuracy.md §2.3「観測者/楕円体/標高 sub-m」「WGS84 で十分」、§2.1（最大食 ±1.5s の幾何バジェット）から:
- `(ξ,η,ζ)` 位置: **sub-m 相当（Re 換算で ≤ 1m/6.378e6 ≈ 1.6e-7 Re、目標 ≪0.1m）**。根拠: 観測者起点誤差は接触時刻・高度方位へ直結し、§2.3 の sub-m を満たすこと。
- 既知点一致: 純幾何変換の丸めのみ（≤ 1e-9 Re 目標）。
- 微分: 数値微分との相対一致 ≤ 1e-6（接触感度・継続時間に効く。§2.1 solver バジェット 0.05″ を侵さない）。
- 根拠: 本層は「計算律速」(accuracy.md §0(a))。射影誤差で食分 ±0.0005・接触 ±2s（§2）の余裕を食わない。

## 実装メモ
- 高度/方位（ISSUE-028）と恒星時・赤緯・時角を共有する。重複評価を避けるため d, μ, 恒星時を一度だけ計算し両者へ渡す設計を doc に記す。
- ζ 符号: 月側（基本面前方）か地球内部かの判定に使う。接触幾何（ISSUE-025）の前提なので符号規約を実装コメントに固定（conventions §10）。
- μ の符号・時角定義（H = μ − λ か μ + λ か）は Explanatory Supplement / Meeus Ch.54 のどちらを正本にするか決め、もう一方で交差検証。要確認事項として明記。
- 単位は Re（conventions §1）。km/m との混在禁止。**微分の時間単位は SI秒基準で 1 箇所に固定（Re/SI秒、conventions §1 整合）**。分単位は内部に持ち込まない（境界で換算）。
- レビュー重点: 時角符号、東経正、地心/測地、Re 単位、μ−λ の連続化、微分の符号と単位。
