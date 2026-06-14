# ISSUE-022: Bessel polynomial fitting（時系列サンプリング→多項式近似→残差 BesselFitError 保持）

- crate: umbra-eclipse
- 依存: ISSUE-037（直接供給源 InstantaneousEvaluator/BesselianSource = サンプリング元）, ISSUE-021（瞬時要素定義）, ISSUE-008（数値解法基盤）, umbra-core（TtInstant, TimeInterval, Polynomial 候補）
- モード(tdd-workflow): strict（多項式は経路/エクスポートの値供給源（公開型 BesselianPolynomial, api-draft §3.3/§3.4）。fit 残差 BesselFitError を必ず保持し許容ガードする契約が精度保証の要。strict）

## 目的
最大食付近の時系列でベッセル要素をサンプリングし、各成分（x, y, d, μ, l1, l2）を時間多項式へ近似して `BesselianPolynomial` を生成する（architecture §6.1, api-draft §3.3）。
- NASA 形式（基準時刻 t0 からの時間 t を変数とする低次多項式）に合わせる。
- **残差 `BesselFitError` を必ず保持**し、許容を超えたら `BesselFitExceededTolerance`（api-draft §3.5, conventions「誤差を隠さない」）。
- `BesselianSource` を実装し、ISSUE-037（直接）と差し替え可能（経路/GeoJSON 用・高速, architecture §6.1）。

## 非目的
- 瞬時要素の定義・1時刻計算（ISSUE-021）・直接供給（ISSUE-037。本 Issue はこれを**サンプリング元**に使う）。
- 経路（中心線・限界線）生成そのもの（umbra-geo。本 Issue は経路が使う多項式を供給）。
- 多項式次数の固定（NASA 低次から開始し fit_error でガード, architecture §6.1）。

## 公開インターフェース
api-draft §3.3 / architecture §6.1 に準拠（公開型 ＋ BesselianSource 実装）。

```rust
#[derive(Clone, Debug)] pub struct Polynomial { pub coefficients: Vec<f64> }  // Horner 評価（architecture §12）
impl Polynomial { pub fn eval(&self, t: f64) -> f64; pub fn derivative(&self) -> Polynomial; }

#[derive(Clone, Copy, Debug)]
pub struct BesselFitError { pub max_x: f64, pub max_y: f64, pub max_l1: f64, pub max_l2: f64 }

#[derive(Clone, Debug)]
pub struct BesselianPolynomial {
    pub epoch_tt: TtInstant,                 // t0（NASA 形式の基準時刻）
    pub x: Polynomial, pub y: Polynomial, pub d: Polynomial, pub mu: Polynomial,
    pub l1: Polynomial, pub l2: Polynomial, pub tan_f1: f64, pub tan_f2: f64, // tan f は定数扱い（NASA 慣習）
    pub fit_interval: TimeInterval<TtInstant>,
    pub fit_error: BesselFitError,           // 残差を必ず保持
}

impl BesselianPolynomial {
    /// 直接供給源からサンプリングして fit。許容超過は Err。
    pub fn fit(
        source: &impl BesselianSource,       // ISSUE-037（直接）を基準にサンプリング
        epoch_tt: TtInstant,
        interval: TimeInterval<TtInstant>,
        degree: usize,                        // NASA 低次から（固定しない）
        tolerance: BesselFitError,            // 許容（accuracy.md §2.1 多項式 <0.10″ 相当）
    ) -> Result<Self, EclipseError>;          // BesselFitExceededTolerance
}

impl BesselianSource for BesselianPolynomial {
    fn at(&self, time: TtInstant) -> Result<InstantaneousBesselianElements, EclipseError>; // Horner 評価
    fn fit_interval(&self) -> TimeInterval<TtInstant>;
}
```

## 数式・アルゴリズムの出典
- **NASA ベッセル多項式形式**: Espenak の Besselian Elements（NASA GSFC / NASA TP-2006-214141, data-sources §4.1）。要素を `t0` からの時間 t（時間単位は時 hour が NASA 慣習, 要確認）の多項式 `x = x0 + x1·t + x2·t² + x3·t³` で表す（通常 3 次）。tan f1, tan f2 は定数（NASA 慣習）。μ は t の1次が主（地球自転）。
- **多項式 fit**: 最小二乗（normal equations or QR）。等間隔サンプル＋低次。出典: 標準的な多項式最小二乗（Numerical Recipes 多項式 fit / polyfit）。Horner 評価（architecture §12）。
- **残差評価**: fit 区間でサンプル点・中間点の直接値（ISSUE-037）と多項式値の最大絶対差を `BesselFitError`（x,y,l1,l2）に保持。出典: accuracy.md §3 実測ガード。
- μ の連続化: μ は赤経基準で `[0,2π)` 折返しがあるため、fit 前に区間内で unwrap（連続化, conventions §2）。

## 単位 / 時刻系 / 座標系
- 要素は ISSUE-021 と同一: **TT 基準**（conventions §6, epoch_tt 起点）、**FundamentalPlane**（conventions §5）、x,y,l1,l2 **Re 無次元**、d,μ ラジアン。
- 多項式変数 t = epoch_tt からの経過時間（単位は NASA 慣習＝時 hour 候補。実装で固定しコメント, conventions §10）。
- 供給源差し替え（ISSUE-037 と）で単位/座標系不変（同一 BesselianSource 契約）。

## アルゴリズム概要
1. `fit_interval` を等間隔サンプリング（点数は次数+余裕。最大食を中心に対称配置）。
2. 各サンプルで `source.at(t)`（ISSUE-037 直接, fit 誤差ゼロ基準）を評価し x,y,d,μ,l1,l2 を取得。μ は unwrap 連続化。
3. 各成分を低次（NASA: 3 次から）多項式へ最小二乗 fit。tan f1/f2 は区間平均 or 定数。
4. fit 区間で密にサンプルし残差最大値を測定 → `BesselFitError`。
5. 残差が `tolerance` 超なら次数を上げて再 fit（上限あり）、それでも超なら `BesselFitExceededTolerance`。
6. `BesselianPolynomial` を返す。`at` は Horner 評価（architecture §12）。
- 数値安定性: μ の折返し unwrap 必須（連続化, conventions §2）。Vandermonde の条件数悪化を避けるため t をスケーリング（区間を [-1,1] 正規化等）。許容を「通すため」に緩めない（fit_error で必ずガード, conventions §11）。

## 受け入れテスト
accuracy.md テストレベル **L4（ベッセル）** ＋ **L7 サブテスト（直接 vs 多項式残差, accuracy.md §3.2, architecture §6.1）**。
- **直接 vs 多項式 残差テスト（最重要, L7）**: ISSUE-037（直接）を基準に、`BesselianPolynomial.at(t)` の x,y,l1,l2 残差を fit 区間で実測。`fit_error` が実残差を正しく報告し、許容超で `BesselFitExceededTolerance` を返すこと。
- **NASA 係数比較（第二義, data-sources §4.1, 品質基準「係数比較＋瞬時値比較の両方」）**: 既知日食で生成多項式係数（x0,x1,x2,x3,...）を NASA 公開ベッセル多項式と比較（k/ΔT 慣習を揃える, accuracy.md §3.1）。係数比較と、評価した瞬時値比較の両方を行う。
- **MockEphemeris（accuracy.md §3.1）**: 人工配置でサンプリング→fit→評価が直接値に一致（fit 誤差が許容内）。
- fit_error 保持テスト: `BesselFitError` が必ず非ゼロで埋まり、結果に同梱されること（誤差を隠さない, conventions §11）。
- 許容超ケース: fit 区間を過大にして残差を悪化させ `BesselFitExceededTolerance` を確認。
- μ 連続化テスト: μ が ±π/2π 境界をまたぐ区間で unwrap が効き fit が破綻しないこと（L1）。

## 許容誤差
- accuracy.md §2.1「多項式 fit（使う場合）**<0.10″**」。fit 区間で残差を実測ガード（§3）。x,y の 0.10″ 相当を Re 換算した値を `tolerance` 既定に（要 Milestone 2 実測, accuracy.md §2.1 注）。
- **直接（ISSUE-037, fit誤差0）vs 多項式の残差 = L7 サブテスト**（accuracy.md §3.2）。profile 毎の採用（**経路/エクスポート=多項式が本務、局地=直接(037)既定**, B2）を残差で最終決定（architecture §6.1）。局地への多項式転用は「**M2 で fit 残差 <0.01″ を実測実証後のバッチ局地（多地点）の任意最適化**・未達はフォールバック」に限定。
- 食分 0.001≈1.9″（accuracy.md §2.2）に対し l1,l2 の fit 残差は十分小さく（<0.10″ 相当）。
- 許容を通すための拡大禁止（conventions §11）。fit_error は常に真の残差を報告。

## 実装メモ
- **責務分担（品質基準・B2）**: ISSUE-021=定義/1時刻 / ISSUE-037=直接供給（**Standard 局地の既定**, fit 誤差ゼロ）/ 本 Issue=多項式供給（**経路/GeoJSON/NASA エクスポートが本務**）。両者とも `BesselianSource`（architecture §6.1）。本 Issue の**局地転用は M2 で fit 残差 <0.01″ を実測実証後のバッチ局地（多地点）の任意最適化に限定**し、未達は直接(037)へフォールバック。実証前は局地で多項式を既定にしない。
- サンプリング元は ISSUE-037（直接, fit誤差ゼロ）を基準にする（誤差を直接 vs 多項式に分離, accuracy.md §4）。暦を直接サンプリングしない（層分解のため）。
- 次数は NASA 低次（3 次）から開始し fit_error でガード、必要時のみ上げる（architecture §6.1, 固定しない）。tan f1/f2 を定数にするか低次にするかは NASA 慣習＋残差で決定（要レビュー）。
- 多項式変数 t の時間単位（時 or 日 or 秒）と t0 起点を固定しコメント（NASA は時単位が多い, 要確認）。
- μ の unwrap・t スケーリング（条件数）を忘れると fit が破綻。レビュー重点。
- レビュー重点: fit_error が真の残差を報告、許容ガード（BesselFitExceededTolerance）、直接 vs 多項式の L7 残差、NASA 係数＋瞬時値の両比較、μ 連続化。
