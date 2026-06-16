//! ベッセル要素の多項式供給源（NASA 形式, `docs/issues/ISSUE-022` スライス2）。
//!
//! [`BesselianSource`](crate::source::BesselianSource) を時間多項式で表す供給源
//! `BesselianPolynomial` を提供する。供給源を区間内で Chebyshev ノードサンプリング →
//! 多項式フィット（`crate::polynomial`）→ 残差ガード（[`BesselFitError`]）する層。
//! NASA 慣習（経路/エクスポート用）の x,y,d,μ,l1,l2 多項式 + 定数 tan f1/f2 を保持し、
//! `at()` 評価時に t = epoch_tt からの経過時間[hour]で各多項式を評価する。
//!
//! フィット手順（`docs/algorithms/07-bessel-polynomial.md`）: ① fit 区間を Chebyshev ノードで
//! サンプリング（`source.at`, ISSUE-037 直接が基準）② μ を時間順 unwrap で連続化 ③ 各成分を
//! Chebyshev 最小二乗 →monomial(t) フィット（`crate::polynomial`）④ 稠密サンプルで残差を実測
//! （[`BesselFitError`]）⑤ 許容超なら次数を上げて再フィット、上限超で
//! [`EclipseError::BesselFitExceededTolerance`]。

// 整数→f64/i32 変換は小さなノード数・サンプル数・次数・テストの多項式指数のみ（精度クリティカルな
// 天文量でない）。
#![allow(clippy::cast_precision_loss, clippy::cast_possible_truncation)]

use umbra_core::{Radians, TimeInterval, TtInstant};

use crate::besselian::InstantaneousBesselianElements;
use crate::error::EclipseError;
use crate::polynomial::{chebyshev_nodes, fit_chebyshev_monomial, Polynomial};
use crate::source::BesselianSource;

/// 自動エスカレーションの上限次数（NASA は 3 次。Runge/過剰適合を避け 6 で打切り, numerical-policy §A4）。
const MAX_AUTO_DEGREE: usize = 6;
/// 残差ガードの稠密サンプル区間数（fit 区間を等分する点数 = この値 + 1）。
const RESIDUAL_SAMPLES: usize = 64;
/// `at()` の区間端許容[hour]（端点を含めるための数値マージン）。
const INTERVAL_EPS_HOURS: f64 = 1.0e-6;

/// epoch からの経過時間[hour]（2要素 JD 差分で桁落ち回避, julian §days_since）。
fn elapsed_hours(time: TtInstant, epoch: TtInstant) -> f64 {
    time.jd2().days_since(epoch.jd2()) * 24.0
}

/// 経過時間[hour]に対応する TT 時刻。
fn tt_at_hours(epoch: TtInstant, hours: f64) -> TtInstant {
    TtInstant::from_jd2(epoch.jd2().add_days(hours / 24.0))
}

/// フィット残差（fit 区間での直接値と多項式値の最大絶対差, Re 無次元）。誤差を隠さず必ず保持。
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BesselFitError {
    /// x の最大絶対残差（Re）。
    pub max_x: f64,
    /// y の最大絶対残差（Re）。
    pub max_y: f64,
    /// l1 の最大絶対残差（Re）。
    pub max_l1: f64,
    /// l2 の最大絶対残差（Re）。
    pub max_l2: f64,
}

impl BesselFitError {
    /// 全成分が `tol` 以下か（許容判定）。
    fn within(&self, tol: &BesselFitError) -> bool {
        self.max_x <= tol.max_x
            && self.max_y <= tol.max_y
            && self.max_l1 <= tol.max_l1
            && self.max_l2 <= tol.max_l2
    }
}

/// ベッセル要素を時間多項式で表す供給源（NASA 形式, 経路/エクスポート用, ISSUE-022）。
///
/// 変数 t = `epoch_tt` からの経過時間[hour]。各成分多項式を Horner 評価して
/// [`InstantaneousBesselianElements`] を構成する。tan f1/f2 は定数（NASA 慣習）。
/// 直接評価（`DirectBesselianSource`, ISSUE-037, fit 誤差ゼロ）と `BesselianSource` 契約で差し替え可能。
#[derive(Clone, Debug)]
pub struct BesselianPolynomial {
    /// 多項式の基準時刻 t0（NASA 形式, 経過時間[hour]の原点）。
    pub epoch_tt: TtInstant,
    /// x(t)（Re）。
    pub x: Polynomial,
    /// y(t)（Re）。
    pub y: Polynomial,
    /// d(t)（影軸赤緯, rad）。
    pub d: Polynomial,
    /// μ(t)（連続化済み, rad。評価時に [0,2π) 正規化）。
    pub mu: Polynomial,
    /// l1(t)（半影錐半径, Re）。
    pub l1: Polynomial,
    /// l2(t)（本影錐半径, Re, 符号付き）。
    pub l2: Polynomial,
    /// tan f1（半影半頂角の正接, 定数）。
    pub tan_f1: f64,
    /// tan f2（本影半頂角の正接, 定数）。
    pub tan_f2: f64,
    /// fit 区間。
    pub fit_interval: TimeInterval<TtInstant>,
    /// fit 残差（必ず保持, 誤差を隠さない, conventions §11）。
    pub fit_error: BesselFitError,
}

impl BesselianPolynomial {
    /// 供給源 `source` を区間 `interval` で Chebyshev ノードサンプリングし、`degree` 次から
    /// 多項式フィットする。残差が `tolerance` 以下になる最小次数を採用（必要なら内部で次数を上げ、
    /// 上限 [`MAX_AUTO_DEGREE`] まで）。許容を満たせなければ
    /// [`EclipseError::BesselFitExceededTolerance`]。区間不正は [`EclipseError::InvalidFitInterval`]。
    pub fn fit(
        source: &impl BesselianSource,
        epoch_tt: TtInstant,
        interval: TimeInterval<TtInstant>,
        degree: usize,
        tolerance: BesselFitError,
    ) -> Result<Self, EclipseError> {
        let t_start = elapsed_hours(interval.start, epoch_tt);
        let t_end = elapsed_hours(interval.end, epoch_tt);
        if !t_start.is_finite() || !t_end.is_finite() || t_end <= t_start {
            return Err(EclipseError::InvalidFitInterval);
        }
        let t_center = 0.5 * (t_start + t_end);
        let half_width = 0.5 * (t_end - t_start);

        // tan f1/f2 は定数（NASA 慣習）: 区間中心の代表値を採る（定数 tan f はそのまま保存）。
        let center = source.at(tt_at_hours(epoch_tt, t_center))?;
        let tan_f1 = center.tan_f1;
        let tan_f2 = center.tan_f2;

        let max_degree = degree.max(MAX_AUTO_DEGREE);
        let mut last_achieved: Option<BesselFitError> = None;

        for deg in degree..=max_degree {
            let m = deg + 3; // ノード数 ≥ degree+1（余裕 +2）。
            let taus = chebyshev_nodes(m);

            // ① Chebyshev ノードでサンプリング。
            let mut node_t = Vec::with_capacity(m);
            let mut xs = Vec::with_capacity(m);
            let mut ys = Vec::with_capacity(m);
            let mut ds = Vec::with_capacity(m);
            let mut l1s = Vec::with_capacity(m);
            let mut l2s = Vec::with_capacity(m);
            let mut mus = Vec::with_capacity(m);
            for &tau in &taus {
                let th = t_center + half_width * tau;
                let elems = source.at(tt_at_hours(epoch_tt, th))?;
                node_t.push(th);
                xs.push(elems.x);
                ys.push(elems.y);
                ds.push(elems.declination.0);
                l1s.push(elems.l1);
                l2s.push(elems.l2);
                mus.push(elems.mu.0);
            }

            // ② μ を時間順に unwrap（連続化）。
            let mu_unwrapped = unwrap_mu(&node_t, &mus);

            // ③ 各成分を Chebyshev 最小二乗 →monomial(t) フィット。
            let x = fit_chebyshev_monomial(&taus, &xs, deg, t_center, half_width);
            let y = fit_chebyshev_monomial(&taus, &ys, deg, t_center, half_width);
            let d = fit_chebyshev_monomial(&taus, &ds, deg, t_center, half_width);
            let l1 = fit_chebyshev_monomial(&taus, &l1s, deg, t_center, half_width);
            let l2 = fit_chebyshev_monomial(&taus, &l2s, deg, t_center, half_width);
            let mu = fit_chebyshev_monomial(&taus, &mu_unwrapped, deg, t_center, half_width);

            // ④ 残差ガード: 稠密サンプルで x,y,l1,l2 の最大絶対残差を実測。
            let mut achieved = BesselFitError {
                max_x: 0.0,
                max_y: 0.0,
                max_l1: 0.0,
                max_l2: 0.0,
            };
            for k in 0..=RESIDUAL_SAMPLES {
                let frac = k as f64 / RESIDUAL_SAMPLES as f64;
                let th = t_start + (t_end - t_start) * frac;
                let direct = source.at(tt_at_hours(epoch_tt, th))?;
                achieved.max_x = achieved.max_x.max((x.eval(th) - direct.x).abs());
                achieved.max_y = achieved.max_y.max((y.eval(th) - direct.y).abs());
                achieved.max_l1 = achieved.max_l1.max((l1.eval(th) - direct.l1).abs());
                achieved.max_l2 = achieved.max_l2.max((l2.eval(th) - direct.l2).abs());
            }

            // ⑤ 許容内なら採用、超なら次数を上げて再フィット。
            if achieved.within(&tolerance) {
                return Ok(BesselianPolynomial {
                    epoch_tt,
                    x,
                    y,
                    d,
                    mu,
                    l1,
                    l2,
                    tan_f1,
                    tan_f2,
                    fit_interval: interval,
                    fit_error: achieved,
                });
            }
            last_achieved = Some(achieved);
        }

        Err(EclipseError::BesselFitExceededTolerance {
            achieved: last_achieved.unwrap_or(BesselFitError {
                max_x: f64::INFINITY,
                max_y: f64::INFINITY,
                max_l1: f64::INFINITY,
                max_l2: f64::INFINITY,
            }),
            tolerance,
        })
    }
}

/// μ サンプルを時間昇順に並べて 2π 折返しを除去し（連続化）、元の順序に戻す。
/// μ は赤経基準で [0,2π) 折返しがあるため、fit 前に区間内で unwrap する（conventions §2）。
fn unwrap_mu(node_t: &[f64], mus: &[f64]) -> Vec<f64> {
    use core::f64::consts::TAU;
    let n = mus.len();
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&a, &b| {
        node_t[a]
            .partial_cmp(&node_t[b])
            .expect("finite node times")
    });

    // 時間昇順の連続列を作る（隣接差の 2π 折返しを round で一括除去 = 隣接差を [−π,π] に収める）。
    let mut seq: Vec<f64> = order.iter().map(|&i| mus[i]).collect();
    for k in 1..seq.len() {
        let diff = seq[k] - seq[k - 1];
        let wraps = (diff / TAU).round(); // 何回 2π を跨いだか（増減両方向）。
        seq[k] -= wraps * TAU;
    }

    // 元の順序へ散らし戻す。
    let mut out = vec![0.0; n];
    for (k, &i) in order.iter().enumerate() {
        out[i] = seq[k];
    }
    out
}

impl BesselianSource for BesselianPolynomial {
    /// t = epoch_tt からの経過時間[hour]で各多項式を Horner 評価して瞬時要素を構成する。
    /// 区間外は [`EclipseError::EvaluationOutsideFitInterval`]。μ は [0,2π) 正規化。
    fn at(&self, time: TtInstant) -> Result<InstantaneousBesselianElements, EclipseError> {
        let t = elapsed_hours(time, self.epoch_tt);
        let t0 = elapsed_hours(self.fit_interval.start, self.epoch_tt);
        let t1 = elapsed_hours(self.fit_interval.end, self.epoch_tt);
        if t < t0 - INTERVAL_EPS_HOURS || t > t1 + INTERVAL_EPS_HOURS {
            return Err(EclipseError::EvaluationOutsideFitInterval);
        }
        let mu = Radians::new(self.mu.eval(t)).normalized_two_pi();
        Ok(InstantaneousBesselianElements {
            x: self.x.eval(t),
            y: self.y.eval(t),
            declination: Radians(self.d.eval(t)),
            mu,
            l1: self.l1.eval(t),
            l2: self.l2.eval(t),
            tan_f1: self.tan_f1,
            tan_f2: self.tan_f2,
            time_tt: time,
        })
    }

    fn fit_interval(&self) -> TimeInterval<TtInstant> {
        self.fit_interval
    }
}

#[cfg(test)]
mod tests {
    //! `BesselianPolynomial` の契約テスト（ISSUE-022 S2, 契約1-10）。
    //!
    //! 公開IF（`fit` / `at` / `fit_interval` / 各フィールド・新 `EclipseError` バリアント）と
    //! 契約のみに基づく外部検証。オラクルは「t の既知多項式を返す合成供給源 `PolySource`」と、
    //! 実 ephemeris（`DirectBesselianSource`, 2017-08-21）の 2 系統。
    //!
    //! 多項式評価はテスト内自前の素朴べき乗和（`naive_eval`）で行い、実装側 `Polynomial::eval` に
    //! 依存しない（オラクル独立性）。期待値はマジック直書きせず、合成源係数 / source.at から導出する。

    #![allow(clippy::excessive_precision)]

    use super::*;

    use crate::besselian::InstantaneousBesselianElements;
    use crate::error::EclipseError;
    use crate::source::{BesselianSource, DirectBesselianSource};

    use umbra_core::constants::{EARTH_EQUATORIAL_RADIUS_M, SOLAR_RADIUS_KM};
    use umbra_core::{EspenakMeeusDeltaT, JulianDate2, Radians, TimeInterval, TtInstant};

    /// 太陽物理半径[km]（実 ephemeris テスト用）。
    const R_SUN: f64 = SOLAR_RADIUS_KM;
    /// 月半径[km]（= k·Re, k=0.2725076）。
    const R_MOON: f64 = 0.2725076 * (EARTH_EQUATORIAL_RADIUS_M / 1000.0);

    // ------------------------------------------------------------------
    // 浮動小数比較ヘルパ
    // ------------------------------------------------------------------

    /// 絶対誤差比較。
    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    /// 角度を [0,2π) 正規化して絶対誤差比較（μ の境界折返しを跨いだ比較に使う）。
    fn close_angle(a: f64, b: f64, tol: f64) -> bool {
        let d = Radians::new(a - b).normalized_two_pi().0;
        d < tol || (std::f64::consts::TAU - d) < tol
    }

    // ------------------------------------------------------------------
    // 時刻ヘルパ
    // ------------------------------------------------------------------

    /// TT 時刻を 2 要素 JD から構築。
    fn tt(jd1: f64, jd2: f64) -> TtInstant {
        TtInstant::from_jd2(JulianDate2::new(jd1, jd2))
    }

    /// J2000.0（TT）を多項式の基準時刻 epoch として使う。
    fn epoch_j2000() -> TtInstant {
        tt(2_451_545.0, 0.0)
    }

    /// epoch から経過時間 `hours`[hour] だけ進んだ TT 時刻。
    /// 合成源 / 区間 / サンプル時刻を「経過時間（hour）」で素直に組むためのヘルパ。
    fn at_hours(epoch: TtInstant, hours: f64) -> TtInstant {
        TtInstant::from_jd2(epoch.jd2().add_days(hours / 24.0))
    }

    /// 契約1の t（= epoch からの経過時間[hour]）をテスト側で独立に再計算するオラクル。
    fn elapsed_hours(time: TtInstant, epoch: TtInstant) -> f64 {
        time.jd2().days_since(epoch.jd2()) * 24.0
    }

    /// epoch を中心に ±`half_hours`[hour] の fit 区間。
    fn interval_hours(epoch: TtInstant, half_hours: f64) -> TimeInterval<TtInstant> {
        TimeInterval {
            start: at_hours(epoch, -half_hours),
            end: at_hours(epoch, half_hours),
        }
    }

    // ------------------------------------------------------------------
    // 自前オラクル: 素朴べき乗和（実装側 Polynomial::eval に依存しない）
    // ------------------------------------------------------------------

    /// p(t) = Σ coeffs[i]·t^i（昇べき, 素朴和）。
    fn naive_eval(coeffs: &[f64], t: f64) -> f64 {
        let mut acc = 0.0;
        for (i, &c) in coeffs.iter().enumerate() {
            acc += c * t.powi(i as i32);
        }
        acc
    }

    // ------------------------------------------------------------------
    // 合成供給源 PolySource（テスト用クリーンオラクル）
    //
    // 各ベッセル成分を t（= epoch からの経過時間[hour]）の既知多項式で返す。
    // μ は連続な生多項式値を [0,2π) に正規化して返す（実供給源と同じ。fit 側が unwrap で復元）。
    // ------------------------------------------------------------------

    /// t の既知多項式で各ベッセル要素を供給する合成源。
    struct PolySource {
        /// この源の t 原点（fit の epoch_tt と一致させて使う）。
        epoch: TtInstant,
        /// この源の妥当区間（`fit_interval()` が返す。fit には別途 interval を渡す）。
        interval: TimeInterval<TtInstant>,
        /// x,y,d の t 多項式係数（昇べき）。
        cx: Vec<f64>,
        cy: Vec<f64>,
        cd: Vec<f64>,
        /// l1,l2 の係数。
        cl1: Vec<f64>,
        cl2: Vec<f64>,
        /// μ(t) の係数（連続な生 μ。0..2π 折返し前）。
        cmu: Vec<f64>,
        /// tan f（定数, NASA 慣習）。
        tan_f1: f64,
        tan_f2: f64,
    }

    impl PolySource {
        /// 典型的な皆既日食レンジを模した 2〜3 次の合成源（t=epoch 原点）。
        /// 係数はテスト側のオラクル（naive_eval）と完全に共有する。
        fn synthetic(epoch: TtInstant, interval: TimeInterval<TtInstant>) -> Self {
            PolySource {
                epoch,
                interval,
                // x,y は緩やかに動く（皆既帯の影軸交点を模す）。3次まで含めて非自明に。
                cx: vec![0.20, 0.45, -0.01, 0.002],
                cy: vec![-0.30, 0.18, 0.02, -0.001],
                // d（赤緯, rad）はほぼ一定 + 微小ドリフト。
                cd: vec![0.2070, 1.0e-4, -2.0e-6],
                // l1,l2（Re）。l2<0=皆既。
                cl1: vec![0.5400, 1.0e-4, 5.0e-6],
                cl2: vec![-0.0090, -3.0e-5, 1.0e-6],
                // μ(t): 連続。区間内では 2π を跨がない素直なケース（別途 unwrap 専用源を用意）。
                cmu: vec![1.2, 0.05, 1.0e-4],
                tan_f1: 0.004_65,
                tan_f2: 0.004_63,
            }
        }

        /// μ が区間内で 2π 境界を跨ぐ合成源（契約5: unwrap 検証用）。
        /// μ(t) = 6.0 + 0.30·t を 2π(≈6.283) 近傍から立ち上げ、区間を十分長く取り折返しを跨がせる。
        fn mu_wrapping(epoch: TtInstant, interval: TimeInterval<TtInstant>) -> Self {
            let mut s = Self::synthetic(epoch, interval);
            s.cmu = vec![6.0, 0.30];
            s
        }

        /// μ が時間とともに **減少**して 0/2π 境界を跨ぐ合成源（unwrap の逆方向折返し検証用）。
        /// μ(t) = 0.3 − 0.30·t は raw が負へ抜け、正規化で +2π 跳躍する（増加版と逆符号の wrap）。
        fn mu_wrapping_decreasing(epoch: TtInstant, interval: TimeInterval<TtInstant>) -> Self {
            let mut s = Self::synthetic(epoch, interval);
            s.cmu = vec![0.3, -0.30];
            s
        }

        /// 全成分が定数（degree=0 フィットの検証用）。
        fn constant(epoch: TtInstant, interval: TimeInterval<TtInstant>) -> Self {
            PolySource {
                epoch,
                interval,
                cx: vec![0.20],
                cy: vec![-0.30],
                cd: vec![0.2070],
                cl1: vec![0.5400],
                cl2: vec![-0.0090],
                cmu: vec![1.2],
                tan_f1: 0.004_65,
                tan_f2: 0.004_63,
            }
        }
    }

    impl BesselianSource for PolySource {
        fn at(&self, time: TtInstant) -> Result<InstantaneousBesselianElements, EclipseError> {
            let t = elapsed_hours(time, self.epoch);
            Ok(InstantaneousBesselianElements {
                x: naive_eval(&self.cx, t),
                y: naive_eval(&self.cy, t),
                declination: Radians(naive_eval(&self.cd, t)),
                // 実供給源と同じく [0,2π) 正規化して返す（fit 側が連続化する）。
                mu: Radians::new(naive_eval(&self.cmu, t)).normalized_two_pi(),
                l1: naive_eval(&self.cl1, t),
                l2: naive_eval(&self.cl2, t),
                tan_f1: self.tan_f1,
                tan_f2: self.tan_f2,
                time_tt: time,
            })
        }

        fn fit_interval(&self) -> TimeInterval<TtInstant> {
            self.interval
        }
    }

    /// 非多項式（三角）源: x=sin(ω t), y=cos(ω t)。サンプリング配置（t_center/half_width/
    /// 残差サンプル位置）が正しいときだけ低残差でフィットでき、配置を崩すと外挿で残差が悪化する。
    /// l1/l2/d/μ は素直（ガード対象の x,y で配置の正しさを縛る）。
    struct TrigSource {
        epoch: TtInstant,
        interval: TimeInterval<TtInstant>,
        omega: f64,
    }

    impl BesselianSource for TrigSource {
        fn at(&self, time: TtInstant) -> Result<InstantaneousBesselianElements, EclipseError> {
            let t = elapsed_hours(time, self.epoch);
            let a = self.omega * t;
            Ok(InstantaneousBesselianElements {
                x: a.sin(),
                y: a.cos(),
                declination: Radians(0.2070),
                mu: Radians::new(1.0 + 0.05 * t).normalized_two_pi(),
                l1: 0.5400,
                l2: -0.0090,
                tan_f1: 0.004_65,
                tan_f2: 0.004_63,
                time_tt: time,
            })
        }

        fn fit_interval(&self) -> TimeInterval<TtInstant> {
            self.interval
        }
    }

    /// 区間内に散らしたサンプル時刻（端点含む）を経過時間[hour]で生成。
    fn samples_in(epoch: TtInstant, half_hours: f64, n: usize) -> Vec<TtInstant> {
        (0..=n)
            .map(|k| {
                let frac = (k as f64) / (n as f64); // 0..1
                let h = -half_hours + 2.0 * half_hours * frac;
                at_hours(epoch, h)
            })
            .collect()
    }

    /// 緩い許容（合成源往復は機械精度級に厳しいが、係数導出残差の余裕を見込む）。
    fn tight_tol() -> BesselFitError {
        BesselFitError {
            max_x: 1.0e-6,
            max_y: 1.0e-6,
            max_l1: 1.0e-6,
            max_l2: 1.0e-6,
        }
    }

    // ==================================================================
    // 契約2 / 10 / 受け入れ「合成多項式源の往復（最重要）」
    // ==================================================================

    /// 契約2: 合成多項式源（各成分が t の低次多項式そのもの）を fit すると、
    /// `fit(..).at(t)` が区間内多数の t で source.at(t) に高精度一致し、fit_error ≈ 0。
    /// オラクルは合成源係数を共有する PolySource.at（= naive_eval）。マジック値直書きしない。
    #[test]
    fn fit_recovers_synthetic_polynomial_source() {
        let epoch = epoch_j2000();
        let half = 2.4; // ±2.4h
        let iv = interval_hours(epoch, half);
        let src = PolySource::synthetic(epoch, iv);

        // degree=3（合成源の最高次に一致）から fit。
        let poly = BesselianPolynomial::fit(&src, epoch, iv, 3, tight_tol())
            .expect("fit of a polynomial source must succeed");

        for time in samples_in(epoch, half, 40) {
            let got = poly.at(time).expect("at() inside interval must succeed");
            let want = src.at(time).expect("oracle source.at must succeed");
            assert!(
                close(got.x, want.x, 1e-9),
                "x: got={} want={}",
                got.x,
                want.x
            );
            assert!(
                close(got.y, want.y, 1e-9),
                "y: got={} want={}",
                got.y,
                want.y
            );
            assert!(
                close(got.declination.0, want.declination.0, 1e-9),
                "d: got={} want={}",
                got.declination.0,
                want.declination.0
            );
            assert!(close(got.l1, want.l1, 1e-9), "l1");
            assert!(close(got.l2, want.l2, 1e-9), "l2");
            assert!(
                close_angle(got.mu.0, want.mu.0, 1e-9),
                "mu: got={} want={}",
                got.mu.0,
                want.mu.0
            );
        }

        // fit_error は機械精度級（誤差を隠さないが、多項式源なので ~0）。
        assert!(
            poly.fit_error.max_x < 1e-9,
            "max_x={}",
            poly.fit_error.max_x
        );
        assert!(
            poly.fit_error.max_y < 1e-9,
            "max_y={}",
            poly.fit_error.max_y
        );
        assert!(
            poly.fit_error.max_l1 < 1e-9,
            "max_l1={}",
            poly.fit_error.max_l1
        );
        assert!(
            poly.fit_error.max_l2 < 1e-9,
            "max_l2={}",
            poly.fit_error.max_l2
        );
    }

    /// 契約1/2 系: epoch と区間中心がずれていても（t_center ≠ epoch）往復が成立する。
    /// 区間を [epoch+1h, epoch+5h]（中心=epoch+3h）に取り、t 原点は epoch のまま。
    #[test]
    fn fit_round_trip_with_offcenter_interval() {
        let epoch = epoch_j2000();
        let iv = TimeInterval {
            start: at_hours(epoch, 1.0),
            end: at_hours(epoch, 5.0),
        };
        let src = PolySource::synthetic(epoch, iv);
        let poly = BesselianPolynomial::fit(&src, epoch, iv, 3, tight_tol())
            .expect("off-center fit must succeed");

        // 区間 [1h,5h] を 30 分割して検証。
        for k in 0..=30 {
            let h = 1.0 + 4.0 * (k as f64) / 30.0;
            let time = at_hours(epoch, h);
            let got = poly.at(time).expect("at() must succeed in interval");
            let want = src.at(time).expect("oracle");
            assert!(close(got.x, want.x, 1e-9), "x at {h}h");
            assert!(close(got.y, want.y, 1e-9), "y at {h}h");
            assert!(close(got.l2, want.l2, 1e-9), "l2 at {h}h");
        }
    }

    /// サンプリング幾何（t_center / half_width / 残差サンプル位置）の健全性。
    /// 非多項式（緩い三角）源は、Chebyshev ノードと残差サンプルが fit 区間を正しく覆うときだけ
    /// 低残差で張れる。中心寄せ・スケール崩し・区間外サンプリングをすると外挿誤差が残差ゲートに現れ、
    /// 許容（1e-3）を満たせず fit が失敗する。つまり t_center/half_width/サンプル位置の演算が壊れると
    /// このテストが落ちる（多項式源では配置不変で検出できない幾何を縛る）。
    #[test]
    fn fit_trig_source_requires_correct_sampling_geometry() {
        let epoch = epoch_j2000();
        let half = 2.4;
        let iv = interval_hours(epoch, half);
        // ω: 区間 ±2.4h で偏角 ±1.5 rad（緩い曲率）。正配置の degree6 で残差 ≪1e-3、
        //    中心寄せ/区間外サンプリングだと外挿で残差 >1e-3。
        let src = TrigSource {
            epoch,
            interval: iv,
            omega: 1.5 / half,
        };

        let tol = BesselFitError {
            max_x: 1.0e-3,
            max_y: 1.0e-3,
            max_l1: 1.0e-3,
            max_l2: 1.0e-3,
        };
        let poly = BesselianPolynomial::fit(&src, epoch, iv, 6, tol)
            .expect("correct sampling geometry must fit a mild sine within 1e-3");

        // 正しい配置なら残差は許容内、かつ at() が源に一致（外挿でないことの確認）。
        for time in samples_in(epoch, half, 40) {
            let got = poly.at(time).unwrap();
            let want = src.at(time).unwrap();
            assert!(
                close(got.x, want.x, 2.0e-3),
                "x: got={} want={}",
                got.x,
                want.x
            );
            assert!(
                close(got.y, want.y, 2.0e-3),
                "y: got={} want={}",
                got.y,
                want.y
            );
        }
        // 残差ゲートが実際に小さい残差を報告している（幾何が正しい証拠）。
        assert!(
            poly.fit_error.max_x < 1.0e-3,
            "max_x={}",
            poly.fit_error.max_x
        );
        assert!(
            poly.fit_error.max_y < 1.0e-3,
            "max_y={}",
            poly.fit_error.max_y
        );
    }

    /// 残差ゲートのサンプル位置 `th = t_start + (t_end − t_start)·frac` の健全性。
    /// **非対称**区間 + 非多項式源では、サンプル span の演算が壊れる（例 `t_end − t_start` →
    /// `t_end + t_start`）とゲートが区間外を測り、正しく張れた多項式でも外挿領域の大残差で fit が失敗する。
    /// 対称区間では `t_end + t_start = 0` で検出できないため、非対称で縛る。
    #[test]
    fn fit_trig_asymmetric_interval_guards_residual_sample_span() {
        let epoch = epoch_j2000();
        // 非対称 [1h, 4h]。span 演算が `−`→`+` だとゲートは [1h, 1+(4+1)=6h] を測り、外挿で残差が悪化。
        let iv = TimeInterval {
            start: at_hours(epoch, 1.0),
            end: at_hours(epoch, 4.0),
        };
        let src = TrigSource {
            epoch,
            interval: iv,
            omega: 0.6,
        };
        let tol = BesselFitError {
            max_x: 1.0e-3,
            max_y: 1.0e-3,
            max_l1: 1.0e-3,
            max_l2: 1.0e-3,
        };
        let poly = BesselianPolynomial::fit(&src, epoch, iv, 6, tol)
            .expect("correct residual-sample span must fit within 1e-3 on an asymmetric interval");

        for k in 0..=40 {
            let h = 1.0 + 3.0 * (k as f64 / 40.0);
            let time = at_hours(epoch, h);
            let got = poly.at(time).unwrap();
            let want = src.at(time).unwrap();
            assert!(close(got.x, want.x, 2.0e-3), "x at {h}h");
            assert!(close(got.y, want.y, 2.0e-3), "y at {h}h");
        }
    }

    /// degree=0（定数フィット）でノード数 `m = deg + 3` が正しく確保されること。
    /// 定数源を degree=0 で fit → at() が定数を返す。`m` の算出が壊れて 0 ノードになると
    /// fit_chebyshev_monomial が NaN 多項式を返し（残差ゲートは `f64::max` が NaN を無視して 0 と誤報し
    /// 素通りするため）、at() が定数を返せず NaN になる → このテストが落ちる。
    #[test]
    fn fit_constant_source_at_degree_zero() {
        let epoch = epoch_j2000();
        let half = 2.4;
        let iv = interval_hours(epoch, half);
        let src = PolySource::constant(epoch, iv);
        let poly = BesselianPolynomial::fit(&src, epoch, iv, 0, tight_tol())
            .expect("degree-0 fit of a constant source must succeed");
        // degree=0 → 係数 1 個（定数多項式）。
        assert_eq!(
            poly.x.coefficients.len(),
            1,
            "degree-0 fit must yield 1 coefficient"
        );
        for time in samples_in(epoch, half, 10) {
            let got = poly.at(time).unwrap();
            let want = src.at(time).unwrap();
            assert!(got.x.is_finite(), "x must be finite (0 ノードだと NaN)");
            assert!(
                close(got.x, want.x, 1.0e-9),
                "x const: got={} want={}",
                got.x,
                want.x
            );
            assert!(close(got.l2, want.l2, 1.0e-9), "l2 const");
        }
    }

    /// `within` は **全成分** が許容以下を要求する（`&&`）。ある成分の許容を満たせない（max_x を
    /// 達成不能に小さく）と、他成分が緩くても fit は失敗する（`||` では誤って成功してしまう）。
    #[test]
    fn fit_requires_all_components_within_tolerance() {
        let epoch = epoch_j2000();
        let half = 2.4;
        let iv = interval_hours(epoch, half);
        // 非多項式源（三角）は x 残差が必ず有限・非ゼロ（~5e-5）。max_x を 1e-9 にすると達成不能。
        let src = TrigSource {
            epoch,
            interval: iv,
            omega: 1.5 / half,
        };
        // x のみ達成不能、他は緩い。全成分要求(&&)なら Err、いずれか(||)なら誤って Ok。
        let tol = BesselFitError {
            max_x: 1.0e-9,
            max_y: 1.0,
            max_l1: 1.0,
            max_l2: 1.0,
        };
        assert!(
            matches!(
                BesselianPolynomial::fit(&src, epoch, iv, 6, tol),
                Err(EclipseError::BesselFitExceededTolerance { .. })
            ),
            "unsatisfiable max_x must fail even if other components are within tolerance (&& not ||)"
        );
    }

    // ==================================================================
    // 契約1: t = epoch からの経過時間[hour]で各多項式を評価
    // ==================================================================

    /// 契約1: at() の time_tt ラベルは入力 time をそのまま保持し、
    /// tan_f1/tan_f2 は source の定数値と一致する（評価時刻に依らない）。
    #[test]
    fn at_preserves_time_label_and_constant_tan_f() {
        let epoch = epoch_j2000();
        let iv = interval_hours(epoch, 2.4);
        let src = PolySource::synthetic(epoch, iv);
        let poly = BesselianPolynomial::fit(&src, epoch, iv, 3, tight_tol()).unwrap();

        let time = at_hours(epoch, 1.1);
        let got = poly.at(time).expect("at() must succeed");
        assert_eq!(got.time_tt, time, "time_tt must echo the query time");
        // tan f は定数（NASA 慣習）として代表値が厳密保存される（close(_,_,0.0) は < 0 で恒偽のため
        // 厳密一致は assert_eq! で表す）。
        assert_eq!(got.tan_f1, src.tan_f1, "tan_f1 constant");
        assert_eq!(got.tan_f2, src.tan_f2, "tan_f2 constant");
    }

    /// 契約1: μ は at() で [0,2π) に正規化されて返る（境界跨ぎでない素直な源でも常に範囲内）。
    #[test]
    fn at_returns_mu_normalized_into_zero_two_pi() {
        let epoch = epoch_j2000();
        let iv = interval_hours(epoch, 2.4);
        let src = PolySource::synthetic(epoch, iv);
        let poly = BesselianPolynomial::fit(&src, epoch, iv, 3, tight_tol()).unwrap();

        for time in samples_in(epoch, 2.4, 20) {
            let mu = poly.at(time).unwrap().mu.0;
            assert!(
                (0.0..std::f64::consts::TAU).contains(&mu),
                "mu out of [0,2π): {mu}"
            );
        }
    }

    // ==================================================================
    // 契約3: 残差ガード（fit_error は誤差を隠さず有限・非負で同梱）
    // ==================================================================

    /// 契約3: 合成多項式源でも fit_error は有限・非負で結果に同梱される（NaN/Inf/負でない）。
    #[test]
    fn fit_error_is_finite_and_nonnegative_for_synthetic() {
        let epoch = epoch_j2000();
        let iv = interval_hours(epoch, 2.4);
        let src = PolySource::synthetic(epoch, iv);
        let poly = BesselianPolynomial::fit(&src, epoch, iv, 3, tight_tol()).unwrap();
        let e = poly.fit_error;
        for (name, v) in [
            ("max_x", e.max_x),
            ("max_y", e.max_y),
            ("max_l1", e.max_l1),
            ("max_l2", e.max_l2),
        ] {
            assert!(v.is_finite(), "{name} not finite: {v}");
            assert!(v >= 0.0, "{name} negative: {v}");
        }
    }

    // ==================================================================
    // 契約4: 許容超過 → Err(BesselFitExceededTolerance { achieved, tolerance })
    // ==================================================================

    /// 契約4: tolerance を厳密 0.0 にすると、実 ephemeris の有限な残差では満たせず
    /// `BesselFitExceededTolerance` を返す。achieved は実測残差を保持し、渡した tolerance と一致する。
    #[test]
    fn fit_exceeding_tolerance_errors_with_achieved() {
        let epoch = tt_2017_max();
        let half = 2.4;
        let iv = interval_hours(epoch, half);
        let dt = EspenakMeeusDeltaT;
        let src = DirectBesselianSource::new(R_SUN, R_MOON, &dt, iv);

        // 全成分 0.0 = どんな有限残差も超過する（実 fit 残差は非ゼロ）。
        let zero_tol = BesselFitError {
            max_x: 0.0,
            max_y: 0.0,
            max_l1: 0.0,
            max_l2: 0.0,
        };
        let err = BesselianPolynomial::fit(&src, epoch, iv, 3, zero_tol)
            .expect_err("zero tolerance must be exceeded by a real source");
        match err {
            EclipseError::BesselFitExceededTolerance {
                achieved,
                tolerance,
            } => {
                assert_eq!(tolerance, zero_tol, "tolerance must be echoed back");
                // achieved は実測残差: 有限・非負・少なくとも 1 成分は厳密 0 を超える。
                for v in [
                    achieved.max_x,
                    achieved.max_y,
                    achieved.max_l1,
                    achieved.max_l2,
                ] {
                    assert!(v.is_finite() && v >= 0.0, "achieved component invalid: {v}");
                }
                assert!(
                    achieved.max_x > 0.0
                        || achieved.max_y > 0.0
                        || achieved.max_l1 > 0.0
                        || achieved.max_l2 > 0.0,
                    "achieved must report a nonzero residual: {achieved:?}"
                );
            }
            other => panic!("expected BesselFitExceededTolerance, got {other:?}"),
        }
    }

    // ==================================================================
    // 契約5: μ 連続化（unwrap）
    // ==================================================================

    /// 契約5: μ が fit 区間内で 0/2π 境界を跨ぐ合成源でも fit が破綻せず、
    /// at() の μ が source（[0,2π) 正規化込み）と一致する。
    /// 区間を十分長く取り μ(t)=6.0+0.30·t が 2π を跨ぐようにする。
    #[test]
    fn fit_handles_mu_wrapping_across_two_pi() {
        let epoch = epoch_j2000();
        // half=6h ⇒ μ は [6.0+0.3·(-6), 6.0+0.3·6] = [4.2, 7.8] で 2π≈6.283 を跨ぐ。
        let half = 6.0;
        let iv = interval_hours(epoch, half);
        let src = PolySource::mu_wrapping(epoch, iv);

        // 区間内で実際に 2π 跨ぎが起きていることを前提として固定（テストの自己検証）。
        let mu_start = src.at(iv.start).unwrap().mu.0;
        let mu_end = src.at(iv.end).unwrap().mu.0;
        // 連続 μ では端で 4.2→7.8 だが、正規化後は折り返して mu_end < mu_start になる。
        assert!(
            mu_end < mu_start,
            "test setup must wrap across 2π: mu_start={mu_start} mu_end={mu_end}"
        );

        let poly = BesselianPolynomial::fit(&src, epoch, iv, 2, tight_tol())
            .expect("fit must not break across the 2π wrap");

        for time in samples_in(epoch, half, 40) {
            let got = poly.at(time).unwrap().mu.0;
            let want = src.at(time).unwrap().mu.0;
            assert!(
                (0.0..std::f64::consts::TAU).contains(&got),
                "mu out of [0,2π): {got}"
            );
            assert!(
                close_angle(got, want, 1e-7),
                "mu mismatch across wrap: got={got} want={want}"
            );
        }
    }

    /// 契約5（逆方向）: μ が時間とともに **減少** して 2π 境界を跨ぐ源でも unwrap が効く。
    /// 増加版と逆符号の wrap（round 補正の符号両方向）を踏ませる。
    #[test]
    fn fit_handles_mu_wrapping_decreasing() {
        let epoch = epoch_j2000();
        let half = 6.0;
        let iv = interval_hours(epoch, half);
        let src = PolySource::mu_wrapping_decreasing(epoch, iv);

        // 連続 μ は減少（0.3−0.3t）。正規化で raw<0 が +2π へ跳ぶので mu_end > mu_start になる
        //（増加版と逆。実際に跨ぎが起きることを自己検証）。
        let mu_start = src.at(iv.start).unwrap().mu.0;
        let mu_end = src.at(iv.end).unwrap().mu.0;
        assert!(
            mu_end > mu_start,
            "test setup must wrap downward across 2π: mu_start={mu_start} mu_end={mu_end}"
        );

        let poly = BesselianPolynomial::fit(&src, epoch, iv, 1, tight_tol())
            .expect("fit must not break across the downward 2π wrap");

        for time in samples_in(epoch, half, 40) {
            let got = poly.at(time).unwrap().mu.0;
            let want = src.at(time).unwrap().mu.0;
            assert!(
                (0.0..std::f64::consts::TAU).contains(&got),
                "mu out of [0,2π): {got}"
            );
            assert!(
                close_angle(got, want, 1e-7),
                "mu mismatch across downward wrap: got={got} want={want}"
            );
        }
    }

    // ==================================================================
    // 契約6: at() の区間チェック
    // ==================================================================

    /// 契約6: fit_interval の外（前後）で at() は EvaluationOutsideFitInterval、
    /// 区間内（端含む）は Ok。
    #[test]
    fn at_rejects_times_outside_fit_interval() {
        let epoch = epoch_j2000();
        let half = 2.4;
        let iv = interval_hours(epoch, half);
        let src = PolySource::synthetic(epoch, iv);
        let poly = BesselianPolynomial::fit(&src, epoch, iv, 3, tight_tol()).unwrap();

        // 端点（含む）は Ok。
        assert!(poly.at(iv.start).is_ok(), "start endpoint must be Ok");
        assert!(poly.at(iv.end).is_ok(), "end endpoint must be Ok");
        // 中央も Ok。
        assert!(poly.at(epoch).is_ok(), "interior must be Ok");

        // 区間直前 / 直後は Err(EvaluationOutsideFitInterval)。
        let before = at_hours(epoch, -(half + 0.5));
        let after = at_hours(epoch, half + 0.5);
        assert!(
            matches!(
                poly.at(before),
                Err(EclipseError::EvaluationOutsideFitInterval)
            ),
            "before interval must error with EvaluationOutsideFitInterval"
        );
        assert!(
            matches!(
                poly.at(after),
                Err(EclipseError::EvaluationOutsideFitInterval)
            ),
            "after interval must error with EvaluationOutsideFitInterval"
        );
    }

    // ==================================================================
    // 契約7: 不正区間 → Err(InvalidFitInterval)
    // ==================================================================

    /// 契約7: interval.start ≥ interval.end（経過時間で）の fit は InvalidFitInterval。
    /// start>end と start==end の両方を確認する。
    #[test]
    fn fit_rejects_invalid_interval() {
        let epoch = epoch_j2000();

        // start > end（経過時間で逆転）。
        let reversed = TimeInterval {
            start: at_hours(epoch, 2.0),
            end: at_hours(epoch, -2.0),
        };
        let src_r = PolySource::synthetic(epoch, reversed);
        assert!(
            matches!(
                BesselianPolynomial::fit(&src_r, epoch, reversed, 3, tight_tol()),
                Err(EclipseError::InvalidFitInterval)
            ),
            "start>end must error with InvalidFitInterval"
        );

        // start == end（退化区間）。
        let degenerate = TimeInterval {
            start: at_hours(epoch, 1.0),
            end: at_hours(epoch, 1.0),
        };
        let src_d = PolySource::synthetic(epoch, degenerate);
        assert!(
            matches!(
                BesselianPolynomial::fit(&src_d, epoch, degenerate, 3, tight_tol()),
                Err(EclipseError::InvalidFitInterval)
            ),
            "start==end must error with InvalidFitInterval"
        );
    }

    /// 契約7（非有限）: 区間端の経過時間が非有限（∞）でも `InvalidFitInterval`。
    /// start を −∞ 経過時間にすると `end ≤ start` では弾けず、`is_finite` ガード（`||` の別項）だけが
    /// 弾く。これで区間検査の各 `||` 項が独立に効くことを縛る。
    #[test]
    fn fit_rejects_nonfinite_interval() {
        let epoch = epoch_j2000();
        // start = −∞ 経過時間（end は有限・end≤start でない）→ is_finite ガードのみが弾く。
        let iv = TimeInterval {
            start: TtInstant::from_jd2(JulianDate2::new(f64::NEG_INFINITY, 0.0)),
            end: at_hours(epoch, 2.4),
        };
        let src = PolySource::synthetic(epoch, iv);
        assert!(
            matches!(
                BesselianPolynomial::fit(&src, epoch, iv, 3, tight_tol()),
                Err(EclipseError::InvalidFitInterval)
            ),
            "non-finite interval start must error with InvalidFitInterval"
        );
    }

    // ==================================================================
    // 契約8: epoch_tt / fit_interval 保持
    // ==================================================================

    /// 契約8: fit に渡した interval が fit_interval() と構造体フィールドに保持され、epoch_tt も保持される。
    #[test]
    fn fit_preserves_epoch_and_interval() {
        let epoch = epoch_j2000();
        let iv = interval_hours(epoch, 2.4);
        let src = PolySource::synthetic(epoch, iv);
        let poly = BesselianPolynomial::fit(&src, epoch, iv, 3, tight_tol()).unwrap();

        assert_eq!(poly.epoch_tt, epoch, "epoch_tt must be preserved");
        assert_eq!(poly.fit_interval, iv, "fit_interval field must equal input");
        assert_eq!(
            BesselianSource::fit_interval(&poly),
            iv,
            "fit_interval() accessor must equal input"
        );
    }

    // ==================================================================
    // 契約9: object-safe / &dyn 経由呼出（DirectBesselianSource と差替可能）
    // ==================================================================

    /// 契約9: `&dyn BesselianSource` 経由で BesselianPolynomial の at()/fit_interval() が呼べ、
    /// DirectBesselianSource と同じ trait オブジェクト型に束ねられる（差し替え可能）。
    #[test]
    fn usable_through_dyn_besselian_source() {
        let epoch = epoch_j2000();
        let iv = interval_hours(epoch, 2.4);
        let src = PolySource::synthetic(epoch, iv);
        let poly = BesselianPolynomial::fit(&src, epoch, iv, 3, tight_tol()).unwrap();

        let dt = EspenakMeeusDeltaT;
        let direct = DirectBesselianSource::new(R_SUN, R_MOON, &dt, iv);

        // 同一の &dyn 型に両供給源を束ねられる（差し替え可能性）。
        let sources: [&dyn BesselianSource; 2] = [&poly, &direct];
        for s in sources {
            // fit_interval() を trait オブジェクト越しに。
            assert_eq!(s.fit_interval(), iv, "fit_interval via &dyn");
        }

        // 多項式版の at() を &dyn 越しに呼べる（区間内）。
        let s: &dyn BesselianSource = &poly;
        let time = at_hours(epoch, 0.7);
        let got = s.at(time).expect("at() via &dyn must succeed");
        let want = src.at(time).expect("oracle");
        assert!(close(got.x, want.x, 1e-9), "x via &dyn");
        assert!(close(got.l2, want.l2, 1e-9), "l2 via &dyn");
    }

    // ==================================================================
    // 契約10: 次数エスカレーション
    // ==================================================================

    /// 契約10: 合成源が 3 次でも、開始 degree=1 から内部で必要な次数まで上げて tolerance を満たす。
    /// fit は成功し fit_error ≤ tolerance, at() は source に一致する。
    #[test]
    fn fit_escalates_degree_until_tolerance_met() {
        let epoch = epoch_j2000();
        let half = 2.4;
        let iv = interval_hours(epoch, half);
        let src = PolySource::synthetic(epoch, iv); // cx は 3 次成分を含む

        let tol = tight_tol();
        // 開始 degree=1（合成源の真の次数 3 より低い）でも成功しなければならない。
        let poly = BesselianPolynomial::fit(&src, epoch, iv, 1, tol)
            .expect("low starting degree must escalate and succeed");

        // fit_error ≤ tolerance（各成分）。
        assert!(
            poly.fit_error.max_x <= tol.max_x,
            "max_x {} > tol {}",
            poly.fit_error.max_x,
            tol.max_x
        );
        assert!(poly.fit_error.max_y <= tol.max_y, "max_y > tol");
        assert!(poly.fit_error.max_l1 <= tol.max_l1, "max_l1 > tol");
        assert!(poly.fit_error.max_l2 <= tol.max_l2, "max_l2 > tol");

        // at() が source に一致。
        for time in samples_in(epoch, half, 30) {
            let got = poly.at(time).unwrap();
            let want = src.at(time).unwrap();
            assert!(close(got.x, want.x, 1e-8), "x after escalation");
            assert!(close(got.y, want.y, 1e-8), "y after escalation");
            assert!(close(got.l2, want.l2, 1e-8), "l2 after escalation");
        }
    }

    // ==================================================================
    // 受け入れ: 実 ephemeris L7 残差（DirectBesselianSource, 2017-08-21）
    // ==================================================================

    /// 2017-08-21 最大食付近の TT epoch（besselian.rs / source.rs と同一）。
    fn tt_2017_max() -> TtInstant {
        tt(2_457_986.5, 0.768_532_2)
    }

    /// 受け入れ: 実 ephemeris（DirectBesselianSource）を ±2.4h で fit すると成功し、
    /// fit_error が緩い許容（1e-3 Re 程度）に収まる。at() が DirectBesselianSource.at() に
    /// fit_error 内で一致し、fit_error は非ゼロ（実源なので残差を隠さない）。
    #[test]
    fn fit_real_ephemeris_2017_has_small_residual() {
        let epoch = tt_2017_max();
        let half = 2.4; // ±2.4h = 皆既日食の窓
        let iv = interval_hours(epoch, half);
        let dt = EspenakMeeusDeltaT;
        let src = DirectBesselianSource::new(R_SUN, R_MOON, &dt, iv);

        // 緩い許容（実 fit 残差の桁を許す）。degree は高め（5）で開始。
        let loose = BesselFitError {
            max_x: 1.0e-3,
            max_y: 1.0e-3,
            max_l1: 1.0e-3,
            max_l2: 1.0e-3,
        };
        let poly = BesselianPolynomial::fit(&src, epoch, iv, 5, loose)
            .expect("real ephemeris fit within loose tolerance must succeed");

        // fit_error は有限・非負・許容内、かつ非ゼロ（実源は残差を隠さない）。
        let e = poly.fit_error;
        for (name, v, t) in [
            ("max_x", e.max_x, loose.max_x),
            ("max_y", e.max_y, loose.max_y),
            ("max_l1", e.max_l1, loose.max_l1),
            ("max_l2", e.max_l2, loose.max_l2),
        ] {
            assert!(v.is_finite() && v >= 0.0, "{name} invalid: {v}");
            assert!(v <= t, "{name} exceeds loose tol: {v} > {t}");
        }
        assert!(
            e.max_x > 0.0 || e.max_y > 0.0 || e.max_l1 > 0.0 || e.max_l2 > 0.0,
            "real-source fit_error must be nonzero: {e:?}"
        );

        // at() が直接供給源と fit_error 内で一致（区間内多数点）。
        // 許容は緩い許容の小さな倍率（多項式 vs 直接の差は fit_error 規模）。
        for time in samples_in(epoch, half, 40) {
            let got = poly.at(time).expect("at() inside interval");
            let want = src.at(time).expect("direct source");
            assert!(
                close(got.x, want.x, 2.0e-3),
                "x: got={} want={}",
                got.x,
                want.x
            );
            assert!(
                close(got.y, want.y, 2.0e-3),
                "y: got={} want={}",
                got.y,
                want.y
            );
            assert!(close(got.l1, want.l1, 2.0e-3), "l1");
            assert!(close(got.l2, want.l2, 2.0e-3), "l2");
            assert!(close_angle(got.mu.0, want.mu.0, 2.0e-3), "mu");
            assert!(close(got.declination.0, want.declination.0, 2.0e-3), "d");
        }
    }
}
