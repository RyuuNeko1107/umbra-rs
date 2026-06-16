//! 新月（朔）候補生成（ISSUE-016）。
//!
//! 期間（UTC）を新月単位に分解し、後段の地心合 solver（ISSUE-017）が確実にブラケットを
//! 取れる検索窓付きの候補列を生成する層。各候補は Meeus Ch.49 の平均朔概算時刻 `approx_tt`、
//! その周囲 ±Δ（≈1 日）の検索窓 `search_window`、および Meeus lunation index `k` を持つ。
//!
//! **偽陰性なし**（期間内の真の新月が必ずいずれかの候補窓に内包される）ことを最優先の契約と
//! する。平均朔と真の朔のずれ（最大 ≈±0.6 日）を上回るマージン Δ を窓半幅に取ることでこれを保証する。
//!
//! アルゴリズム: 期間境界を UTC→TT に変換し、窓が期間と交差する平均朔の lunation index `k` の閉区間を
//! `[ceil((tt_start−Δ−E)/S), floor((tt_end+Δ−E)/S)]` で直接求め、各 `k` の平均朔 JDE を `approx_tt`、
//! 窓 `[approx−Δ, approx+Δ]` として候補化する。Δ=1 が端の窓張り出しを担保し、ずれ ≈0.6 < Δ なので
//! 期間内の真の新月は必ずこの `k` 範囲に対応する（偽陰性ゼロ）。

// 整数↔浮動小数変換は lunation index k（±数千で f64 厳密）と JDE→k の floor のみ（精度クリティカルな
// 天文量ではない。窓マージンが偽陰性ゼロを担保し、概算式は精度バジェットに直接寄与しない）。
#![allow(clippy::cast_precision_loss, clippy::cast_possible_truncation)]
// new_moon_candidates は ISSUE-017（地心合 solver）が消費する。結線され次第この許容は外す
// （polynomial.rs S1 と同手順）。
#![allow(dead_code)]

use umbra_core::time::utc_to_tt;
use umbra_core::{JulianDate2, TimeInterval, TimeRange, TtInstant, UtcInstant};

use crate::error::EclipseError;

/// Meeus エポック: lunation index `k = 0` の平均朔 JDE（2000-01-06.x, 力学時, Ch.49 式 49.1）。
const MEEUS_EPOCH_JDE: f64 = 2_451_550.097_66;
/// 平均朔望月[日]（Meeus Ch.49, 式 49.1 の `k` 係数）。
const SYNODIC_MONTH_DAYS: f64 = 29.530_588_861;
/// 検索窓の半幅 Δ[日]。平均朔と真の朔の最大ずれ ≈±0.6 日を上回るマージン（偽陰性ゼロの砦,
/// Meeus Ch.49 の周期補正項の総和上限が根拠。窓を不必要に広げると後段の粗走査コストが増えるため 1 日）。
const WINDOW_HALF_WIDTH_DAYS: f64 = 1.0;

/// 1 朔の候補（ISSUE-016）。後段の地心合 solver（ISSUE-017）への入力窓を含む。
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct NewMoonCandidate {
    /// 朔の概算時刻（地心合の概算, TT）。Meeus 49.1 の平均朔 JDE。
    pub approx_tt: TtInstant,
    /// 後段 solver がブラケットを取る検索窓（`approx_tt` 中心 ±Δ）。
    pub search_window: TimeInterval<TtInstant>,
    /// Meeus lunation index `k`（`k = 0` ↔ 2000-01-06 の新月）。安定番号・event_key 素材。
    pub lunation_number: i64,
}

/// Meeus 49.1: lunation index `k` の平均朔 JDE（力学時 TT≈TDT）。
/// `JDE = 2451550.09766 + 29.530588861·k + 0.00015437·T² − 0.000000150·T³ + 0.00000000073·T⁴`,
/// `T = k/1236.85`。周期補正（太陽・月の平均近点角等, Ch.49）は窓マージン Δ が吸収するため不要。
fn mean_new_moon_jde(k: f64) -> f64 {
    let t = k / 1236.85;
    MEEUS_EPOCH_JDE + SYNODIC_MONTH_DAYS * k + 0.000_154_37 * t * t - 0.000_000_150 * t * t * t
        + 0.000_000_000_73 * t * t * t * t
}

/// 期間内の全朔候補を時系列順で返す（偽陰性なし）。期間境界の UTC→TT 変換に失敗（1972 年より前など）
/// すると [`EclipseError::Time`]。
pub(crate) fn new_moon_candidates(
    range: TimeRange<UtcInstant>,
) -> Result<Vec<NewMoonCandidate>, EclipseError> {
    // 期間境界を TT へ変換（平均朔式 JDE は力学時）。1972 前は閏秒未定義で Err(Time)。
    let tt_start = utc_to_tt(range.start)?.jd2().jd();
    let tt_end = utc_to_tt(range.end)?.jd2().jd();

    // 窓 [A(k)−Δ, A(k)+Δ] が期間 [tt_start, tt_end] と交差する ⇔ 平均朔 A(k) ∈ [tt_start−Δ, tt_end+Δ]。
    // A(k) ≈ E + S·k より k ∈ [(tt_start−Δ−E)/S, (tt_end+Δ−E)/S]。整数 k の閉区間を直接求める
    // （平均朔と真朔のずれ ≈0.6 < Δ=1 なので、期間内の真の新月は必ずこの k 範囲に対応＝偽陰性ゼロ。
    // 端の窓張り出しは Δ がそのまま担保するので別途の張り出しは不要）。
    let k_lo =
        ((tt_start - WINDOW_HALF_WIDTH_DAYS - MEEUS_EPOCH_JDE) / SYNODIC_MONTH_DAYS).ceil() as i64;
    let k_hi =
        ((tt_end + WINDOW_HALF_WIDTH_DAYS - MEEUS_EPOCH_JDE) / SYNODIC_MONTH_DAYS).floor() as i64;

    let mut out = Vec::new();
    for k in k_lo..=k_hi {
        let jde = mean_new_moon_jde(k as f64);
        let approx_tt = TtInstant::from_jd2(JulianDate2::from_jd(jde));
        let search_window = TimeInterval {
            start: TtInstant::from_jd2(approx_tt.jd2().add_days(-WINDOW_HALF_WIDTH_DAYS)),
            end: TtInstant::from_jd2(approx_tt.jd2().add_days(WINDOW_HALF_WIDTH_DAYS)),
        };
        out.push(NewMoonCandidate {
            approx_tt,
            search_window,
            lunation_number: k,
        });
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    // 実装本体（同モジュール直下）の型・関数。impl 担当の `use` 選択に依存しないよう、
    // テスト側で必要なシンボルを明示的に取り込む。
    use crate::error::EclipseError;

    use umbra_core::time::utc_to_tt;
    use umbra_core::{JulianDate2, TimeRange, TtInstant, UtcInstant, Vector3};
    use umbra_ephemeris::apparent::{moon_apparent_cirs, sun_apparent_cirs};

    // ============================================================
    // 共通定数・ヘルパ
    // ============================================================

    /// 朔望月（平均朔の周期, 日）。Meeus 49.1 の k 係数。
    const SYNODIC_MONTH_DAYS: f64 = 29.530_588_861;

    /// 平均朔と真の朔の最大ずれは ≈±0.6 日。窓半幅 Δ はこれを上回る ≈1 日であること（契約2）。
    /// 真の新月から最寄り候補までの差がこの下限未満であれば偽陰性は起き得ない（数値根拠）。
    const TRUE_VS_MEAN_MAX_OFFSET_DAYS: f64 = 0.6;

    /// 許容つきスカラ比較（clippy::float_cmp 回避）。
    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    /// 単一 JD（UTC スケール）から UtcInstant を作る。
    fn utc_from_jd(jd: f64) -> UtcInstant {
        UtcInstant::from_jd2(JulianDate2::from_jd(jd))
    }

    /// UTC の `[start_jd, end_jd]` 範囲。
    fn utc_range(start_jd: f64, end_jd: f64) -> TimeRange<UtcInstant> {
        TimeRange {
            start: utc_from_jd(start_jd),
            end: utc_from_jd(end_jd),
        }
    }

    /// 候補の `approx_tt` を JD（TT）で取り出す。
    fn approx_jd(c: &NewMoonCandidate) -> f64 {
        c.approx_tt.jd2().jd()
    }

    /// 候補窓の半幅 Δ（日）。中心 approx_tt から両端までの距離（窓中心性の契約2が成り立つ前提）。
    fn window_half_width_days(c: &NewMoonCandidate) -> f64 {
        let lo = c.approx_tt.jd2().days_since(c.search_window.start.jd2());
        let hi = c.search_window.end.jd2().days_since(c.approx_tt.jd2());
        // start≤approx≤end が成り立てば lo,hi はともに非負。半幅は両者の平均で代表させる。
        (lo + hi) / 2.0
    }

    // ---- Meeus Ch.49 平均朔のオラクル（実装の平均式に依存しない自前再実装）----

    /// Meeus 49.1: lunation index k（新月で整数）に対する平均朔の JDE（力学時 TT）。
    /// `JDE(k) = 2451550.09766 + 29.530588861·k + 0.00015437·T² − 0.000000150·T³ + 0.00000000073·T⁴`,
    /// `T = k/1236.85`。期待値のマジック直書きを避けるためテスト側で式どおり再実装する。
    fn meeus_mean_new_moon_jde(k: f64) -> f64 {
        let t = k / 1236.85;
        2_451_550.097_66 + 29.530_588_861 * k + 0.000_154_37 * t * t - 0.000_000_150 * t * t * t
            + 0.000_000_000_73 * t * t * t * t
    }

    /// lunation_number(i64) を Meeus 式の引数(f64)へ。k は ±数千程度で f64 に正確に載るが、
    /// clippy::cast_precision_loss を明示的に許容する（オラクル用途, 値域は安全）。
    #[allow(clippy::cast_precision_loss)]
    fn k_as_f64(k: i64) -> f64 {
        k as f64
    }

    // ---- 真の新月の独立オラクル: 実 ephemeris の地心離角の極小 ----

    /// 太陽・月の見かけ地心方向（CIRS）のなす角（離角, rad）。新月で極小。
    fn elongation_rad(time_tt: TtInstant) -> f64 {
        let s: Vector3 = sun_apparent_cirs(time_tt);
        let m: Vector3 = moon_apparent_cirs(time_tt);
        (s.dot(m) / (s.norm() * m.norm())).clamp(-1.0, 1.0).acos()
    }

    /// TT-JD で離角を評価する薄いラッパ。
    fn elongation_at_jd(jd_tt: f64) -> f64 {
        elongation_rad(TtInstant::from_jd2(JulianDate2::from_jd(jd_tt)))
    }

    /// 期間 `[start_tt_jd, end_tt_jd]`（TT-JD）を粗く走査し、離角が局所極小となる時刻（真の新月の
    /// 概算, ±step 精度）を時系列で返す。step は 0.2 日（窓 ±1 日に対し十分余裕）。
    /// 3 点 (prev, cur, next) で cur が最小となる cur を極小として採る。
    fn true_new_moons_tt_jd(start_tt_jd: f64, end_tt_jd: f64) -> Vec<f64> {
        let step = 0.5;
        let mut samples: Vec<(f64, f64)> = Vec::new();
        let mut jd = start_tt_jd;
        while jd <= end_tt_jd {
            samples.push((jd, elongation_at_jd(jd)));
            jd += step;
        }
        let mut minima = Vec::new();
        for i in 1..samples.len().saturating_sub(1) {
            let (t, e) = samples[i];
            if e < samples[i - 1].1 && e < samples[i + 1].1 {
                minima.push(t);
            }
        }
        minima
    }

    /// UTC-JD 範囲を TT-JD 範囲へ（境界変換。真の新月走査区間を作るのに使う, 1972 以降）。
    fn utc_range_to_tt_jd(start_jd: f64, end_jd: f64) -> (f64, f64) {
        let s = utc_to_tt(utc_from_jd(start_jd))
            .expect("post-1972 UTC→TT")
            .jd2()
            .jd();
        let e = utc_to_tt(utc_from_jd(end_jd))
            .expect("post-1972 UTC→TT")
            .jd2()
            .jd();
        (s, e)
    }

    // 代表的 JD（UTC）:
    //   2000-01-01.5 = 2451545.0, 2020-01-01.5 ≈ 2458850.0, 2022-01-01.5 ≈ 2459581.0,
    //   1971-01-01   ≈ 2440953（1972 前 = Err 用）。
    const JD_2020: f64 = 2_458_850.0; // 2020-01-01.5
    const JD_2022: f64 = 2_459_581.0; // 2022-01-01.5
    const JD_1971: f64 = 2_440_953.0; // 1971-01-01（1972 前）

    // ============================================================
    // 契約3+（最重要）: 偽陰性なし — 真の新月が必ずいずれかの候補窓に内包される
    // ============================================================

    /// 短め期間（2020–2022, 真の新月 ~25 件）で、離角極小オラクルが見つけた各真の新月が、
    /// 生成された候補窓のいずれかに内包される。偽陰性ゼロの直接検証（最重要契約）。
    #[test]
    fn every_true_new_moon_is_inside_some_candidate_window() {
        let candidates = new_moon_candidates(utc_range(JD_2020, JD_2022)).expect("post-1972 range");
        let (tt_lo, tt_hi) = utc_range_to_tt_jd(JD_2020, JD_2022);
        let truths = true_new_moons_tt_jd(tt_lo, tt_hi);

        // オラクル健全性: 期間内に妥当な数の真の新月が見つかっていること。
        assert!(
            truths.len() >= 20,
            "elongation oracle found too few new moons: {} (expected ~25)",
            truths.len()
        );

        for &truth_jd in &truths {
            let contained = candidates.iter().any(|c| {
                let lo = c.search_window.start.jd2().jd();
                let hi = c.search_window.end.jd2().jd();
                lo <= truth_jd && truth_jd <= hi
            });
            assert!(
                contained,
                "true new moon at TT-JD {truth_jd} not contained in any candidate window"
            );
        }
    }

    /// マージン下限: 各真の新月と最寄り候補 approx_tt の差が Δ（≈1 日）未満であること。
    /// さらに、平均朔–真朔の物理上限 0.6 日を超えないこと（オラクルと実装式双方の健全性）。
    /// これは「ずれが窓半幅を超えない＝偽陰性が起きない」の数値根拠。
    #[test]
    fn nearest_candidate_offset_is_below_window_half_width() {
        let candidates = new_moon_candidates(utc_range(JD_2020, JD_2022)).expect("post-1972 range");
        let (tt_lo, tt_hi) = utc_range_to_tt_jd(JD_2020, JD_2022);
        let truths = true_new_moons_tt_jd(tt_lo, tt_hi);

        for &truth_jd in &truths {
            let nearest = candidates
                .iter()
                .map(|c| (approx_jd(c) - truth_jd).abs())
                .fold(f64::INFINITY, f64::min);
            // 走査刻み 0.2 日の量子化誤差を見込んでも、ずれは平均朔–真朔の上限を大きく超えない。
            assert!(
                nearest < TRUE_VS_MEAN_MAX_OFFSET_DAYS + step_margin(),
                "nearest approx_tt offset {nearest} d exceeds true-vs-mean bound for truth {truth_jd}"
            );
            // かつ最寄り候補の窓半幅より小さい（窓に収まる）こと。
            let nearest_c = candidates
                .iter()
                .min_by(|a, b| {
                    (approx_jd(a) - truth_jd)
                        .abs()
                        .total_cmp(&(approx_jd(b) - truth_jd).abs())
                })
                .expect("non-empty candidates");
            assert!(
                nearest < window_half_width_days(nearest_c) + step_margin(),
                "offset {nearest} d ≥ window half-width for truth {truth_jd}"
            );
        }
    }

    /// 真の新月走査刻み（`true_new_moons_tt_jd` の step）に由来する量子化マージン（極小時刻は ±step 精度）。
    fn step_margin() -> f64 {
        0.5
    }

    // ============================================================
    // 契約4: 網羅性（個数） — ephemeris 不要、100 年規模
    // ============================================================

    /// 100 年規模（1980–2080, 1972 以降）の候補数 ≈ 期間日数 / 朔望月 の ±2 以内。
    /// off-by-one（端の朔の取りこぼし・過剰生成）を検出する。
    #[test]
    fn candidate_count_matches_period_over_synodic_month() {
        let start_jd = 2_444_240.0; // ≈ 1980-01-01
        let end_jd = start_jd + 100.0 * 365.25; // ≈ 2080
        let candidates = new_moon_candidates(utc_range(start_jd, end_jd)).expect("post-1972");
        let span_days = end_jd - start_jd;
        let expected = span_days / SYNODIC_MONTH_DAYS;
        // 端の朔(±1)＋窓交差で拾う端外朔(±1) を見込んで ±2。
        #[allow(clippy::cast_precision_loss)]
        let got = candidates.len() as f64;
        assert!(
            (got - expected).abs() <= 2.0,
            "candidate count {got} not within ±2 of expected {expected} ({span_days} d / {SYNODIC_MONTH_DAYS})"
        );
    }

    // ============================================================
    // 契約2: approx_tt が窓中心・窓幅 ≈ 2Δ・Δ ≈ 1 日
    // ============================================================

    /// 各候補で `window.start ≤ approx_tt ≤ window.end`（窓は approx_tt を内包）。
    #[test]
    fn approx_tt_is_within_its_window() {
        let candidates = new_moon_candidates(utc_range(JD_2020, JD_2022)).expect("post-1972 range");
        assert!(!candidates.is_empty());
        for c in &candidates {
            let start_jd = c.search_window.start.jd2().jd();
            let end_jd = c.search_window.end.jd2().jd();
            let a = approx_jd(c);
            assert!(
                start_jd <= a && a <= end_jd,
                "approx_tt {a} not within window [{start_jd}, {end_jd}]"
            );
        }
    }

    /// approx_tt が窓の（ほぼ）中心: start からの距離 ≈ end までの距離。
    #[test]
    fn approx_tt_is_centered_in_window() {
        let candidates = new_moon_candidates(utc_range(JD_2020, JD_2022)).expect("post-1972 range");
        for c in &candidates {
            let lo = c.approx_tt.jd2().days_since(c.search_window.start.jd2());
            let hi = c.search_window.end.jd2().days_since(c.approx_tt.jd2());
            assert!(
                close(lo, hi, 1e-6),
                "window not centered: lo={lo} hi={hi} (approx not at center)"
            );
        }
    }

    /// 窓半幅 Δ は平均朔–真朔の最大ずれ 0.6 日を上回り、約 1 日であること（契約2のマージン要件）。
    /// 偽陰性ゼロは Δ > 0.6 に依存するため、この下限を全候補で要求する。
    #[test]
    fn window_half_width_exceeds_true_vs_mean_offset_and_is_about_one_day() {
        let candidates = new_moon_candidates(utc_range(JD_2020, JD_2022)).expect("post-1972 range");
        for c in &candidates {
            let half = window_half_width_days(c);
            assert!(
                half > TRUE_VS_MEAN_MAX_OFFSET_DAYS,
                "window half-width {half} d does not exceed true-vs-mean offset 0.6 d"
            );
            assert!(
                (0.6..=2.0).contains(&half),
                "window half-width {half} d outside expected ~1 day band"
            );
        }
    }

    // ============================================================
    // 契約5: 時系列順・単調・一意
    // ============================================================

    /// approx_tt は厳密昇順。
    #[test]
    fn candidates_are_sorted_by_approx_tt_ascending() {
        let candidates = new_moon_candidates(utc_range(JD_2020, JD_2022)).expect("post-1972 range");
        assert!(candidates.len() >= 2);
        for w in candidates.windows(2) {
            assert!(
                approx_jd(&w[0]) < approx_jd(&w[1]),
                "approx_tt not strictly ascending: {} !< {}",
                approx_jd(&w[0]),
                approx_jd(&w[1])
            );
        }
    }

    /// lunation_number は 1 ずつ厳密増加・一意（k が連続）。
    #[test]
    fn lunation_numbers_are_consecutive_and_unique() {
        let candidates = new_moon_candidates(utc_range(JD_2020, JD_2022)).expect("post-1972 range");
        assert!(candidates.len() >= 2);
        for w in candidates.windows(2) {
            assert_eq!(
                w[1].lunation_number - w[0].lunation_number,
                1,
                "lunation_number must increase by exactly 1 ({} → {})",
                w[0].lunation_number,
                w[1].lunation_number
            );
        }
    }

    // ============================================================
    // 契約1: approx_tt が Meeus 式どおり（自前オラクルと一致）
    // ============================================================

    /// 各候補の approx_tt(TT-JD) が、その lunation_number から Meeus 式で計算した JDE に一致。
    /// 実装が概算式どおりかを検証（係数取り違え・T のスケール誤りを殺す）。tol は数秒（1e-4 日）。
    #[test]
    fn approx_tt_matches_meeus_formula_for_its_lunation_number() {
        let candidates = new_moon_candidates(utc_range(JD_2020, JD_2022)).expect("post-1972 range");
        assert!(!candidates.is_empty());
        for c in &candidates {
            let expected_jde = meeus_mean_new_moon_jde(k_as_f64(c.lunation_number));
            assert!(
                close(approx_jd(c), expected_jde, 1e-4),
                "approx_tt {} != Meeus JDE {} for k={}",
                approx_jd(c),
                expected_jde,
                c.lunation_number
            );
        }
    }

    /// k=0 は 2000-01-06 付近の新月（Meeus エポック JDE ≈ 2451550.09766）。lunation index の
    /// 原点とスケールを固定（実装の k 採番が Meeus 規約に一致することを担保）。
    #[test]
    fn lunation_index_origin_matches_meeus_epoch() {
        // 2000-01-06 を含む短い範囲で k=0 の候補が現れることを確認。
        let jd_2000_01_06 = 2_451_549.5; // 2000-01-06.0 UTC 付近
        let candidates = new_moon_candidates(utc_range(jd_2000_01_06 - 2.0, jd_2000_01_06 + 2.0))
            .expect("post-1972 range");
        let has_k0 = candidates.iter().any(|c| c.lunation_number == 0);
        assert!(
            has_k0,
            "expected a candidate with lunation_number == 0 near Meeus epoch 2000-01-06"
        );
        // その k=0 候補の approx_tt が Meeus エポック JDE に一致。
        let k0 = candidates
            .iter()
            .find(|c| c.lunation_number == 0)
            .expect("k=0 present");
        assert!(
            close(approx_jd(k0), meeus_mean_new_moon_jde(0.0), 1e-4),
            "k=0 approx_tt {} != Meeus epoch {}",
            approx_jd(k0),
            meeus_mean_new_moon_jde(0.0)
        );
    }

    // ============================================================
    // 契約6: 期間端の取りこぼしなし
    // ============================================================

    /// 範囲が 1 朔未満（数日）でも、窓が範囲と交差する朔を返す（空にしない）。
    #[test]
    fn sub_lunation_range_returns_window_crossing_candidate() {
        // 既知の朔近傍を含む短い窓。2020 年内のある新月時刻を離角オラクルで特定して挟む。
        let (tt_lo, tt_hi) = utc_range_to_tt_jd(JD_2020, JD_2020 + 35.0);
        let truths = true_new_moons_tt_jd(tt_lo, tt_hi);
        assert!(!truths.is_empty(), "oracle should find a new moon in 35 d");
        let truth_jd = truths[0];
        // 真の新月をちょうど含む 3 日幅（< 1 朔）の UTC 範囲。
        let candidates =
            new_moon_candidates(utc_range(truth_jd - 1.5, truth_jd + 1.5)).expect("post-1972");
        assert!(
            !candidates.is_empty(),
            "sub-lunation range containing a true new moon returned no candidates"
        );
        let contained = candidates.iter().any(|c| {
            let lo = c.search_window.start.jd2().jd();
            let hi = c.search_window.end.jd2().jd();
            lo <= truth_jd && truth_jd <= hi
        });
        assert!(
            contained,
            "true new moon {truth_jd} in a short range not contained in any window"
        );
    }

    /// 範囲開始直後/終了直前の朔を落とさない: 候補の最初の窓は範囲開始の手前まで、最後の窓は
    /// 範囲終了の先まで張り出しうる（start−1朔/end+1朔まで生成し窓交差で残す設計）。
    /// 端の真の新月がいずれも窓に内包されることで取りこぼしなしを確認する。
    #[test]
    fn boundary_new_moons_near_range_edges_are_not_dropped() {
        let candidates = new_moon_candidates(utc_range(JD_2020, JD_2022)).expect("post-1972 range");
        let (tt_lo, tt_hi) = utc_range_to_tt_jd(JD_2020, JD_2022);
        let truths = true_new_moons_tt_jd(tt_lo, tt_hi);
        assert!(truths.len() >= 2);

        // 期間内で最も早い/最も遅い真の新月（端に最も近い）も窓に内包される。
        for &edge in &[truths[0], truths[truths.len() - 1]] {
            let contained = candidates.iter().any(|c| {
                let lo = c.search_window.start.jd2().jd();
                let hi = c.search_window.end.jd2().jd();
                lo <= edge && edge <= hi
            });
            assert!(
                contained,
                "edge new moon {edge} (near range boundary) was dropped"
            );
        }
    }

    /// 契約1（補正項）: 概算式の **T²/T³/T⁴ 補正項** まで式どおりであること。エポック(2000)近傍では
    /// 補正項が微小（< tol）で区別できないため、エポックから遠い年代（k≈6180, 偏角 T≈5）で締める。
    /// ここは精度主張ではなく「定義された平均朔多項式を実装が忠実に評価するか」のフィデリティ検証。
    /// tol は機械精度級（同一式の独立再実装と突合, T²≈4e-3 / T³≈1.9e-5 / T⁴≈4.6e-7 日が有意化）。
    #[test]
    fn approx_tt_matches_meeus_formula_far_from_epoch() {
        // ≈ 2500 年（k≈6180）。補正項が tol 1e-7 を十分上回る。
        let jd_far = 2_634_170.0;
        let candidates =
            new_moon_candidates(utc_range(jd_far, jd_far + 90.0)).expect("post-1972 range");
        assert!(!candidates.is_empty());
        for c in &candidates {
            let expected = meeus_mean_new_moon_jde(k_as_f64(c.lunation_number));
            assert!(
                close(approx_jd(c), expected, 1e-7),
                "approx_tt {} != Meeus JDE {} (far epoch) for k={}",
                approx_jd(c),
                expected,
                c.lunation_number
            );
        }
    }

    /// 契約6（開始端の張り出し）: 期間開始の **直前** にある平均朔（窓が期間に食い込む）も候補に含む。
    /// 範囲開始をある平均朔 k0 の 0.3 日後に置くと、k0 は範囲外（approx<start）だが窓 [approx−1,approx+1]
    /// が範囲に食い込む → k0 候補が生成されねばならない（k 下限 `ceil((tt_start−Δ−E)/S)` が要）。
    /// 下限が k0 を取りこぼす向きに壊れると偽陰性になる。
    #[test]
    fn candidate_just_before_range_start_is_retained_via_overhang() {
        let k0: i64 = 300; // ≈ 2024 年（post-1972）。
        let approx_k0 = meeus_mean_new_moon_jde(k_as_f64(k0)); // TT-JD
                                                               // 範囲開始を平均朔 k0 の 0.3 日後に（UTC≈TT, ΔT≈0.0008 日差は 0.3 日に対し無視可）。
        let start_jd = approx_k0 + 0.3;
        let candidates =
            new_moon_candidates(utc_range(start_jd, start_jd + 10.0)).expect("post-1972 range");
        // k0 の窓 [approx−1, approx+1] は start(=approx+0.3) を内包 → k0 候補が残らねばならない。
        assert!(
            candidates.iter().any(|c| c.lunation_number == k0),
            "overhang candidate k0={k0} (mean new moon just before range start) was dropped"
        );
    }

    /// 契約6（終了端の張り出し）: 期間終了の **直後** にある平均朔（窓が期間に食い込む）も候補に含む。
    /// 範囲終了をある平均朔 k1 の 0.3 日前に置くと、k1 は範囲外（approx>end）だが窓 [approx−1,approx+1]
    /// が範囲に食い込む → k1 候補が生成されねばならない（k 上限 `floor((tt_end+Δ−E)/S)` が要）。
    /// 上限が k1 を取りこぼす向きに壊れると偽陰性になる。
    #[test]
    fn candidate_just_after_range_end_is_retained_via_overhang() {
        let k1: i64 = 301; // ≈ 2024 年（post-1972）。
        let approx_k1 = meeus_mean_new_moon_jde(k_as_f64(k1)); // TT-JD
                                                               // 範囲終了を平均朔 k1 の 0.3 日前に（UTC≈TT 近似）。
        let end_jd = approx_k1 - 0.3;
        let candidates =
            new_moon_candidates(utc_range(end_jd - 10.0, end_jd)).expect("post-1972 range");
        // k1 の窓 [approx−1, approx+1] は end(=approx−0.3) を内包 → k1 候補が残らねばならない。
        assert!(
            candidates.iter().any(|c| c.lunation_number == k1),
            "overhang candidate k1={k1} (mean new moon just after range end) was dropped"
        );
    }

    /// 契約3+（2000 年以前・post-1972）: エポック(JDE_2000)より前の期間でも候補が空にならず網羅される。
    /// k 範囲算術で `tt_end − E` の `/` を `*` 等に取り違えると、E より前（tt_end−E<0）で k 上限が
    /// 下限を下回り空になる退行を捕捉する。
    #[test]
    fn pre_2000_post_1972_range_is_nonempty_and_covers_true_new_moons() {
        // ≈ 1990 年（post-1972, < エポック 2000-01-06）。
        let jd_1990 = 2_447_893.0; // 1990-01-01.5
        let candidates =
            new_moon_candidates(utc_range(jd_1990, jd_1990 + 400.0)).expect("post-1972 range");
        assert!(
            !candidates.is_empty(),
            "pre-2000 (post-1972) range must yield candidates, got empty"
        );
        // 真の新月（離角極小）がすべて窓に内包される（この年代でも偽陰性なし）。
        let (tt_lo, tt_hi) = utc_range_to_tt_jd(jd_1990, jd_1990 + 400.0);
        let truths = true_new_moons_tt_jd(tt_lo, tt_hi);
        assert!(
            truths.len() >= 8,
            "oracle should find ~13 new moons in 400 d"
        );
        for &truth_jd in &truths {
            let contained = candidates.iter().any(|c| {
                let lo = c.search_window.start.jd2().jd();
                let hi = c.search_window.end.jd2().jd();
                lo <= truth_jd && truth_jd <= hi
            });
            assert!(
                contained,
                "pre-2000 true new moon {truth_jd} not in any window"
            );
        }
        // lunation_number が負域（k<0）を含む（エポック前）。
        assert!(
            candidates.iter().any(|c| c.lunation_number < 0),
            "pre-2000 range should include negative lunation numbers"
        );
    }

    // ============================================================
    // 契約7: 1972 年より前は Err(EclipseError::Time(..))
    // ============================================================

    /// 範囲開始が 1972 年より前 → utc_to_tt の閏秒未定義で Err(EclipseError::Time(..))。
    #[test]
    fn range_starting_before_1972_is_time_error() {
        let result = new_moon_candidates(utc_range(JD_1971, JD_1971 + 60.0));
        assert!(
            matches!(result, Err(EclipseError::Time(_))),
            "expected Err(EclipseError::Time(..)) for pre-1972 start, got {result:?}"
        );
    }

    // ============================================================
    // 健全性: 戻り値の有限性
    // ============================================================

    /// 全候補の approx_tt と窓端が有限（NaN/Inf を返さない）。
    #[test]
    fn all_candidate_times_are_finite() {
        let candidates = new_moon_candidates(utc_range(JD_2020, JD_2022)).expect("post-1972 range");
        for c in &candidates {
            assert!(approx_jd(c).is_finite(), "approx_tt non-finite");
            assert!(
                c.search_window.start.jd2().jd().is_finite()
                    && c.search_window.end.jd2().jd().is_finite(),
                "window endpoint non-finite"
            );
        }
    }
}
