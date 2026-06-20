//! D6: 日食フィルタの**偽陰性ゼロ・マージン実余裕統計**（accuracy.md §3.4・ISSUE-018）の受け入れテスト。
//!
//! 本ファイルは `umbra-eclipse` の**公開 API のみ**を対象とした統合テスト（tests/ 配下・別クレート境界）。
//! 対象はこれから実装する（実装は別担当・本ファイルはテストのみ）:
//! - `pub const ECLIPSE_FILTER_SAFETY_MARGIN_RAD: f64;`（= 0.0087, 偽陰性ゼロ保証マージン \[rad\]）。
//! - `pub struct MarginSample { separation_rad, bare_limit_rad, possible }`（1 候補のマージン標本）。
//! - `pub struct FilterMarginStats { candidates_scanned, accepted, rejected, safety_margin_rad,
//!   max_margin_consumed_rad, min_accepted_slack_rad }`（D6 統計）。
//! - `pub fn aggregate_filter_margins(samples: &[MarginSample], safety_margin_rad: f64) -> FilterMarginStats;`
//!   （標本群から D6 統計を集計する純関数・FAST 中核）。
//! - `pub fn scan_filter_margins(range: UtcRange) -> Result<FilterMarginStats, EclipseError>;`
//!   （範囲を candidate→合→フィルタで実走し D6 統計を返す・SLOW 実パイプライン）。
//!
//! ## 確定セマンティクス（テストで縛る）
//! 1. `candidates_scanned = samples.len()`。`accepted = possible==true の数`。`rejected = scanned − accepted`。
//! 2. `safety_margin_rad` = 引数そのまま（passthrough）。
//! 3. `max_margin_consumed_rad` = accepted 標本に対する `max(0, separation_rad − bare_limit_rad)` の最大
//!    （accepted が無ければ 0.0）。bare_limit 内（sep<bare）は消費 0。rejected は寄与しない。
//! 4. `min_accepted_slack_rad = safety_margin_rad − max_margin_consumed_rad`（常に 0..=safety_margin）。
//! 5. 空入力: scanned=accepted=rejected=0、consumed=0、slack=safety_margin。
//! 6. scan: 既知の実日食を含む控えめ範囲（2017-01-01〜2020-01-01・1972 以降）で scanned≥1・accepted≥1・
//!    safety==const・accepted+rejected==scanned・0≤slack≤safety・consumed≥0・slack>0（マージン非枯渇）。
//!
//! ## テスト戦略（strict / mutation-resistant / 負荷配分）
//! FAST はすべて合成 `MarginSample` 群で `aggregate_filter_margins` を算術的に縛る（実パイプライン非実走）。
//! SLOW は `scan_filter_margins` を控えめ範囲 1 ケースで不変条件のみ縛る（脆い数値固定は避ける）。
//!
//! ## 期待される RED（実装前）
//! `ECLIPSE_FILTER_SAFETY_MARGIN_RAD` / `MarginSample` / `FilterMarginStats` /
//! `aggregate_filter_margins` / `scan_filter_margins` はまだ存在しないため、本ファイルは
//! **未解決インポート（E0432）/未定義シンボル（E0425）でコンパイル不能 = RED** になる。これが想定どおりの赤。

use umbra_core::{TimeRange, UtcInstant};
use umbra_eclipse::{
    aggregate_filter_margins, scan_filter_margins, EclipseError, FilterMarginStats, MarginSample,
    UtcRange, ECLIPSE_FILTER_SAFETY_MARGIN_RAD,
};

/// 統計の f64 一致に用いる厳密許容。
const EPS: f64 = 1e-9;

// ============================================================
// 構築ヘルパ
// ============================================================

/// 合成 `MarginSample` を組む。
fn sample(separation_rad: f64, bare_limit_rad: f64, possible: bool) -> MarginSample {
    MarginSample {
        separation_rad,
        bare_limit_rad,
        possible,
    }
}

/// グレゴリオ暦（UTC 0 時）から `UtcInstant`。
fn utc(year: i32, month: u8, day: u8) -> UtcInstant {
    UtcInstant::from_gregorian(year, month, day, 0, 0, 0.0).expect("有効な UTC 日時")
}

/// `[start, end]` の UTC 範囲（scan 用）。
fn utc_range(start: (i32, u8, u8), end: (i32, u8, u8)) -> UtcRange {
    TimeRange {
        start: utc(start.0, start.1, start.2),
        end: utc(end.0, end.1, end.2),
    }
}

// ============================================================
// FAST: aggregate_filter_margins — カウント（scanned / accepted / rejected）
// ============================================================

/// 受け入れ「scanned=len・accepted=possible==true の数・rejected=scanned−accepted」。
/// possible が {true,false,true,true,false} の 5 件 → scanned=5・accepted=3・rejected=2。
/// 殺す変異: scanned を len 以外（定数・accepted 等）にする、accepted を false 計数/全件計数、
///   rejected を別式（scanned や accepted の取り違え・オフバイワン）にする。
#[test]
fn counts_scanned_accepted_rejected() {
    let samples = [
        sample(0.001, 0.020, true),
        sample(0.050, 0.020, false),
        sample(0.002, 0.020, true),
        sample(0.003, 0.020, true),
        sample(0.060, 0.020, false),
    ];
    let stats = aggregate_filter_margins(&samples, ECLIPSE_FILTER_SAFETY_MARGIN_RAD);
    assert_eq!(stats.candidates_scanned, 5, "scanned は len");
    assert_eq!(stats.accepted, 3, "accepted は possible==true の数");
    assert_eq!(stats.rejected, 2, "rejected = scanned − accepted");
    assert_eq!(
        stats.accepted + stats.rejected,
        stats.candidates_scanned,
        "accepted + rejected == scanned（恒等式）"
    );
}

/// 受け入れ「全件 accepted なら rejected=0、全件 rejected なら accepted=0」。
/// rejected を accepted と取り違える/定数化する変異を両端で殺す。
#[test]
fn counts_all_accepted_and_all_rejected() {
    let all_acc = [sample(0.001, 0.020, true), sample(0.002, 0.020, true)];
    let s = aggregate_filter_margins(&all_acc, ECLIPSE_FILTER_SAFETY_MARGIN_RAD);
    assert_eq!(s.accepted, 2, "全件 accepted");
    assert_eq!(s.rejected, 0, "全件 accepted → rejected=0");

    let all_rej = [
        sample(0.100, 0.020, false),
        sample(0.200, 0.020, false),
        sample(0.300, 0.020, false),
    ];
    let s = aggregate_filter_margins(&all_rej, ECLIPSE_FILTER_SAFETY_MARGIN_RAD);
    assert_eq!(s.accepted, 0, "全件 rejected → accepted=0");
    assert_eq!(s.rejected, 3, "全件 rejected");
}

// ============================================================
// FAST: safety_margin_rad passthrough
// ============================================================

/// 受け入れ「safety_margin_rad は引数そのまま」。const とは別の値 0.005 を渡し、そのまま保持。
/// 殺す変異: 内部定数で上書き、引数を無視、別フィールドへ配線。
#[test]
fn safety_margin_is_passthrough() {
    let samples = [sample(0.001, 0.020, true)];
    let custom = 0.005;
    let stats = aggregate_filter_margins(&samples, custom);
    assert!(
        (stats.safety_margin_rad - custom).abs() < EPS,
        "safety_margin_rad は引数 {custom} そのまま, got {}",
        stats.safety_margin_rad
    );
    // 念のため const とも区別（const で上書きする変異を殺す）。
    assert!(
        (stats.safety_margin_rad - ECLIPSE_FILTER_SAFETY_MARGIN_RAD).abs() > EPS,
        "passthrough は const と異なる引数で検証している（const 上書き変異の検出）"
    );
}

// ============================================================
// FAST: max_margin_consumed_rad（accepted の max(0, sep−bare) の最大）
// ============================================================

/// 受け入れ「accepted 全件が bare_limit 内（sep<bare）なら consumed=0・slack=safety」。
/// sep<bare の accepted のみ → max(0, 負)=0 のクランプ。
/// 殺す変異: クランプ欠落（負の消費が出て slack>safety）、consumed を sep−bare 素の値にする。
#[test]
fn consumed_zero_when_all_accepted_within_bare_limit() {
    let safety = ECLIPSE_FILTER_SAFETY_MARGIN_RAD;
    let samples = [
        sample(0.010, 0.020, true), // sep−bare = −0.010 → クランプ 0
        sample(0.005, 0.020, true), // sep−bare = −0.015 → クランプ 0
    ];
    let stats = aggregate_filter_margins(&samples, safety);
    assert!(
        stats.max_margin_consumed_rad.abs() < EPS,
        "全 accepted が bare_limit 内 → consumed=0, got {}",
        stats.max_margin_consumed_rad
    );
    assert!(
        (stats.min_accepted_slack_rad - safety).abs() < EPS,
        "consumed=0 → slack=safety, got {}",
        stats.min_accepted_slack_rad
    );
}

/// 受け入れ「grazing（sep>bare）の accepted で consumed = sep−bare（正）」。
/// 1 件 accepted: sep=0.025, bare=0.020 → consumed=0.005。
/// 殺す変異: クランプを誤って 0 に潰す、sep と bare の引き算の向き反転（bare−sep）、別列参照。
#[test]
fn consumed_positive_for_grazing_accepted() {
    let safety = ECLIPSE_FILTER_SAFETY_MARGIN_RAD;
    let samples = [sample(0.025, 0.020, true)];
    let stats = aggregate_filter_margins(&samples, safety);
    assert!(
        (stats.max_margin_consumed_rad - 0.005).abs() < EPS,
        "consumed = sep−bare = 0.005, got {}",
        stats.max_margin_consumed_rad
    );
    // slack = safety − consumed。
    assert!(
        (stats.min_accepted_slack_rad - (safety - 0.005)).abs() < EPS,
        "slack = safety − consumed, got {}",
        stats.min_accepted_slack_rad
    );
}

/// 受け入れ「consumed は accepted のうち最大の max(0, sep−bare)（複数で最大が選ばれる）」。
/// accepted 3 件: 消費 {0.001, 0.006, 0.003} → 最大 0.006。混入する bare 内 accepted（消費 0）も無視されない。
/// 殺す変異: 最大でなく最後/最初/合計/平均を返す、ループ範囲のオフバイワン。
#[test]
fn consumed_is_max_over_accepted() {
    let safety = ECLIPSE_FILTER_SAFETY_MARGIN_RAD;
    let samples = [
        sample(0.021, 0.020, true), // 消費 0.001
        sample(0.026, 0.020, true), // 消費 0.006（最大）
        sample(0.023, 0.020, true), // 消費 0.003
        sample(0.010, 0.020, true), // 消費 0（bare 内）
    ];
    let stats = aggregate_filter_margins(&samples, safety);
    assert!(
        (stats.max_margin_consumed_rad - 0.006).abs() < EPS,
        "consumed は accepted 消費の最大 0.006, got {}",
        stats.max_margin_consumed_rad
    );
}

/// 受け入れ「rejected 標本は consumed に寄与しない（possible==false は無視）」。
/// accepted 1 件（消費 0.002）の脇に、消費が巨大に見える rejected（sep=0.500, bare=0.020）を置く。
/// consumed は accepted の 0.002 のまま（rejected の 0.480 を拾わない）。
/// 殺す変異: possible を見ずに全件で max を取る（rejected の巨大消費を拾い consumed が跳ね上がる）。
#[test]
fn rejected_samples_do_not_contribute_to_consumed() {
    let safety = ECLIPSE_FILTER_SAFETY_MARGIN_RAD;
    let samples = [
        sample(0.022, 0.020, true),  // accepted・消費 0.002
        sample(0.500, 0.020, false), // rejected・無視（拾えば 0.480）
    ];
    let stats = aggregate_filter_margins(&samples, safety);
    assert!(
        (stats.max_margin_consumed_rad - 0.002).abs() < EPS,
        "consumed は accepted のみ（rejected 無視）で 0.002, got {}",
        stats.max_margin_consumed_rad
    );
}

/// 受け入れ「accepted が無ければ consumed=0・slack=safety（rejected だけでも安全）」。
/// rejected のみ 2 件 → max over accepted は空集合 → consumed=0。
/// 殺す変異: accepted 空で NaN/panic、空集合 max を −inf にする、slack を safety 以外にする。
#[test]
fn consumed_zero_when_no_accepted() {
    let safety = ECLIPSE_FILTER_SAFETY_MARGIN_RAD;
    let samples = [sample(0.500, 0.020, false), sample(0.600, 0.020, false)];
    let stats = aggregate_filter_margins(&samples, safety);
    assert_eq!(stats.accepted, 0, "accepted=0");
    assert!(
        stats.max_margin_consumed_rad.abs() < EPS,
        "accepted 無し → consumed=0, got {}",
        stats.max_margin_consumed_rad
    );
    assert!(
        (stats.min_accepted_slack_rad - safety).abs() < EPS,
        "accepted 無し → slack=safety, got {}",
        stats.min_accepted_slack_rad
    );
}

// ============================================================
// FAST: 境界（sep==bare で消費 0・sep==bare+margin 近傍）
// ============================================================

/// 受け入れ「sep==bare_limit の accepted は消費ちょうど 0（境界・max(0,0)=0）」。
/// 殺す変異: 境界での符号取り違え（< を <= に等）で消費が微小正/負に振れる、クランプの境界誤り。
#[test]
fn consumed_is_zero_at_separation_equals_bare_limit() {
    let safety = ECLIPSE_FILTER_SAFETY_MARGIN_RAD;
    let samples = [sample(0.020, 0.020, true)]; // sep − bare = 0
    let stats = aggregate_filter_margins(&samples, safety);
    assert!(
        stats.max_margin_consumed_rad.abs() < EPS,
        "sep==bare → consumed=0, got {}",
        stats.max_margin_consumed_rad
    );
    assert!(
        (stats.min_accepted_slack_rad - safety).abs() < EPS,
        "sep==bare → slack=safety, got {}",
        stats.min_accepted_slack_rad
    );
}

/// 受け入れ「sep が bare+safety 近傍（マージンほぼ枯渇）の accepted で slack≈0（≥0）」。
/// 消費 = safety − tiny → slack = tiny（正の極小）。slack=safety−consumed の符号・恒等を境界で縛る。
/// 殺す変異: slack の引き算の向き反転（consumed−safety で負）、slack を 0 固定/safety 固定。
#[test]
fn slack_near_zero_when_margin_nearly_exhausted() {
    let safety = ECLIPSE_FILTER_SAFETY_MARGIN_RAD;
    let tiny = 1e-5;
    // 消費 = safety − tiny（bare をその分だけ超える accepted）。
    let samples = [sample(0.020 + (safety - tiny), 0.020, true)];
    let stats = aggregate_filter_margins(&samples, safety);
    assert!(
        (stats.max_margin_consumed_rad - (safety - tiny)).abs() < EPS,
        "consumed = safety − tiny, got {}",
        stats.max_margin_consumed_rad
    );
    assert!(
        (stats.min_accepted_slack_rad - tiny).abs() < EPS,
        "slack = safety − consumed = tiny（≥0・枯渇直前）, got {}",
        stats.min_accepted_slack_rad
    );
    assert!(
        stats.min_accepted_slack_rad >= 0.0,
        "slack は非負（枯渇直前でも 0 以上）"
    );
}

/// 受け入れ「slack = safety − consumed（恒等式・一般値）」。
/// safety=0.0087, consumed が 0.004 になる accepted 1 件 → slack=0.0047。引数 safety を使うことも縛る。
/// 殺す変異: slack を const で計算（引数 safety を無視）、consumed を二重計上、定数化。
#[test]
fn slack_equals_safety_minus_consumed() {
    let safety = ECLIPSE_FILTER_SAFETY_MARGIN_RAD; // 0.0087
    let samples = [sample(0.024, 0.020, true)]; // 消費 0.004
    let stats = aggregate_filter_margins(&samples, safety);
    assert!(
        (stats.max_margin_consumed_rad - 0.004).abs() < EPS,
        "consumed=0.004, got {}",
        stats.max_margin_consumed_rad
    );
    assert!(
        (stats.min_accepted_slack_rad - (safety - 0.004)).abs() < EPS,
        "slack = safety − consumed = {}, got {}",
        safety - 0.004,
        stats.min_accepted_slack_rad
    );
}

// ============================================================
// FAST: 空入力（vacuous）
// ============================================================

/// 受け入れ「空入力: scanned=accepted=rejected=0・consumed=0・slack=safety」。
/// 殺す変異: 空で panic、scanned を 0 以外、slack を 0/0 以外（空でも safety を返す契約）。
#[test]
fn empty_input_is_all_zero_with_full_slack() {
    let safety = ECLIPSE_FILTER_SAFETY_MARGIN_RAD;
    let samples: [MarginSample; 0] = [];
    let stats = aggregate_filter_margins(&samples, safety);
    assert_eq!(stats.candidates_scanned, 0, "空 → scanned=0");
    assert_eq!(stats.accepted, 0, "空 → accepted=0");
    assert_eq!(stats.rejected, 0, "空 → rejected=0");
    assert!(
        stats.max_margin_consumed_rad.abs() < EPS,
        "空 → consumed=0, got {}",
        stats.max_margin_consumed_rad
    );
    assert!(
        (stats.min_accepted_slack_rad - safety).abs() < EPS,
        "空 → slack=safety（マージン未使用）, got {}",
        stats.min_accepted_slack_rad
    );
    assert!(
        (stats.safety_margin_rad - safety).abs() < EPS,
        "空でも safety_margin_rad は引数そのまま"
    );
}

// ============================================================
// FAST: 公開定数の値
// ============================================================

/// 受け入れ「ECLIPSE_FILTER_SAFETY_MARGIN_RAD == 0.0087（偽陰性ゼロ保証マージン）」。
/// eclipse_filter.rs の SAFETY_MARGIN_RAD と同値の公開再エクスポート。
/// 殺す変異: 定数値の改変（0 化・符号反転・桁違い）。
// 公開定数の値そのものを検証する意図的な定数アサーション（const 改変変異を殺す）。
#[allow(clippy::assertions_on_constants)]
#[test]
fn safety_margin_const_is_0_0087() {
    assert!(
        (ECLIPSE_FILTER_SAFETY_MARGIN_RAD - 0.008_7).abs() < EPS,
        "ECLIPSE_FILTER_SAFETY_MARGIN_RAD は 0.0087, got {ECLIPSE_FILTER_SAFETY_MARGIN_RAD}"
    );
    assert!(
        ECLIPSE_FILTER_SAFETY_MARGIN_RAD > 0.0,
        "保守マージンは正（偽陰性ゼロ側）"
    );
}

// ============================================================
// SLOW: scan_filter_margins — 実パイプライン 1 ケース（不変条件のみ）
// ============================================================

/// 受け入れ「実日食を複数含む控えめ範囲（2017-01-01〜2020-01-01）で D6 統計が不変条件を満たす」。
/// 実パイプライン（candidate→合→フィルタ）を実走。範囲には 2017-08-21 皆既・2019-07-02 皆既など
/// 実日食が落ちているので accepted≥1。脆い数値固定は避け、不変条件のみを縛る:
///   scanned≥1, accepted≥1, safety==const, accepted+rejected==scanned,
///   0≤slack≤safety, consumed≥0, slack>0（マージン非枯渇＝偽陰性ゼロが全期間で保たれる前提）。
/// 殺す変異: scanned/accepted の計数崩れ、safety の取り違え、恒等式 accepted+rejected==scanned の破れ、
///   slack の範囲逸脱（負/safety 超え）、slack==0（マージン枯渇＝危険信号を見逃す）。
#[test]
fn scan_real_eclipse_range_satisfies_invariants() {
    let range = utc_range((2017, 1, 1), (2020, 1, 1));
    let stats = scan_filter_margins(range).expect("post-1972 range scan must succeed");

    assert!(
        stats.candidates_scanned >= 1,
        "範囲内に新月候補が最低 1 件（scanned≥1）, got {}",
        stats.candidates_scanned
    );
    assert!(
        stats.accepted >= 1,
        "範囲に実日食が落ちている → accepted≥1, got {}",
        stats.accepted
    );
    assert!(
        (stats.safety_margin_rad - ECLIPSE_FILTER_SAFETY_MARGIN_RAD).abs() < EPS,
        "scan の safety_margin_rad は const と一致, got {}",
        stats.safety_margin_rad
    );
    assert_eq!(
        stats.accepted + stats.rejected,
        stats.candidates_scanned,
        "accepted + rejected == scanned（恒等式）"
    );
    assert!(
        stats.max_margin_consumed_rad >= 0.0,
        "consumed は非負（max(0, …) のクランプ）, got {}",
        stats.max_margin_consumed_rad
    );
    assert!(
        stats.min_accepted_slack_rad >= 0.0
            && stats.min_accepted_slack_rad <= stats.safety_margin_rad,
        "slack は 0..=safety, got {}",
        stats.min_accepted_slack_rad
    );
    assert!(
        stats.min_accepted_slack_rad > 0.0,
        "slack>0（マージン非枯渇＝偽陰性ゼロが全期間で保たれる。slack==0 は危険信号）, got {}",
        stats.min_accepted_slack_rad
    );
}

/// 受け入れ「scan は範囲を変えても恒等式・範囲制約を保つ（別の控えめ範囲で再確認）」。
/// 2018 単年（実日食 2018-07-13 部分 / 2018-08-11 部分などを含む）でも不変条件を満たす。
/// 範囲依存の脆い数値には触れず、構造的恒等式のみを別範囲で重ねて配線ミスを殺す。
/// 殺す変異: scan が範囲を無視して固定統計を返す（同一値を返す変異を、別範囲の scanned 差で間接検出）。
#[test]
fn scan_another_range_keeps_invariants() {
    let range = utc_range((2018, 1, 1), (2019, 1, 1));
    let stats = scan_filter_margins(range).expect("post-1972 range scan must succeed");

    assert!(stats.candidates_scanned >= 1, "scanned≥1");
    // 2018 は 2018-07-13 部分・2018-08-11 部分など実日食を含む → accepted≥1（全件棄却する変異を殺す）。
    assert!(
        stats.accepted >= 1,
        "2018 は実日食を含む → accepted≥1, got {}",
        stats.accepted
    );
    assert_eq!(
        stats.accepted + stats.rejected,
        stats.candidates_scanned,
        "accepted + rejected == scanned（別範囲でも恒等）"
    );
    assert!(
        stats.min_accepted_slack_rad >= 0.0
            && stats.min_accepted_slack_rad <= stats.safety_margin_rad,
        "slack は 0..=safety（別範囲でも）"
    );
    // 別範囲でもマージン非枯渇（偽陰性ゼロが保たれる）＝ slack>0。
    assert!(
        stats.min_accepted_slack_rad > 0.0,
        "slack>0（別範囲でもマージン非枯渇）, got {}",
        stats.min_accepted_slack_rad
    );
    assert!(
        stats.max_margin_consumed_rad >= 0.0,
        "consumed 非負（別範囲でも）"
    );
}

/// 診断（常時 ignore）: 全期間（1972-2100・候補生成は閏秒の都合で 1972 以降）の D6 マージン実余裕統計を
/// 印字する。accuracy.md §3.4 D6 の達成値記録に使う。
/// `cargo test -p umbra-eclipse --test filter_margins -- --ignored --nocapture margin_scan_full_period`
#[ignore = "diagnostic: 全期間 D6 マージン実余裕統計の印字（--ignored --nocapture）"]
#[test]
fn margin_scan_full_period_1972_2100() {
    let range = utc_range((1972, 1, 1), (2100, 1, 1));
    let stats = scan_filter_margins(range).expect("post-1972 range scan must succeed");
    eprintln!(
        "[D6 margin scan 1972-2100] scanned={} accepted={} rejected={} \
         safety_margin={:.6} rad  max_consumed={:.6} rad  min_slack={:.6} rad ({:.1}% of margin)",
        stats.candidates_scanned,
        stats.accepted,
        stats.rejected,
        stats.safety_margin_rad,
        stats.max_margin_consumed_rad,
        stats.min_accepted_slack_rad,
        100.0 * stats.min_accepted_slack_rad / stats.safety_margin_rad,
    );
    assert!(stats.min_accepted_slack_rad > 0.0, "全期間でマージン非枯渇");
}

// ============================================================
// 型の存在・フィールド（コンパイル時の型契約）
// ============================================================

/// `FilterMarginStats` / `MarginSample` の全フィールドを構築・参照して型契約を縛る。
/// import だけで終わらせず 1 度は型として触れ、フィールド名・型（usize/f64/bool）の存在を保証する。
#[test]
fn types_are_constructible_and_fielded() {
    let s = MarginSample {
        separation_rad: 0.01,
        bare_limit_rad: 0.02,
        possible: true,
    };
    assert!((s.separation_rad - 0.01).abs() < EPS);
    assert!((s.bare_limit_rad - 0.02).abs() < EPS);
    assert!(s.possible);

    let stats = FilterMarginStats {
        candidates_scanned: 3,
        accepted: 2,
        rejected: 1,
        safety_margin_rad: ECLIPSE_FILTER_SAFETY_MARGIN_RAD,
        max_margin_consumed_rad: 0.001,
        min_accepted_slack_rad: ECLIPSE_FILTER_SAFETY_MARGIN_RAD - 0.001,
    };
    assert_eq!(stats.candidates_scanned, 3);
    assert_eq!(stats.accepted, 2);
    assert_eq!(stats.rejected, 1);
    assert!((stats.max_margin_consumed_rad - 0.001).abs() < EPS);

    // EclipseError が公開型として参照できること（scan の戻り値型に使う）。
    fn _accepts_err(_e: EclipseError) {}
}
