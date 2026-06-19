//! ISSUE-030 S30a 受け入れテスト（strict / 純プリミティブ: ErrorStats + ToleranceProfile）。
//!
//! 本ファイルは `umbra-fixtures` の **公開 API のみ**を対象とした統合テスト（tests/ 配下）。
//! 対象は本スライスの純プリミティブのみ:
//! - `ErrorStats`（絶対誤差ベースの記述統計: n / max_abs / mean_abs / p95 / units）
//! - `ToleranceProfile`（standard / reference の許容誤差プロファイル定数）
//!
//! ## スコープ外（後続スライス・本ファイルでは検証しない）
//! ゴールデン比較・レイヤ分解・JPL-DE 差分・1900–2100 スイープ・JSON・CLI 統合は後続。
//!
//! ## パーセンタイル規約（FIXED — R-7 / 線形補間, NumPy 既定）
//! 昇順ソート済み絶対誤差 `a[0..n]`、分数 `p`（本テストでは p = 0.95）に対し:
//! - `n == 0` → 0.0
//! - `n == 1` → `a[0]`
//! - それ以外: `h = (n-1) as f64 * p`; `lo = h.floor() as usize`;
//!   `lo + 1 >= n` なら `a[n-1]`、さもなくば `a[lo] + (h - lo) * (a[lo+1] - a[lo])`。
//!
//! 期待 p95 は上記規約に従い **本ファイル内で手計算した literal** を用いる（実装をミラーしない）。
//!
//! ## 期待される RED（実装前）
//! `report` モジュールと公開再エクスポート（`umbra_fixtures::ErrorStats` /
//! `umbra_fixtures::ToleranceProfile`）はまだ存在しないため、本ファイルは
//! **未解決インポート（E0432）でコンパイル不能 = RED** になる。これが想定どおりの赤である。

use umbra_fixtures::{ErrorStats, ToleranceProfile};

/// 計算統計の f64 一致に用いる厳密許容。
const EPS: f64 = 1e-12;

// ============================================================
// ErrorStats::from_errors — 基本（符号付き入力に abs 適用）
// ============================================================

/// 受け入れ「from_errors 基本」: `&[-3.0, 1.0, -2.0]` units "s"。
/// |e| = [3,1,2] → n=3, max_abs=3.0, mean_abs=(3+1+2)/3=2.0。
/// p95（R-7, 昇順 [1,2,3]）: h=(3-1)*0.95=1.9, lo=1, a[1]=2.0, a[2]=3.0
///   → 2.0 + 0.9*(3.0-2.0) = 2.9。
/// Kills: abs 忘れ（max/mean が変わる）, パーセンタイル規約違い, フィールド入れ替え。
#[test]
fn from_errors_basic_applies_abs_and_r7_p95() {
    let stats = ErrorStats::from_errors(&[-3.0, 1.0, -2.0], "s");
    assert_eq!(stats.n, 3, "n must equal errors.len()");
    assert!(
        (stats.max_abs - 3.0).abs() < EPS,
        "max_abs must be max|e|=3.0, got {}",
        stats.max_abs
    );
    assert!(
        (stats.mean_abs - 2.0).abs() < EPS,
        "mean_abs must be mean|e|=2.0, got {}",
        stats.mean_abs
    );
    assert!(
        (stats.p95 - 2.9).abs() < EPS,
        "p95 (R-7 over [1,2,3]) must be 2.9, got {}",
        stats.p95
    );
    assert_eq!(stats.units, "s", "units must be stored verbatim");
}

/// 受け入れ「from_errors: abs 適用（符号が効く）」: `&[-5.0, -5.0]`。
/// |e| = [5,5] → max_abs=5.0, mean_abs=5.0。
/// p95: n=2, h=(2-1)*0.95=0.95, lo=0, a[0]=5.0,a[1]=5.0 → 5.0。
/// Kills: abs を落とすと mean_abs が -5.0 になる。
#[test]
fn from_errors_negative_inputs_use_absolute_values() {
    let stats = ErrorStats::from_errors(&[-5.0, -5.0], "s");
    assert_eq!(stats.n, 2, "n must equal 2");
    assert!(
        (stats.max_abs - 5.0).abs() < EPS,
        "max_abs must be 5.0 (abs of -5.0), got {}",
        stats.max_abs
    );
    assert!(
        (stats.mean_abs - 5.0).abs() < EPS,
        "mean_abs must be 5.0 (abs applied), got {}",
        stats.mean_abs
    );
    assert!(
        (stats.p95 - 5.0).abs() < EPS,
        "p95 over [5,5] must be 5.0, got {}",
        stats.p95
    );
    assert_eq!(stats.units, "s", "units must be stored verbatim");
}

/// 受け入れ「from_errors n==1」: `&[4.2]`。
/// n=1, max_abs=4.2, mean_abs=4.2, p95=4.2（n==1 分岐 → a[0]）。
/// Kills: パーセンタイル索引パニック / 単一要素誤処理。
#[test]
fn from_errors_single_element_n1_branch() {
    let stats = ErrorStats::from_errors(&[4.2], "s");
    assert_eq!(stats.n, 1, "n must equal 1");
    assert!(
        (stats.max_abs - 4.2).abs() < EPS,
        "max_abs must be 4.2, got {}",
        stats.max_abs
    );
    assert!(
        (stats.mean_abs - 4.2).abs() < EPS,
        "mean_abs must be 4.2, got {}",
        stats.mean_abs
    );
    assert!(
        (stats.p95 - 4.2).abs() < EPS,
        "p95 must be a[0]=4.2 for n==1, got {}",
        stats.p95
    );
    assert_eq!(stats.units, "s", "units must be stored verbatim");
}

/// 受け入れ「from_errors 空列」: `&[]`。
/// n=0, max_abs=0.0, mean_abs=0.0, p95=0.0, units 保持。
/// Kills: 空列パニック, 0/0 による mean の NaN。
#[test]
fn from_errors_empty_slice_yields_zeros() {
    let stats = ErrorStats::from_errors(&[], "deg");
    assert_eq!(stats.n, 0, "n must be 0 for empty slice");
    assert!(
        stats.max_abs == 0.0,
        "max_abs must be exactly 0.0 for empty, got {}",
        stats.max_abs
    );
    assert!(
        stats.mean_abs == 0.0 && !stats.mean_abs.is_nan(),
        "mean_abs must be exactly 0.0 (not NaN) for empty, got {}",
        stats.mean_abs
    );
    assert!(
        stats.p95 == 0.0,
        "p95 must be exactly 0.0 for empty, got {}",
        stats.p95
    );
    assert_eq!(
        stats.units, "deg",
        "units must be stored verbatim even when empty"
    );
}

/// 受け入れ「from_errors p95（大きめの既知列・未ソート入力）」:
/// |e| 集合 = {0..10}（n=11）を **未ソート順**で与える。
/// h=(11-1)*0.95=9.5, lo=9, a[9]=9.0, a[10]=10.0 → 9.0 + 0.5*(10.0-9.0)=9.5。
/// max_abs=10.0, mean_abs=(0+1+...+10)/11=55/11=5.0。
/// Kills: 入力ソート前提, パーセンタイル規約違い, max/mean 取り違え。
#[test]
fn from_errors_p95_on_unsorted_known_list() {
    // 未ソート（負号も混ぜて abs 適用も同時に縛る）。|e| は {0,1,...,10}。
    let errors = [7.0, -3.0, 10.0, 0.0, 5.0, -1.0, 9.0, 2.0, -8.0, 4.0, -6.0];
    let stats = ErrorStats::from_errors(&errors, "s");
    assert_eq!(stats.n, 11, "n must equal 11");
    assert!(
        (stats.max_abs - 10.0).abs() < EPS,
        "max_abs must be 10.0, got {}",
        stats.max_abs
    );
    assert!(
        (stats.mean_abs - 5.0).abs() < EPS,
        "mean_abs must be 55/11=5.0, got {}",
        stats.mean_abs
    );
    assert!(
        (stats.p95 - 9.5).abs() < EPS,
        "p95 (R-7 over 0..10) must be 9.5, got {}",
        stats.p95
    );
}

// ============================================================
// ErrorStats::within — 境界（包含的 <=）
// ============================================================

/// 受け入れ「within 境界」: max_abs ちょうど == tolerance で true（包含的）, 直下で false。
/// `from_errors(&[2.0, -1.0], "s")` → |e|=[2,1] → max_abs=2.0。
/// `.within(2.0)` == true（<=）, `.within(1.9999)` == false。
/// Kills: `<` vs `<=` の変異, max ではなく mean を使う変異（mean_abs=1.5 なら 1.9999 で true になってしまう）。
#[test]
fn within_is_inclusive_on_max_abs() {
    let stats = ErrorStats::from_errors(&[2.0, -1.0], "s");
    assert!(
        (stats.max_abs - 2.0).abs() < EPS,
        "precondition: max_abs must be 2.0, got {}",
        stats.max_abs
    );
    assert!(
        stats.within(2.0),
        "within must be inclusive: max_abs(2.0) <= tolerance(2.0) -> true"
    );
    assert!(
        !stats.within(1.9999),
        "within must be false when max_abs(2.0) > tolerance(1.9999)"
    );
}

// ============================================================
// ToleranceProfile::standard / reference — 定数（フィールド毎に厳密）
// ============================================================

/// 受け入れ「ToleranceProfile::standard() 定数」: 6 フィールドを厳密に検証。
/// contact_seconds=2.0, maximum_seconds=1.5, magnitude=0.0005, obscuration=0.0005,
/// altitude_degrees=0.1, note_utc_is_delta_t_limited=true。
/// Kills: 任意のフィールド入れ替え / 定数違い。
#[test]
fn tolerance_profile_standard_exact_fields() {
    let p = ToleranceProfile::standard();
    assert_eq!(
        p.contact_seconds, 2.0,
        "standard.contact_seconds must be 2.0"
    );
    assert_eq!(
        p.maximum_seconds, 1.5,
        "standard.maximum_seconds must be 1.5"
    );
    assert_eq!(p.magnitude, 0.0005, "standard.magnitude must be 0.0005");
    assert_eq!(p.obscuration, 0.0005, "standard.obscuration must be 0.0005");
    assert_eq!(
        p.altitude_degrees, 0.1,
        "standard.altitude_degrees must be 0.1"
    );
    assert!(
        p.note_utc_is_delta_t_limited,
        "standard.note_utc_is_delta_t_limited must be true"
    );
}

/// 受け入れ「ToleranceProfile::reference() 定数」: 6 フィールドを厳密に検証。
/// contact_seconds=1.0, maximum_seconds=1.0, magnitude=0.0002, obscuration=0.0002,
/// altitude_degrees=0.05, note_utc_is_delta_t_limited=true。
/// Kills: フィールド入れ替え / 定数違い。
#[test]
fn tolerance_profile_reference_exact_fields() {
    let p = ToleranceProfile::reference();
    assert_eq!(
        p.contact_seconds, 1.0,
        "reference.contact_seconds must be 1.0"
    );
    assert_eq!(
        p.maximum_seconds, 1.0,
        "reference.maximum_seconds must be 1.0"
    );
    assert_eq!(p.magnitude, 0.0002, "reference.magnitude must be 0.0002");
    assert_eq!(
        p.obscuration, 0.0002,
        "reference.obscuration must be 0.0002"
    );
    assert_eq!(
        p.altitude_degrees, 0.05,
        "reference.altitude_degrees must be 0.05"
    );
    assert!(
        p.note_utc_is_delta_t_limited,
        "reference.note_utc_is_delta_t_limited must be true"
    );
}

/// 受け入れ「reference は standard より厳しい」: 時刻/食分/食面積/高度の各許容で
/// reference < standard が成立すること（contact/maximum/magnitude/obscuration/altitude）。
/// Kills: "reference == standard" 変異 / standard と reference の取り違え。
#[test]
fn reference_profile_is_strictly_tighter_than_standard() {
    let s = ToleranceProfile::standard();
    let r = ToleranceProfile::reference();
    assert!(
        r.contact_seconds < s.contact_seconds,
        "reference.contact_seconds ({}) must be < standard ({})",
        r.contact_seconds,
        s.contact_seconds
    );
    assert!(
        r.maximum_seconds < s.maximum_seconds,
        "reference.maximum_seconds ({}) must be < standard ({})",
        r.maximum_seconds,
        s.maximum_seconds
    );
    assert!(
        r.magnitude < s.magnitude,
        "reference.magnitude ({}) must be < standard ({})",
        r.magnitude,
        s.magnitude
    );
    assert!(
        r.obscuration < s.obscuration,
        "reference.obscuration ({}) must be < standard ({})",
        r.obscuration,
        s.obscuration
    );
    assert!(
        r.altitude_degrees < s.altitude_degrees,
        "reference.altitude_degrees ({}) must be < standard ({})",
        r.altitude_degrees,
        s.altitude_degrees
    );
}

// ============================================================
// units 伝播（ハードコード禁止）
// ============================================================

/// 受け入れ「units 伝播」: units が "deg" / "s" でそのまま保持される。
/// Kills: units のハードコード（特定文字列固定）。
#[test]
fn units_are_propagated_verbatim() {
    assert_eq!(
        ErrorStats::from_errors(&[1.0], "deg").units,
        "deg",
        "units 'deg' must be stored verbatim"
    );
    assert_eq!(
        ErrorStats::from_errors(&[1.0], "s").units,
        "s",
        "units 's' must be stored verbatim"
    );
}
