//! ISSUE-030 S30c 受け入れテスト（strict / 純地点別比較: compare_local + aggregate_local）。
//!
//! 本ファイルは `umbra-fixtures` の **公開 API のみ**を対象とした統合テスト（tests/ 配下）。
//! 対象は本スライスの **純地点別比較**のみ:
//! - `compare_local(computed, golden) -> LocalErrors`（符号付き = computed − golden）。
//! - `aggregate_local(errors, profile) -> LocalReport`（metric 別 `ErrorStats` ＋ 合否）。
//!
//! ## スコープ外（後続スライス・本ファイルでは検証しない）
//! エンジン駆動の `report_against_golden`・レイヤ分解・JPL-DE 差分・JSON/CLI 出力は後続スライス。
//! 本スライスは **エンジンを一切走らせず**、合成（手構築）値のみで純比較を縛る。
//!
//! ## オラクル戦略（mutation-resistant）
//! 期待値はすべてテスト側で手計算した literal。各フィールド・各分岐・符号・境界・各 gate を独立に縛る。
//! `compare_local` の符号は **computed − golden** に固定（被減数/減数の入れ替えを撃破）。
//! 接触時刻誤差は **TT 優先・なければ UTC** ＋ 2 部安全な `days_since` ×86400 を厳守する。
//! `ErrorStats::from_errors` は内部で abs を取るため、集計統計は |誤差| 上の値（手計算 R-7 p95）。
//! `LocalErrors` の各フィールドは**符号付きのまま**保持される。
//!
//! ## 期待される RED（実装前）
//! `LocalErrors` / `compare_local` / `LocalReport` / `aggregate_local` はまだ存在しないため、
//! 本ファイルは **未解決インポート（E0432/E0425）でコンパイル不能 = RED** になる。これが想定どおりの赤。

#![allow(clippy::excessive_precision)]

use umbra_core::{Degrees, JulianDate2, TtInstant, UtcInstant};
use umbra_eclipse::{
    AccuracyProfile, CalculationMetadata, EclipseMagnitude, LocalCircumstances, LocalContact,
    LocalContactSet, Obscuration, Visibility,
};

use umbra_fixtures::{
    aggregate_local, compare_local, ErrorStats, GoldenContact, GoldenLocation, LocalErrors,
    LocalReport, LocationClass, ToleranceProfile,
};

/// 計算統計の f64 一致に用いる厳密許容。
const EPS: f64 = 1e-12;
/// TT/UTC 秒差（×86400 換算）の許容。`Δs/86400` の往復で浮動小数の丸めが出るため緩め。
const SEC_EPS: f64 = 1e-6;

// ============================================================
// 構築ヘルパ（results.rs / types.rs のフィールドをミラー）
// ============================================================

/// UTC 瞬時を整数引数で組む。
fn utc(year: i32, month: u8, day: u8, hour: u8, minute: u8, second: f64) -> UtcInstant {
    UtcInstant::from_gregorian(year, month, day, hour, minute, second).expect("有効な UTC 日時")
}

/// TT 瞬時を 2 要素 JD で組む。`Δs/86400` を part2 に足して既知秒差を作る。
fn tt(jd1: f64, jd2: f64) -> TtInstant {
    TtInstant::from_jd2(JulianDate2::new(jd1, jd2))
}

/// 固定のメタデータ（フィラー。compare_local は読まない）。
fn metadata() -> CalculationMetadata {
    CalculationMetadata {
        library_version: "0.1.0".to_string(),
        ephemeris_model: "ELP/MPP02+VSOP87D".to_string(),
        ephemeris_version: "2024a".to_string(),
        delta_t_model: "EspenakMeeus".to_string(),
        delta_t_uncertainty_seconds: 0.5,
        earth_model: "WGS84".to_string(),
        lunar_radius_model: "IauMean".to_string(),
        accuracy_profile: AccuracyProfile::Standard,
        generated_at: utc(2026, 6, 18, 0, 0, 0.0),
    }
}

/// 計算側の局地接触。時刻 2 値以外（高度/方位/位置角/visible）は固定フィラー（compare_local は読まない）。
fn local_contact(time_utc: UtcInstant, time_tt: TtInstant) -> LocalContact {
    LocalContact {
        time_utc,
        time_tt,
        sun_altitude: Degrees(40.0),
        sun_azimuth: Degrees(200.0),
        position_angle: Degrees(300.0),
        visible: true,
    }
}

/// 計算側の局地接触で **visible フラグを明示**したもの（地平下接触の合成用）。
/// 時刻 2 値と visible 以外は固定フィラー。USNO オラクルは地平下接触を None として省くため、
/// `visible:false` の computed 接触 vs golden None は「不一致でない」べき（本スライスの是正対象）。
fn local_contact_vis(time_utc: UtcInstant, time_tt: TtInstant, visible: bool) -> LocalContact {
    LocalContact {
        time_utc,
        time_tt,
        sun_altitude: Degrees(40.0),
        sun_azimuth: Degrees(200.0),
        position_angle: Degrees(300.0),
        visible,
    }
}

/// 合成 `computed: LocalCircumstances`。compare_local が読む値だけを引数で指定し、metadata はフィラー。
fn computed_local(
    contacts: LocalContactSet,
    magnitude: f64,
    obscuration: f64,
    max_alt: f64,
    visibility: Visibility,
) -> LocalCircumstances {
    LocalCircumstances {
        contacts,
        magnitude: EclipseMagnitude(magnitude),
        obscuration: Obscuration(obscuration),
        maximum_altitude: Degrees(max_alt),
        visibility,
        metadata: metadata(),
    }
}

/// ゴールデン接触（UTC・任意 TT・高度）。compare_local は time_tt / time_utc のみ読む。
fn golden_contact(
    time_utc: UtcInstant,
    time_tt: Option<TtInstant>,
    altitude_deg: f64,
) -> GoldenContact {
    GoldenContact {
        time_utc,
        time_tt,
        altitude_deg,
    }
}

/// 合成 `golden: GoldenLocation`。compare_local が読む値（c1..c4 / maximum / magnitude /
/// obscuration / max_altitude_deg / visibility_expected）だけを引数で指定し、残りはフィラー。
#[allow(clippy::too_many_arguments)]
fn golden_location(
    c1: Option<GoldenContact>,
    c2: Option<GoldenContact>,
    maximum: GoldenContact,
    c3: Option<GoldenContact>,
    c4: Option<GoldenContact>,
    magnitude: f64,
    obscuration: f64,
    max_altitude_deg: f64,
    visibility_expected: Visibility,
) -> GoldenLocation {
    GoldenLocation {
        name: "x".to_string(),
        latitude_deg: 35.0,
        east_longitude_deg: 139.0,
        elevation_m: 0.0,
        location_class: LocationClass::Centerline,
        c1,
        c2,
        maximum,
        c3,
        c4,
        magnitude,
        obscuration,
        max_altitude_deg,
        max_azimuth_deg: 200.0,
        visibility_expected,
    }
}

/// 既定の computed contacts（全 Some・全時刻 golden と一致）。各テストが上書きする土台。
/// part2=0 を基準に、golden 側も同じ JD を使うことで「差 0」を作る。
fn aligned_set(base_tt: TtInstant, base_utc: UtcInstant) -> LocalContactSet {
    LocalContactSet {
        c1: Some(local_contact(base_utc, base_tt)),
        c2: Some(local_contact(base_utc, base_tt)),
        maximum: local_contact(base_utc, base_tt),
        c3: Some(local_contact(base_utc, base_tt)),
        c4: Some(local_contact(base_utc, base_tt)),
    }
}

// ============================================================
// compare_local — maximum_seconds（TT 優先・UTC フォールバック・符号）
// ============================================================

/// 受け入れ「maximum_seconds: golden.maximum に TT があれば TT 優先・符号は computed − golden」。
/// golden maximum TT = base、computed maximum TT = base + 2.0 s → +2.0。computed の UTC は
/// **わざと遠い別日**にして TT 分岐が UTC を読んでいないことを証明する。
/// 別ケース: golden.maximum.time_tt = None → UTC フォールバックで computed UTC が golden UTC より
/// −1.0 s 前 → −1.0（computed TT は無関係な別値で「常に TT を読む」変異を撃破）。
/// 殺す変異: 符号反転（golden − computed）, TT があるのに UTC を使う / 常に TT を読む, ×86400 欠落,
/// `.jd()` 桁落ち（days_since の精度を 1e-6 で要求）。
#[test]
fn compare_local_maximum_seconds_tt_preferred_and_utc_fallback_signed() {
    let base_tt = tt(2_451_545.0, 0.0);
    let g_utc = utc(2024, 4, 8, 18, 0, 0.0);

    // --- TT 優先ケース: golden maximum に TT あり。computed TT = base + 2.0 s。---
    let golden_tt = golden_location(
        None,
        None,
        golden_contact(g_utc, Some(base_tt), 50.0),
        None,
        None,
        1.0,
        1.0,
        50.0,
        Visibility::FullyVisible,
    );
    let computed_tt = computed_local(
        LocalContactSet {
            c1: None,
            c2: None,
            // computed TT は +2.0 s 後。UTC は遠い別日（TT 優先の証明）。
            maximum: local_contact(utc(2000, 1, 1, 0, 0, 0.0), tt(2_451_545.0, 2.0 / 86400.0)),
            c3: None,
            c4: None,
        },
        1.0,
        1.0,
        50.0,
        Visibility::FullyVisible,
    );
    let e_tt = compare_local(&computed_tt, &golden_tt);
    assert!(
        (e_tt.maximum_seconds - 2.0).abs() < SEC_EPS,
        "TT 優先・computed が +2.0 s 後 → maximum_seconds ≈ +2.0（computed − golden）, got {}",
        e_tt.maximum_seconds
    );

    // --- UTC フォールバックケース: golden maximum の TT = None。computed UTC = golden UTC − 1.0 s。---
    let golden_utc = golden_location(
        None,
        None,
        golden_contact(utc(2024, 4, 8, 18, 0, 5.0), None, 50.0),
        None,
        None,
        1.0,
        1.0,
        50.0,
        Visibility::FullyVisible,
    );
    let computed_utc = computed_local(
        LocalContactSet {
            c1: None,
            c2: None,
            // computed UTC は golden UTC(…5.0s) より −1.0 s 前。computed TT は無関係（読まれないはず）。
            maximum: local_contact(utc(2024, 4, 8, 18, 0, 4.0), tt(2_400_000.0, 0.0)),
            c3: None,
            c4: None,
        },
        1.0,
        1.0,
        50.0,
        Visibility::FullyVisible,
    );
    let e_utc = compare_local(&computed_utc, &golden_utc);
    assert!(
        (e_utc.maximum_seconds + 1.0).abs() < SEC_EPS,
        "TT=None → UTC フォールバックで −1.0 s（computed − golden）, got {}",
        e_utc.maximum_seconds
    );
}

// ============================================================
// compare_local — contact_seconds（時系列順・両 Some 収集）
// ============================================================

/// 受け入れ「c1..c4 が両方 Some なら時刻誤差を c1→c4 の順で contact_seconds に積む」。
/// c1=+0.5, c2=+1.0, c3=+1.5, c4=+2.0 s の **互いに異なる**既知誤差を与え、順序・取りこぼし・符号を縛る。
/// 全 golden 接触に TT を付与（TT 分岐）。computed maximum/その他フィールドは差 0 にして単離。
/// 殺す変異: 順序入替（c1↔c4 等）, 接触の取りこぼし, 符号反転, ×86400 欠落, presence 誤カウント。
#[test]
fn compare_local_contacts_both_some_collected_in_order() {
    let g_utc = utc(2024, 4, 8, 18, 0, 0.0);
    let base = tt(2_451_545.0, 0.0);

    // golden の 4 接触はすべて同一基準 TT（part2=0）。
    let gc = |alt: f64| golden_contact(g_utc, Some(base), alt);
    let golden = golden_location(
        Some(gc(10.0)),
        Some(gc(20.0)),
        gc(50.0), // maximum: 差 0
        Some(gc(30.0)),
        Some(gc(40.0)),
        1.0,
        1.0,
        50.0,
        Visibility::FullyVisible,
    );

    // computed の 4 接触は base + {0.5,1.0,1.5,2.0} s。maximum は base（差 0）。
    let off = |s: f64| local_contact(g_utc, tt(2_451_545.0, s / 86400.0));
    let computed = computed_local(
        LocalContactSet {
            c1: Some(off(0.5)),
            c2: Some(off(1.0)),
            maximum: off(0.0),
            c3: Some(off(1.5)),
            c4: Some(off(2.0)),
        },
        1.0,
        1.0,
        50.0,
        Visibility::FullyVisible,
    );
    let e = compare_local(&computed, &golden);
    assert_eq!(
        e.contact_seconds.len(),
        4,
        "両 Some の 4 接触すべてが収集される, got {:?}",
        e.contact_seconds
    );
    let expected = [0.5, 1.0, 1.5, 2.0];
    for (i, &exp) in expected.iter().enumerate() {
        assert!(
            (e.contact_seconds[i] - exp).abs() < SEC_EPS,
            "contact_seconds[{i}] は {exp}（c1→c4 順, computed − golden）, got {}",
            e.contact_seconds[i]
        );
    }
    assert_eq!(
        e.contact_presence_mismatches, 0,
        "全接触が両 Some → presence mismatch は 0"
    );
}

/// 受け入れ「presence 不一致（片側のみ Some）は contact_presence_mismatches に計上し、
/// contact_seconds には積まない」。computed c2=Some/golden c2=None、computed c3=None/golden c3=Some
/// の 2 件が不一致。c1/c4 は両 Some（既知誤差 +0.5 / +2.0）。
/// → contact_presence_mismatches==2、contact_seconds は c1,c4 のみで len 2（順序保持）。
/// 殺す変異: presence 不一致を数えない, 不一致接触を contact_seconds に積む, カウント値の誤り。
#[test]
fn compare_local_contact_presence_mismatch_counted() {
    let g_utc = utc(2024, 4, 8, 18, 0, 0.0);
    let base = tt(2_451_545.0, 0.0);
    let gc = |alt: f64| golden_contact(g_utc, Some(base), alt);
    let off = |s: f64| local_contact(g_utc, tt(2_451_545.0, s / 86400.0));

    let golden = golden_location(
        Some(gc(10.0)), // c1: 両 Some
        None,           // c2: golden None（computed Some → 不一致）
        gc(50.0),
        Some(gc(30.0)), // c3: golden Some（computed None → 不一致）
        Some(gc(40.0)), // c4: 両 Some
        1.0,
        1.0,
        50.0,
        Visibility::FullyVisible,
    );
    let computed = computed_local(
        LocalContactSet {
            c1: Some(off(0.5)), // 両 Some → +0.5
            c2: Some(off(9.0)), // computed Some / golden None → 不一致（積まれない）
            maximum: off(0.0),
            c3: None,           // computed None / golden Some → 不一致
            c4: Some(off(2.0)), // 両 Some → +2.0
        },
        1.0,
        1.0,
        50.0,
        Visibility::FullyVisible,
    );
    let e = compare_local(&computed, &golden);
    assert_eq!(
        e.contact_presence_mismatches, 2,
        "c2(computed のみ Some)+c3(golden のみ Some)= 2 件の presence 不一致"
    );
    assert_eq!(
        e.contact_seconds.len(),
        2,
        "両 Some は c1,c4 のみ → len 2（不一致接触は積まない）, got {:?}",
        e.contact_seconds
    );
    assert!(
        (e.contact_seconds[0] - 0.5).abs() < SEC_EPS,
        "contact_seconds[0] は c1 の +0.5, got {}",
        e.contact_seconds[0]
    );
    assert!(
        (e.contact_seconds[1] - 2.0).abs() < SEC_EPS,
        "contact_seconds[1] は c4 の +2.0, got {}",
        e.contact_seconds[1]
    );
}

/// 受け入れ「両側 None の接触は不一致でも収集でもなく、ただ無視される」。
/// c2,c3 を **両側 None**、c1,c4 は両 Some（+0.5 / +2.0）。
/// → contact_presence_mismatches==0、contact_seconds は c1,c4 のみで len 2。
/// 殺す変異: 両 None を presence 不一致として数える, 両 None を contact_seconds に積む。
#[test]
fn compare_local_both_none_contacts_skipped() {
    let g_utc = utc(2024, 4, 8, 18, 0, 0.0);
    let base = tt(2_451_545.0, 0.0);
    let gc = |alt: f64| golden_contact(g_utc, Some(base), alt);
    let off = |s: f64| local_contact(g_utc, tt(2_451_545.0, s / 86400.0));

    let golden = golden_location(
        Some(gc(10.0)), // c1: 両 Some
        None,           // c2: 両 None
        gc(50.0),
        None,           // c3: 両 None
        Some(gc(40.0)), // c4: 両 Some
        1.0,
        1.0,
        50.0,
        Visibility::FullyVisible,
    );
    let computed = computed_local(
        LocalContactSet {
            c1: Some(off(0.5)),
            c2: None, // 両 None
            maximum: off(0.0),
            c3: None, // 両 None
            c4: Some(off(2.0)),
        },
        1.0,
        1.0,
        50.0,
        Visibility::FullyVisible,
    );
    let e = compare_local(&computed, &golden);
    assert_eq!(
        e.contact_presence_mismatches, 0,
        "両 None は presence 不一致ではない → 0"
    );
    assert_eq!(
        e.contact_seconds.len(),
        2,
        "両 None は積まれない → c1,c4 のみ len 2, got {:?}",
        e.contact_seconds
    );
    assert!(
        (e.contact_seconds[0] - 0.5).abs() < SEC_EPS,
        "contact_seconds[0]=c1 の +0.5, got {}",
        e.contact_seconds[0]
    );
    assert!(
        (e.contact_seconds[1] - 2.0).abs() < SEC_EPS,
        "contact_seconds[1]=c4 の +2.0, got {}",
        e.contact_seconds[1]
    );
}

// ============================================================
// compare_local — 地平下接触の presence 規約（USNO は地平下接触を None で省く）
//   computed が Some{visible:false} ＝ 地平下接触のとき、golden None は「不一致でない」。
//   visible:true の Some vs golden None は従来どおり不一致。
// ============================================================

/// 受け入れ「地平下の computed 接触（Some{visible:false}）vs golden None は presence 不一致でない」。
/// 日没皆既（例 2002-12-04 の Ceduna/Lyndhurst/Adelaide）の C4 部分食終了は地平下で起こり、
/// エンジンは `Some{visible:false}` を返すが USNO オラクルは省いて `None` を格納する。両表現は
/// 「観測不能」で一致しているので不一致にしてはならない。
/// 構成: computed C4 = Some(visible:false)、golden C4 = None。他接触は両 None で単離。
/// → contact_presence_mismatches == 0（現行コードは 1 と数える＝本テストが PRIMARY RED）。
/// 殺す変異/バグ: `(Some,None)` を無条件で不一致計上する現行実装、`visible` ガードの脱落・反転。
#[test]
fn compare_local_below_horizon_contact_vs_golden_none_not_mismatch() {
    let g_utc = utc(2002, 12, 4, 7, 0, 0.0);
    let base = tt(2_452_612.0, 0.0);

    // golden: C4 は None（USNO が地平下接触を省く）。他接触も None。maximum のみ存在。
    let golden = golden_location(
        None,
        None,
        golden_contact(g_utc, Some(base), 1.0),
        None,
        None, // c4: golden None
        1.0,
        1.0,
        1.0,
        Visibility::FullyVisible,
    );
    // computed: C4 は Some だが地平下（visible:false）。
    let computed = computed_local(
        LocalContactSet {
            c1: None,
            c2: None,
            maximum: local_contact_vis(g_utc, base, true),
            c3: None,
            c4: Some(local_contact_vis(g_utc, base, false)), // 地平下接触
        },
        1.0,
        1.0,
        1.0,
        Visibility::FullyVisible,
    );
    let e = compare_local(&computed, &golden);
    assert_eq!(
        e.contact_presence_mismatches, 0,
        "地平下 computed C4(visible:false) vs golden None は不一致でない → 0, got {}",
        e.contact_presence_mismatches
    );
    // 他フィールドが地平下規約の影響を受けないことの単離確認。
    assert!(
        e.contact_seconds.is_empty(),
        "両 Some の接触は無いので contact_seconds は空, got {:?}",
        e.contact_seconds
    );
    assert!(
        (e.maximum_seconds).abs() < SEC_EPS,
        "maximum は差 0（地平下規約と独立）, got {}",
        e.maximum_seconds
    );
}

/// 受け入れ「地平**上**の computed 接触（Some{visible:true}）vs golden None は従来どおり不一致」。
/// 可視接触をエンジンが返したのに USNO が省くのは真の食い違い（過剰抑制ガード）。
/// 構成: computed C4 = Some(visible:true)、golden C4 = None。他は両 None。
/// → contact_presence_mismatches == 1（現行でも通る。地平下抑制が visible:true まで漏れないことを縛る）。
/// 殺す変異/バグ: `visible` を見ずに常に抑制する、ガード条件を `visible==false` でなく恒真にする実装。
#[test]
fn compare_local_above_horizon_contact_vs_golden_none_is_mismatch() {
    let g_utc = utc(2024, 4, 8, 18, 0, 0.0);
    let base = tt(2_451_545.0, 0.0);

    let golden = golden_location(
        None,
        None,
        golden_contact(g_utc, Some(base), 50.0),
        None,
        None, // c4: golden None
        1.0,
        1.0,
        50.0,
        Visibility::FullyVisible,
    );
    let computed = computed_local(
        LocalContactSet {
            c1: None,
            c2: None,
            maximum: local_contact_vis(g_utc, base, true),
            c3: None,
            c4: Some(local_contact_vis(g_utc, base, true)), // 地平上の可視接触
        },
        1.0,
        1.0,
        50.0,
        Visibility::FullyVisible,
    );
    let e = compare_local(&computed, &golden);
    assert_eq!(
        e.contact_presence_mismatches, 1,
        "地平上 computed C4(visible:true) vs golden None は真の不一致 → 1, got {}",
        e.contact_presence_mismatches
    );
}

/// 受け入れ「computed None vs golden Some は（visible 規約に関わらず）従来どおり不一致」。
/// golden だけが接触を持つ＝エンジンの取りこぼしであり、地平下抑制の対象外。
/// 構成: computed C4 = None、golden C4 = Some。他は両 None。
/// → contact_presence_mismatches == 1（不変）。
/// 殺す変異/バグ: 地平下抑制を `(None,Some)` 側へ誤って広げる実装。
#[test]
fn compare_local_computed_none_vs_golden_some_is_mismatch_unchanged() {
    let g_utc = utc(2024, 4, 8, 18, 0, 0.0);
    let base = tt(2_451_545.0, 0.0);

    let golden = golden_location(
        None,
        None,
        golden_contact(g_utc, Some(base), 50.0),
        None,
        Some(golden_contact(g_utc, Some(base), 50.0)), // c4: golden Some
        1.0,
        1.0,
        50.0,
        Visibility::FullyVisible,
    );
    let computed = computed_local(
        LocalContactSet {
            c1: None,
            c2: None,
            maximum: local_contact_vis(g_utc, base, true),
            c3: None,
            c4: None, // computed None
        },
        1.0,
        1.0,
        50.0,
        Visibility::FullyVisible,
    );
    let e = compare_local(&computed, &golden);
    assert_eq!(
        e.contact_presence_mismatches, 1,
        "computed None vs golden Some は不一致（不変）→ 1, got {}",
        e.contact_presence_mismatches
    );
}

/// 受け入れ「両 Some の接触は visible に関係なく時刻比較され、presence 不一致にならない」。
/// 地平下抑制は **golden None のときだけ**効く。両 Some なら（地平下でも）時刻誤差を積む（不変）。
/// 構成: computed C4 = Some(visible:false, +2.0 s)、golden C4 = Some。他は両 None。
/// → contact_presence_mismatches == 0、contact_seconds == [+2.0]（時刻比較が生きている）。
/// 殺す変異/バグ: 地平下接触を両 Some でも抑制してしまい時刻誤差を取りこぼす実装。
#[test]
fn compare_local_below_horizon_both_some_still_time_compared() {
    let g_utc = utc(2024, 4, 8, 18, 0, 0.0);
    let base = tt(2_451_545.0, 0.0);

    let golden = golden_location(
        None,
        None,
        golden_contact(g_utc, Some(base), 50.0),
        None,
        Some(golden_contact(g_utc, Some(base), -1.0)), // c4: golden Some（地平下高度でも存在）
        1.0,
        1.0,
        50.0,
        Visibility::FullyVisible,
    );
    let computed = computed_local(
        LocalContactSet {
            c1: None,
            c2: None,
            maximum: local_contact_vis(g_utc, base, true),
            c3: None,
            // computed C4: 地平下（visible:false）だが golden も Some → 時刻比較される。+2.0 s。
            c4: Some(local_contact_vis(
                g_utc,
                tt(2_451_545.0, 2.0 / 86400.0),
                false,
            )),
        },
        1.0,
        1.0,
        50.0,
        Visibility::FullyVisible,
    );
    let e = compare_local(&computed, &golden);
    assert_eq!(
        e.contact_presence_mismatches, 0,
        "両 Some は presence 不一致でない（visible 無関係）→ 0, got {}",
        e.contact_presence_mismatches
    );
    assert_eq!(
        e.contact_seconds.len(),
        1,
        "両 Some の C4 は時刻誤差として積まれる, got {:?}",
        e.contact_seconds
    );
    assert!(
        (e.contact_seconds[0] - 2.0).abs() < SEC_EPS,
        "C4 の時刻誤差 = +2.0 s（地平下でも両 Some なら比較される）, got {}",
        e.contact_seconds[0]
    );
}

/// 受け入れ「C1..C4 混在: 地平下 Some vs None は抑制、地平上 Some vs None と None vs Some は計上」。
/// 構成（全 4 接触を使う）:
///   C1 = computed Some(visible:false) / golden None → 抑制（不一致でない）。
///   C2 = computed Some(visible:true)  / golden None → 不一致。
///   C3 = computed None / golden Some                → 不一致。
///   C4 = computed Some / golden Some（+1.5 s）       → 両 Some・時刻比較（不一致でない）。
/// → contact_presence_mismatches == 2（C2,C3 のみ）、contact_seconds == [+1.5]（C4 のみ）。
/// 殺す変異/バグ: 地平下ガードの取りこぼし/過剰適用、presence と time の取り違え、件数誤り。
#[test]
fn compare_local_mixed_contacts_only_genuine_mismatches_counted() {
    let g_utc = utc(2002, 12, 4, 7, 0, 0.0);
    let base = tt(2_452_612.0, 0.0);

    let golden = golden_location(
        None,                                         // c1: golden None
        None,                                         // c2: golden None
        golden_contact(g_utc, Some(base), 5.0),       // maximum
        Some(golden_contact(g_utc, Some(base), 5.0)), // c3: golden Some
        Some(golden_contact(g_utc, Some(base), 5.0)), // c4: golden Some
        1.0,
        1.0,
        5.0,
        Visibility::FullyVisible,
    );
    let computed = computed_local(
        LocalContactSet {
            c1: Some(local_contact_vis(g_utc, base, false)), // 地平下 vs None → 抑制
            c2: Some(local_contact_vis(g_utc, base, true)),  // 地平上 vs None → 不一致
            maximum: local_contact_vis(g_utc, base, true),
            c3: None, // None vs Some → 不一致
            // c4: 両 Some。+1.5 s。
            c4: Some(local_contact_vis(
                g_utc,
                tt(2_452_612.0, 1.5 / 86400.0),
                true,
            )),
        },
        1.0,
        1.0,
        5.0,
        Visibility::FullyVisible,
    );
    let e = compare_local(&computed, &golden);
    assert_eq!(
        e.contact_presence_mismatches, 2,
        "真の不一致は C2(地平上 vs None) と C3(None vs Some) の 2 件のみ, got {}",
        e.contact_presence_mismatches
    );
    assert_eq!(
        e.contact_seconds.len(),
        1,
        "両 Some は C4 のみ → contact_seconds len 1, got {:?}",
        e.contact_seconds
    );
    assert!(
        (e.contact_seconds[0] - 1.5).abs() < SEC_EPS,
        "C4 の時刻誤差 = +1.5 s, got {}",
        e.contact_seconds[0]
    );
}

// ============================================================
// compare_local — magnitude / obscuration / max_altitude（符号付き・独立）
// ============================================================

/// 受け入れ「magnitude / obscuration / max_altitude_deg は computed − golden の符号付き差で、
/// 互いに異なる非ゼロ値を取り違えずに各フィールドへ配線する」。
/// magnitude: 1.01 − 1.00 = +0.01、obscuration: 0.98 − 0.99 = −0.01、max_alt: 60.5 − 60.0 = +0.5。
/// 接触はすべて差 0（時刻フィールドと混線しないことも単離）。各フィールドを独立 assert。
/// 殺す変異: 符号反転, フィールド間ミスルーティング（magnitude↔obscuration↔altitude）,
/// EclipseMagnitude(.0)/Obscuration(.0)/Degrees(.0) の取り違え。
#[test]
fn compare_local_magnitude_obscuration_altitude_signed() {
    let g_utc = utc(2024, 4, 8, 18, 0, 0.0);
    let base = tt(2_451_545.0, 0.0);

    let golden = golden_location(
        None,
        None,
        golden_contact(g_utc, Some(base), 50.0),
        None,
        None,
        1.00, // golden magnitude
        0.99, // golden obscuration
        60.0, // golden max_altitude_deg
        Visibility::FullyVisible,
    );
    let computed = computed_local(
        aligned_set(base, g_utc),
        1.01, // computed magnitude → +0.01
        0.98, // computed obscuration → −0.01
        60.5, // computed max_alt → +0.5
        Visibility::FullyVisible,
    );
    let e = compare_local(&computed, &golden);
    assert!(
        (e.magnitude - 0.01).abs() < EPS,
        "magnitude = 1.01 − 1.00 = +0.01（obscuration/altitude と混線しない）, got {}",
        e.magnitude
    );
    assert!(
        (e.obscuration + 0.01).abs() < EPS,
        "obscuration = 0.98 − 0.99 = −0.01（符号・他フィールドと混線しない）, got {}",
        e.obscuration
    );
    assert!(
        (e.max_altitude_deg - 0.5).abs() < EPS,
        "max_altitude_deg = 60.5 − 60.0 = +0.5（他フィールドと混線しない）, got {}",
        e.max_altitude_deg
    );
}

/// 受け入れ「visibility_matches は computed.visibility == golden.visibility_expected」。
/// 一致（両 FullyVisible）→ true、不一致（computed FullyVisible / golden PartialVisible）→ false。
/// 殺す変異: 比較の反転（!=）, 常に true/false を返す固定化。
#[test]
fn compare_local_visibility_matches_true_and_false() {
    let g_utc = utc(2024, 4, 8, 18, 0, 0.0);
    let base = tt(2_451_545.0, 0.0);

    // 一致ケース。
    let golden_match = golden_location(
        None,
        None,
        golden_contact(g_utc, Some(base), 50.0),
        None,
        None,
        1.0,
        1.0,
        50.0,
        Visibility::FullyVisible,
    );
    let computed_match = computed_local(
        aligned_set(base, g_utc),
        1.0,
        1.0,
        50.0,
        Visibility::FullyVisible,
    );
    assert!(
        compare_local(&computed_match, &golden_match).visibility_matches,
        "computed.visibility == golden.visibility_expected → true"
    );

    // 不一致ケース。
    let golden_mismatch = golden_location(
        None,
        None,
        golden_contact(g_utc, Some(base), 50.0),
        None,
        None,
        1.0,
        1.0,
        50.0,
        Visibility::PartialVisible,
    );
    let computed_mismatch = computed_local(
        aligned_set(base, g_utc),
        1.0,
        1.0,
        50.0,
        Visibility::FullyVisible,
    );
    assert!(
        !compare_local(&computed_mismatch, &golden_mismatch).visibility_matches,
        "computed(FullyVisible) != golden(PartialVisible) → false"
    );
}

// ============================================================
// aggregate_local — metric 別 ErrorStats ＋ 件数 ＋ 合否
//   （LocalErrors を直接構築して単離する）
// ============================================================

/// `LocalErrors` を必要フィールドだけ指定して組むビルダ（残りは中立値）。
#[allow(clippy::too_many_arguments)]
fn errs(
    maximum_seconds: f64,
    contact_seconds: Vec<f64>,
    contact_presence_mismatches: usize,
    magnitude: f64,
    obscuration: f64,
    max_altitude_deg: f64,
    visibility_matches: bool,
) -> LocalErrors {
    LocalErrors {
        maximum_seconds,
        contact_seconds,
        contact_presence_mismatches,
        magnitude,
        obscuration,
        max_altitude_deg,
        visibility_matches,
    }
}

/// 受け入れ「contacts ErrorStats は全地点・全接触をフラットに集計する」。
/// 2 地点: contact_seconds=[1.0,-2.0] と [0.5]。|e|=[1,2,0.5] → n=3, max_abs=2.0,
/// mean_abs=(1+2+0.5)/3, 昇順[0.5,1,2] の R-7 p95: h=(3-1)*0.95=1.9, lo=1 → 1.0+0.9*(2.0-1.0)=1.9。
/// units="s"。
/// 殺す変異: フラット化しない（地点ごと別集計）, units 取り違え, abs 忘れ（max が 2.0 でなくなる）。
#[test]
fn aggregate_local_flattens_contacts_across_locations() {
    let errors = [
        errs(0.0, vec![1.0, -2.0], 0, 0.0, 0.0, 0.0, true),
        errs(0.0, vec![0.5], 0, 0.0, 0.0, 0.0, true),
    ];
    let report: LocalReport = aggregate_local(&errors, &ToleranceProfile::standard());
    assert_eq!(report.contacts.n, 3, "全接触フラットで n=3");
    assert!(
        (report.contacts.max_abs - 2.0).abs() < EPS,
        "contacts.max_abs = |−2.0| = 2.0, got {}",
        report.contacts.max_abs
    );
    let mean = (1.0 + 2.0 + 0.5) / 3.0;
    assert!(
        (report.contacts.mean_abs - mean).abs() < EPS,
        "contacts.mean_abs = (1+2+0.5)/3, got {}",
        report.contacts.mean_abs
    );
    assert!(
        (report.contacts.p95 - 1.9).abs() < EPS,
        "contacts.p95 (R-7 over [0.5,1,2]) = 1.9, got {}",
        report.contacts.p95
    );
    assert_eq!(report.contacts.units, "s", "contacts.units は 's'");
}

/// 受け入れ「metric→ErrorStats の配線と単位（maximum 's' / magnitude '' / obscuration '' /
/// max_altitude 'deg'）」。3 地点を異なる既知値で構築し、metric ごとに max_abs と units を縛る。
/// maximum_seconds=[1.0,-3.0,2.0] → max_abs=3.0, units "s"。
/// magnitude=[0.0003,-0.0002,0.0001] → max_abs=0.0003, units ""。
/// obscuration=[-0.0004,0.0001,0.0002] → max_abs=0.0004, units ""。
/// max_altitude_deg=[0.2,-0.5,0.1] → max_abs=0.5, units "deg"。
/// 殺す変異: metric↔stats の取り違え（maximum→magnitude 等）, 単位取り違え, abs 忘れ。
#[test]
fn aggregate_local_per_metric_stats_and_units() {
    let errors = [
        errs(1.0, vec![], 0, 0.0003, -0.0004, 0.2, true),
        errs(-3.0, vec![], 0, -0.0002, 0.0001, -0.5, true),
        errs(2.0, vec![], 0, 0.0001, 0.0002, 0.1, true),
    ];
    let report = aggregate_local(&errors, &ToleranceProfile::standard());

    // maximum（units "s"）。
    assert_eq!(report.maximum.n, 3, "maximum.n=3");
    assert!(
        (report.maximum.max_abs - 3.0).abs() < EPS,
        "maximum.max_abs = |−3.0| = 3.0, got {}",
        report.maximum.max_abs
    );
    assert_eq!(report.maximum.units, "s", "maximum.units は 's'");

    // magnitude（units ""）。
    assert_eq!(report.magnitude.n, 3, "magnitude.n=3");
    assert!(
        (report.magnitude.max_abs - 0.0003).abs() < EPS,
        "magnitude.max_abs = 0.0003, got {}",
        report.magnitude.max_abs
    );
    assert_eq!(report.magnitude.units, "", "magnitude.units は ''");

    // obscuration（units ""）。
    assert_eq!(report.obscuration.n, 3, "obscuration.n=3");
    assert!(
        (report.obscuration.max_abs - 0.0004).abs() < EPS,
        "obscuration.max_abs = |−0.0004| = 0.0004, got {}",
        report.obscuration.max_abs
    );
    assert_eq!(report.obscuration.units, "", "obscuration.units は ''");

    // max_altitude（units "deg"）。
    assert_eq!(report.max_altitude.n, 3, "max_altitude.n=3");
    assert!(
        (report.max_altitude.max_abs - 0.5).abs() < EPS,
        "max_altitude.max_abs = |−0.5| = 0.5, got {}",
        report.max_altitude.max_abs
    );
    assert_eq!(
        report.max_altitude.units, "deg",
        "max_altitude.units は 'deg'"
    );
}

/// 受け入れ「visibility_mismatches は !visibility_matches の地点数、contact_presence_mismatches は
/// 全地点の合計」。4 地点: visibility_matches=[true,false,false,false]（false 3 件・true 1 件で**非対称**）、
/// contact_presence_mismatches=[1,0,2,1]（合計 4）。
/// 殺す変異: true を数える（mismatch でなく match をカウント＝非対称ゆえ 1≠3 で撃破）,
///   presence の合計を取らない / 件数誤り。
#[test]
fn aggregate_local_counts_visibility_and_presence_mismatches() {
    let errors = [
        errs(0.0, vec![], 1, 0.0, 0.0, 0.0, true),
        errs(0.0, vec![], 0, 0.0, 0.0, 0.0, false),
        errs(0.0, vec![], 2, 0.0, 0.0, 0.0, false),
        errs(0.0, vec![], 1, 0.0, 0.0, 0.0, false),
    ];
    let report = aggregate_local(&errors, &ToleranceProfile::standard());
    assert_eq!(
        report.visibility_mismatches, 3,
        "!visibility_matches は 3 件（false の数・match カウントなら 1 で撃破）"
    );
    assert_eq!(
        report.contact_presence_mismatches, 4,
        "presence 不一致の合計 = 1+0+2+1 = 4"
    );
}

/// 受け入れ「全 metric が standard 内・可視一致・presence 不一致 0 なら pass==true。
/// かつ pass でも統計は populated（誤差を隠さない）」。
/// 2 地点・全 metric は standard 許容内（maximum ≤1.5, contacts ≤2.0, magnitude ≤0.0005,
/// obscuration ≤0.0005, altitude ≤0.1）、visibility 全一致、presence 0。
/// 殺す変異: pass が常に false, pass が統計をゼロ化/空化。
#[test]
fn aggregate_local_pass_true_all_within_zero_mismatches() {
    let errors = [
        errs(1.0, vec![1.5, -0.5], 0, 0.0003, -0.0002, 0.05, true),
        errs(-0.8, vec![0.5], 0, 0.0001, 0.0004, -0.08, true),
    ];
    let report = aggregate_local(&errors, &ToleranceProfile::standard());
    assert!(report.pass, "全条件 OK → pass==true");
    // 統計は populated（誤差を隠さない）。
    assert_eq!(report.maximum.n, 2, "maximum.n=2（統計は出す）");
    assert!(
        report.maximum.max_abs > 0.0,
        "maximum.max_abs は非ゼロ, got {}",
        report.maximum.max_abs
    );
    assert_eq!(
        report.contacts.n, 3,
        "contacts.n=3（フラット集計・統計は出す）"
    );
    assert!(
        report.contacts.max_abs > 0.0,
        "contacts.max_abs は非ゼロ, got {}",
        report.contacts.max_abs
    );
    assert_eq!(report.visibility_mismatches, 0, "visibility_mismatches=0");
    assert_eq!(
        report.contact_presence_mismatches, 0,
        "presence_mismatches=0"
    );
}

// ------------------------------------------------------------
// per-gate pass==false（各 gate を 1 件ずつ単離して撃破）
//   各テストで「対象 gate 以外はすべて pass 側」に固定する。
// ------------------------------------------------------------

/// gate「maximum」: maximum_seconds=2.0 > standard.maximum_seconds(1.5) → pass==false。
/// 他 metric/可視/presence はすべて OK。殺す変異: maximum gate の脱落。
#[test]
fn aggregate_local_pass_false_maximum_exceeds() {
    let errors = [errs(2.0, vec![0.5], 0, 0.0001, 0.0001, 0.05, true)];
    let report = aggregate_local(&errors, &ToleranceProfile::standard());
    assert!(!report.pass, "maximum 2.0 > 1.5 → pass==false");
}

/// gate「contacts」: contact_seconds に 3.0 > standard.contact_seconds(2.0) → pass==false。
/// 他はすべて OK。殺す変異: contacts gate の脱落。
#[test]
fn aggregate_local_pass_false_contacts_exceed() {
    let errors = [errs(0.5, vec![3.0], 0, 0.0001, 0.0001, 0.05, true)];
    let report = aggregate_local(&errors, &ToleranceProfile::standard());
    assert!(!report.pass, "contact 3.0 > 2.0 → pass==false");
}

/// gate「magnitude」: magnitude=0.001 > standard.magnitude(0.0005) → pass==false。
/// 他はすべて OK。殺す変異: magnitude gate の脱落。
#[test]
fn aggregate_local_pass_false_magnitude_exceeds() {
    let errors = [errs(0.5, vec![0.5], 0, 0.001, 0.0001, 0.05, true)];
    let report = aggregate_local(&errors, &ToleranceProfile::standard());
    assert!(!report.pass, "magnitude 0.001 > 0.0005 → pass==false");
}

/// gate「obscuration」: obscuration=0.001 > standard.obscuration(0.0005) → pass==false。
/// 他はすべて OK。殺す変異: obscuration gate の脱落。
#[test]
fn aggregate_local_pass_false_obscuration_exceeds() {
    let errors = [errs(0.5, vec![0.5], 0, 0.0001, 0.001, 0.05, true)];
    let report = aggregate_local(&errors, &ToleranceProfile::standard());
    assert!(!report.pass, "obscuration 0.001 > 0.0005 → pass==false");
}

/// gate「altitude」: max_altitude_deg=0.2 > standard.altitude_degrees(0.1) → pass==false。
/// 他はすべて OK。殺す変異: altitude gate の脱落。
#[test]
fn aggregate_local_pass_false_altitude_exceeds() {
    let errors = [errs(0.5, vec![0.5], 0, 0.0001, 0.0001, 0.2, true)];
    let report = aggregate_local(&errors, &ToleranceProfile::standard());
    assert!(!report.pass, "altitude 0.2 > 0.1 → pass==false");
}

/// gate「visibility_mismatches==0」: 1 件の visibility_matches=false（他 metric/presence は OK）
/// → pass==false。殺す変異: visibility gate の脱落。
#[test]
fn aggregate_local_pass_false_visibility_mismatch() {
    let errors = [errs(0.5, vec![0.5], 0, 0.0001, 0.0001, 0.05, false)];
    let report = aggregate_local(&errors, &ToleranceProfile::standard());
    assert!(!report.pass, "visibility 不一致 1 件 → pass==false");
}

/// gate「contact_presence_mismatches==0」: 1 件の presence 不一致（他 metric/可視は OK）
/// → pass==false。殺す変異: presence gate の脱落。
#[test]
fn aggregate_local_pass_false_presence_mismatch() {
    let errors = [errs(0.5, vec![0.5], 1, 0.0001, 0.0001, 0.05, true)];
    let report = aggregate_local(&errors, &ToleranceProfile::standard());
    assert!(!report.pass, "presence 不一致 1 件 → pass==false");
}

/// 受け入れ「空列は vacuous pass・5 ErrorStats すべてゼロ・件数 0」。
/// aggregate_local(&[], standard) → 5 ErrorStats n=0 かつ全 0.0、visibility_mismatches=0、
/// contact_presence_mismatches=0、pass==true。
/// 殺す変異: 空でパニック, 空で pass==false, 空で件数が 0 にならない。
#[test]
fn aggregate_local_empty_is_vacuous_pass_zero_stats() {
    let report = aggregate_local(&[], &ToleranceProfile::standard());
    assert!(report.pass, "空列は vacuous pass==true");
    assert_eq!(
        report.visibility_mismatches, 0,
        "空で visibility_mismatches=0"
    );
    assert_eq!(
        report.contact_presence_mismatches, 0,
        "空で contact_presence_mismatches=0"
    );

    let zero = |s: &ErrorStats, label: &str| {
        assert_eq!(s.n, 0, "{label}.n must be 0 for empty");
        assert!(
            s.max_abs == 0.0,
            "{label}.max_abs must be 0.0, got {}",
            s.max_abs
        );
        assert!(
            s.mean_abs == 0.0 && !s.mean_abs.is_nan(),
            "{label}.mean_abs must be 0.0 (not NaN), got {}",
            s.mean_abs
        );
        assert!(s.p95 == 0.0, "{label}.p95 must be 0.0, got {}", s.p95);
    };
    zero(&report.maximum, "maximum");
    zero(&report.contacts, "contacts");
    zero(&report.magnitude, "magnitude");
    zero(&report.obscuration, "obscuration");
    zero(&report.max_altitude, "max_altitude");
}
