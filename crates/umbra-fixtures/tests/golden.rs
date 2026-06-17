//! ISSUE-029 ゴールデンフィクスチャ受け入れテスト（strict / 構造・サニティ・被覆・出典完全性・固定性）。
//!
//! 本ファイルは `umbra-fixtures` の **公開 API のみ**を対象とした統合テスト（tests/ 配下）。
//! 検証対象は ISSUE-029 の受け入れ基準（issue §受け入れテスト, accuracy.md §3.4, conventions §3/§7/§9）。
//!
//! ## 規律（ハードコード禁止 / 数値事実の検証ではない）
//! 本テストは **構造・サニティ・被覆・出典完全性・チェックサム固定**のみを縛る。
//! 特定の日食の gamma/接触時刻などのオラクル数値は **一切ハードコードしない**
//! （実装出力とオラクル値の一致は ISSUE-030 = `ToleranceProfile` の領分。conventions §11）。
//!
//! ## seed フェーズの分割（ACTIVE / DEFERRED）
//! 本スライスは「インフラ＋実データ seed」であり、ゴールデン20 全件は後続拡張。
//! - **ACTIVE**: `golden_eclipses()` が返した内容に対し、現時点の seed で必ず成立すべき不変条件。
//! - **DEFERRED**: `#[ignore = "..."]` でゴールデン20 全件の受け入れ基準を可視化のまま残す。
//!
//! 各テスト関数の doc にカバーする受け入れ箇条を明記する。
//!
//! ## 期待される RED（実装前）
//! `umbra_fixtures` の公開 API（`GoldenEclipse` 他・`golden_eclipses()`・`fixtures_checksum()`・
//! `FIXTURES_CHECKSUM`）と `umbra-eclipse` への依存はまだ存在しないため、本ファイルは
//! **コンパイル不能（API 不在）で RED** になる。これが想定どおりの赤である。

#![allow(clippy::manual_range_contains)]

use umbra_core::{TtInstant, UtcInstant};
use umbra_eclipse::{SolarEclipseKind, Visibility};
use umbra_fixtures::{
    fixtures_checksum, golden_eclipses, GoldenContact, GoldenEclipse, GoldenLocation,
    LocationClass, OracleSource, FIXTURES_CHECKSUM,
};

// ============================================================
// 共通ヘルパ
// ============================================================

/// 文字列が trim 後に非空であること。
fn non_empty(s: &str) -> bool {
    !s.trim().is_empty()
}

/// `retrieved` が YYYY-MM-DD 形（10 文字・idx 4,7 がダッシュ・他は数字）であること。
fn is_yyyy_mm_dd(s: &str) -> bool {
    let b = s.as_bytes();
    if b.len() != 10 {
        return false;
    }
    for (i, &c) in b.iter().enumerate() {
        match i {
            4 | 7 => {
                if c != b'-' {
                    return false;
                }
            }
            _ => {
                if !c.is_ascii_digit() {
                    return false;
                }
            }
        }
    }
    true
}

/// k 慣習として許容する文書化済みの言及（conventions §9 / issue §慣習）。
/// Espenak 2 値慣習に言及するか、2 つの k 値を明示すること（いずれもケースインセンシティブ）。
fn k_convention_is_documented(k: &str) -> bool {
    let lower = k.to_lowercase();
    let mentions_espenak = lower.contains("espenak");
    // 2 値明示: umbral/penumbral の両語、または 2 つの数値風トークンを含む。
    let mentions_two_named = lower.contains("umbral") && lower.contains("penumbral");
    let mentions_iau_mean = lower.contains("iaumean") || lower.contains("iau mean");
    // 数値 2 個（"0.2722" と "0.2725" のような）を見つける緩いチェック。
    let numeric_tokens = lower
        .split(|c: char| !(c.is_ascii_digit() || c == '.'))
        .filter(|t| t.contains('.') && t.len() >= 4)
        .count();
    mentions_espenak || mentions_two_named || mentions_iau_mean || numeric_tokens >= 2
}

/// 接触の**論理順序**で present(`Some`)のみを集める: c1, c2, maximum, c3, c4。
/// maximum は常に present。
fn present_contact_times(loc: &GoldenLocation) -> Vec<UtcInstant> {
    let mut v: Vec<UtcInstant> = Vec::new();
    if let Some(c) = &loc.c1 {
        v.push(c.time_utc);
    }
    if let Some(c) = &loc.c2 {
        v.push(c.time_utc);
    }
    v.push(loc.maximum.time_utc);
    if let Some(c) = &loc.c3 {
        v.push(c.time_utc);
    }
    if let Some(c) = &loc.c4 {
        v.push(c.time_utc);
    }
    v
}

/// 全地点をフラットに走査。
fn all_locations(eclipses: &[GoldenEclipse]) -> impl Iterator<Item = &GoldenLocation> {
    eclipses.iter().flat_map(|e| e.locations.iter())
}

// ============================================================
// ACTIVE 1: ローダ非空・整合性（受け入れ: ローダ §構造検証）
// ============================================================

/// ACTIVE / 受け入れ「ローダ: golden_eclipses() が件数・地点・OracleSource 充足を返す」。
/// seed フェーズの下限: 日食 ≥4・総地点 ≥15・各日食 地点 ≥3。
#[test]
fn loader_returns_nonempty_with_minimum_seed_integrity() {
    let eclipses = golden_eclipses();
    assert!(
        eclipses.len() >= 4,
        "seed must bundle at least 4 eclipses, got {}",
        eclipses.len()
    );
    let total_locations: usize = eclipses.iter().map(|e| e.locations.len()).sum();
    assert!(
        total_locations >= 15,
        "seed must bundle at least 15 locations total, got {total_locations}"
    );
    for e in &eclipses {
        assert!(
            e.locations.len() >= 3,
            "each eclipse must have >=3 locations, '{}' has {}",
            e.event_key,
            e.locations.len()
        );
    }
}

// ============================================================
// ACTIVE 2: 食種被覆（受け入け: 被覆検証 §食種）
// ============================================================

/// ACTIVE / 受け入れ「被覆検証: 食種（皆既/金環/部分/ハイブリッド）が漏れなく含まれる」。
/// seed の kind_expected 集合 ⊇ {Total, Annular, Hybrid, Partial}。
#[test]
fn kind_coverage_includes_total_annular_hybrid_partial() {
    let eclipses = golden_eclipses();
    let kinds: Vec<SolarEclipseKind> = eclipses.iter().map(|e| e.kind_expected).collect();
    for required in [
        SolarEclipseKind::Total,
        SolarEclipseKind::Annular,
        SolarEclipseKind::Hybrid,
        SolarEclipseKind::Partial,
    ] {
        assert!(
            kinds.contains(&required),
            "kind coverage must include {required:?}; present kinds = {kinds:?}"
        );
    }
}

// ============================================================
// ACTIVE 3: 地点分類被覆（受け入れ: 被覆検証 §地点分類）
// ============================================================

/// ACTIVE / 受け入れ「被覆検証: 地点分類（中心線上〜限界〜可視域外〜日の出日没）が漏れなく含まれる」。
/// classes ⊇ {Centerline, NearLimit, PartialZone}、かつ {Sunrise, Sunset} の少なくとも一方、
/// かつ {BelowHorizon, HighElevation} の少なくとも一方。
#[test]
fn location_class_coverage_spans_core_and_edge_classes() {
    let eclipses = golden_eclipses();
    let classes: Vec<LocationClass> = all_locations(&eclipses).map(|l| l.location_class).collect();
    let has = |c: LocationClass| classes.contains(&c);

    for required in [
        LocationClass::Centerline,
        LocationClass::NearLimit,
        LocationClass::PartialZone,
    ] {
        assert!(
            has(required),
            "location class coverage must include {required:?}; present = {classes:?}"
        );
    }
    assert!(
        has(LocationClass::Sunrise) || has(LocationClass::Sunset),
        "coverage must include at least one of Sunrise/Sunset; present = {classes:?}"
    );
    assert!(
        has(LocationClass::BelowHorizon) || has(LocationClass::HighElevation),
        "coverage must include at least one of BelowHorizon/HighElevation; present = {classes:?}"
    );
}

// ============================================================
// ACTIVE 4: 出典完全性（受け入れ: 出典完全性 §data-sources §4）
// ============================================================

/// ACTIVE / 受け入れ「出典完全性: 各 OracleSource に name/url/retrieved/慣習/ライセンス注記が非空」。
/// retrieved は YYYY-MM-DD としてパース可能であること。
#[test]
fn source_completeness_all_fields_nonempty_and_retrieved_parses() {
    let eclipses = golden_eclipses();
    for e in &eclipses {
        let s: &OracleSource = &e.source;
        assert!(
            non_empty(&s.name),
            "source.name empty for '{}'",
            e.event_key
        );
        assert!(non_empty(&s.url), "source.url empty for '{}'", e.event_key);
        assert!(
            non_empty(&s.retrieved),
            "source.retrieved empty for '{}'",
            e.event_key
        );
        assert!(
            non_empty(&s.delta_t_convention),
            "source.delta_t_convention empty for '{}'",
            e.event_key
        );
        assert!(
            non_empty(&s.k_convention),
            "source.k_convention empty for '{}'",
            e.event_key
        );
        assert!(
            non_empty(&s.license_note),
            "source.license_note empty for '{}'",
            e.event_key
        );
        assert!(
            is_yyyy_mm_dd(&s.retrieved),
            "source.retrieved '{}' must be YYYY-MM-DD for '{}'",
            s.retrieved,
            e.event_key
        );
    }
}

// ============================================================
// ACTIVE 5: k 慣習の記録（受け入れ: 慣習整合 §conventions §9）
// ============================================================

/// ACTIVE / 受け入れ「慣習整合: k_convention が Espenak 2 値 or 明示」。
/// 非空 + 文書化済みの許容集合（Espenak / umbral+penumbral 2 語 / IauMean / 2 数値）に該当すること。
#[test]
fn k_convention_recorded_and_references_known_convention() {
    let eclipses = golden_eclipses();
    for e in &eclipses {
        let k = &e.source.k_convention;
        assert!(
            non_empty(k),
            "k_convention must be non-empty for '{}'",
            e.event_key
        );
        assert!(
            k_convention_is_documented(k),
            "k_convention '{}' for '{}' must reference a documented convention \
             (Espenak two-value, umbral+penumbral, IauMean, or two explicit k values)",
            k,
            e.event_key
        );
    }
}

// ============================================================
// ACTIVE 6: 日食ごとのサニティ（受け入れ: 値妥当性 §サニティ）
// ============================================================

/// ACTIVE / 受け入れ「値妥当性: gamma ∈ [-1.6,1.6]、食分 >0、皆既で食分 ≥1」。
/// 加えて金環では食分 <1（issue §サニティ・accuracy.md §3.4）。
#[test]
fn per_eclipse_sanity_gamma_and_magnitude() {
    let eclipses = golden_eclipses();
    for e in &eclipses {
        assert!(
            e.gamma >= -1.6 && e.gamma <= 1.6,
            "gamma {} out of [-1.6,1.6] for '{}'",
            e.gamma,
            e.event_key
        );
        assert!(
            e.magnitude > 0.0,
            "global magnitude {} must be > 0 for '{}'",
            e.magnitude,
            e.event_key
        );
        match e.kind_expected {
            SolarEclipseKind::Total | SolarEclipseKind::NonCentralTotal => assert!(
                e.magnitude >= 1.0,
                "Total eclipse '{}' must have magnitude >= 1.0, got {}",
                e.event_key,
                e.magnitude
            ),
            SolarEclipseKind::Annular | SolarEclipseKind::NonCentralAnnular => assert!(
                e.magnitude < 1.0,
                "Annular eclipse '{}' must have magnitude < 1.0, got {}",
                e.event_key,
                e.magnitude
            ),
            _ => {}
        }
    }
}

// ============================================================
// ACTIVE 7: 地点ごとのサニティ（受け入れ: 値妥当性 §サニティ・接触順序・規約）
// ============================================================

/// ACTIVE / 受け入れ「値妥当性（地点）」: 緯度経度標高・食分/食面積・高度方位のレンジ。
/// conventions §3（東経正・範囲）, §7（方位北0東回り [0,360)）。
#[test]
fn per_location_sanity_ranges() {
    let eclipses = golden_eclipses();
    for (e, l) in eclipses
        .iter()
        .flat_map(|e| e.locations.iter().map(move |l| (e, l)))
    {
        let id = format!("{} / {}", e.event_key, l.name);
        assert!(
            l.latitude_deg >= -90.0 && l.latitude_deg <= 90.0,
            "latitude_deg {} out of [-90,90] ({id})",
            l.latitude_deg
        );
        // 東経正・範囲 (-180,180]（conventions §3）。
        assert!(
            l.east_longitude_deg > -180.0 && l.east_longitude_deg <= 180.0,
            "east_longitude_deg {} out of (-180,180] ({id})",
            l.east_longitude_deg
        );
        assert!(
            l.elevation_m >= -500.0 && l.elevation_m <= 9000.0,
            "elevation_m {} out of [-500,9000] ({id})",
            l.elevation_m
        );
        assert!(
            l.obscuration >= 0.0 && l.obscuration <= 1.0,
            "obscuration {} out of [0,1] ({id})",
            l.obscuration
        );
        assert!(
            l.magnitude >= 0.0,
            "location magnitude {} must be >= 0 ({id})",
            l.magnitude
        );
        // 方位北0東回り [0,360)（conventions §7）。
        assert!(
            l.max_azimuth_deg >= 0.0 && l.max_azimuth_deg < 360.0,
            "max_azimuth_deg {} out of [0,360) ({id})",
            l.max_azimuth_deg
        );
        assert!(
            l.max_altitude_deg >= -90.0 && l.max_altitude_deg <= 90.0,
            "max_altitude_deg {} out of [-90,90] ({id})",
            l.max_altitude_deg
        );
    }
}

/// ACTIVE / 受け入れ「値妥当性: 接触順序 c1<max<c4」。
/// present(`Some`)の (c1,c2,maximum,c3,c4) を論理順に並べた time_utc が非減少。
/// c1/c4 が present なら c1 < maximum < c4 が厳密に成立。
#[test]
fn per_location_contact_ordering_is_monotonic() {
    let eclipses = golden_eclipses();
    for (e, l) in eclipses
        .iter()
        .flat_map(|e| e.locations.iter().map(move |l| (e, l)))
    {
        let id = format!("{} / {}", e.event_key, l.name);
        let times = present_contact_times(l);
        // 論理順で非減少（PartialOrd 由来の time_utc 比較）。
        for w in times.windows(2) {
            assert!(
                w[0] <= w[1],
                "present contact times must be non-decreasing in logical order ({id})"
            );
        }
        // c1/c4 が present なら maximum を厳密に挟む。
        if let (Some(c1), Some(c4)) = (&l.c1, &l.c4) {
            assert!(
                c1.time_utc < l.maximum.time_utc && l.maximum.time_utc < c4.time_utc,
                "c1 < maximum < c4 must hold strictly when both present ({id})"
            );
        }
    }
}

/// ACTIVE / 受け入れ「可視性整合」: visibility_expected と max_altitude_deg の一貫性。
/// BelowHorizon → max_altitude_deg < 0; FullyVisible/PartialVisible/SunriseEclipse/
/// SunsetEclipse → max_altitude_deg >= 0; NotVisible は高度チェックをスキップ。
#[test]
fn per_location_visibility_consistent_with_max_altitude() {
    let eclipses = golden_eclipses();
    for (e, l) in eclipses
        .iter()
        .flat_map(|e| e.locations.iter().map(move |l| (e, l)))
    {
        let id = format!("{} / {}", e.event_key, l.name);
        match l.visibility_expected {
            Visibility::BelowHorizon => assert!(
                l.max_altitude_deg < 0.0,
                "BelowHorizon requires max_altitude_deg < 0, got {} ({id})",
                l.max_altitude_deg
            ),
            Visibility::FullyVisible
            | Visibility::PartialVisible
            | Visibility::SunriseEclipse
            | Visibility::SunsetEclipse => assert!(
                l.max_altitude_deg >= 0.0,
                "{:?} requires max_altitude_deg >= 0, got {} ({id})",
                l.visibility_expected,
                l.max_altitude_deg
            ),
            Visibility::NotVisible => { /* altitude check skipped */ }
            // Visibility は #[non_exhaustive]: 将来 variant は高度チェックをスキップ（前方互換）。
            _ => { /* unknown future visibility: skip */ }
        }
    }
}

/// ACTIVE / 受け入れ「皆既性/金環性」: Total 日食かつ当該地点が c2/c3 ともに present(中心相)なら
/// obscuration ≈ 1.0（>= 0.999）かつ location magnitude >= 1.0。
#[test]
fn per_location_totality_implies_full_obscuration() {
    let eclipses = golden_eclipses();
    for e in &eclipses {
        if e.kind_expected != SolarEclipseKind::Total {
            continue;
        }
        for l in &e.locations {
            if l.c2.is_some() && l.c3.is_some() {
                let id = format!("{} / {}", e.event_key, l.name);
                assert!(
                    l.obscuration >= 0.999,
                    "Total central location must have obscuration ~1.0 (>=0.999), got {} ({id})",
                    l.obscuration
                );
                assert!(
                    l.magnitude >= 1.0,
                    "Total central location must have magnitude >= 1.0, got {} ({id})",
                    l.magnitude
                );
            }
        }
    }
}

/// ACTIVE / 受け入れ「出典完全性: TT 捏造禁止」(accuracy.md §0)。
/// USNO オラクルは局所接触を **UT のみ**で提供する。UT-only の接触に TT を捏造することは禁止。
/// よって全日食・全地点で、present(`Some`)な各接触 (c1,c2,maximum,c3,c4) の `time_tt` は
/// すべて `None` でなければならない。
/// (日食レベルの `greatest_time_tt` は NASA TD 由来で Some を許容するため、ここでは検査しない。)
#[test]
fn per_contact_time_tt_is_none_no_fabricated_tt() {
    let eclipses = golden_eclipses();
    for (e, l) in eclipses
        .iter()
        .flat_map(|e| e.locations.iter().map(move |l| (e, l)))
    {
        let id = format!("{} / {}", e.event_key, l.name);
        let check = |label: &str, c: &Option<GoldenContact>| {
            if let Some(c) = c {
                assert!(
                    c.time_tt.is_none(),
                    "{label} time_tt must be None (no fabricated TT for UT-only USNO contact) ({id})"
                );
            }
        };
        check("c1", &l.c1);
        check("c2", &l.c2);
        check("c3", &l.c3);
        check("c4", &l.c4);
        // maximum は常に present。
        assert!(
            l.maximum.time_tt.is_none(),
            "maximum time_tt must be None (no fabricated TT for UT-only USNO contact) ({id})"
        );
    }
}

/// ACTIVE / 受け入れ「金環性」: Annular 日食かつ当該地点が c2/c3 ともに present(金環相/中心線)なら
/// obscuration ∈ (0,1) かつ location magnitude ∈ (0,1)。
/// 金環中心の地点はリング（金環）を見るため、皆既のような全被覆には決してならない
/// (per_location_totality_implies_full_obscuration の対称ガード)。
#[test]
fn per_location_annular_phase_is_a_ring() {
    let eclipses = golden_eclipses();
    for e in &eclipses {
        if e.kind_expected != SolarEclipseKind::Annular {
            continue;
        }
        for l in &e.locations {
            if l.c2.is_some() && l.c3.is_some() {
                let id = format!("{} / {}", e.event_key, l.name);
                assert!(
                    l.obscuration < 1.0 && l.obscuration > 0.0,
                    "Annular central location must see a ring: obscuration in (0,1), got {} ({id})",
                    l.obscuration
                );
                assert!(
                    l.magnitude < 1.0 && l.magnitude > 0.0,
                    "Annular central location must see a ring: magnitude in (0,1), got {} ({id})",
                    l.magnitude
                );
            }
        }
    }
}

// ============================================================
// ACTIVE 8: チェックサム固定（受け入れ: 固定性 §変更検知）
// ============================================================

/// ACTIVE / 受け入れ「固定性: checksum 一致（変更検知）」。
/// fixtures_checksum() == FIXTURES_CHECKSUM、かつ 64 桁の小文字 16 進。
#[test]
fn checksum_matches_recorded_and_is_64_lower_hex() {
    let actual = fixtures_checksum();
    assert_eq!(
        actual, FIXTURES_CHECKSUM,
        "fixtures_checksum() must equal recorded FIXTURES_CHECKSUM (fixture drift detector)"
    );
    assert_eq!(
        actual.len(),
        64,
        "checksum must be 64 hex chars (SHA-256), got {}",
        actual.len()
    );
    assert!(
        actual
            .chars()
            .all(|c| c.is_ascii_digit() || ('a'..='f').contains(&c)),
        "checksum must be lowercase hex, got '{actual}'"
    );
}

// ============================================================
// ACTIVE 9: event_key 一意性（受け入れ: ローダ §構造検証）
// ============================================================

/// ACTIVE / 受け入れ「ローダ構造: event_key が全て非空・一意」。
#[test]
fn event_keys_are_unique_and_nonempty() {
    let eclipses = golden_eclipses();
    let mut seen: Vec<&str> = Vec::new();
    for e in &eclipses {
        assert!(non_empty(&e.event_key), "event_key must be non-empty");
        assert!(
            !seen.contains(&e.event_key.as_str()),
            "event_key '{}' is duplicated",
            e.event_key
        );
        seen.push(e.event_key.as_str());
    }
}

// ============================================================
// DEFERRED: ゴールデン20 全件の受け入れ基準（後続拡張で外す）
// ============================================================

/// DEFERRED / 受け入れ「ローダ: 20 件」。seed を golden-twenty に拡張したら有効化する。
#[test]
#[ignore = "full-set acceptance (ISSUE-029): expand seed to golden-twenty"]
fn deferred_golden_twenty_has_exactly_20_eclipses() {
    let eclipses = golden_eclipses();
    assert_eq!(
        eclipses.len(),
        20,
        "golden-twenty must contain exactly 20 eclipses, got {}",
        eclipses.len()
    );
}

/// DEFERRED / 受け入れ「各日食に地点 5〜10」。
#[test]
#[ignore = "full-set acceptance (ISSUE-029): expand seed to golden-twenty"]
fn deferred_each_eclipse_has_5_to_10_locations() {
    let eclipses = golden_eclipses();
    for e in &eclipses {
        let n = e.locations.len();
        assert!(
            (5..=10).contains(&n),
            "each eclipse must have 5..=10 locations, '{}' has {n}",
            e.event_key
        );
    }
}

/// DEFERRED / 受け入れ「accuracy.md §3.4 の食種別件数」:
/// 皆既≥5 / 金環≥5 / 部分≥3 / ハイブリッド≥2 + 境界/日の出日没/極域 地点 ≥5。
#[test]
#[ignore = "full-set acceptance (ISSUE-029): expand seed to golden-twenty"]
fn deferred_full_category_counts_match_accuracy_3_4() {
    let eclipses = golden_eclipses();
    let count_kind = |k: SolarEclipseKind| eclipses.iter().filter(|e| e.kind_expected == k).count();
    assert!(
        count_kind(SolarEclipseKind::Total) >= 5,
        "need >=5 total eclipses"
    );
    assert!(
        count_kind(SolarEclipseKind::Annular) >= 5,
        "need >=5 annular eclipses"
    );
    assert!(
        count_kind(SolarEclipseKind::Partial) >= 3,
        "need >=3 partial eclipses"
    );
    assert!(
        count_kind(SolarEclipseKind::Hybrid) >= 2,
        "need >=2 hybrid eclipses"
    );

    // 境界/日の出日没/極域に該当する地点（Sunrise/Sunset/BelowHorizon/NearLimit/HighElevation）が
    // セット全体で >=5。
    let boundary_like = all_locations(&eclipses)
        .filter(|l| {
            matches!(
                l.location_class,
                LocationClass::Sunrise
                    | LocationClass::Sunset
                    | LocationClass::BelowHorizon
                    | LocationClass::NearLimit
                    | LocationClass::HighElevation
            )
        })
        .count();
    assert!(
        boundary_like >= 5,
        "need >=5 boundary/sunrise-sunset/polar locations across the set, got {boundary_like}"
    );
}

// 未使用 import 警告を避けるための型参照アンカー（実装後は各テストが直接使う）。
const _: fn() -> Option<(GoldenContact, TtInstant)> = || None;
