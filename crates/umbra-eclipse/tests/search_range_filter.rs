//! `EclipseEngine::search(range)` の**範囲絞り込み**統合テスト（ISSUE-050）。
//!
//! `umbra-eclipse` の**公開 API のみ**を対象とした統合テスト（tests/ 配下・別クレート境界）。
//! 対象は `EclipseEngine::search(UtcRange) -> Result<Vec<SolarEclipse>, EclipseError>`。
//!
//! ## 縛る確定仕様（ISSUE-050 §確定仕様）
//! 1. **判定時刻＝最大食 `global.greatest.time_utc`**。P1〜P4 の区間交差では判定しない
//!    （範囲に部分食だけ掛かる日食は**返らない**）。
//! 2. **閉区間 `[range.start, range.end]`**。最大食がちょうど端点の日食は**返る**。
//! 3. **候補生成（`new_moon_candidates` の ±`WINDOW_HALF_WIDTH_DAYS`=1 日の広げ）は不変**。
//!    絞り込みは組み立て後に `search` が行う ⇒ 範囲の外側 ±1 日に最大食がある日食は返らない。
//!
//! ## オラクル戦略（外部表のハードコード禁止・conventions §11）
//! NASA/USNO 等の外部数値は一切書かない。各テストはまず**エンジン自身**に 2024-04-08 皆既を
//! 探させ、その `greatest.time_utc` / `partial_begin.time_utc` / `partial_end.time_utc` を読み、
//! そこから相対に要求範囲を組み立てて件数を縛る。したがって「2024-04-08 に皆既日食がある」
//! （独立事実・kind は暦から出る）以外の外部値に依存しない。
//!
//! ## 負荷（SLOW）
//! 本ファイルは全て **SLOW**（実エンジン `standard_engine(bundled_time_data())` の実 search）。
//! `search` 呼び出し回数を抑えるため、各テストは 1 回の「発見用」探索（1 日窓）で最大食時刻を
//! 取り、その後に狭い probe を最小限だけ回す。de440s 不要（解析暦）。
//!
//! ## 期待される RED（修正前）
//! 現状 `search` は候補窓（要求範囲 ±1 日）の日食をそのまま返すため、範囲外を 0 件と要求する
//! テスト（4/7・4/9 の日単位・just-outside・部分食のみ重なる範囲）が「1 件返る」で落ちる。
//! これが想定どおりの赤。

use umbra_core::UtcInstant;
use umbra_eclipse::{standard_engine, SolarEclipse, SolarEclipseKind, StandardEngine, UtcRange};
use umbra_ephemeris::bundled_time_data;

const SECONDS_PER_DAY: f64 = 86_400.0;

// ============================================================
// ヘルパ
// ============================================================

/// UTC 瞬時を暦日時から構築。
fn utc(year: i32, month: u8, day: u8, hour: u8, minute: u8, second: f64) -> UtcInstant {
    UtcInstant::from_gregorian(year, month, day, hour, minute, second).expect("有効な UTC 日時")
}

/// `t` から `seconds` 秒ずらした UTC 瞬時（負値可）。
fn shift_seconds(t: UtcInstant, seconds: f64) -> UtcInstant {
    UtcInstant::from_jd2(t.jd2().add_days(seconds / SECONDS_PER_DAY))
}

/// UTC 範囲を組み立てる。
fn range(start: UtcInstant, end: UtcInstant) -> UtcRange {
    UtcRange { start, end }
}

/// 比較用の UTC 日数（JD）。
fn jd(t: UtcInstant) -> f64 {
    t.jd2().jd()
}

/// 2024-04-08 の皆既日食を**エンジン自身に**探させて返す（発見用の 1 日窓・SLOW）。
///
/// 外部表を引かずに最大食 UTC・P1・P4 を得るための足場。窓 `[2024-04-08, 2024-04-09]` は
/// 修正の前後いずれでもこの日食を含む（修正後は最大食が窓内、修正前は候補窓がさらに広い）。
fn find_2024_total(engine: &StandardEngine) -> SolarEclipse {
    let found = engine
        .search(range(
            utc(2024, 4, 8, 0, 0, 0.0),
            utc(2024, 4, 9, 0, 0, 0.0),
        ))
        .expect("2024-04-08 窓の search は成功する");
    found
        .into_iter()
        .find(|e| matches!(e.kind, SolarEclipseKind::Total))
        .expect("2024-04-08 の皆既日食が見つかる")
}

/// 結果の最大食 UTC を JD で列挙（メッセージ用）。
fn greatest_jds(results: &[SolarEclipse]) -> Vec<f64> {
    results
        .iter()
        .map(|e| jd(e.global.greatest.time_utc))
        .collect()
}

// ============================================================
// 1. 再現（回帰ガード）: 日単位の範囲で 4/7・4/9 は 0 件、4/8 は 1 件
// ============================================================

/// **ISSUE-050 の再現そのもの**（§背景の実測表・§確定仕様 1・2 を固定）。
/// 2024-04-08 の皆既（最大食は 4/8 の日中）に対し、日単位の要求範囲
/// `[4/7, 4/8]` と `[4/9, 4/10]` は **0 件**、`[4/8, 4/9]` は **ちょうど 1 件**で、
/// その 1 件の `event_key` 日付部と最大食 UTC 日付が 2024-04-08 であることを縛る。
///
/// 殺す実装: 候補窓（要求範囲 ±1 日）の日食をそのまま返す現行実装（＝4/7・4/9 でも 1 件返る）。
/// また「範囲を絞るが片側だけ」「絞り込みに最大食以外の時刻（P1 や P4）を使う」実装も撃破する
/// （P1/P4 は最大食の数時間前後なので 4/7・4/9 の日単位範囲には掛からないが、日付を跨ぐ実装
/// ずれがあれば件数が変わる）。
#[test]
fn search_day_ranges_around_2024_04_08_exclude_neighbouring_days() {
    let engine = standard_engine(bundled_time_data());
    let eclipse = find_2024_total(&engine);
    let (gy, gm, gd, _, _, _) = eclipse.global.greatest.time_utc.to_gregorian();
    assert_eq!(
        (gy, gm, gd),
        (2024, 4, 8),
        "発見した皆既の最大食 UTC 日付は 2024-04-08"
    );

    // 前日 `[4/7, 4/8]`: 最大食（4/8 の日中）は範囲外 ⇒ 0 件。
    let prev = engine
        .search(range(
            utc(2024, 4, 7, 0, 0, 0.0),
            utc(2024, 4, 8, 0, 0, 0.0),
        ))
        .expect("4/7 窓の search は成功する");
    assert!(
        prev.is_empty(),
        "[2024-04-07, 04-08] は 0 件（最大食は 4/8 の日中）: {} 件 greatest_jd={:?}",
        prev.len(),
        greatest_jds(&prev)
    );

    // 翌日 `[4/9, 4/10]`: 同様に 0 件。
    let next = engine
        .search(range(
            utc(2024, 4, 9, 0, 0, 0.0),
            utc(2024, 4, 10, 0, 0, 0.0),
        ))
        .expect("4/9 窓の search は成功する");
    assert!(
        next.is_empty(),
        "[2024-04-09, 04-10] は 0 件（最大食は 4/8 の日中）: {} 件 greatest_jd={:?}",
        next.len(),
        greatest_jds(&next)
    );

    // 当日 `[4/8, 4/9]`: ちょうど 1 件・同一性（event_key）まで確認。
    let same = engine
        .search(range(
            utc(2024, 4, 8, 0, 0, 0.0),
            utc(2024, 4, 9, 0, 0, 0.0),
        ))
        .expect("4/8 窓の search は成功する");
    assert_eq!(
        same.len(),
        1,
        "[2024-04-08, 04-09] はちょうど 1 件: greatest_jd={:?}",
        greatest_jds(&same)
    );
    assert_eq!(
        same[0].event_key, eclipse.event_key,
        "返る 1 件は 2024-04-08 の皆既そのもの（event_key 一致）"
    );
    assert!(
        same[0].event_key.starts_with("2024-04-08#"),
        "event_key の日付部は最大食 UTC 日付 2024-04-08: {:?}",
        same[0].event_key
    );
    assert_eq!(
        jd(same[0].global.greatest.time_utc),
        jd(eclipse.global.greatest.time_utc),
        "最大食 UTC が発見用探索と一致（同一イベント）"
    );
}

// ============================================================
// 2. 閉区間の端点（ISSUE-050 §確定仕様 2）
// ============================================================

/// **閉区間 `[start, end]`**: 最大食の瞬間ちょうどを `start` に置いた範囲でも、
/// `end` に置いた範囲でも、その日食は**返る**。
/// 最大食時刻はエンジン自身の出力から取り（外部表を使わない）、両側 2 時間の幅を与える。
///
/// 殺す実装: 開区間 `(start, end)` / 半開区間 `[start, end)` あるいは `(start, end]` での絞り込み
/// （どちらか片方の端点で 0 件に落ちる）。
#[test]
fn search_closed_interval_includes_greatest_instant_at_both_endpoints() {
    let engine = standard_engine(bundled_time_data());
    let eclipse = find_2024_total(&engine);
    let g = eclipse.global.greatest.time_utc;
    let two_hours = 2.0 * 3600.0;

    // start == 最大食ちょうど。
    let at_start = engine
        .search(range(g, shift_seconds(g, two_hours)))
        .expect("start=最大食の search は成功する");
    assert_eq!(
        at_start.len(),
        1,
        "start が最大食ちょうど ⇒ 閉区間なので返る: {} 件",
        at_start.len()
    );
    assert_eq!(
        at_start[0].event_key, eclipse.event_key,
        "返るのは当該日食（start 端点）"
    );

    // end == 最大食ちょうど。
    let at_end = engine
        .search(range(shift_seconds(g, -two_hours), g))
        .expect("end=最大食の search は成功する");
    assert_eq!(
        at_end.len(),
        1,
        "end が最大食ちょうど ⇒ 閉区間なので返る: {} 件",
        at_end.len()
    );
    assert_eq!(
        at_end[0].event_key, eclipse.event_key,
        "返るのは当該日食（end 端点）"
    );
}

// ============================================================
// 3. 端点のすぐ外側（ISSUE-050 §確定仕様 2 の裏側）
// ============================================================

/// **端点のすぐ外側は返らない**: `end` を最大食の 1 秒前に置いた範囲、`start` を最大食の
/// 1 秒後に置いた範囲は、どちらも **0 件**。
///
/// マージンに 1 秒を選んだ理由: (a) JD 表現の分解能（~1e-5 秒）より十分大きく、比較の丸めで
/// 揺れない。(b) 日食の継続時間（数時間）より遥かに小さいので、「範囲内の日食を誤って落とす」
/// 側の余地を残さない＝純粋に境界の向きだけを縛る。
/// 最大食時刻はエンジン自身の出力から取る（外部表を使わない）。
///
/// 殺す実装: 比較の向き反転・`<=`/`<` の取り違えの**逆側**（端点を跨いで 1 秒外まで拾う実装）、
/// および絞り込みを入れ忘れた現行実装（候補窓のまま返る）。
#[test]
fn search_excludes_eclipse_one_second_outside_the_range() {
    let engine = standard_engine(bundled_time_data());
    let eclipse = find_2024_total(&engine);
    let g = eclipse.global.greatest.time_utc;
    let two_hours = 2.0 * 3600.0;

    // 範囲の終端が最大食の 1 秒前 ⇒ 範囲外。
    let before = engine
        .search(range(shift_seconds(g, -two_hours), shift_seconds(g, -1.0)))
        .expect("end=最大食-1s の search は成功する");
    assert!(
        before.is_empty(),
        "end が最大食の 1 秒前 ⇒ 0 件: {} 件 greatest_jd={:?}",
        before.len(),
        greatest_jds(&before)
    );

    // 範囲の始端が最大食の 1 秒後 ⇒ 範囲外。
    let after = engine
        .search(range(shift_seconds(g, 1.0), shift_seconds(g, two_hours)))
        .expect("start=最大食+1s の search は成功する");
    assert!(
        after.is_empty(),
        "start が最大食の 1 秒後 ⇒ 0 件: {} 件 greatest_jd={:?}",
        after.len(),
        greatest_jds(&after)
    );
}

// ============================================================
// 4. 部分食だけ重なる範囲（ISSUE-050 §確定仕様 1）
// ============================================================

/// **判定は最大食時刻であり、P1〜P4 の区間交差ではない**（§確定仕様 1）。
/// エンジン自身が返した P1・最大食・P4 から範囲を組み立て、
/// (a) `[P1-10分, 最大食-1分]`（部分食の前半だけが範囲に掛かる）
/// (b) `[最大食+1分, P4+10分]`（後半だけが掛かる）
/// のいずれでも当該日食は**返らない**ことを縛る。
///
/// 殺す実装: 絞り込み条件を「`[P1,P4]` と `[start,end]` の区間交差」にした実装
/// （このテストでは 1 件返って落ちる）。§確定仕様 1 が退けた案そのもの。
/// 範囲が本当に部分食に重なっていること（P1 < 最大食 < P4 で、範囲が P1〜P4 に食い込む）も
/// 独立に assert する＝「たまたま重なっていなかった」偽の合格を防ぐ。
#[test]
fn search_excludes_eclipse_whose_partial_phase_overlaps_but_greatest_is_outside() {
    let engine = standard_engine(bundled_time_data());
    let eclipse = find_2024_total(&engine);
    let g = eclipse.global.greatest.time_utc;
    let p1 = eclipse
        .global
        .partial_begin
        .expect("皆既なので P1=Some")
        .time_utc;
    let p4 = eclipse
        .global
        .partial_end
        .expect("皆既なので P4=Some")
        .time_utc;
    assert!(
        jd(p1) < jd(g) && jd(g) < jd(p4),
        "P1 < 最大食 < P4: p1={} g={} p4={}",
        jd(p1),
        jd(g),
        jd(p4)
    );

    let ten_min = 600.0;
    let one_min = 60.0;

    // (a) 前半だけ掛かる範囲: [P1-10分, 最大食-1分]。
    let leading = range(shift_seconds(p1, -ten_min), shift_seconds(g, -one_min));
    assert!(
        jd(leading.end) > jd(p1) && jd(leading.end) < jd(g),
        "範囲の終端は P1 より後・最大食より前（部分食に確かに重なる）"
    );
    let leading_results = engine
        .search(leading)
        .expect("部分食前半に重なる範囲の search は成功する");
    assert!(
        leading_results.is_empty(),
        "最大食が範囲外なら P1〜P4 が重なっても返らない（前半側）: {} 件 greatest_jd={:?}",
        leading_results.len(),
        greatest_jds(&leading_results)
    );

    // (b) 後半だけ掛かる範囲: [最大食+1分, P4+10分]。
    let trailing = range(shift_seconds(g, one_min), shift_seconds(p4, ten_min));
    assert!(
        jd(trailing.start) > jd(g) && jd(trailing.start) < jd(p4),
        "範囲の始端は最大食より後・P4 より前（部分食に確かに重なる）"
    );
    let trailing_results = engine
        .search(trailing)
        .expect("部分食後半に重なる範囲の search は成功する");
    assert!(
        trailing_results.is_empty(),
        "最大食が範囲外なら P1〜P4 が重なっても返らない（後半側）: {} 件 greatest_jd={:?}",
        trailing_results.len(),
        greatest_jds(&trailing_results)
    );
}

// ============================================================
// 5. 広範囲の不変条件（絞り込みの過剰を検出）
// ============================================================

/// **広い範囲でも複数件が昇順で返り、全件の最大食が範囲内**（一般不変条件）。
/// 2023-01-01 〜 2024-12-31 の 2 年窓（日食は年 2〜4 回あるので複数件）。件数の外部表固定は
/// しない（≥4 件の下限のみ）。昇順（最大食 UTC の単調増加）と、
/// `start <= greatest <= end` を全件について縛る。
///
/// 殺す実装: 絞り込みが過剰（例: 常に空を返す・先頭/末尾を落とす・範囲を狭く取り違える）、
/// 絞り込み時に並び順を壊す（HashSet 経由など）、絞り込みを入れ忘れる（範囲外が混じる）。
#[test]
fn search_wide_range_returns_multiple_eclipses_in_range_and_in_order() {
    let engine = standard_engine(bundled_time_data());
    let r = range(utc(2023, 1, 1, 0, 0, 0.0), utc(2024, 12, 31, 23, 59, 59.0));
    let results = engine.search(r).expect("2 年窓の search は成功する");

    assert!(
        results.len() >= 4,
        "2 年窓には日食が複数（≥4 件）: {} 件",
        results.len()
    );

    let mut prev: Option<f64> = None;
    for e in &results {
        let g = jd(e.global.greatest.time_utc);
        assert!(
            jd(r.start) <= g && g <= jd(r.end),
            "全件の最大食が要求範囲内: event_key={:?} greatest_jd={} range=[{}, {}]",
            e.event_key,
            g,
            jd(r.start),
            jd(r.end)
        );
        if let Some(p) = prev {
            assert!(
                p < g,
                "最大食 UTC の昇順（狭義単調増加）: prev={p} cur={g} event_key={:?}",
                e.event_key
            );
        }
        prev = Some(g);
    }
}
