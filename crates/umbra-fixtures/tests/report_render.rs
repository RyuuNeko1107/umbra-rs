//! ISSUE-030 S30e 受け入れテスト（strict / ゴールデンレポートのレンダリング）。
//!
//! 本ファイルは `umbra-fixtures` の **公開 API のみ**を対象とした統合テスト（tests/ 配下）。
//! 対象は本スライスが追加する 2 つのレンダリング関数:
//! - `render_text(&GoldenReport) -> String`（人間可読の複数行サマリ）
//! - `render_json(&GoldenReport) -> Result<String, serde_json::Error>`
//!   （serde 派生で `GoldenReport` 全体を pretty-print した機械可読 JSON＋末尾改行）
//!
//! ## テスト戦略（mutation-resistant）
//! `sample_report()` は **フィールド毎に区別可能な値**で完全に埋めた `GoldenReport` を作る。
//! これにより「どのフィールドがどの出力位置に配線されたか」を独立に縛り、フィールド取り違え
//! （例 global↔local の pass 入れ替え, max_abs を別 metric から取る）変異を殺せる。
//! JSON は実装の整形に依存しないよう `serde_json::Value` にパースして構造を検査し、テキストは
//! CLI text-format テストと同じく **部分文字列**で検査する（厳密レイアウトは固定しない）。
//!
//! ## 期待される RED（実装前）
//! `render_text` / `render_json` はまだ存在しないため、本ファイルは **未解決インポート（E0432）
//! でコンパイル不能 = RED** になる。これが想定どおりの赤である（`GoldenReport` 等の型は既存・
//! テストはそれらを直接構築するのでインポートは解決する）。

use serde_json::Value;

use umbra_fixtures::{
    render_json, render_text, ErrorStats, GlobalReport, GoldenReport, LocalReport,
};

/// JSON 数値の一致に用いる厳密許容。
const EPS: f64 = 1e-9;

// ============================================================
// 構築ヘルパ
// ============================================================

/// フィールド毎に区別可能な値で完全に埋めた `GoldenReport`。
///
/// 各 metric / count / pass を一意に判別できる値にし、JSON でもテキストでも「どのフィールドが
/// どこへ流れたか」を個別に縛れるようにする（フィールド取り違え変異を殺す）。
fn sample_report() -> GoldenReport {
    let global = GlobalReport {
        greatest: ErrorStats {
            n: 5,
            max_abs: 1.5,
            mean_abs: 0.6,
            p95: 1.3,
            units: "s",
        },
        gamma: ErrorStats {
            n: 5,
            max_abs: 0.0030,
            mean_abs: 0.0011,
            p95: 0.0025,
            units: "Re",
        },
        magnitude: ErrorStats {
            n: 5,
            max_abs: 0.00040,
            mean_abs: 0.00012,
            p95: 0.00035,
            units: "",
        },
        pass: true,
    };
    let local = LocalReport {
        maximum: ErrorStats {
            n: 23,
            max_abs: 1.9,
            mean_abs: 0.7,
            p95: 1.7,
            units: "s",
        },
        contacts: ErrorStats {
            n: 80,
            max_abs: 2.4,
            mean_abs: 0.9,
            p95: 2.2,
            units: "s",
        },
        magnitude: ErrorStats {
            n: 23,
            max_abs: 0.00030,
            mean_abs: 0.00010,
            p95: 0.00028,
            units: "",
        },
        obscuration: ErrorStats {
            n: 23,
            max_abs: 0.00020,
            mean_abs: 0.00008,
            p95: 0.00018,
            units: "",
        },
        max_altitude: ErrorStats {
            n: 23,
            max_abs: 0.07,
            mean_abs: 0.03,
            p95: 0.06,
            units: "deg",
        },
        visibility_mismatches: 2,
        contact_presence_mismatches: 1,
        pass: false,
    };
    GoldenReport {
        global,
        local,
        eclipses_found: 5,
        eclipses_missing: 1,
        locations_compared: 23,
    }
}

/// 全フィールドが 0 / 空の `GoldenReport`（非パニック検査用）。
fn zero_report() -> GoldenReport {
    let zero = |units: &'static str| ErrorStats {
        n: 0,
        max_abs: 0.0,
        mean_abs: 0.0,
        p95: 0.0,
        units,
    };
    GoldenReport {
        global: GlobalReport {
            greatest: zero("s"),
            gamma: zero("Re"),
            magnitude: zero(""),
            pass: true,
        },
        local: LocalReport {
            maximum: zero("s"),
            contacts: zero("s"),
            magnitude: zero(""),
            obscuration: zero(""),
            max_altitude: zero("deg"),
            visibility_mismatches: 0,
            contact_presence_mismatches: 0,
            pass: true,
        },
        eclipses_found: 0,
        eclipses_missing: 0,
        locations_compared: 0,
    }
}

/// `render_json` の出力を `serde_json::Value` にパースする（末尾改行は trim せず as-is で parse 可能）。
fn parse_json(report: &GoldenReport) -> Value {
    let s = render_json(report).expect("render_json は有効な JSON を返す");
    serde_json::from_str::<Value>(&s).expect("render_json の出力は有効な JSON")
}

// ============================================================
// render_json — 構造・キー・末尾改行
// ============================================================

/// 受け入れ「render_json は有効な JSON でトップレベルキーと件数を持つ」。
/// 出力をパースし、トップレベルに `global`/`local`/`eclipses_found`/`eclipses_missing`/
/// `locations_compared` が存在し、件数 5/1/23 を整数で反映する。出力は '\n' で終わる。
/// 殺す変異: キーの欠落・改名、件数を JSON に出さない・取り違える、末尾改行の欠落。
#[test]
fn render_json_is_valid_and_has_top_level_keys() {
    let report = sample_report();
    let s = render_json(&report).expect("render_json は有効な JSON を返す");
    assert!(s.ends_with('\n'), "render_json の出力は末尾に改行を持つ");

    let v: Value = serde_json::from_str(&s).expect("出力は有効な JSON");
    let obj = v.as_object().expect("トップレベルは JSON オブジェクト");
    for key in [
        "global",
        "local",
        "eclipses_found",
        "eclipses_missing",
        "locations_compared",
    ] {
        assert!(obj.contains_key(key), "トップレベルにキー {key} がある");
    }
    assert_eq!(v["eclipses_found"].as_u64(), Some(5), "eclipses_found は 5");
    assert_eq!(
        v["eclipses_missing"].as_u64(),
        Some(1),
        "eclipses_missing は 1"
    );
    assert_eq!(
        v["locations_compared"].as_u64(),
        Some(23),
        "locations_compared は 23"
    );
}

/// 受け入れ「render_json は global/local の metric 統計を正しく入れ子で出す」。
/// `global.greatest` は `{n,max_abs,mean_abs,p95,units}`。max_abs≈1.5・units=="s"・n==5。
/// gamma.units=="Re"・magnitude.units==""・local.max_altitude.units=="deg"・
/// local.maximum.max_abs≈1.9。
/// 殺す変異: metric→JSON の取り違え（max_abs を別 metric から）、units 取り違え、n の誤配線。
#[test]
fn render_json_global_and_local_stats() {
    let report = sample_report();
    let v = parse_json(&report);

    // global.greatest
    let greatest = &v["global"]["greatest"];
    assert!(
        (greatest["max_abs"].as_f64().expect("max_abs は数値") - 1.5).abs() < EPS,
        "global.greatest.max_abs は 1.5, got {greatest:?}"
    );
    assert_eq!(
        greatest["units"].as_str(),
        Some("s"),
        "global.greatest.units は \"s\""
    );
    assert_eq!(greatest["n"].as_u64(), Some(5), "global.greatest.n は 5");
    // mean_abs / p95 も値検査（max_abs=1.5 と互いに区別できる 0.6 / 1.3）。
    // stats_line / serde の mean↔p95 取り違えを撃破する。
    assert!(
        (greatest["mean_abs"].as_f64().expect("mean_abs は数値") - 0.6).abs() < EPS,
        "global.greatest.mean_abs は 0.6, got {greatest:?}"
    );
    assert!(
        (greatest["p95"].as_f64().expect("p95 は数値") - 1.3).abs() < EPS,
        "global.greatest.p95 は 1.3, got {greatest:?}"
    );

    // units の取り違えを縛る（gamma=Re, magnitude="", max_altitude=deg）。
    assert_eq!(
        v["global"]["gamma"]["units"].as_str(),
        Some("Re"),
        "global.gamma.units は \"Re\""
    );
    assert_eq!(
        v["global"]["magnitude"]["units"].as_str(),
        Some(""),
        "global.magnitude.units は空文字列"
    );
    assert_eq!(
        v["local"]["max_altitude"]["units"].as_str(),
        Some("deg"),
        "local.max_altitude.units は \"deg\""
    );

    // local.maximum.max_abs（global.greatest.max_abs=1.5 と区別できる 1.9）。
    assert!(
        (v["local"]["maximum"]["max_abs"]
            .as_f64()
            .expect("max_abs は数値")
            - 1.9)
            .abs()
            < EPS,
        "local.maximum.max_abs は 1.9（global.greatest と混ざらない）"
    );
}

/// 受け入れ「render_json は pass 真偽と mismatch カウントを出す」。
/// global.pass==true・local.pass==false（true/false の差で global↔local 入れ替えも殺す）、
/// local.visibility_mismatches==2・local.contact_presence_mismatches==1。
/// 殺す変異: pass を落とす、global/local の pass 入れ替え、mismatch カウントを落とす・取り違える。
#[test]
fn render_json_pass_and_counts() {
    let report = sample_report();
    let v = parse_json(&report);

    assert_eq!(
        v["global"]["pass"].as_bool(),
        Some(true),
        "global.pass は true"
    );
    assert_eq!(
        v["local"]["pass"].as_bool(),
        Some(false),
        "local.pass は false（global と取り違えていない）"
    );
    assert_eq!(
        v["local"]["visibility_mismatches"].as_u64(),
        Some(2),
        "local.visibility_mismatches は 2"
    );
    assert_eq!(
        v["local"]["contact_presence_mismatches"].as_u64(),
        Some(1),
        "local.contact_presence_mismatches は 1"
    );
}

/// 受け入れ「render_json は入れ子の数値を正しくシリアライズする（round-trip）」。
/// パース後、`local.contacts.max_abs`≈2.4 と `global.gamma.max_abs`≈0.0030 を再確認する。
/// 殺す変異: 入れ子 struct（contacts/gamma）のシリアライズ漏れ・取り違え。
#[test]
fn render_json_round_trips_value() {
    let report = sample_report();
    let v = parse_json(&report);

    assert!(
        (v["local"]["contacts"]["max_abs"]
            .as_f64()
            .expect("max_abs は数値")
            - 2.4)
            .abs()
            < EPS,
        "local.contacts.max_abs は 2.4"
    );
    assert!(
        (v["global"]["gamma"]["max_abs"]
            .as_f64()
            .expect("max_abs は数値")
            - 0.0030)
            .abs()
            < EPS,
        "global.gamma.max_abs は 0.0030"
    );
}

// ============================================================
// render_text — 部分文字列（厳密レイアウトは固定しない）
// ============================================================

/// 受け入れ「render_text は件数・単位・pass の可読表現を含む」。
/// 非空・複数行。判別的な件数 23（locations_compared）/ 2（visibility_mismatches）/
/// 1（contact_presence_mismatches）と、単位 "s"/"Re"/"deg" を含む。global pass(true)・
/// local pass(false) の双方が読めること（global/local の節マーカ＋ pass 真偽の両表現を要求）。
/// 殺す変異: 件数を出さない、単位を出さない、pass を出さない・global/local 片方しか出さない。
#[test]
fn render_text_contains_counts_and_pass() {
    let report = sample_report();
    let text = render_text(&report);

    assert!(!text.is_empty(), "render_text は非空文字列を返す");
    assert!(
        text.lines().count() >= 2,
        "render_text は複数行サマリ, got {text:?}"
    );

    // 判別的な件数（"5"/"1" は他の数にも現れうるので、より一意な 23/2 を主検査に使う）。
    assert!(
        text.contains("23"),
        "locations_compared=23 を含む, got {text:?}"
    );
    assert!(
        text.contains('2'),
        "visibility_mismatches=2 を含む, got {text:?}"
    );

    // 単位（取り違え防止のため 3 種すべて）。
    assert!(text.contains('s'), "単位 \"s\" を含む");
    assert!(text.contains("Re"), "単位 \"Re\" を含む");
    assert!(text.contains("deg"), "単位 \"deg\" を含む");

    // global/local 双方の節マーカと pass 真偽の可読表現。
    let lower = text.to_lowercase();
    assert!(
        lower.contains("global"),
        "global 節マーカを含む, got {text:?}"
    );
    assert!(
        lower.contains("local"),
        "local 節マーカを含む, got {text:?}"
    );
    assert!(
        text.contains("true"),
        "global pass(true) の可読表現を含む, got {text:?}"
    );
    assert!(
        text.contains("false"),
        "local pass(false) の可読表現を含む, got {text:?}"
    );
}

/// 受け入れ「render_text は各 metric の大きさを含む」。
/// global.greatest.max_abs=1.5 / local.maximum.max_abs=1.9 / contacts.max_abs=2.4 /
/// max_altitude.max_abs=0.07 を部分文字列で含む（4 桁整形でも "1.5" は "1.5000" に含まれる）。
/// 殺す変異: metric 値を出さない・別 metric の値に取り違える。
#[test]
fn render_text_contains_metric_magnitudes() {
    let report = sample_report();
    let text = render_text(&report);

    assert!(
        text.contains("1.5"),
        "global.greatest.max_abs=1.5 を含む, got {text:?}"
    );
    assert!(
        text.contains("1.9"),
        "local.maximum.max_abs=1.9 を含む, got {text:?}"
    );
    assert!(
        text.contains("2.4"),
        "local.contacts.max_abs=2.4 を含む, got {text:?}"
    );
    assert!(
        text.contains("0.07"),
        "local.max_altitude.max_abs=0.07 を含む, got {text:?}"
    );
}

/// 受け入れ「render_text は全ゼロ・空レポートでもパニックしない」。
/// 0 件・全空統計の `zero_report()` で非空文字列を返し、panic しない。
/// 殺す変異: 空/0 で除算や index で panic、空で空文字列を返す。
#[test]
fn render_text_non_panicking_on_zero_report() {
    let text = render_text(&zero_report());
    assert!(
        !text.is_empty(),
        "全ゼロレポートでも render_text は非空文字列を返す"
    );
}
