//! M9.3 経路生成スライス（南北限界線・**厳密な錐接線解**）の統合テスト（確定 path）。
//!
//! `umbra-eclipse` の**公開 API のみ**を対象とした統合テスト（tests/ 配下・別クレート境界）。
//! 対象は `EclipseEngine::path(&SolarEclipse, PathOptions) -> Result<EclipsePath, EclipseError>`。
//! M9.1（中心線）・M9.2（GeoJSON）は実装済み。本スライスは中心食の本影帯の北/南縁
//! （`northern_limit` / `southern_limit`）を **厳密に**（自己整合ζ＋相対速度包絡）生成する
//! （ISSUE-045 残(5)・現状の geometric 近似を置換）。
//!
//! ## 確定セマンティクス（テストで縛る）
//! 1. `greatest_point`・`center_line`・`samples`・`partial_limit` は M9.1/現状どおり
//!    （center_line は中心食で Some・samples 空・partial_limit None）。
//! 2. **中心食（central_begin/central_end 両方 Some）かつ `include_limits == true`**:
//!    `northern_limit = Some(GeoLine)` かつ `southern_limit = Some(GeoLine)`。中心線と同じサンプル
//!    時刻列で、各時刻に本影帯の北縁/南縁の地表点を結ぶ。**高緯度側が北限・低緯度側が南限**。
//! 3. **`include_limits == false`**: 限界線は両方 None（center_line は include_limits に依らず
//!    中心食なら Some）。
//! 4. **非中心**（central_begin か central_end が None）: 限界線も center_line も None。
//! 5. 限界線の GeoJSON 出力は本スライス（M9.3）では扱わない。`to_geojson` への限界線 feature 化は
//!    M9.5（`path_geojson.rs` で検証）。
//!
//! ## 厳密性のオラクル（追認回避）
//! 各限界点 P を **検証済み前方射影** `project_observer_to_fundamental`（ISSUE-024・公開）へ通して
//! 自身の基本面座標 (ξ,η,ζ) を独立復元し、次の 2 条件を絶対値で表明する:
//!   条件1（錐exact・自己整合ζ）: `hypot(ξ−x, η−y) == |l2 − ζ·tan_f2|`（ζ は点自身の値）。
//!   条件2（包絡）: `(ξ−x)·rel_vx + (η−y)·rel_vy == 0`,
//!     `rel_vx = x' − μ'·(ζ·cos d − η·sin d)`, `rel_vy = y' − μ'·ξ·sin d`。
//! 真値 x,y,d,μ,l2,tan_f2 は `BesselianSource::at(t)`（公開・検証済）、微分 x',y',μ' は
//! `Polynomial::derivative()`（公開・検証済）。被テスト関数の戻りを期待値生成に流用しない。
//! **μ' を非零**にした合成中心食を使い、影速度のみに垂直な geometric 近似が条件2 を満たせないことで
//! 厳密化前の RED を保証する。
//!
//! ## テスト戦略（strict / mutation-resistant / 負荷配分）
//! - FAST（実 search 非実走・合成 SolarEclipse）: 厳密 2 条件（μ'≠0 合成中心食）・南北分離・
//!   高緯度側=北限・include_limits=false / 非中心は None。
//! - SLOW（実エンジン 1 件・2024-04-08 皆既）: search → path()。限界線 Some・中心線が南北の間・
//!   最大食付近で厳密 2 条件・帯幅が NASA 公表 ~197 km の妥当域。
//!
//! ## 期待される RED（実装前）
//! 現状 path() は `northern_limit = None` / `southern_limit = None`（M9.3 近似が未マージ、または
//! 近似が μ' を無視）を返すため、厳密 2 条件・帯幅域 assert が落ちる。コンパイルは通る。

use umbra_core::ellipsoid::{observer_geocentric, Ellipsoid};
use umbra_core::{JulianDate2, Radians, TimeInterval, TtInstant, UtcInstant};
use umbra_eclipse::{
    project_observer_to_fundamental, standard_engine, AccuracyProfile, BesselFitError,
    BesselianPolynomial, BesselianSource, CalculationMetadata, EclipseMagnitude,
    GlobalCircumstances, GlobalContact, GreatestEclipse, InstantaneousBesselianElements,
    Obscuration, ObserverFundamental, PathOptions, Polynomial, SolarEclipse, SolarEclipseKind,
};
use umbra_ephemeris::bundled_time_data;

// ============================================================
// 時刻 / 地理ヘルパ（path_center_line.rs ミラー）
// ============================================================

/// TT 時刻を 2 要素 JD から構築。
fn tt(jd1: f64, jd2: f64) -> TtInstant {
    TtInstant::from_jd2(JulianDate2::new(jd1, jd2))
}

/// epoch から経過時間 `hours`[hour] だけ進んだ TT 時刻。
fn tt_at_hours(epoch: TtInstant, hours: f64) -> TtInstant {
    TtInstant::from_jd2(epoch.jd2().add_days(hours / 24.0))
}

/// UTC 瞬時（合成日食の時刻ラベル用・幾何には無関係）。
fn utc(year: i32, month: u8, day: u8, hour: u8, minute: u8, second: f64) -> UtcInstant {
    UtcInstant::from_gregorian(year, month, day, hour, minute, second).expect("有効な UTC 日時")
}

/// 地表点（度）。
fn geo(lat: f64, lon: f64) -> umbra_geo::GeoPoint {
    umbra_geo::GeoPoint::from_degrees(lat, lon).expect("有効な地表点")
}

/// 緯度（度）取り出し。
fn lat_deg(p: &umbra_geo::GeoPoint) -> f64 {
    p.lat.degrees().0
}

/// 経度（度）取り出し。
fn lon_deg(p: &umbra_geo::GeoPoint) -> f64 {
    p.lon.degrees().0
}

/// 緯度・経度が妥当な範囲か（限界線各点のサニティ）。
fn lat_lon_in_range(p: &umbra_geo::GeoPoint) -> bool {
    let lat = lat_deg(p);
    let lon = lon_deg(p);
    (-90.0..=90.0).contains(&lat) && (-180.0..=180.0).contains(&lon)
}

/// 2 点間の概算大円距離 [km]（haversine・WGS84 平均半径 6371 km）。
fn great_circle_km(a: &umbra_geo::GeoPoint, b: &umbra_geo::GeoPoint) -> f64 {
    let r = 6371.0_f64;
    let (lat1, lon1) = (lat_deg(a).to_radians(), lon_deg(a).to_radians());
    let (lat2, lon2) = (lat_deg(b).to_radians(), lon_deg(b).to_radians());
    let dlat = lat2 - lat1;
    let dlon = lon2 - lon1;
    let h = (dlat / 2.0).sin().powi(2) + lat1.cos() * lat2.cos() * (dlon / 2.0).sin().powi(2);
    2.0 * r * h.sqrt().asin()
}

// ============================================================
// 合成 BesselianPolynomial / SolarEclipse（path_center_line.rs ミラー）
// ============================================================

/// 合成日食の epoch（窓の中心 TT）。J2000 を借用（解析暦の評価可能域内・幾何のみ使用）。
fn synth_epoch() -> TtInstant {
    tt(2_451_545.0, 0.0)
}

/// 中心食を模す合成 bessel poly（影軸が地表に当たる小さい x,y）。
/// path_center_line.rs の central_bessel と同値。
fn central_bessel() -> BesselianPolynomial {
    let epoch = synth_epoch();
    let c = |v: f64| Polynomial {
        coefficients: vec![v],
    };
    BesselianPolynomial {
        epoch_tt: epoch,
        x: Polynomial {
            coefficients: vec![0.0, 0.05],
        },
        y: c(0.10),
        d: c(0.20),
        mu: c(1.2),
        l1: c(0.54),
        l2: c(-0.009),
        tan_f1: 0.004_65,
        tan_f2: 0.004_63,
        fit_interval: TimeInterval {
            start: tt_at_hours(epoch, -2.0),
            end: tt_at_hours(epoch, 2.0),
        },
        fit_error: BesselFitError {
            max_x: 1.0e-7,
            max_y: 1.0e-7,
            max_l1: 1.0e-7,
            max_l2: 1.0e-7,
        },
    }
}

/// 厳密限界線テスト用の合成 bessel poly。**μ を一次（μ'≠0）**にして相対速度包絡を幾何近似から分離する。
/// l2<0（皆既）・gamma≪1（中心・両縁とも地表に当たる）。
///
/// **x を二次**（x''≠0）にして速度 x'(t)=0.45+2·0.02·t_hours を **t_hours 依存**にする。これにより
/// 実装の `t_hours = days_since*24` のスケール（`*24`→`+24`/`/24`）を変える変異が x_deriv.eval(t_hours)
/// を変え、相対速度包絡⊥（条件2）テストが拾える（x が一次だと x' が定数で t_hours に依らず変異が生存する）。
/// 二次係数 0.02 は微小で、epoch±1h・60s 刻み（t_hours∈[-1,1]）でも x∈[-0.40,0.52]・gamma≪1 を維持し
/// 中心軸・南北両縁が全サンプルで地表に当たる（既存の非空・lockstep 同点数 assert を壊さない）。
fn rigorous_bessel() -> BesselianPolynomial {
    let epoch = synth_epoch();
    let p = |coeffs: Vec<f64>| Polynomial {
        coefficients: coeffs,
    };
    BesselianPolynomial {
        epoch_tt: epoch,
        // x(t)=0.05 + 0.45 t + 0.02 t² — x'(t)=0.45+0.04 t_hours（東進・t_hours 依存で変異を露出）。
        x: p(vec![0.05, 0.45, 0.02]),
        y: p(vec![0.02, 0.06]),
        d: p(vec![0.20]),
        // μ'=0.26 rad/hour ≠ 0（地球自転）。rel 速度に μ' が効く。
        mu: p(vec![1.2, 0.26]),
        l1: p(vec![0.54]),
        l2: p(vec![-0.009]),
        tan_f1: 0.004_65,
        tan_f2: 0.004_63,
        fit_interval: TimeInterval {
            start: tt_at_hours(epoch, -2.0),
            end: tt_at_hours(epoch, 2.0),
        },
        fit_error: BesselFitError {
            max_x: 1.0e-7,
            max_y: 1.0e-7,
            max_l1: 1.0e-7,
            max_l2: 1.0e-7,
        },
    }
}

/// M9 残(3) 3c-iii: **limb bulge（terminator 連結）**を発火させる合成 bessel。
///
/// `rigorous_bessel` と同形（μ'≠0・x 二次で t_hours 依存）だが、影軸を**地球の縁寄り**に置く
/// （x ベース 1.15・y ベース −0.08）ことで、[P1,P4]＝epoch±1.5h の**時間端で半影縁が昼面（ζ>0）を
/// 外れる**＝`solve_limit_edge(l1)` が `Ok(None)`（`RootNotBracketed`／未収束）になる端区間を作る。
/// その端では昼面包絡が無く terminator 交点（円∩terminator 楕円・ζ=0）で連結すべき領域になる。
///
/// 中央付近のサンプルでは南北とも昼面包絡が解ける（v1 リボンの種＝`partial_limit=Some`・非退化を保つ）一方、
/// 端区間は terminator 連結が必要。**v1（連結未実装）の外環は昼面包絡頂点のみ＝ζ≈0 terminator 頂点を含まない**
/// → 3c-iii の FAST red の発火源。連結実装後は端が terminator まで張り出し ζ≈0 頂点が現れる。
///
/// 値は `rigorous_bessel` と非対称（x/y/μ/l1/l2/tan_f が別値）で取り違え変異を撃つ。影軸位置オフセット
/// （x 定数 1.15・y 定数 −0.08）以外は `rigorous_bessel` を踏襲し、fit_interval（epoch±2h）に [P1,P4]±1.5h を収める。
fn limb_continuation_bessel() -> BesselianPolynomial {
    let epoch = synth_epoch();
    let p = |coeffs: Vec<f64>| Polynomial {
        coefficients: coeffs,
    };
    BesselianPolynomial {
        epoch_tt: epoch,
        // x(t)=1.15 + 0.45 t + 0.02 t²（影軸を縁寄りに・東進・t_hours 依存で変異露出）。
        x: p(vec![1.15, 0.45, 0.02]),
        // y(t)=−0.08 + 0.06 t（縁寄り・端で半影縁が昼面を外す）。
        y: p(vec![-0.08, 0.06]),
        d: p(vec![0.20]),
        mu: p(vec![1.2, 0.26]),
        l1: p(vec![0.54]),
        l2: p(vec![-0.009]),
        tan_f1: 0.004_65,
        tan_f2: 0.004_63,
        fit_interval: TimeInterval {
            start: tt_at_hours(epoch, -2.0),
            end: tt_at_hours(epoch, 2.0),
        },
        fit_error: BesselFitError {
            max_x: 1.0e-7,
            max_y: 1.0e-7,
            max_l1: 1.0e-7,
            max_l2: 1.0e-7,
        },
    }
}

/// 与えた bessel で中心食 SolarEclipse を構築する（central_begin/end=Some, ±span_hours）。
fn central_eclipse_with_bessel(bessel: BesselianPolynomial, span_hours: f64) -> SolarEclipse {
    let epoch = bessel.epoch_tt;
    let u1 = contact(tt_at_hours(epoch, -span_hours));
    let u4 = contact(tt_at_hours(epoch, span_hours));
    let global = GlobalCircumstances {
        kind: SolarEclipseKind::Total,
        partial_begin: None,
        central_begin: Some(u1),
        greatest: greatest_at(geo(0.0, 0.0)),
        central_end: Some(u4),
        partial_end: None,
        gamma: 0.05,
    };
    SolarEclipse {
        event_key: "synthetic-rigorous#0".to_string(),
        kind: SolarEclipseKind::Total,
        global,
        bessel,
        metadata: metadata(),
    }
}

/// 限界点 P を **検証済み前方射影**へ通し、自身の (ξ,η,ζ) を返す独立オラクル
/// （axis_intercept.rs `assert_forward_roundtrip` と同パターン・逆射影の内部式は再実装しない）。
fn forward_project(
    p: &umbra_geo::GeoPoint,
    e: &InstantaneousBesselianElements,
) -> ObserverFundamental {
    let phi = p.lat.radians().0;
    let lam = p.lon.radians().0;
    let obs = observer_geocentric(&Ellipsoid::WGS84, phi, 0.0);
    project_observer_to_fundamental(&obs, Radians::new(lam), e)
}

/// path() の lockstep サンプル時刻列を独立再構成（始点・終点を必ず含み終点へクランプ）。
fn lockstep_sample_times(
    start: TtInstant,
    end: TtInstant,
    interval_seconds: f64,
) -> Vec<TtInstant> {
    let span = end.jd2().days_since(start.jd2()) * 86_400.0;
    let mut out = Vec::new();
    let mut t_sec = 0.0_f64;
    loop {
        out.push(TtInstant::from_jd2(start.jd2().add_days(t_sec / 86_400.0)));
        if t_sec >= span || interval_seconds <= 0.0 {
            break;
        }
        t_sec = (t_sec + interval_seconds).min(span);
    }
    out
}

/// 各限界点で厳密 2 条件（自己整合ζの錐exact＋相対速度包絡⊥）を表明する共有チェック。
/// `cone_tol`/`dot_tol` は呼び側が（合成＝厳密 / 実日食＝中心線位置律速で緩め）で与える。
fn assert_exact_limit_conditions(
    north: &umbra_geo::GeoLine,
    south: &umbra_geo::GeoLine,
    bessel: &BesselianPolynomial,
    sample_times: &[TtInstant],
    cone_tol: f64,
    dot_tol: f64,
) {
    let x_deriv = bessel.x.derivative();
    let y_deriv = bessel.y.derivative();
    let mu_deriv = bessel.mu.derivative();
    let epoch = bessel.epoch_tt;
    assert_eq!(north.points.len(), south.points.len(), "北縁・南縁は同点数");
    assert!(!north.points.is_empty(), "限界線は非空");
    for (i, t) in sample_times.iter().enumerate() {
        if i >= north.points.len() {
            break;
        }
        let e = bessel.at(*t).expect("区間内サンプルは評価成功");
        let t_hours = t.jd2().days_since(epoch.jd2()) * 24.0;
        let vx = x_deriv.eval(t_hours);
        let vy = y_deriv.eval(t_hours);
        let mu_rate = mu_deriv.eval(t_hours);
        let (sin_d, cos_d) = e.declination.0.sin_cos();
        for p in [&north.points[i], &south.points[i]] {
            let of = forward_project(p, &e);
            let off_x = of.xi - e.x;
            let off_y = of.eta - e.y;
            // 条件1: 面内距離 = ζ補正本影半径（自己整合ζ）。
            let in_plane = off_x.hypot(off_y);
            let umbral = (e.l2 - of.zeta * e.tan_f2).abs();
            assert!(
                (in_plane - umbral).abs() < cone_tol,
                "サンプル{i}: 面内距離 {in_plane} = |L2'| {umbral}（自己整合ζ={}）でない",
                of.zeta
            );
            // 条件2: offset ⊥ rel 速度（μ' 項込み）。
            let rel_vx = vx - mu_rate * (of.zeta * cos_d - of.eta * sin_d);
            let rel_vy = vy - mu_rate * of.xi * sin_d;
            let dot = off_x * rel_vx + off_y * rel_vy;
            assert!(
                dot.abs() < dot_tol,
                "サンプル{i}: offset·rel = {dot}（≈0 でない＝包絡条件違反）"
            );
        }
    }
}

/// 実エンジン `path()` は影軸が地表を外す/掠めるサンプル（grazing＝二重根や RootNotBracketed）を
/// **スキップ**するため、`north.points[i]` が `lockstep_sample_times` の `times[i]` に対応するとは
/// 限らない（実 2024 では U1 付近で先頭サンプルがスキップされ index↔時刻が 1〜数サンプルずれる）。
///
/// そこで各 kept index i について、**実サンプル時刻 t_i を中心線点から復元**してから厳密 2 条件を
/// 検証する。経路は東進で x が単調なので `bessel.x.eval(t_hours) == ξ_C` の根は `[U1,U4]` 区間で一意。
///
/// オラクル独立性（strict）: t_i の復元には path() 出力（中心線点）を使うが、これは「その点に対応する
/// サンプル時刻を引く」ためだけであり、cone/envelope の**期待値は bessel 多項式（path とは独立な入力）
/// から組む**。被テスト関数 path() の戻りを期待値生成に流用しないので追認にはならない。
///
/// `cone_tol`/`dot_tol` は実日食＝中心線位置律速ゆえ合成より緩く呼び側が与える。
#[allow(clippy::too_many_arguments)]
fn assert_exact_limit_conditions_real(
    center: &umbra_geo::GeoLine,
    north: &umbra_geo::GeoLine,
    south: &umbra_geo::GeoLine,
    bessel: &BesselianPolynomial,
    u1: TtInstant,
    u4: TtInstant,
    cone_tol: f64,
    dot_tol: f64,
) {
    let x_deriv = bessel.x.derivative();
    let y_deriv = bessel.y.derivative();
    let mu_deriv = bessel.mu.derivative();
    let epoch = bessel.epoch_tt;
    assert_eq!(north.points.len(), south.points.len(), "北縁・南縁は同点数");
    assert_eq!(
        north.points.len(),
        center.points.len(),
        "限界線と中心線は同点数"
    );
    assert!(!north.points.is_empty(), "限界線は非空");

    // [U1,U4] を hours で表す（x 単調の根を二分法で挟む区間）。
    let t0_hours = u1.jd2().days_since(epoch.jd2()) * 24.0;
    let t1_hours = u4.jd2().days_since(epoch.jd2()) * 24.0;

    // 時刻復元の根条件: 中心線点 P_c は時刻 t_i の影軸が地表に当たる点なので、その時刻の瞬時要素で
    // 前方射影すると基本面で (ξ,η) = (x(t_i), y(t_i)) に閉じる（P_c が軸の足＝gamma≈0）。
    //
    // 単一変数で挟むには **ξ の自己整合残差** `g(t) = forward_project(P_c, e(t)).ξ − x(t)` を使う。
    // 射影 ξ = ρcosφ′·sin(μ(t)+λ) は μ'>0（地球自転・実 2024 で +0.25 rad/hr 級）で h が単調増加、
    // x(t) は緩慢な東進。よって g は区間内で単調＝根は一意で二分法が挟める（x 単調前提のコメント）。
    for i in 0..center.points.len() {
        let p_c = &center.points[i];
        // g(t_hours) = 射影ξ(P_c, e(t)) − x(t)。t_i で 0。
        let g = |th: f64| -> f64 {
            let t = TtInstant::from_jd2(epoch.jd2().add_days(th / 24.0));
            let e = bessel.at(t).expect("区間内サンプルは評価成功");
            forward_project(p_c, &e).xi - bessel.x.eval(th)
        };

        let (mut a, mut b) = (t0_hours, t1_hours);
        let (mut ga, gb) = (g(a), g(b));
        // 端で根を挟めない（グレージング端でわずかに外れる）場合は端へクランプ。
        let t_i_hours = if ga.signum() == gb.signum() {
            if ga.abs() <= gb.abs() {
                a
            } else {
                b
            }
        } else {
            // g(a),g(b) が異符号＝根を挟む。80 反復の二分法（十分収束）。
            for _ in 0..80 {
                let m = 0.5 * (a + b);
                let gm = g(m);
                if ga.signum() == gm.signum() {
                    a = m;
                    ga = gm;
                } else {
                    b = m;
                }
            }
            0.5 * (a + b)
        };

        // 復元時刻で瞬時要素を構成し、軸 (ξ_C, η_C) は復元時刻の前方射影で一貫させる。
        let t_i = TtInstant::from_jd2(epoch.jd2().add_days(t_i_hours / 24.0));
        let e = bessel.at(t_i).expect("復元サンプル時刻は区間内");
        let center_of = forward_project(p_c, &e);
        let xi_c = center_of.xi;
        let eta_c = center_of.eta;
        // 健全性: 復元時刻で射影した (ξ_C, η_C) が (x(t_i), y(t_i)) に整合する＝時刻復元成功。
        assert!(
            (xi_c - e.x).abs() < 1e-4 && (eta_c - e.y).abs() < 1e-4,
            "サンプル{i}: 復元 (ξ_C {xi_c}, η_C {eta_c}) が (x {}, y {}) に整合しない（時刻復元失敗）",
            e.x,
            e.y
        );

        let vx = x_deriv.eval(t_i_hours);
        let vy = y_deriv.eval(t_i_hours);
        let mu_rate = mu_deriv.eval(t_i_hours);
        let (sin_d, cos_d) = e.declination.0.sin_cos();

        for p in [&north.points[i], &south.points[i]] {
            let of = forward_project(p, &e);
            let off_x = of.xi - xi_c;
            let off_y = of.eta - eta_c;
            // 条件1: 面内距離 = ζ補正本影半径（自己整合ζ）。期待値は bessel 多項式から独立に組む。
            let in_plane = off_x.hypot(off_y);
            let umbral = (e.l2 - of.zeta * e.tan_f2).abs();
            assert!(
                (in_plane - umbral).abs() < cone_tol,
                "サンプル{i}: 面内距離 {in_plane} = |L2'| {umbral}（自己整合ζ={}）でない",
                of.zeta
            );
            // 条件2: offset ⊥ rel 速度（μ' 項込み）。
            let rel_vx = vx - mu_rate * (of.zeta * cos_d - of.eta * sin_d);
            let rel_vy = vy - mu_rate * of.xi * sin_d;
            let dot = off_x * rel_vx + off_y * rel_vy;
            assert!(
                dot.abs() < dot_tol,
                "サンプル{i}: offset·rel = {dot}（≈0 でない＝包絡条件違反）"
            );
        }
    }
}

/// 代表メタデータ（幾何に無関係）。
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

/// U1/U4 の全球接触点（時刻のみ幾何に効く）。
fn contact(time_tt: TtInstant) -> GlobalContact {
    GlobalContact {
        time_utc: utc(2024, 4, 8, 18, 0, 0.0),
        time_tt,
        position: geo(20.0, -100.0),
    }
}

/// 中心食の最大食。
fn greatest_at(position: umbra_geo::GeoPoint) -> GreatestEclipse {
    GreatestEclipse {
        time_utc: utc(2024, 4, 8, 18, 17, 0.0),
        time_tt: synth_epoch(),
        position,
        magnitude: EclipseMagnitude(1.05),
        obscuration: Obscuration(1.0),
        path_width: Some(umbra_core::Kilometers(180.0)),
        central_duration: Some(200.0),
        sun_altitude: umbra_core::Degrees(70.0),
    }
}

/// 合成「中心食」SolarEclipse（central_begin/central_end=Some・kind=Total）。
fn central_eclipse(greatest_position: umbra_geo::GeoPoint) -> SolarEclipse {
    let epoch = synth_epoch();
    let u1 = contact(tt_at_hours(epoch, -1.0));
    let u4 = contact(tt_at_hours(epoch, 1.0));
    let global = GlobalCircumstances {
        kind: SolarEclipseKind::Total,
        partial_begin: None,
        central_begin: Some(u1),
        greatest: greatest_at(greatest_position),
        central_end: Some(u4),
        partial_end: None,
        gamma: 0.10,
    };
    SolarEclipse {
        event_key: "synthetic-central#0".to_string(),
        kind: SolarEclipseKind::Total,
        global,
        bessel: central_bessel(),
        metadata: metadata(),
    }
}

/// 合成「非中心/部分食」SolarEclipse（central_begin か central_end が None）。
fn noncentral_eclipse(
    greatest_position: umbra_geo::GeoPoint,
    with_begin: bool,
    with_end: bool,
) -> SolarEclipse {
    let epoch = synth_epoch();
    let central_begin = with_begin.then(|| contact(tt_at_hours(epoch, -1.0)));
    let central_end = with_end.then(|| contact(tt_at_hours(epoch, 1.0)));
    let global = GlobalCircumstances {
        kind: SolarEclipseKind::Partial,
        partial_begin: Some(contact(tt_at_hours(epoch, -1.5))),
        central_begin,
        greatest: greatest_at(greatest_position),
        central_end,
        partial_end: Some(contact(tt_at_hours(epoch, 1.5))),
        gamma: 0.10,
    };
    SolarEclipse {
        event_key: "synthetic-noncentral#0".to_string(),
        kind: SolarEclipseKind::Partial,
        global,
        bessel: central_bessel(),
        metadata: metadata(),
    }
}

// ============================================================
// FAST: 中心食＋include_limits=true → 北/南限界線 Some・非空・各点妥当
// ============================================================

/// FAST / 新規: 中心食＋include_limits 既定(true) で northern/southern_limit=Some・点列非空（≥2）・
/// 各点が妥当な緯度経度。center_line も Some（限界線と取り違えていない＝3 本とも独立に存在）。
///
/// 殺す変異: 限界線を常に None にする・空点列を返す・限界線を center_line の別名にして 1 本しか作らない・
///   include_limits を無視して常に None にする。
#[test]
fn central_eclipse_with_limits_produces_nonempty_north_and_south() {
    let engine = standard_engine(bundled_time_data());
    let eclipse = central_eclipse(geo(12.5, 77.5));

    let path = engine
        .path(&eclipse, PathOptions::default())
        .expect("中心食の path() は成功する");

    let north = path
        .northern_limit
        .as_ref()
        .expect("中心食＋include_limits=true では northern_limit=Some");
    let south = path
        .southern_limit
        .as_ref()
        .expect("中心食＋include_limits=true では southern_limit=Some");
    assert!(
        path.center_line.is_some(),
        "中心食では center_line も Some（限界線と独立）"
    );

    assert!(
        north.points.len() >= 2,
        "北限界線は ≥2 点, got {}",
        north.points.len()
    );
    assert!(
        south.points.len() >= 2,
        "南限界線は ≥2 点, got {}",
        south.points.len()
    );

    for p in north.points.iter().chain(south.points.iter()) {
        assert!(
            lat_lon_in_range(p),
            "限界線の点が妥当な緯度経度域にない: lat={} lon={}",
            lat_deg(p),
            lon_deg(p)
        );
    }
}

/// FAST / 新規（**厳密性の主検証**）: μ'≠0 の合成中心食で、北/南限各点が厳密 2 条件
/// （自己整合ζの錐exact＋相対速度包絡⊥）を満たす。前方射影で各点を基本面へ戻して独立に検証する。
///
/// 殺す変異（厳密化前 RED の主因）:
/// - |L2'| を中心軸 ζ₀ で計算する（点自身の ζ でない）→ 条件1 で距離不一致。
/// - rel に μ' 項を含めず影速度 (x',y') のみに垂直とする geometric 近似 → 条件2 で dot≠0（μ'≠0 ゆえ）。
/// - rel_vx/rel_vy の μ' 項・cos d↔sin d・ξ↔η・ζ の取り違え、tan_f2 の符号反転、l1↔l2 取り違え。
#[test]
fn synthetic_limits_satisfy_exact_cone_and_envelope_conditions() {
    let engine = standard_engine(bundled_time_data());
    let bessel = rigorous_bessel();
    let eclipse = central_eclipse_with_bessel(bessel.clone(), 1.0);

    let path = engine
        .path(&eclipse, PathOptions::default())
        .expect("中心食の path() は成功する");
    let north = path.northern_limit.as_ref().expect("northern_limit=Some");
    let south = path.southern_limit.as_ref().expect("southern_limit=Some");

    // μ'≠0 を独立に確認（このテストの分離力の前提）。
    assert!(
        bessel.mu.derivative().eval(0.0).abs() > 1e-6,
        "μ'≠0 構成（rel 速度に μ' が効く）"
    );

    let u1 = eclipse.global.central_begin.as_ref().unwrap().time_tt;
    let u4 = eclipse.global.central_end.as_ref().unwrap().time_tt;
    let times = lockstep_sample_times(u1, u4, PathOptions::default().sample_interval_seconds);
    // 合成は前方射影が厳密に閉じるので強い許容で締める。dot のスケールは |off|·|rel|~0.01·0.5。
    assert_exact_limit_conditions(north, south, &bessel, &times, 1e-7, 1e-9);
}

/// FAST / 新規: 対応サンプルで **北限緯度 ≥ 中心線緯度 ≥ 南限緯度**（高緯度側=北限・低緯度側=南限）。
/// 近似ゆえ等号許容（小マージン）。3 本が同じサンプル時刻列＝同点数で並ぶことも縛る。
///
/// 殺す変異: 北限と南限を入れ替える・中心線の外に両方とも同じ側へずらす・限界線を中心線のコピーにする
///   （緯度差ゼロ）・サンプル列を北/南で食い違わせる。
#[test]
fn northern_limit_is_north_of_center_is_north_of_southern() {
    let engine = standard_engine(bundled_time_data());
    let eclipse = central_eclipse(geo(0.0, 0.0));

    let path = engine
        .path(&eclipse, PathOptions::default())
        .expect("中心食の path() は成功する");

    let center = path.center_line.as_ref().expect("center_line=Some");
    let north = path.northern_limit.as_ref().expect("northern_limit=Some");
    let south = path.southern_limit.as_ref().expect("southern_limit=Some");

    // 同一サンプル時刻列＝3 本が同点数（取り違え・食い違いを撃破）。
    assert_eq!(
        north.points.len(),
        center.points.len(),
        "北限と中心線は同点数（同サンプル列）"
    );
    assert_eq!(
        south.points.len(),
        center.points.len(),
        "南限と中心線は同点数（同サンプル列）"
    );

    // 各対応サンプルで 北 ≥ 中心 ≥ 南（近似ゆえ微小マージン許容）。
    const EPS: f64 = 1.0e-6;
    for i in 0..center.points.len() {
        let n = lat_deg(&north.points[i]);
        let c = lat_deg(&center.points[i]);
        let s = lat_deg(&south.points[i]);
        assert!(
            n >= c - EPS,
            "サンプル{i}: 北限緯度 {n} ≥ 中心線緯度 {c}（高緯度側=北限）"
        );
        assert!(
            c >= s - EPS,
            "サンプル{i}: 中心線緯度 {c} ≥ 南限緯度 {s}（低緯度側=南限）"
        );
    }
}

/// FAST / 新規: 北限と南限が**分離**している（帯幅が正で過大でない）。代表サンプルで北限・南限の
/// 緯度差が正、かつ各対応点の概算距離が 0 でなく数百 km オーダー（数度未満＝過大でない）に収まる。
/// 脆い絶対値固定は避け、下限>0・上限を緩い帯（< ~5°/~600 km）で縛る。
///
/// 殺す変異: 北限=南限（幅ゼロ）にする・帯幅を桁違いに大きく（地球規模）/小さく（数値誤差）する・
///   ±オフセットの符号を片側に倒して幅を消す。
#[test]
fn limits_are_separated_with_plausible_band_width() {
    let engine = standard_engine(bundled_time_data());
    let eclipse = central_eclipse(geo(0.0, 0.0));

    let path = engine
        .path(&eclipse, PathOptions::default())
        .expect("中心食の path() は成功する");
    let north = path.northern_limit.as_ref().expect("northern_limit=Some");
    let south = path.southern_limit.as_ref().expect("southern_limit=Some");
    let center = path.center_line.as_ref().expect("center_line=Some");

    // 代表サンプル（中央付近）で北限・南限が分離している。
    let mid = center.points.len() / 2;
    let n = &north.points[mid];
    let s = &south.points[mid];

    let dlat = lat_deg(n) - lat_deg(s);
    assert!(
        dlat > 0.0,
        "中央サンプルで北限緯度 - 南限緯度 = {dlat} > 0（帯が分離）"
    );
    assert!(
        dlat < 5.0,
        "帯幅（緯度差 {dlat}°）が過大でない（< 5°＝数百 km オーダー）"
    );

    let band_km = great_circle_km(n, s);
    assert!(
        band_km > 1.0,
        "北限-南限の概算距離 {band_km} km が 0 でない（>1 km）"
    );
    assert!(
        band_km < 600.0,
        "北限-南限の概算距離 {band_km} km が過大でない（< 600 km）"
    );

    // 中心点が北限と南限の緯度で挟まれる（代表サンプル）。
    let c_lat = lat_deg(&center.points[mid]);
    assert!(
        lat_deg(s) <= c_lat && c_lat <= lat_deg(n),
        "中心点緯度 {c_lat} が南限 {} ≤ ・北限 {} ≥ で挟まれる",
        lat_deg(s),
        lat_deg(n)
    );
}

// ============================================================
// FAST: include_limits=false → 限界線 None（center_line は Some のまま）
// ============================================================

/// FAST / 新規: 中心食でも include_limits=false なら northern/southern_limit=None。center_line は
/// include_limits に依らず中心食なら Some（限界線フラグが中心線生成を巻き込まない）。
///
/// 殺す変異: include_limits を無視して常に限界線を作る・include_limits=false で center_line まで
///   None にする・フラグを反転して解釈する。
#[test]
fn include_limits_false_yields_no_limits_but_keeps_center_line() {
    let engine = standard_engine(bundled_time_data());
    let eclipse = central_eclipse(geo(0.0, 0.0));

    let path = engine
        .path(
            &eclipse,
            PathOptions {
                include_limits: false,
                ..PathOptions::default()
            },
        )
        .expect("中心食の path() は成功する");

    assert!(
        path.northern_limit.is_none(),
        "include_limits=false では northern_limit=None"
    );
    assert!(
        path.southern_limit.is_none(),
        "include_limits=false では southern_limit=None"
    );
    assert!(
        path.center_line.is_some(),
        "include_limits=false でも中心食なら center_line=Some"
    );
    // partial_limit / samples は本スライスでも未生成。
    assert!(path.partial_limit.is_none(), "partial_limit None");
    assert!(
        path.samples.is_empty(),
        "samples 空, got {}",
        path.samples.len()
    );
}

// ============================================================
// FAST: 非中心 → 限界線も中心線も None（include_limits=true でも）
// ============================================================

/// FAST / 新規: 非中心（central_begin か central_end が None）では include_limits=true でも
/// northern/southern_limit=None・center_line=None。「両方 Some」のときだけ限界線を作る（&& 条件）。
///
/// 殺す変異: 非中心でも限界線を作る・central_begin/end の片方だけ見て限界線を作る（|| 化）・
///   include_limits=true なら無条件で限界線を出す。
#[test]
fn noncentral_eclipse_has_no_limits() {
    let engine = standard_engine(bundled_time_data());
    let greatest_position = geo(-33.0, 151.0);

    for (with_begin, with_end) in [(false, false), (true, false), (false, true)] {
        let eclipse = noncentral_eclipse(greatest_position, with_begin, with_end);
        let path = engine
            .path(&eclipse, PathOptions::default())
            .expect("非中心の path() も成功する");

        assert!(
            path.northern_limit.is_none(),
            "非中心(begin={with_begin}, end={with_end}) では northern_limit=None"
        );
        assert!(
            path.southern_limit.is_none(),
            "非中心(begin={with_begin}, end={with_end}) では southern_limit=None"
        );
        assert!(
            path.center_line.is_none(),
            "非中心(begin={with_begin}, end={with_end}) では center_line=None"
        );
    }
}

// ============================================================
// SLOW: 実 2024-04-08 皆既を search → path() の限界線が NASA 帯幅域
// ============================================================

/// SLOW / 新規（**厳密化に伴う狭帯＋実日食での 2 条件**）: 実エンジンで 2024-04-08 皆既を
/// search → path()。北/南限界線が Some・各点妥当・中心線が南北の間にあり、(a) 最大食付近の帯幅が
/// NASA 公表 197.5 km の妥当域 [185, 215] km、(b) 最大食付近サンプルで限界点が厳密 2 条件
/// （自己整合ζの錐exact＋相対速度包絡⊥）を満たすことを縛る。de440s 不要（解析暦）。
///
/// 帯幅域 [185, 215] km の根拠: NASA 公表 197.5 km（18:16/18:18 で 197–198 km）に対し ±~7%。
/// 残差源は k 値（IAU mean lunar radius vs NASA 限界用 k=0.2725076 で l2 が ~1–2% スケール）・ΔT・
/// 解析暦差・最大食に最も近いサンプルが ≤30 s ズレること。厳密化前の geometric 近似は影速度のみに
/// 垂直で |L2'| を中心軸 ζ₀ で測るため帯幅がこの狭域から外れる（過去の緩い [100,350] は近似ゆえ）。
///
/// (b) の実日食許容: cone_tol/dot_tol は合成（厳密に閉じる）より緩く取る。前方射影自体は厳密だが、
/// l2/tan_f2/d/μ の実暦評価と中心線位置律速で微小残差が乗るため。NASA 緯度経度の直接一致は
/// 中心線位置精度律速ゆえ縛らない（帯幅と 2 条件で締める）。
///
/// 殺す変異: 限界線を捏造/空にする・南北を取り違える・帯幅を桁違いにする・中心線が帯の外に出る・
///   実日食で限界線を生成しない・rel に μ' を含めない近似のまま（条件2 が実日食でも崩れる）。
#[test]
fn real_2024_eclipse_limits_match_nasa_band_width() {
    let engine = standard_engine(bundled_time_data());
    let range = umbra_core::TimeRange {
        start: utc(2024, 4, 8, 0, 0, 0.0),
        end: utc(2024, 4, 9, 0, 0, 0.0),
    };
    let eclipses = engine
        .search(range)
        .expect("2024-04-08 範囲の search は成功する");
    let eclipse = eclipses
        .iter()
        .find(|e| matches!(e.kind, SolarEclipseKind::Total))
        .expect("2024-04-08 皆既が見つかる");

    let path = engine
        .path(eclipse, PathOptions::default())
        .expect("実皆既の path() は成功する");

    let center = path
        .center_line
        .as_ref()
        .expect("皆既なので center_line=Some");
    let north = path
        .northern_limit
        .as_ref()
        .expect("皆既なので northern_limit=Some");
    let south = path
        .southern_limit
        .as_ref()
        .expect("皆既なので southern_limit=Some");

    assert!(
        north.points.len() >= 2,
        "実北限界線は ≥2 点, got {}",
        north.points.len()
    );
    assert_eq!(
        north.points.len(),
        center.points.len(),
        "北限と中心線は同サンプル列"
    );
    assert_eq!(
        south.points.len(),
        center.points.len(),
        "南限と中心線は同サンプル列"
    );

    for p in north.points.iter().chain(south.points.iter()) {
        assert!(
            lat_lon_in_range(p),
            "実限界線の点が妥当な緯度経度域にない: lat={} lon={}",
            lat_deg(p),
            lon_deg(p)
        );
    }

    // 各対応サンプルで 北 ≥ 中心 ≥ 南（近似ゆえ微小マージン許容）。
    const EPS: f64 = 1.0e-3;
    for i in 0..center.points.len() {
        let n = lat_deg(&north.points[i]);
        let c = lat_deg(&center.points[i]);
        let s = lat_deg(&south.points[i]);
        assert!(n >= c - EPS, "サンプル{i}: 北限 {n} ≥ 中心 {c}");
        assert!(c >= s - EPS, "サンプル{i}: 中心 {c} ≥ 南限 {s}");
    }

    // 最大食点に最も近いサンプル付近で帯幅 ~197 km の妥当域（geometric 近似ゆえ広め 100–350 km）。
    let g_lat = lat_deg(&path.greatest_point);
    let g_lon = lon_deg(&path.greatest_point);
    let mid = (0..center.points.len())
        .min_by(|&a, &b| {
            let da = {
                let dlat = lat_deg(&center.points[a]) - g_lat;
                let dlon = lon_deg(&center.points[a]) - g_lon;
                dlat * dlat + dlon * dlon
            };
            let db = {
                let dlat = lat_deg(&center.points[b]) - g_lat;
                let dlon = lon_deg(&center.points[b]) - g_lon;
                dlat * dlat + dlon * dlon
            };
            da.partial_cmp(&db).expect("有限距離")
        })
        .expect("中心線は非空");

    let band_km = great_circle_km(&north.points[mid], &south.points[mid]);
    assert!(
        (185.0..=215.0).contains(&band_km),
        "最大食付近の帯幅 {band_km} km が NASA 197.5 km の妥当域 [185, 215] に入る（厳密錐接線解）"
    );

    // (b) 実日食でも厳密 2 条件を縛る（最大食付近の数サンプル）。実暦・中心線位置律速ゆえ
    // 合成より緩い許容（cone_tol 5e-4 Re ≈ 3 km, dot_tol は |off|·|rel| スケールに対し緩め）。
    // 全サンプルで回すと SLOW がさらに重くなるため、最大食付近の窓に限定する。
    //
    // 実 path() は grazing/RootNotBracketed のサンプルをスキップするため index↔時刻がずれる。
    // よって lockstep の times[i] を信用せず、各 kept index について中心線点から実サンプル時刻を
    // 二分法で復元する `assert_exact_limit_conditions_real` を使う（同サンプル列の center を渡す）。
    let u1 = eclipse.global.central_begin.as_ref().unwrap().time_tt;
    let u4 = eclipse.global.central_end.as_ref().unwrap().time_tt;
    let lo = mid.saturating_sub(2);
    let hi = (mid + 3).min(north.points.len());
    let win_center = umbra_geo::GeoLine::new(center.points[lo..hi].to_vec());
    let win_north = umbra_geo::GeoLine::new(north.points[lo..hi].to_vec());
    let win_south = umbra_geo::GeoLine::new(south.points[lo..hi].to_vec());
    assert_exact_limit_conditions_real(
        &win_center,
        &win_north,
        &win_south,
        &eclipse.bessel,
        u1,
        u4,
        5.0e-4,
        5.0e-6,
    );
}

// ============================================================
// M9.7: 経路サンプル列 samples（中心食で充足・center_line/南北限界線と lockstep）
//
// 確定仕様（観測可能な契約）:
//  1. 中心食（central_begin/central_end 両方 Some）かつ include_limits=true で
//     samples.len() == center_line.len() == northern_limit.len() == southern_limit.len()、
//     samples[i].center == center_line.points[i]（完全 lockstep）。
//  2. 中心食でも include_limits=false なら samples は空（限界線 None と整合）。
//  3. 非中心/部分食では samples 空（center_line=None と整合）。
//  各 PathSample フィールド:
//   - time_utc = tt_to_utc(サンプル時刻 TT)。U1〜U4 内で単調増加。
//   - center = center_line.points[i]（影軸地表点）。
//   - duration_seconds = 2|L2'|/|rel|×3600（中心軸 ζ で評価。M9.6 と同定義）。
//   - sun_altitude = その時刻・中心点の幾何高度（RefractionModel::None）。
//   - path_width = 南北本影縁点間の大圏距離（M9.6 と同定義 = great_circle(north[i],south[i])）。
//   - kind = L2'=l2−ζ·tan f2 の符号（<0=Total / それ以外=Annular）。
// ============================================================

/// サンプル時刻 TT を、その中心点を前方射影した ζ を使って M9.6 と同方式で
/// `duration_seconds` の期待値を独立に組む（被テスト関数の戻りは流用しない）。
/// rel は中心軸 (ξ=x, η=y, ζ) の地表相対速度（μ' 項込み）。
fn expected_duration_seconds(
    center_point: &umbra_geo::GeoPoint,
    bessel: &BesselianPolynomial,
    t: TtInstant,
) -> f64 {
    let e = bessel.at(t).expect("区間内サンプルは評価成功");
    let zeta = forward_project(center_point, &e).zeta;
    let epoch = bessel.epoch_tt;
    let t_hours = t.jd2().days_since(epoch.jd2()) * 24.0;
    let vx = bessel.x.derivative().eval(t_hours);
    let vy = bessel.y.derivative().eval(t_hours);
    let mu_rate = bessel.mu.derivative().eval(t_hours);
    let (sin_d, cos_d) = e.declination.0.sin_cos();
    let rel_x = vx - mu_rate * (zeta * cos_d - e.y * sin_d);
    let rel_y = vy - mu_rate * e.x * sin_d;
    let rel_speed = rel_x.hypot(rel_y);
    let l2p_abs = (e.l2 - zeta * e.tan_f2).abs();
    2.0 * l2p_abs / rel_speed * 3600.0
}

/// FAST / 新規（**lockstep の主検証**）: 中心食＋include_limits=true で samples が
/// center_line・北限・南限と完全に同一サンプル列（同点数）になり、samples[i].center が
/// center_line.points[i] と一致する。samples は非空（≥2）。
///
/// 殺す変異: samples を常に空にする（M9.6 以前の挙動）・samples 長さを center/限界線とズラす
///   （off-by-one・別ループ上限）・samples[i].center に限界線点や別 index の点を入れる
///   （center↔north/south 取り違え、index ズレ）。
#[test]
fn samples_are_lockstep_with_center_and_limit_lines() {
    let engine = standard_engine(bundled_time_data());
    let eclipse = central_eclipse(geo(0.0, 0.0));

    let path = engine
        .path(&eclipse, PathOptions::default())
        .expect("中心食の path() は成功する");

    let center = path.center_line.as_ref().expect("center_line=Some");
    let north = path.northern_limit.as_ref().expect("northern_limit=Some");
    let south = path.southern_limit.as_ref().expect("southern_limit=Some");

    assert!(
        path.samples.len() >= 2,
        "中心食では samples は非空（≥2）, got {}",
        path.samples.len()
    );
    // 4 本が完全に同点数（lockstep）。
    assert_eq!(
        path.samples.len(),
        center.points.len(),
        "samples.len() == center_line.len()"
    );
    assert_eq!(
        path.samples.len(),
        north.points.len(),
        "samples.len() == northern_limit.len()"
    );
    assert_eq!(
        path.samples.len(),
        south.points.len(),
        "samples.len() == southern_limit.len()"
    );

    // samples[i].center == center_line.points[i]（影軸地表点・index 一致）。
    for (i, s) in path.samples.iter().enumerate() {
        assert_eq!(
            s.center, center.points[i],
            "samples[{i}].center が center_line.points[{i}] と一致しない（center 取り違え/index ズレ）"
        );
    }
}

/// FAST / 新規（**フィールド・オラクルの主検証**）: μ'≠0 の合成中心食で、各サンプルの
/// duration_seconds・path_width・kind を独立オラクルで縛る。
///   - duration_seconds = 2|L2'|/|rel|×3600（中心点 ζ・μ' 項込み・M9.6 同方式）。
///   - path_width = great_circle(north[i], south[i])（M9.6 同方式）。
///   - kind = L2'<0 → Total（合成は l2=−0.009<0 ゆえ Total）。
///
/// duration（秒・~200 s）と path_width（km・<1000）は桁が異なる非対称値なので、両者を
/// 取り違える変異は両方のオラクルを同時に外す。
///
/// 殺す変異:
/// - duration↔path_width フィールド取り違え（秒 vs km で両域同時に外れる）。
/// - duration の 2× / ½ / |rel| の逆数誤り・×3600 脱落（秒域外）。
/// - |L2'| を中心軸 ζ₀=0 で測る（点自身の ζ でない）→ duration ズレ（μ' で ζ≠0）。
/// - rel に μ' 項を含めない近似 → duration ズレ（μ'≠0）。
/// - path_width を北限・南限以外（中心線等）から測る・南北片側だけ → 距離ズレ。
/// - kind の符号反転（L2'<0 を Annular にする）。
#[test]
fn samples_field_values_match_independent_oracles() {
    let engine = standard_engine(bundled_time_data());
    let bessel = rigorous_bessel();
    let eclipse = central_eclipse_with_bessel(bessel.clone(), 1.0);

    let path = engine
        .path(&eclipse, PathOptions::default())
        .expect("中心食の path() は成功する");
    let north = path.northern_limit.as_ref().expect("northern_limit=Some");
    let south = path.southern_limit.as_ref().expect("southern_limit=Some");

    // μ'≠0 を独立確認（ζ・rel オラクルの分離力の前提）。
    assert!(
        bessel.mu.derivative().eval(0.0).abs() > 1e-6,
        "μ'≠0 構成（rel 速度に μ' が効く）"
    );

    let u1 = eclipse.global.central_begin.as_ref().unwrap().time_tt;
    let u4 = eclipse.global.central_end.as_ref().unwrap().time_tt;
    let times = lockstep_sample_times(u1, u4, PathOptions::default().sample_interval_seconds);

    assert_eq!(
        path.samples.len(),
        times.len(),
        "samples 列が lockstep 時刻列と同点数（時刻復元の前提）"
    );

    for (i, s) in path.samples.iter().enumerate() {
        let t = times[i];

        // duration_seconds: 中心点 ζ・μ' 項込みの M9.6 式で独立に組む。
        let want_dur = expected_duration_seconds(&s.center, &bessel, t);
        assert!(
            (s.duration_seconds - want_dur).abs() < 1e-6,
            "samples[{i}].duration_seconds {} != 2|L2'|/|rel|×3600 期待 {want_dur}",
            s.duration_seconds
        );

        // path_width: 南北本影縁点間の大圏距離（M9.6 同方式）。
        let want_width = great_circle_km(&north.points[i], &south.points[i]);
        // great_circle_km は haversine（実装側と同一近似でなくてもよい）ゆえ相対 1% 許容。
        assert!(
            (s.path_width.0 - want_width).abs() <= 1.0e-2 * want_width.max(1.0),
            "samples[{i}].path_width {} km != great_circle(north,south) 期待 {want_width} km",
            s.path_width.0
        );
        // 桁の独立性: 帯幅(km) と継続(秒) は別物（取り違え検出の補強）。
        assert!(
            (s.path_width.0 - s.duration_seconds).abs() > 1.0,
            "samples[{i}]: path_width と duration_seconds が同値（取り違えの疑い）"
        );

        // kind: 合成は l2=−0.009<0 ⇒ L2'<0 ⇒ Total。
        assert_eq!(
            s.kind,
            SolarEclipseKind::Total,
            "samples[{i}].kind は L2'<0 ゆえ Total（符号規約 l2<0=皆既）"
        );

        // 各フィールドが有限・妥当域。
        assert!(
            s.duration_seconds.is_finite() && s.duration_seconds > 0.0,
            "samples[{i}].duration_seconds {} は正・有限",
            s.duration_seconds
        );
        assert!(
            s.path_width.0.is_finite() && s.path_width.0 > 0.0,
            "samples[{i}].path_width {} は正・有限",
            s.path_width.0
        );
        assert!(
            (-90.0..=90.0).contains(&s.sun_altitude.0) && s.sun_altitude.0.is_finite(),
            "samples[{i}].sun_altitude {}° は [-90,90] で有限",
            s.sun_altitude.0
        );
        assert!(
            lat_lon_in_range(&s.center),
            "samples[{i}].center が妥当な緯度経度域にない"
        );
    }
}

/// FAST / 新規: 各サンプルの time_utc が tt_to_utc(サンプル時刻 TT) と一致し、列全体で
/// 単調増加する。サンプル時刻 TT は lockstep 時刻列から独立再構成する。
///
/// 殺す変異: time_utc に TT をそのまま入れる（UTC 変換脱落・ΔT 分ズレ）・別 index の時刻を
///   入れる（時刻↔index 取り違え）・時刻列を逆順/定数にする（単調増加が崩れる）。
#[test]
fn samples_time_utc_equals_tt_to_utc_and_is_monotonic() {
    let engine = standard_engine(bundled_time_data());
    let eclipse = central_eclipse(geo(0.0, 0.0));

    let path = engine
        .path(&eclipse, PathOptions::default())
        .expect("中心食の path() は成功する");

    let u1 = eclipse.global.central_begin.as_ref().unwrap().time_tt;
    let u4 = eclipse.global.central_end.as_ref().unwrap().time_tt;
    let times = lockstep_sample_times(u1, u4, PathOptions::default().sample_interval_seconds);
    assert_eq!(
        path.samples.len(),
        times.len(),
        "samples 列が lockstep 時刻列と同点数（時刻復元の前提）"
    );

    let mut prev_jd = f64::NEG_INFINITY;
    for (i, s) in path.samples.iter().enumerate() {
        // time_utc == tt_to_utc(その TT)。同一瞬時（< 1ms 相当）。
        let want_utc = umbra_core::time::tt_to_utc(times[i])
            .expect("サンプル TT は post-1972 で UTC 変換可能");
        let got_jd = s.time_utc.jd2().jd();
        let want_jd = want_utc.jd2().jd();
        assert!(
            (got_jd - want_jd).abs() < 1.0 / 86_400.0,
            "samples[{i}].time_utc == tt_to_utc(time_tt): got_jd={got_jd} want_jd={want_jd}"
        );
        // 単調増加（厳密に増加・等間隔サンプル）。
        assert!(
            got_jd > prev_jd,
            "samples[{i}].time_utc が単調増加でない: got_jd={got_jd} prev_jd={prev_jd}"
        );
        prev_jd = got_jd;
    }
    // U1〜U4 の範囲内（始点 ≥ U1 相当・終点 ≤ U4 相当を UTC で確認）。
    let first = path
        .samples
        .first()
        .expect("samples 非空")
        .time_utc
        .jd2()
        .jd();
    let last = path
        .samples
        .last()
        .expect("samples 非空")
        .time_utc
        .jd2()
        .jd();
    let u1_utc = umbra_core::time::tt_to_utc(u1).unwrap().jd2().jd();
    let u4_utc = umbra_core::time::tt_to_utc(u4).unwrap().jd2().jd();
    assert!(
        first >= u1_utc - 1.0 / 86_400.0 && last <= u4_utc + 1.0 / 86_400.0,
        "samples の時刻が [U1,U4] 内: first={first} last={last} U1={u1_utc} U4={u4_utc}"
    );
}

/// FAST / 新規: 中心食でも include_limits=false なら samples は空（限界線が None になるのと整合）。
/// center_line は include_limits に依らず Some のまま（samples 空が center_line 生成を巻き込まない）。
///
/// 殺す変異: include_limits を無視して常に samples を作る・include_limits=false で
///   center_line まで None にする・samples 充足を限界線フラグから切り離して常時充足する。
#[test]
fn include_limits_false_yields_empty_samples_but_keeps_center_line() {
    let engine = standard_engine(bundled_time_data());
    let eclipse = central_eclipse(geo(0.0, 0.0));

    let path = engine
        .path(
            &eclipse,
            PathOptions {
                include_limits: false,
                ..PathOptions::default()
            },
        )
        .expect("中心食の path() は成功する");

    assert!(
        path.samples.is_empty(),
        "include_limits=false では samples 空, got {}",
        path.samples.len()
    );
    // 限界線も None（整合）。center_line は Some のまま。
    assert!(
        path.northern_limit.is_none(),
        "include_limits=false で northern_limit=None"
    );
    assert!(
        path.southern_limit.is_none(),
        "include_limits=false で southern_limit=None"
    );
    assert!(
        path.center_line.is_some(),
        "include_limits=false でも中心食なら center_line=Some"
    );
}

/// FAST / 新規: 非中心（central_begin か central_end が None）では include_limits=true でも
/// samples 空（center_line=None と整合）。「両方 Some」のときだけ samples を充足する（&& 条件）。
///
/// 殺す変異: 非中心でも samples を作る・central_begin/end の片方だけ見て充足（|| 化）・
///   include_limits=true なら無条件で samples を出す。
#[test]
fn noncentral_eclipse_has_empty_samples() {
    let engine = standard_engine(bundled_time_data());
    let greatest_position = geo(-33.0, 151.0);

    for (with_begin, with_end) in [(false, false), (true, false), (false, true)] {
        let eclipse = noncentral_eclipse(greatest_position, with_begin, with_end);
        let path = engine
            .path(&eclipse, PathOptions::default())
            .expect("非中心の path() も成功する");

        assert!(
            path.samples.is_empty(),
            "非中心(begin={with_begin}, end={with_end}) では samples 空, got {}",
            path.samples.len()
        );
        assert!(
            path.center_line.is_none(),
            "非中心(begin={with_begin}, end={with_end}) では center_line=None（samples と整合）"
        );
    }
}

// ============================================================
// SLOW: 実 2024-04-08 皆既の samples（lockstep・各フィールド NASA 域・全 Total・単調 UTC）
// ============================================================

/// SLOW / 新規: 実エンジンで 2024-04-08 皆既を search → path()。samples が center_line/北限/南限と
/// 同点数（lockstep）かつ samples[i].center == center_line.points[i]。最大食付近のサンプルで
/// path_width ∈ [185,215] km・duration_seconds ∈ [250,286] s（NASA ~197.5 km / ~268 s の妥当域、
/// 既存 SLOW テストと同域）。全サンプルの kind=Total（2024 は皆既）。time_utc は U1〜U4 内で単調増加。
/// de440s 不要（解析暦）。
///
/// 殺す変異: 実日食で samples を空/捏造にする・lockstep 長さズレ・center に限界線点を入れる・
///   width↔duration 取り違え（km vs 秒で両域同時外し）・duration の 2×/½・kind を Annular にする
///   （皆既で金環）・time_utc に TT を入れる（ΔT 分ズレ）・時刻を非単調にする。
#[test]
fn real_2024_eclipse_samples_lockstep_and_field_domains() {
    let engine = standard_engine(bundled_time_data());
    let range = umbra_core::TimeRange {
        start: utc(2024, 4, 8, 0, 0, 0.0),
        end: utc(2024, 4, 9, 0, 0, 0.0),
    };
    let eclipses = engine
        .search(range)
        .expect("2024-04-08 範囲の search は成功する");
    let eclipse = eclipses
        .iter()
        .find(|e| matches!(e.kind, SolarEclipseKind::Total))
        .expect("2024-04-08 皆既が見つかる");

    let path = engine
        .path(eclipse, PathOptions::default())
        .expect("実皆既の path() は成功する");

    let center = path
        .center_line
        .as_ref()
        .expect("皆既なので center_line=Some");
    let north = path
        .northern_limit
        .as_ref()
        .expect("皆既なので northern_limit=Some");
    let south = path
        .southern_limit
        .as_ref()
        .expect("皆既なので southern_limit=Some");

    assert!(
        path.samples.len() >= 2,
        "実皆既の samples は非空（≥2）, got {}",
        path.samples.len()
    );
    // lockstep: samples == center == north == south（同点数）。
    assert_eq!(
        path.samples.len(),
        center.points.len(),
        "samples==center_line 同点数"
    );
    assert_eq!(
        path.samples.len(),
        north.points.len(),
        "samples==northern_limit 同点数"
    );
    assert_eq!(
        path.samples.len(),
        south.points.len(),
        "samples==southern_limit 同点数"
    );

    // samples[i].center == center_line.points[i]、全 kind=Total、time_utc 単調増加。
    let u1 = eclipse.global.central_begin.as_ref().unwrap().time_tt;
    let u4 = eclipse.global.central_end.as_ref().unwrap().time_tt;
    let u1_utc = umbra_core::time::tt_to_utc(u1).unwrap().jd2().jd();
    let u4_utc = umbra_core::time::tt_to_utc(u4).unwrap().jd2().jd();

    let mut prev_jd = f64::NEG_INFINITY;
    for (i, s) in path.samples.iter().enumerate() {
        assert_eq!(
            s.center, center.points[i],
            "samples[{i}].center が center_line.points[{i}] と一致しない"
        );
        assert_eq!(
            s.kind,
            SolarEclipseKind::Total,
            "samples[{i}].kind は皆既なので Total"
        );
        let jd = s.time_utc.jd2().jd();
        assert!(jd > prev_jd, "samples[{i}].time_utc が単調増加でない");
        prev_jd = jd;
        assert!(
            jd >= u1_utc - 1.0 / 86_400.0 && jd <= u4_utc + 1.0 / 86_400.0,
            "samples[{i}].time_utc が [U1,U4] 内でない: jd={jd}"
        );
        assert!(
            s.duration_seconds.is_finite() && s.duration_seconds > 0.0,
            "samples[{i}].duration_seconds は正・有限"
        );
        assert!(
            s.path_width.0.is_finite() && s.path_width.0 > 0.0,
            "samples[{i}].path_width は正・有限"
        );
    }

    // 最大食点に最も近いサンプルで path_width / duration_seconds が NASA 妥当域。
    let g_lat = lat_deg(&path.greatest_point);
    let g_lon = lon_deg(&path.greatest_point);
    let mid = (0..center.points.len())
        .min_by(|&a, &b| {
            let da = {
                let dlat = lat_deg(&center.points[a]) - g_lat;
                let dlon = lon_deg(&center.points[a]) - g_lon;
                dlat * dlat + dlon * dlon
            };
            let db = {
                let dlat = lat_deg(&center.points[b]) - g_lat;
                let dlon = lon_deg(&center.points[b]) - g_lon;
                dlat * dlat + dlon * dlon
            };
            da.partial_cmp(&db).expect("有限距離")
        })
        .expect("中心線は非空");

    let width = path.samples[mid].path_width.0;
    let duration = path.samples[mid].duration_seconds;
    assert!(
        (185.0..=215.0).contains(&width),
        "最大食付近 samples[{mid}].path_width {width} km が NASA 妥当域 [185,215] に入る（NASA≈197.5 km）"
    );
    assert!(
        (250.0..=286.0).contains(&duration),
        "最大食付近 samples[{mid}].duration_seconds {duration} s が NASA 妥当域 [250,286] に入る（NASA≈268.1 s）"
    );
}

// ============================================================
// SLOW: 実 2024-04-08 皆既の GreatestEclipse.path_width / central_duration が NASA 値
// （M9.6 — 中心食の帯幅 path_width と中心食継続 central_duration を Some・NASA ballpark に縛る）
// ============================================================

/// SLOW / 新規（**M9.6 強オラクル**）: 実エンジンで 2024-04-08 皆既を search し、
/// `eclipse.global.greatest.path_width` ≈ 197.5 km・`central_duration` ≈ 268.1 s（4m28.1s）の
/// NASA 公表値の妥当域に入ることを縛る（中心食ゆえともに Some）。de440s 不要（解析暦）。
///
/// 量の定義（オラクル根拠・実装式は写経しない）:
/// - path_width [km] = 最大食時刻の本影帯の北縁-南縁の地表点間 大圏距離（M9.4 限界線・相対速度包絡⊥）。
/// - central_duration [s] = 2·|L2'|/|rel|（umbra 直径 ÷ 影の地表相対速度）。
///
/// 帯幅域 [185, 215] km の根拠: NASA 公表 197.5 km（18:16/18:18 で 197–198 km）に ±~7%
/// （`real_2024_eclipse_limits_match_nasa_band_width` と同じ・k 値/ΔT/解析暦差/最大食サンプルズレ）。
/// 継続域 [250, 286] s の根拠: NASA 公表 268.1 s（4m28.1s）に ±~7%。NASA 秒/km の等値ハードコードは
/// 禁止（conventions §11）＝範囲 check に限定。`umbra_core::Kilometers` から `.0` で値を取り出す。
///
/// 殺す変異: 中心食で path_width/central_duration を None にする（None↔Some 分岐）・width↔duration の
///   取り違え（km vs 秒で桁が違い両域を同時に外す）・2 倍/半分（範囲外）・|rel| の逆数誤り（duration 域外）。
#[test]
fn real_2024_greatest_path_width_and_central_duration_match_nasa() {
    let engine = standard_engine(bundled_time_data());
    let range = umbra_core::TimeRange {
        start: utc(2024, 4, 8, 0, 0, 0.0),
        end: utc(2024, 4, 9, 0, 0, 0.0),
    };
    let eclipses = engine
        .search(range)
        .expect("2024-04-08 範囲の search は成功する");
    let eclipse = eclipses
        .iter()
        .find(|e| matches!(e.kind, SolarEclipseKind::Total))
        .expect("2024-04-08 皆既が見つかる");

    let greatest = &eclipse.global.greatest;

    // 中心食（皆既）ゆえ path_width・central_duration はともに Some。
    let width = greatest
        .path_width
        .expect("2024 皆既は中心食ゆえ path_width=Some（M9.6）");
    let duration = greatest
        .central_duration
        .expect("2024 皆既は中心食ゆえ central_duration=Some（M9.6）");

    // 帯幅 [km]: NASA 公表 197.5 km の妥当域 ±~7%。Kilometers から .0 で取り出す。
    assert!(
        (185.0..=215.0).contains(&width.0),
        "2024 greatest path_width {} km not in NASA ballpark [185,215] (NASA≈197.5 km)",
        width.0
    );
    // 継続 [s]: NASA 公表 268.1 s（4m28.1s）の妥当域 ±~7%。
    assert!(
        (250.0..=286.0).contains(&duration),
        "2024 greatest central_duration {duration} s not in NASA ballpark [250,286] (NASA≈268.1 s)"
    );
}

// ============================================================
// M9 残(3) (3f): 部分食域 partial_limit ＝ **3 領域の多角形ユニオン**
//
// 確定仕様（docs/algorithms/11-path-partial-domain.md §11.6 (3f) 確定仕様・受け入れテスト戦略）:
//   部分食 phase（global.partial_begin / partial_end が両方 Some）かつ include_limits=true の日食で
//   path().partial_limit = Some(GeoPolygon)。部分食域は
//     D = 帯(band) ∪ morning lune ∪ evening lune
//   の **3 領域の和**であり、`partial_limit` はそのユニオン結果の**最大面積成分ひとつ**
//   （`rings[0]`=外環・`rings[1..]`=穴。複数成分は落とす＝単一 GeoPolygon）。
//   include_limits=false／部分 phase 無し（P1 or P4 None）では None。center_line（中心食のみ）と独立。
//   `sample_interval_seconds=0.0`（1 サンプル）は退化して None。
//
// **撤廃された契約（(3f) §11.6(e)）**: 「外環＝北 n ++ 南 n 逆順」のリボン位相（偶数頂点・前半/後半の
//   同時刻北南対・経度単調）と、「**全**頂点が半影縁条件を**厳密に**満たす」はもはや成り立たない
//   （ユニオンは入力リングの断片を繋ぎ直し、交点頂点を導入する）。
//
// オラクル戦略（strict・観測契約 §11.6(d) / 受け入れテスト戦略）:
//   - 存在: 上記の Some/None 分岐を縛る（維持）。
//   - **頂点の正当性（緩和後）**: 各外環頂点 P を **検証済み前方射影**（forward_project）で基本面へ戻し、
//       ある [P1,P4] サンプル時刻 t で **閉半影内**（面内距離 ≤ |l1 − ζ·tan f1| + tol）**かつ昼面側**
//       （ζ ≥ −tol）であること。半影縁点は等号・交点頂点は折れ線離散化の範囲内で満たす。
//       期待値は bessel 多項式（path とは独立な入力）から組む＝被テスト関数の戻りを流用しない。
//   - **単純性**: 外環は自己交差しない（隣接しない辺同士が交差しない）。
//   - **穴の扱い**: `rings[0]`=外環・`rings[1..]`=穴。穴は外環の内側・外環より小面積・互いに素。
//   - 非退化: 外環 ≥3 頂点。
//   - **包含（partial ⊃ umbral path）**: 中心線各点が **平面 (lon,lat) ray-casting point-in-polygon** で
//       外環の内側 **かつ どの穴の内側でもない**。
//   - SLOW headline（(3f) 昇格）: 実 2024-04-08 で **中心線の全点が例外なく内包**される
//       （従来の「最大食まわり ±10 サンプル窓」から昇格）。加えて南北端が中心線より外・terminator 頂点 ≥1。
// ============================================================

/// 部分食 phase（partial_begin/partial_end=Some・P1/P4）を持つ合成中心食を `rigorous_bessel` で組む。
/// 部分食区間 [P1,P4] は中心食区間 [U1,U4]（±span_hours）より広く（±partial_span_hours）、
/// `rigorous_bessel` の fit_interval（epoch±2h）に収める。これにより build_partial_limit が
/// [P1,P4] で半影限界・limb 点を収集でき partial_limit=Some になる。
fn partial_eclipse_with_bessel(
    bessel: BesselianPolynomial,
    span_hours: f64,
    partial_span_hours: f64,
) -> SolarEclipse {
    let epoch = bessel.epoch_tt;
    let p1 = contact(tt_at_hours(epoch, -partial_span_hours));
    let p4 = contact(tt_at_hours(epoch, partial_span_hours));
    let u1 = contact(tt_at_hours(epoch, -span_hours));
    let u4 = contact(tt_at_hours(epoch, span_hours));
    let global = GlobalCircumstances {
        kind: SolarEclipseKind::Total,
        partial_begin: Some(p1),
        central_begin: Some(u1),
        greatest: greatest_at(geo(0.0, 0.0)),
        central_end: Some(u4),
        partial_end: Some(p4),
        gamma: 0.05,
    };
    SolarEclipse {
        event_key: "synthetic-partial#0".to_string(),
        kind: SolarEclipseKind::Total,
        global,
        bessel,
        metadata: metadata(),
    }
}

/// 部分食 phase を持つ「部分食 only」（central_begin/end=None）の合成日食。
/// 中心線は出ない（center_line=None）が、部分食域 partial_limit は Some になりうる（中心食と独立）。
fn partial_only_eclipse(bessel: BesselianPolynomial, partial_span_hours: f64) -> SolarEclipse {
    let epoch = bessel.epoch_tt;
    let p1 = contact(tt_at_hours(epoch, -partial_span_hours));
    let p4 = contact(tt_at_hours(epoch, partial_span_hours));
    let global = GlobalCircumstances {
        kind: SolarEclipseKind::Partial,
        partial_begin: Some(p1),
        central_begin: None,
        greatest: greatest_at(geo(0.0, 0.0)),
        central_end: None,
        partial_end: Some(p4),
        gamma: 0.05,
    };
    SolarEclipse {
        event_key: "synthetic-partial-only#0".to_string(),
        kind: SolarEclipseKind::Partial,
        global,
        bessel,
        metadata: metadata(),
    }
}

/// (3f) §11.6(d): 外環頂点 P が **閉半影内かつ昼面側**か（緩和後の「頂点の正当性」）。
/// ある [P1,P4] のサンプル時刻 t で
///   (i) **昼面側**: `ζ ≥ −zeta_tol`（前方射影した P 自身の ζ）
///   (ii) **閉半影内**: 影軸からの基本面内距離 ≤ `|l1 − ζ·tan f1| + dist_tol`（自己整合ζ）
/// を同時に満たせば妥当とする。半影縁点（昼面包絡・terminator 交点）は (ii) を等号で、ユニオンが導入する
/// **交点頂点**は折れ線離散化の範囲内で満たす。どのサンプル時刻由来か不定なので全時刻で試す。
/// 期待値 l1 / tan f1 は bessel 多項式（path とは独立な入力）から組む＝被テスト関数の戻りを流用しない。
///
/// 非対称性: (ii) の半径を |l2|（本影）で測る・l1↔l2 / tan_f1↔tan_f2 取り違えは半影が本影の ~60 倍ゆえ
/// 大多数の頂点が閉半影外に落ちて偽になる。(i) を外すと夜面（ζ<0）の捏造点を見逃す。
fn vertex_is_in_closed_penumbra(
    p: &umbra_geo::GeoPoint,
    bessel: &BesselianPolynomial,
    sample_times: &[TtInstant],
    dist_tol: f64,
    zeta_tol: f64,
) -> bool {
    for t in sample_times {
        let e = match bessel.at(*t) {
            Ok(e) => e,
            Err(_) => continue,
        };
        let of = forward_project(p, &e);
        if of.zeta < -zeta_tol {
            continue; // 夜面側は不可（昼面側の契約）
        }
        let in_plane = (of.xi - e.x).hypot(of.eta - e.y);
        let penumbral = (e.l1 - of.zeta * e.tan_f1).abs();
        if in_plane <= penumbral + dist_tol {
            return true;
        }
    }
    false
}

/// 線分 (a,b) と (c,d) が **真に交差**するか（端点で触れるだけ・共線は false）。
/// 平面 (lon,lat) 度。外環の単純性（自己交差なし）判定に使う。
fn segments_properly_intersect(
    a: &umbra_geo::GeoPoint,
    b: &umbra_geo::GeoPoint,
    c: &umbra_geo::GeoPoint,
    d: &umbra_geo::GeoPoint,
) -> bool {
    let cross = |ox: f64, oy: f64, px: f64, py: f64, qx: f64, qy: f64| {
        (px - ox) * (qy - oy) - (py - oy) * (qx - ox)
    };
    let (ax, ay) = (lon_deg(a), lat_deg(a));
    let (bx, by) = (lon_deg(b), lat_deg(b));
    let (cx, cy) = (lon_deg(c), lat_deg(c));
    let (dx, dy) = (lon_deg(d), lat_deg(d));
    let d1 = cross(ax, ay, bx, by, cx, cy);
    let d2 = cross(ax, ay, bx, by, dx, dy);
    let d3 = cross(cx, cy, dx, dy, ax, ay);
    let d4 = cross(cx, cy, dx, dy, bx, by);
    // 真の交差（両線分が互いを厳密に跨ぐ）のみ。端点接触・共線重なりは false（退化は別契約）。
    ((d1 > 0.0 && d2 < 0.0) || (d1 < 0.0 && d2 > 0.0))
        && ((d3 > 0.0 && d4 < 0.0) || (d3 < 0.0 && d4 > 0.0))
}

/// 閉リングが**単純**（隣接しない辺同士が交差しない）か。交差が見つかれば辺 index 対を返す。
/// 辺 i = (ring[i], ring[(i+1)%n])。隣接辺（index 差 1 / 巻き戻りの先頭-末尾）は端点共有ゆえ除外。
fn ring_self_intersection(ring: &[umbra_geo::GeoPoint]) -> Option<(usize, usize)> {
    let n = ring.len();
    if n < 4 {
        return None;
    }
    for i in 0..n {
        for j in (i + 1)..n {
            // 隣接（端点共有）はスキップ。
            if j == i + 1 || (i == 0 && j == n - 1) {
                continue;
            }
            if segments_properly_intersect(
                &ring[i],
                &ring[(i + 1) % n],
                &ring[j],
                &ring[(j + 1) % n],
            ) {
                return Some((i, j));
            }
        }
    }
    None
}

/// 平面 (lon,lat) の shoelace 符号付き面積 [deg²]（向きに依存・絶対値が面積）。
fn signed_area_deg2(ring: &[umbra_geo::GeoPoint]) -> f64 {
    let n = ring.len();
    let mut s = 0.0;
    let mut j = n - 1;
    for i in 0..n {
        s += (lon_deg(&ring[j]) + lon_deg(&ring[i])) * (lat_deg(&ring[j]) - lat_deg(&ring[i]));
        j = i;
    }
    s / 2.0
}

/// (3f) §11.6(c): `GeoPolygon`（rings[0]=外環・rings[1..]=穴）への包含判定。
/// 点 q は **外環の内側 かつ どの穴の内側でもない**とき内包される。
fn polygon_contains(poly: &umbra_geo::GeoPolygon, q: &umbra_geo::GeoPoint) -> bool {
    if poly.rings.is_empty() || !point_in_polygon(&poly.rings[0], q) {
        return false;
    }
    !poly.rings[1..].iter().any(|hole| point_in_polygon(hole, q))
}

/// M9 残(3) 3c-iii: 外環頂点 P が **terminator（日の出入り・ζ=0）上の半影縁点**かを判定する独立オラクル。
/// terminator 点は ζ=0 ゆえ、前方射影で **(i) |ζ|≈0** かつ **(ii) 軸からの面内距離 ≈ l1**
/// （ζ=0 で `|l1 − ζ·tan f1| = l1`）を同時に満たす。どのサンプル時刻由来か不定なので [P1,P4] の各時刻で試し、
/// いずれか 1 時刻で両条件を満たせば terminator 頂点とみなす。期待値 l1 は bessel 多項式（path とは独立な
/// 入力）から組む＝被テスト関数の戻りを流用しない。
///
/// 非対称性: 条件(i) の ζ≈0 を外すと昼面包絡点（ζ>0）が誤って terminator 扱いになるので、ζ_tol は十分小さく。
/// 条件(ii) を |l2| で測る変異・l1↔l2 取り違えは面内距離が桁違いに外れて偽になる。
fn vertex_is_terminator(
    p: &umbra_geo::GeoPoint,
    bessel: &BesselianPolynomial,
    sample_times: &[TtInstant],
    zeta_tol: f64,
    plane_tol: f64,
) -> bool {
    for t in sample_times {
        let e = match bessel.at(*t) {
            Ok(e) => e,
            Err(_) => continue,
        };
        let of = forward_project(p, &e);
        let in_plane = (of.xi - e.x).hypot(of.eta - e.y);
        // (i) ζ≈0（terminator 上）かつ (ii) 面内距離 ≈ l1（ζ=0 の半影半径）。
        if of.zeta.abs() < zeta_tol && (in_plane - e.l1).abs() < plane_tol {
            return true;
        }
    }
    false
}

/// 平面 (lon,lat) ray-casting による point-in-polygon（標準アルゴリズム・star-shaped を仮定しない）。
/// 点 q から +経度方向へ無限に伸ばした半直線が外環辺と交差する回数の偶奇で内外を判定する。
/// 経度は度・**反子午線非跨ぎ・非極**の合成/実 2024 を前提（§11.5・(3d) までの制約）。
/// リボン外環の位相（北南）を裏返す変異は包含を破るので非対称オラクルになる。
fn point_in_polygon(ring: &[umbra_geo::GeoPoint], q: &umbra_geo::GeoPoint) -> bool {
    let n = ring.len();
    let (qx, qy) = (lon_deg(q), lat_deg(q));
    let mut inside = false;
    let mut j = n - 1;
    for i in 0..n {
        let (xi, yi) = (lon_deg(&ring[i]), lat_deg(&ring[i]));
        let (xj, yj) = (lon_deg(&ring[j]), lat_deg(&ring[j]));
        // 辺 (i,j) が q の緯度 qy を跨ぎ、その交点経度が q の東（>qx）にあれば交差 1 回。
        let crosses = (yi > qy) != (yj > qy);
        if crosses {
            let x_at = xi + (qy - yi) / (yj - yi) * (xj - xi);
            if x_at > qx {
                inside = !inside;
            }
        }
        j = i;
    }
    inside
}

// ------------------------------------------------------------
// FAST: 存在（Some/None 分岐）
// ------------------------------------------------------------

/// FAST / 改訂((3f)): 部分食 phase（P1/P4=Some）＋include_limits=true で partial_limit=Some(GeoPolygon)・
/// `rings` は非空で `rings[0]`=外環が ≥3 頂点（非退化）。(3f) では**穴（rings[1..]）が存在しうる**ので
/// `rings.len()==1` は縛らず、外環の非退化と（穴があれば）各穴も ≥3 頂点であることを縛る。
///
/// 殺す変異: partial_limit を常に None にする・include_limits を無視・外環を空/2 点未満にする・
///   rings を 0 本にする・穴に退化リング（<3 頂点）を混ぜる。
#[test]
fn partial_phase_with_limits_produces_some_polygon() {
    let engine = standard_engine(bundled_time_data());
    let bessel = rigorous_bessel();
    let eclipse = partial_eclipse_with_bessel(bessel, 1.0, 1.5);

    let path = engine
        .path(&eclipse, PathOptions::default())
        .expect("部分食 phase の path() は成功する");

    let poly = path
        .partial_limit
        .as_ref()
        .expect("部分食 phase＋include_limits=true では partial_limit=Some");
    assert!(
        !poly.rings.is_empty(),
        "partial_limit は少なくとも外環 rings[0] を持つ"
    );
    assert!(
        poly.rings[0].len() >= 3,
        "外環は ≥3 頂点（非退化）, got {}",
        poly.rings[0].len()
    );
    for (h, hole) in poly.rings[1..].iter().enumerate() {
        assert!(
            hole.len() >= 3,
            "穴[{h}] も ≥3 頂点（退化リングを捏造しない）, got {}",
            hole.len()
        );
    }
}

/// FAST / 新規: include_limits=false なら partial_limit=None（部分食 phase があっても）。
///
/// 殺す変異: include_limits を無視して常に partial_limit を作る・フラグを反転して解釈する。
#[test]
fn partial_limit_none_when_include_limits_false() {
    let engine = standard_engine(bundled_time_data());
    let bessel = rigorous_bessel();
    let eclipse = partial_eclipse_with_bessel(bessel, 1.0, 1.5);

    let path = engine
        .path(
            &eclipse,
            PathOptions {
                include_limits: false,
                ..PathOptions::default()
            },
        )
        .expect("path() は成功する");
    assert!(
        path.partial_limit.is_none(),
        "include_limits=false では partial_limit=None"
    );
}

/// FAST / 新規: 部分食 phase 無し（partial_begin か partial_end が None）なら partial_limit=None。
/// 中心食でも P1/P4 が無ければ部分食域は組まない（partial_begin/end の && 条件）。
/// `central_eclipse`（partial_begin/end=None）で確認。
///
/// 殺す変異: P1/P4 を見ずに中心食で常に partial_limit を作る・片方だけ見て作る（|| 化）。
#[test]
fn partial_limit_none_without_partial_phase() {
    let engine = standard_engine(bundled_time_data());
    // central_eclipse は partial_begin/partial_end=None（中心食のみ・部分 phase 無し）。
    let eclipse = central_eclipse(geo(0.0, 0.0));
    // 前提を独立確認（部分 phase が無いこと）。
    assert!(
        eclipse.global.partial_begin.is_none() && eclipse.global.partial_end.is_none(),
        "central_eclipse は部分 phase 無し（P1/P4=None）"
    );

    let path = engine
        .path(&eclipse, PathOptions::default())
        .expect("中心食の path() は成功する");
    assert!(
        path.partial_limit.is_none(),
        "部分 phase 無し（P1/P4=None）では partial_limit=None"
    );
    // 中心線・限界線は中心食ゆえ Some のまま（partial_limit None が他を巻き込まない）。
    assert!(path.center_line.is_some(), "中心食なので center_line=Some");
}

/// FAST / 新規: 部分食 only（central_begin/end=None・P1/P4=Some）でも partial_limit=Some。
/// center_line は None（中心食でない）。部分食域は中心食と**独立**（center_line とは別経路）。
///
/// 殺す変異: partial_limit を center_line（中心食）と連動させる・部分食 only で partial_limit を
///   None にする・center_line=None のとき partial_limit も無条件 None にする。
#[test]
fn partial_only_eclipse_has_partial_limit_but_no_center_line() {
    let engine = standard_engine(bundled_time_data());
    let bessel = rigorous_bessel();
    let eclipse = partial_only_eclipse(bessel, 1.5);

    let path = engine
        .path(&eclipse, PathOptions::default())
        .expect("部分食 only の path() は成功する");

    assert!(
        path.partial_limit.is_some(),
        "部分食 only でも P1/P4=Some なら partial_limit=Some（中心食と独立）"
    );
    assert!(
        path.center_line.is_none(),
        "部分食 only では center_line=None（中心食でない）"
    );
    assert!(
        path.northern_limit.is_none() && path.southern_limit.is_none(),
        "部分食 only では本影南北限界線も None（中心食でない）"
    );
}

/// FAST / 新規: `sample_interval_seconds = 0.0` では `trace_penumbral_limits` が 1 サンプルのみ
/// 評価し、外環は北 1 点＋南 1 点 = 2 頂点（または 0）＝多角形を成さない（退化ガード `ring.len()<3`）。
/// よって partial_limit=None。
///
/// (3f) 更新: 3 領域ユニオンでも、1 サンプルでは帯も lune も頂点 <3 に退化し、ユニオン入力が全て捨てられる
/// （§11.6(b)「頂点 3 未満の領域は捨てる」・§11.7 入力正規化）ので結果は空＝`None`。
///
/// 殺す変異: 退化ガード（`ring.len() < 3` / ユニオンの「頂点 <3 のリングを捨てる」）を外す・
///   退化入力から面積 0 の多角形を捏造して `Some` を返す・1 サンプルでも 3 領域のどれかを無理に閉じる。
#[test]
fn partial_limit_none_when_single_sample_degenerate() {
    let engine = standard_engine(bundled_time_data());
    let bessel = rigorous_bessel();
    // partial_phase_with_limits_produces_some_polygon と同じ fixture（P1/P4=±1.5h ⊆ fit_interval ±2h）。
    let eclipse = partial_eclipse_with_bessel(bessel, 1.0, 1.5);

    let path = engine
        .path(
            &eclipse,
            PathOptions {
                sample_interval_seconds: 0.0,
                include_limits: true,
                // split_antimeridian は既定値（partial_limit=None の判定には無関係）。
                split_antimeridian: PathOptions::default().split_antimeridian,
            },
        )
        .expect("interval=0 でも path() は成功する");

    // interval=0 では center_line 等も 1 点になるが、本テストは partial_limit=None のみを縛る
    //（他フィールドは別テストの責務）。
    assert!(
        path.partial_limit.is_none(),
        "interval=0（1 サンプル）では外環頂点 < 3＝帯を成さず partial_limit=None。`ring.len()==3` 変異を撃つ"
    );
}

// ------------------------------------------------------------
// FAST: 頂点の正当性・方位ソート・非退化
// ------------------------------------------------------------

/// FAST / 改訂((3f) §11.6(d)・**頂点の正当性（緩和後）の主検証**): 外環の各頂点は、ある [P1,P4] サンプル
/// 時刻 t で **閉半影内**（面内距離 ≤ |l1 − ζ·tan f1| + tol・自己整合ζ）**かつ昼面側**（ζ ≥ −tol）にある。
/// 半影縁点は等号で、ユニオンが導入する**交点頂点**は折れ線離散化の範囲内でこれを満たす。
/// 期待値は bessel 多項式から独立に組む（被テスト関数の戻りを流用しない）。
///
/// **撤廃**: 「全頂点が半影縁条件を**厳密に**（等号で）満たす」は (3f) で成り立たない（§11.6(d)）ので
/// 等号契約はここで緩和した。代わりに「閉半影内＋昼面側」を全頂点に課す（捏造点・夜面点を撃つ）。
///
/// 殺す変異: 外環に閉半影の外の捏造点を入れる・夜面（ζ<0）の点を混ぜる・面内距離を |l2|（本影）で測る
///   （半影は本影の ~60 倍ゆえ大多数が閉半影外に落ちる）・|L1'| を中心軸 ζ₀ で測る（点自身の ζ でない）・
///   l1↔l2 や tan_f1↔tan_f2 の取り違え。
#[test]
fn partial_limit_vertices_are_inside_closed_penumbra_on_day_side() {
    let engine = standard_engine(bundled_time_data());
    let bessel = rigorous_bessel();
    let eclipse = partial_eclipse_with_bessel(bessel.clone(), 1.0, 1.5);

    let path = engine
        .path(&eclipse, PathOptions::default())
        .expect("部分食 phase の path() は成功する");
    let poly = path.partial_limit.as_ref().expect("partial_limit=Some");
    let ring = &poly.rings[0];

    // [P1,P4] のサンプル時刻列（build_partial_limit が収集に使う区間）。
    let p1 = eclipse.global.partial_begin.as_ref().unwrap().time_tt;
    let p4 = eclipse.global.partial_end.as_ref().unwrap().time_tt;
    let times = lockstep_sample_times(p1, p4, PathOptions::default().sample_interval_seconds);

    // 前方射影は厳密に閉じるが、頂点がどのサンプル時刻由来か不定なので各時刻で試す。
    // dist_tol=1e-6 [Re]（≈6 m）は半影縁点の等号一致（厳密 ~1e-7）＋交点頂点の離散化残差を包む合成域。
    // zeta_tol=1e-6 は「昼面側」の境界（terminator ζ=0 を許容しつつ夜面の捏造点を落とす）。
    for (j, p) in ring.iter().enumerate() {
        assert!(
            vertex_is_in_closed_penumbra(p, &bessel, &times, 1e-6, 1e-6),
            "外環頂点[{j}] (lat={}, lon={}) がどのサンプル時刻でも「閉半影内かつ昼面側」でない（捏造点/夜面点）",
            lat_deg(p),
            lon_deg(p)
        );
    }
    // 穴の頂点も同じ契約（ユニオンは入力境界上の点しか作らない）。
    for (h, hole) in poly.rings[1..].iter().enumerate() {
        for (j, p) in hole.iter().enumerate() {
            assert!(
                vertex_is_in_closed_penumbra(p, &bessel, &times, 1e-6, 1e-6),
                "穴[{h}] の頂点[{j}] (lat={}, lon={}) が「閉半影内かつ昼面側」でない（捏造点）",
                lat_deg(p),
                lon_deg(p)
            );
        }
    }
}

// ------------------------------------------------------------
// FAST: 単純性・穴構造（(3f) §11.6(c)(e) — リボン位相の撤廃に代わる位相契約）
// ------------------------------------------------------------

/// FAST / 新規((3f) §11.6(e)・**単純性の主検証**): 外環は**自己交差しない**（隣接しない辺同士が
/// 真に交差しない）。リボン位相（偶数頂点・前半北/後半南逆順・経度単調）は (3f) で撤廃されたので、
/// 位相の正しさはここ（単純性）・包含・頂点正当性で縛る。
///
/// 判定は平面 (lon,lat) の全辺ペア総当たり（端点共有の隣接辺は除外・端点接触/共線は交差としない）。
/// 反子午線非跨ぎ・非極の合成/実 2024 前提（(3d) までの制約）。
///
/// 殺す変異: ユニオンを行わず自己交差する帯をそのまま返す（(3f) 以前の帯は実データで自己交差する）・
///   環の再結合で断片を誤った順に繋ぐ（八の字）・穴を外環に連結して一本の自己交差リングにする。
#[test]
fn partial_limit_outer_ring_is_simple() {
    let engine = standard_engine(bundled_time_data());
    let bessel = rigorous_bessel();
    let eclipse = partial_eclipse_with_bessel(bessel, 1.0, 1.5);

    let path = engine
        .path(&eclipse, PathOptions::default())
        .expect("部分食 phase の path() は成功する");
    let poly = path.partial_limit.as_ref().expect("partial_limit=Some");

    if let Some((i, j)) = ring_self_intersection(&poly.rings[0]) {
        panic!(
            "外環が自己交差している（辺{i} × 辺{j}）: 単純多角形でない（ユニオン未実施/環再結合の誤り）。\
             ring={} 頂点",
            poly.rings[0].len()
        );
    }
    for (h, hole) in poly.rings[1..].iter().enumerate() {
        if let Some((i, j)) = ring_self_intersection(hole) {
            panic!("穴[{h}] が自己交差している（辺{i} × 辺{j}）");
        }
    }
}

/// FAST / 新規((3f) §11.6(c)・**穴構造の契約**): `partial_limit` は
/// `rings[0]`=外環・`rings[1..]`=穴の**単一多角形**（最大面積成分ひとつ）。穴は
///   (a) 外環より小さい面積、(b) 代表頂点が外環の内側、(c) 互いに入れ子でない（他の穴の内側にない）
/// を満たす。穴が 0 本でも契約を満たす（穴の存在は強制しない）。
///
/// 殺す変異: 2 番目以降の成分（別の外環）を穴として rings に混ぜる（外環の外に出る→(b) 破れ）・
///   外環と穴を取り違える（面積の大小が反転→(a) 破れ）・穴を重複して二重登録する（(c) 破れ）。
#[test]
fn partial_limit_holes_are_nested_inside_the_single_outer_ring() {
    let engine = standard_engine(bundled_time_data());
    let bessel = rigorous_bessel();
    let eclipse = partial_eclipse_with_bessel(bessel, 1.0, 1.5);

    let path = engine
        .path(&eclipse, PathOptions::default())
        .expect("部分食 phase の path() は成功する");
    let poly = path.partial_limit.as_ref().expect("partial_limit=Some");

    let outer_area = signed_area_deg2(&poly.rings[0]).abs();
    assert!(
        outer_area > 0.0,
        "外環の面積が 0（退化多角形を捏造している）"
    );
    for (h, hole) in poly.rings[1..].iter().enumerate() {
        let ha = signed_area_deg2(hole).abs();
        assert!(
            ha > 0.0 && ha < outer_area,
            "穴[{h}] の面積 {ha} は 0 より大きく外環 {outer_area} より小さい（外環/穴の取り違え）"
        );
        for (v, p) in hole.iter().enumerate() {
            assert!(
                point_in_polygon(&poly.rings[0], p),
                "穴[{h}] の頂点[{v}] (lat={}, lon={}) が外環の内側にない（別成分を穴として混入）",
                lat_deg(p),
                lon_deg(p)
            );
        }
        for (k, other) in poly.rings[1..].iter().enumerate() {
            if k == h {
                continue;
            }
            assert!(
                !hole.iter().all(|p| point_in_polygon(other, p)),
                "穴[{h}] が穴[{k}] の内側に入れ子になっている（穴の重複登録/割当誤り）"
            );
        }
    }
}

// ------------------------------------------------------------
// FAST: 包含（中心線点が外環内側かつ穴の外）— 平面 point-in-polygon
// ------------------------------------------------------------

/// FAST / 改訂((3f) §11.6(c)・**包含の主検証**): 中心線の各点が、**平面 (lon,lat) ray-casting
/// point-in-polygon**（star-shaped を仮定しない）で **外環の内側 かつ どの穴の内側でもない**
/// （partial ⊃ umbral path）。中心軸は半影域の内側を通るので、本影中心線は部分食域に内包される。
///
/// 注: `path.greatest_point` は合成メタデータの便宜値（geo(0,0)）で実際の半影帯（≈30–49°N）上に無いため
/// 包含判定の対象にしない。包含の本質は中心線（実 bessel 由来）が半影帯に入ること。
///
/// 殺す変異: リボンの北南を取り違える/逆順を欠いて自己交差させる（包含が崩れる）・外環を中心線より
///   内側に縮める・南北を取り違えて中心線が外に出る・外環頂点を捏造して領域が中心線を含まなくなる。
#[test]
fn partial_limit_contains_center_line() {
    let engine = standard_engine(bundled_time_data());
    let bessel = rigorous_bessel();
    let eclipse = partial_eclipse_with_bessel(bessel, 1.0, 1.5);

    let path = engine
        .path(&eclipse, PathOptions::default())
        .expect("部分食 phase の path() は成功する");
    let poly = path.partial_limit.as_ref().expect("partial_limit=Some");

    // 中心線の各点が「外環の内側 かつ どの穴の内側でもない」（partial ⊃ umbral path）。
    let center = path
        .center_line
        .as_ref()
        .expect("中心食なので center_line=Some");
    for (i, c) in center.points.iter().enumerate() {
        assert!(
            polygon_contains(poly, c),
            "中心線点[{i}] (lat={}, lon={}) が部分食域の外（外環の外 or 穴の中）＝partial ⊅ umbral path",
            lat_deg(c),
            lon_deg(c)
        );
    }
}

// ============================================================
// M9 残(3) 3c-iii: limb bulge 精緻化（terminator 連結）
//
// 確定仕様（docs/algorithms/11-path-partial-domain.md §11.4 (3c-iii)）:
//   各サンプル時刻で 北縁=昼面包絡 `solve_limit_edge(l1,+1)` が解ければそれ、ELSE terminator 交点
//   （`cone_terminator_intersections`・円∩terminator 楕円・ζ=0）の**高緯度側**。南縁=昼面包絡 `(−1)` ELSE
//   terminator 交点の**低緯度側**。北[i]/南[i] は lockstep。これで外環が [P1,P4] の時間端で terminator まで
//   張り出し、v1 リボンの limb 方向過小被覆を解消する。
//
// 観測契約（要確認#4 解決・§11.5）:
//   - FAST: 連結が**発火する合成 fixture**（`limb_continuation_bessel`＝影軸が縁寄りで端区間の昼面包絡が
//       無い）で、外環に **ζ≈0（terminator）かつ面内距離 ≈ l1** の頂点が現れる（前方射影で機械精度）。
//       v1（連結未実装）の外環は昼面包絡頂点のみ＝ζ≈0 頂点を含まず → RED。リボン不変条件（単一リング・
//       偶数頂点・北≥南 lockstep）も保つ。
//   - SLOW: 実 2024 で **terminator bulge が発火**（外環に ζ≈0 ＆ 面内距離 ≈ l1 の terminator 頂点が ≥1）し、
//       かつ**最大食まわりの核（中心線 ±10 サンプル窓）**が平面 point-in-polygon で内包される。
//       **中心線全点の内包は本スライスでは未達**（早朝端は deferred な 4 曲線 terminator-limb 境界が必要・要確認#4）。
// ============================================================

/// FAST / 新規（**3c-iii の主検証**・terminator 連結の発火）: 影軸が縁寄りの合成
/// （`limb_continuation_bessel`）では [P1,P4] の時間端で半影縁が昼面を外れ、外環は terminator まで
/// 張り出す。外環に **前方射影で ζ≈0（terminator 上）かつ軸からの面内距離 ≈ l1** の頂点が
/// **少なくとも 1 つ**現れる（terminator 連結が組み込まれた証拠・捏造でない＝半影縁条件を厳密に満たす）。
///
/// red（実装前）: v1 リボンは昼面の南北半影限界点（ζ>0）のみで terminator まで張り出さないため、ζ≈0 頂点が
///   存在せず本 assert が落ちる。
///
/// 殺す変異: 端区間で連結を行わず昼面包絡のみで外環を組む（ζ≈0 頂点が出ない）・terminator 交点を ζ≠0 で
///   捏造する（前方射影 ζ≈0 を外す）・面内距離を |l2| で測る/連結半径を本影 l2 にする（面内距離 ≈ l1 が破れる）。
#[test]
fn partial_limit_ring_includes_terminator_vertices_on_limb() {
    let engine = standard_engine(bundled_time_data());
    let bessel = limb_continuation_bessel();
    // [P1,P4]=epoch±1.5h。中央付近は昼面包絡が解ける（partial_limit=Some・非退化）一方、端区間は
    // 半影縁が昼面を外れ terminator 連結が必要。
    let eclipse = partial_eclipse_with_bessel(bessel.clone(), 1.0, 1.5);

    let path = engine
        .path(&eclipse, PathOptions::default())
        .expect("部分食 phase の path() は成功する");
    let poly = path
        .partial_limit
        .as_ref()
        .expect("limb 連結 fixture でも partial_limit=Some（中央は昼面包絡で種が出る）");
    let ring = &poly.rings[0];

    let p1 = eclipse.global.partial_begin.as_ref().unwrap().time_tt;
    let p4 = eclipse.global.partial_end.as_ref().unwrap().time_tt;
    let times = lockstep_sample_times(p1, p4, PathOptions::default().sample_interval_seconds);

    // 外環に terminator 頂点（ζ≈0・面内距離 ≈ l1）が ≥1 つある。前方射影は厳密に閉じる（ζ_tol 小・
    // plane_tol は半影縁の面内一致）。
    let terminator_vertices = ring
        .iter()
        .filter(|p| vertex_is_terminator(p, &bessel, &times, 1e-6, 1e-6))
        .count();
    assert!(
        terminator_vertices >= 1,
        "外環に terminator 頂点（前方射影 ζ≈0 ＆ 面内距離 ≈ l1）が無い＝limb 端で terminator 連結未発火 \
         （v1 リボンの limb 過小被覆・3c-iii 未実装）。ring={} 頂点",
        ring.len()
    );
}

/// FAST / 改訂((3f) §11.6(d)・**頂点正当性はユニオン後も維持**): limb fixture でも外環の各頂点が
/// 「閉半影内（面内距離 ≤ |l1 − ζ·tan f1| + tol）かつ昼面側（ζ ≥ −tol）」を満たす。terminator 頂点（ζ=0）も
/// 半影縁点も等号で、ユニオンの交点頂点は離散化の範囲内で満たす（捏造点が無い）。
///
/// **撤廃**: 厳密な等号（半影縁条件を全頂点が満たす）は (3f) で撤廃（交点頂点が入るため）。
///
/// 殺す変異: 連結/交点頂点を閉半影の外（半径 > |l1−ζ·tan f1|）に置く・terminator 交点を本影半径 l2 で解く・
///   夜面（ζ<0）の点を境界に混ぜる・昼面包絡頂点と terminator 頂点で半径式を取り違える。
#[test]
fn partial_limit_limb_vertices_are_inside_closed_penumbra_on_day_side() {
    let engine = standard_engine(bundled_time_data());
    let bessel = limb_continuation_bessel();
    let eclipse = partial_eclipse_with_bessel(bessel.clone(), 1.0, 1.5);

    let path = engine
        .path(&eclipse, PathOptions::default())
        .expect("部分食 phase の path() は成功する");
    let poly = path.partial_limit.as_ref().expect("partial_limit=Some");
    let ring = &poly.rings[0];

    let p1 = eclipse.global.partial_begin.as_ref().unwrap().time_tt;
    let p4 = eclipse.global.partial_end.as_ref().unwrap().time_tt;
    let times = lockstep_sample_times(p1, p4, PathOptions::default().sample_interval_seconds);

    for (j, p) in ring.iter().enumerate() {
        assert!(
            vertex_is_in_closed_penumbra(p, &bessel, &times, 1e-6, 1e-6),
            "limb fixture の外環頂点[{j}] (lat={}, lon={}) が「閉半影内かつ昼面側」でない（捏造点/夜面点）",
            lat_deg(p),
            lon_deg(p)
        );
    }
}

/// FAST / 新規((3f) §11.6(e) 置換・limb fixture の位相契約): limb fixture（terminator 連結が発火し
/// 帯が自己交差しうる fixture）でも、ユニオン後の外環は**単純**（自己交差なし）で、穴があれば
/// 外環の内側に入れ子になる。撤廃された「単一リング・偶数頂点・北南 lockstep 対」の代わりの契約。
///
/// 殺す変異: 自己交差する帯をユニオンせずそのまま返す・環再結合で断片を誤って繋ぐ（八の字）・
///   別成分を穴として混ぜる（外環の外に出る）。
#[test]
fn partial_limit_limb_polygon_is_simple_with_nested_holes() {
    let engine = standard_engine(bundled_time_data());
    let bessel = limb_continuation_bessel();
    let eclipse = partial_eclipse_with_bessel(bessel, 1.0, 1.5);

    let path = engine
        .path(&eclipse, PathOptions::default())
        .expect("部分食 phase の path() は成功する");
    let poly = path.partial_limit.as_ref().expect("partial_limit=Some");

    assert!(
        poly.rings[0].len() >= 3,
        "外環は ≥3 頂点, got {}",
        poly.rings[0].len()
    );
    if let Some((i, j)) = ring_self_intersection(&poly.rings[0]) {
        panic!(
            "limb fixture の外環が自己交差している（辺{i} × 辺{j}）: ユニオン未実施/環再結合の誤り。\
             ring={} 頂点",
            poly.rings[0].len()
        );
    }
    let outer_area = signed_area_deg2(&poly.rings[0]).abs();
    for (h, hole) in poly.rings[1..].iter().enumerate() {
        assert!(hole.len() >= 3, "穴[{h}] は ≥3 頂点, got {}", hole.len());
        if let Some((i, j)) = ring_self_intersection(hole) {
            panic!("穴[{h}] が自己交差している（辺{i} × 辺{j}）");
        }
        assert!(
            signed_area_deg2(hole).abs() < outer_area,
            "穴[{h}] の面積が外環以上（外環/穴の取り違え）"
        );
        for (v, p) in hole.iter().enumerate() {
            assert!(
                point_in_polygon(&poly.rings[0], p),
                "穴[{h}] の頂点[{v}] が外環の内側にない（別成分の混入）"
            );
        }
    }
}

// ------------------------------------------------------------
// SLOW: 実 2024-04-08 — partial_limit ballpark
// ------------------------------------------------------------

/// SLOW / 改訂((3f)・**headline acceptance＝中心線全点内包**・§11.6 / 受け入れテスト戦略):
/// 実エンジンで 2024-04-08 皆既を search → path()。partial_limit=Some・外環 ≥3 頂点・各頂点が妥当な緯度経度・
/// (a) 部分食域が皆既帯より緯度方向に広い（スパン > 中心線スパン・北端が中心線北端より外）・
/// (3) **実データで terminator 頂点が現れる**（外環に前方射影 ζ≈0 ＆ 軸からの面内距離 ≈ l1 の頂点が ≥1＝
///   terminator limb が境界に織り込まれた証拠）・
/// (4) **headline: 中心線の全点（194 点規模）が例外なく内包**される（外環の内側 かつ どの穴の内側でもない・
///   平面 (lon,lat) ray-casting point-in-polygon）・(5) 外環が**自己交差しない**（単純）。
/// NASA 緯度経度の直接一致は中心線位置精度律速ゆえ縛らず、桁の整合（広さ）＋位相（全点内包・単純性）で締める。
/// de440s 不要（解析暦）。
///
/// **(3f) での昇格**: 従来（3c-iii）は「最大食まわり ±10 サンプル窓だけ内包」という正直な妥協だった
/// （帯単独では U1/U4 近傍の中心線 4 点が帯の外に落ちる）。(3f) は部分食域を
/// **帯 ∪ morning lune ∪ evening lune** の多角形ユニオンとして組むため、morning/evening limb が境界に入り
/// **全点内包が達成されるべき契約**になる。よって窓限定を廃し、全点内包を headline として表明する。
///
/// 殺す変異: 実日食で partial_limit を None/捏造にする・外環を皆既帯より狭く縮める（(a) 破れ）・
///   morning lune / evening lune をユニオン入力から落とす（U1/U4 近傍の中心線点が外に落ち (4) 破れ）・
///   ユニオンを行わず帯だけを返す（(4)(5) 同時破れ＝実 2024 の帯は自己交差する）・
///   ユニオン結果の最大面積成分でなく別成分を返す（(4) 破れ）・穴を外環と取り違える（(4) 破れ）。
#[test]
fn real_2024_eclipse_partial_limit_is_plausible() {
    let engine = standard_engine(bundled_time_data());
    let range = umbra_core::TimeRange {
        start: utc(2024, 4, 8, 0, 0, 0.0),
        end: utc(2024, 4, 9, 0, 0, 0.0),
    };
    let eclipses = engine
        .search(range)
        .expect("2024-04-08 範囲の search は成功する");
    let eclipse = eclipses
        .iter()
        .find(|e| matches!(e.kind, SolarEclipseKind::Total))
        .expect("2024-04-08 皆既が見つかる");

    let path = engine
        .path(eclipse, PathOptions::default())
        .expect("実皆既の path() は成功する");

    let poly = path
        .partial_limit
        .as_ref()
        .expect("実 2024 は部分食 phase を持つので partial_limit=Some");
    let ring = &poly.rings[0];
    assert!(ring.len() >= 3, "実外環は ≥3 頂点, got {}", ring.len());
    for p in ring {
        assert!(
            lat_lon_in_range(p),
            "外環頂点が妥当な緯度経度域にない: lat={} lon={}",
            lat_deg(p),
            lon_deg(p)
        );
    }

    let center = path
        .center_line
        .as_ref()
        .expect("皆既なので center_line=Some");
    let greatest = &path.greatest_point;

    // (a) 部分食域は本影帯より緯度方向に広い（半影帯 ⊃ 皆既帯）。リボンの緯度スパンが中心線より大きく、
    //     北端は中心線北端より外（半影縁は皆既帯の外側へ張り出す）。v1 リボンは limb 方向に過小被覆
    //     （§11.4・要確認3）なので南端の厳密外側化までは縛らず、スパン優位と北端外側で「広い」を締める。
    let ring_max_lat = ring.iter().map(lat_deg).fold(f64::NEG_INFINITY, f64::max);
    let ring_min_lat = ring.iter().map(lat_deg).fold(f64::INFINITY, f64::min);
    let center_max_lat = center
        .points
        .iter()
        .map(lat_deg)
        .fold(f64::NEG_INFINITY, f64::max);
    let center_min_lat = center
        .points
        .iter()
        .map(lat_deg)
        .fold(f64::INFINITY, f64::min);
    assert!(
        ring_max_lat > center_max_lat,
        "部分食域の北端緯度 {ring_max_lat} が中心線北端 {center_max_lat} より北（半影帯は皆既帯より広い）"
    );
    assert!(
        (ring_max_lat - ring_min_lat) > (center_max_lat - center_min_lat),
        "部分食域の緯度スパン {} が中心線スパン {} より広い（半影帯 ⊃ 皆既帯）",
        ring_max_lat - ring_min_lat,
        center_max_lat - center_min_lat
    );

    // (3) **実データで terminator bulge が発火**: 外環に前方射影で ζ≈0（terminator 上）かつ軸からの面内距離 ≈ l1
    //     （ζ=0 の半影半径）の terminator 頂点が ≥1 つある。これは 3c-iii の端区間 terminator 連結が**合成 fixture
    //     だけでなく実 2024 でも発火した**証拠（v1 リボンは昼面包絡頂点 ζ>0 のみで terminator 頂点を持たない）。
    //     [P1,P4] のサンプル時刻列は FAST テストと同じ `lockstep_sample_times` で、P1/P4=partial_begin/end・
    //     bessel=eclipse.bessel から独立に組む（被テスト関数の戻りを流用しない）。
    //
    //     許容（zeta_tol/plane_tol）: FAST 合成は前方射影が厳密に閉じるので 1e-6 だが、実 2024 は l2/tan_f/d/μ の
    //     実暦評価＋中心線位置律速で round-trip がやや緩む。
    //       - plane_tol = 5e-4 [Re]（≈3 km）: 同ファイルの実日食用厳密 2 条件 `assert_exact_limit_conditions_real`
    //         で採用済みの面内距離許容と同値（実暦の半影縁面内一致のロバスト下限）。
    //       - zeta_tol = 1e-3: ζ=sin(太陽高度) なので 1e-3 は太陽高度 ≈0.057°。terminator（ζ=0）近傍だけを拾い、
    //         **昼面包絡頂点を誤分類しない**。実 2024 の昼面半影帯の典型 ζ は高度 数十° ＝ ζ~O(0.3–0.9) で、
    //         1e-3 とは 2–3 桁の隔たりがある（zeta_tol ≪ 帯の典型 ζ）。よって day-side ζ>0 頂点が (i) を満たすことはない。
    let p1 = eclipse.global.partial_begin.as_ref().unwrap().time_tt;
    let p4 = eclipse.global.partial_end.as_ref().unwrap().time_tt;
    let times = lockstep_sample_times(p1, p4, PathOptions::default().sample_interval_seconds);
    let terminator_vertices = ring
        .iter()
        .filter(|p| vertex_is_terminator(p, &eclipse.bessel, &times, 1e-3, 5e-4))
        .count();
    assert!(
        terminator_vertices >= 1,
        "実 2024 の外環に terminator 頂点（前方射影 ζ≈0 ＆ 面内距離 ≈ l1）が無い＝limb bulge が実データで未発火 \
         （3c-iii の terminator 連結が実 2024 で効いていない）。ring={} 頂点",
        ring.len()
    );

    // (4) **headline acceptance（(3f) 昇格）**: 中心線の**全点**が部分食域に内包される
    //     （外環の内側 かつ どの穴の内側でもない）。従来の「最大食まわり ±10 サンプル窓」から昇格。
    //     窓限定は帯単独（(3c-iii)）の限界に合わせた妥協であり、3 領域ユニオン後は例外を許さない。
    assert!(
        center.points.len() >= 100,
        "実 2024 の中心線は 194 点規模（サンプル列が痩せていると全点内包の意味が薄れる）, got {}",
        center.points.len()
    );
    let outside: Vec<usize> = center
        .points
        .iter()
        .enumerate()
        .filter(|(_, c)| !polygon_contains(poly, c))
        .map(|(i, _)| i)
        .collect();
    assert!(
        outside.is_empty(),
        "実 2024: 中心線 {} 点のうち {} 点が部分食域の外（partial ⊅ umbral path）。外れた index={:?} \
         先頭の点 (lat={}, lon={})",
        center.points.len(),
        outside.len(),
        outside,
        lat_deg(&center.points[outside[0]]),
        lon_deg(&center.points[outside[0]])
    );

    // (5) 外環は単純（自己交差なし）。実 2024 の帯は単独では自己交差する（§11.6(b)）ので、
    //     ユニオンを経ていない出力はここで落ちる。
    if let Some((i, j)) = ring_self_intersection(ring) {
        panic!(
            "実 2024: 外環が自己交差している（辺{i} × 辺{j}）＝ユニオン未実施/環再結合の誤り。ring={} 頂点",
            ring.len()
        );
    }

    // greatest_point は実データでは中心線上にあるので、同じ契約で内包される。
    assert!(
        polygon_contains(poly, greatest),
        "実 2024: 最大食点 (lat={}, lon={}) が部分食域の外",
        lat_deg(greatest),
        lon_deg(greatest)
    );
}

// ------------------------------------------------------------
// SLOW: 実 2016-03-09 皆既 — partial_limit の GeoJSON が反子午線で MultiPolygon に割れる
// ------------------------------------------------------------

/// `serde_json::Value` の座標ペア `[lon, lat]` を取り出す（型・長さも検証）。
fn geojson_coord_pair(v: &serde_json::Value) -> (f64, f64) {
    let arr = v.as_array().expect("座標ペアは配列");
    assert_eq!(arr.len(), 2, "座標ペアは [lon, lat] の長さ 2, got {arr:?}");
    (
        arr[0].as_f64().expect("lon は数値"),
        arr[1].as_f64().expect("lat は数値"),
    )
}

/// (lon, lat) 度平面の符号付き面積（shoelace・独立計算オラクル）。
fn geojson_ring_signed_area(coords: &[(f64, f64)]) -> f64 {
    let n = coords.len();
    if n < 3 {
        return 0.0;
    }
    let mut s = 0.0;
    for i in 0..n {
        let (x1, y1) = coords[i];
        let (x2, y2) = coords[(i + 1) % n];
        s += x1 * y2 - x2 * y1;
    }
    s / 2.0
}

/// SLOW / 新規（(3g)・§11.8 反子午線分割の**実データ結線**）: 実エンジンで **2016-03-09 皆既**（インドネシア〜
/// 太平洋）を search → `path(PathOptions::default())` し、`EclipsePath::to_geojson()` の
/// `role="partial_limit"` feature の geometry を検証する。
///
/// この日食を選ぶ理由（実測・2026-09-14）: 部分食域の外環が反子午線を **2 回**（＝偶数・真の跨ぎ）横切る
/// ＝ §11.8(b) の前提を満たす**極を囲まない**跨ぎ領域であり、(3g) の分割が実データ経路で発火する唯一の
/// 検証対象になる（実 2024-04-08 は跨がないので分割経路を通らない）。
/// 奇数回跨ぎ（＝極を囲む・§11.8(d) 未対応）の 2021-12-04 / 2028-07-22 / 2037-07-13 は対象外。
///
/// 縛るもの（**構造性質のみ**・外部座標表は使わない＝本プロジェクトのオラクル ハードコード禁止に従う）:
/// 1. geometry の `type` が **`"MultiPolygon"`**（跨ぐ領域が地球を一周する不正な単一 Polygon にならない）。
/// 2. 断片は **2 枚以上**、各多角形は **≥1 リング**。
/// 3. 全リングが **閉じている**（先頭==末尾）。
/// 4. 全リングの座標数が **≥4**（閉じた面のある環の最小＝3 頂点 + 先頭複製。退化断片を出さない）。
/// 5. 全頂点の経度が **[−180, 180]**・緯度が [−90, 90]（0..360 への付け替え・跨ぎ残りの検出）。
/// 6. 断片の **|符号付き面積| の総和 > 0**（分割で面積を失っていない）。
/// 7. 各断片の外環（rings[0]）が **CCW**（面積 > 0・RFC 7946 §3.1.6・§11.8(c) の分割後正規化）。
///
/// 殺す実装: 跨ぎリングを分割せず単一 Polygon のまま出す（1 破れ）・弧のペアリングに失敗して
///   片半球の弧を黙って捨てる（2 か 6 破れ）・子午線で閉じずに開いた弧を出す（3 破れ）・
///   経度を 0..360 へ付け替えて跨ぎを誤魔化す（5 破れ）・分割後に環向き正規化をしない（7 破れ）・
///   空/2 点の偽断片を生やす（4 破れ）・実データの unpaired 弧で panic する（テスト全体が落ちる）。
///
/// de440s 不要（解析暦）。
#[test]
fn real_2016_eclipse_partial_limit_geojson_splits_at_antimeridian() {
    let engine = standard_engine(bundled_time_data());
    let range = umbra_core::TimeRange {
        start: utc(2016, 3, 9, 0, 0, 0.0),
        end: utc(2016, 3, 10, 0, 0, 0.0),
    };
    let eclipses = engine
        .search(range)
        .expect("2016-03-09 範囲の search は成功する");
    let eclipse = eclipses
        .iter()
        .find(|e| matches!(e.kind, SolarEclipseKind::Total))
        .expect("2016-03-09 皆既が見つかる");

    let path = engine
        .path(eclipse, PathOptions::default())
        .expect("実皆既の path() は成功する");
    assert!(
        path.partial_limit.is_some(),
        "実 2016 は部分食 phase を持つので partial_limit=Some"
    );

    let json = path.to_geojson().expect("to_geojson は成功する");
    let fc: serde_json::Value = serde_json::from_str(&json).expect("to_geojson は妥当な JSON");
    let features = fc["features"].as_array().expect("features は配列");
    let geom = features
        .iter()
        .find(|f| f["properties"]["role"] == serde_json::Value::String("partial_limit".into()))
        .map(|f| &f["geometry"])
        .expect("role=partial_limit の feature がある");

    // 1. 反子午線を跨ぐので MultiPolygon。
    assert_eq!(
        geom["type"],
        serde_json::Value::String("MultiPolygon".into()),
        "実 2016 の部分食域は反子午線を跨ぐので MultiPolygon になる, got {}",
        geom["type"]
    );

    let polys = geom["coordinates"].as_array().expect("coordinates は配列");
    // 2. 断片は 2 枚以上。
    assert!(
        polys.len() >= 2,
        "跨ぎ分割の断片は 2 枚以上, got {}",
        polys.len()
    );

    let mut total_abs_area = 0.0_f64;
    for (pi, p) in polys.iter().enumerate() {
        let rings = p.as_array().expect("多角形はリング配列");
        assert!(
            !rings.is_empty(),
            "断片{pi} にリングが無い（空の多角形を捏造している）"
        );
        for (ri, r) in rings.iter().enumerate() {
            let coords: Vec<(f64, f64)> = r
                .as_array()
                .expect("リングは座標配列")
                .iter()
                .map(geojson_coord_pair)
                .collect();
            // 4. 退化しない（閉じた面のある環は ≥4 座標）。
            assert!(
                coords.len() >= 4,
                "断片{pi} リング{ri} の座標数が {} < 4（退化断片）",
                coords.len()
            );
            // 3. 閉じている。
            let (f_lon, f_lat) = coords[0];
            let (l_lon, l_lat) = coords[coords.len() - 1];
            assert!(
                (f_lon - l_lon).abs() < 1e-9 && (f_lat - l_lat).abs() < 1e-9,
                "断片{pi} リング{ri} が閉じていない: 先頭=[{f_lon},{f_lat}] 末尾=[{l_lon},{l_lat}]"
            );
            // 5. 座標域。
            for (ci, &(lon, lat)) in coords.iter().enumerate() {
                assert!(
                    (-180.0..=180.0).contains(&lon),
                    "断片{pi} リング{ri} 頂点{ci} の経度 {lon} が [−180,180] の外（0..360 付け替え/跨ぎ残り）"
                );
                assert!(
                    (-90.0..=90.0).contains(&lat),
                    "断片{pi} リング{ri} 頂点{ci} の緯度 {lat} が [−90,90] の外"
                );
            }
            let area = geojson_ring_signed_area(&coords);
            if ri == 0 {
                // 7. 外環は CCW（面積 > 0）。
                assert!(
                    area > 0.0,
                    "断片{pi} の外環が CCW でない（分割後の環向き正規化の欠落）: signed_area={area}"
                );
                total_abs_area += area.abs();
            }
        }
    }
    // 6. 面積が失われていない。
    assert!(
        total_abs_area > 0.0,
        "断片の外環面積の総和が 0（分割で領域を失っている）, got {total_abs_area}"
    );
}
