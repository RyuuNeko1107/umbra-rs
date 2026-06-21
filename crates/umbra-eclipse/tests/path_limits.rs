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
