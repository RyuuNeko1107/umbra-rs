//! `umbra_geo::union_rings` — 多角形ユニオン（§11.6 / §11.7 の**確定契約**）の受け入れテスト。
//!
//! 対象は `pub fn union_rings(rings: &[Vec<GeoPoint>]) -> Vec<GeoPolygon>`。
//! 内部アルゴリズム（アレンジメント方式・分割/量子化/再結合の手順）には**一切依存しない**。
//! 外から観測できる契約だけを縛る。
//!
//! ## 確定セマンティクス（テストで縛る正本）
//! - 座標系: `(lon, lat)` 度の**平面**。反子午線跨ぎ・極は未対応（テストも小さな領域のみ使う）。
//! - 充填: 各入力リングは **even-odd**。出力 U =「いずれかの入力リングの even-odd 内部」の境界。
//! - 入力: 閉/非閉どちらでも同結果。向き（CW/CCW）は結果に影響しない。頂点 3 未満は無視。
//! - 出力: `Vec<GeoPolygon>`、`rings[0]`=外環（CCW・shoelace>0）、`rings[1..]`=その多角形の穴（CW・<0）。
//!   各リングは**非閉**（先頭 != 末尾）。並びは**外環の面積降順**。ほぼ面積 0 の環は返さない。
//! - 空入力 / 全リング頂点 3 未満 / 全て退化 → 空 `Vec`。
//! - 退化: 自己交差（八の字）→分離 2 多角形 / 共有辺・共線重なり→1 本に畳む / 頂点上の交点→特別扱いなし /
//!   点接触→分離 2 多角形 / 完全重複→元と同じ 1 多角形 / 空洞を囲めば穴が付く。
//! - 数値許容: 座標一致は 1e-9 度程度（内部 1e-9 度スナップ）。
//!
//! ## テスト計画（工程1）
//! 1. **正常系**: 離れた 2 四角形 / 重なる 2 四角形（L 字・頂点と面積を厳密に）/ 3 リング /
//!    完全内包。
//! 2. **境界値**: 辺をちょうど共有（長方形に結合）/ 点接触 / 完全重複 / 一方の頂点が他方の辺上。
//! 3. **退化・異常系**: 空入力 / 頂点 0・1・2 のリング / 自己交差（八の字）単体 / 一直線（面積 0）。
//! 4. **不変条件**: 向き反転で同結果 / 閉・非閉で同結果 / 入力順入替で同結果（面積降順＝決定的なので
//!    出力を丸ごと比較できる）。
//! 5. **出力構造**: 外環 CCW・穴 CW・非閉・面積降順・連続重複頂点なし（共通ヘルパで全ケースに適用）。
//! 6. **点を捏造しない**: 出力の全頂点が、いずれかの入力辺（閉リングとして解釈）上に機械精度で乗る。
//! 7. **穴の生成**: 4 本のバーで正方形アニュラスを作り、穴 1 つを厳密に縛る。
//!
//! ## 期待される RED（実装前）
//! `union_rings` が未定義のため、`use umbra_geo::union_rings;` が解決できずコンパイル不能（E0432）。
//! これが想定どおりの赤であり、テスト側では解消しない。
//!
//! 全テストは FAST（実エンジン・IO・乱数・時刻に依存しない）。テスト間で状態を共有しない。

use umbra_geo::{union_rings, GeoPoint, GeoPolygon};

/// 座標一致の許容（内部スナップ 1e-9 度と同オーダー）。
const TOL: f64 = 1e-9;
/// 面積比較の許容（座標許容の伝播ぶんを見込む）。
const AREA_TOL: f64 = 1e-7;
/// 「入力辺の上に乗っている」判定の許容（スナップ 1e-9 度の伝播ぶん）。
const ON_EDGE_TOL: f64 = 1e-8;

// ============================================================
// ヘルパ（全て (lon, lat) 順で扱う。GeoPoint::from_degrees は (lat, lon) 順なので必ずここを通す）
// ============================================================

/// `(lon, lat)` 度から `GeoPoint` を作る。
fn p(lon: f64, lat: f64) -> GeoPoint {
    GeoPoint::from_degrees(lat, lon).expect("テストの緯度は全て有効範囲内")
}

/// `(lon, lat)` の並びからリング（点列）を作る。
fn ring(pts: &[(f64, f64)]) -> Vec<GeoPoint> {
    pts.iter().map(|&(lon, lat)| p(lon, lat)).collect()
}

/// 軸平行な長方形リング（CCW・非閉）。`(lon0,lat0)` が左下、`(lon1,lat1)` が右上。
fn rect(lon0: f64, lat0: f64, lon1: f64, lat1: f64) -> Vec<GeoPoint> {
    ring(&[(lon0, lat0), (lon1, lat0), (lon1, lat1), (lon0, lat1)])
}

/// `GeoPoint` を `(lon, lat)` の度タプルへ。
fn lonlat(pt: &GeoPoint) -> (f64, f64) {
    (pt.lon.degrees().0, pt.lat.degrees().0)
}

/// 非閉リングの符号付き面積（shoelace・(lon,lat) 平面・CCW>0）。末尾→先頭の辺も含める。
fn signed_area(r: &[GeoPoint]) -> f64 {
    let n = r.len();
    if n < 3 {
        return 0.0;
    }
    let mut sum = 0.0;
    for i in 0..n {
        let (x1, y1) = lonlat(&r[i]);
        let (x2, y2) = lonlat(&r[(i + 1) % n]);
        sum += x1 * y2 - x2 * y1;
    }
    sum / 2.0
}

/// 多角形の正味面積（外環 − 穴）。
fn net_area(poly: &GeoPolygon) -> f64 {
    let outer = signed_area(&poly.rings[0]).abs();
    let holes: f64 = poly.rings[1..].iter().map(|h| signed_area(h).abs()).sum();
    outer - holes
}

fn close(a: f64, b: f64, tol: f64) -> bool {
    (a - b).abs() < tol
}

/// 出力多角形の**構造契約**を全て検証する（全ケースで呼ぶ共通ガード）。
/// - 外環 CCW（面積>0）・穴 CW（面積<0）
/// - 各リングは非閉（先頭 != 末尾）・頂点 3 以上・連続重複頂点なし
/// - 面積ゼロの環を返さない
fn assert_polygon_structure(poly: &GeoPolygon) {
    assert!(!poly.rings.is_empty(), "多角形は最低 1 つの外環を持つ");
    for (index, r) in poly.rings.iter().enumerate() {
        assert!(
            r.len() >= 3,
            "リング{index}の頂点数は 3 以上（実際 {}）",
            r.len()
        );
        let (first, last) = (lonlat(&r[0]), lonlat(&r[r.len() - 1]));
        assert!(
            !(close(first.0, last.0, TOL) && close(first.1, last.1, TOL)),
            "リング{index}は非閉表現（先頭 {first:?} と末尾 {last:?} が一致してはならない）"
        );
        for i in 0..r.len() {
            let a = lonlat(&r[i]);
            let b = lonlat(&r[(i + 1) % r.len()]);
            assert!(
                !(close(a.0, b.0, TOL) && close(a.1, b.1, TOL)),
                "リング{index}に連続重複頂点 {a:?}"
            );
        }
        let area = signed_area(r);
        assert!(
            area.abs() > AREA_TOL,
            "リング{index}は面積ゼロの環であってはならない（面積 {area}）"
        );
        if index == 0 {
            assert!(area > 0.0, "外環は CCW（符号付き面積>0）だが {area}");
        } else {
            assert!(area < 0.0, "穴は CW（符号付き面積<0）だが {area}");
        }
    }
}

/// 出力全体の構造＋外環面積降順を検証する。
fn assert_output_structure(polys: &[GeoPolygon]) {
    for poly in polys {
        assert_polygon_structure(poly);
    }
    for w in polys.windows(2) {
        let a = signed_area(&w[0].rings[0]).abs();
        let b = signed_area(&w[1].rings[0]).abs();
        assert!(a >= b - AREA_TOL, "外環の面積降順が崩れている（{a} < {b}）");
    }
}

/// 点 `q` が線分 `a-b` 上にあるか（許容 `ON_EDGE_TOL`・端点含む）。
fn point_on_segment(q: (f64, f64), a: (f64, f64), b: (f64, f64)) -> bool {
    let (vx, vy) = (b.0 - a.0, b.1 - a.1);
    let (wx, wy) = (q.0 - a.0, q.1 - a.1);
    let len2 = vx * vx + vy * vy;
    if len2 <= 0.0 {
        return wx.hypot(wy) <= ON_EDGE_TOL;
    }
    // 線分上へ射影して [0,1] にクランプし、距離を測る（端点外は端点距離になる）。
    let t = ((wx * vx + wy * vy) / len2).clamp(0.0, 1.0);
    let (px, py) = (a.0 + t * vx, a.1 + t * vy);
    (q.0 - px).hypot(q.1 - py) <= ON_EDGE_TOL
}

/// **点を捏造しない**契約: 出力の全頂点が、いずれかの入力リングの辺（閉リングとして解釈）上に乗る。
/// 入力頂点も、入力 2 辺の交点も、必ず「いずれかの入力辺の上」にあるので、これで両方を包含する。
fn assert_no_fabricated_vertices(inputs: &[Vec<GeoPoint>], polys: &[GeoPolygon]) {
    // 入力辺（閉リング。頂点 3 未満のリングは union の入力として無視されるが、辺集合にも寄与しない）。
    let mut edges: Vec<((f64, f64), (f64, f64))> = Vec::new();
    for r in inputs {
        if r.len() < 3 {
            continue;
        }
        let coords: Vec<(f64, f64)> = r.iter().map(lonlat).collect();
        for i in 0..coords.len() {
            let a = coords[i];
            let b = coords[(i + 1) % coords.len()];
            edges.push((a, b));
        }
    }
    for poly in polys {
        for r in &poly.rings {
            for pt in r {
                let q = lonlat(pt);
                assert!(
                    edges.iter().any(|&(a, b)| point_on_segment(q, a, b)),
                    "出力頂点 {q:?} はどの入力辺の上にも乗っていない（点の捏造）"
                );
            }
        }
    }
}

/// リングが期待頂点列（`(lon,lat)` 順）と**巡回一致**するか（回転は許すが順序と向きは固定）。
/// 出力の開始頂点を実装詳細として縛らないための比較。
fn assert_ring_matches(actual: &[GeoPoint], expected: &[(f64, f64)]) {
    let got: Vec<(f64, f64)> = actual.iter().map(lonlat).collect();
    assert_eq!(
        got.len(),
        expected.len(),
        "頂点数が違う: got {got:?} / expected {expected:?}"
    );
    let n = got.len();
    let matched = (0..n).any(|shift| {
        (0..n).all(|i| {
            let g = got[(i + shift) % n];
            let e = expected[i];
            close(g.0, e.0, TOL) && close(g.1, e.1, TOL)
        })
    });
    assert!(
        matched,
        "リングが期待列と巡回一致しない: got {got:?} / expected {expected:?}"
    );
}

/// 出力を比較可能な正規形（各リングを最小頂点から回転した `(lon,lat)` 列）へ。
/// 不変条件テスト（向き反転・閉/非閉・入力順入替）で結果同一性を丸ごと比較するために使う。
fn canonical(polys: &[GeoPolygon]) -> Vec<Vec<Vec<(i64, i64)>>> {
    // 1e-9 度グリッドの整数化で、許容つき比較を「完全一致比較」に落とす。
    // `round()` 済みの値に対する変換で、テスト座標は |値| < 180 ゆえ 1.8e11 と i64 の範囲内。
    #[allow(clippy::cast_possible_truncation)]
    fn q(v: f64) -> i64 {
        (v / TOL).round() as i64
    }
    polys
        .iter()
        .map(|poly| {
            poly.rings
                .iter()
                .map(|r| {
                    let coords: Vec<(i64, i64)> = r
                        .iter()
                        .map(|pt| {
                            let (lon, lat) = lonlat(pt);
                            (q(lon), q(lat))
                        })
                        .collect();
                    // 最小頂点を先頭に回転（開始頂点の違いを吸収。向き・順序は保つ）。
                    let n = coords.len();
                    let start = (0..n).min_by_key(|&i| coords[i]).unwrap_or(0);
                    (0..n).map(|i| coords[(i + start) % n]).collect()
                })
                .collect()
        })
        .collect()
}

// ============================================================
// 正常系
// ============================================================

/// 離れた 2 つの四角形は、2 つの独立多角形として**面積降順**で返る。
/// 小: (0,0)-(1,1)（面積 1）、大: (3,0)-(5,2)（面積 4）。
///
/// 殺す変異: 常に 1 多角形へ畳む（成分分離の脱落）・面積昇順や入力順のまま返す（ソート脱落/反転）・
/// 外環頂点の取り違え・穴を捏造する。
#[test]
fn union_disjoint_squares_returns_two_components_sorted_by_area_desc() {
    let small = rect(0.0, 0.0, 1.0, 1.0);
    let big = rect(3.0, 0.0, 5.0, 2.0);
    let out = union_rings(&[small.clone(), big.clone()]);

    assert_eq!(out.len(), 2, "離れた 2 領域は 2 多角形");
    assert_output_structure(&out);
    assert_no_fabricated_vertices(&[small, big], &out);

    // 先頭＝大きい方（面積 4）。
    assert_eq!(out[0].rings.len(), 1, "穴は無い");
    assert!(
        close(net_area(&out[0]), 4.0, AREA_TOL),
        "{}",
        net_area(&out[0])
    );
    assert_ring_matches(
        &out[0].rings[0],
        &[(3.0, 0.0), (5.0, 0.0), (5.0, 2.0), (3.0, 2.0)],
    );

    // 2 番目＝小さい方（面積 1）。
    assert_eq!(out[1].rings.len(), 1, "穴は無い");
    assert!(
        close(net_area(&out[1]), 1.0, AREA_TOL),
        "{}",
        net_area(&out[1])
    );
    assert_ring_matches(
        &out[1].rings[0],
        &[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)],
    );
}

/// 角で重なる 2 四角形 A=(0,0)-(2,2)・B=(1,1)-(3,3) の和は L 字の 8 頂点多角形。
/// 頂点列と面積（4+4−1=7）を**機械精度で厳密に**縛る。
///
/// 殺す変異: 交点 (2,1)/(1,2) を落とす・重なり分を二重計上（面積 8）や欠損（面積 6）・
/// 重なり領域を差し引いてしまう・頂点順序を逆転（CW 出力）・余計な頂点を挿入する。
#[test]
fn union_overlapping_squares_forms_exact_l_shape() {
    let a = rect(0.0, 0.0, 2.0, 2.0);
    let b = rect(1.0, 1.0, 3.0, 3.0);
    let out = union_rings(&[a.clone(), b.clone()]);

    assert_eq!(out.len(), 1, "重なるので単一成分");
    assert_output_structure(&out);
    assert_no_fabricated_vertices(&[a, b], &out);
    assert_eq!(out[0].rings.len(), 1, "穴は無い");
    assert!(
        close(net_area(&out[0]), 7.0, AREA_TOL),
        "{}",
        net_area(&out[0])
    );
    // CCW の L 字（交点は (2,1) と (1,2)）。
    assert_ring_matches(
        &out[0].rings[0],
        &[
            (0.0, 0.0),
            (2.0, 0.0),
            (2.0, 1.0),
            (3.0, 1.0),
            (3.0, 3.0),
            (1.0, 3.0),
            (1.0, 2.0),
            (0.0, 2.0),
        ],
    );
}

/// 3 つの四角形が鎖状に重なる場合、単一の長方形 (0,0)-(4,2) に融合する。
/// A=(0,0)-(2,2), B=(1,0)-(3,2), C=(2,0)-(4,2)。面積 8・バウンディングボックス・
/// 全頂点が長方形の境界上にあることを縛る（共線頂点が残るか否かは仕様外なので数は縛らない）。
///
/// 殺す変異: 2 リングだけ処理して 3 本目を無視する・重なりの二重計上（面積>8）・
/// 内部に埋もれた辺を境界として残す（境界外の頂点が出る）。
#[test]
fn union_three_overlapping_rings_merges_into_one_rectangle() {
    let a = rect(0.0, 0.0, 2.0, 2.0);
    let b = rect(1.0, 0.0, 3.0, 2.0);
    let c = rect(2.0, 0.0, 4.0, 2.0);
    let inputs = vec![a, b, c];
    let out = union_rings(&inputs);

    assert_eq!(out.len(), 1, "鎖状に重なるので単一成分");
    assert_output_structure(&out);
    assert_no_fabricated_vertices(&inputs, &out);
    assert_eq!(out[0].rings.len(), 1, "穴は無い");
    assert!(
        close(net_area(&out[0]), 8.0, AREA_TOL),
        "{}",
        net_area(&out[0])
    );

    // 全頂点が長方形 (0,0)-(4,2) の境界上（内部に取り残された頂点が無い）。
    for pt in &out[0].rings[0] {
        let (lon, lat) = lonlat(pt);
        let on_vertical = close(lon, 0.0, TOL) || close(lon, 4.0, TOL);
        let on_horizontal = close(lat, 0.0, TOL) || close(lat, 2.0, TOL);
        assert!(
            (on_vertical && (-TOL..=2.0 + TOL).contains(&lat))
                || (on_horizontal && (-TOL..=4.0 + TOL).contains(&lon)),
            "頂点 ({lon},{lat}) が長方形 (0,0)-(4,2) の境界上にない"
        );
    }
    // 4 隅は必ず頂点として存在する。
    for corner in [(0.0, 0.0), (4.0, 0.0), (4.0, 2.0), (0.0, 2.0)] {
        assert!(
            out[0].rings[0]
                .iter()
                .any(|pt| close(lonlat(pt).0, corner.0, TOL) && close(lonlat(pt).1, corner.1, TOL)),
            "隅 {corner:?} が出力頂点に無い"
        );
    }
}

/// 一方が他方を完全に内包する場合、結果は外側のリングそのもの（内側は消える）。
/// 外 (0,0)-(10,10)・内 (2,2)-(4,4)。穴は**できない**（和なので内側は埋まる）。
///
/// 殺す変異: 内包リングを穴として扱う（rings.len()==2）・交差（内側だけ返す）にすり替える・
/// 内側の頂点を外環に混ぜる・面積 100−4=96 になる（差集合化）。
#[test]
fn union_fully_contained_ring_yields_outer_only_without_hole() {
    let outer = rect(0.0, 0.0, 10.0, 10.0);
    let inner = rect(2.0, 2.0, 4.0, 4.0);
    let inputs = vec![outer, inner];
    let out = union_rings(&inputs);

    assert_eq!(out.len(), 1);
    assert_output_structure(&out);
    assert_no_fabricated_vertices(&inputs, &out);
    assert_eq!(out[0].rings.len(), 1, "内包リングは穴にならない");
    assert!(
        close(net_area(&out[0]), 100.0, AREA_TOL),
        "{}",
        net_area(&out[0])
    );
    assert_ring_matches(
        &out[0].rings[0],
        &[(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)],
    );
}

// ============================================================
// 境界値
// ============================================================

/// 辺をちょうど共有する 2 四角形 (0,0)-(1,1) と (1,0)-(2,1) は、1 つの長方形 (0,0)-(2,1) に結合する。
/// 共有辺 lon=1 は結果に**現れない**（重複辺が畳まれる）。
///
/// 殺す変異: 共有辺を境界として残す（2 成分・面積 2 だが自己接触リング）・
/// 共有辺の分だけ内部に線が残り頂点 (1,0)/(1,1) が内部辺として二重に出る・
/// 接触を「重なり無し」と見なして 2 多角形を返す。
#[test]
fn union_edge_sharing_squares_merge_into_single_rectangle() {
    let left = rect(0.0, 0.0, 1.0, 1.0);
    let right = rect(1.0, 0.0, 2.0, 1.0);
    let inputs = vec![left, right];
    let out = union_rings(&inputs);

    assert_eq!(out.len(), 1, "辺共有は内部が繋がるので単一多角形");
    assert_output_structure(&out);
    assert_no_fabricated_vertices(&inputs, &out);
    assert_eq!(out[0].rings.len(), 1, "穴は無い");
    assert!(
        close(net_area(&out[0]), 2.0, AREA_TOL),
        "{}",
        net_area(&out[0])
    );

    // 全頂点が長方形 (0,0)-(2,1) の境界上、かつ同一頂点の重複出現なし（重複辺が畳まれた証拠）。
    let coords: Vec<(f64, f64)> = out[0].rings[0].iter().map(lonlat).collect();
    for &(lon, lat) in &coords {
        let on_v = close(lon, 0.0, TOL) || close(lon, 2.0, TOL);
        let on_h = close(lat, 0.0, TOL) || close(lat, 1.0, TOL);
        assert!(
            (on_v && (-TOL..=1.0 + TOL).contains(&lat))
                || (on_h && (-TOL..=2.0 + TOL).contains(&lon)),
            "頂点 ({lon},{lat}) が長方形 (0,0)-(2,1) の境界上にない（共有辺が残っている疑い）"
        );
    }
    for i in 0..coords.len() {
        for j in (i + 1)..coords.len() {
            assert!(
                !(close(coords[i].0, coords[j].0, TOL) && close(coords[i].1, coords[j].1, TOL)),
                "同一頂点 {:?} が重複して出現（重複辺が畳まれていない）",
                coords[i]
            );
        }
    }
}

/// 1 頂点だけを共有する 2 四角形 (0,0)-(1,1) と (1,1)-(2,2) は、**分離した 2 多角形**になる
/// （内部が繋がっていない）。各面積 1。
///
/// 殺す変異: 点接触を「繋がっている」と見なして 1 つの自己接触多角形にする・
/// 一方を捨てる（1 多角形）・接触点を交点として余計な頂点を増やす。
#[test]
fn union_point_touching_squares_stay_two_polygons() {
    let a = rect(0.0, 0.0, 1.0, 1.0);
    let b = rect(1.0, 1.0, 2.0, 2.0);
    let inputs = vec![a, b];
    let out = union_rings(&inputs);

    assert_eq!(out.len(), 2, "点接触は内部が繋がらないので 2 多角形");
    assert_output_structure(&out);
    assert_no_fabricated_vertices(&inputs, &out);
    for poly in &out {
        assert_eq!(poly.rings.len(), 1);
        assert_eq!(poly.rings[0].len(), 4, "各成分は 4 頂点の四角形");
        assert!(close(net_area(poly), 1.0, AREA_TOL), "{}", net_area(poly));
    }
    // 面積が同値なので、2 成分が別位置であることを外接ボックスで確認する。
    let lon_max_0 = out[0].rings[0]
        .iter()
        .map(|pt| lonlat(pt).0)
        .fold(f64::NEG_INFINITY, f64::max);
    let lon_max_1 = out[1].rings[0]
        .iter()
        .map(|pt| lonlat(pt).0)
        .fold(f64::NEG_INFINITY, f64::max);
    assert!(
        close(lon_max_0.min(lon_max_1), 1.0, TOL) && close(lon_max_0.max(lon_max_1), 2.0, TOL),
        "2 成分は (0,0)-(1,1) と (1,1)-(2,2) のはず（lon_max = {lon_max_0}, {lon_max_1}）"
    );
}

/// 完全に同じ四角形を 2 枚重ねた入力は、元と同じ 1 多角形を返す（重複辺が畳まれる）。
///
/// 殺す変異: 2 枚を別成分として返す（len==2）・even-odd を入力全体に適用して空にする
/// （2 枚重ね＝偶数 → 穴、は誤り。合併の membership は「いずれかのリングの内部」）・
/// 面積が 2 倍/0 になる・頂点が 8 個に増える。
#[test]
fn union_duplicate_rings_returns_the_same_single_polygon() {
    let square = rect(0.0, 0.0, 1.0, 1.0);
    let inputs = vec![square.clone(), square.clone()];
    let out = union_rings(&inputs);

    assert_eq!(out.len(), 1, "完全重複は 1 多角形");
    assert_output_structure(&out);
    assert_no_fabricated_vertices(&inputs, &out);
    assert_eq!(out[0].rings.len(), 1, "穴は無い");
    assert!(
        close(net_area(&out[0]), 1.0, AREA_TOL),
        "{}",
        net_area(&out[0])
    );
    assert_ring_matches(
        &out[0].rings[0],
        &[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)],
    );
}

/// 一方の頂点が他方の辺の**途中**に乗り、かつ辺の一部が共線で重なるケース。
/// A=(0,0)-(2,2)、B=(1,2)-(3,4)。B の下辺 lat=2 は A の上辺 lat=2 と区間 [1,2] で共線重複し、
/// B の頂点 (1,2) は A の上辺の内部に乗る。和は 8 頂点・面積 8。
///
/// 殺す変異: 頂点上の交点を特別扱いして落とす（(1,2)/(2,2) の欠落）・共線重なりを分割点として
/// 登録せず内部に線を残す・接触のみと誤判定して 2 多角形にする・面積が 4 や 7 になる。
#[test]
fn union_vertex_on_other_edge_with_collinear_overlap() {
    let a = rect(0.0, 0.0, 2.0, 2.0);
    let b = rect(1.0, 2.0, 3.0, 4.0);
    let inputs = vec![a, b];
    let out = union_rings(&inputs);

    assert_eq!(out.len(), 1, "共線区間で内部が繋がるので単一多角形");
    assert_output_structure(&out);
    assert_no_fabricated_vertices(&inputs, &out);
    assert_eq!(out[0].rings.len(), 1, "穴は無い");
    assert!(
        close(net_area(&out[0]), 8.0, AREA_TOL),
        "{}",
        net_area(&out[0])
    );
    assert_ring_matches(
        &out[0].rings[0],
        &[
            (0.0, 0.0),
            (2.0, 0.0),
            (2.0, 2.0),
            (3.0, 2.0),
            (3.0, 4.0),
            (1.0, 4.0),
            (1.0, 2.0),
            (0.0, 2.0),
        ],
    );
}

// ============================================================
// 退化・異常系
// ============================================================

/// 空入力は空 `Vec`（panic しない・多角形を捏造しない）。
///
/// 殺す変異: 空入力で panic / unwrap する・空の `GeoPolygon` を 1 つ返す。
#[test]
fn union_empty_input_returns_empty_vec() {
    let out = union_rings(&[]);
    assert!(out.is_empty(), "空入力は空 Vec（実際 {} 件）", out.len());
}

/// 頂点 0・1・2 のリングだけの入力は無視され、空 `Vec` になる。
///
/// 殺す変異: 頂点 <3 の足切り条件を `<2` / `<=3` にずらす・退行リングで index out of bounds panic・
/// 2 点リングから面積ゼロの環を返す。
#[test]
fn union_rings_with_fewer_than_three_vertices_are_ignored() {
    let zero: Vec<GeoPoint> = Vec::new();
    let one = ring(&[(0.0, 0.0)]);
    let two = ring(&[(0.0, 0.0), (1.0, 1.0)]);
    let out = union_rings(&[zero, one, two]);
    assert!(
        out.is_empty(),
        "頂点 3 未満のみの入力は空 Vec（実際 {} 件）",
        out.len()
    );
}

/// 頂点 3 未満のリングが混ざっても、有効なリングの結果は変わらない。
///
/// 殺す変異: 退行リングに引きずられて全体を空にする・退行リングの点を出力に混ぜる。
#[test]
fn union_ignores_degenerate_rings_mixed_with_a_valid_ring() {
    let square = rect(0.0, 0.0, 1.0, 1.0);
    let out = union_rings(&[
        ring(&[(5.0, 5.0)]),
        square.clone(),
        ring(&[(7.0, 7.0), (8.0, 8.0)]),
        Vec::new(),
    ]);

    assert_eq!(out.len(), 1);
    assert_output_structure(&out);
    assert_eq!(out[0].rings.len(), 1);
    assert_ring_matches(
        &out[0].rings[0],
        &[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)],
    );
    // 退行リングの座標は一切現れない。
    for pt in &out[0].rings[0] {
        let (lon, _) = lonlat(pt);
        assert!(lon < 2.0, "退行リング由来の頂点が混入している（lon={lon}）");
    }
}

/// 自己交差リング（八の字／蝶ネクタイ）単体は、even-odd により**分離した 2 多角形**になる。
/// リング `[(0,0),(2,0),(0,2),(2,2)]` の 2 本の対角線は (1,1) で交差し、
/// 下三角 {(0,0),(2,0),(1,1)}（面積 1）と上三角 {(0,2),(2,2),(1,1)}（面積 1）に分かれる。
///
/// 殺す変異: 自己交差を分割せずそのまま 1 リングで返す（面積 0 や自己交差多角形）・
/// 交点 (1,1) を打たない・片方の三角形だけ返す・even-odd でなく nonzero 充填にする。
#[test]
fn union_self_intersecting_figure_eight_splits_into_two_triangles() {
    let bowtie = ring(&[(0.0, 0.0), (2.0, 0.0), (0.0, 2.0), (2.0, 2.0)]);
    let inputs = vec![bowtie];
    let out = union_rings(&inputs);

    assert_eq!(out.len(), 2, "八の字は even-odd で 2 領域");
    assert_output_structure(&out);
    assert_no_fabricated_vertices(&inputs, &out);
    for poly in &out {
        assert_eq!(poly.rings.len(), 1, "穴は無い");
        assert_eq!(poly.rings[0].len(), 3, "各成分は三角形");
        assert!(close(net_area(poly), 1.0, AREA_TOL), "{}", net_area(poly));
        // どの成分も交点 (1,1) を頂点に持つ。
        assert!(
            poly.rings[0]
                .iter()
                .any(|pt| close(lonlat(pt).0, 1.0, TOL) && close(lonlat(pt).1, 1.0, TOL)),
            "自己交点 (1,1) が頂点に無い"
        );
    }
    // 2 成分は上下に分かれる（下三角の lat 最大 = 1、上三角の lat 最小 = 1）。
    let lat_mins: Vec<f64> = out
        .iter()
        .map(|poly| {
            poly.rings[0]
                .iter()
                .map(|pt| lonlat(pt).1)
                .fold(f64::INFINITY, f64::min)
        })
        .collect();
    assert!(
        close(
            lat_mins.iter().cloned().fold(f64::INFINITY, f64::min),
            0.0,
            TOL
        ) && close(
            lat_mins.iter().cloned().fold(f64::NEG_INFINITY, f64::max),
            1.0,
            TOL
        ),
        "2 成分は下三角(lat_min=0)・上三角(lat_min=1) のはず（{lat_mins:?}）"
    );
}

/// 一直線に並んだ面積ゼロのリングは、環を返さない（捏造しない）。
///
/// 殺す変異: 面積ゼロ足切りの削除（退化した環が返る）・|A| 比較の符号ミス・
/// 一直線入力で panic する。
#[test]
fn union_zero_area_collinear_ring_returns_empty_vec() {
    let line = ring(&[(0.0, 0.0), (1.0, 0.0), (2.0, 0.0), (1.0, 0.0)]);
    let out = union_rings(&[line]);
    assert!(
        out.is_empty(),
        "面積ゼロのリングは環を返さない（実際 {} 件）",
        out.len()
    );
}

// ============================================================
// 不変条件
// ============================================================

/// 入力リングの**向き（CW/CCW）を反転**しても結果は同一。
/// L 字ケース（重なる 2 四角形）で、両方 CCW / 両方 CW / 片方ずつ の 4 通りを比較する。
///
/// 殺す変異: 入力の向きで membership 判定を変える（nonzero 的な扱い）・
/// 出力の向き正規化を入力向きに委ねる・CW 入力で空を返す。
#[test]
fn union_is_invariant_to_input_ring_orientation() {
    let a_ccw = rect(0.0, 0.0, 2.0, 2.0);
    let b_ccw = rect(1.0, 1.0, 3.0, 3.0);
    let mut a_cw = a_ccw.clone();
    a_cw.reverse();
    let mut b_cw = b_ccw.clone();
    b_cw.reverse();
    // 前提の確認（ヘルパが本当に向きを変えていること）。
    assert!(signed_area(&a_ccw) > 0.0 && signed_area(&a_cw) < 0.0);

    let base = canonical(&union_rings(&[a_ccw.clone(), b_ccw.clone()]));
    assert_eq!(
        base,
        canonical(&union_rings(&[a_cw.clone(), b_ccw.clone()]))
    );
    assert_eq!(
        base,
        canonical(&union_rings(&[a_ccw.clone(), b_cw.clone()]))
    );
    assert_eq!(base, canonical(&union_rings(&[a_cw, b_cw])));
}

/// 入力リングを**閉じても閉じなくても**結果は同一（末尾==先頭の重複は落とされる）。
///
/// 殺す変異: 閉リングの重複点を落とさず連続重複頂点や面積ゼロ辺として残す・
/// 非閉リングの「末尾→先頭」辺を張り忘れる（開いた領域として扱う）。
#[test]
fn union_is_invariant_to_closed_or_open_input_rings() {
    let a_open = rect(0.0, 0.0, 2.0, 2.0);
    let b_open = rect(1.0, 1.0, 3.0, 3.0);
    let mut a_closed = a_open.clone();
    a_closed.push(a_open[0]);
    let mut b_closed = b_open.clone();
    b_closed.push(b_open[0]);

    let base = canonical(&union_rings(&[a_open.clone(), b_open.clone()]));
    assert!(!base.is_empty(), "前提: 非閉入力でも結果は非空");
    assert_eq!(
        base,
        canonical(&union_rings(&[a_closed.clone(), b_closed.clone()]))
    );
    assert_eq!(base, canonical(&union_rings(&[a_closed, b_open])));
    assert_eq!(base, canonical(&union_rings(&[a_open, b_closed])));
}

/// 入力リングの**順序**を入れ替えても結果は同一（出力は面積降順で決定的）。
/// 3 リング（面積がすべて異なる離散＋重なり混在）で全 6 順列を比較する。
///
/// 殺す変異: 出力順が入力順に依存する（ソート脱落）・最初のリングを基準に据える実装・
/// 同面積タイでの非決定な並び。
#[test]
fn union_is_invariant_to_input_ring_order() {
    let a = rect(0.0, 0.0, 2.0, 2.0); // 面積 4（b と重なる）
    let b = rect(1.0, 1.0, 3.0, 3.0); // 面積 4
    let c = rect(10.0, 10.0, 11.0, 11.0); // 面積 1（離れている）
    let base = canonical(&union_rings(&[a.clone(), b.clone(), c.clone()]));
    assert_eq!(base.len(), 2, "前提: L 字 + 離れた小四角 = 2 成分");

    let perms = [
        [&a, &c, &b],
        [&b, &a, &c],
        [&b, &c, &a],
        [&c, &a, &b],
        [&c, &b, &a],
    ];
    for perm in perms {
        let input: Vec<Vec<GeoPoint>> = perm.iter().map(|r| (*r).clone()).collect();
        assert_eq!(
            base,
            canonical(&union_rings(&input)),
            "入力順で結果が変わる"
        );
    }
}

// ============================================================
// 穴の生成
// ============================================================

/// 4 本の細長いバーで正方形のアニュラスを作ると、穴が 1 つできる。
/// 下 (0,0)-(10,1)・上 (0,9)-(10,10)・左 (0,0)-(1,10)・右 (9,0)-(10,10)。
/// 外環 = (0,0)-(10,10) の 4 頂点（面積 100）、穴 = (1,1)-(9,9) の 4 頂点（面積 64・CW）、
/// 正味面積 36。
///
/// 殺す変異: 穴を返さない（rings.len()==1・正味 100）・穴を独立した外環として返す（out.len()==2）・
/// 穴の向きを CCW のまま返す・穴を誤った外環へ割り当てる・外環に内部のバー境界を残す。
#[test]
fn union_four_bars_form_annulus_with_exactly_one_hole() {
    let bottom = rect(0.0, 0.0, 10.0, 1.0);
    let top = rect(0.0, 9.0, 10.0, 10.0);
    let left = rect(0.0, 0.0, 1.0, 10.0);
    let right = rect(9.0, 0.0, 10.0, 10.0);
    let inputs = vec![bottom, top, left, right];
    let out = union_rings(&inputs);

    assert_eq!(out.len(), 1, "アニュラスは 1 多角形");
    assert_output_structure(&out);
    assert_no_fabricated_vertices(&inputs, &out);
    assert_eq!(out[0].rings.len(), 2, "外環 1 + 穴 1");

    // 外環: (0,0)-(10,10)・CCW・面積 100。
    assert!(
        close(signed_area(&out[0].rings[0]), 100.0, AREA_TOL),
        "外環面積 {}",
        signed_area(&out[0].rings[0])
    );
    assert_ring_matches(
        &out[0].rings[0],
        &[(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)],
    );

    // 穴: (1,1)-(9,9)・CW（符号付き面積 −64）。CW なので頂点は逆回り。
    assert!(
        close(signed_area(&out[0].rings[1]), -64.0, AREA_TOL),
        "穴の符号付き面積 {}",
        signed_area(&out[0].rings[1])
    );
    assert_ring_matches(
        &out[0].rings[1],
        &[(1.0, 1.0), (1.0, 9.0), (9.0, 9.0), (9.0, 1.0)],
    );

    // 正味面積 = 100 − 64 = 36。
    assert!(
        close(net_area(&out[0]), 36.0, AREA_TOL),
        "{}",
        net_area(&out[0])
    );
}

/// アニュラスの穴の中に独立した島（小四角形）がある場合、
/// 穴付き多角形（面積降順で先頭）と島（2 番目）の 2 成分になる。
/// 島は穴の内部にあり、外環に囲まれているが**和領域としては別成分**。
///
/// 殺す変異: 島を穴の rings に混ぜる・島を落とす・穴を島の外環に割り当てる・
/// 面積降順が崩れて島が先頭に来る。
#[test]
fn union_island_inside_hole_is_returned_as_separate_component() {
    let bottom = rect(0.0, 0.0, 10.0, 1.0);
    let top = rect(0.0, 9.0, 10.0, 10.0);
    let left = rect(0.0, 0.0, 1.0, 10.0);
    let right = rect(9.0, 0.0, 10.0, 10.0);
    let island = rect(4.0, 4.0, 6.0, 6.0);
    let inputs = vec![bottom, top, left, right, island];
    let out = union_rings(&inputs);

    assert_eq!(out.len(), 2, "アニュラス + 島 = 2 成分");
    assert_output_structure(&out);
    assert_no_fabricated_vertices(&inputs, &out);

    // 先頭 = 外環面積 100 のアニュラス（穴 1 つ・正味 36）。
    assert_eq!(out[0].rings.len(), 2, "先頭は穴 1 つを持つ");
    assert!(
        close(net_area(&out[0]), 36.0, AREA_TOL),
        "{}",
        net_area(&out[0])
    );
    assert!(
        close(signed_area(&out[0].rings[1]), -64.0, AREA_TOL),
        "穴の符号付き面積 {}",
        signed_area(&out[0].rings[1])
    );

    // 2 番目 = 島（面積 4・穴なし）。
    assert_eq!(out[1].rings.len(), 1, "島に穴は無い");
    assert!(
        close(net_area(&out[1]), 4.0, AREA_TOL),
        "{}",
        net_area(&out[1])
    );
    assert_ring_matches(
        &out[1].rings[0],
        &[(4.0, 4.0), (6.0, 4.0), (6.0, 6.0), (4.0, 6.0)],
    );
}

// ============================================================
// 出力の一般契約（構造・捏造しない）
// ============================================================

/// 複雑な入力（重なり・共線・自己交差・離散を同時に含む）でも、出力の構造契約が保たれ、
/// 全頂点が入力辺の上に乗る（点を捏造しない）。
///
/// 殺す変異: 交点計算で座標をずらす（入力辺から外れた頂点が出る）・出力リングを閉じて返す・
/// 外環を CW で返す・面積降順を破る・面積ゼロの環を返す。
#[test]
fn union_complex_input_satisfies_structure_and_no_fabricated_vertices() {
    let inputs = vec![
        rect(0.0, 0.0, 4.0, 4.0),
        rect(2.0, 2.0, 6.0, 6.0),
        rect(4.0, 0.0, 8.0, 2.0),
        ring(&[(20.0, 0.0), (22.0, 0.0), (20.0, 2.0), (22.0, 2.0)]), // 八の字
        rect(30.0, 30.0, 30.5, 30.5),
    ];
    let out = union_rings(&inputs);

    assert!(!out.is_empty(), "有効な領域があるので非空");
    assert_output_structure(&out);
    assert_no_fabricated_vertices(&inputs, &out);

    // 面積降順が厳密に成り立つ（隣接ペアの比較を明示）。
    let areas: Vec<f64> = out
        .iter()
        .map(|poly| signed_area(&poly.rings[0]).abs())
        .collect();
    for w in areas.windows(2) {
        assert!(w[0] >= w[1] - AREA_TOL, "面積降順が崩れている（{areas:?}）");
    }
    // 最小の成分（面積 0.25 の四角）が消えていない＝小領域を丸め落としていない。
    assert!(
        areas.iter().any(|&a| close(a, 0.25, AREA_TOL)),
        "小さな独立成分（面積 0.25）が失われている（{areas:?}）"
    );
}

// ============================================================
// 回帰（レビュー指摘）
// ============================================================

/// **回帰**: 交点が量子化グリッド（1e-9 度）の境界付近に落ちても、環が閉じて出力が消えない。
///
/// `a` は底辺を `slope` 度だけ傾けた単位正方形 `[(0,0),(1,slope),(1,1),(0,1)]`、
/// `b` は下から突き刺さる細長い長方形 `lon∈[0.5,0.6] × lat∈[-1,0.5]`。
/// 2 領域は重なるので**単一多角形**が返らねばならない。
///
/// 期待面積は傾き 0 の実測ではなく**幾何から独立に**算出する:
/// - `area(a) = 1 − slope/2`（shoelace）
/// - `area(b) = 0.1 × 1.5 = 0.15`
/// - 重なり = `∫_{0.5}^{0.6} (0.5 − slope·x) dx = 0.05 − 0.055·slope`
/// - `union = 1.1 − 0.445·slope`
///
/// 殺す変異: 交点を辺ごとに独立に計算して別々に量子化する（同一交点が 1 ulp ずれて 2 点になり、
/// 環が繋がらない）・環が閉じなかったときに黙って捨てる（空 Vec / 成分欠落）・
/// 重なりを二重計上する（面積 1.15 付近）。
#[test]
fn union_regression_intersection_snapping_keeps_ring_closed() {
    for &slope in &[0.0, 5e-10, 1e-9, 1e-8, 1e-7] {
        let a = vec![p(0.0, 0.0), p(1.0, slope), p(1.0, 1.0), p(0.0, 1.0)];
        let b = vec![p(0.5, -1.0), p(0.6, -1.0), p(0.6, 0.5), p(0.5, 0.5)];
        let inputs = vec![a, b];
        let out = union_rings(&inputs);

        assert!(
            !out.is_empty(),
            "slope={slope}: 重なる 2 領域の和が空 Vec になった（環が閉じていない）"
        );
        assert_eq!(
            out.len(),
            1,
            "slope={slope}: 重なるので単一成分（実際 {} 件）",
            out.len()
        );
        assert_output_structure(&out);
        assert_no_fabricated_vertices(&inputs, &out);
        assert_eq!(out[0].rings.len(), 1, "slope={slope}: 穴は無い");

        let expected = 1.1 - 0.445 * slope;
        let got = net_area(&out[0]);
        assert!(
            close(got, expected, 1e-6),
            "slope={slope}: 面積 {got} が期待 {expected} と一致しない"
        );
        // 傾き 0 の面積 1.1 からも 1e-6 以内（傾きの微小変化で結果が飛ばない）。
        assert!(
            close(got, 1.1, 1e-6),
            "slope={slope}: 面積 {got} が傾き 0 の 1.1 から離れている"
        );
    }
}

/// **回帰**: 交差しないまま極めて近接して並走する 2 領域は、**2 つの多角形**として残る。
///
/// `a = (0,0)-(10,1)`・`b = (0, 1+gap)-(10,2)`。隙間 `gap` は 5e-8〜1e-5 度。
/// 対向する境界辺（`lat=1` と `lat=1+gap`）が「内部に埋もれた辺」と誤判定されて捨てられると、
/// 出力が空になったり 1 成分に融合したりする。
///
/// 殺す変異: 境界判定の法線オフセット幅を入力スケールと無関係な固定値にする（隙間より広い
/// オフセットで対向辺を内部と誤判定 → 空 Vec / 融合）・近接成分を 1 つに畳む・
/// 片方の成分を落とす。
#[test]
fn union_regression_near_parallel_rings_stay_separate() {
    for &gap in &[5e-8, 1e-7, 1e-6, 1e-5] {
        let a = rect(0.0, 0.0, 10.0, 1.0);
        let b = rect(0.0, 1.0 + gap, 10.0, 2.0);
        let inputs = vec![a, b];
        let out = union_rings(&inputs);

        assert!(
            !out.is_empty(),
            "gap={gap}: 近接する 2 領域の和が空 Vec になった（境界辺が捨てられている）"
        );
        assert_eq!(
            out.len(),
            2,
            "gap={gap}: 重なっていないので 2 成分（実際 {} 件）",
            out.len()
        );
        assert_output_structure(&out);
        assert_no_fabricated_vertices(&inputs, &out);
        for poly in &out {
            assert_eq!(poly.rings.len(), 1, "gap={gap}: 穴は無い");
        }

        // 面積降順なので先頭＝下の帯（10.0）、2 番目＝上の帯（10·(1−gap)）。
        let expected_upper = 10.0 * (1.0 - gap);
        assert!(
            close(net_area(&out[0]), 10.0, 1e-9),
            "gap={gap}: 下の帯の面積 {}",
            net_area(&out[0])
        );
        assert!(
            close(net_area(&out[1]), expected_upper, 1e-9),
            "gap={gap}: 上の帯の面積 {} が期待 {expected_upper} と一致しない",
            net_area(&out[1])
        );
    }
}

// ============================================================
// 回帰（ミューテーション）
// ============================================================

/// **A**: 頂点ちょうど 3（三角形）のリング 1 枚は、そのまま 1 多角形として返る。
/// `[(0,0),(4,0),(0,3)]`・面積 6・CCW・頂点 3。
///
/// 殺す変異: 入力正規化や point-in-polygon の「頂点 3 未満は無効」判定を `<= 3` / `== 3` に
/// ずらす（三角形が丸ごと捨てられて空 Vec になる）・三角形に余計な頂点を足す。
#[test]
fn union_regression_single_triangle_is_returned_unchanged() {
    let tri = ring(&[(0.0, 0.0), (4.0, 0.0), (0.0, 3.0)]);
    let inputs = vec![tri];
    let out = union_rings(&inputs);

    assert_eq!(out.len(), 1, "三角形 1 枚は 1 多角形");
    assert_output_structure(&out);
    assert_no_fabricated_vertices(&inputs, &out);
    assert_eq!(out[0].rings.len(), 1, "穴は無い");
    assert_eq!(out[0].rings[0].len(), 3, "頂点はちょうど 3");
    assert!(
        close(net_area(&out[0]), 6.0, AREA_TOL),
        "{}",
        net_area(&out[0])
    );
    assert_ring_matches(&out[0].rings[0], &[(0.0, 0.0), (4.0, 0.0), (0.0, 3.0)]);
}

/// **A**: 三角形 2 枚が重なる和。
/// `T1=[(0,0),(4,0),(0,4)]`（面積 8）と `T2=[(1,1),(5,1),(1,5)]`（面積 8）。
/// 重なりは `{x>=1, y>=1, x+y<=4}` の三角形 `(1,1),(3,1),(1,3)`（面積 2）なので
/// **和の面積は 8+8-2 = 14**（幾何から独立に算出）。外環は 7 頂点で、
/// 交点 `(3,1)`（T1 斜辺 ∩ T2 底辺）と `(1,3)`（T1 斜辺 ∩ T2 左辺）を含む。
///
/// 殺す変異: 三角形入力を捨てる（空 Vec）・重なりを二重計上（面積 16）・交点を打たない・
/// 斜辺と水平/垂直辺の交点計算の取り違え。
#[test]
fn union_regression_two_overlapping_triangles_exact_outline() {
    let t1 = ring(&[(0.0, 0.0), (4.0, 0.0), (0.0, 4.0)]);
    let t2 = ring(&[(1.0, 1.0), (5.0, 1.0), (1.0, 5.0)]);
    let inputs = vec![t1, t2];
    let out = union_rings(&inputs);

    assert_eq!(out.len(), 1, "重なるので単一成分");
    assert_output_structure(&out);
    assert_no_fabricated_vertices(&inputs, &out);
    assert_eq!(out[0].rings.len(), 1, "穴は無い");
    assert!(
        close(net_area(&out[0]), 14.0, AREA_TOL),
        "{}",
        net_area(&out[0])
    );
    assert_ring_matches(
        &out[0].rings[0],
        &[
            (0.0, 0.0),
            (4.0, 0.0),
            (3.0, 1.0),
            (5.0, 1.0),
            (1.0, 5.0),
            (1.0, 3.0),
            (0.0, 4.0),
        ],
    );
}

/// **A**: 三角形と四角形が重なる和。
/// `T=[(0,0),(4,0),(0,4)]`（面積 8）と `R=(1,1)-(5,5)`（面積 16）。
/// 重なりは三角形 `(1,1),(3,1),(1,3)`（面積 2）なので **和は 8+16-2 = 22**。外環は 8 頂点。
///
/// 殺す変異: 三角形側だけ/四角形側だけを採る・重なりの二重計上（面積 24）・
/// 斜辺との交点 `(3,1)`/`(1,3)` の欠落。
#[test]
fn union_regression_triangle_and_rectangle_overlap_exact_outline() {
    let tri = ring(&[(0.0, 0.0), (4.0, 0.0), (0.0, 4.0)]);
    let r = rect(1.0, 1.0, 5.0, 5.0);
    let inputs = vec![tri, r];
    let out = union_rings(&inputs);

    assert_eq!(out.len(), 1, "重なるので単一成分");
    assert_output_structure(&out);
    assert_no_fabricated_vertices(&inputs, &out);
    assert_eq!(out[0].rings.len(), 1, "穴は無い");
    assert!(
        close(net_area(&out[0]), 22.0, AREA_TOL),
        "{}",
        net_area(&out[0])
    );
    assert_ring_matches(
        &out[0].rings[0],
        &[
            (0.0, 0.0),
            (4.0, 0.0),
            (3.0, 1.0),
            (5.0, 1.0),
            (5.0, 5.0),
            (1.0, 5.0),
            (1.0, 3.0),
            (0.0, 4.0),
        ],
    );
}

/// **A**: 3 本の帯で直角三角形のアニュラスを作ると、**穴が三角形（頂点 3・CW）**になる。
/// 3 帯はいずれも大三角形 `x>=0, y>=0, x+y<=12` の内側に収まる台形:
/// 底帯 `[(0,0),(12,0),(10,2),(0,2)]`・左帯 `[(0,0),(2,0),(2,10),(0,12)]`・
/// 斜辺帯 `[(10,0),(12,0),(0,12),(0,10)]`。
/// 外環は大三角形 `(0,0),(12,0),(0,12)`（面積 72）、穴は `(2,2),(8,2),(2,8)`（面積 18・
/// 斜辺の内側は `x+y=10`）、正味 54。
///
/// 殺す変異: 頂点 3 の穴を「頂点 3 未満」と誤判定して捨てる（正味 72）・穴を外環として返す
/// （out.len()==2）・穴の向きを CCW のまま返す・穴の頂点を 4 個に増やす・
/// 斜辺帯（台形）を無視して穴が開かない。
#[test]
fn union_regression_three_strips_form_triangular_hole() {
    let bottom = ring(&[(0.0, 0.0), (12.0, 0.0), (10.0, 2.0), (0.0, 2.0)]);
    let left = ring(&[(0.0, 0.0), (2.0, 0.0), (2.0, 10.0), (0.0, 12.0)]);
    let hypotenuse = ring(&[(10.0, 0.0), (12.0, 0.0), (0.0, 12.0), (0.0, 10.0)]);
    let inputs = vec![bottom, left, hypotenuse];
    let out = union_rings(&inputs);

    assert_eq!(out.len(), 1, "アニュラスは 1 多角形");
    assert_output_structure(&out);
    assert_no_fabricated_vertices(&inputs, &out);
    assert_eq!(out[0].rings.len(), 2, "外環 1 + 穴 1");

    // 外環: 大三角形（面積 72）。斜辺帯の端点による共線頂点が残りうるので頂点数は縛らず、
    // 三角形の境界上にあることと 3 隅の存在で縛る。
    assert!(
        close(signed_area(&out[0].rings[0]), 72.0, AREA_TOL),
        "外環面積 {}",
        signed_area(&out[0].rings[0])
    );
    for pt in &out[0].rings[0] {
        let (lon, lat) = lonlat(pt);
        let on_bottom = close(lat, 0.0, TOL) && (-TOL..=12.0 + TOL).contains(&lon);
        let on_left = close(lon, 0.0, TOL) && (-TOL..=12.0 + TOL).contains(&lat);
        let on_hyp = close(lon + lat, 12.0, 1e-8);
        assert!(
            on_bottom || on_left || on_hyp,
            "外環頂点 ({lon},{lat}) が大三角形の境界上にない"
        );
    }
    for corner in [(0.0, 0.0), (12.0, 0.0), (0.0, 12.0)] {
        assert!(
            out[0].rings[0]
                .iter()
                .any(|pt| close(lonlat(pt).0, corner.0, TOL) && close(lonlat(pt).1, corner.1, TOL)),
            "隅 {corner:?} が外環頂点に無い"
        );
    }

    // 穴: 三角形 (2,2),(8,2),(2,8)・CW（符号付き面積 -18）・頂点ちょうど 3。
    assert_eq!(out[0].rings[1].len(), 3, "穴は三角形（頂点 3）");
    assert!(
        close(signed_area(&out[0].rings[1]), -18.0, AREA_TOL),
        "穴の符号付き面積 {}",
        signed_area(&out[0].rings[1])
    );
    assert_ring_matches(&out[0].rings[1], &[(2.0, 2.0), (2.0, 8.0), (8.0, 2.0)]);

    // 正味 72 - 18 = 54。
    assert!(
        close(net_area(&out[0]), 54.0, AREA_TOL),
        "{}",
        net_area(&out[0])
    );
}

/// **B**: 共線の**部分**重なりが、第三のリングによって**重なり区間の内部で**分断されるケース。
/// `A=(0,0)-(3,1)` の上辺と `B=(1,1)-(4,2)` の下辺は `lon in [1,3]` で共線重複する。
/// さらに `C=(1.5,0.5)-(2.5,1.5)` が `lat=1` を跨いで A・B 内部に収まり、その垂直辺が
/// 共線区間の**内部** `lon=1.5, 2.5`（A・B どちらの頂点でもない）で分割点を作る。
///
/// A と B は面積 3 ずつで重なりゼロ（接するのみ）、C は完全に A∪B の内部なので、
/// **和の面積は 3+3 = 6**（幾何から独立に算出）。外環は Z 字の 8 頂点で、
/// 分割点 `(1.5,1)` `(2.5,1)` は内部なので**出力に現れない**。
///
/// 殺す変異: 共線重なりの分割点登録を無効化する（内部の線が境界として残り、頂点 (1.5,1)/(2.5,1)
/// が外環に現れる／成分が割れる）・共線辺を畳まず 2 本残す（面積の二重計上・自己接触リング）・
/// C を独立成分として返す（out.len()>1）。
#[test]
fn union_regression_collinear_partial_overlap_split_in_its_interior() {
    let a = rect(0.0, 0.0, 3.0, 1.0);
    let b = rect(1.0, 1.0, 4.0, 2.0);
    let c = rect(1.5, 0.5, 2.5, 1.5);
    let inputs = vec![a, b, c];
    let out = union_rings(&inputs);

    assert_eq!(
        out.len(),
        1,
        "A∪B は lat=1 で繋がり、C は内部なので単一成分"
    );
    assert_output_structure(&out);
    assert_no_fabricated_vertices(&inputs, &out);
    assert_eq!(out[0].rings.len(), 1, "穴は無い");
    assert!(
        close(net_area(&out[0]), 6.0, AREA_TOL),
        "{}",
        net_area(&out[0])
    );
    assert_ring_matches(
        &out[0].rings[0],
        &[
            (0.0, 0.0),
            (3.0, 0.0),
            (3.0, 1.0),
            (4.0, 1.0),
            (4.0, 2.0),
            (1.0, 2.0),
            (1.0, 1.0),
            (0.0, 1.0),
        ],
    );
    // 共線区間の内部分割点は境界に現れない（内部に埋もれた辺として捨てられる）。
    for pt in &out[0].rings[0] {
        let (lon, lat) = lonlat(pt);
        assert!(
            !(close(lat, 1.0, TOL) && (close(lon, 1.5, TOL) || close(lon, 2.5, TOL))),
            "共線区間内部の分割点 ({lon},{lat}) が境界に残っている"
        );
    }
}

/// **B**: 完全に重なる共線辺（上下の長方形が `lat=1` を全長で共有）を、第三のリングが
/// その**内部**で分断するケース。`A=(0,0)-(4,1)`・`B=(0,1)-(4,2)`・
/// `C=(1,0.5)-(2,1.5)`（`lat=1` を跨ぎ、垂直辺が `lon=1,2` で共線辺を分断）。
/// 和は長方形 `(0,0)-(4,2)`＝**面積 8**。
///
/// 殺す変異: 分断された共線断片の一部を境界として残す（`lat=1` の内部線が出る・面積が 8 から
/// ずれる）・重複除去を分断後に行わず二重辺を残す・C を別成分にする。
#[test]
fn union_regression_duplicate_collinear_edge_split_in_its_interior() {
    let a = rect(0.0, 0.0, 4.0, 1.0);
    let b = rect(0.0, 1.0, 4.0, 2.0);
    let c = rect(1.0, 0.5, 2.0, 1.5);
    let inputs = vec![a, b, c];
    let out = union_rings(&inputs);

    assert_eq!(out.len(), 1, "上下が辺を共有するので単一成分");
    assert_output_structure(&out);
    assert_no_fabricated_vertices(&inputs, &out);
    assert_eq!(out[0].rings.len(), 1, "穴は無い");
    assert!(
        close(net_area(&out[0]), 8.0, AREA_TOL),
        "{}",
        net_area(&out[0])
    );

    // 全頂点が長方形 (0,0)-(4,2) の境界上（lat=1 の内部線が残っていない）。
    for pt in &out[0].rings[0] {
        let (lon, lat) = lonlat(pt);
        let on_v = close(lon, 0.0, TOL) || close(lon, 4.0, TOL);
        let on_h = close(lat, 0.0, TOL) || close(lat, 2.0, TOL);
        assert!(
            (on_v && (-TOL..=2.0 + TOL).contains(&lat))
                || (on_h && (-TOL..=4.0 + TOL).contains(&lon)),
            "頂点 ({lon},{lat}) が長方形 (0,0)-(4,2) の境界上にない（共線辺が残っている）"
        );
    }
    for corner in [(0.0, 0.0), (4.0, 0.0), (4.0, 2.0), (0.0, 2.0)] {
        assert!(
            out[0].rings[0]
                .iter()
                .any(|pt| close(lonlat(pt).0, corner.0, TOL) && close(lonlat(pt).1, corner.1, TOL)),
            "隅 {corner:?} が出力頂点に無い"
        );
    }
}

/// **C**: 1 点に 3 領域（＝境界辺 6 本）が集まり、各辺の角度が**すべて非対称**な配置。
/// 原点で接する 3 つの三角形:
/// `T1=[(0,0),(4,0),(4,1)]`（面積 2）・`T2=[(0,0),(3,3),(0,4)]`（面積 6）・
/// `T3=[(0,0),(-2,-1),(-1,-3)]`（面積 2.5）。互いに原点以外では交わらない。
/// 期待: **3 つの独立多角形**が面積降順 6 → 2.5 → 2 で返る。
///
/// 殺す変異: 環再結合の角度選択の符号を反転する（本来 3 つに分かれる成分が繋がって 1〜2 個に
/// なる・頂点列の順序が変わる）・点接触を融合する・成分を落とす・面積降順の崩れ。
#[test]
fn union_regression_three_regions_meeting_at_one_point_stay_separate() {
    let t1 = ring(&[(0.0, 0.0), (4.0, 0.0), (4.0, 1.0)]);
    let t2 = ring(&[(0.0, 0.0), (3.0, 3.0), (0.0, 4.0)]);
    let t3 = ring(&[(0.0, 0.0), (-2.0, -1.0), (-1.0, -3.0)]);
    let inputs = vec![t1, t2, t3];
    let out = union_rings(&inputs);

    assert_eq!(out.len(), 3, "1 点で出会う 3 領域は分離した 3 多角形");
    assert_output_structure(&out);
    assert_no_fabricated_vertices(&inputs, &out);
    for poly in &out {
        assert_eq!(poly.rings.len(), 1, "穴は無い");
        assert_eq!(poly.rings[0].len(), 3, "各成分は三角形（頂点 3）");
    }
    // 面積降順: T2(6) -> T3(2.5) -> T1(2)。頂点列も厳密に縛る。
    assert!(
        close(net_area(&out[0]), 6.0, AREA_TOL),
        "{}",
        net_area(&out[0])
    );
    assert_ring_matches(&out[0].rings[0], &[(0.0, 0.0), (3.0, 3.0), (0.0, 4.0)]);
    assert!(
        close(net_area(&out[1]), 2.5, AREA_TOL),
        "{}",
        net_area(&out[1])
    );
    assert_ring_matches(&out[1].rings[0], &[(0.0, 0.0), (-2.0, -1.0), (-1.0, -3.0)]);
    assert!(
        close(net_area(&out[2]), 2.0, AREA_TOL),
        "{}",
        net_area(&out[2])
    );
    assert_ring_matches(&out[2].rings[0], &[(0.0, 0.0), (4.0, 0.0), (4.0, 1.0)]);
}

/// **C**: 長方形の角と斜めの三角形が 1 点 `(3,2)` だけで接し、その点に集まる 4 本の辺の角度が
/// すべて異なる（180 度・270 度・約 14 度・約 63 度）非対称な分岐頂点。
/// `R=(0,0)-(3,2)`（面積 6）・`T=[(3,2),(7,3),(5,6)]`（面積 7）。
/// 期待: **2 つの多角形**が面積降順 7 → 6 で返る（内部は繋がっていない）。
///
/// 殺す変異: 角度選択の符号反転で 2 領域が 1 つの自己接触多角形に繋がる・接触点を交点として
/// 余計な頂点を増やす・順序が入力順のまま（面積降順の脱落）・片方を落とす。
#[test]
fn union_regression_asymmetric_point_contact_rectangle_and_triangle() {
    let r = rect(0.0, 0.0, 3.0, 2.0);
    let t = ring(&[(3.0, 2.0), (7.0, 3.0), (5.0, 6.0)]);
    let inputs = vec![r, t];
    let out = union_rings(&inputs);

    assert_eq!(out.len(), 2, "点接触なので 2 多角形");
    assert_output_structure(&out);
    assert_no_fabricated_vertices(&inputs, &out);

    // 先頭 = 三角形（面積 7）。
    assert_eq!(out[0].rings.len(), 1);
    assert_eq!(out[0].rings[0].len(), 3, "三角形は頂点 3");
    assert!(
        close(net_area(&out[0]), 7.0, AREA_TOL),
        "{}",
        net_area(&out[0])
    );
    assert_ring_matches(&out[0].rings[0], &[(3.0, 2.0), (7.0, 3.0), (5.0, 6.0)]);

    // 2 番目 = 長方形（面積 6）。
    assert_eq!(out[1].rings.len(), 1);
    assert_eq!(out[1].rings[0].len(), 4, "長方形は頂点 4");
    assert!(
        close(net_area(&out[1]), 6.0, AREA_TOL),
        "{}",
        net_area(&out[1])
    );
    assert_ring_matches(
        &out[1].rings[0],
        &[(0.0, 0.0), (3.0, 0.0), (3.0, 2.0), (0.0, 2.0)],
    );
}
