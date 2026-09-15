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
/// - `rings[1..]`（穴）は `rings[0]`（外環）の**内部**にある
///   （穴の全頂点が外環の閉内部・穴の面積 < 外環の面積）— ISSUE-053 §確定仕様 2
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
    assert_holes_inside_outer(poly);
}

/// **ISSUE-053 §確定仕様 2（構造契約）**: 穴は外環の内部にある。
/// - 穴の符号なし面積 < 外環の符号なし面積
/// - 穴の全頂点が外環の**閉内部**（even-odd 内部、または外環の辺の上）
fn assert_holes_inside_outer(poly: &GeoPolygon) {
    if poly.rings.len() < 2 {
        return;
    }
    let outer = ring_coords(&poly.rings[0]);
    let outer_area = signed_area(&poly.rings[0]).abs();
    for (index, hole) in poly.rings[1..].iter().enumerate() {
        let hole_area = signed_area(hole).abs();
        assert!(
            hole_area < outer_area - AREA_TOL,
            "穴{index}の面積 {hole_area} が外環の面積 {outer_area} 以上（穴が外環の内部に無い）\
             / 外環 {outer:?} / 穴 {:?}",
            ring_coords(hole)
        );
        for pt in hole {
            let q = lonlat(pt);
            let inside = point_in_ring_evenodd(&outer, q)
                || edges_of(&outer).any(|(a, b)| point_on_segment(q, a, b));
            assert!(
                inside,
                "穴{index}の頂点 {q:?} が外環の閉内部に無い（穴が誤った外環に割り当たっている）\
                 / 外環 {outer:?} / 穴 {:?}",
                ring_coords(hole)
            );
        }
    }
}

// ============================================================
// 点オラクル（ISSUE-053 §確定仕様 1 を縛るための独立実装）
//
// 期待値は**入力リングの even-odd 内外を点ごとに直接評価した値**から導出する。
// `union_rings` の出力は一切参照しない。
// ============================================================

/// `GeoPoint` 列を `(lon, lat)` 列へ。
fn ring_coords(r: &[GeoPoint]) -> Vec<(f64, f64)> {
    r.iter().map(lonlat).collect()
}

/// 閉リングとして解釈した辺の列（末尾→先頭を含む）。
fn edges_of(r: &[(f64, f64)]) -> impl Iterator<Item = ((f64, f64), (f64, f64))> + '_ {
    (0..r.len()).map(move |i| (r[i], r[(i + 1) % r.len()]))
}

/// 点 `q` がリング `r`（閉リングとして解釈）の **even-odd 内部**か。
/// 水平レイキャスト（+lon 方向）の交差本数の偶奇で判定する。境界上の点の結果は**不定**なので、
/// 呼び出し側で境界近傍の点を除外すること。
fn point_in_ring_evenodd(r: &[(f64, f64)], q: (f64, f64)) -> bool {
    if r.len() < 3 {
        return false;
    }
    let mut inside = false;
    for ((x1, y1), (x2, y2)) in edges_of(r) {
        // 半開区間 [min, max) で辺を数えることで、頂点を通るレイの二重計上を避ける。
        if (y1 > q.1) != (y2 > q.1) {
            let t = (q.1 - y1) / (y2 - y1);
            if x1 + t * (x2 - x1) > q.0 {
                inside = !inside;
            }
        }
    }
    inside
}

/// **点オラクル**: 「いずれかの入力リングの even-odd 内部」か。
fn oracle_inside(inputs: &[Vec<(f64, f64)>], q: (f64, f64)) -> bool {
    inputs
        .iter()
        .any(|r| r.len() >= 3 && point_in_ring_evenodd(r, q))
}

/// 出力の内外: 「外環の内部 ∧ どの穴の内部でもない」の OR。
fn output_inside(polys: &[GeoPolygon], q: (f64, f64)) -> bool {
    polys.iter().any(|poly| {
        let outer = ring_coords(&poly.rings[0]);
        point_in_ring_evenodd(&outer, q)
            && !poly.rings[1..]
                .iter()
                .any(|h| point_in_ring_evenodd(&ring_coords(h), q))
    })
}

/// 点 `q` から線分 `a-b` までの距離。
fn point_segment_distance(q: (f64, f64), a: (f64, f64), b: (f64, f64)) -> f64 {
    let (vx, vy) = (b.0 - a.0, b.1 - a.1);
    let (wx, wy) = (q.0 - a.0, q.1 - a.1);
    let len2 = vx * vx + vy * vy;
    if len2 <= 0.0 {
        return wx.hypot(wy);
    }
    let t = ((wx * vx + wy * vy) / len2).clamp(0.0, 1.0);
    (q.0 - (a.0 + t * vx)).hypot(q.1 - (a.1 + t * vy))
}

/// 境界近傍の除外幅（even-odd レイキャストは境界上で不定なので、この距離未満の点は評価しない）。
const ORACLE_CLEARANCE: f64 = 1e-6;

/// 格子上の点で、出力の内外が**点オラクル**と一致することを検証する
/// （ISSUE-053 §確定仕様 1）。入力辺から `ORACLE_CLEARANCE` 未満の点は判定不定なので飛ばす。
///
/// `label` は失敗時に配置を再現するための文脈（シード・リング座標など）。
/// 戻り値は実際に評価した点の数（オラクルが全点を除外していないことの確認用）。
fn check_point_oracle(
    inputs: &[Vec<GeoPoint>],
    polys: &[GeoPolygon],
    lo: f64,
    hi: f64,
    step: f64,
    label: &str,
) -> usize {
    let coords: Vec<Vec<(f64, f64)>> = inputs.iter().map(|r| ring_coords(r)).collect();
    let input_edges: Vec<((f64, f64), (f64, f64))> = coords
        .iter()
        .filter(|r| r.len() >= 3)
        .flat_map(|r| edges_of(r).collect::<Vec<_>>())
        .collect();
    // 格子は整数座標・半整数座標を避けるようオフセットする（入力頂点は整数格子上）。
    let mut evaluated = 0;
    let mut y = lo + 0.2417;
    while y < hi {
        let mut x = lo + 0.1731;
        while x < hi {
            let q = (x, y);
            let clear = input_edges
                .iter()
                .all(|&(a, b)| point_segment_distance(q, a, b) >= ORACLE_CLEARANCE);
            if clear {
                evaluated += 1;
                let want = oracle_inside(&coords, q);
                let got = output_inside(polys, q);
                assert_eq!(
                    got, want,
                    "{label}: 点 {q:?} の内外が点オラクルと違う（want={want} / got={got}）\
                     ＝入力の even-odd 和と一致していない"
                );
            }
            x += step;
        }
        y += step;
    }
    evaluated
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

// ============================================================
// 回帰（§11.7 手順 7・8 と退化方針「共線重なり」を直接観測する）
// ============================================================

/// 穴を**符号付き面積の絶対値**で探す（穴の並び順は契約外なので順序に依存しない）。
fn find_hole_by_area(poly: &GeoPolygon, area: f64) -> &[GeoPoint] {
    poly.rings[1..]
        .iter()
        .find(|h| close(signed_area(h).abs(), area, AREA_TOL))
        .unwrap_or_else(|| {
            let got: Vec<f64> = poly.rings[1..].iter().map(|h| signed_area(h)).collect();
            panic!("面積 {area} の穴が無い（穴の符号付き面積 {got:?}）")
        })
}

/// 出力リングに同一頂点が 2 度現れない（自己接触環・重複辺の残留を検出する）。
fn assert_no_repeated_vertex(r: &[GeoPoint]) {
    let coords: Vec<(f64, f64)> = r.iter().map(lonlat).collect();
    for i in 0..coords.len() {
        for j in (i + 1)..coords.len() {
            assert!(
                !(close(coords[i].0, coords[j].0, TOL) && close(coords[i].1, coords[j].1, TOL)),
                "同一頂点 {:?} がリングに重複して出現",
                coords[i]
            );
        }
    }
}

/// **§11.7 手順 8「各穴を、それを含む最小面積の外環に割り当てる」**（二重入れ子）。
///
/// 4 本のバーで大アニュラス（外環 `(0,0)-(20,16)`＝320・穴 `(1,1)-(19,15)`＝252）を作り、
/// その穴の中に 4 本の細いバーで**島アニュラス**（外環 `(5,4)-(13,11)`＝56・穴 `(6,5)-(12,10)`＝30）
/// を置く。島の穴は幾何的に**2 つの外環**（大外環・島外環）に含まれるので、
/// 「含む外環のうち最小のもの」＝島外環に割り当てられねばならない。
///
/// 期待: 2 多角形（面積降順 320 → 56）で、**それぞれちょうど 1 つの穴**を持つ。
/// 正味面積は 320−252=68 と 56−30=26。
///
/// 殺す変異: 穴を「最初に見つかった含む外環」や「最大の含む外環」へ割り当てる（先頭多角形が
/// 穴 2 つ・島が穴なし＝正味 38 と 56）・穴の包含判定を外環ではなく穴で行う・島の穴を
/// 独立多角形として返す（out.len()==3）・面積降順の崩れ。
#[test]
fn union_regression_hole_assigned_to_smallest_containing_outer_ring() {
    let big = [
        rect(0.0, 0.0, 20.0, 1.0),
        rect(0.0, 15.0, 20.0, 16.0),
        rect(0.0, 0.0, 1.0, 16.0),
        rect(19.0, 0.0, 20.0, 16.0),
    ];
    let island = [
        rect(5.0, 4.0, 13.0, 5.0),
        rect(5.0, 10.0, 13.0, 11.0),
        rect(5.0, 4.0, 6.0, 11.0),
        rect(12.0, 4.0, 13.0, 11.0),
    ];
    let inputs: Vec<Vec<GeoPoint>> = big.iter().chain(island.iter()).cloned().collect();
    let out = union_rings(&inputs);

    assert_eq!(out.len(), 2, "大アニュラス + 島アニュラス = 2 多角形");
    assert_output_structure(&out);
    assert_no_fabricated_vertices(&inputs, &out);

    // 先頭 = 大アニュラス。穴はちょうど 1 つ（島の穴を奪ってはならない）。
    assert_eq!(out[0].rings.len(), 2, "大アニュラスの穴はちょうど 1 つ");
    assert_ring_matches(
        &out[0].rings[0],
        &[(0.0, 0.0), (20.0, 0.0), (20.0, 16.0), (0.0, 16.0)],
    );
    assert!(
        close(signed_area(&out[0].rings[1]), -252.0, AREA_TOL),
        "大アニュラスの穴の符号付き面積 {}",
        signed_area(&out[0].rings[1])
    );
    assert_ring_matches(
        &out[0].rings[1],
        &[(1.0, 1.0), (1.0, 15.0), (19.0, 15.0), (19.0, 1.0)],
    );
    assert!(
        close(net_area(&out[0]), 68.0, AREA_TOL),
        "{}",
        net_area(&out[0])
    );

    // 2 番目 = 島アニュラス。自分の穴をちょうど 1 つ持つ。
    assert_eq!(out[1].rings.len(), 2, "島アニュラスの穴はちょうど 1 つ");
    assert_ring_matches(
        &out[1].rings[0],
        &[(5.0, 4.0), (13.0, 4.0), (13.0, 11.0), (5.0, 11.0)],
    );
    assert!(
        close(signed_area(&out[1].rings[1]), -30.0, AREA_TOL),
        "島の穴の符号付き面積 {}",
        signed_area(&out[1].rings[1])
    );
    assert_ring_matches(
        &out[1].rings[1],
        &[(6.0, 5.0), (6.0, 10.0), (12.0, 10.0), (12.0, 5.0)],
    );
    assert!(
        close(net_area(&out[1]), 26.0, AREA_TOL),
        "{}",
        net_area(&out[1])
    );
}

/// **§11.7 手順 7「各頂点の出辺を入射辺の逆向きからの角度で選ぶ（最も右回り側）」**:
/// 穴が**内部の 1 点でくびれて 2 つに分かれる**配置。
///
/// 下バー `(0,0)-(10,1)`・上バー `(0,9)-(10,10)`・左バー `(0,0)-(1,10)` と、右側を塞ぐ
/// くさび三角形 `[(10,0),(10,10),(1,4)]`。くさびの頂点 `(1,4)` は左バーの右辺 `lon=1` の
/// 途中に乗る。頂点 `(1,4)` には境界半辺が 4 本（`lon=1` の上下 2 本・くさびの 2 辺、
/// 傾き −4/9 と 2/3 で非対称）集まる。
///
/// 期待: 1 多角形、外環は正方形 `(0,0)-(10,10)`（面積 100）、穴は**2 つ**:
/// 下穴 `(1,1),(7.75,1),(1,4)`（面積 10.125）・上穴 `(1,4),(8.5,9),(1,9)`（面積 18.75）。
/// 正味 100 − 28.875 = 71.125。
///
/// 殺す変異: 分岐頂点で最も右回りでなく最も左回り（または入射順・任意）の出辺を選ぶ
/// （2 穴が `(1,4)` を 2 度通る 1 本の自己接触環になり rings.len()==2・穴の頂点 6）・
/// 角度比較の基準を入射辺の逆向きにしない・くさびの端点接触 `(1,4)` を分割点にしない
/// （穴が 1 つになり面積が変わる）。
#[test]
fn union_regression_branch_vertex_pinches_hole_into_two_holes() {
    let bottom = rect(0.0, 0.0, 10.0, 1.0);
    let top = rect(0.0, 9.0, 10.0, 10.0);
    let left = rect(0.0, 0.0, 1.0, 10.0);
    let wedge = ring(&[(10.0, 0.0), (10.0, 10.0), (1.0, 4.0)]);
    let inputs = vec![bottom, top, left, wedge];
    let out = union_rings(&inputs);

    assert_eq!(out.len(), 1, "くびれた穴は同じ外環の 2 穴（多角形は 1 つ）");
    assert_output_structure(&out);
    assert_no_fabricated_vertices(&inputs, &out);
    assert_eq!(
        out[0].rings.len(),
        3,
        "外環 1 + 穴 2（自己接触した 1 穴に融合してはならない）"
    );
    assert_ring_matches(
        &out[0].rings[0],
        &[(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)],
    );

    let lower = find_hole_by_area(&out[0], 10.125);
    assert_eq!(lower.len(), 3, "下穴は三角形");
    assert_ring_matches(lower, &[(1.0, 1.0), (1.0, 4.0), (7.75, 1.0)]);
    let upper = find_hole_by_area(&out[0], 18.75);
    assert_eq!(upper.len(), 3, "上穴は三角形");
    assert_ring_matches(upper, &[(1.0, 4.0), (1.0, 9.0), (8.5, 9.0)]);

    assert!(
        close(net_area(&out[0]), 71.125, AREA_TOL),
        "{}",
        net_area(&out[0])
    );
}

/// **§11.7 手順 7（最も右回りの出辺）**: 穴の角が**外環に 1 点で接する**配置。
///
/// 下バー `(0,0)-(10,1)`・上バー `(0,9)-(10,10)`・右バー `(9,0)-(10,10)` と、左側を塞ぐ
/// 2 つの三角形 `[(0,0),(2.5,0),(0,5)]`・`[(0,5),(0,10),(5,10)]`。2 三角形は `(0,5)` だけを
/// 共有し、そこで外環（`lon=0` の上下）と穴の 2 辺（傾き 5/−2.5 と 5/5、非対称）の
/// 計 4 本の境界半辺が集まる。
///
/// 期待: 1 多角形、外環は正方形 `(0,0)-(10,10)`（面積 100・全頂点が正方形の境界上）、
/// 穴はちょうど 1 つ `(0,5),(4,9),(9,9),(9,1),(2,1)`（CW・面積 60）。正味 40。
///
/// 殺す変異: 分岐頂点 `(0,5)` で外環の歩行が穴の辺へ曲がる（外環と穴が 1 本の自己接触環に
/// 融合して rings.len()==1・正味面積が 40 と一致しない）・穴側の歩行が外環へ抜ける・
/// 接触点 `(0,5)` を穴の頂点から落とす（穴が 4 頂点・面積 60 以外になる）。
#[test]
fn union_regression_branch_vertex_hole_touching_outer_ring_stays_hole() {
    let bottom = rect(0.0, 0.0, 10.0, 1.0);
    let top = rect(0.0, 9.0, 10.0, 10.0);
    let right = rect(9.0, 0.0, 10.0, 10.0);
    let lower_wedge = ring(&[(0.0, 0.0), (2.5, 0.0), (0.0, 5.0)]);
    let upper_wedge = ring(&[(0.0, 5.0), (0.0, 10.0), (5.0, 10.0)]);
    let inputs = vec![bottom, top, right, lower_wedge, upper_wedge];
    let out = union_rings(&inputs);

    assert_eq!(out.len(), 1, "穴が外環に接しても多角形は 1 つ");
    assert_output_structure(&out);
    assert_no_fabricated_vertices(&inputs, &out);
    assert_eq!(
        out[0].rings.len(),
        2,
        "外環 1 + 穴 1（融合して 1 環になってはならない）"
    );

    // 外環: 正方形 (0,0)-(10,10)。接触点 (0,5) は lon=0 上の共線点なので頂点数は縛らず、
    // 面積 100・全頂点が正方形の境界上・4 隅の存在で縛る。
    assert!(
        close(signed_area(&out[0].rings[0]), 100.0, AREA_TOL),
        "外環面積 {}",
        signed_area(&out[0].rings[0])
    );
    for pt in &out[0].rings[0] {
        let (lon, lat) = lonlat(pt);
        let on_v = close(lon, 0.0, TOL) || close(lon, 10.0, TOL);
        let on_h = close(lat, 0.0, TOL) || close(lat, 10.0, TOL);
        assert!(
            (on_v && (-TOL..=10.0 + TOL).contains(&lat))
                || (on_h && (-TOL..=10.0 + TOL).contains(&lon)),
            "外環頂点 ({lon},{lat}) が正方形 (0,0)-(10,10) の境界上にない（穴と融合している疑い）"
        );
    }
    for corner in [(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)] {
        assert!(
            out[0].rings[0]
                .iter()
                .any(|pt| close(lonlat(pt).0, corner.0, TOL) && close(lonlat(pt).1, corner.1, TOL)),
            "隅 {corner:?} が外環頂点に無い"
        );
    }

    // 穴: 5 頂点・CW・面積 60。接触点 (0,5) を頂点として含む。
    assert_eq!(out[0].rings[1].len(), 5, "穴は 5 頂点");
    assert!(
        close(signed_area(&out[0].rings[1]), -60.0, AREA_TOL),
        "穴の符号付き面積 {}",
        signed_area(&out[0].rings[1])
    );
    assert_ring_matches(
        &out[0].rings[1],
        &[(0.0, 5.0), (4.0, 9.0), (9.0, 9.0), (9.0, 1.0), (2.0, 1.0)],
    );
    assert!(
        close(net_area(&out[0]), 40.0, AREA_TOL),
        "{}",
        net_area(&out[0])
    );
}

/// **§11.7 手順 3「共線重なりも分割点として登録する」＋手順 5「無向重複除去」**:
/// 共線重なりの分割点が、共線登録**以外**からは得られない配置（外周側・水平）。
///
/// `A=(0,0)-(10,2)` と `B=[(3,0),(6,0),(12,0),(12,5),(3,5)]`。両者は下辺 `lat=0` を
/// `lon∈[3,10]` で共線重複し、この区間は**両方から見て外側**（＝和の境界に残る）。
/// B の頂点 `(6,0)` は A の下辺の内部に乗るが、その両隣の辺 `(3,0)-(6,0)`・`(6,0)-(12,0)` は
/// ともに A の辺と共線で、`(6,0)` を通る横断辺は存在しない。したがって A の下辺の `lon=6` での
/// 分割は共線重なりの登録からしか生じない。分割されないと A の断片 `(3,0)-(10,0)` と
/// B の断片 `(3,0)-(6,0)`・`(6,0)-(10,0)` が別キーとなって畳まれず、境界辺が二重に残る。
///
/// 期待: 1 多角形・穴なし・外環は `(0,0),(12,0),(12,5),(3,5),(3,2),(0,2)`（共線点は落ちる）・
/// 面積 20+45−14 = **51**・同一頂点の重複なし。
///
/// 殺す変異: 共線重なりの分割点登録を無効化する／端点だけ登録し内部の頂点を落とす
/// （境界辺の二重残留＝頂点の重複出現・余計な頂点・面積不一致）・重複除去を分割前に行う。
#[test]
fn union_regression_collinear_overlap_is_sole_split_source_horizontal() {
    let a = rect(0.0, 0.0, 10.0, 2.0);
    let b = ring(&[(3.0, 0.0), (6.0, 0.0), (12.0, 0.0), (12.0, 5.0), (3.0, 5.0)]);
    let inputs = vec![a, b];
    let out = union_rings(&inputs);

    assert_eq!(out.len(), 1, "重なるので単一成分");
    assert_output_structure(&out);
    assert_no_fabricated_vertices(&inputs, &out);
    assert_eq!(out[0].rings.len(), 1, "穴は無い");
    assert_no_repeated_vertex(&out[0].rings[0]);
    assert!(
        close(net_area(&out[0]), 51.0, AREA_TOL),
        "{}",
        net_area(&out[0])
    );
    assert_ring_matches(
        &out[0].rings[0],
        &[
            (0.0, 0.0),
            (12.0, 0.0),
            (12.0, 5.0),
            (3.0, 5.0),
            (3.0, 2.0),
            (0.0, 2.0),
        ],
    );
}

/// **§11.7 手順 3・5（共線重なりの登録と重複除去）**: 同上を**垂直辺・共線頂点が反対側のリング**
/// にある配置で縛る（座標軸の取り違え・「短い辺側だけ登録」の欠陥を検出）。
///
/// `A=[(0,0),(2,0),(2,10),(0,10),(0,7)]`（左辺 `lon=0` に共線頂点 `(0,7)`）と
/// `B=(0,3)-(5,14)`。両者は左辺 `lon=0` を `lat∈[3,10]` で共線重複し、この区間は外周。
/// `(0,7)` の両隣の辺は B の左辺と共線で、横断辺は無い。
///
/// 期待: 1 多角形・穴なし・外環 `(0,0),(2,0),(2,3),(5,3),(5,14),(0,14)`・
/// 面積 20+55−14 = **61**・同一頂点の重複なし。
///
/// 殺す変異: 共線重なりの登録を水平辺（`lat` 一致）にしか適用しない・重なり端点だけを
/// 登録し相手辺の内部頂点を落とす・共線判定に `lon` と `lat` を取り違える。
#[test]
fn union_regression_collinear_overlap_is_sole_split_source_vertical() {
    let a = ring(&[(0.0, 0.0), (2.0, 0.0), (2.0, 10.0), (0.0, 10.0), (0.0, 7.0)]);
    let b = rect(0.0, 3.0, 5.0, 14.0);
    let inputs = vec![a, b];
    let out = union_rings(&inputs);

    assert_eq!(out.len(), 1, "重なるので単一成分");
    assert_output_structure(&out);
    assert_no_fabricated_vertices(&inputs, &out);
    assert_eq!(out[0].rings.len(), 1, "穴は無い");
    assert_no_repeated_vertex(&out[0].rings[0]);
    assert!(
        close(net_area(&out[0]), 61.0, AREA_TOL),
        "{}",
        net_area(&out[0])
    );
    assert_ring_matches(
        &out[0].rings[0],
        &[
            (0.0, 0.0),
            (2.0, 0.0),
            (2.0, 3.0),
            (5.0, 3.0),
            (5.0, 14.0),
            (0.0, 14.0),
        ],
    );
}

/// **§11.7 手順 3・5（共線重なりの登録と重複除去）**: **両リングとも**共線区間の内部に
/// 共線頂点を持ち、互いに相手の辺の内部で分割点を要求する配置。
///
/// `A=[(0,0),(8,0),(10,0),(10,2),(0,2)]`（下辺に共線頂点 `(8,0)`）と
/// `B=[(3,0),(6,0),(12,0),(12,5),(3,5)]`（下辺に共線頂点 `(6,0)`）。共線重複区間は
/// `lon∈[3,10]`。`(6,0)` は A の辺 `(0,0)-(8,0)` の内部、`(8,0)` は B の辺 `(6,0)-(12,0)` の
/// 内部にあり、いずれも横断辺を持たない。両方向の登録が無いと畳めない断片が残る。
///
/// 期待: 1 多角形・穴なし・外環 `(0,0),(12,0),(12,5),(3,5),(3,2),(0,2)`・面積 **51**・
/// 同一頂点の重複なし。
///
/// 殺す変異: 辺ペアの片側（第 1 辺）にだけ分割点を登録する・重なり区間の端点のみ登録する・
/// 共線重なりの `t` を一方の辺のパラメータで両辺に流用する。
#[test]
fn union_regression_collinear_overlap_split_required_on_both_rings() {
    let a = ring(&[(0.0, 0.0), (8.0, 0.0), (10.0, 0.0), (10.0, 2.0), (0.0, 2.0)]);
    let b = ring(&[(3.0, 0.0), (6.0, 0.0), (12.0, 0.0), (12.0, 5.0), (3.0, 5.0)]);
    let inputs = vec![a, b];
    let out = union_rings(&inputs);

    assert_eq!(out.len(), 1, "重なるので単一成分");
    assert_output_structure(&out);
    assert_no_fabricated_vertices(&inputs, &out);
    assert_eq!(out[0].rings.len(), 1, "穴は無い");
    assert_no_repeated_vertex(&out[0].rings[0]);
    assert!(
        close(net_area(&out[0]), 51.0, AREA_TOL),
        "{}",
        net_area(&out[0])
    );
    assert_ring_matches(
        &out[0].rings[0],
        &[
            (0.0, 0.0),
            (12.0, 0.0),
            (12.0, 5.0),
            (3.0, 5.0),
            (3.0, 2.0),
            (0.0, 2.0),
        ],
    );
}

/// **§11.7 手順 7（最も右回りの出辺）＋退化方針「点接触」**: 2 領域が**2 つの異なる頂点だけ**で接する配置
/// （軸平行）。
///
/// `A=(0,2)-(10,4)` と `B=[(0,0),(10,0),(8,2),(5,1),(2,2)]`。B の頂点 `(2,2)`・`(8,2)` は A の下辺の
/// 内部に乗るが、B の上縁 `(2,2)-(5,1)-(8,2)` は A の下辺から離れるので、A と B は**辺を共有せず**
/// 2 点でのみ接する。内部は繋がっていないので §11.7「点接触は分離した多角形」により **2 多角形**
/// （A と B がそのまま）で、両者に囲まれた三角形 `(2,2),(8,2),(5,1)` は**どちらの穴にもならない**。
/// 両接触点では境界半辺が 4 本ずつ集まる。
///
/// 殺す変異: 分岐頂点で右回りでなく左回りの出辺を選ぶ（(2,2)/(8,2) で A のチェーンと B のチェーンが
/// 対にされ、三角形を囲む 1 環＝面積 36 の外環＋穴、または向きの逆転した環になる。2 点接触は同一頂点を
/// 2 度通らないので頂点重複の分割では修復できない）・接触点の落とし。
#[test]
fn union_regression_two_regions_touching_at_two_vertices_stay_separate_axis_aligned() {
    let a = rect(0.0, 2.0, 10.0, 4.0);
    let b = ring(&[(0.0, 0.0), (10.0, 0.0), (8.0, 2.0), (5.0, 1.0), (2.0, 2.0)]);
    let inputs = vec![a, b];
    let out = union_rings(&inputs);

    assert_eq!(out.len(), 2, "2 点接触のみ＝内部が繋がらないので 2 多角形");
    assert_output_structure(&out);
    assert_no_fabricated_vertices(&inputs, &out);
    assert_eq!(out[0].rings.len(), 1, "A に穴は無い");
    assert_eq!(out[1].rings.len(), 1, "B に穴は無い");
    assert_ring_matches(
        &out[0].rings[0],
        &[(0.0, 2.0), (10.0, 2.0), (10.0, 4.0), (0.0, 4.0)],
    );
    assert!(
        close(signed_area(&out[0].rings[0]), 20.0, AREA_TOL),
        "{}",
        signed_area(&out[0].rings[0])
    );
    assert_ring_matches(
        &out[1].rings[0],
        &[(0.0, 0.0), (10.0, 0.0), (8.0, 2.0), (5.0, 1.0), (2.0, 2.0)],
    );
    assert!(
        close(signed_area(&out[1].rings[0]), 13.0, AREA_TOL),
        "{}",
        signed_area(&out[1].rings[0])
    );
}

/// **§11.7 手順 7（最も右回りの出辺）＋退化方針「点接触」**: 同上を**非軸平行**な配置で縛る
/// （角度の鏡映ミスも検出）。
///
/// `A=[(0,2),(10,3),(10,6),(1,5)]`（下辺は `lat = 2 + 0.1·lon` の斜線）と
/// `B=[(0,0),(10,0),(8,2.8),(5,1),(2,2.2)]`。B の頂点 `(2,2.2)`・`(8,2.8)` は A の斜めの下辺の
/// 内部に乗り、A と B は 2 点でのみ接する。各接触点に集まる 4 半辺の角度はすべて異なる。
/// 期待: **2 多角形**（A: 面積 28・B: 面積 15.5）・穴なし。
///
/// 殺す変異: 出辺選択の角度符号の反転／鏡映（軸平行では偶然通る）・入射辺の逆向きを基準にしない・
/// 接触点の落とし・A と B の融合。
#[test]
fn union_regression_two_regions_touching_at_two_vertices_stay_separate_skewed() {
    let a = ring(&[(0.0, 2.0), (10.0, 3.0), (10.0, 6.0), (1.0, 5.0)]);
    let b = ring(&[(0.0, 0.0), (10.0, 0.0), (8.0, 2.8), (5.0, 1.0), (2.0, 2.2)]);
    let inputs = vec![a, b];
    let out = union_rings(&inputs);

    assert_eq!(out.len(), 2, "2 点接触のみ＝内部が繋がらないので 2 多角形");
    assert_output_structure(&out);
    assert_no_fabricated_vertices(&inputs, &out);
    assert_eq!(out[0].rings.len(), 1, "A に穴は無い");
    assert_eq!(out[1].rings.len(), 1, "B に穴は無い");
    assert_ring_matches(
        &out[0].rings[0],
        &[(0.0, 2.0), (10.0, 3.0), (10.0, 6.0), (1.0, 5.0)],
    );
    assert!(
        close(signed_area(&out[0].rings[0]), 28.0, AREA_TOL),
        "{}",
        signed_area(&out[0].rings[0])
    );
    assert_ring_matches(
        &out[1].rings[0],
        &[(0.0, 0.0), (10.0, 0.0), (8.0, 2.8), (5.0, 1.0), (2.0, 2.2)],
    );
    assert!(
        close(signed_area(&out[1].rings[0]), 15.5, AREA_TOL),
        "{}",
        signed_area(&out[1].rings[0])
    );
}

/// 共線重なり＋点接触の共通検証: 主多角形（面積 51・6 頂点）と点接触する三角形 C（面積 3）が
/// **別々の 2 多角形**として、どちらも失われずに返る。
fn assert_collinear_overlap_plus_touching_triangle(
    inputs: &[Vec<GeoPoint>],
    c_expected: &[(f64, f64)],
    label: &str,
) {
    let out = union_rings(inputs);
    assert_eq!(
        out.len(),
        2,
        "{label}: 点接触の C は別多角形（主 + C = 2、実際 {} 件）",
        out.len()
    );
    assert_output_structure(&out);
    assert_no_fabricated_vertices(inputs, &out);

    // 先頭 = 主多角形（面積 51・穴なし・頂点の重複なし）。
    assert_eq!(out[0].rings.len(), 1, "{label}: 主多角形に穴は無い");
    assert_no_repeated_vertex(&out[0].rings[0]);
    assert!(
        close(net_area(&out[0]), 51.0, AREA_TOL),
        "{label}: 主多角形の面積 {}",
        net_area(&out[0])
    );
    assert_ring_matches(
        &out[0].rings[0],
        &[
            (0.0, 0.0),
            (12.0, 0.0),
            (12.0, 5.0),
            (3.0, 5.0),
            (3.0, 2.0),
            (0.0, 2.0),
        ],
    );

    // 2 番目 = C（面積 3・三角形そのもの）。
    assert_eq!(out[1].rings.len(), 1, "{label}: C に穴は無い");
    assert_eq!(out[1].rings[0].len(), 3, "{label}: C は三角形");
    assert!(
        close(net_area(&out[1]), 3.0, AREA_TOL),
        "{label}: C の面積 {}",
        net_area(&out[1])
    );
    assert_ring_matches(&out[1].rings[0], c_expected);
}

/// **§11.7 手順 3・5（共線重なりの登録と重複除去）＋点接触方針**: 共線重なりの**遠端**
/// `(10,0)`（A の頂点・B の辺 `(6,0)-(12,0)` の内部）に、第三の三角形 C が 1 点で接する。
///
/// `A=(0,0)-(10,2)`・`B=[(3,0),(6,0),(12,0),(12,5),(3,5)]`・`C=[(10,0),(9,-3),(11,-3)]`。
/// 共線重複区間 `lon∈[3,10]` が正しく分割・畳み込みされないと、二重に残った断片が
/// `(10,0)` で C の境界へ継ぎ足され、開いた歩行の切り捨てでは消えない誤環（C と主多角形の
/// 融合・C の欠落・面積のずれ）が生じる。入力順は **C が先頭／末尾** の両方を検証する
/// （歩行順で切り捨て挙動が変わりうる）。
///
/// 期待: 2 多角形（主 51 → C 3）。主は 6 頂点、C は三角形 `(9,-3),(11,-3),(10,0)` そのもの。
///
/// 殺す変異: 共線重なりの分割点登録を無効化する・重複除去を分割前に行う・点接触で C を主へ
/// 融合する・歩行の打ち切りで C を落とす。
#[test]
fn union_regression_collinear_overlap_with_triangle_touching_far_endpoint() {
    let a = rect(0.0, 0.0, 10.0, 2.0);
    let b = ring(&[(3.0, 0.0), (6.0, 0.0), (12.0, 0.0), (12.0, 5.0), (3.0, 5.0)]);
    let c = ring(&[(10.0, 0.0), (9.0, -3.0), (11.0, -3.0)]);
    let c_expected = [(9.0, -3.0), (11.0, -3.0), (10.0, 0.0)];

    assert_collinear_overlap_plus_touching_triangle(
        &[c.clone(), a.clone(), b.clone()],
        &c_expected,
        "C 先頭",
    );
    assert_collinear_overlap_plus_touching_triangle(&[a, b, c], &c_expected, "C 末尾");
}

/// **§11.7 手順 3・5（共線重なりの登録と重複除去）＋点接触方針**: 共線重なりの**直進頂点**
/// `(6,0)`（B の共線頂点・A の下辺の内部）に、第三の三角形 C が 1 点で接する。
///
/// `A=(0,0)-(10,2)`・`B=[(3,0),(6,0),(12,0),(12,5),(3,5)]`・`C=[(6,0),(5,-3),(7,-3)]`。
/// `(6,0)` は共線登録が無ければ A の辺上に分割点として存在しない頂点で、C の接触により
/// そこに境界半辺が集まる。A の断片が `(6,0)` で割れていないと、`(6,0)` に集まる半辺の
/// 集合が不整合になり、C の環が主多角形へ継ぎ足されるか失われる。
/// 入力順は **C が先頭／末尾** の両方を検証する。
///
/// 期待: 2 多角形（主 51 → C 3）。主は 6 頂点、C は三角形 `(5,-3),(7,-3),(6,0)` そのもの。
///
/// 殺す変異: 共線重なりの分割点登録の無効化（相手辺の内部頂点での分割欠落）・点接触の融合・
/// 開いた歩行の切り捨てで C を落とす。
#[test]
fn union_regression_collinear_overlap_with_triangle_touching_straight_through_vertex() {
    let a = rect(0.0, 0.0, 10.0, 2.0);
    let b = ring(&[(3.0, 0.0), (6.0, 0.0), (12.0, 0.0), (12.0, 5.0), (3.0, 5.0)]);
    let c = ring(&[(6.0, 0.0), (5.0, -3.0), (7.0, -3.0)]);
    let c_expected = [(5.0, -3.0), (7.0, -3.0), (6.0, 0.0)];

    assert_collinear_overlap_plus_touching_triangle(
        &[c.clone(), a.clone(), b.clone()],
        &c_expected,
        "C 先頭",
    );
    assert_collinear_overlap_plus_touching_triangle(&[a, b, c], &c_expected, "C 末尾");
}

/// **§11.7 手順 3・5（共線重なりの分割点登録）— 縦向き**: 直進頂点に接する三角形の配置
/// （`union_regression_collinear_overlap_with_triangle_touching_straight_through_vertex`）を
/// lon/lat を入れ替え（さらに lon を +5 平行移動して負の経度を避ける）**縦の共線重なり**にする。重なりの分割点登録がパラメータ計算を
/// 横方向でしか正しく行わない（`lon` 成分だけで t を求める・縦線分で 0 除算になる）変異は、
/// 横向きの配置では通るが縦向きでは分割が失われて主多角形が壊れる。
///
/// 期待は横向き配置の転置（外環は転置で向きが反転するので CCW へ戻す）: 主多角形の外環
/// `(5,0),(7,0),(7,3),(10,3),(10,12),(5,12)`（面積 51・穴なし）、C は `(5,6),(2,7),(2,5)`（直進頂点に接する場合）または `(5,10),(2,11),(2,9)`（遠端に接する場合）（面積 3）。
#[test]
fn union_regression_collinear_overlap_with_triangle_touching_straight_through_vertex_vertical() {
    // 転置＋平行移動: (lon,lat) -> (lat+5, lon)。入力リングは転置で CW になるので反転して CCW に戻す。
    let t = |pts: &[(f64, f64)]| -> Vec<GeoPoint> {
        let mut v: Vec<(f64, f64)> = pts.iter().map(|&(x, y)| (y + 5.0, x)).collect();
        v.reverse();
        ring(&v)
    };
    let a = t(&[(0.0, 0.0), (10.0, 0.0), (10.0, 2.0), (0.0, 2.0)]);
    let b = t(&[(3.0, 0.0), (6.0, 0.0), (12.0, 0.0), (12.0, 5.0), (3.0, 5.0)]);
    // C は直進頂点 (6,0) に接する三角形と、重なりの遠端 (10,0) に接する三角形の 2 種
    // （横向き配置の 2 テストに対応）。遠端側では C の横断辺が (6,0) の分割を供給しないので、
    // 縦線分の共線重なりの分割点登録が唯一の供給源になる。
    let c_mid = t(&[(6.0, 0.0), (5.0, -3.0), (7.0, -3.0)]);
    let c_far = t(&[(10.0, 0.0), (9.0, -3.0), (11.0, -3.0)]);
    let c_mid_expected = [(5.0, 6.0), (2.0, 7.0), (2.0, 5.0)];
    let c_far_expected = [(5.0, 10.0), (2.0, 11.0), (2.0, 9.0)];

    for (inputs, c_expected, label) in [
        (
            vec![c_mid.clone(), a.clone(), b.clone()],
            &c_mid_expected,
            "C(直進頂点) 先頭",
        ),
        (
            vec![a.clone(), b.clone(), c_mid],
            &c_mid_expected,
            "C(直進頂点) 末尾",
        ),
        (
            vec![c_far.clone(), a.clone(), b.clone()],
            &c_far_expected,
            "C(遠端) 先頭",
        ),
        (vec![a, b, c_far], &c_far_expected, "C(遠端) 末尾"),
    ] {
        let out = union_rings(&inputs);
        assert_eq!(out.len(), 2, "{label}: 主 + 点接触の C = 2 多角形");
        assert_output_structure(&out);
        assert_no_fabricated_vertices(&inputs, &out);
        assert_eq!(out[0].rings.len(), 1, "{label}: 主多角形に穴は無い");
        assert_no_repeated_vertex(&out[0].rings[0]);
        assert!(
            close(net_area(&out[0]), 51.0, AREA_TOL),
            "{label}: 主多角形の面積 {}",
            net_area(&out[0])
        );
        assert_ring_matches(
            &out[0].rings[0],
            &[
                (5.0, 0.0),
                (7.0, 0.0),
                (7.0, 3.0),
                (10.0, 3.0),
                (10.0, 12.0),
                (5.0, 12.0),
            ],
        );
        assert_eq!(out[1].rings.len(), 1, "{label}: C に穴は無い");
        assert!(
            close(net_area(&out[1]), 3.0, AREA_TOL),
            "{label}: C の面積 {}",
            net_area(&out[1])
        );
        assert_ring_matches(&out[1].rings[0], c_expected);
    }
}

// ============================================================
// ISSUE-052: 経度フレーム非依存（反子午線跨ぎ）
// ============================================================
//
// ISSUE-052 §確定仕様:
//   1. 正規化は `union_rings` の内部（公開シグネチャは不変）。
//   2. `|Δlon| > 180` の辺が 1 本も無ければ回転量 0（従来出力を変えない）。
//   3. 回転量は「全入力頂点の経度の**最大の空き区間**の中央に ±180 の継ぎ目が来る」ように決める。
//   4. 回転後も `|Δlon| > 180` の辺が残れば（＝全経度を覆う）**回転せず**従来どおり解く。
//   5. 結果は逆回転で戻す。緯度は不変、経度は `wrap` の丸めで最下位桁が動きうる。
//
// ここでのテストは「球面上で同一の領域を、複数の経度フレームで与えても同じ結果が返る」ことを
// 縛る。比較は**参照フレーム**（その領域が継ぎ目を跨がない経度表現）へ出力を回してから行う。
// 跨ぐフレームでの出力は経度が ±180 を挟むので、平面 shoelace（`signed_area`）を直接当てても
// 意味を持たない（それ自体は実装の誤りではない）ため、必ず参照フレームで測る。
//
// 許容について: 回転は `wrap(lon + offset)` の加減算なので、|lon| < 180 の範囲で最下位桁
// （~3e-14 度）が動きうる。座標比較は既存の `TOL`（1e-9）、面積比較は `AREA_TOL`（1e-7・
// 面積 200〜3400 に対して相対 1e-9 以下）で、いずれも丸めの伝播より数桁厳しい。

/// 経度を [-180, 180) へ畳む（`EastLongitude::from_degrees` と同じ約束）。
fn wrap180(deg: f64) -> f64 {
    let mut v = (deg + 180.0) % 360.0;
    if v < 0.0 {
        v += 360.0;
    }
    v - 180.0
}

/// 「真の経度」で書いた点列を、`offset` だけ回した経度表現のリングにする。
/// `offset = 0` が ISSUE-052 の実測配置（反子午線跨ぎ）に当たる。
fn framed(pts: &[(f64, f64)], offset: f64) -> Vec<GeoPoint> {
    pts.iter()
        .map(|&(lon, lat)| p(wrap180(lon + offset), lat))
        .collect()
}

/// 真の経度で書いた軸平行長方形を、`offset` だけ回した表現のリングにする。
fn framed_rect(lon0: f64, lat0: f64, lon1: f64, lat1: f64, offset: f64) -> Vec<GeoPoint> {
    framed(
        &[(lon0, lat0), (lon1, lat0), (lon1, lat1), (lon0, lat1)],
        offset,
    )
}

/// 出力リングを `delta` だけ回して比較用フレームへ移す。
fn rotate_ring(r: &[GeoPoint], delta: f64) -> Vec<GeoPoint> {
    r.iter()
        .map(|pt| {
            let (lon, lat) = lonlat(pt);
            p(wrap180(lon + delta), lat)
        })
        .collect()
}

/// `delta` だけ回したフレームで測った正味面積（外環 − 穴）。
fn rotated_net_area(poly: &GeoPolygon, delta: f64) -> f64 {
    let outer = signed_area(&rotate_ring(&poly.rings[0], delta)).abs();
    let holes: f64 = poly.rings[1..]
        .iter()
        .map(|h| signed_area(&rotate_ring(h, delta)).abs())
        .sum();
    outer - holes
}

/// リングの外接ボックス `(lon_min, lon_max, lat_min, lat_max)`。
fn bbox(r: &[GeoPoint]) -> (f64, f64, f64, f64) {
    let mut b = (
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::INFINITY,
        f64::NEG_INFINITY,
    );
    for pt in r {
        let (lon, lat) = lonlat(pt);
        b.0 = b.0.min(lon);
        b.1 = b.1.max(lon);
        b.2 = b.2.min(lat);
        b.3 = b.3.max(lat);
    }
    b
}

/// **ISSUE-052 §確定仕様 1・3・5（headline）**: 球面上で同一の領域は、どの経度表現で与えても
/// 同じ成分数・同じ面積・同じ輪郭を返す。
///
/// 領域は ISSUE-052「背景（実測）」の配置そのもの: 真の経度 `170..190`・緯度 `0..10` を
/// `170..179` と `179..190` の **2 枚**に割って与える（ユニオンが継ぎ目を跨いで融合せねばならない）。
/// 正しい面積は**幾何から独立に**求まる: 経度幅 20 × 緯度幅 10 = **200**（ticket の (b) と一致）。
///
/// フレームは `offset = 0`（= 跨ぐ表現。`179 → -170` の辺が平面上 −349°）・`-180`（ticket の (b)・
/// 参照フレーム）・`+90`・`-90`・`+45`（いずれも跨がない）。比較は参照フレーム `-180` で行い、
/// 領域は `-10..10` に写る。
///
/// 期待 RED（修正前）: `offset = 0` のフレームだけが、面積 **3490**・外環の経度幅 349° を返す
/// （成分数は 1 のままなので、赤くなるのは面積と輪郭の assert）。他のフレームは緑。
///
/// 殺す変異: 経度を単なる平面座標として扱う（＝現状）・回転を掛けて戻し忘れる（輪郭が別フレームに
/// 出る）・回転量の符号を逆にする（別の辺が跨ぐ）・最大の空き区間でなく最小の空き区間に継ぎ目を置く。
#[test]
fn union_frame_invariant_area_for_seam_crossing_region() {
    const REF: f64 = -180.0;
    let west = [(170.0, 0.0), (179.0, 0.0), (179.0, 10.0), (170.0, 10.0)];
    let east = [(179.0, 0.0), (190.0, 0.0), (190.0, 10.0), (179.0, 10.0)];

    for &offset in &[0.0_f64, -180.0, 90.0, -90.0, 45.0] {
        let inputs = vec![framed(&west, offset), framed(&east, offset)];
        let out = union_rings(&inputs);

        assert_eq!(
            out.len(),
            1,
            "offset={offset}: 2 枚は辺を共有して繋がるので単一成分（実際 {} 件）",
            out.len()
        );
        assert_eq!(out[0].rings.len(), 1, "offset={offset}: 穴は無い");

        // 参照フレーム（-180）へ戻して測る。
        let delta = REF - offset;
        let outer = rotate_ring(&out[0].rings[0], delta);
        assert!(
            signed_area(&outer) > 0.0,
            "offset={offset}: 参照フレームの外環は CCW（符号付き面積 {}）",
            signed_area(&outer)
        );
        let area = rotated_net_area(&out[0], delta);
        assert!(
            close(area, 200.0, AREA_TOL),
            "offset={offset}: 面積 {area} が正しい 200 と一致しない（跨ぎ表現で 3490 になる欠陥）"
        );
        let (lon_min, lon_max, lat_min, lat_max) = bbox(&outer);
        assert!(
            close(lon_min, -10.0, TOL)
                && close(lon_max, 10.0, TOL)
                && close(lat_min, 0.0, TOL)
                && close(lat_max, 10.0, TOL),
            "offset={offset}: 参照フレームの外接ボックスが (-10..10, 0..10) でない \
             （{lon_min}..{lon_max}, {lat_min}..{lat_max}）"
        );
        assert_ring_matches(
            &outer,
            &[(-10.0, 0.0), (10.0, 0.0), (10.0, 10.0), (-10.0, 10.0)],
        );
    }
}

/// **ISSUE-052 §確定仕様 1・5**: 球面上でだけ重なる 2 枚が、継ぎ目を跨いで**融合**する。
///
/// 真の経度で `A = 160..185`・`B = 180..200`（緯度 0..10）。重なりは `180..185` で、
/// `offset = 0` の表現では A が `160 → -175`（跨ぐ）、B が `-180 → -160` になる。
/// 平面として読むと両者は重ならない（A が地球を逆走する）ので、欠陥実装は融合できない。
/// 正しい和は経度幅 40 × 緯度幅 10 = **400**、参照フレーム `-180` では `-20..20`。
///
/// さらに **確定仕様 5「結果は逆回転で戻す」** を直接縛る: `offset = 0` の生の出力は
/// 与えたフレームの値、すなわち真の経度 160 と 200 に対応する `160` と `-160` を頂点に持つ
/// （回転したフレームのまま返してはならない）。
///
/// 期待 RED（修正前）: `offset = 0` で成分数・面積・輪郭のいずれも一致しない。
///
/// 殺す変異: 跨ぎ検出を辺単位で行わない・回転後の結果を戻さない・
/// 重なり判定を回転前の座標で行う。
#[test]
fn union_seam_crossing_overlapping_rects_merge_into_one_component() {
    const REF: f64 = -180.0;

    for &offset in &[0.0_f64, -180.0, 90.0, -90.0, 45.0] {
        let a = framed_rect(160.0, 0.0, 185.0, 10.0, offset);
        let b = framed_rect(180.0, 0.0, 200.0, 10.0, offset);
        let out = union_rings(&[a, b]);

        assert_eq!(
            out.len(),
            1,
            "offset={offset}: 球面上では 180..185 で重なるので単一成分（実際 {} 件）",
            out.len()
        );
        assert_eq!(out[0].rings.len(), 1, "offset={offset}: 穴は無い");

        let delta = REF - offset;
        let outer = rotate_ring(&out[0].rings[0], delta);
        let area = rotated_net_area(&out[0], delta);
        assert!(
            close(area, 400.0, AREA_TOL),
            "offset={offset}: 面積 {area} が正しい 400 と一致しない"
        );
        assert_ring_matches(
            &outer,
            &[(-20.0, 0.0), (20.0, 0.0), (20.0, 10.0), (-20.0, 10.0)],
        );
    }

    // 確定仕様 5: 跨ぐフレームで与えたら、跨ぐフレームのまま返る（内部回転を戻す）。
    let a = framed_rect(160.0, 0.0, 185.0, 10.0, 0.0);
    let b = framed_rect(180.0, 0.0, 200.0, 10.0, 0.0);
    let out = union_rings(&[a, b]);
    assert_eq!(out.len(), 1);
    assert_ring_matches(
        &out[0].rings[0],
        &[(160.0, 0.0), (-160.0, 0.0), (-160.0, 10.0), (160.0, 10.0)],
    );
}

/// **ISSUE-052 §確定仕様 5「回転は経度の平行移動のみで緯度に触れない」**:
/// 回転が発生するフレームでも、出力頂点の緯度は入力の緯度と**ビット単位で**一致する。
///
/// 入力の緯度は 0 と 10 の 2 種類しかないので、出力の全頂点の緯度はそのいずれかに厳密一致
/// せねばならない（比較対象は同じ構築経路 `p()` を通した値なので、丸めの差は生じない）。
/// 経度側は `wrap` の丸めで最下位桁が動きうるため、ここでは縛らない（確定仕様 5 の明記どおり）。
///
/// 期待: 修正の前後どちらでも緑になりうる（現状は回転そのものが無いので自明に緑）。
/// 緯度にも `wrap`/正規化を掛けてしまう変異・経度と緯度を取り違えて回す変異を殺す。
#[test]
fn union_rotation_leaves_vertex_latitudes_bit_exact() {
    let west = [(170.0, 0.0), (179.0, 0.0), (179.0, 10.0), (170.0, 10.0)];
    let east = [(179.0, 0.0), (190.0, 0.0), (190.0, 10.0), (179.0, 10.0)];
    let allowed = [lonlat(&p(0.0, 0.0)).1, lonlat(&p(0.0, 10.0)).1];

    for &offset in &[0.0_f64, 45.0] {
        let inputs = vec![framed(&west, offset), framed(&east, offset)];
        let out = union_rings(&inputs);
        assert!(!out.is_empty(), "offset={offset}: 出力が空");
        for poly in &out {
            for r in &poly.rings {
                for pt in r {
                    let lat = lonlat(pt).1;
                    assert!(
                        allowed.contains(&lat),
                        "offset={offset}: 出力頂点の緯度 {lat} が入力の緯度 {allowed:?} と厳密一致しない（回転が緯度に触れている）"
                    );
                }
            }
        }
    }
}

/// **ISSUE-052 §確定仕様 2「跨ぐ辺が 1 本も無ければ何もしない」**: 跨ぎ判定の閾値が
/// `|Δlon| > 180`（**厳密な超過**）であることを、閾値ぎりぎりの 2 配置で縛る。
///
/// - `|Δlon| = 179`: `A = 0..179`・`B = 100..179`（B は A に内包）→ 和は A そのもの・面積 1790。
/// - `|Δlon| = 180`（ちょうど・跨がない）: `A = -90..90`・`B = 0..90` → 和は A・面積 1800。
///
/// どちらも回転してはならない配置で、既存の 40 本（すべて経度幅 30 度以下の非跨ぎ入力）が
/// 担保していない**閾値近傍**を埋める。バイト不変そのもの（出力の最下位桁が動かないこと）は
/// 公開 API からは観測できないため、ここでは「形が変わらない」ことまでを縛る。
///
/// 期待: 修正の前後どちらでも緑（回帰テスト）。
///
/// 殺す変異: 跨ぎ判定の閾値を 180 以外（90・170 等）にする・
/// 経度差に `abs()` を付け忘れて片側だけ検出する。
#[test]
fn union_non_crossing_input_near_180_threshold_is_unchanged() {
    // |Δlon| = 179
    let a = rect(0.0, 0.0, 179.0, 10.0);
    let b = rect(100.0, 0.0, 179.0, 10.0);
    let inputs = vec![a, b];
    let out = union_rings(&inputs);
    assert_eq!(out.len(), 1, "B は A に内包されるので単一成分");
    assert_output_structure(&out);
    assert_no_fabricated_vertices(&inputs, &out);
    assert_eq!(out[0].rings.len(), 1, "穴は無い");
    assert!(
        close(net_area(&out[0]), 1790.0, AREA_TOL),
        "{}",
        net_area(&out[0])
    );
    assert_ring_matches(
        &out[0].rings[0],
        &[(0.0, 0.0), (179.0, 0.0), (179.0, 10.0), (0.0, 10.0)],
    );

    // |Δlon| = 180 ちょうど（「> 180」ではないので回転しない）
    let a = rect(-90.0, 0.0, 90.0, 10.0);
    let b = rect(0.0, 0.0, 90.0, 10.0);
    let inputs = vec![a, b];
    let out = union_rings(&inputs);
    assert_eq!(out.len(), 1, "B は A に内包されるので単一成分");
    assert_output_structure(&out);
    assert_no_fabricated_vertices(&inputs, &out);
    assert_eq!(out[0].rings.len(), 1, "穴は無い");
    assert!(
        close(net_area(&out[0]), 1800.0, AREA_TOL),
        "{}",
        net_area(&out[0])
    );
    assert_ring_matches(
        &out[0].rings[0],
        &[(-90.0, 0.0), (90.0, 0.0), (90.0, 10.0), (-90.0, 10.0)],
    );
}

/// **ISSUE-052 §確定仕様 4「検証して退避する」**: 全経度を覆う（＝極を囲む）入力では、
/// どこに継ぎ目を置いても `|Δlon| > 180` の辺が残るので、回転は**失敗**し、従来どおり解く。
///
/// 入力は経度 −180..150 を 30 度刻みで単調に一周するリング（緯度は 80/85 を交互に振って
/// 平面面積を持たせる）。末尾 `150` から先頭 `-180` へ戻る辺が平面上 −330° で、頂点の
/// 空き区間（どれも 30°）の中央に継ぎ目を置いても、その区間を跨ぐ辺が必ず残る。
///
/// **幾何的な正しさは意図的に検証しない**: ISSUE-052 §非目的が「極を囲む領域は本 issue の
/// 対象外・退避する（結果は従来と同じく未保証）」と明記しているため。ここで縛るのは
/// 「パニックせず・結果を返し・同じ入力に対して決定的である」ことだけで、退避の結果に
/// 期待値を置くと据え置きの仕様を勝手に確定させてしまう。
///
/// 殺す変異: 回転後の再検証を省いて壊れた回転結果を返す（＝別の答えを捏造する）・
/// 回転失敗時に panic / unwrap する・退避経路を非決定にする。
#[test]
fn union_all_longitude_coverage_falls_back_without_panicking() {
    let pts: Vec<(f64, f64)> = (0..12)
        .map(|k| {
            let lon = -180.0 + f64::from(k) * 30.0;
            let lat = if k % 2 == 0 { 80.0 } else { 85.0 };
            (lon, lat)
        })
        .collect();
    let ring_all_lons = framed(&pts, 0.0);

    // 退避しても panic しない。
    let first = union_rings(std::slice::from_ref(&ring_all_lons));
    // 決定的（同じ入力から同じ出力）。
    let second = union_rings(&[ring_all_lons]);
    assert_eq!(
        canonical(&first),
        canonical(&second),
        "全経度を覆う入力での退避結果が決定的でない"
    );
}

/// **ISSUE-052 §確定仕様 3「最大の空き区間の中央に継ぎ目を置く」**: 空き区間が**狭い**
/// （＝領域がほぼ全経度に広がる）配置でも、回転が自分で新しい継ぎ目を作らずフレーム不変を保つ。
///
/// 領域は真の経度 `0..340`・緯度 `0..10`。頂点は 10 度刻みに密に置く（`0,10,…,340` の 35 本）ので、
/// 頂点間の空き区間は領域内部ではどれも 10°、領域外の `340..360` だけが **20°** で唯一の最大となる。
/// 継ぎ目はその中央 350° に置かれ、回転量は −170° になる。`0..340` は参照フレーム `-170` で
/// `-170..170` に写り、どの辺も `|Δlon| = 10` で跨がない。
/// 入力は緯度で重なる 2 本の帯（`lat 0..6` と `lat 4..10`）にして、ユニオンが実際に融合を要求する。
/// 正しい面積は経度幅 340 × 緯度幅 10 = **3400**。
///
/// フレームは `0`・`90`・`-90`（いずれも跨ぐ）と `-165`・`-170`・`-175`（継ぎ目が真の空き区間
/// `340..360` に入るので跨がない）。
///
/// 期待 RED（修正前）: 跨ぐ 3 フレームで面積・外接ボックスが一致しない。
///
/// 殺す変異: 空き区間を「最大」でなく「最初」や「最小」で選ぶ（継ぎ目が領域内部に落ち、
/// 回転後も跨ぐ辺が残って退避 → 跨ぎフレームの結果が壊れたまま）・空き区間の中央でなく端に
/// 継ぎ目を置く（境界の頂点がちょうど ±180 に乗って跨ぎ判定が揺れる）・
/// 円周上の並べ替えで `340 → 0` の折り返し区間を数え落とす。
#[test]
fn union_frame_invariant_when_largest_longitude_gap_is_narrow() {
    const REF: f64 = -170.0;
    /// 経度 0..=340（10 度刻み・35 頂点）を往復する帯リング（真の経度）。
    fn band(lat0: f64, lat1: f64) -> Vec<(f64, f64)> {
        let mut pts: Vec<(f64, f64)> = (0..=34).map(|k| (f64::from(k) * 10.0, lat0)).collect();
        pts.extend((0..=34).rev().map(|k| (f64::from(k) * 10.0, lat1)));
        pts
    }

    for &offset in &[0.0_f64, 90.0, -90.0, -165.0, -170.0, -175.0] {
        let lower = framed(&band(0.0, 6.0), offset);
        let upper = framed(&band(4.0, 10.0), offset);
        let out = union_rings(&[lower, upper]);

        assert_eq!(
            out.len(),
            1,
            "offset={offset}: 緯度で重なる 2 帯は単一成分（実際 {} 件）",
            out.len()
        );
        assert_eq!(out[0].rings.len(), 1, "offset={offset}: 穴は無い");

        let delta = REF - offset;
        let outer = rotate_ring(&out[0].rings[0], delta);
        assert!(
            signed_area(&outer) > 0.0,
            "offset={offset}: 参照フレームの外環は CCW（符号付き面積 {}）",
            signed_area(&outer)
        );
        let area = rotated_net_area(&out[0], delta);
        assert!(
            close(area, 3400.0, AREA_TOL),
            "offset={offset}: 面積 {area} が正しい 3400 と一致しない"
        );
        // 共線頂点が残るか否かは仕様外なので頂点列は縛らず、外接ボックスで縛る。
        let (lon_min, lon_max, lat_min, lat_max) = bbox(&outer);
        assert!(
            close(lon_min, -170.0, TOL)
                && close(lon_max, 170.0, TOL)
                && close(lat_min, 0.0, TOL)
                && close(lat_max, 10.0, TOL),
            "offset={offset}: 参照フレームの外接ボックスが (-170..170, 0..10) でない （{lon_min}..{lon_max}, {lat_min}..{lat_max}）"
        );
    }
}

/// **ISSUE-052 §確定仕様 3・4（候補の再試行）**: 最大の空き区間は**頂点**から求めるので、
/// その区間を**無関係な領域の辺**が頂点なしで跨いでいることがある。そこに継ぎ目を置くと
/// その領域が跨ぐようになり検証（確定仕様 4）が失敗する。候補を 1 つしか試さない実装は
/// そこで回転を諦めて offset 0 に落ち、**本来直すべき跨ぎ領域を壊れたまま返す**
/// （＝本 issue が直したはずのバグの再現）。空き区間を広い順に試し、跨ぐ辺が残らない
/// 最初の候補を採る実装だけが通る。
///
/// 配置（レビュー指摘の入力）:
/// - `A`: 真の経度 `170..190`・緯度 `0..5`（＝ `(170,0),(-170,0),(-170,5),(170,5)`）。**跨ぐ**。
///   面積は幾何から独立に 20 × 5 = **100**。
/// - `B`: 経度 `-80..90`・緯度 `20..30`。跨がないが**幅 170°**。面積 170 × 10 = **1700**。
/// - 両者は緯度で離れているので **2 成分**（面積降順に B → A）。
///
/// 頂点経度は `-170, -80, 90, 170`。空き区間は順に `-80→90`（**170°**・最大）・`-170→-80`（90°）・
/// `90→170`（80°）・`170→-170`（20°）。最大の 170° の区間は **B 自身の辺**が走っているので、
/// 継ぎ目を中央（真の経度 5°）に置くと B が跨ぎ、検証が失敗する。次点の 90° の区間
/// （真の経度 190..280・中央 235°・回転量 −55°）なら A も B も跨がず、回転が成立する。
///
/// フレームは `0`（A が跨ぐ）・`90`・`-140`・`180`（いずれも B が跨ぐ）・`-70`（参照フレーム・
/// どの辺も跨がない）。比較は参照フレーム `-70` へ戻して行い、A は `100..120`、B は `-150..20` に写る。
///
/// 期待 RED: 候補を 1 つしか試さない実装では、`-70` 以外の全フレームで退避が起き、
/// 跨ぐ側の成分が平面上の巨大な環（A なら経度幅 340°）になって面積 100 の assert が落ちる。
///
/// 殺す変異: 候補の再試行を止めて最初の 1 つで諦める・候補を幅の降順でなく昇順/入力順で試す・
/// 検証（跨ぐ辺が残っていないか）を回転**前**の座標で行う・最初に成立した候補でなく最後の候補を採る。
#[test]
fn union_retries_gap_candidates_when_widest_gap_is_occupied_by_another_edge() {
    const REF: f64 = -70.0;
    let a = [(170.0, 0.0), (190.0, 0.0), (190.0, 5.0), (170.0, 5.0)];
    let b = [(-80.0, 20.0), (90.0, 20.0), (90.0, 30.0), (-80.0, 30.0)];

    for &offset in &[0.0_f64, 90.0, -140.0, 180.0, -70.0] {
        let inputs = vec![framed(&a, offset), framed(&b, offset)];
        let out = union_rings(&inputs);

        assert_eq!(
            out.len(),
            2,
            "offset={offset}: A と B は緯度で離れているので 2 成分（実際 {} 件）",
            out.len()
        );
        for poly in &out {
            assert_eq!(poly.rings.len(), 1, "offset={offset}: 穴は無い");
        }

        let delta = REF - offset;
        // 先頭 = B（面積 1700）。
        let b_outer = rotate_ring(&out[0].rings[0], delta);
        assert!(
            signed_area(&b_outer) > 0.0,
            "offset={offset}: 参照フレームで B の外環は CCW（符号付き面積 {}）",
            signed_area(&b_outer)
        );
        let b_area = rotated_net_area(&out[0], delta);
        assert!(
            close(b_area, 1700.0, AREA_TOL),
            "offset={offset}: B の面積 {b_area} が正しい 1700 と一致しない"
        );
        assert_ring_matches(
            &b_outer,
            &[(-150.0, 20.0), (20.0, 20.0), (20.0, 30.0), (-150.0, 30.0)],
        );

        // 2 番目 = A（面積 100）。跨ぎが直っていないと平面上の経度幅 340° の環になる。
        let a_outer = rotate_ring(&out[1].rings[0], delta);
        assert!(
            signed_area(&a_outer) > 0.0,
            "offset={offset}: 参照フレームで A の外環は CCW（符号付き面積 {}）",
            signed_area(&a_outer)
        );
        let a_area = rotated_net_area(&out[1], delta);
        assert!(
            close(a_area, 100.0, AREA_TOL),
            "offset={offset}: A の面積 {a_area} が正しい 100 と一致しない\
             （最大の空き区間が B の辺に塞がれ、候補 1 つで諦めて跨ぎを直し損ねた疑い）"
        );
        assert_ring_matches(
            &a_outer,
            &[(100.0, 0.0), (120.0, 0.0), (120.0, 5.0), (100.0, 5.0)],
        );
    }
}

/// **ISSUE-052 §確定仕様 4「検証して退避する」**: **すべての**空き区間候補が他の領域の辺に
/// 塞がれている（＝入力が全経度を覆う）場合は、再試行しても回転は成立せず、従来どおり解く。
/// 単一リングで一周する `union_all_longitude_coverage_falls_back_without_panicking` に対し、
/// こちらは**複数の跨がない/跨ぐ帯の重ね合わせ**で全経度が覆われる場合を縛る
/// （再試行ループが候補を使い切った先で無限ループ・panic しないこと）。
///
/// 配置: 緯度で離れた 3 本の帯。`R1` 真の経度 `-170..-10`・`R2` `-50..110`・`R3` `90..250`（跨ぐ）。
/// 合わせて経度 `-170..250`（420°）＝全経度を覆う。頂点経度は `-170,-110,-50,-10,90,110` で、
/// 6 つの空き区間はいずれかの帯の辺が走っているので、どの候補に継ぎ目を置いても跨ぐ辺が残る。
///
/// **幾何的な正しさは検証しない**: ISSUE-052 §非目的が全経度被覆を対象外（結果は未保証）と
/// 明記しているため。縛るのは「パニックせず・結果を返し・決定的である」ことだけ。
///
/// 殺す変異: 候補の再試行を終端しない（無限ループ・候補リストの添字外参照で panic）・
/// 候補を使い切ったときに退避せず最後の壊れた回転結果を返す・退避経路を非決定にする。
#[test]
fn union_all_gap_candidates_occupied_falls_back_without_panicking() {
    let r1 = framed_rect(-170.0, 0.0, -10.0, 10.0, 0.0);
    let r2 = framed_rect(-50.0, 20.0, 110.0, 30.0, 0.0);
    let r3 = framed_rect(90.0, 40.0, 250.0, 50.0, 0.0);
    let inputs = vec![r1, r2, r3];

    let first = union_rings(&inputs);
    let second = union_rings(&inputs);
    assert_eq!(
        canonical(&first),
        canonical(&second),
        "全候補が塞がれた入力での退避結果が決定的でない"
    );
}

// ============================================================
// 回帰（ISSUE-053: 退化配置での穴の割り当て）
// ============================================================

/// **ISSUE-053 反例 C**（4 リング・複数リングが同一頂点で会合する退化配置）。
///
/// 期待値は**点オラクル**（入力リングの even-odd 和）から導出し、実装の出力からは取らない。
/// 縛るのは §確定仕様 1（点ごとの一致）と §確定仕様 2（穴は外環の内部）。
///
/// 既知の症状: 微小な外環にそれより大きい穴が割り当たり、点 `(1.6731, 3.7417)` /
/// `(2.1731, 3.7417)` が「どの入力リングの even-odd 内部でもない」のに出力では内部になる。
///
/// 殺す変異: 穴を面積や包含でなく別の基準（最初に見つかった外環・入力順）で割り当てる・
/// 穴の帰属先を探す包含判定を落とす・自己接触頂点で環の後継選択を誤る。
#[test]
fn union_regression_issue053_counterexample_c_matches_point_oracle() {
    let inputs = vec![
        ring(&[(1.0, 4.0), (3.0, 4.0), (3.0, 5.0)]),
        ring(&[(3.0, 3.0), (5.0, 5.0), (3.0, 5.0)]),
        ring(&[(2.0, 2.0), (3.0, 3.0), (1.0, 4.0)]),
        ring(&[(3.0, 0.0), (3.0, 4.0), (1.0, 3.0), (2.0, 2.0), (2.0, 4.0)]),
    ];
    let out = union_rings(&inputs);

    // 構造契約（穴が外環の内部にあること）は共通ヘルパで縛る。
    assert_output_structure(&out);
    assert_no_fabricated_vertices(&inputs, &out);

    // 点ごとの一致。issue 記載の反例点 2 つは、この格子（オフセット 0.1731/0.2417・刻み 0.5）に乗る。
    let evaluated = check_point_oracle(&inputs, &out, 0.0, 5.0, 0.5, "ISSUE-053 反例 C");
    assert!(evaluated > 50, "評価点が少なすぎる（{evaluated} 点）");
}

/// **ISSUE-053 反例 D**（3 リング・1 本目が `(4,0)` を 2 度通る**自己接触**）。
///
/// 期待値は**点オラクル**から導出する。既知の症状は点 `(1.6731, 1.2417)` /
/// `(2.1731, 0.7417)` で `want=false / got=true`（微小な外環にそれより広い穴が割り当たる）。
///
/// 殺す変異: 自己接触頂点での環の後継選択を誤る（`back` の符号・角度差の比較）・
/// 穴を包含関係でなく面積順だけで割り当てる・同一頂点を通る複数の環を 1 本に潰す。
#[test]
fn union_regression_issue053_counterexample_d_matches_point_oracle() {
    let inputs = vec![
        ring(&[
            (4.0, 1.0),
            (4.0, 0.0),
            (2.0, 0.0),
            (3.0, 2.0),
            (4.0, 0.0),
            (0.0, 1.0),
            (4.0, 3.0),
            (0.0, 0.0),
        ]),
        ring(&[(4.0, 0.0), (3.0, 0.0), (0.0, 4.0)]),
        ring(&[(2.0, 0.0), (2.0, 1.0), (0.0, 1.0)]),
    ];
    let out = union_rings(&inputs);

    assert_output_structure(&out);
    assert_no_fabricated_vertices(&inputs, &out);

    let evaluated = check_point_oracle(&inputs, &out, 0.0, 5.0, 0.5, "ISSUE-053 反例 D");
    assert!(evaluated > 50, "評価点が少なすぎる（{evaluated} 点）");
}

/// **ISSUE-053 反例 E**（3 リング・2 本目の 3 頂点リングと 3 本目の自己接触リングが
/// 同一頂点 `(2,4)` / `(3,3)` で会合する退化配置）。
///
/// これは**環再結合の後継選択規則（角度規則）の判別テスト**である。期待値は
/// **点オラクル**（入力リング群の even-odd 和）から導出し、**実装の出力からは一切取らない**。
/// 後継選択規則の一部が誤っていると、出力の内外が点オラクルと点ごとに食い違う。
/// 縛るのは §確定仕様 1（点ごとの一致）と §確定仕様 2（穴は外環の内部＝構造契約）。
///
/// 3 本目のリングは `(2,4)` と `(3,3)` をそれぞれ 2 度通る自己接触リング。契約上、入力リングは
/// even-odd・自己交差可なので有効な入力である。代表的な不一致点は `(1.1731, 2.7417)`。
///
/// 殺す変異: 自己接触頂点での後継選択（角度差の比較・`back` 方向の符号・同角度のタイ処理）を誤る・
/// 複数リングが会合する頂点で半辺を 1 本に潰す・穴を包含でなく面積順だけで割り当てる。
#[test]
fn union_regression_issue053_counterexample_e_matches_point_oracle() {
    let inputs = vec![
        ring(&[(3.0, 3.0), (6.0, 6.0), (3.0, 6.0)]),
        ring(&[(2.0, 4.0), (5.0, 7.0), (2.0, 7.0)]),
        ring(&[
            (2.0, 4.0),
            (3.0, 3.0),
            (0.0, 2.0),
            (2.0, 4.0),
            (3.0, 2.0),
            (3.0, 3.0),
            (0.0, 1.0),
        ]),
    ];
    let out = union_rings(&inputs);

    assert_output_structure(&out);
    assert_no_fabricated_vertices(&inputs, &out);

    // 代表点 `(1.1731, 2.7417)` は格子（オフセット 0.1731/0.2417・刻み 0.5）に乗り、
    // かつ境界近傍の除外に掛からない（＝評価対象である）ことを明示的に確かめる。
    let coords: Vec<Vec<(f64, f64)>> = inputs.iter().map(|r| ring_coords(r)).collect();
    let probe = (1.1731, 2.7417);
    let clearance = coords
        .iter()
        .flat_map(|r| edges_of(r).collect::<Vec<_>>())
        .map(|(a, b)| point_segment_distance(probe, a, b))
        .fold(f64::INFINITY, f64::min);
    assert!(
        clearance >= ORACLE_CLEARANCE,
        "代表点 {probe:?} が境界近傍として除外されている（clearance={clearance}）"
    );
    assert_eq!(
        output_inside(&out, probe),
        oracle_inside(&coords, probe),
        "代表点 {probe:?} の内外が点オラクルと違う"
    );

    // 格子全体（0..7 を覆う）でも点ごとに一致する。
    let evaluated = check_point_oracle(&inputs, &out, 0.0, 7.0, 0.5, "ISSUE-053 反例 E");
    assert!(evaluated > 50, "評価点が少なすぎる（{evaluated} 点）");
}

/// 決定的な擬似乱数（xorshift64*）。外部 crate に依存せず、**固定シードで完全に再現可能**。
struct Rng(u64);

impl Rng {
    fn new(seed: u64) -> Self {
        // 0 を避ける（xorshift は 0 で固定点になる）。
        Self(seed | 1)
    }

    fn next_u64(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x >> 12;
        x ^= x << 25;
        x ^= x >> 27;
        self.0 = x;
        x.wrapping_mul(0x2545_F491_4F6C_DD1D)
    }

    /// `[0, n)` の一様整数。
    fn below(&mut self, n: u64) -> u64 {
        self.next_u64() % n
    }
}

/// ランダムなリングを 1 本作る（軸平行矩形 / 三角形 / 自己交差・自己接触しうる自由多角形）。
/// 返すのは格子 `0..=5` 上の `(lon, lat)` 整数座標列（ISSUE-053 の探索空間と同じ）。
/// 自由多角形を 3/5 の比率で引く（頂点の一致・辺の共線重なりが起きやすく、退化配置に当たりやすい）。
fn random_ring_coords(rng: &mut Rng) -> Vec<(f64, f64)> {
    #[allow(clippy::cast_precision_loss)]
    fn c(v: u64) -> f64 {
        v as f64
    }
    match rng.below(5) {
        0 => {
            // 軸平行矩形（退化しないよう幅・高さを 1 以上にする）。
            let x0 = rng.below(4);
            let y0 = rng.below(4);
            let x1 = x0 + 1 + rng.below(5 - x0);
            let y1 = y0 + 1 + rng.below(5 - y0);
            vec![
                (c(x0), c(y0)),
                (c(x1), c(y0)),
                (c(x1), c(y1)),
                (c(x0), c(y1)),
            ]
        }
        1 => (0..3).map(|_| (c(rng.below(6)), c(rng.below(6)))).collect(),
        _ => {
            // 頂点 3〜8 の自由多角形（自己交差・自己接触・重複頂点を許す）。
            let n = 3 + rng.below(6);
            (0..n).map(|_| (c(rng.below(6)), c(rng.below(6)))).collect()
        }
    }
}

/// **ISSUE-053 §確定仕様 1・2 のランダム化差分テスト**（固定シード・有界反復・外部 crate 非依存）。
///
/// 格子 `0..=5` 上にリング 2〜5 本（頂点 3〜8・軸平行矩形/三角形/自己交差多角形の混合）を置き、
/// - 出力の内外（外環 ∧ ¬穴 の OR）が**点オラクル**（入力リングの even-odd 和）と一致すること、
/// - 出力の構造契約（穴が外環の内部・外環 CCW・穴 CW・非閉・面積降順）が保たれること、
/// - 出力頂点が入力辺の上に乗る（点を捏造しない）こと
///
/// を縛る。失敗時はシードとリング座標をメッセージに出すので、その配置を決定的テストへ昇格できる。
///
/// 殺す変異: 退化配置での穴の割り当て・環の後継選択・交点の量子化に入る、
/// 決定的テストが個別には捕まえきれない欠陥全般。
#[test]
fn union_randomized_differential_against_point_oracle() {
    const ITERATIONS: u64 = 20000;
    const BASE_SEED: u64 = 0x5155_1ED0_0053;

    for iter in 0..ITERATIONS {
        let seed = BASE_SEED ^ iter.wrapping_mul(0x9E37_79B9_7F4A_7C15);
        let mut rng = Rng::new(seed);
        let count = 2 + rng.below(4);
        let coords: Vec<Vec<(f64, f64)>> =
            (0..count).map(|_| random_ring_coords(&mut rng)).collect();
        let inputs: Vec<Vec<GeoPoint>> = coords.iter().map(|r| ring(r)).collect();

        let out = union_rings(&inputs);
        // 失敗時に配置を再現できるよう、シードと全リング座標を文脈に載せる。
        let label = format!("iter={iter} seed={seed:#x} rings={coords:?}");

        // 共通ヘルパ（構造契約）は label を受け取れないので、捕捉して配置を添えて投げ直す。
        let checked = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            for poly in &out {
                assert_polygon_structure(poly);
            }
            for w in out.windows(2) {
                let a = signed_area(&w[0].rings[0]).abs();
                let b = signed_area(&w[1].rings[0]).abs();
                assert!(a >= b - AREA_TOL, "外環の面積降順が崩れている（{a} < {b}）");
            }
            assert_no_fabricated_vertices(&inputs, &out);
            check_point_oracle(&inputs, &out, 0.0, 5.0, 0.5, &label);
        }));
        if let Err(payload) = checked {
            let detail = payload
                .downcast_ref::<String>()
                .cloned()
                .or_else(|| payload.downcast_ref::<&str>().map(|s| (*s).to_string()))
                .unwrap_or_else(|| "（メッセージ不明）".to_string());
            panic!("{label}\n{detail}");
        }
    }
}
