//! 平面 `(lon, lat)` 度の**多角形ユニオン**（部分食域の 3 領域合成・`docs/algorithms/11-path-partial-domain.md`
//! §11.6 / §11.7）。
//!
//! [`union_rings`] は入力リング集合の **even-odd 内部の和** U の境界を、外環＋穴の多角形列として返す。
//!
//! **アルゴリズムと選定理由（§11.7）**: 平面アレンジメント方式（辺分割 → 両側 membership による境界判定 →
//! 角度優先の環再結合。Weiler–Atherton 系）。Greiner–Hormann は退化（頂点上の交点・共有辺）で破綻し、
//! Martinez–Rueda（sweep line）は**入力が単純多角形であること**を前提にする。本 crate の用途（部分食域の
//! 「帯」）は**実データで自己交差する**ため、自己交差入力を構造的に許容する本方式を採る（自己交点も他リングとの
//! 交点と同列の分割点として扱えば同じ機構で吸収できる）。精度規約は ballpark（数十 km・§許容誤差）なので
//! sweep line の漸近性能は不要で、素朴な O(n²) 交点列挙で足りる。
//!
//! **退化ケースの方針（§11.7 の表）**:
//! - 自己交差リング → even-odd 充填により分離した領域になる（八の字は 2 多角形）。
//! - 共有辺・共線重なり → 量子化＋無向重複除去で 1 本に畳む。
//! - 頂点上の交点 → [`SNAP`] 量子化で既存頂点と一致し、特別扱いしない。
//! - 点接触（1 頂点のみ共有） → 角度優先の環再結合で分離した 2 多角形になる。
//! - 面積ゼロの環・頂点 3 未満のリング → 捨てる（捏造しない）。
//!
//! **残る近似（conventions §11「近似は明記」）**:
//! - 座標は `(lon, lat)` 度の**平面**として扱う（球面ではない）。**反子午線（±180）跨ぎ・極を含む領域は未対応**
//!   （[`crate::GeoPolygon::geojson_geometry`] の v1 制約と同じ）。
//! - 端点は [`SNAP`]（1e-9 度 ≈ 0.1 mm）グリッドへ量子化する。これより細かい特徴（極端に鋭い刺・1e-9 度
//!   未満の重なり）は解像しない。部分食域の許容誤差に対し 8 桁以上の余裕がある。
//! - 境界判定は中点からの法線オフセットによる両側 membership。オフセット幅は断片長と「中点から他の辺までの
//!   最短距離」に**適応**させ、近接して並走する別領域を誤って融合／消失させない。
//! - 交点は**辺ペアにつき 1 度だけ**計算して両辺で同じ量子化座標を共有する。辺ごとに独立に計算すると
//!   同一交点が 1 ulp ずれて別キーになり、環が閉じずに結果が黙って消える。
//! - 交点列挙は **O(n²)**。1 万辺規模を超える用途では sweep line への差し替えが要る。

use std::collections::{HashMap, HashSet};
use std::f64::consts::TAU;

use crate::geometry::{GeoPoint, GeoPolygon};

/// 端点量子化グリッド \[度\]（≈0.1 mm）。環再結合の端点一致を整数比較で厳密にする。
const SNAP: f64 = 1.0e-9;
/// ゼロ面積環の閾値（shoelace 2 倍値・\[度²\]）。これ未満の環は捏造せず捨てる。
const AREA2_EPS: f64 = 1.0e-14;
/// 平行判定の外積閾値。入力は [`SNAP`] 量子化済みゆえ、真に非平行な辺対の外積はこれを超える。
const CROSS_EPS: f64 = 1.0e-20;
/// 境界判定の法線オフセット係数。実効値は「断片長 × 本値」と「中点から他の辺までの最短距離 × 本値 × 2」の
/// 小さい方＝**入力スケールに適応**する（固定値だと、交差せずに近接して並走する別領域の境界辺を
/// 誤って内部と判定して捨てる）。
const PROBE_FRACTION: f64 = 0.25;
/// 共線頂点除去の判定閾値（隣接辺の正規化外積 = sin θ）。
const COLLINEAR_SIN_EPS: f64 = 1.0e-12;

/// 平面点 `[lon, lat]`（度）。
type P = [f64; 2];
/// [`SNAP`] グリッド上の整数キー（端点一致の厳密判定用）。
type Key = (i64, i64);

/// [`SNAP`] グリッドの整数キー。
///
/// `as i64` は切り捨てではなく `round()` 済みの値に対する変換で、入力は緯度経度（|値| ≤ 180）ゆえ
/// `180 / 1e-9 = 1.8e11` と i64 の範囲内に収まる（`clippy::cast_possible_truncation` を局所許可）。
#[allow(clippy::cast_possible_truncation)]
fn key(p: P) -> Key {
    ((p[0] / SNAP).round() as i64, (p[1] / SNAP).round() as i64)
}

/// [`SNAP`] グリッドへ丸めた点。
fn snap(p: P) -> P {
    [(p[0] / SNAP).round() * SNAP, (p[1] / SNAP).round() * SNAP]
}

fn cross(a: P, b: P) -> f64 {
    a[0] * b[1] - a[1] * b[0]
}

fn sub(a: P, b: P) -> P {
    [a[0] - b[0], a[1] - b[1]]
}

/// 非閉リングの符号付き面積 ×2（shoelace・CCW で正）。
fn signed_area2(ring: &[P]) -> f64 {
    let n = ring.len();
    let mut s = 0.0;
    for i in 0..n {
        let a = ring[i];
        let b = ring[(i + 1) % n];
        s += a[0] * b[1] - b[0] * a[1];
    }
    s
}

/// even-odd ray-casting point-in-polygon（+lon 方向への半直線との交差回数の偶奇）。
fn pip(ring: &[P], q: P) -> bool {
    let n = ring.len();
    if n < 3 {
        return false;
    }
    let mut inside = false;
    let mut j = n - 1;
    for i in 0..n {
        let (xi, yi) = (ring[i][0], ring[i][1]);
        let (xj, yj) = (ring[j][0], ring[j][1]);
        if (yi > q[1]) != (yj > q[1]) {
            let x_at = xi + (q[1] - yi) / (yj - yi) * (xj - xi);
            if x_at > q[0] {
                inside = !inside;
            }
        }
        j = i;
    }
    inside
}

/// 入力リングを量子化・連続重複点除去・非閉表現へ正規化する。頂点 3 未満は `None`。
fn normalize(ring: &[GeoPoint]) -> Option<Vec<P>> {
    let mut v: Vec<P> = ring
        .iter()
        .map(|g| snap([g.lon.degrees().0, g.lat.degrees().0]))
        .collect();
    v.dedup_by(|a, b| key(*a) == key(*b));
    while v.len() > 1 && key(v[0]) == key(v[v.len() - 1]) {
        v.pop();
    }
    if v.len() < 3 {
        None
    } else {
        Some(v)
    }
}

/// 共線と分かっている線分 `seg` 上に点 `p` が載るか（パラメータ範囲のみ判定）。
fn param_on_segment(p: P, seg: (P, P)) -> bool {
    let r = sub(seg.1, seg.0);
    let rr = r[0] * r[0] + r[1] * r[1];
    if rr == 0.0 {
        return false;
    }
    let d = sub(p, seg.0);
    let t = (d[0] * r[0] + d[1] * r[1]) / rr;
    (0.0..=1.0).contains(&t)
}

/// 線分 `a` と `b` の交点（[`SNAP`] 量子化済み）を返す。
///
/// **辺ペアにつき 1 度だけ呼び、返った点を両辺の分割点として共有すること。** 辺ごとに独立に計算すると、
/// 数学的には同一の交点が浮動小数で 1 ulp ずれ、[`SNAP`] の半グリッド境界をまたぐと別キーになり、
/// 交差点で鎖が繋がらず環が閉じない（結果が黙って消える）。
///
/// 共線重なりでは、重なり区間の端点（各線分の端点のうち他方に載るもの）を返す。これらは入力頂点そのもので
/// 既に量子化済みゆえ、両辺で自動的に同一キーになる。
fn intersection_points(a: (P, P), b: (P, P)) -> Vec<P> {
    let r = sub(a.1, a.0);
    let s = sub(b.1, b.0);
    let denom = cross(r, s);
    let qp = sub(b.0, a.0);
    if denom.abs() > CROSS_EPS {
        let t = cross(qp, s) / denom;
        let u = cross(qp, r) / denom;
        if (0.0..=1.0).contains(&t) && (0.0..=1.0).contains(&u) {
            return vec![snap([a.0[0] + r[0] * t, a.0[1] + r[1] * t])];
        }
        return Vec::new();
    }
    // 平行。共線でなければ交点なし。
    if cross(qp, r).abs() > CROSS_EPS {
        return Vec::new();
    }
    let mut out = Vec::new();
    for (p, seg) in [(b.0, a), (b.1, a), (a.0, b), (a.1, b)] {
        if param_on_segment(p, seg) {
            out.push(p);
        }
    }
    out
}

/// 点 `p` から線分 `a`–`b` までの最短距離。境界判定オフセットの適応幅に使う。
fn point_segment_distance(p: P, a: P, b: P) -> f64 {
    let r = sub(b, a);
    let rr = r[0] * r[0] + r[1] * r[1];
    let t = if rr > 0.0 {
        let d = sub(p, a);
        ((d[0] * r[0] + d[1] * r[1]) / rr).clamp(0.0, 1.0)
    } else {
        0.0
    };
    (p[0] - (a[0] + r[0] * t)).hypot(p[1] - (a[1] + r[1] * t))
}

/// 環から**共線の中間頂点**を落とす（直進の継続のみ・180° 折返しは残す）。安定するまで反復する。
/// 点を捏造せず、入力境界上の冗長な頂点だけを取り除く（辺共有・共線重なりで生じる）。
fn drop_collinear(mut ring: Vec<P>) -> Vec<P> {
    loop {
        let n = ring.len();
        if n <= 3 {
            // 三角形はこれ以上減らせない（共線な三角形は面積ゼロとして後段で捨てる）。
            return ring;
        }
        let mut removed = false;
        let mut out: Vec<P> = Vec::with_capacity(n);
        for i in 0..n {
            let prev = ring[(i + n - 1) % n];
            let cur = ring[i];
            let next = ring[(i + 1) % n];
            let v1 = sub(cur, prev);
            let v2 = sub(next, cur);
            let l1 = v1[0].hypot(v1[1]);
            let l2 = v2[0].hypot(v2[1]);
            let straight = l1 > 0.0
                && l2 > 0.0
                && cross(v1, v2).abs() <= COLLINEAR_SIN_EPS * l1 * l2
                && (v1[0] * v2[0] + v1[1] * v2[1]) > 0.0;
            if straight {
                removed = true;
            } else {
                out.push(cur);
            }
        }
        ring = out;
        if !removed {
            return ring;
        }
    }
}

/// 辺 `a`–`b` の中点から、その左右へ**入力スケールに適応した幅**でずらした 2 点 `(left, right)` を返す。
///
/// オフセット幅は「辺長 × [`PROBE_FRACTION`]」と「中点から**他の辺**までの最短距離 × [`PROBE_FRACTION`] × 2」の
/// 小さい方。固定幅にすると、交差せずに幅より近接して並走する別領域の境界辺を「内部に埋もれた辺」と
/// 誤判定して捨ててしまう（＝出力が黙って消える／誤って融合する）。中点からの距離が [`SNAP`] 以下の辺は
/// 「中点を含む辺」とみなして距離の候補から外す。
fn probe_sides(a: P, b: P, edges: &[(P, P)]) -> Option<(P, P)> {
    let d = sub(b, a);
    let len = d[0].hypot(d[1]);
    if len == 0.0 {
        return None;
    }
    let n = [-d[1] / len, d[0] / len];
    let mid = [(a[0] + b[0]) * 0.5, (a[1] + b[1]) * 0.5];
    let mut eps = len * PROBE_FRACTION;
    for e in edges {
        let dist = point_segment_distance(mid, e.0, e.1);
        if dist > SNAP {
            eps = eps.min(dist * PROBE_FRACTION * 2.0);
        }
    }
    Some((
        [mid[0] + n[0] * eps, mid[1] + n[1] * eps],
        [mid[0] - n[0] * eps, mid[1] - n[1] * eps],
    ))
}

/// 平面 `(lon, lat)` 度でリング集合の**ユニオン（和領域）**を計算し、外環＋穴の多角形列を返す。
///
/// - 各入力リングの充填規則は **even-odd**。`U` =「いずれかのリングの even-odd 内部」。
/// - 入力リングは閉じていても非閉でもよく、**向き（CW/CCW）・リングの並び順は結果に影響しない**。
///   頂点 3 未満・面積ゼロのリングは無視する。連続重複点は無視する。
/// - 返る各 [`GeoPolygon`] は `rings[0]` = 外環（**CCW**・shoelace > 0）・`rings[1..]` = その多角形に属する
///   穴（**CW**）。いずれも**非閉表現**（先頭 != 末尾）で、共線の中間頂点は落としてある。
///   並びは**外環の面積降順**（決定的）。
/// - 入力が空／全て退化なら**空 `Vec`**（多角形を捏造しない）。出力頂点は必ず入力リングの頂点か
///   入力辺同士の交点であり、入力境界の外に点を作らない。
///
/// 前提・退化ケースの方針・残る近似は本モジュールの doc および
/// `docs/algorithms/11-path-partial-domain.md` §11.7 を参照。
pub fn union_rings(rings: &[Vec<GeoPoint>]) -> Vec<GeoPolygon> {
    // 面積ゼロ（一直線）のリングはここでは落とさない: **自己交差リングは shoelace が相殺して 0 になりうる**
    // （対称な八の字など）。退化リングは even-odd 内部を持たないので、境界判定（手順 4）で自然に消える。
    let regions: Vec<Vec<P>> = rings.iter().filter_map(|r| normalize(r)).collect();
    if regions.is_empty() {
        return Vec::new();
    }

    // 1. 全リングの全辺。
    let mut edges: Vec<(P, P)> = Vec::new();
    for region in &regions {
        let n = region.len();
        for i in 0..n {
            edges.push((region[i], region[(i + 1) % n]));
        }
    }

    // 2. 交点を**辺ペアにつき 1 度だけ**求め、両辺の分割点として同一座標を共有する（同一リング内の
    //    ペアも含む＝自己交差を吸収）。辺ごとに独立計算すると 1 ulp ずれて環が閉じない。
    let mut splits: Vec<Vec<P>> = vec![Vec::new(); edges.len()];
    for i in 0..edges.len() {
        for j in (i + 1)..edges.len() {
            for p in intersection_points(edges[i], edges[j]) {
                splits[i].push(p);
                splits[j].push(p);
            }
        }
    }

    // 分割点を辺上のパラメータ順に並べ、**点そのもの**（再計算しない）を繋いで断片にする。
    let mut fragments: Vec<(P, P)> = Vec::new();
    for (i, e) in edges.iter().enumerate() {
        let r = sub(e.1, e.0);
        let rr = r[0] * r[0] + r[1] * r[1];
        let mut ordered: Vec<(f64, P)> = vec![(0.0, e.0), (1.0, e.1)];
        for p in &splits[i] {
            let t = if rr > 0.0 {
                let d = sub(*p, e.0);
                ((d[0] * r[0] + d[1] * r[1]) / rr).clamp(0.0, 1.0)
            } else {
                0.0
            };
            ordered.push((t, *p));
        }
        ordered.sort_by(|a, b| a.0.total_cmp(&b.0));
        let mut prev = ordered[0].1;
        for (_, p) in ordered.into_iter().skip(1) {
            if key(p) != key(prev) {
                fragments.push((prev, p));
                prev = p;
            }
        }
    }

    // 3. 無向端点対で重複除去（共有辺・共線重なり・完全重複リング）。
    let mut seen: HashSet<(Key, Key)> = HashSet::new();
    let mut unique: Vec<(P, P)> = Vec::new();
    for (a, b) in fragments {
        let (ka, kb) = (key(a), key(b));
        let k = if ka <= kb { (ka, kb) } else { (kb, ka) };
        if seen.insert(k) {
            unique.push((a, b));
        }
    }

    // 4. 境界判定（中点の両側 membership）＋「内部が左」になるよう向き付け。
    let inside = |q: P| regions.iter().any(|r| pip(r, q));
    let mut half: Vec<(P, P)> = Vec::new();
    for (a, b) in unique {
        let Some((left_probe, right_probe)) = probe_sides(a, b, &edges) else {
            continue;
        };
        let left = inside(left_probe);
        let right = inside(right_probe);
        if left == right {
            continue; // 内部に埋もれた辺／領域外の辺。
        }
        if left {
            half.push((a, b));
        } else {
            half.push((b, a));
        }
    }
    if half.is_empty() {
        return Vec::new();
    }

    // 5. 角度優先で環を再結合。頂点 v へ半辺 e で入ったとき、ē（e の逆向き）から**時計回り**に
    //    最初に現れる未使用の出辺を後継とする（面の後継半辺の標準規則）。点接触はこれで 2 環に分かれる。
    let mut out_edges: HashMap<Key, Vec<usize>> = HashMap::new();
    for (idx, (a, _)) in half.iter().enumerate() {
        out_edges.entry(key(*a)).or_default().push(idx);
    }
    let mut used = vec![false; half.len()];
    let mut cycles: Vec<Vec<P>> = Vec::new();
    for start in 0..half.len() {
        if used[start] {
            continue;
        }
        let start_key = key(half[start].0);
        let mut ring: Vec<P> = Vec::new();
        let mut cur = start;
        let mut closed = false;
        loop {
            used[cur] = true;
            ring.push(half[cur].0);
            let v = half[cur].1;
            if key(v) == start_key {
                closed = true;
                break;
            }
            let din = sub(v, half[cur].0);
            let back = (-din[1]).atan2(-din[0]);
            let Some(candidates) = out_edges.get(&key(v)) else {
                break;
            };
            let mut best: Option<(f64, usize)> = None;
            for &c in candidates {
                if used[c] {
                    continue;
                }
                let dout = sub(half[c].1, half[c].0);
                // ē から出辺までの**時計回り**の回転角を `(0, TAU]` へ正規化する。`back` も `atan2` も
                // `(-π, π]` なので差は `(-2π, 2π)`。TAU を足す繰り返しで必ず `(0, TAU]` に収まる
                // （上側への折返しは起こり得ないので不要）。
                let mut diff = back - dout[1].atan2(dout[0]);
                while diff <= 0.0 {
                    diff += TAU;
                }
                let better = match best {
                    None => true,
                    Some((d, _)) => diff < d,
                };
                if better {
                    best = Some((diff, c));
                }
            }
            match best {
                Some((_, next)) => cur = next,
                None => break,
            }
        }
        if closed {
            let simplified = drop_collinear(ring);
            if simplified.len() >= 3 {
                cycles.push(simplified);
            }
        }
    }

    // 6. 外環（面積>0・CCW）と穴（面積<0・CW）へ分離。
    let mut outers: Vec<Vec<P>> = Vec::new();
    let mut holes: Vec<Vec<P>> = Vec::new();
    for c in cycles {
        let a2 = signed_area2(&c);
        if a2.abs() <= AREA2_EPS {
            continue;
        }
        if a2 > 0.0 {
            outers.push(c);
        } else {
            holes.push(c);
        }
    }
    if outers.is_empty() {
        return Vec::new();
    }
    outers.sort_by(|a, b| signed_area2(b).abs().total_cmp(&signed_area2(a).abs()));

    // 7. 穴を、その空洞側代表点を含む**最小面積**の外環へ割り当てる（穴の中の島の入れ子に対応）。
    let mut assigned: Vec<Vec<Vec<P>>> = outers.iter().map(|_| Vec::new()).collect();
    for h in holes {
        // 穴は CW（内部が左）なので、空洞は走査方向の**右**側。
        let Some((_, probe)) = probe_sides(h[0], h[1], &edges) else {
            continue;
        };
        let mut best: Option<usize> = None;
        for (i, o) in outers.iter().enumerate() {
            if !pip(o, probe) {
                continue;
            }
            let better = match best {
                None => true,
                Some(j) => signed_area2(o).abs() < signed_area2(&outers[j]).abs(),
            };
            if better {
                best = Some(i);
            }
        }
        if let Some(i) = best {
            assigned[i].push(h);
        }
    }

    // 8. GeoPolygon へ（緯度経度域外は起こらないが、変換できない環は捏造せず落とす）。
    let to_geo = |ring: &Vec<P>| -> Option<Vec<GeoPoint>> {
        ring.iter()
            .map(|p| GeoPoint::from_degrees(p[1], p[0]).ok())
            .collect()
    };
    let mut out = Vec::new();
    for (i, o) in outers.iter().enumerate() {
        let Some(outer) = to_geo(o) else {
            continue;
        };
        let mut poly_rings = vec![outer];
        for h in &assigned[i] {
            if let Some(hole) = to_geo(h) {
                poly_rings.push(hole);
            }
        }
        out.push(GeoPolygon::new(poly_rings));
    }
    out
}
