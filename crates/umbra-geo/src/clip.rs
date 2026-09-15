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

use crate::geometry::{EnclosedPole, GeoPoint, GeoPolygon};

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

/// 点 `q` がリング `ring` の**閉内部**（even-odd 内部または境界上）にあるか。
///
/// 穴の帰属判定（手順 7）で使う。穴の頂点は外環の頂点・辺の上に**乗る**ことが普通にあり
/// （点接触・共有辺）、そこで [`pip`] の開内部判定は境界上の点を落としてしまう。
fn point_in_or_on(ring: &[P], q: P) -> bool {
    let n = ring.len();
    if n < 3 {
        return false;
    }
    for i in 0..n {
        let seg = (ring[i], ring[(i + 1) % n]);
        let r = sub(seg.1, seg.0);
        let d = sub(q, seg.0);
        // 辺上（共線かつパラメータ範囲内）なら境界上。
        if (r[0] * d[1] - r[1] * d[0]).abs() <= SNAP * (r[0].abs() + r[1].abs() + 1.0)
            && param_on_segment(q, seg)
        {
            return true;
        }
    }
    pip(ring, q)
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
/// 共線ゆえパラメータは**成分が大きい側の軸 1 本**で決まる（縦線分でも 0 除算にならない）。
fn param_on_segment(p: P, seg: (P, P)) -> bool {
    let r = sub(seg.1, seg.0);
    let d = sub(p, seg.0);
    let axis = if r[0].abs() >= r[1].abs() { 0 } else { 1 };
    if r[axis] == 0.0 {
        return false;
    }
    let t = d[axis] / r[axis];
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

/// 閉歩行 `ring` を、同一（量子化）頂点を 2 度通る箇所で単純環に分割する（pinch 分割・§11.7 手順 7）。
/// 頂点 v を 2 度通る歩行 `[.., v, X.., v, Y..]` は、部分歩行 `[v, X..]` が閉環をなすので切り出し、残り
/// `[.., v, Y..]` に対して繰り返す。重複が無ければ `ring` をそのまま 1 要素で返す。点を捏造しない。
fn split_at_repeated_vertices(ring: Vec<P>) -> Vec<Vec<P>> {
    let mut out = Vec::new();
    let mut stack: Vec<P> = Vec::with_capacity(ring.len());
    let mut seen: HashMap<Key, usize> = HashMap::new();
    for p in ring {
        let k = key(p);
        if let Some(&start) = seen.get(&k) {
            // v の前回出現位置から現在までが閉環。切り出し、v は残りの歩行の頂点として残す。
            let sub: Vec<P> = stack.drain(start..).collect();
            for q in &sub {
                seen.remove(&key(*q));
            }
            out.push(sub);
        }
        seen.insert(k, stack.len());
        stack.push(p);
    }
    if !stack.is_empty() {
        out.push(stack);
    }
    out
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
/// - **経度フレーム（ISSUE-052）**: 反子午線を跨ぐ入力は、継ぎ目が空き経度区間に来るよう内部で経度を
///   回してから解き、結果を逆回転で戻す（跨がない入力は回さない＝出力はバイト不変）。
///   **回しても跨ぎが消えない入力（領域が全経度を覆う＝極を囲む）だけは従来どおり平面で解く**ので、
///   その場合に限り返る外環は**呼び出し側のフレームで反子午線を跨ぎうる**。跨ぐ外環の
///   **平面 shoelace は領域の面積ではない**（実測で真値の約 17 倍になる）ので、面積として使わないこと。
///   `CCW`・面積降順の保証は**内部の回転後フレームで**成立したものを保持している。
///   跨ぐ可能性の判定と GeoJSON 化は [`crate::GeoPolygon::geojson_geometry_with_pole`] が行う。
/// - 入力が空／全て退化なら**空 `Vec`**（多角形を捏造しない）。出力頂点は必ず入力リングの頂点か
///   入力辺同士の交点であり、入力境界の外に点を作らない。
///
/// 前提・退化ケースの方針・残る近似は本モジュールの doc および
/// `docs/algorithms/11-path-partial-domain.md` §11.7 を参照。
pub fn union_rings(rings: &[Vec<GeoPoint>]) -> Vec<GeoPolygon> {
    // 面積ゼロ（一直線）のリングはここでは落とさない: **自己交差リングは shoelace が相殺して 0 になりうる**
    // （対称な八の字など）。退化リングは even-odd 内部を持たないので、境界判定（手順 4）で自然に消える。
    let mut regions: Vec<Vec<P>> = rings.iter().filter_map(|r| normalize(r)).collect();
    if regions.is_empty() {
        return Vec::new();
    }

    // 0. **経度フレームの正規化**（ISSUE-052）。`lon` を単なる平面座標として扱うため、反子午線を跨ぐ
    //    入力は「地球を逆走する辺」を持ち、**別の図形**になる（実測: 同一領域が面積 3490 対 正 200）。
    //    跨ぐ辺があれば、全頂点の経度の**最大の空き区間**の中央へ ±180 の継ぎ目が来るよう回してから解き、
    //    最後に逆回転で戻す。跨ぐ辺が無ければ回さない（従来の出力はバイト不変）。
    //    回しても跨ぐ辺が残る（＝領域が全経度を覆う＝極を囲む）場合は**回さずに従来どおり**解く
    //    （結果は従来と同じく未保証だが、勝手に別の答えを作らない）。
    let offset = longitude_frame_offset(&regions);
    if offset != 0.0 {
        for region in &mut regions {
            for p in region.iter_mut() {
                p[0] = wrap_longitude(p[0] + offset);
            }
        }
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
            // 5'. pinch 分割: 面の境界としての閉歩行は、同じ面に属する 2 つの穴が 1 点で接する場合や
            //     穴が外環に 1 点で接する場合に、その頂点を 2 度通る**自己接触環**になる（面は 1 つでも
            //     境界成分は複数）。同一頂点を 2 度通る環は、その頂点で単純環に分割する（§11.7 手順 7）。
            for loop_ in split_at_repeated_vertices(ring) {
                let simplified = drop_collinear(loop_);
                if simplified.len() >= 3 {
                    cycles.push(simplified);
                }
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

    // 7. 穴を、**穴リング全体を含む**最小面積の外環へ割り当てる（穴の中の島の入れ子に対応）。
    //
    //    空洞側の代表点 1 点だけで判定すると、**代表点をたまたま含む微小な外環**が最小面積として
    //    選ばれ、穴が外環より大きいまま割り当たる（ISSUE-053 の反例 C/D＝自己接触配置で実測）。
    //    穴 H が外環 O の穴であるためには H ⊂ O が必要なので、**H の全頂点が O の閉内部**
    //    （内部または境界上）であることを要求する。面積の必要条件 |O| > |H| も併せて課す。
    //
    //    含む外環が 1 つも無い CW 環は**穴として採用しない**（非有界面の境界成分＝偽の環）。
    //    幾何的には H ⊂ O なら必ず |O| > |H| なので、正しい外環が面積条件で落ちることはない。
    //    量子化 SNAP のぶんの誤差は `point_in_or_on` の境界許容が吸収する。
    let mut assigned: Vec<Vec<Vec<P>>> = outers.iter().map(|_| Vec::new()).collect();
    for h in holes {
        let h_area2 = signed_area2(&h).abs();
        let mut best: Option<usize> = None;
        for (i, o) in outers.iter().enumerate() {
            if signed_area2(o).abs() <= h_area2 {
                continue; // 穴以下の面積の外環は H を含み得ない。
            }
            if !h.iter().all(|&q| point_in_or_on(o, q)) {
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

    // 8'. 経度フレームを戻す（手順 0 で回した場合のみ）。緯度には触れない（ISSUE-052 §確定仕様 5）。
    if offset != 0.0 {
        for ring in outers.iter_mut() {
            for p in ring.iter_mut() {
                p[0] = wrap_longitude(p[0] - offset);
            }
        }
        for holes in assigned.iter_mut() {
            for ring in holes.iter_mut() {
                for p in ring.iter_mut() {
                    p[0] = wrap_longitude(p[0] - offset);
                }
            }
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

// ============================================================
// 反子午線分割（(3g)・§11.8）
// ============================================================

/// 跨ぎ判定の閾値（度）。`|Δlon| = 180` ちょうどは跨ぎとしない（測度ゼロ境界・`GeoLine` と同一規約）。
const ANTIMERIDIAN_DELTA: f64 = 180.0;

/// 閉リング（`[経度, 緯度]` 列・**非閉表現**＝末尾に先頭を重複させない）を反子午線で分割し、
/// 東半球側・西半球側の閉環に切り分ける（§11.8(b)）。跨ぎが無ければ入力をそのまま 1 要素で返す。
///
/// 跨ぎ点の緯度は子午線上に線形補間し（東進 `t=(180−lon1)/(360+Δlon)`・西進 `t=(lon1+180)/(360−Δlon)`、
/// `lat_c = lat1 + t·(lat2−lat1)`）、弧を子午線上で結んで閉じる。結線の向きは**リングの実際の向きが
/// 内部を左に保つ**ように取る（実 CCW なら東側は北向き・西側は南向き。実 CW なら反転）。出力の環向き
/// 正規化は**呼び出し側が分割の後に**行う（§11.8(c)）。点を捏造しない（子午線上の頂点は全て補間点＝境界上）。
pub(crate) fn split_ring_at_antimeridian(ring: &[P], pole: Option<EnclosedPole>) -> Vec<Vec<P>> {
    if ring.len() < 3 {
        return vec![ring.to_vec()];
    }
    // 極を囲むリング（§11.8(d) / ISSUE-051）。閉リングの跨ぎ回数は経度の巻き数に等しく、
    // **奇数回＝経度が一周する＝極を囲む**。このとき弧は対を成さず子午線だけでは閉じられない。
    // `pole` が与えられていれば**その極を通して閉じる**（ISSUE-051）。無ければ従来どおり
    // **元のリングをそのまま返す**（面積を失わず・捏造もしない）。
    if antimeridian_crossings(ring) % 2 == 1 && pole.is_none() {
        return vec![ring.to_vec()];
    }
    // **リングの向きは結線に不要**（ISSUE-051）。子午線の結線はパリティ（極から数えた内外の反転）で
    // 決まり、リングをどちら回りに辿っていても同じ対応になる。従来はここで経度を連続化して向きを
    // 求めていたが、パリティ規則への一般化で不要になったので削除した（死コードを残さない）。
    // 1. 跨ぎ位置で弧に切る。閉曲線なので最後の辺（末尾→先頭）も走査する。
    //    弧は `(始点が子午線上か, 点列)` を持ち、最初の弧は途中から始まりうるので最後に連結する。
    let n = ring.len();
    let mut arcs: Vec<Vec<P>> = Vec::new();
    let mut current: Vec<P> = vec![ring[0]];
    let mut crossed = false;
    for i in 0..n {
        let a = ring[i];
        let b = ring[(i + 1) % n];
        let delta = b[0] - a[0];
        if delta < -ANTIMERIDIAN_DELTA {
            // 東進（+180 を越える）: a → +180 → −180 → b。
            let t = (180.0 - a[0]) / (360.0 + delta);
            let lat_c = a[1] + t * (b[1] - a[1]);
            current.push([180.0, lat_c]);
            arcs.push(std::mem::take(&mut current));
            current.push([-180.0, lat_c]);
            crossed = true;
        } else if delta > ANTIMERIDIAN_DELTA {
            // 西進（−180 を越える）: a → −180 → +180 → b。
            let t = (a[0] + 180.0) / (360.0 - delta);
            let lat_c = a[1] + t * (b[1] - a[1]);
            current.push([-180.0, lat_c]);
            arcs.push(std::mem::take(&mut current));
            current.push([180.0, lat_c]);
            crossed = true;
        }
        // b は次の辺の始点。最後の辺では先頭に戻るので積まない。
        if i + 1 < n {
            current.push(b);
        }
    }
    if !crossed {
        return vec![ring.to_vec()];
    }
    // 最後の弧は先頭の弧の続き（閉曲線の起点は任意）。先頭へ連結する。
    if !current.is_empty() {
        if arcs.is_empty() {
            return vec![ring.to_vec()];
        }
        let head = arcs.remove(0);
        current.extend(head);
        arcs.insert(0, current);
    }

    // 2. 子午線上の端点を**緯度順に並べて隣どうしを対にする**（§11.8(b)・ISSUE-051 で一般化）。
    //    どの子午線区間が領域の内部かは**パリティ**で決まる: 極から緯度を上げていくと、境界と交わる
    //    たびに内外が反転する。極を囲まない領域は `lat=−90` が外部なので下から (1,2),(3,4),… が内部。
    //    極を囲む領域はその極が内部なので、**極を端点リストの極側に挿入**してから同じ規則で対にする。
    //    これで偶数跨ぎ（従来）と奇数跨ぎ（極を囲む）が同じ規則で扱える。
    //
    //    **「両端が同じ子午線に載る弧は自己閉じ」という規則は誤り**（ISSUE-051 で判明）: 子午線の
    //    結線が許されるのはその区間が**内部**のときだけで、弧の帳簿ではなくパリティが決める。
    //    誤った規則は外部の帯を内部として塗り、同時に内部の帯を落とす。
    // 極が効くのは**奇数跨ぎ（経度が一周する＝極を囲む）ときだけ**。偶数跨ぎでは極の指定は
    // 結果に影響しない（§確定仕様 4）。
    let encircles_pole = antimeridian_crossings(ring) % 2 == 1;
    let pole_lat = pole.filter(|_| encircles_pole).map(|p| match p {
        EnclosedPole::North => 90.0_f64,
        EnclosedPole::South => -90.0_f64,
    });
    let Some(links) = pair_meridian_endpoints(&arcs, pole_lat) else {
        return vec![ring.to_vec()];
    };
    let Some(out) = walk_arc_cycles(&arcs, &links, pole_lat) else {
        return vec![ring.to_vec()];
    };
    if out.is_empty() {
        vec![ring.to_vec()]
    } else {
        out
    }
}

/// 弧の端点の識別子。`(弧の番号, 終端か)`。
type EndpointId = (usize, bool);

/// 子午線上の端点を緯度順に並べ、隣どうしを対にする（ISSUE-051 §確定仕様 3 の一般化規則）。
///
/// 極を囲む場合（`pole_lat` が `Some`）は、**その極側の端**に極を表す番兵を挿入してから対にする。
/// 対にできない（端点数の偶奇が合わない）場合は `None`＝呼び出し側が分割を諦める。
/// 返り値は端点 → 相方の対応表で、相方が `None` の端点は**極へ接続する**。
#[allow(clippy::type_complexity)]
fn pair_meridian_endpoints(
    arcs: &[Vec<P>],
    pole_lat: Option<f64>,
) -> Option<HashMap<EndpointId, Option<EndpointId>>> {
    let mut links: HashMap<EndpointId, Option<EndpointId>> = HashMap::new();
    for east_side in [true, false] {
        // この子午線に載る端点を集める（緯度・識別子）。
        let mut points: Vec<(f64, EndpointId)> = Vec::new();
        for (i, arc) in arcs.iter().enumerate() {
            for is_end in [false, true] {
                let p = if is_end { arc[arc.len() - 1] } else { arc[0] };
                if p[0].abs() < ANTIMERIDIAN_DELTA {
                    continue;
                }
                if (p[0] > 0.0) == east_side {
                    points.push((p[1], (i, is_end)));
                }
            }
        }
        if points.is_empty() {
            continue;
        }
        points.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
        // 極を囲むなら、極側の端に「極へ抜ける」枠を 1 つ足す（南極なら先頭・北極なら末尾）。
        let mut slots: Vec<Option<EndpointId>> = points.iter().map(|&(_, id)| Some(id)).collect();
        match pole_lat {
            Some(lat) if lat < 0.0 => slots.insert(0, None),
            Some(_) => slots.push(None),
            None => {}
        }
        if slots.len() % 2 != 0 {
            return None; // 対にならない（不正入力）。捏造せず諦める。
        }
        for pair in slots.chunks(2) {
            match (pair[0], pair[1]) {
                (Some(a), Some(b)) => {
                    links.insert(a, Some(b));
                    links.insert(b, Some(a));
                }
                // 極に接する枠。相方が極であることを `None` で表す。
                (Some(a), None) | (None, Some(a)) => {
                    links.insert(a, None);
                }
                (None, None) => return None,
            }
        }
    }
    Some(links)
}

/// 対応表に従って弧を辿り、閉環を取り出す（ISSUE-051）。
///
/// 弧の終点から相方の端点へ渡り、その弧を順方向に辿る。相方が極（`None`）なら**極の 2 頂点**
/// （`(±180, ±90)`）を挿入して反対側の子午線の極接続端点へ渡る。全弧をちょうど 1 回ずつ使い切れない
/// 場合は `None`＝分割を諦める（面積を黙って失わない・捏造しない）。
fn walk_arc_cycles(
    arcs: &[Vec<P>],
    links: &HashMap<EndpointId, Option<EndpointId>>,
    pole_lat: Option<f64>,
) -> Option<Vec<Vec<P>>> {
    // 極へ抜ける端点は両子午線に 1 つずつ。互いの相手として繋ぐ。
    let pole_ends: Vec<EndpointId> = {
        let mut v: Vec<EndpointId> = links
            .iter()
            .filter_map(|(&id, &to)| if to.is_none() { Some(id) } else { None })
            .collect();
        v.sort_unstable();
        v
    };
    if pole_lat.is_some() && pole_ends.len() != 2 {
        return None;
    }

    let mut used = vec![false; arcs.len()];
    let mut rings = Vec::new();
    for start in 0..arcs.len() {
        if used[start] {
            continue;
        }
        let mut ring: Vec<P> = Vec::new();
        let mut cur = start;
        let mut closed = false;
        for _ in 0..=arcs.len() {
            if used[cur] {
                return None; // 同じ弧を 2 度使う＝対応表が壊れている。
            }
            used[cur] = true;
            ring.extend(arcs[cur].iter().copied());
            // 終点の相方へ渡る。
            let partner = links.get(&(cur, true)).copied();
            let next_start = match partner {
                Some(Some(id)) => id,
                Some(None) => {
                    // 極を経由して反対側の子午線の極接続端点へ。
                    let lat = pole_lat?;
                    let here = arcs[cur][arcs[cur].len() - 1][0];
                    let other = *pole_ends.iter().find(|&&id| id != (cur, true))?;
                    let other_lon = arc_endpoint(arcs, other)[0];
                    ring.push([here, lat]);
                    ring.push([other_lon, lat]);
                    other
                }
                None => return None,
            };
            // 相方が弧の始点なら順方向に続く。終点なら向きが不整合。
            if next_start.1 {
                return None;
            }
            if next_start.0 == start {
                closed = true;
                break;
            }
            cur = next_start.0;
        }
        if !closed || ring.len() < 3 {
            return None;
        }
        rings.push(ring);
    }
    Some(rings)
}

/// 端点識別子から実際の座標を引く。
fn arc_endpoint(arcs: &[Vec<P>], id: EndpointId) -> P {
    let arc = &arcs[id.0];
    if id.1 {
        arc[arc.len() - 1]
    } else {
        arc[0]
    }
}

/// 多角形（外環＋穴・`[経度, 緯度]` の**非閉**列）を反子午線で分割し、GeoJSON 用の多角形列
/// （各要素は `[外環, 穴…]`）を返す（§11.8(c)）。
///
/// 手順: (1) 各リングを [`split_ring_at_antimeridian`] で分割（跨がないリングはそのまま）、
/// (2) 各断片の環向きを正規化（外環由来 CCW・穴由来 CW。**分割の後**に行う＝跨ぐリングの平面 shoelace は
/// 意味を持たないため）、(3) 穴の断片を、それを含む外環断片のうち**最小面積**のものへ割り当てる
/// （§11.7 手順 8 と同一規則・含む外環が無ければ捨てる）。多角形は外環の面積降順で返す。
pub(crate) fn split_polygon_at_antimeridian(
    rings: &[Vec<P>],
    pole: Option<EnclosedPole>,
) -> Vec<Vec<Vec<P>>> {
    let Some(outer_ring) = rings.first() else {
        return Vec::new();
    };
    // (1)(2) 外環: 分割して CCW へ正規化。
    let mut outers: Vec<Vec<P>> = split_ring_at_antimeridian(outer_ring, pole)
        .into_iter()
        .map(|mut r| {
            if signed_area2(&r) < 0.0 {
                r.reverse();
            }
            r
        })
        .collect();
    // (1)(2) 穴: 分割して CW へ正規化。
    let holes: Vec<Vec<P>> = rings[1..]
        .iter()
        .flat_map(|h| split_ring_at_antimeridian(h, None))
        .map(|mut r| {
            if signed_area2(&r) > 0.0 {
                r.reverse();
            }
            r
        })
        .collect();

    // 外環を面積降順に並べる（§11.7 手順 8 と同一・バイト安定な出力のため）。
    outers.sort_by(|a, b| {
        signed_area2(b)
            .abs()
            .partial_cmp(&signed_area2(a).abs())
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    let mut polygons: Vec<Vec<Vec<P>>> = outers.into_iter().map(|o| vec![o]).collect();
    // (3) 穴を、それを含む最小面積の外環断片へ割り当てる。
    for hole in holes {
        // 判定点は**子午線から外れた穴の頂点**を使う。分割後の穴は端点が子午線上（±180）に載り、外環断片の
        // 境界と一致するため、その頂点では `pip` が不定になる。子午線外の頂点は（穴が外環に内包される限り）
        // 外環断片の内部にある。頂点重心は外環が凹（C 字）のとき切り欠きに落ちて誤判定しうるので使わない。
        // 全頂点が子午線上（退行）の場合のみ重心へ退避する。
        // 全頂点が子午線上の穴は面積ゼロの退行なので、判定せず捨てる（捏造しない）。
        let Some(probe) = hole
            .iter()
            .copied()
            .find(|p| p[0].abs() < ANTIMERIDIAN_DELTA)
        else {
            continue;
        };
        let mut best: Option<usize> = None;
        for (i, poly) in polygons.iter().enumerate() {
            if !pip(&poly[0], probe) {
                continue;
            }
            let take = match best {
                None => true,
                Some(j) => signed_area2(&poly[0]).abs() < signed_area2(&polygons[j][0]).abs(),
            };
            if take {
                best = Some(i);
            }
        }
        if let Some(i) = best {
            polygons[i].push(hole);
        }
    }
    polygons
}

/// 閉リングが反子午線を跨ぐ回数（§11.8(b) の跨ぎ判定・末尾→先頭の辺も含む）。
/// 偶奇が経度の巻き数の偶奇に一致するので、極を囲むかの判別（§11.8(d)）にも使う。
fn antimeridian_crossings(ring: &[P]) -> usize {
    let n = ring.len();
    if n < 2 {
        return 0;
    }
    (0..n)
        .filter(|&i| {
            let delta = ring[(i + 1) % n][0] - ring[i][0];
            !(-ANTIMERIDIAN_DELTA..=ANTIMERIDIAN_DELTA).contains(&delta)
        })
        .count()
}

/// リング列のいずれかが反子午線を跨ぐか。跨がない入力では分割経路に入らず、従来の出力を
/// **バイト不変**に保つために使う。
pub(crate) fn crosses_antimeridian(rings: &[Vec<P>]) -> bool {
    rings.iter().any(|ring| antimeridian_crossings(ring) > 0)
}

/// 経度を `[-180, 180)` へ折り返す（フレーム回転用・ISSUE-052）。
fn wrap_longitude(lon: f64) -> f64 {
    let mut v = (lon + 180.0) % 360.0;
    if v < 0.0 {
        v += 360.0;
    }
    v - 180.0
}

/// 入力が反子午線を跨ぐとき、継ぎ目を**空き経度区間**の中央へ移す回転量を返す（ISSUE-052）。
///
/// 跨ぐ辺（`|Δlon| > 180`）が 1 本も無ければ `0.0`＝回さない（出力はバイト不変）。
///
/// **候補は幅の広い空き区間から順に試し、回した結果に跨ぐ辺が 1 本も残らない最初のものを採る。**
/// 空き区間は**頂点**の並びから求めるので、「頂点は無いが**辺が横切っている**」区間が最大になりうる
/// （実装レビュー指摘）。そこへ継ぎ目を置くとその辺が跨ぎに変わって検証に失敗するので、1 候補で
/// 諦めると**本当に直すべき領域まで未修正のまま**になる。よって候補を順に試す。
/// どの候補でも跨ぎが残る場合（領域が全経度を覆う＝極を囲む）は `0.0`＝**従来どおり**解かせる
/// （別の答えを作らない）。同幅の区間は経度の小さい側から試す（決定的）。
fn longitude_frame_offset(regions: &[Vec<P>]) -> f64 {
    if !crosses_antimeridian(regions) {
        return 0.0;
    }
    let mut lons: Vec<f64> = regions.iter().flatten().map(|p| p[0]).collect();
    lons.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    lons.dedup();
    if lons.len() < 2 {
        return 0.0;
    }
    // 円周上の隣接する空き区間（末尾→先頭の巻き戻りも 1 区間）を、**幅の降順**に候補化する。
    let mut candidates: Vec<(f64, f64)> = (0..lons.len())
        .map(|i| {
            let (a, b) = (lons[i], lons[(i + 1) % lons.len()]);
            let gap = if i + 1 == lons.len() {
                b + 360.0 - a
            } else {
                b - a
            };
            // 区間の中央が継ぎ目（±180）に来る回転量。
            (gap, wrap_longitude(180.0 - wrap_longitude(a + gap / 2.0)))
        })
        .collect();
    candidates.sort_by(|x, y| y.0.partial_cmp(&x.0).unwrap_or(std::cmp::Ordering::Equal));

    // 候補を順に試し、回した結果に跨ぐ辺が 1 本も残らない最初のものを採る。
    for (_, offset) in candidates {
        let rotated: Vec<Vec<P>> = regions
            .iter()
            .map(|r| {
                r.iter()
                    .map(|p| [wrap_longitude(p[0] + offset), p[1]])
                    .collect()
            })
            .collect();
        if !crosses_antimeridian(&rotated) {
            return offset;
        }
    }
    0.0
}
