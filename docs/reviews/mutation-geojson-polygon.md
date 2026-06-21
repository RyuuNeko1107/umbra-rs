# mutation レビュー: 部分食域 GeoJSON（`GeoPolygon::geojson_geometry` / `signed_area_lonlat` / `to_geojson`）

対象: ISSUE-045 M9 残(3) サブスライス 3d（部分食域の GeoJSON 化）。`crates/umbra-geo/src/geometry.rs` の
`GeoPolygon::geojson_geometry`（RFC 7946 Polygon・閉リング・環向き正規化）と `signed_area_lonlat`（shoelace 向き判定）、
`crates/umbra-eclipse/src/path.rs` の `EclipsePath::to_geojson` への `partial_limit` feature 追加。

## 実行
```
cargo mutants -p umbra-geo --re 'geojson_geometry|signed_area_lonlat' --no-shuffle -- --test geojson
cargo mutants -p umbra-eclipse --re 'to_geojson' --no-shuffle -- --test path_geojson
```
killer は umbra-geo geojson 統合テスト（Polygon 構造・[lon,lat]順・閉リング・環向き CCW/CW・退行・反子午線単一 Polygon）＋
umbra-eclipse path_geojson（partial_limit feature の有無・順序・geometry 一致）。

## 結果（2026-06-22）
- **umbra-geo `geojson_geometry`/`signed_area_lonlat`: 64 mutants・64 caught・0 missed。**
- **umbra-eclipse `to_geojson`: 2 mutants・2 caught・0 missed。**
除外（mutation.yml 追記）は不要。

## 初回 3 missed → 弁別テスト追加で 0 missed
初回走で環向き正規化の orientation ロジックに 3 件生存。いずれも**弁別テスト**で撃破（等価ではなく test-gap だった）:
| 変異 | 撃ち方（追加テスト） |
|---|---|
| `geometry.rs:164` 外環 `area < 0.0`→`<=` | `geo_polygon_geojson_zero_area_outer_not_reversed`: 共線（面積 0）外環で original は非反転・`<=` は反転＝出力経度順 `[0,2,1]` vs `[0,1,2]` で弁別。 |
| `geometry.rs:164` 穴 `area > 0.0`→`>=` | `geo_polygon_geojson_zero_area_hole_not_reversed`: 共線（面積 0）穴で同様に順序弁別。 |
| `geometry.rs:181` shoelace `x1*y2`→`x1+y2` | `geo_polygon_geojson_shoelace_product_sign_decides_winding`: `(0,0),(0,1),(1,0)` 閉リングで original 式=−1（CW→反転して CCW 出力）/ `+` 式=+1（非反転で CW のまま）＝出力面積符号・座標順で弁別。 |

面積 0 境界（`<`/`>`→`<=`/`>=`）も「退行リングを反転しない」契約を共線リングで明示的に縛ることで等価扱いを避け caught 化した
（degenerate でも決定的出力＝順序保存を保証）。

実回帰ガード: 通常 CI の `cargo test -p umbra-geo`/`-p umbra-eclipse`（Polygon 構造・環向き・退行・反子午線単一 Polygon・
partial_limit feature 順序）が GeoJSON 出力の正しさを縛る。

## v1 の限界（明記）
反子午線（±180 跨ぎ）の MultiPolygon 分割は未対応＝跨ぐリングは単一 Polygon（ポリゴンクリッピングは後続精緻化）。
テスト `geo_polygon_geojson_antimeridian_stays_single_polygon` がこの v1 仕様を固定。
