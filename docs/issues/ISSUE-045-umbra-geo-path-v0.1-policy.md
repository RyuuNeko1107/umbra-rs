# ISSUE-045: umbra-geo / path の v0.1 方針（型・境界のみ定義／本実装は Milestone 9）

- crate: umbra-geo
- 依存: ISSUE-022（`BesselianPolynomial`＝経路計算の供給源・公開型）, ISSUE-023（`GlobalCircumstances`・種別・最大食地点）, ISSUE-021（瞬時ベッセル要素 x,y,d,μ,l1,l2,tan f）, ISSUE-002（角度・緯度経度 newtype）, ISSUE-010（WGS84・測地/地心緯度）, ISSUE-044（`EclipseError::NotImplemented` 等のエラー集約）, ISSUE-001（規約）
- milestone: M9 経路（v0.1 では型・境界のみ。本実装は Milestone 9）
- モード(tdd-workflow): standard（v0.1 では公開型シグネチャと未実装スタブの契約のみを固定する。`SolarEclipse.bessel`（`BesselianPolynomial`＝ISSUE-022）は v0.1 で必須だが、`path()`/中心線/限界線/GeoJSON の数式本体は Milestone 9 のため、型境界の前方互換性確保が要点で standard。本実装時に strict へ昇格）

## M9 実装状況
- **M9.1 中心線トラック**（2026-06-21・strict）: `EclipseEngine::path()` を実装（旧 `Err(NotImplemented)` スタブから昇格）。中心食（全球 U1/U4 接触が両方 Some）で `center_line` を生成＝`[U1,U4]` を `PathOptions::sample_interval_seconds` 刻みでサンプルし、各時刻のベッセル要素（`BesselianSource::at`）から影軸地表貫通点（`axis_intercept::shadow_axis_surface_point`・WGS84）を結んだ `GeoLine`。軸が地球を外す端（`RootNotBracketed`）はスキップ。非中心は `center_line=None`。`greatest_point` は `global.greatest.position` passthrough。**北/南限界線・部分食域・`samples`（帯幅/継続）・GeoJSON は未実装**（後続スライス）。5 テスト（FAST 4＋SLOW 1: 実 2017-08-21 皆既で中心線が太平洋〜大西洋を延び北米を横断・最大食点近傍を通ることを実証）。mutation 12 中 8 caught・2 unviable・2 timeout（ループ終端変異＝ハング検出）・生存0。
- **M9.2 GeoJSON 出力**（2026-06-21・strict）: `GeoPoint::geojson_geometry()`（Point・[経度,緯度]順）/ `GeoLine::geojson_geometry()`（LineString・日付変更線 |Δlon|>180 で MultiLineString 分割）を umbra-geo に、`EclipsePath::to_geojson() -> Result<String, serde_json::Error>`（FeatureCollection: greatest_point の Point＋center_line〔Some 時〕の折れ線・pretty＋末尾改行）を umbra-eclipse に実装。両 crate に `serde_json` を production 依存追加。15 テスト（umbra-geo 9＋umbra-eclipse 6）/ mutation 12 中 12 caught・生存0。北/南限界線・部分食域・samples は未出力（後続）。交点の ±180 補間は後続改良。
- **M9.3 北/南限界線**（2026-06-21・strict・幾何近似）: 中心食で `include_limits` 時、各サンプルで ζ補正本影半径 `|L2'| = |l2 − ζ₀·tan f2|`（ζ₀＝影軸地表交点の基本面 ζ）を影の運動方向 (x′,y′) に**垂直**へ ±オフセットした 2 基本面点を地表へ射影し、高緯度側＝北限・低緯度側＝南限とする `northern_limit`/`southern_limit` を生成。中心線と**同一サンプル列**＝軸/縁が地表を外す（`RootNotBracketed`）/ 影速度ゼロのサンプルは 3 本ともスキップ（lockstep）。`axis_intercept::surface_point_for_fundamental` を新設＝任意の基本面 (ξ,η) から地表点＋ζ を返す（中心線 `shadow_axis_surface_point` も本関数へ委譲）。**幾何近似**（経路ほぼ東西走行・垂直オフセット前提。厳密な本影錐∩地表楕円体の接線解 ExplSup §11.3.5 は後続）である旨を accuracy.md §4.2 に明記（ユーザー指示「近似明記」）。テスト: path_limits 6（FAST 5＋SLOW 1: 実 2024-04-08 皆既の帯幅が NASA 公表 ~197 km と整合する [100,350] km）＋オラクル単体 4（`surface_point_for_fundamental` の二次閉形式オラクル / `sample_central_point` の独立オラクルで本影半径・時間尺・垂直法線・縁オフセットの各演算子を縛る FAST）＋path_center_line/path_geojson 回帰。mutation: 2 核（`surface_point_for_fundamental`/`sample_central_point`）65 中 **60 caught・5 unviable・生存0**（FAST オラクル単体で全算術変異を捕捉）、`trace_central` は path 統合テストで caught＋ループ終端 2 timeout（無限ループ＝ハング検出）・生存0。**併せて M9.1 で陳腐化した stub 契約 lib テスト 2 本**（`path` が `NotImplemented` を返す前提）を実装契約（非中心食は `Ok` ＋ `center_line=None`、greatest passthrough）へ修正＝M9.1〜M9.2 の検証が統合テストのみで `--lib` 未実行だったため見落とし、今回 `--lib` mutation ゲートで捕捉・修正。
- **M9.4 北/南限界線の厳密錐接線解**（2026-06-21・strict・M9.3 幾何近似を置換）: 中心食の南北限界線を**本影錐∩WGS84 楕円体の厳密接線解（路限界＝移動本影の包絡）**で生成（`engine::solve_limit_edge`・不動点反復）。各縁点 P が 2 条件を同時に満たす: (1) 錐exact＝影軸からの基本面内距離 = ζ補正本影半径 `|l2−ζ·tan f2|`（ζ は**P 自身**の値・M9.3 は中心軸 ζ₀ を流用していた誤差を是正）、(2) 包絡＝オフセットが影の**地表に対する相対速度** rel に直交（`rel=(x′−μ′(ζcosd−ηsind), y′−μ′ξsind)`＝影軸運動−自転運搬 ω×P。M9.3 は地球自転 μ′≈0.26rad/h を欠いていた）。WGS84 扁平は `surface_point_for_fundamental` の楕円体 root-find で厳密に処理（Almanac の ρ1/d1 近似不要）。d′ は無視（明記）。方式は第一原理（剛体回転＋錐∩楕円体）から導出、概念出典 ExplSup §11 / NASA Espenak（**式番号は一次資料未確認のため転記せず**）、数値オラクルは NASA 2024-04-08 公開 path table。テスト: path_limits SLOW の帯幅域を NASA 197.5km の **[185,215]km**（M9.3 の緩い [100,350] から狭帯化）に締め＋実日食で 2 条件（時刻復元でサンプル skip にロバスト）、FAST 合成（μ′≠0）で前方射影による 2 条件を機械精度（cone 1e-7・dot 1e-9）検証＋南北割当・include_limits=false・非中心 None。accuracy.md §4.2 を厳密版へ更新。
- **M9.5 限界線 GeoJSON 化＋日付変更線 ±180 補間**（2026-06-21・strict）: (A) `GeoLine::geojson_geometry`（umbra-geo）に日付変更線交点の **±180 線形補間**を実装＝跨ぎ点で交点緯度を子午線上に補間し前/次セグメント端へ補う（東進 Δlon<−180: 末尾 +180／先頭 −180・`t=(180−lon1)/(360+Δlon)`、西進 Δlon>180: 逆・`t=(lon1+180)/(360−Δlon)`、`lat_c=lat1+t·(lat2−lat1)`。RFC 7946 §3.1.9・M9.2 の「隙間が残る」分割を改良）。(B) `EclipsePath::to_geojson`（umbra-eclipse）に `northern_limit`/`southern_limit` の Feature を決定的順序（greatest→center_line→northern_limit→southern_limit）で追加（`role` プロパティ・Some 時のみ）。**`partial_limit` の GeoJSON 化は対象外**（常に None・GeoPolygon の GeoJSON 化は (3) と同時）。テスト: umbra-geo geojson 12（東進/西進/非対称t/二重跨ぎ/ちょうど±180 両側）＋umbra-eclipse path_geojson 8（4 feature・順序・geometry 一致・補間端点）。mutation `geojson_geometry` 40 中 **40 caught・生存0**（西進境界 `>→>=` の生存をレビュー指摘→西進ちょうど180テスト追加で撃破）。
- **M9.6 最大食点の帯幅・中心線継続**（2026-06-21・strict）: 中心食で `GreatestEclipse.path_width`/`central_duration` を Some に（従来 None）。**帯幅** = 相対速度包絡の南北本影縁点（M9.4 `solve_limit_edge` を t_max で適用）の大圏距離。**中心線継続** = `2·|L2'| / |rel|`（本影直径÷影の地表相対速度・初等運動学）× 3600 秒。`solve_limit_edge` と新 `great_circle_distance_km` を **engine.rs から axis_intercept.rs へ移設・共有**（global/engine の二層が消費）。供給源は直接評価ゆえ x'/y'/μ' を**数値中心差分**（±0.1h・per hour、μ' は ±2π 折返し補正 `wrap_to_pi`）。中心食のみ Some・部分/非中心は None（`shadow_axis_surface_point` の Ok/Err 分岐）。**算法 §8.11/8.12 の「要一次資料確認」は解消**＝帯幅は限界線の定義量・継続は初等運動学（直径÷速度）で第一原理導出、数値オラクルは NASA 2024-04-08（width 197.5km・duration 268.1s=4m28.1s）。テスト: global.rs 合成 source で `2|L2'|/|rel|` 厳密一致＋2017/2023 ballpark＋部分/非中心 None、path_limits SLOW で実 2024 が width∈[185,215]km・duration∈[250,286]s。**`samples`（PathSample 列）は後続**（per-sample width/duration/sun_altitude/kind）。
- **M9.7 経路サンプル列 samples**（2026-06-21・strict）: 中心食で `EclipsePath.samples`（`Vec<PathSample>`）を充足（従来空）。`include_limits=true` のとき各サンプルで中心線・北限・南限と**完全に同一サンプル列（lockstep）**＝`samples.len()==center_line.len()==northern_limit.len()==southern_limit.len()`・`samples[i].center==center_line.points[i]`。各 `PathSample`: `time_utc=tt_to_utc(t)`、`center`＝影軸地表点、`duration_seconds=2|L2'|/|rel|×3600`（中心軸 ζ₀・rel に μ' 項込み、M9.6 と同定義）、`sun_altitude`＝center 点の太陽幾何高度（`RefractionModel::None`）、`path_width`＝南北本影縁点間 大圏距離、`kind`＝`L2'=l2−ζ₀·tan f2` 符号別（`<0`→Total／それ以外 Annular。l2<0=皆既規約・hybrid 経路で切替）。`sample_central_point` を 4 本 lockstep（軸/縁/相対速度/継続定義不能で同期 skip）に拡張、`trace_central`/`path` は `delta_t` を太陽高度へ供給。`include_limits=false`・非中心では samples 空。`tt_to_utc` 失敗は伝播（lockstep を壊さず Err 終了）。テスト: path_limits FAST 5（lockstep・独立オラクル〔duration の 2|L2'|/|rel|・width の大圏・kind 符号〕・time_utc=tt_to_utc＆単調・include_limits=false 空・非中心空）＋SLOW 1（実 2024 で lockstep・全 Total・width∈[185,215]・duration∈[250,286]・UTC 単調）。M9.1 で陳腐化した path_center_line の「samples 空」契約 2 本を lockstep 充足契約へ更新。mutation: `docs/reviews/mutation-path-samples.md`。
- **M9残(3) 部分食域**（`docs/algorithms/11-path-partial-domain.md`・第一原理＋NASA 2024-04-08 オラクル）。設計ドラフト→厳密化済み（要確認 #1 rise/set は外周に 2 根曲線の両方使用・4 分類は出力ラベル／#2 terminator は WGS84 楕円 `ξ²+k·η²=1`（k=sin²d+cos²d/(1−f)²）で厳密化／#5 solve_limit_edge 一般化の M9.4 mutation 影響なしを確認、で解決）。サブスライス:
  - **(3a) 半影限界の半径引数化** ✅（2026-06-21・strict）: `solve_limit_edge` に錐半径 `(cone_l,cone_tan_f)` を引数化＝本影 `(l2,tan f2)`／半影 `(l1,tan f1)` 両対応。本影 M9.4 退化で完全回帰。mutation 52/51 caught・0 missed。併せて mutation.yml の未エスケープ `||`（空 alternation で全除外＝CI mutation no-op バグ）を `\|\|` へ修正。
  - **(3b) rise/set limb 点** ✅（2026-06-21・strict）: `cone_terminator_intersections`＝錐縁 ∩ WGS84 terminator 楕円（θ 媒介の円残差を粗走査＋Brent、機構は `scan_periodic_sign_change_roots` に分離）→ `fundamental_to_geodetic(ζ=0)`。WGS84 前方射影往復（ζ≈0・面内距離=cone_l）＋d=π/2 二円閉形式で検証。mutation 40/39 caught・0 missed（機構は wholesale 除外）。API 露出なし。
  - **(3c-i) 南北半影限界の曲線化** ✅（2026-06-21・strict）: `trace_penumbral_limits`＝[P1,P4] を sample し `solve_limit_edge(l1,tan f1)` で南北半影限界の lockstep 2 GeoLine。`trace_central` と同型・ζ₀=0（部分食で軸が地球を外しうる）。FAST 4（2 条件・lockstep・半影>本影距離・ループ規約）。mutation 14/11 caught・2 timeout（ループ制御＝ハング検出）・0 missed。API 露出なし。
  - **(3c-ii) 外環組立（リボン法）** ✅（2026-06-22・strict）: `build_partial_limit`＝南北半影限界（lockstep）を `北(P1→P4)++南(P4→P1 逆順)` で帯状単純多角形にし `path()` の `partial_limit` を Some に（部分食 phase〔P1/P4 両 Some〕＋include_limits 時・中心食と独立）。**初の観測可能な部分食域 API 成果物**。検証: 平面 point-in-polygon で partial⊃umbral・リボン位相・退化 interval=0 で None・実 2024。mutation 4/3 caught・0 missed（`<=` は ring 長偶数ゆえ等価除外）。**当初の方位ソート設計は実 2024〔太平洋〜欧州の巨大領域〕で star-shaped 破綻→SLOW オラクルが捕捉→リボン法へ是正**。`initial_bearing`（geo ユーティリティ）追加。docs/reviews/mutation-partial-domain.md。
  - **(3d) GeoPolygon GeoJSON** ✅（2026-06-22・strict）: `GeoPolygon::geojson_geometry`（RFC 7946 Polygon・閉リング・環向き正規化〔外環 CCW/穴 CW・`signed_area_lonlat`〕・退行非捏造・**v1 反子午線非分割**）＋ `to_geojson` に `partial_limit` feature（`role="partial_limit"`・southern_limit の後）。これで `to_geojson` が EclipsePath の全主要要素（greatest/center_line/北南限界線/partial_limit）を出力。mutation 64/64＋2/2 caught・0 missed（docs/reviews/mutation-geojson-polygon.md）。
  - **(3c-iii) limb bulge 精緻化（端区間 terminator 連結・genuine bulge）** ✅（2026-06-22・strict）: `trace_penumbral_limits` の各サンプルで、昼面包絡 `solve_limit_edge(l1,tan f1)` が**解けない**端区間を terminator 交点（`cone_terminator_intersections`・円∩terminator 楕円・ζ=0・半影半径 l1）の高/低緯度側（新ヘルパ `highest_latitude`/`lowest_latitude`）で連結（最終対は緯度ソートで北≥南・lockstep 維持・両昼面包絡が解ければ terminator 非計算）。外環が昼面包絡の欠ける端区間で terminator まで張り出す（v1 より limb 方向に広い）。公開IF（`path()`/`partial_limit`）・リボン構成不変。terminator 頂点も半影縁条件（ζ=0 で面内距離=l1）を満たし「頂点正当性」維持。**到達範囲の確定**: 端区間連結は v1 より広いが **full center-line containment は未達**（実 2024 早期端 ~6.7°S sunrise が帯西外＝真の west/east 境界は morning/evening terminator limb の [P1,P4] 全域追跡＝4 曲線境界で、設計が v1 で「脆い」と回避・**後続へ繰延**）。ユーザー判断で本スライスは genuine bulge まで・受入を「bulge 発火・帯が皆既帯より広い・最大食付近の中心線内包」へ較正。検証: FAST 3（`limb_continuation_bessel` で端 terminator 連結発火・ζ≈0 頂点・半影縁条件・リボン不変）＋ヘルパ単体 2（極値選択を緯度ソートに masking されず固定）＋SLOW 実 2024（partial=Some・帯>皆既帯・実データでも ζ≈0 terminator 頂点・最大食付近包含）。mutation 19 中 13 caught・2 timeout（ループ＝ハング）・3 unviable・1 survivor（混在分岐 `&&→\|\|`＝loose 近似・実 2024 含む全受入で不可分ゆえ除外＝過仕様回避）・実質0 missed（docs/reviews/mutation-limb-bulge.md）。`initial_bearing`（(3c-ii) で「(3c-iii) で消費予定」とした geo ユーティリティ）は本実装が緯度ベースのため**未消費**＝test-only のまま（後続 4 曲線境界での消費は未定）。
  - **(3f) full containment＝3 領域の多角形ユニオン** ✅（2026-09-14・strict）: 要確認 6 の実測で「4 曲線の連結では外周が出ない」ことが確定（terminator limb の 2 枝が帯と重なる独立ループ＝継ぎ目の長い弦が中心線点を外に落とす）。人間決定で**多角形ユニオン方式**を採用。`umbra_geo::union_rings`（新規 `clip.rs`・平面アレンジメント方式＝辺分割→量子化 1e-9°→無向重複除去→両側 membership 境界判定→角度優先の環再結合→穴割当。**自己交差入力を許容**・even-odd 充填）を追加し、`build_partial_limit` を **帯 ∪ morning lune ∪ evening lune** の 3 領域ユニオン（最大面積成分）に置換。morning/evening は基本面 ξ の符号（`dζ/dt=−μ′cos d·ξ`）で厳密分類（`cone_terminator_intersections_detailed`）。lune は lockstep ribbon（交点無しの時刻は両枝同時に落とす・接点は hi=lo）。**中心線全点内包を headline acceptance に昇格**＝実 2024-04-08 で 194 点全て内包（従来 190/194）。リボン位相契約（偶数頂点・同時刻対）は撤廃（§11.6(e)）。検証: umbra-geo `tests/union.rs` FAST 29（合成多角形の厳密頂点列・面積・環数・向き・三角形/穴/点接触/共線重なり/自己交差/閉開入力不変）＋path_limits 27（合成 bessel で外環頂点が閉半影内かつ昼面側・中心線内包、SLOW 実 2024 全点内包）。mutation（clip.rs）: 275 変異・234 caught・5 timeout（ループ＝ハング検出）・1 unviable・**35 missed＝31 等価＋3 系統未解決**（分岐頂点の後継選択・二重入れ子の穴割当・共線部分重なりの分割点＝現行オラクルが到達しない。要確認 7 として登録）。**追試（同日・strict）**: 判別テスト 11 本追加（union 40 本）で穴割当・共線重なりを撃破、**実装欠陥を発見・修正**（穴同士／穴と外環の 1 点接触が自己接触環に融合→pinch 分割 `split_at_repeated_vertices`）、`param_on_segment` を主軸パラメータ化。最終 268 変異・236 caught・26 missed＝21 等価＋後継選択 5 件（7 配置で不可分・許容・要確認 7 縮小）。docs/reviews/mutation-polygon-union.md。umbra-eclipse 側（lunes/detailed/build_partial_limit）: 58 変異・39 caught・3 timeout（ループ＝ハング）・7 unviable・9 missed＝5 は detailed の楕円算術で `--lib` 単体が撃つ（40/36 caught・0 missed）＋4 は分類式 `(ξ<0)==is_morning` のラベル入替/測度ゼロ/受入不感（ユニオンが対称ゆえ公開 IF 不変・mutation.yml 除外）。docs/reviews/mutation-terminator-lunes.md。近似レジスタ（conventions §12 新設・accuracy §4.3）: 平面 (lon,lat) ユニオン＝**反子午線/極 未対応**・複数成分時は最大面積成分のみ（要確認 8）・同一 limb 交点 3 点以上で中間点を捨てる・1e-9° 量子化。mutation.yml の `build_partial_limit` 除外（ring 長偶数前提）は構造変更で撤去。
  - **(3g) 反子午線 MultiPolygon 分割** ✅（2026-09-14・strict・§11.8）: `GeoPolygon::geojson_geometry` が ±180 を跨ぐリングを **MultiPolygon へ分割**するように（跨ぎ判定・緯度の線形補間は `GeoLine` (3d) と**同一規約**＝`|Δlon|=180` ちょうどは跨ぎとしない）。**適用層は出力層のみ**＝`EclipsePath::partial_limit`（`Option<GeoPolygon>`）の型・値は不変（公開型拡張を避ける・要確認 8 の決定と整合）。
    実装（`umbra-geo::clip`）: `split_ring_at_antimeridian`（跨ぎで弧に切り、子午線上で結んで閉じる。結線の向きは**リング自身の実際の向き**が内部を左に保つ向き＝跨ぐリングの平面 shoelace は無意味なので `unwrap_longitudes` で経度を連続化してから判定）＋ `close_arcs_along_meridian`（緯度順の後継表で巡回を取り出す）＋ `split_polygon_at_antimeridian`（**分割→環向き正規化→穴の再割当**の順。穴は含む外環断片のうち最小面積へ・判定点は境界に載らない**頂点重心**）。跨ぎが無ければ分割経路に入らず**出力はバイト不変**。
    **極を囲む領域は未対応**を厳密化: 閉リングの跨ぎ回数は経度の巻き数に等しく**奇数回＝極を囲む**ので、この場合は**分割せず元のリングを返す**（弧が対を成さない。面積を失わず捏造もしない）。実測で 2021-12-04（1 回）・2028-07-22／2037-07-13（3 回）が該当し、2016-03-09（2 回）は通常の跨ぎ。
    実装レビュー（別エージェント）が**テスト未到達の欠陥 2 件**を指摘・修正: ①弧の対応付けが部分失敗すると**面積を黙って失う**（自己交差入力）→ 全弧を使い切れなければ分割を諦めて元のリングを返す全か無かのガード、②穴の内包判定に**頂点重心**を使うと凹（C 字）外環断片の切り欠きに落ちて穴が消える→**子午線外の穴頂点**を判定点に（`vertex_centroid` は死コードとなり削除）。
    mutation（clip.rs＋geometry.rs 496 変異）: 445 caught・7 timeout・2 unviable・42 missed。うち (3g) 由来は **32→15 件**（判別テスト 17 件撃破）で、残りは測度ゼロ境界（`|Δlon|=180` ちょうど・面積厳密ゼロ）・同値 tie-break・到達不能な防御ガードのみ。docs/reviews/mutation-antimeridian-split.md。
    検証: geojson FAST 38（+18・跨がない場合のバイト不変／Δlon=180 の測度ゼロ境界／東西両方向の分割と厳密座標・補間緯度／穴も跨ぐ場合の再割当／片半球のみの穴／4 回跨ぎ＝3 断片／極を囲む奇数跨ぎ 2 本／退行）＋ path_limits SLOW 実 2016-03-09（partial_limit の GeoJSON が MultiPolygon・全リング閉・経度域・外環 CCW・面積正）。
  - **(3h) 極を囲む領域** ✅（ISSUE-051・commit 67b3e59/90dbd8a）: 跨ぎ回数が奇数のとき「分割せず元のリングを返す」退避が経度幅 358.5° の不正多角形になっていた欠陥を、**リングを該当極で閉じる**方式で解消（エンジン側で ζ から極の昼面判定）。
  - **(3i) union の経度フレーム依存** ✅（ISSUE-052・commit 4a2d23c）: 反子午線を跨ぐ入力で平面 even-odd/pip/shoelace が壊れ、同一領域が約 17 倍の別図形になっていた欠陥を、**跨がない経度フレームへ回してから解く**方式で解消（極を囲む領域のみ回転で解消できず従来どおり）。
  - **(3j) 退化配置の穴割当** ✅（ISSUE-053）: 自己接触を含む配置で穴が外環より大きくなる構造契約の破れ＝
    領域の捏造を、穴の帰属判定を「穴リング全頂点の閉内部包含＋面積必要条件」へ厳密化して解消。
  - **残**: **要確認 7**（union 後継選択の反例未構成・許容済み）。ISSUE-053 のランダム化差分探索で
    **5 件→4 件へ縮小**（`back - atan2(…)` → `+` の判別配置 E を決定的テスト化）。修正前に見つかった
    判別配置は後継選択ではなく ISSUE-053 の穴割当欠陥の露出だった。要確認 8 は「v1 は最大面積成分のみ・MultiPolygon 化は反子午線分割と一括」で決定済み。**v0.1 完成条件に path/GeoJSON は非含**（search/local/next_visible が v0.1）ゆえ M9 path 群はここで一区切り可。

## 目的
`umbra-geo` の経路 API（中心線・限界線・部分食域・GeoJSON）の **公開型と境界のみを v0.1 で確定**し、**本実装を Milestone 9 へ明示的に後回し**する方針を文書化する（レビュー minor 確定事項 / milestone0-review §Minor「045 umbra-geo/path はv0.1スコープ外だが結果型が bessel多項式(022)必須 → v0.1は path未実装方針を明文化」）。
- v0.1 完成条件（search・種別・最大食時刻・C1/最大/C4・食分食面積・50地点誤差レポート）に **`path()` は含まれない**。一方、`SolarEclipse.bessel: BesselianPolynomial`（api-draft §3.4）は v0.1 でも必須フィールドであり、ISSUE-022 が供給する。
- 本 Issue では `umbra-geo` の公開型（`GeoPoint`/`GeoLine`/`GeoPolygon`/`EclipsePath`/`PathSample`/`PathOptions`、api-draft §4）を **型として定義**し、`EclipseEngine::path()` は **v0.1 では未実装スタブ＝`Err(EclipseError::NotImplemented)`**（PATH 確定）とする。`EclipsePath::to_geojson()` も v0.1 未実装（Milestone 9）。
- 型と境界（フレーム規約・単位・日付変更線/極域の扱い方針）だけを固定し、中心線/限界線/GeoJSON の**数式本体は Milestone 9** であることを明文化する。

## 非目的
- 中心線・北限/南限・部分食域・GeoJSON の**実計算**（Milestone 9。本 Issue はスタブと型のみ）。
- 経路サンプリングの数式・日付変更線分割・極域特異点処理の実装（Milestone 9、accuracy.md / algorithms.md で別途）。
- `BesselianPolynomial`（ISSUE-022）・全球分類（ISSUE-023）の実装。本 Issue はそれらを**消費する境界**を置くのみ。
- v0.1 CLI の `path` サブコマンド本体（umbra-cli。スタブ呼出しで「未実装」を明示する整形のみ許容）。

## 公開インターフェース
api-draft §4 をそのまま型として確定（実装は Milestone 9）。
```rust
#[derive(Clone, Copy, Debug)] pub struct GeoPoint { pub lat: GeodeticLatitude, pub lon: EastLongitude }
#[derive(Clone, Debug)] pub struct GeoLine { pub points: Vec<GeoPoint> }
#[derive(Clone, Debug)] pub struct GeoPolygon { pub rings: Vec<Vec<GeoPoint>> }

#[derive(Clone, Debug)]
pub struct EclipsePath {
    pub center_line: Option<GeoLine>,
    pub northern_limit: Option<GeoLine>,
    pub southern_limit: Option<GeoLine>,
    pub partial_limit: Option<GeoPolygon>,
    pub greatest_point: GeoPoint,
    pub samples: Vec<PathSample>,
}
#[derive(Clone, Copy, Debug)]
pub struct PathSample {
    pub time_utc: UtcInstant, pub center: GeoPoint,
    pub duration_seconds: f64, pub sun_altitude: Degrees,
    pub path_width: Kilometers, pub kind: SolarEclipseKind,
}
#[derive(Clone, Copy, Debug)]
pub struct PathOptions { pub sample_interval_seconds: f64, pub include_limits: bool, pub split_antimeridian: bool }

impl EclipsePath {
    /// v0.1 未実装。Milestone 9 で実装。
    #[cfg(feature = "geojson")] pub fn to_geojson(&self) -> String;   // v0.1: 未実装（戻り型が String のため呼出経路に乗せない。CLI は「未実装」整形表示）
}
```
- `EclipseEngine::path(&self, eclipse: &SolarEclipse, options: PathOptions) -> Result<EclipsePath, EclipseError>`（api-draft §3.2）は v0.1 では**未実装スタブ**。
- **v0.1 スタブの戻り方（統一規則・確定 PATH）**: 「対応年代外」ではなく「機能未提供」を表すため、**`Err(EclipseError::NotImplemented)` を返す**（panic/`unimplemented!` は採用しない）。`UnsupportedTimeRange` は「対応年代外」専用語義に保ち、**未実装の意味に流用しない**。`NotImplemented` variant は ISSUE-044 で追加。CLI など実行経路に乗っても `Err(NotImplemented)` を「Milestone 9 で対応予定」と整形表示し、空 `EclipsePath` を成功として返さない。

## 数式・アルゴリズムの出典
- 本 Issue は**型・境界の確定のみ**で数式を持たない（数式本体は Milestone 9）。
- 参照（Milestone 9 で使う出典の予約・本 Issue では実装しない）:
  - 中心線・限界線・部分食域: ベッセル要素からの地上投影（Explanatory Supplement to the Astronomical Almanac、Espenak/NASA の経路生成手順）。**要確認**（一次資料の式番号は Milestone 9 で確定）。
  - GeoJSON: RFC 7946（日付変更線をまたぐ線分の分割規約・§3.1.9）。**要確認**。

## 単位 / 時刻系 / 座標系
- 角度: 公開は度（`GeodeticLatitude`/`EastLongitude`、conventions §3）。経度は**東経正** `[-180°,180°)`。
- 時刻: `PathSample.time_utc` は UTC（accuracy.md §0。TT 併記が必要なら Milestone 9 で `PathSample` を拡張、`#[non_exhaustive]` 検討）。
- 座標系: 地上点は ITRS→測地座標（WGS84、conventions §4/§5）。フレーム連鎖は ISSUE-035（GCRS→CIRS→TIRS→ITRS）に従う。
- 距離: 食帯幅は km（`Kilometers`、conventions §1）。継続時間は秒。

## アルゴリズム概要
v0.1（本 Issue のスコープ）:
1. api-draft §4 の公開型を `umbra-geo` に定義（フィールド・単位・フレーム規約を確定）。
2. `EclipseEngine::path()` を**未実装スタブ＝`Err(EclipseError::NotImplemented)`**として置く（前項「戻り方」規則、PATH）。`EclipsePath::to_geojson()` も v0.1 未実装。
3. `SolarEclipse.bessel`（`BesselianPolynomial`、ISSUE-022）は v0.1 で必須のため、**型として参照可能**にする（umbra-eclipse 側で生成。本 Issue は経路側で消費する境界の型整合のみ）。
4. ドキュメント（本 Issue・README・api-draft §4 注記）に「path は Milestone 9」と明記。

Milestone 9（本 Issue の非目的・予約）: ベッセル多項式から中心線/限界線/部分食域をサンプリングし、`EclipsePath` を構築。日付変更線分割・極域特異点処理・GeoJSON 出力。

## 受け入れテスト
v0.1（本 Issue）:
- 型整合: `EclipsePath`/`PathSample`/`PathOptions`/`GeoPoint`/`GeoLine`/`GeoPolygon` が api-draft §4 のフィールド・単位で定義され、`SolarEclipse.bessel: BesselianPolynomial` を含む `SolarEclipse` がコンパイル可能（型レベル検証）。
- スタブ契約: `path()` 呼出しが「未実装」を表す＝`Err(EclipseError::NotImplemented)` を返す（panic でなく Result、`UnsupportedTimeRange` を流用しない）。`assert!(matches!(.., Err(EclipseError::NotImplemented)))` で固定。**v0.1 の通常経路（search/local/next_visible）から path が呼ばれないこと**もテストで保証。
- CLI 整合: `umbra path`（あれば）は「Milestone 9 で対応予定」を表示し、誤った経路（空 `EclipsePath` を成功として返す等）を作らない。
- 前方互換: 列挙・設定型は `#[non_exhaustive]`/`Default` で Milestone 9 拡張時に破壊的変更を避けられる（api-draft §0）。
- 二段オラクルゲート（ISSUE-047 連動）: **本 Issue は v0.1 でスタブのため数値ゲート対象外**。Milestone 9 実装時に「M2 暫定ゲート（Mock+SOFA+NASA 経路値）」と「M10 最終ゲート（DE 差分）」を付す（ISSUE-047 の二段方針を継承）。

## 許容誤差
- v0.1（本 Issue）: 数値計算を行わないため**許容誤差なし**。
- Milestone 9（予約・本 Issue では保証しない）: 中心線位置 sub-km（≲0.5 km、幾何分、accuracy.md §1 Standard）。fit 残差は `BesselianPolynomial.fit_error`（ISSUE-022）でガード。

## 実装メモ
- 本 Issue は **milestone0-review §Minor の確定事項**の反映: 「v0.1 は path 未実装方針を明文化」。型と境界だけ定義し、本実装は Milestone 9 へ後回しを明示する。
- `SolarEclipse.bessel`（`BesselianPolynomial`）は **v0.1 で必須**（api-draft §3.4）。これは ISSUE-022 が供給し、本 Issue（umbra-geo）はそれを消費する経路側の型のみを持つ。両者の責務を混同しない。
- スタブの戻り方は ISSUE-044（`EclipseError` 集約）と整合し、**`Err(EclipseError::NotImplemented)` に統一**（確定 PATH）。`UnsupportedTimeRange` は「対応年代外」専用語義に保ち未実装に流用しない。`unimplemented!`（panic）は撤回。
- `umbra-geo` は v0.1 では実質スケルトン。ただし公開型は SemVer 境界なので、Milestone 9 で破壊しないよう `#[non_exhaustive]` とフィールド追加余地を意識する。
- レビュー重点: 「v0.1 で path を呼ばせない」保証、型の前方互換、`bessel` 必須と path 未実装の責務分離、スタブ語義の一貫性（**`Err(EclipseError::NotImplemented)` に統一**、panic/`unimplemented!` と `UnsupportedTimeRange` 流用を排除、PATH）。
