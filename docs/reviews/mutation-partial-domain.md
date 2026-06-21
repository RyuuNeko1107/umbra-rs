# mutation レビュー: 部分食域 partial_limit（`build_partial_limit` / `initial_bearing`）

対象: ISSUE-045 M9 残(3) サブスライス 3c-ii（部分食域 `EclipsePath::partial_limit` の外環組立・**リボン法**）。
`crates/umbra-eclipse/src/engine.rs` の `build_partial_limit`（南北半影限界を `北(P1→P4)++南(P4→P1 逆順)` で帯状
単純多角形に）と、`crates/umbra-eclipse/src/axis_intercept.rs` の `initial_bearing`（大圏初期方位・geo ユーティリティ）。

## 実行
```
cargo mutants -p umbra-eclipse --re 'build_partial_limit|initial_bearing' --no-shuffle \
  --exclude-re 'replace < with <= in.*build_partial_limit' \
  -- --test path_limits --lib -- \
     partial_phase partial_limit_none partial_only partial_limit_vertices partial_limit_ring partial_limit_contains bearing
```
killer は path_limits.rs の FAST partial 統合テスト（存在 Some/None・頂点が半影縁条件・リボン位相〔前半北/後半南逆順〕・
平面 point-in-polygon で partial⊃umbral・退化 interval=0 で None）＋ axis_intercept.rs の initial_bearing 単体（北/東/南/西の既知方位）。

## 結果（2026-06-22）
**build_partial_limit: 4 mutants（`<=` 除外後）・3 caught・1 unviable・0 missed。initial_bearing: 全 caught。**

## 生存 → 是正/許容の経緯
初回走（除外・新テスト前）で `build_partial_limit` の退化ガード `if ring.len() < 3` に 2 件生存:
| 変異 | 区分 | 対応 |
|---|---|---|
| `< → ==`（`ring.len() == 3`） | 撃破 | production では到達しない（`PathOptions` 既定 interval 60s・[P1,P4] は数時間＝多数サンプル）が、**`interval=0` で 1 サンプル＝外環 2 頂点**になり `<3`→None / `==3`→Some と分岐する。`partial_limit_none_when_single_sample_degenerate`（interval=0 → partial_limit None）を追加して撃破。 |
| `< → <=`（`ring.len() <= 3`） | 等価（許容） | 外環 = `北 n 点 ++ 南 n 点`（lockstep ゆえ north.len()==south.len()）で **長さは常に偶数**＝3 になり得ない。よって全到達長で `<3` と `<=3` は同一分岐＝真の等価。CI で `--exclude-re 'replace < with <= in.*build_partial_limit'` 除外（`build_partial_limit` 内の `<` は当該 1 箇所のみ）。 |

`initial_bearing` は production 未使用（リボン法は方位ソートを使わない・(3c-iii) limb 精緻化の順序付けで消費予定）だが、
単体テスト（既知方位）で atan2 引数順・sinΔλ/cosφ 各項・正規化を全捕捉。

## 方位ソート是正の記録（設計の轍）
当初 (3c-ii) は全境界点を最大食点まわりの大圏方位ソートで外環化したが、**実 2024-04-08（太平洋〜欧州の巨大領域）で
star-shaped 仮定が破綻**し中心線北東端が外環の外に落ちた（SLOW `real_2024_eclipse_partial_limit_is_plausible` が検出）。
位相保存のリボン法（限界線の時系列順を保つ帯）へ是正し、検証も radial（star-shaped 前提）から平面 point-in-polygon へ。
mutation だけでなく**実日食 SLOW オラクルが設計欠陥を捕捉**した好例。

実回帰ガード: 通常 CI の `cargo test -p umbra-eclipse`（FAST partial＋実 2024 SLOW）がリボン位相・包含・退化ガードを縛る。
