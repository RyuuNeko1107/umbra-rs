# mutation レビュー: 部分食域 limb bulge / terminator 連結（`trace_penumbral_limits` 3c-iii ＋ `highest_latitude`/`lowest_latitude`）

対象: ISSUE-045 M9 残(3) サブスライス 3c-iii（部分食域 `EclipsePath::partial_limit` の **limb bulge＝端区間 terminator 連結**）。
`crates/umbra-eclipse/src/engine.rs` の `trace_penumbral_limits`（昼面南北半影限界に、昼面包絡が欠ける端区間で
terminator 交点〔円∩terminator 楕円・ζ=0・半影半径 l1〕を連結する分岐を追加）と、新ヘルパ `highest_latitude`/
`lowest_latitude`（terminator 交点群の緯度極値選択・北縁/南縁連結）。

## 実行
```
cargo mutants -p umbra-eclipse --re 'trace_penumbral_limits|highest_latitude|lowest_latitude' --no-shuffle \
  -- --test path_limits --lib -- \
     penumbral latitude partial_limit_ring partial_limit_limb partial_phase partial_only partial_limit_none partial_limit_contains partial_limit_vertices
```
killer:
- **新ヘルパ単体**（engine.rs `highest_latitude_picks_maximum_latitude_point` / `lowest_latitude_picks_minimum_latitude_point`）
  ＝極値選択を `trace_penumbral_limits` の緯度ソート（北≥南 再整列）に masking されない純関数レベルで固定（`max_by↔min_by`・
  比較方向・最初/最後の要素を撃つ。非対称緯度 +40/−7/−30）。
- **FAST limb fixture**（path_limits.rs `limb_continuation_bessel` の `partial_limit_ring_includes_terminator_vertices_on_limb` /
  `partial_limit_limb_vertices_satisfy_penumbral_conditions` / `partial_limit_limb_ring_preserves_ribbon_invariants`）＝端区間で
  terminator 連結が発火し ζ≈0 頂点が現れる・半影縁条件・リボン不変条件。
- **(3c-i) `penumbral_*` 単体**（trace_penumbral_limits の昼面包絡・lockstep・ループ規約。回帰）。

## 結果（2026-06-22）
**19 mutants: 13 caught・3 unviable・2 timeout（ループ制御＝ハング検出）・1 survivor（除外＝下記）。実質 0 missed。**

| 区分 | 件数 | 備考 |
|---|---|---|
| caught | 13 | ヘルパ極値・terminator 連結発火・半影半径・lockstep の算術/分岐 |
| timeout（caught） | 2 | line 766 `\|\| → &&`（break 条件が成立せず無限ループ）／ line 769 `+ → *`（`t_sec*interval` が 0 固定で無限ループ）＝ハング検出 |
| unviable | 3 | 型不整合でコンパイル不能 |
| survivor（除外） | 1 | line 747 `&& → \|\|`（下記・mutation.yml 除外） |

## survivor → 除外の経緯
| 変異 | 区分 | 対応 |
|---|---|---|
| line 747 `north_day.is_some() && south_day.is_some()` の `&& → \|\|` | 除外（loose 近似分岐・過仕様回避） | この変異は**混在サンプル**（昼面包絡が片側のみ解ける）でのみ挙動が変わる: `&&`＝欠側を terminator で補いサンプルを保持／`\|\|`＝if 枝に入り片側 None のまま lockstep skip（drop）。混在分岐の正確な頂点は**意図的に loose な近似**（limb 整合を取らず off-limb partner になりうる＝実装レビュー Finding 1・§11.4 (3c-iii) で明記。full containment を実現する 4 曲線 terminator-limb 境界は後続へ繰延）。**SLOW 実 2024 を含む全受入オラクルで両挙動が不可分**＝手動で `&&→\|\|` を適用し `real_2024_eclipse_partial_limit_is_plausible`＋FAST limb 3 本が全 pass することを確認（実 2024 の bulge は both-None の clean cap サンプル由来で、混在サンプルは幾何的に希少/受入不感）。loose と宣言した内部挙動を test で固定すると過仕様になるため、`--exclude-re 'replace && with \|\| in.*trace_penumbral_limits'` で除外（関数内の `&&` は当該 1 箇所のみ）。 |

## 設計到達範囲の記録（scope 判断の轍）
3c-iii の当初狙いは「中心線全点が partial_limit に平面 point-in-polygon で内包」（full containment・§11.5）。しかし端区間
terminator 連結（**昼面包絡が無いサンプルだけ terminator で補う**）では、実 2024 の中心線早期端（~6.7°S・sunrise）が帯の西外
に落ち full containment に**未達**と判明（probe: 北限界が lon≈−144 から東進開始・中心線 lon≈−152 が帯西外）。真の west/east 境界＝
morning/evening rise/set **terminator limb を [P1,P4] 全域追跡する 4 曲線境界**が必要だが、設計が v1 で「脆い」として回避した
ステッチであり、ユーザー判断で **本スライスは genuine bulge までを成果物**とし full containment は後続スライスへ繰延（受入を
「bulge 発火・帯が皆既帯より広い・最大食付近の中心線内包」へ較正）。**SLOW オラクルが scope の限界を捕捉**した好例。

実回帰ガード: 通常 CI の `cargo test -p umbra-eclipse`（FAST limb＋実 2024 SLOW partial）が連結発火・リボン不変・包含核を縛る。
