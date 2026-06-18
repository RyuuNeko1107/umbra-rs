# Mutation review — `umbra-eclipse::engine` 局地条件（ISSUE-043 S7b-i）

`EclipseEngine::local_circumstances`（局地最大食の結線）と補助 `build_local_contact` の
`cargo mutants` 生存変異の列挙・許容判断、および退行ガード `mutation.yml` での扱い。

## `build_local_contact`（新規ヘルパ・ミューテーション実行済）

```
cargo mutants --package umbra-eclipse --jobs 2 --re 'build_local_contact' -- local_circumstances
→ 2 mutants tested: 1 caught, 1 unviable（0 missed）
```

- **caught**: `visible: altitude_geometric.0 >= 0.0` の比較・フィールド代入系。
  局地条件テスト（`local_circumstances_central_site_is_total` の `visible==true`、
  `local_circumstances_angles_in_valid_ranges` の alt/az/PA 値域）が撃破。
- **unviable**: 関数本体を `Default::default()` 等へ置換する変異。`LocalContact` は `Default`
  非実装のためコンパイル不能（＝到達不能・生存ではない）。

生存ゼロ。`build_local_contact` のフィールド配線（alt/az/PA/visible）は局地条件テストで縛られている。

## `local_circumstances` 本体はオーケストレーション

`local_circumstances()` は確定済みの下位プリミティブを配線して `LocalCircumstances` を組み立てる
工程であり、固有の数値ロジックは ζ 補正半径による食分/食面積のみ。各下位工程は**それ自身のテストで
独立にミューテーション検証済み**:

- `project_observer_to_fundamental`（ISSUE-024）/ `solve_local_maximum`（ISSUE-026,
  `mutation-local-maximum.md`）/ `sun_horizontal`・`classify_visibility`（ISSUE-028,
  `mutation-apparent.md` 系）/ `eclipse_magnitude`・`eclipse_obscuration`（ISSUE-027）/
  `contact_position_angle`（ISSUE-043 S7a, `position_angle.rs` ミューテーション 0 missed）。
- ζ 補正半径 `L1'=l1−ζ·tanf1 / L2'=l2−ζ·tanf2` と `radius_ratio`/`separation` の式は
  **`global.rs::solve_greatest_eclipse` と同一**で、そちらは S6a-ii で
  ミューテーション検証済み（`mutation-local-maximum.md` 系・16caught/1unviable/0missed）。

配線と分岐は統合テスト（合成 `standard_engine` で 2017-08-21 を全パイプライン実走）が縛る:

| テスト | 縛る対象 |
|---|---|
| `local_circumstances_central_site_is_total` | 中心食地点: magnitude>1・obscuration≈1・PartialVisible（S7b-i 中間）・観測者配線 |
| `local_circumstances_partial_site_is_partial` | 部分食地点: magnitude<1・中心食との差>0.05（observer lat/lon 配線） |
| `local_circumstances_invisible_site_is_not_visible` | 非可視（Ok 分岐・min_sep≥L1）: NotVisible・magnitude/obscuration 0・in_eclipse 判定 |
| `local_circumstances_unbracketable_window_anchors_at_global_greatest` | **錨分岐（RootNotBracketed）**: 退化窓 `[t,t]` で機構的に励起。全球最大食 TT/UTC 錨・magnitude/obscuration 0・NotVisible。**search 非依存で FAST** |
| `local_circumstances_angles_in_valid_ranges` | alt/az/PA 値域（build_local_contact 配線） |
| `local_circumstances_maximum_utc_tt_consistent` | maximum.time_utc == tt_to_utc(time_tt) |
| `local_circumstances_metadata_recipe` | metadata レシピ転記 |

これらのうち錨分岐テスト以外は内部で `engine.search`（≈300s）を呼ぶため各々が重い。

## `mutation.yml`（週次ガード）での扱い

`local_circumstances` 系テスト（錨分岐テストを除く 6 件）は内部で `search` を実走するため各々
≈300s。これらを毎ミューテーションの `cargo test` に含めると umbra-eclipse の**全**変異が
300s 超になりガードが非現実的になる（`search` と同事情）。よって週次ガードでは:

- `-- --skip search_finds_2017_08_21_total_eclipse --skip local_circumstances` で
  search 統合テストと局地条件テスト群を**スキップ**（高速化。錨分岐テストも `local_circumstances`
  接頭辞で一緒にスキップされるが、対象変異も下記で除外するため問題なし）。
- 上記スキップにより `local_circumstances` / `build_local_contact` の変異は週次ガードでは
  撃破されず MISSED になるため、`--exclude-re 'EclipseEngine<.*>::local_circumstances'` および
  `--exclude-re 'EclipseEngine<.*>::build_local_contact'` で**両者の全変異を除外**。

実回帰ガードは、通常 CI / tdd-workflow 工程6 の `cargo test -p umbra-eclipse`（局地条件
統合テストを**含む**）が担う。週次ミューテーションガードは下位ソルバ群のテスト有効性を見るものとし、
オーケストレーション層（`search` / `local_circumstances`）は統合テストへ委譲する（`mutation-search.md`
と同方針）。
