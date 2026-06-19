# Mutation review — `umbra-eclipse::engine` next_visible_eclipse（ISSUE-043 S8）

`EclipseEngine::next_visible_eclipse`（可視日食の遅延窓走査）と純ヘルパ
`next_visible_is_observable` の `cargo mutants` 生存変異の扱い。

## `next_visible_is_observable`（純ヘルパ・高速・ミューテーション実行済）

```
cargo mutants --package umbra-eclipse --re 'next_visible_is_observable' -- next_visible_is_observable
→ 2 mutants tested: 2 caught（0 missed）。baseline+全体 8s。
```

可視判定 `matches!(FullyVisible|PartialVisible|SunriseEclipse|SunsetEclipse)` の true/false 置換は、
高速ユニットテスト（6 値真理値表 ＋ 各バリアント個別）が撃破。**search を呼ばないので高速**＝
週次ガードに残す（観点: 可視判定はオーケストレーションでなく純ロジック）。

## `next_visible_eclipse` 本体はオーケストレーション

`next_visible_eclipse()` は `after` 以降を窓刻みで [`search`] 走査し、各日食を
[`local_circumstances`] 評価して最初の可視を返す遅延ループ。固有ロジックは
窓タイリング・`after`/重複フィルタ・可視判定（純ヘルパに分離済）・早期 return のみで、
重い計算は下位（search/local）に委譲。下位は独立にミューテーション検証済み／統合テストで縛る。

配線は統合テスト 2 件が縛る:

| テスト | 縛る対象 |
|---|---|
| `next_visible_eclipse_central_site_returns_2017_total` | happy: 可視日食を `Some` で返す・event_key・可視種別・magnitude>1・local 整合 |
| `next_visible_eclipse_skips_invisible_eclipse` | **skip**: 不可視日食(2017-08-21 from −40°S)を飛ばし後続の可視を返す（「最初の日食を無条件返却」バグの唯一のガード） |

両テストとも内部で `search`（≈300s/件）を実走し、skip は複数件解くため非常に遅い（2 件並行で ≈756s）。

## `mutation.yml`（週次ガード）での扱い

統合テスト（`next_visible_eclipse_*`）は数百秒のため、毎ミューテーションの `cargo test` に含めると
ガードが非現実的になる（search / local_circumstances と同事情）。よって週次ガードでは:

- `-- --skip next_visible_eclipse_` で 2 統合テストを**スキップ**（高速化）。`next_visible_is_observable_*`
  の高速ヘルパテストは接頭辞が異なる（`next_visible_is_observable`）ため**スキップされず**、ガードに残る。
- スキップにより撃破できなくなる本体変異を `--exclude-re 'EclipseEngine<.*>::next_visible_eclipse'`
  で除外。可視判定の純ヘルパ `next_visible_is_observable` は除外しない（高速ガード対象・上記 0 missed）。

実回帰ガードは通常 CI / tdd-workflow 工程6 の `cargo test -p umbra-eclipse`（統合テスト込み）が担う
（`mutation-search.md` / `mutation-local-circumstances.md` と同方針）。
