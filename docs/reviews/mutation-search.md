# Mutation review — `umbra-eclipse::engine::EclipseEngine::search`（ISSUE-043 S6c-ii）

`cargo mutants` の生存変異の列挙と許容判断、および退行ガード `mutation.yml` での扱い。

## search() はオーケストレーション

`search()` は確定済みの下位ソルバを配線して `Vec<SolarEclipse>` を組み立てる工程であり、
固有の数値ロジックを持たない。各下位工程は**それ自身のテストで独立にミューテーション検証済み**:

- `new_moon_candidates`（ISSUE-016）/ `solve_conjunction`（ISSUE-017）/
  `assess_eclipse_possibility`（ISSUE-018, 早期棄却）
- `classify_global_kind` / `l2_changes_sign`（S6b-iii, `mutation-global-kind.md`）
- `solve_greatest_eclipse`（S6a/S6c-i, `mutation-local-maximum.md`）
- `solve_global_contact_set`（S6b-i/ii, `mutation-global-contacts.md`）
- `BesselianPolynomial::fit`（ISSUE-022, `mutation-bessel-poly.md`）

`search()` 自身の配線は統合テスト `engine::tests::search_finds_2017_08_21_total_eclipse`
（2017-08-21 皆既を合成源で全パイプライン実走）が縛る。同テストは全球解＋ベッセル fit を
通すため約 300s 要する。

## 生存変異の列挙（`--re 'engine\.rs:13' -- --lib engine::tests::search` 限定実行）

合計 4 生存。内訳は **等価 2 / 可殺 2**。可殺分は統合テストにピン追加で撃破済み。

### 可殺（統合テストのピンで撃破）

- **`engine.rs:135:60` `*`→`+`**: `r_moon_km = k() * Re` の月半径配線取り違え。
  gamma は半径非依存ゆえ gamma 判定では捕捉できない（S5b と同じ盲点）。
  `search_finds_...` に**半径配線ピン**を追加して撃破:
  最大食 TT でベッセル多項式を評価した `l1`/`l2` が、**未変異**の
  `instantaneous_elements`（同一 config・別経路の半径計算 = 独立オラクル）と
  `1e-3` 以内で一致することを要求（fit 残差 `1e-4 ≪ 1e-3 ≪ 半径取り違え誤差`）。
- **`engine.rs:170:17` match arm `(Some(p1), Some(p4))` の削除**: fit 区間が
  全球部分食 `[P1,P4]` から候補窓フォールバックへ落ちる。候補窓も最大食を bracket するため
  「最大食 bracket」判定だけでは見逃す。`search_finds_...` に
  **`fit_interval.start == P1` / `fit_interval.end == P4` の厳密一致ピン**を追加して撃破。

### 生存（許容）= 粗走査局在による根許容差変換の等価変異

- **`engine.rs:130:66` `/`→`%` / `/`→`*`**: `RootConfig.x_tolerance_days =
  root_tolerance_seconds / SECONDS_PER_DAY`（秒→日換算）。
  この許容差は最大食・接触の求根（Brent）の停止幅にのみ効く。最大食 TT は
  `solve_local_maximum` の**粗走査が ±60s 程度に局在**させた後で Brent 精解に渡るため、
  Brent の停止幅を変えても返る根は粗走査分解能の範囲で不変＝挙動不変。
  下位ソルバの精度はそれぞれのテストで検証済み。`solver.rs` /
  `solve_zero_in_window` / `scan_point_count` の求根機構除外と同方針の等価変異。

## `mutation.yml`（週次ガード）での扱い

`search_finds_...` は約 300s で、これを毎ミューテーションの `cargo test` に含めると
umbra-eclipse の**全**変異が 300s 超になりガードが非現実的になる。よって週次ガードでは:

- `-- --skip search_finds_2017_08_21_total_eclipse` で当該統合テストを**スキップ**（高速化）。
- 上記スキップにより `search` の可殺 2 変異は週次ガードでは撃破されず MISSED になるため、
  `--exclude-re 'EclipseEngine<.*>::search'` で **`search` の全変異を除外**（等価 2 + 可殺 2）。

可殺 2 変異の実回帰ガードは、通常 CI / tdd-workflow 工程6 の
`cargo test -p umbra-eclipse`（統合テストを**含む**）が担う。週次ミューテーションガードは
あくまで下位ソルバ群のテスト有効性を見るものとし、オーケストレーション層は統合テストへ委譲する。
