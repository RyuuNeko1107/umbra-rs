# ミューテーション生存変異の許容判断 — `xtask::validate`（ISSUE-030 S30f）

`cargo mutants -p xtask -f validate.rs -- -- --skip real_engine`
（cargo-mutants 27.1.0, Docker 内）。**19 mutants: 8 caught / 5 unviable / 6 missed**。

対象: `cargo xtask validate`（ゴールデン照合を実エンジンで実走しレポート出力）。

## caught（純関数・引数解釈・結線）

- `flag_value` / `parse_format` / `parse_accuracy`（既定値・未知値→`InvalidArgument{flag,value}`・
  値欠落→`MissingArgument`。InvalidArgument は flag/value まで値検査）。
- `tolerance_profile`（Standard→standard() / Reference→reference() 写像。非公開だが単体テストで直接縛る）。
- `validate_report`（report_against_golden→render_text/render_json の format 分岐・エラー伝播）。

これら純関数は高速モックテスト（`tests/validate.rs` 7 件＋ lib 内 1 件）で生存変異ゼロ。

## 生存 6 件（実エンジン・オーケストレーション・許容）

いずれも **実エンジンを走らせる経路**で、SLOW 統合テスト
`validate_report_real_engine_one_golden_json`（実 `search`＋`local_circumstances`・≈2 分）でのみ実効検証される。
各変異を撃破するには変異ごとに実エンジンを実走する必要があり（19 変異 × 数分 = 非現実的）、
mutants 実行では `-- -- --skip real_engine` で当該 SLOW テストを除外しているため生存する。
solver.rs / EclipseEngine::search / local_circumstances を mutation 対象外とする既存方針
（mutation-search.md / mutation-local-circumstances.md）と同種の**オーケストレーション等価扱い**:

- `EngineGoldenComputer::eclipse_on` の本体置換 `-> Ok(None)`（line 116）。
- `eclipse_on` の探索窓算術 `center_jd - 0.5` の `-`→`+`/`/`、`center_jd + 0.5` の `+`→`-`/`*`（line 117-118, 計 4）。
- `run_validate` の本体置換 `-> Ok(())`（line 159）。

窓算術 4 件・本体 2 件は、窓を壊す/食を返さない/レポートを出さない変異であり、SLOW 統合テスト
（`eclipses_found==1`・`locations_compared==g.locations.len()` を assert）が確実に撃破する。
実回帰ガードは通常 CI の SLOW 統合テストが担う。

## 結論

純関数は生存変異ゼロ。残る 6 件は実エンジン経路の等価/SLOW 専用変異で、SLOW 統合テストが担保
（tdd-workflow 工程7・許容判断）。
