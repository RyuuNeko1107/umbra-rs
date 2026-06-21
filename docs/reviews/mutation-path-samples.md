# mutation レビュー: 経路サンプル列 samples（`sample_central_point` / `trace_central`）

対象: ISSUE-045 M9.7（中心食の経路サンプル列 `EclipsePath.samples`）の核心関数
`sample_central_point`（PathSample 構築＝per-sample の duration/path_width/kind/sun_altitude/time_utc）と
`trace_central`（lockstep サンプル収集ループ）。いずれも `crates/umbra-eclipse/src/engine.rs`。

## 実行
```
cargo install cargo-mutants --locked
cargo mutants -p umbra-eclipse \
  --re 'sample_central_point|trace_central' --no-shuffle \
  -- --test path_limits -- \
     samples_are samples_field samples_time include_limits_false noncentral_eclipse \
     synthetic_limits northern_limit limits_are central_eclipse_with_limits
```
killer は path_limits.rs の **FAST 統合テスト**（samples の lockstep・独立オラクル〔duration の
`2|L2'|/|rel|`・width の大圏距離・kind 符号〕・time_utc=tt_to_utc＆単調・include_limits=false 空・
非中心空、＋既存の南北限界線 FAST）。SLOW（実 2024・各 ~130s）は mutation の再走に不適ゆえ
テスト名フィルタで除外（`real_2024` を含まない FAST 名のみ列挙）。

## 結果（2026-06-21・工程7）
**55 mutants: 44 caught・3 missed・6 unviable・2 timeouts。**

### timeout 2 件（＝ハング検出・caught 相当・許容）
| 行 | 変異 | 区分 |
|---|---|---|
| engine.rs:656 `\|\|` → `&&`（`t_sec >= span_seconds \|\| interval_seconds <= 0.0`） | ループ終端 | 終端条件破壊で無限ループ＝**20s timeout でハング検出**（M9.1/M9.3 と同方針）。 |
| engine.rs:659 `+` → `*`（`t_sec = (t_sec + interval_seconds).min(span_seconds)`） | ループ前進 | 前進が壊れ進まない/発散＝無限ループ＝timeout で検出。 |

`trace_central` の刻み前進・終端は **解像度/終端機構**（サンプル点の総数を決めるが真の経路位置は
`sample_central_point` が決める）。これらは値ではなくハングで検出されるため CI 除外は付けない
（M9.1 以来この扱い・trace_central は mutation.yml で未除外のまま timeout 検出に委ねる）。

### 生存 3 件（すべて到達不能境界の等価変異・許容）
| # | 行 | 変異 | 理由 |
|---|---|---|---|
| 1 | engine.rs:741:20 | `rel_speed > 0.0` の `>` → `>=` | `rel_speed = rel_x.hypot(rel_y)` が**厳密 0.0** になるのは rel_x・rel_y がともに厳密 0 のときのみ。実/合成中心食では rel≈0.2 Re/h で `>` と `>=` は同分岐＝到達不能。**M9.6 `central_width_and_duration:249` と同型の既許容等価**。 |
| 2 | engine.rs:741:26 | `rel_speed > 0.0 && rel_speed.is_finite()` の `&&` → `\|\|` | `&&`→`\|\|` が分岐を変えるのは「`rel_speed==0.0`（有限）」または「`rel_speed=+∞`」のときのみ。前者は #1 と同じく到達不能、後者は rel 成分が無限大＝瞬時要素が非有限となる構成で実/合成中心食とも生成不能。**M9.6:249 と同型の既許容等価**。 |
| 3 | engine.rs:749:23 | `l2p < 0.0`（kind） の `<` → `<=` | `l2p = l2 − ζ₀·tan f2` が**厳密 0.0**（皆既↔金環の hybrid 遷移点を bit-exact に踏む）ときのみ Total/Annular が変わる。ζ₀ は root-finder 出力ゆえ `l2 == ζ₀·tan f2` を bit 一致で踏む入力は測度ゼロ＝オラクルで到達不能。境界の語義（`l2p<0`→Total・`l2p==0`→Annular）は既存 `global::classify`（`l2_exactly_zero_is_annular_not_total`）と一致しコードで固定済み。**到達不能境界の等価**。 |

注: #3 を狙う「金環中心食（l2>0）」テストでも `l2p>0` では `<` と `<=` は同分岐ゆえ kill 不能
（差は `l2p==0` の一点のみ）。よって追加オラクルでの撃破は原理的に不可能で、等価と確定。

## CI 除外提案（mutation.yml）
等価 3 件のみを**演算子・関数を厳密指定**して除外（`sample_central_point` は比較演算子が複数行に
散在＝南北割当の `edge_a.lat >= edge_b.lat` など killable な比較も持つため、`with <= in.*sample_central_point`
のような関数レベル粗除外は使わず、`replace X with Y` 形で当該等価変異のみを的確に除外）:
```
--exclude-re 'replace > with >= in.*sample_central_point'
--exclude-re 'replace && with \|\| in.*sample_central_point'
--exclude-re 'replace < with <= in.*sample_central_point'
```
実回帰ガード: 通常 CI の `cargo test -p umbra-eclipse`（FAST 独立オラクル＋実 2024 SLOW）が
duration/width/kind/time_utc・lockstep を縛り続ける。
