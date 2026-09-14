# mutation レビュー: `umbra path` CLI（ISSUE-048）

対象: `crates/umbra-cli/src/lib.rs` の `run_path` / `format_path_text` / `format_limit_line` /
`check_path_sample_count`。正本は `docs/issues/ISSUE-048-cli-path.md`。

## 実行

```
docker compose -p umbra-rs run --rm rust bash -c \
  "cargo install cargo-mutants --locked && \
   cargo mutants --package umbra-cli \
     --re 'run_path|format_path_text|format_limit_line|check_path_sample_count' -j 2 --timeout 600"
```

## 結果（部分実行・2026-09-14）

**38 変異中 6 件を実行した時点で中断**（下記「中断の理由」）。実行分は **全て caught**（うち 5 件は TIMEOUT）。

| 変異 | 結果 |
|---|---|
| `check_path_sample_count -> Ok(())`（検査を無効化） | **TIMEOUT** |
| `check_path_sample_count` の `span / interval` → `%` / `*` | **TIMEOUT** |
| `check_path_sample_count` の `>` → `==` / `<` | **TIMEOUT** |

`check_path_sample_count` が生む変異は上記 5 件で**全て**であり、**いずれも TIMEOUT＝ハング検出で caught**。
これは本ガードの存在意義そのものの証明になっている: **ガードを壊すと `run_path` が実際に停止しなくなる**
（`run_path_tiny_interval_is_interval_too_small` が返らなくなる）。本プロジェクトで timeout を caught として
扱う既存方針（ループ制御の変異＝`mutation-limb-bulge.md` / `mutation-polygon-union.md`）と同じ扱い。

## 追試（2026-09-14・FAST 分離後）

**中断の原因はテスト構成**だった。text 整形は**純関数**（`format_path_text` / `format_limit_line`）なので、
合成 `SolarEclipse` / `EclipsePath` を直接渡す **FAST 単体テスト 13 本**を追加し、mutation の killer を
FAST のみに絞れるようにした（実エンジンを一切起動しない）。

```
cargo mutants --package umbra-cli --re 'format_path_text|format_limit_line'   -j 2 --timeout 120 -- --lib -- format_path_text format_limit_line
```

| 総数 | caught | missed | 所要 |
|---|---|---|---|
| 19 | 17 | **2** | **44 秒** |

**3〜4 時間 → 44 秒**。整形系の変異は**全て caught**（ラベルの脱落・重複・順序、`None`/空の
「なし」分岐 4 種、点数の取り違え、始点終点の入れ替え、環数と外環頂点数の転置、`samples` の 0 表示、
最大食の時刻・座標の出所と lat/lon 順、`format_limit_line` の `n > 0` 境界）。

missed 2 件は `run_path` の `PathOptions` 構築（`sample_interval_seconds` / `include_limits` の
フィールド削除＝`Default` 値へフォールバック）で、**本走では killer に含めていない**
（`run_path` のテストは実エンジン依存＝SLOW）。別走で `run_path` 側を確認した結果は下記。

## `run_path` 側の走（SLOW killer・15 変異）

```
cargo mutants --package umbra-cli --re 'run_path' -j 2 --timeout 600   -- --lib -- run_path check_path_sample_count
```

| 総数 | caught | missed | 所要 |
|---|---|---|---|
| 15 | 13 | **2** | 35 分 |

| missed | 判定 |
|---|---|
| `PathOptions` から `sample_interval_seconds` を削除（`--interval` が無視され既定 60 s になる） | **非等価・撃破済み**。`run_path_interval_is_wired_into_path_options_finer_gives_more_points` を追加（120 s と 360 s で中心線の点数を比較。変異時は両方とも既定 60 s になり点数が一致するので厳密不等号が破れる） |
| `date.jd2().jd() + 1.0` → `* 1.0`（探索範囲が幅ゼロに潰れる） | **等価ではないが撃てない**。原因は `run_path` ではなく **`search` の範囲セマンティクス**（下記） |

## 派生して発見した `search` の範囲セマンティクス（別件・要対応）

`+ 1.0` → `* 1.0` が撃てないのを追ったところ、**`search` は結果を要求範囲で絞り込んでいない**ことが判明した。
`new_moon_candidates` が平均朔の窓（`WINDOW_HALF_WIDTH_DAYS = 1` 日）で範囲を**両側に広げて**候補を採り、
`search` はその候補から得た日食を**範囲で再フィルタせずに**返す。

実測（`standard_engine`・2024-04-08 の皆既）:

| 要求範囲 | 返った件数 |
|---|---|
| `[2024-04-06, 04-07)` | 0 |
| `[2024-04-07, 04-08)` | **1**（`2024-04-08#300`） |
| `[2024-04-08, 04-09)` | 1（`2024-04-08#300`） |
| `[2024-04-09, 04-10)` | **1**（`2024-04-08#300`） |

つまり `umbra path --date 2024-04-07` / `--date 2024-04-09` は **4/8 の日食を「その日の日食」として返す**。
`run_local`（ISSUE-032）も同じ日付解決なので同様。`umbra search --from/--to` も範囲外の日食を返しうる。
**CLI の「指定日に起こる日食」という契約と食い違う**ので、別 issue として起票して対応する。


## 当初の中断とその理由（記録）

`umbra-cli` のテストスイートは**実エンジンの SLOW テスト**（1 日範囲の `search` を複数回）を含むため
baseline だけで約 236 秒かかる。cargo-mutants は変異ごとに全スイートを回すので、38 変異の完走には
**3〜4 時間**かかる見積りだった（timeout 到達分は 1 件 600 秒）。マシン負荷を優先して中断した
（`docker-light` の運用方針）。

**当時未実行だった 32 件**は `run_path` / `format_path_text` / `format_limit_line` の変異で、
**text 整形のラベル・件数・分岐**が中心。これらは ISSUE-048 節の 20 テスト
（`--no-limits` で「なし」・`None` 要素の明示・geojson の `role` 集合一致・実 2024 の種別/点数）が
**通常の `cargo test` で**縛っている。ただし当時は **mutation で撃てたことを確認できておらず**、
「テストは存在するが判別力は未証明」の状態だった（→ 上の「追試」「`run_path` 側の走」で解消）。

**残作業だったもの**（上の「追試」で解消済み）: CLI のテストを FAST（合成 fixture）と SLOW（実エンジン）に
分離し、mutation は FAST のみを killer にして完走させる。**整形系は完了**（19 変異 44 秒・整形の生存 0）。
`run_path` のオーケストレーション（日付解決・エンジン構築・no-eclipse 分岐）は実エンジン依存のままで、
ここだけは SLOW が要る。
