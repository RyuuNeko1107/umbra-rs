# mutation 運用手順（cargo-mutants）

tdd-workflow 工程7 の実行手順。**killer の選び方を間違えると baseline が timeout して
「unmutated tree で失敗」となり、1 件も走らない**。本書はその轍を残す。

## 原則: killer は「変異箇所を覆う最小のテスト集合」に絞る

cargo-mutants は**変異ごとに全テストを回す**ので、所要は `変異数 × baseline 時間`。
baseline が長いと二重に効く（timeout 値も上げる必要があり、timeout 到達分がさらに高くつく）。

## 実測（2026-09-15・umbra-eclipse / umbra-cli）

| 対象 | baseline | 備考 |
|---|---|---|
| `umbra-eclipse --lib`（全部） | **330 s** | `engine` モジュールだけで 267 s |
| `umbra-eclipse --lib` ＋ 下記 skip | **99 s** | 438/453 本が走る |
| `umbra-eclipse --test path_limits` | 103 s | path 系の killer |
| `umbra-eclipse --test search_range_filter` | 410 s | 実エンジンの search を 5 本 |
| `umbra-cli`（全部） | 236 s | 実エンジンの SLOW を含む |
| `umbra-cli --lib -- format_path_text format_limit_line` | **< 1 s** | 純関数の合成データ単体 |

**重い 3 グループ（`umbra-eclipse` lib）**:

| テスト | 本数 | 所要 |
|---|---|---|
| `next_visible_eclipse_` | 2 | 243 s |
| `local_circumstances` | 12 | 187 s |
| `search_finds_2017_08_21_total_eclipse` | 1 | 102 s |

`.github/workflows/mutation.yml` は**既にこの 3 つを `--skip` している**。ローカルで回すときも同じ
skip を付ける（付け忘れると baseline timeout で全滅する＝2026-09-15 に実際に踏んだ）。

## レシピ

```bash
# (1) lib 全体を killer にする（重い 3 グループを除外）
cargo mutants --package umbra-eclipse --re '<対象関数>' -j 2 --timeout 400 \
  -- --lib -- --skip search_finds_2017_08_21_total_eclipse \
                --skip local_circumstances --skip next_visible_eclipse_

# (2) 統合テスト 1 本だけを killer にする（その関数を縛るテストが分かっているとき）
cargo mutants --package umbra-eclipse --re '<対象関数>' -j 2 --timeout 500 \
  -- --test path_limits

# (3) 純関数は合成データの単体だけを killer にする（最速・最優先）
cargo mutants --package umbra-cli --re 'format_path_text|format_limit_line' \
  -j 2 --timeout 120 -- --lib -- format_path_text format_limit_line
```

`--timeout` は **baseline の 3〜4 倍**を目安にする（baseline 自体が timeout に掛かると 1 件も走らない）。

## timeout した変異は「caught」

ループ制御やガードの変異は、壊れると**停止しなくなる**。本プロジェクトはこれを検出とみなす
（`mutation-limb-bulge.md` / `mutation-polygon-union.md` / `mutation-path-sample-limit.md` の前例）。
`check_path_sample_limit` は全 12 変異中 2 件が timeout で、これは「ガードを外すと実際にハングする」
ことの証明になっている。

## 絞った走の結果の書き方

killer を絞ると、**絞った範囲の外の変異は当然 missed になる**。レビュー文書では
「テストが弱い」のか「killer の範囲外」なのかを**必ず区別して書く**
（例: `mutation-search-range-filter.md`）。区別せずに missed 件数だけ書くと、後任が誤読する。

## 純関数化が効く

`umbra-cli` の text 整形は実エンジン経由でしか叩けず、38 変異の完走に 3〜4 時間の見積りだった。
整形を純関数のまま**合成データで直接叩く単体テスト**を足したら、**19 変異 44 秒**になった
（`mutation-cli-path.md` の「追試」）。**mutation が重いときは、まず対象を純関数として叩けないかを疑う。**
