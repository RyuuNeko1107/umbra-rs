# Mutation review — `umbra-eclipse::eclipse_filter`（ISSUE-018 日食候補フィルタ）

`cargo mutants` の生存変異の列挙と許容判断。退行ガード `mutation.yml` は eclipse を含む。

## 結果

- eclipse_filter.rs: 初回 29 mutants → 25 caught / 2 unviable / 2 missed。
- うち 1 件（margin の符号 `+ → -`）は **grazing 部分日食フィクスチャ追加**で kill。残り 1 件は境界の等価変異。

### kill: `88 + → -`（margin の符号反転）
`limit = sum_semidiameters + parallax_moon + SAFETY_MARGIN_RAD` の最後の `+` を `-` にする変異。既知日食 4 件が
すべて**中心食**（separation 小）で、margin 無し（むしろ −margin）でも possible=true のため初回は生存。
**grazing 部分日食（2011-07-01 gamma≈1.49 / 2018-07-13 gamma≈1.35）** は合の角距離が食限（マージン抜き
`s_sun+s_moon+π_moon`）を上回り、**保守マージンが効いて初めて possible=true** になる。これらを possible=true で
縛るテスト `known_grazing_partial_eclipses_are_all_possible` を追加し、`−margin` 退行（偽陰性）を捕捉＝kill。
（偽陰性ゼロの最難ケースの検証も兼ねる。）

### 生存（許容）= 1 件・境界の等価変異

| 変異 | 理由（等価） |
|---|---|
| `eclipse_filter.rs:89 < → <=`（`separation < limit`） | 差が出るのは `separation == limit` ちょうどの一点のみ。separation（acos 由来 f64）と limit（視半径和+視差+margin の f64）が厳密一致するのは測度ゼロで、いかなるテストも踏めず無意味。BesselianPolynomial::at の EPS 境界・delta_t 区分境界と同 category の**境界等価変異**。 |

`mutation.yml` で `--exclude-re 'with <= in.*assess_eclipse_possibility'` により除外する。本フィルタの実契約
（偽陰性ゼロ＝実日食・grazing 部分日食を全採用、明確な非日食を棄却、視半径/視差/gamma/食限の独立式一致）は
通常 CI の `cargo test`（eclipse_filter テスト群）が担保する。
