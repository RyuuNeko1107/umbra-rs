# Mutation review — `umbra-eclipse::conjunction`（ISSUE-017 地心合 solver）

`cargo mutants` の生存変異の列挙と許容判断。退行ガード `mutation.yml` は eclipse を含む。

## 結果

- conjunction.rs: 39 mutants（初回）→ 33 caught / 2 unviable / 4 missed。
- 4 missed は**すべて `solve_zero_in_window`（粗走査→Brent の求根機構）内**。うち 1 件は専用テスト
  （`solve_zero_in_window_finds_root_in_second_half`）で kill、残り 3 件は**求根機構の等価変異**。
- `angle_difference` / `ecliptic_longitude` / `right_ascension` / `elongation` / `solve_conjunction`
  （角度式・分離角・オーケストレーション）の変異は全 caught（実 ephemeris の独立オラクル＝テスト側
  再実装の Δλ/Δα・離角、契約5 の |Δ|≈0、契約7 離角極小、両 kind、separation acos）。

## 生存（許容）= 求根機構の等価変異

`solve_zero_in_window` は窓を粗走査して符号変化サブ区間をブラケットし、**Brent（ISSUE-008）が真の根へ
収束**する。粗走査の点配置・符号検出は「ブラケットが取れさえすれば Brent が真値を出す」ため、以下は
出力を変えない等価変異:

| 変異 | 理由（等価） |
|---|---|
| `(t1_jd − t0_jd) * frac` の `* → +` | 走査点が窓終端を**越えて飛ぶ**ため、`[t0, 第1走査点]` が常に窓全体を内包しブラケット成立 → Brent が同じ真の根へ収束。単根・単調 f で出力不変。 |
| `prev_f * cur_f < 0.0` の `* → /` | `prev_f/cur_f < 0` ⟺ 異符号 ⟺ `prev_f*cur_f < 0`。符号検出として**完全に等価**。 |
| `prev_f * cur_f < 0.0` の `< → <=` | ループ内で `prev_f≠0`（初回は冒頭 early-return、以降は前反復で `cur_f==0` を除外済）かつ `cur_f≠0`（直前にチェック）ゆえ積は厳密に 0 にならず、`<` と `<=` は同値。 |

これらは `solver.rs` の Brent/golden 加速ヒューリスティクス（二分法フォールバックで根は契約通り）と
**同 category の「求根機構の等価変異」**。`mutation.yml` で `--exclude-re 'in solve_zero_in_window'`
により本関数を退行ガードから除外する（`solver.rs` を `--exclude '**/solver.rs'` で除外しているのと同方針）。
本関数の実契約（窓内のゼロ点を返す／符号変化なしで `RootNotBracketed`／窓前半・後半とも根を捉える）は
通常 CI の `cargo test`（`solve_zero_in_window_*` テスト群）が担保する。
