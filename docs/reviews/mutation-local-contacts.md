# Mutation review — `umbra-eclipse::local_contacts`（ISSUE-025 局地接触 C1/C4・C2/C3 solver）

`cargo mutants` の生存変異の列挙と許容判断。退行ガード `mutation.yml` は eclipse を含む。

## 結果

- `local_contacts.rs`: 全 mutants のうち、**`scan_sign_change_roots` の挙動ロジック（符号変化検出・
  Brent 求根・接触種別割当・C1/C4/C2/C3 分類・None 分岐・TT/UTC）は全 caught**。
- 生存（許容）= **`scan_point_count`（粗走査の分割数 n の計算）内の算術変異のみ**。1 missed（等価）＋
  数件 timeout（巨大 n）。これらは `mutation.yml` で `--exclude-re 'in scan_point_count'` 除外する。

> 注（履歴）: 本 issue の初回コミット（6b1cff9）時のミューテーション報告「0 missed」は、実際には
> 並行実行された別ランと mutants.out を共有して**未完了の中間状態を誤読**したものだった。完走させた
> 結果、下記の `scan_point_count` 等価変異が生存と判明したため、解像度計算を専用ヘルパへ切り出して
> 除外・文書化する形に是正した（コードの振る舞いは不変）。

## 生存（許容）= 走査解像度の等価/timeout 変異

`scan_point_count(t0_jd, t1_jd, step_seconds) -> usize` は粗走査の分割数 `n = ceil(span/step).max(2)`
を返すだけの関数。`n` は**走査解像度のみ**を決め、接触検出の正否には影響しない（偽陰性回避は刻みの
細かさ＝十分大きい n で担保され、n がさらに大きくても／窓全体を別の倍率で刻んでも、符号変化区間は
同じく捉えられ Brent が同じ真根へ収束する）。したがって `span`/`n` の算術には観測可能な振る舞い契約が
無い:

| 変異（`scan_point_count` 内） | 区分 | 理由（等価 / 検出） |
|---|---|---|
| `(t1_jd − t0_jd)` の `− → /` | missed（等価） | `t1/t0 ≈ 1`（JD≈2.46e6）→ `span ≈ 86400` → `n ≈ 2880`。実窓 [t0,t1] は `jd_at` 側で正しく使われるため、**より細かい走査**になるだけで全接触を同じく検出。出力不変。 |
| `(t1_jd − t0_jd)` の `− → +` | timeout | `span ≈ 4.9e6 × 86400` → `n` 巨大 → 走査ループが事実上停止 → 制限時間で検出（振る舞いは壊れている）。 |
| `* SECONDS_PER_DAY` の `* → +` / `* → /` | timeout/等価 | n を巨大化（timeout）または別倍率（細かい走査＝等価）。いずれも解像度のみ。 |
| `/ step_seconds` の `/ → *` | timeout | n を巨大化 → timeout 検出。 |

これらは `solver.rs`（Brent/golden の数値機構）や `conjunction::solve_zero_in_window` と**同 category の
「走査解像度／求根機構の等価変異」**。`scan_point_count` は解像度を決めるだけで挙動契約を持たないため、
`mutation.yml` で `--exclude-re 'in scan_point_count'` により退行ガードから除外する。

接触検出の実契約（C1/C4/C2/C3 の符号方向・順序・None 分岐・grazing 偽陰性ガード・窓幅不変・TT/UTC）は
通常 CI の `cargo test`（`local_contacts::tests` 群）と、`scan_sign_change_roots` 本体の全 caught 変異が
担保する。
