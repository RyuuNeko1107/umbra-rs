# Mutation review — `umbra-eclipse::local_maximum`（ISSUE-026 局地最大食 solver）

`cargo mutants` の生存変異の列挙と許容判断。退行ガード `mutation.yml` は eclipse を含む。

## 結果

- `local_maximum.rs`: 最大食の挙動ロジック（`m²` 合成・粗ブラケット・最小サンプル探索・d(m²)/dt=0 の
  Brent 求根・最小性・m_min・TT/UTC・単調/平底→`RootNotBracketed`）は全 caught。
- 生存（許容）= **`scan_point_count`（粗ブラケットの分割数 n の計算）内の算術変異のみ**。`local_contacts`
  と同じ走査解像度カテゴリで、`mutation.yml` の `--exclude-re 'in scan_point_count'` で除外する。

## 設計で消した等価変異（除外不要にした分）

初回の試走で `d(m²)/dt` の中心差分を `(m²(jd+h) − m²(jd−h)) / (2h)` と正規化していたため、`/(2h)` の
算術変異（`/ → *`, `/ → %`）が **missed（等価）**になっていた。Brent は根の**符号・ゼロ点しか使わず**、
正の定数 `1/(2h)` で割っても根は移動しない（スケール不変）ためである。

これを **非正規化分子 `dm2_sign(jd) = m²(jd+h) − m²(jd−h)`** に変更し、`/(2h)` を除去した。結果:

- `/(2h)` の等価変異は**そもそも生成されない**（除外不要）。
- 残る差分の `−`（`m²(jd+h) − m²(jd−h)`）は **load-bearing**: `− → +` は常に正（二乗和の和）で符号変化が
  消え Brent がブラケットできず挙動が壊れる → caught。`jd + h` / `jd − h` の符号取り違えは差分を 0 化
  または逆符号化し、最小判定が壊れる → caught。

## 生存（許容）= 走査解像度の等価/timeout 変異

`scan_point_count(t0_jd, t1_jd, step_seconds) -> usize`（`local_contacts` と同設計）は粗ブラケットの
分割数 `n = ceil(span/step).max(2)` を返すだけ。`n` は**解像度のみ**を決め、最大食の検出正否には
影響しない（n が十分大きければ単一谷を 3 点で括れる）。`span`/`n` の算術変異は等価（細かい n でも同じ
最大食時刻に収束）か timeout（巨大 n で走査停止＝検出）になる。詳細は
[`mutation-local-contacts.md`](mutation-local-contacts.md) の表と同じ。

`mutation.yml` で `--exclude-re 'in scan_point_count'` 除外（`solver.rs` 除外と同方針）。最大食の実契約
（最小性・接触の内側 c1<max<c4・窓幅不変・単調/平底→RootNotBracketed・TT/UTC）は通常 CI の
`cargo test`（`local_maximum::tests` 群）と本体の全 caught 変異が担保する。
