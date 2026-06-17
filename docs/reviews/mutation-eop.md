# Mutation review — `xtask::eop`（ISSUE-007-EOP part P2a / IERS EOP C04 取り込み）

`cargo mutants --package xtask --file crates/xtask/src/eop.rs` の生存変異の列挙と許容判断。
退行ガード `mutation.yml` は xtask を含む（生存ゼロ要求）。本 crate の純粋ロジック（`parse_eop_c04` /
`pack_records` / `unpack_records` / `build_artifact` / `mjd_to_iso` / `render_metadata`）は捕捉済み。

## 捕捉した変異（テストで殺した）

初回試走の生存 10 件のうち、以下 5 件は**狙い撃ちのテストで捕捉**（実装ロジックは load-bearing）:

- `parse_eop_c04` `tokens.len() < 7` → `<= 7`: ちょうど 7 トークンの有効データ行
  （`parse_seven_token_data_row`）。`<=` だと 7 列の正当行を取りこぼす → 捕捉。
- `parse_eop_c04` 暦範囲 `||` ×2（`!(year) || !(month) || !(day)`）: 月13/日99/年3000 の範囲外行を
  正当行に挟んで「skip されること」を検証（`parse_skips_out_of_calendar_range_rows`）。`&&` 化すると
  範囲外行が data として混入し件数が狂う → 捕捉。
- `unpack_records` `header < 0.0` → `<= 0.0`: 0 件 packed（ヘッダ 0.0）の受理を検証
  （`unpack_accepts_zero_record_buffer`）。`<=` だと正当な 0 件を誤拒否 → 捕捉。
- `unpack_records` 検証 `||`（`fract()!=0 || > MAX`）: 長さ整合バッファに小数ヘッダ 1.5 を入れ、長さ検査
  ではなく**小数検査で拒否されること**を検証（`unpack_rejects_fractional_header_with_consistent_length`）。

## 生存（許容）= IO ラッパ ＋ sanity backstop 境界

残り 5 件は**ユニットでは原理的に殺せない**等価/IO 変異。`mutation.yml` の `--exclude-re` で除外する。

1. **`unpack_records` の `header > MAX_EOP_RECORDS` 境界**（`> → ==` / `> → >=`, 2 件）:
   `MAX_EOP_RECORDS = 10_000_000` は破損データ防御の sanity 上限。境界（==/>=）を区別するには
   **ちょうど 1e7 件のヘッダ＋長さ整合（≈320 MB のバッファ）**が必要で、現実的に作れない。上限の厳密な
   等号挙動は正当性に影響しない（実データは ~2.3 万件、上限は桁違いの backstop）。**等価変異として許容**。

2. **`read_source` の本体置換**（`Ok(String::new())` / `Ok("xyzzy")`, 2 件）:
   `fs::read_to_string` の薄い IO ラッパ。`cargo test` の CWD はパッケージ dir（`crates/xtask`）で
   リポジトリルート相対パスを解決できないため**ユニットテスト不能**。パイプライン本体（parse→pack）は
   `committed_bin_matches_regeneration_from_source`（manifest 相対 `include_str!`/`include_bytes!`）で
   実効検証済み。read_source は CLI（CWD=ルート）の統合経路で行使される。nutation/vsop/elp と同方針。

3. **`verify_against_disk` の本体置換**（`Ok(())`, 1 件）:
   同上の IO ラッパ（原データ再読込→再生成→`compare_checksum`）。CWD 依存でユニット不能。比較ロジックの
   実体（再生成 == コミット済み bin）は `committed_bin_matches_regeneration_from_source` が固定。
   CLI `verify-generated --dataset eop-c04` の統合経路で行使される。nutation と同方針。

## 結論

純粋ロジック（パース・パック・アンパック検証・メタデータ・MJD→ISO）は全 caught。生存 5 件は
(a) 現実に作れない sanity backstop 境界の等価変異、(b) CWD 依存で unit 不能だが統合経路と manifest 相対
回帰テストで実効検証済みの IO ラッパ、のいずれかであり許容する（`mutation.yml` で除外）。
