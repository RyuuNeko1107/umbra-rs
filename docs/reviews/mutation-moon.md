# ミューテーション生存変異の許容判断 — `umbra-ephemeris::moon`（ISSUE-014 S2）

`cargo mutants --package umbra-ephemeris --file crates/umbra-ephemeris/src/moon.rs`（cargo-mutants
27.1.0, Docker 内）。**392 mutants: 378 caught / 3 unviable / 1 timeout / 10 missed**。

オラクル = 著者 Fortran `elp82b_1` の独立 Python 移植（公開リファレンス値 JD 2451555.5 に 5e-10 km
一致）。`moon_geocentric_j2000` の XYZ を 1400–2100 の 7 エポックで照合、許容 **2e-6 km**。
この許容は Rust↔移植の f64 演算順差の実測最大 **6.65e-7 km**（1900, t≈−1.4。37872 項総和の固有丸め）
の ~3 倍。月モデル精度 0.40″≈750m の ~9 桁下。

## 生存 10 件 ＋ timeout 1 件（いずれも許容）

### A. `% TAU → + TAU`（sin 周期性で**等価**）: lines 302, 328, 339
`(y % TAU).sin()` の `%`→`+`。`sin(y + 2π) ≡ sin(y)`（数学的恒等、f64 でも ~1e-16 相対）。
`% TAU` は大引数の精度のためのレンジ縮約で、1 周期加算は sin を変えない。**真の等価変異**、殺せない。

### B. DTASM 経由 tgv の微小補正（**効果が f64 床以下**）: lines 173（×4）, 288
`DTASM = 2·ALFA/(3·AM)` の定数演算（`/`→`%`,`/`→`*`,`*`→`+`,`*`→`/`）と `tgv = c[1] + DTASM·c[5]`
の `+`→`−`。DTASM は主問題振幅補正の `tgv` にのみ入り、`tgv` は `tgv·(DELNP − AM·DELNU)` にのみ
効く。`DELNP − AM·DELNU ≈ 6×10⁻¹¹ rad` のため、DTASM を 88 倍に変えても最終 XYZ への寄与は
< 6.65e-7 km（補正の補正＝サブ mm）。許容床以下で原理的に XYZ テストでは捕捉不可。

### C. Laskar 回転正規化 `ra` の 2 次項（**効果が f64 床以下**）: lines 358:25, 358:30
`ra = 2·√(1 − pw² − qw²)` の `−`→`+` / `pw·pw`→`pw+pw`。1900–2100 で `pw,qw ≈ 1e-5..1e-4`、
よって `pw²,qw² ≈ 1e-10..1e-8`。これらを乱しても `ra` は ~1e-8 しか動かず、`pw_final = pw·ra` は
~1e-13、XYZ への寄与 < 6.65e-7 km。サブ mm で許容床以下。

### D. 復号インデックスの病的変異（**timeout で検出**）: line 248
`values[idx..idx + n_mult/n_coeff]` の `+`→`*`。`idx*n` で系列長/項数が破壊され巨大確保 or 範囲外 →
20s 以内に test が通らない（timeout＝検出）。サイレントな生存ではない。

## 判断
A は等価、B/C は効果が f64 評価床（6.65e-7 km）以下かつ月モデル精度（750m）の 9 桁下で**物理的に
無意味**、D は timeout で検出。構造・論理（系列群の引数組合せ・距離 cos 位相・×t/×t² スケール・
Laskar 回転構造・単位・DE200 補正 1 次・packed 復号）は 378 caught で全て捕捉済み。よって**全件許容**。
小振幅係数の転記精度は ISSUE-034 の round-trip + checksum が別途保証する。
