# §1 時刻系変換（Time scales）

本節は `umbra-rs` の時刻系変換（Gregorian↔JD、UTC↔TAI↔TT↔TDB↔UT1、ΔT 供給）の数式・手順・符号規約を定める。
正本は `algorithms.md §0`（記号表）。本節はそれに厳密準拠し、本節固有の補助記号のみ追加定義する。
関連 Issue: ISSUE-004（JulianDate2）/ 005（暦変換）/ 006（UTC/TAI/TT）/ 007（UT1/ΔT/EOP）/ 042（TimeData/TimeScales 束ね）。

> 状態: ドラフト（Milestone 0）。一次資料で式番号を確認できない箇所は「**要確認**」を残し、推測で式番号を書かない。

---

## 目的と入出力（単位・時刻系・座標系・フレーム）

- **目的**: 公開入力（UTC のグレゴリオ暦）から、天体計算が要求する各時刻系（TT, TDB, UT1, TAI）を、桁落ちを避けた 2 要素ユリウス日 `JulianDate2`（conventions §6）で生成する。あわせて ΔT = TT − UT1 とその不確実性帯（accuracy.md §0）を供給する。
- **入力**: グレゴリオ暦 (y, mo, d, h, mi, s)（UTC として解釈、conventions §3/§6）。または各 `*Instant`。
- **出力**: `JulianDate2`、各 `*Instant`（`UtcInstant`/`TaiInstant`/`TtInstant`/`Ut1Instant`）、ΔT（秒）、不確実性帯（秒）。J2000 からの TT ユリウス世紀 `t`（algorithms.md §0、`t = ((JD−2451545.0))/36525`）の生成もここに属する。
- **単位**: 日（`JulianDate2`）、秒（SI、時刻差・ΔT・閏秒）。角度・距離は本節に現れない。
- **座標系・フレーム**: 該当なし（時刻層）。極運動（xp, yp）の供給は ISSUE-007 だが ITRS フレーム連鎖（§2/§10）で消費される。
- **数値方針**: 桁落ち対策は numerical-policy §A1（`t` は必ず `JulianDate2` から生成）。生 `f64` の JD・時刻を関数間で渡さない（conventions §6/§11）。

---

## 記号（algorithms.md §0 を参照。本節固有の補助記号のみ定義）

algorithms.md §0「時刻系」表（ΔT, ΔAT, TT−TAI, ERA, GAST）を正本とする。本節で補助的に用いる記号:

| 記号 | 意味 | 出典/備考 |
|---|---|---|
| JD | ユリウス日（日）。`JulianDate2` の合算値 = part1 + part2 | conventions §6 |
| part1, part2 | `JulianDate2` の 2 要素。part1 = 整数寄り基準、part2 = 微小差 | ISSUE-004。正規化規約 §B1 |
| MJD | 修正ユリウス日 = JD − 2400000.5 | SOFA `iauCal2jd` の基準（part1=2400000.5） |
| (UT1−UTC) | EOP 由来。秒。δUT1 とも | ISSUE-007、IERS EOP C04 |
| (TAI−UTC) | 閏秒積算（= ΔAT）。整数秒。1972– | ISSUE-006、IERS Bulletin C |
| TDB−TT | 周期項（地心、Fairhead & Bretagnon 1990）。±約 1.7 ms 振幅 | SOFA `iauDtdb`、§E |
| u, v | `iauDtdb` の地上項用 観測者成分（東距離・赤道面距離, km）。**本節では地心なので 0** | SOFA `iauDtdb`。algorithms.md §0 の観測者 u,v とは別物（衝突注意、§E の注記参照） |
| τ | 光行時間オフセット（秒）。part2 に加える微小量 | numerical-policy §A3。本節では加算規約のみ |

> **記号衝突の注意**: `iauDtdb` の引数 `u, v` は algorithms.md §0 観測者欄の `u = x − ξ, v = y − η` とは別概念。本節では地心 TDB のため両者とも 0 を渡し、混同を避けるため §E で明示する。記号表への追記候補は末尾「報告」に記載。

---

## 数式（番号付き。各式に出典）

### A. グレゴリオ暦 ↔ JD（ISSUE-005、Meeus Ch.7）

プロレプティック・グレゴリオ暦で統一（1582 改暦は扱わない。ユリウス分岐に落とさない）。天文学的年番号（1 BC = 0 年, 2 BC = −1 年, ISO 8601）。

**(1.A1) 暦 → JD（Meeus AA 2nd ed., Ch.7「Julian Day」, 式 (7.1) のグレゴリオ分岐）**

mo ≤ 2 のとき y ← y − 1, mo ← mo + 12 とし、

```
A = floor(y / 100)
B = 2 − A + floor(A / 4)              # グレゴリオ補正（常にこの分岐を採用）
JD_int = floor(365.25·(y + 4716)) + floor(30.6001·(mo + 1)) + d + B − 1524.5
```

時分秒は日小数として分離（§B2 参照）:

```
day_frac = (h·3600 + mi·60 + s) / 86400          # 86400 = SI 日（定義値, conventions §4.1 系）
```

`JulianDate2` は part1 = JD_int（0.5 起算＝正午基準を含む整数+0.5）, part2 = day_frac で構築し正規化（§B1）。

- 出典: Meeus AA 2nd ed., Ch.7, 式 (7.1)。グレゴリオ判定（Meeus の改暦境界分岐）は使わず**常にグレゴリオ**（ISSUE-005 非目的）。
- 照合: SOFA `iauCal2jd`（MJD 基準・グレゴリオ・範囲チェックあり）。SOFA と Meeus の年番号規約差（負年・0 年）を実装コメントに明記（ISSUE-005）。

**(1.A2) JD → 暦（Meeus AA 2nd ed., Ch.7 の逆変換手順）**

```
JD' = JD + 0.5
Z = floor(JD')                       # 整数部
F = JD' − Z                          # 小数部（日内）
# グレゴリオ前提（α 補正を常に適用）:
α = floor((Z − 1867216.25) / 36524.25)
A = Z + 1 + α − floor(α / 4)
B = A + 1524
C = floor((B − 122.1) / 365.25)
D = floor(365.25·C)
E = floor((B − D) / 30.6001)
d   = B − D − floor(30.6001·E) + F   # 日（小数含む）
mo  = E − 1   (E < 14) ／ E − 13  (E = 14 or 15)
y   = C − 4716 (mo > 2) ／ C − 4715 (mo ≤ 2)
```

時分秒は F（small part は part2 由来）から復元。小数秒の丸めで 60 秒に達しないようガード（§F）。

- 出典: Meeus AA 2nd ed., Ch.7（Z, F, α, A, B, C, D, E からの復元手順）。**要確認**: 逆変換の各中間式は Meeus 本文記載の手順だが、一次資料で個別の式番号が付されているか未確認（式番号は推測しない）。
- 照合: SOFA `iauJd2cal`。

### B. JulianDate2（2 要素表現、ISSUE-004、SOFA 二引数規約）

**(1.B1) 正規化規約**: part2 ∈ [−0.5, 0.5)（または [0,1)。conventions §6 と numerical-policy §A1 に一致させる。**要確認**: プロジェクト既定を [−0.5,0.5) に固定するか [0,1) かをレビューで一意化）。part1 は整数寄り基準（J2000 相対整数日でもよい）。

**(1.B2) 秒の加算（桁分離、numerical-policy §A1/§A3）**

```
add_seconds(s):  days = s / 86400.0;  part2 ← part2 + days;  再正規化（part1 は直接いじらない）
```

光行時間 τ 等の微小オフセットは必ず part2 へ加え、part2 が肥大したら正規化で part1 へ繰り上げる。

- 出典: SOFA `(date1, date2)` 二引数規約（`iauCal2jd`/`iauJd2cal` 等が `date1+date2 = JD` の可搬分割を採用）。桁分離加算は two-sum / compensated summation（Dekker 1971, Knuth）の考え方（ISSUE-004）。

**(1.B3) 高精度差分（日）**

```
diff_days(a, b) = (a.part1 − b.part1) + (a.part2 − b.part2)
```

大きい part1 同士を先に引いて桁落ちを抑える順序を固定（ISSUE-004）。

**(1.B4) TT ユリウス世紀（algorithms.md §0、numerical-policy §A1）**

```
t = (jd_tt.diff_days(J2000)) / 36525.0           # centuries（J2000 = JulianDate2(2451545.0, 0.0)）
T_millennia = (jd_tt.diff_days(J2000)) / 365250.0
```

エポック 2451545.0 の減算は **part1 側（整数部）で厳密に**行う。単一 f64 JD（≈2.45e6）での減算は ulp ≈ 5.4e-10 day ≈ 4.6e-5 s を失い ±1s 目標を侵食する（numerical-policy §A1）。**`t` は必ず `JulianDate2` から生成**する（生 f64 JD 経由を禁止）。

### C. UTC ↔ TAI ↔ TT（ISSUE-006）

**(1.C1) TT − TAI = 32.184 s（定義定数）**

```
TT = TAI + 32.184 s          # conventions §4.1 物理定数表（IAU 1991 決議, 定義値）
TAI = TT − 32.184 s
```

magic number ではなく定義定数として参照（conventions §4.1 / §11）。

- 出典: IAU 1991 決議、IERS Conventions 2010 (TN 36) §1、SOFA `iauTaitt`/`iauTttai`。

**(1.C2) TAI − UTC = ΔAT（閏秒積算、整数秒、1972–）**

```
TAI = UTC + (TAI−UTC)        # (TAI−UTC) は閏秒テーブルから区間引き（SOFA iauDat 相当）
```

UTC の閏秒挿入日は 86401 秒を持つ（23:59:60 が存在）。SOFA の "quasi-JD" 規約に従い、閏秒日の非線形性を扱う。

- 出典: 閏秒テーブル = IERS Bulletin C / IANA `leap-seconds.list`（data-sources §3.2）。UTC↔TAI 境界処理 = SOFA `iauUtctai`/`iauTaiutc`。閏秒テーブルは SOFA `iauDat` に対応。

**(1.C3) 合成**

```
TT  = UTC + (TAI−UTC) + 32.184 s            # utc_to_tt
UTC = TT − 32.184 s − (TAI−UTC)             # tt_to_utc（閏秒は該当日で逆引き）
```

全加減算は `JulianDate2::add_seconds`（§B2）で µs を保つ。

### D. UT1 と ΔT（ISSUE-007）

**(1.D1) UT1 の取得**

```
UT1 = UTC + (UT1−UTC)        # (UT1−UTC) は EOP（IERS C04）から補間。秒。
```

- 出典: IERS Conventions 2010 ch.5（時刻系定義）。EOP C04 = IERS Earth Orientation Center（data-sources §3.1）。系列 = **EOP 14 C04** または **EOP 20 C04**（採用版を `series_version()` と `CalculationMetadata` に固定、ISSUE-007）。

**(1.D2) ΔT の定義と EOP 由来導出**

```
ΔT = TT − UT1
   = (TAI−UTC) + 32.184 − (UT1−UTC)         # EOP coverage 内はこの恒等式で高精度
```

- 出典: conventions §6 / IERS Conventions 2010。EOP coverage 内の ΔT 履歴は EOP 由来（accuracy.md §3.3）。

**(1.D3) 長期 ΔT（Espenak–Meeus 区間別多項式）**

EOP coverage 外（1972 以前・将来外挿）は Espenak–Meeus の区間別多項式を Horner 評価。各区間は `year` の補助変数（例 `u = (year − 1820)/100` 等、区間ごとに定義）に対する多項式。区間境界の連続性に注意。

- 出典: **Espenak & Meeus, "Five Millennium Canon of Solar Eclipses: −1999 to +3000"（NASA TP-2006-214141）** の ΔT polynomial expressions（NASA `eclipse.gsfc.nasa.gov/SEcat5/deltatpoly.html`、区間別 −1999〜+3000）。各区間係数は出典明記で取り込む（data-sources §3.3）。年 −500 以前は `ΔT = −20 + 32·u²`, `u = (year − 1820)/100`。
- **要確認**: 採用する**不確実性式**（外挿の標準誤差）の一次出典（NASA 記載の uncertainty 式 / Morrison & Stephenson 2004 ベース）。式そのものは推測しない。

**(1.D4) ΔT の供給合成（推奨・別型、ISSUE-007/042）**

```
EOP coverage 内: (1.D2) の恒等式（高精度、不確実性 <0.1 s）
EOP coverage 外: (1.D3) Espenak–Meeus 外挿（不確実性は年数に応じ増大、~数秒/2100）
切替境界の不連続を許容内に収める。
```

不確実性帯は `CalculationMetadata.delta_t_uncertainty_seconds` に必ず出力（accuracy.md §0、欠落させない）。

### E. TDB（Reference バックエンド用、SOFA `iauDtdb`）

**(1.E1) TT ≈ TDB（地心、周期項）**

```
TDB = TT + (TDB−TT)
TDB − TT = iauDtdb(jd_tt.part1, jd_tt.part2, ut, elong, u_obs, v_obs)
```

地心では地上項（観測者依存）を無効化するため `u_obs = v_obs = 0`（および `elong = 0`）を渡し、Fairhead & Bretagnon (1990) の地心周期項のみを得る。振幅は約 ±1.7 ms（algorithms.md §0 注記「TT ≈ TDB ± 2 ms」と整合）。

- 出典: SOFA `iauDtdb`。地心モデル = **Fairhead & Bretagnon (1990)** フル形（約 800 項、現代で数 ns 精度）。地上項 = Moyer (1981) / Murray (1983)、基本引数 = Simon et al. (1994)。
- **適用範囲**: TDB は **Reference プロファイル（JPL DE）専用**（accuracy.md §1）。Standard（VSOP87D/ELP-MPP02、TT 基準）では TDB を経由しない。周期項を無視（TDB≈TT）してよい局面は Reference 以外（algorithms.md §0 の「TT≈TDB±2ms」は近似許容の根拠）。
- **u, v 記号衝突**: ここでの `u_obs, v_obs` は `iauDtdb` の観測者成分であり、algorithms.md §0 観測者欄の `u = x − ξ, v = y − η`（ベッセル相対位置）とは別物。地心では 0。

---

## 手順（実装順・数値注意）

1. **暦 → UTC `JulianDate2`**: (1.A1) で整数日 part1、時分秒 part2。入力検証（§F）。生 f64 JD を作らない（numerical-policy §A1）。
2. **UTC → TAI → TT**: (1.C2)(1.C1)(1.C3)。閏秒テーブルを versioned + checksum でロード（実行時ネットワーク禁止、accuracy.md §5）。`add_seconds`（§B2）で µs 保持。valid_to 超過は `MissingLeapSecondData`（ISSUE-042 確定A2: metadata 記録）。
3. **UTC → UT1**: (1.D1)。EOP coverage 外は `MissingEarthOrientationData`（metadata 記録）。
4. **ΔT 供給**: (1.D4)。coverage 内は (1.D2) 恒等式、外は (1.D3) Espenak–Meeus。不確実性帯を必ず metadata 出力（accuracy.md §0）。
5. **TT 世紀 `t`**: (1.B4) で `JulianDate2` から直接生成。エポック減算を part1 で厳密に（numerical-policy §A1）。**級数引数の事前 mod 2π は禁止**（numerical-policy §A1: 丸め誤差を注入する）。
6. **TDB（Reference のみ）**: (1.E1)。地心は `u_obs=v_obs=0`。

> **接触時刻は TT と UTC の両方を返す**（conventions §6）。将来日食の UTC は ΔT/UT1 予測律速（accuracy.md §0(b)）のため、TT 基準を一級保持し UTC を併記する。

---

## 境界・特異・異常系

- **負の年・0 年**: 天文学的年番号で往復一致（1 BC = 0, 2 BC = −1）。Meeus/SOFA と一致（ISSUE-005）。
- **JD 0.5 境界（正午起算）**: 暦日の真夜中/正午で part1/part2 の繰り上げを厳密に。
- **閏秒境界（23:59:60）**: 双方向で一意写像。23:59:59 → 23:59:60 → 00:00:00 の TAI 差が各 1 s で連続（ISSUE-006）。閏秒日以外の秒=60 は `InvalidDate`。
- **1972 以前の UTC**: UTC が階段閏秒でない年代。ΔT（Espenak–Meeus）経由で UT1/TT へ橋渡し（ISSUE-006/007、**要確認**: `UnsupportedTimeRange` とするかの境界仕様はレビュー確定）。
- **EOP / 閏秒 valid_to 超過**: 沈黙外挿せず `MissingEarthOrientationData` / `MissingLeapSecondData` を返し metadata 記録（ISSUE-042 確定A2）。不確実性 0 の誤報を禁止。
- **NaN / inf**: `JulianDate2` 演算での伝播仕様を固定（ISSUE-004）。
- **異常入力**: 2 月 30 日・13 月・秒=60（非閏秒日）→ `TimeError::InvalidDate`（ISSUE-005）。
- **Espenak–Meeus 区間境界**: 多項式切替で連続性を確認（ISSUE-007）。

---

## 検証（受け入れ基準・基準値の出典）

> オラクル値は実装へコピーしない（conventions §11）。SOFA 同等手順 / Meeus worked example / IERS 公開値から独立に算出する。

- **暦 ↔ JD（L2）**:
  - 2000-01-01 12:00:00 → JD 2451545.0（J2000.0、Meeus Ch.7 記載値）。
  - 1987-01-27 00:00 → JD 2446822.5（Meeus 例）。1957-10-04.81（Sputnik、Meeus 例）。
  - 往復: 1900–2100 多数で `jd_to_gregorian(gregorian_to_jd(..))` が現代日付で ≤ 1 µs（≈1.16e-11 日）一致。整数日部は厳密一致。
  - 負年: 0 年 / −1 年の往復（Meeus/SOFA 一致）。
- **JulianDate2（L1/L8）**: 現代 JD に 1 µs を `add_seconds` し `diff_days` で復元（単一 f64 では不可能な対比テスト、ISSUE-004）。`add_seconds(86400) ≡ add_days(1)`。往復 `add_seconds(s).add_seconds(−s) ≈ jd`。
- **UTC/TAI/TT（L2）**: 2017-01-01 00:00 UTC で TAI−UTC = 37 s。`tt = utc + (TAI−UTC) + 32.184`。閏秒境界 2016-12-31 23:59:60 が表現可、前後 1 s で TAI 連続。決定論部は ≤ 1 µs、閏秒整数値・32.184 は厳密一致（オフバイワン無し）。オラクル = IERS `leap-seconds.list` + SOFA `iauUtctai`/`iauTaitt`。
- **UT1/ΔT（L2）**: 近傍年（2020）の ΔT が EOP 由来で ±0.1 s 内（恒等式経由）。Espenak–Meeus が 1900/2000/2050 の NASA 公表多項式値と一致（区間別）。不確実性 < 0.1 s（近傍）、~数秒（2100）。将来ほど不確実の単調性（プロパティ）。オラクル = IERS EOP C04 / Espenak–Meeus 公表値。
- **TDB（L2）**: `iauDtdb` 地心周期項の振幅が ±約 1.7 ms 規模。オラクル = SOFA `iauDtdb` 同等手順。
- **TimeData/TimeScales（ISSUE-042）**: `bundled()` と `from_path()` の同一性。valid_to 超過で Missing*Data + metadata 記録。不確実性帯伝播。

---

## 出典

- **Meeus, J. "Astronomical Algorithms", 2nd ed.** Chapter 7「Julian Day」式 (7.1)（暦→JD、グレゴリオ分岐）および JD→暦 逆変換手順（式番号は本文手順、個別番号は**要確認**）。
- **IAU SOFA**: `iauCal2jd`, `iauJd2cal`（暦↔JD・二引数 (date1,date2) 規約）、`iauDat`（閏秒）、`iauUtctai`/`iauTaiutc`（UTC↔TAI 境界）、`iauTaitt`/`iauTttai`（TT−TAI = 32.184 s）、`iauDtdb`（TDB−TT, Fairhead & Bretagnon 1990）。
- **IERS Conventions (2010) (TN 36)** §1（TT−TAI 定義）, ch.5（時刻系・地球回転）。
- **IERS EOP C04**（IERS Earth Orientation Center, `datacenter.iers.org`）: UT1−UTC, 極運動。系列 EOP 14 C04 / 20 C04（data-sources §3.1）。
- **IERS Bulletin C / IANA `leap-seconds.list`**（data-sources §3.2）: 閏秒積算。
- **Espenak, F. & Meeus, J. "Five Millennium Canon of Solar Eclipses: −1999 to +3000"（NASA TP-2006-214141）**: ΔT polynomial expressions（`eclipse.gsfc.nasa.gov/SEcat5/deltatpoly.html`、Morrison & Stephenson 2004 ベース）。不確実性式は**要確認**。
- **Fairhead, L. & Bretagnon, P. (1990)**: TDB−TT 地心モデル（SOFA `iauDtdb` 採用）。地上項 Moyer (1981) / Murray (1983)、基本引数 Simon et al. (1994)。
- **Dekker (1971) / Knuth**: two-sum / compensated summation（JulianDate2 桁分離加算）。
- プロジェクト内: conventions §4.1/§6/§11、accuracy.md §0/§2.3/§3.3/§5、numerical-policy §A1/§A3、ISSUE-004/005/006/007/042。
