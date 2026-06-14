# §10 観測者・楕円体（Observer & ellipsoid）

本節は `umbra-rs` の観測者位置変換（測地緯度→地心緯度、WGS84 楕円体での ρsinφ′/ρcosφ′、geodetic→ITRS/ECEF ベクトル）の数式・手順・符号規約を定める。
正本は `algorithms.md §0`（記号表）。本節はそれに厳密準拠し、本節固有の補助記号のみ追加定義する。
関連 Issue: ISSUE-010（WGS84・測地⇔地心緯度）/ 011（観測者→ITRS/ECEF）。

> 状態: ドラフト（Milestone 0）。一次資料で式番号を確認できない箇所は「**要確認**」を残し、推測で式番号を書かない。

---

## 目的と入出力（単位・時刻系・座標系・フレーム）

- **目的**: 観測者（測地緯度 φ・東経 λ・楕円体高 h）から、(1) 地心緯度 φ′ と扁平補正済み動径成分 ρsinφ′・ρcosφ′（ベッセル局地投影用）、(2) 地球固定 ITRS/ECEF 直交ベクトルを生成する。
- **入力**: 測地緯度 φ（`GeodeticLatitude`, rad）、東経 λ（`EastLongitude`, rad, 東経正）、楕円体高 h（`Meters`）。地球モデル `EarthModel::Wgs84`。
- **出力**: 地心緯度 φ′（`Radians`）、ρsinφ′・ρcosφ′（無次元、Re 単位の動径成分）、`Position<Itrs>`（km 既定、Re 無次元化版も供給）。
- **単位**: 角度 = rad（conventions §2）。半径・h = m（楕円体高、conventions §1/§4）。幾何ベクトル = km（conventions §1）。ベッセル投影は Re 無次元化（Re = WGS84 a, conventions §4.1）。
- **時刻系**: なし（地球固定。時刻依存の回転 ITRS→TIRS→CIRS→GCRS は §2 フレーム連鎖の責務）。
- **フレーム**: ITRS（≈WGS84 軸、右手系、Z = 自転軸、X = グリニッジ子午線方向、conventions §5）。**極運動（xp, yp）適用前**＝本節は ITRS（WGS84 軸）ベクトル生成まで。極運動補正は §2 ephemeris フレーム連鎖（ISSUE-007 が値供給）。
- **数値方針**: `√(1−e²sin²φ)` は全緯度で安定。極・赤道は `atan2` で安全に。単位変換（m↔km↔Re）は境界 1 箇所に集約（numerical-policy 横断 / conventions §11）。

---

## 記号（algorithms.md §0 を参照。本節固有の補助記号のみ定義）

algorithms.md §0「観測者（局地）」表（φ, λ, h, φ′, ρsinφ′/ρcosφ′）を正本とする。本節で補助的に用いる記号:

| 記号 | 意味 | 出典/備考 |
|---|---|---|
| a | WGS84 長半径 = Re = 6 378 137.0 m（定義値） | conventions §4.1。ベッセル無次元化基準 |
| 1/f | WGS84 逆扁平率 = 298.257223563（定義値） | conventions §4.1 |
| f | 扁平率 = 1/(1/f) | 導出（式 10.1） |
| b | 極半径 = a(1−f) | 導出（式 10.2） |
| e² | 第一離心率の二乗 = 2f − f² | 導出（式 10.3） |
| N | 卯酉線曲率半径（prime vertical radius of curvature） = a/√(1−e²sin²φ) | 式 10.6 |
| ρ | 地心動径（地表点）／ a（無次元）。ρ² = (ρsinφ′)² + (ρcosφ′)² | Meeus Ch.11 |
| (X, Y, Z) | ITRS/ECEF 直交成分（km、Re 版は無次元） | 式 10.7–10.9 |

> 本節は新規の独自記号を増やさない。a, b, f, e², N, ρ は標準楕円体記号で algorithms.md §0 観測者欄を補完するもの。記号表への追記候補は末尾「報告」に記載。

---

## 数式（番号付き。各式に出典）

### A. WGS84 楕円体パラメータ（ISSUE-010、conventions §4.1）

定義値（conventions §4.1 物理定数表、出典 NIMA TR8350.2 / NGA WGS84）:

```
a   = 6 378 137.0 m          # 定義値（= Re, ベッセル無次元化基準, conventions §4.1）
1/f = 298.257223563          # 定義値
```

導出値（magic number 禁止: a, 1/f から自前計算。conventions §11）:

**(10.1)** `f = 1 / 298.257223563`
**(10.2)** `b = a·(1 − f)`（極半径）
**(10.3)** `e² = 2f − f²`（標準楕円体関係、Meeus Ch.11）

- 出典: WGS84 定義（NIMA TR8350.2 / NGA）。e² = 2f − f² は標準楕円体関係（Meeus AA 2nd ed., Ch.11「The Earth's Globe」）。
- 注記: Meeus Ch.11 は `b/a = 0.99664719`（旧値）を例示。**WGS84 の f から自前計算した b/a を正本とし、Meeus 例値は検証参照に留める**（系統差をコメント、ISSUE-010）。

### B. 測地緯度 → 地心緯度・動径成分（ISSUE-010、Meeus Ch.11）

Meeus AA 2nd ed., Ch.11「The Earth's Globe」の補助角 u と動径成分。h は楕円体高（m）。

**(10.4)** 補助角（reduced latitude）:
```
u = atan( (b/a)·tan φ )      # 極で発散しないよう atan2 表現も可
```

**(10.5)** 扁平補正済み動径成分（標高 h 込み、Meeus Ch.11 式 11.x）:
```
ρ sinφ′ = (b/a)·sin u + (h/a)·sin φ
ρ cosφ′ = cos u       + (h/a)·cos φ
```

**(10.5a)** 地心緯度:
```
φ′ = atan2( ρ sinφ′ , ρ cosφ′ )
```

- 出典: **Meeus AA 2nd ed., Ch.11「The Earth's Globe」**（補助角 u、`ρ·sinφ′`, `ρ·cosφ′`、標高項 h/a を含む）。地心緯度 φ′ = atan2(ρsinφ′, ρcosφ′)。
- **要確認**: ISSUE-010 は式 11.1〜11.4 を参照と記すが、Meeus 本文での `u` / `ρsinφ′` / `ρcosφ′` の**個別の式番号**（11.1 か 11.2 か等）は一次資料（書籍）で再確認が必要。式番号は推測で確定しない。Meeus の標高記号は `H`（本プロジェクトは h、conventions §3/§4）。
- 注記: Meeus の原式は標高項を `(H/a)·sinφ` 等としており、a は楕円体長半径。本プロジェクトでは a = Re = WGS84 a（conventions §4.1）。

### C. 測地座標 → ITRS/ECEF 直交ベクトル（ISSUE-011、IERS Conventions 2010 §4）

標準楕円体公式（卯酉線曲率半径 N 経由）。λ は東経（東経正、conventions §3）。

**(10.6)** 卯酉線曲率半径:
```
N = a / √(1 − e²·sin²φ)
```

**(10.7)** `X = (N + h)·cos φ·cos λ`
**(10.8)** `Y = (N + h)·cos φ·sin λ`
**(10.9)** `Z = (N·(1 − e²) + h)·sin φ`

単位は m → km へ境界変換（conventions §1）。

- 出典: **IERS Conventions 2010 (TN 36) §4**、および測地学標準（Torge "Geodesy" / Hofmann-Wellenhof "GNSS"）。
- **要確認**: ITRS の X 軸（経度 0 = グリニッジ子午線）方向定義を IERS と一致させる。極運動適用前の TIRS との差は §2/ISSUE-007 で扱う（本節は WGS84 軸の ITRS まで）。IERS Conventions §4 の該当式番号は一次資料で確認のこと（推測しない）。

**(10.10)** 等価表現（§B との整合、検証用）:
```
Z          = a · ρ sinφ′
√(X² + Y²) = a · ρ cosφ′
```

§B の Meeus 形（ρsinφ′/ρcosφ′）と §C の N 形は等価。**実装はどちらか一方を正本とし、もう一方を検証に使う**（ISSUE-011）。

### D. Re 無次元化（ベッセル投影用、conventions §4.1）

**(10.11)**:
```
(X, Y, Z)_Re = (X, Y, Z)_m / Re ,   Re = a = 6 378 137.0 m
```

ベッセル要素 x, y, l1, l2 は Re 単位（conventions §1/§4.1）。局地投影（§9 観測者基本面座標 ξ, η, ζ）はこの Re 版を消費する（algorithms.md §0 観測者欄 / §10 索引）。

---

## 手順（実装順・数値注意）

1. **楕円体パラメータ**: `EarthModel::Wgs84` から a, 1/f（定義値）を供給し、f, b, e² を (10.1)–(10.3) で導出。すべて定数化＋出典コメント（magic number 禁止、conventions §11）。Re = a を `reference_equatorial_radius()` で供給。
2. **地心緯度・動径成分**（局地投影用）: (10.4)(10.5)(10.5a)。`tan φ` の発散を避け `atan2` で表現。h は楕円体高。
3. **ITRS ベクトル**: (10.6)–(10.9)。λ = 東経（西経入力は `EastLongitude::from_signed_degrees` で境界吸収、内部に西経正を持ち込まない、conventions §3）。m → km は境界 1 箇所で変換。
4. **Re 版**: (10.11)。ベッセル投影用。km 版と混在させない（ISSUE-011）。
5. **正本/検証の分離**: §B（Meeus ρsinφ′/ρcosφ′）と §C（IERS N 形）の一方を正本、(10.10) で相互検証。
6. **責務分離**: 極運動（xp, yp）適用は行わない（§2 ephemeris フレーム連鎖の責務、ISSUE-011 doc 明記）。地心/測地の取り違え禁止（変数名 `geodetic_lat`/`geocentric_lat` で区別、conventions §3）。

---

## 境界・特異・異常系

- **赤道（φ = 0）**: 地心緯度 φ′ = 0。Z = h·0 ＝ 0（h=0 で）。`(a, 0, 0)`（λ=0, h=0）。
- **極（φ = ±90°）**: 地心緯度 = ±90°。`cos φ = 0` を `atan2`/直接式で安全に。北極 h=0 → `(0, 0, b)`。`tan φ` 発散を (10.4) の atan2 表現で回避。
- **日付変更線（λ = ±180°）**: 東経正で連続。λ→λ+2π でベクトル不変（プロパティ）。
- **緯度符号反転**: φ→−φ で Z 符号反転（プロパティ）。
- **高標高/負標高**: h = 8000 m（高山）、h = −400 m（死海）等で破綻しない。
- **西経入力**: `EastLongitude::from_signed_degrees(−120)` が東経 240° 相当ベクトルに（conventions §3 西経吸収）。
- **単位混在禁止**: m と km、km と Re を混ぜない。変換は境界に集約（conventions §11）。
- **Meeus 旧 b/a との系統差**: 検証で Meeus 例値（b/a=0.99664719）を使う際は WGS84 自前 b/a との系統差をコメント（ISSUE-010）。

---

## 検証（受け入れ基準・基準値の出典）

> オラクル値は実装へコピーしない（conventions §11）。WGS84 定義値・Meeus worked example・独立公式・既知測地点 ECEF から独立に算出する。

- **楕円体パラメータ（L1）**:
  - a = 6 378 137.0 m, 1/f = 298.257223563 が WGS84 定義値と**厳密一致**。
  - e² = 2f − f², b = a(1−f) の導出一致。
- **地心緯度・動径（L1, L6 前提）**:
  - 赤道（φ=0）→ φ′ = 0。極（φ=±90°）→ φ′ = ±90°。
  - Meeus Ch.11 Palomar 例（φ=33°21′22″, h=1706 m）の ρsinφ′, ρcosφ′ 近似一致（旧 b/a 系統差をコメント）。オラクル = Meeus AA 2nd ed. Ch.11 worked example。
  - 中緯度（φ=45°, h=0）: φ′ < φ（差 最大 ~11.5′ 付近）。
  - プロパティ（L8）: |φ′| ≤ |φ|（h=0 で常に成立、符号同じ）。
  - 許容: 地心緯度 ≤ 1e-9 rad（≈0.2 mas、地表 ~6.4 mm。accuracy.md §2.3 sub-m に余裕）。ρ 相対 ≤ 数 ULP。
- **ITRS/ECEF（L1, L6 前提）**:
  - 赤道・λ=0・h=0 → `(a, 0, 0)` km。北極・h=0 → `(0, 0, b)` km。λ=90°E・赤道 → `(0, a, 0)`（東経正符号確認）。
  - 既知地点（api-draft 例 岡山 34.507°N, 133.508°E, 10 m）の ECEF を独立公式で照合。
  - ノルム `√(X²+Y²+Z²)` ≈ 地心距離（§B の ρ·a と一致、(10.10) 検証）。
  - 許容: ECEF 位置 sub-m（≤ 1 m、目標 ≪ 0.1 m。accuracy.md §2.3）。既知点 ≤ 0.01 m 目標（丸めのみ）。Re 無次元化 相対 ≤ 数 ULP。
  - 境界: 極（±90°）、日付変更線（±180°）、高/負標高。西経入力吸収。

---

## 出典

- **conventions §4.1**（物理定数表）: a = Re = 6 378 137.0 m, 1/f = 298.257223563（WGS84 定義値、NIMA TR8350.2 / NGA）。
- **Meeus, J. "Astronomical Algorithms", 2nd ed.** Chapter 11「The Earth's Globe」: 補助角 u、ρ·sinφ′ / ρ·cosφ′（標高項含む）、e² = 2f − f²。地心緯度 φ′ = atan2(ρsinφ′, ρcosφ′)。個別の式番号（11.1〜11.4 の対応）は**要確認**。
- **IERS Conventions (2010) (TN 36) §4**: 測地座標 → ECEF（N = a/√(1−e²sin²φ) 経由）。該当式番号は**要確認**。
- **IAU SOFA `iauGd2gc`**（geodetic → geocentric Cartesian）/ `iauGc2gd`: ECEF 変換の照合基準。
- 測地学標準: Torge "Geodesy" / Hofmann-Wellenhof "GNSS"（楕円体 → ECEF 標準公式）。
- プロジェクト内: conventions §1/§3/§4/§5/§11、accuracy.md §2.3、ISSUE-010/011、api-draft §1.4/§1.5。
