# §9 局地条件（Local circumstances, Step10）

本節は `umbra-rs` の**観測地点における日食状況**（観測者の基本面射影 ξ,η,ζ・局地接触 C1–C4・最大食・食分 magnitude・食面積 obscuration・太陽高度方位・可視性）の数式・手順・符号規約を定める。各 solver が時刻関数として呼ぶ最下層プリミティブ（観測者射影）から、可視性判定までを通す。
正本は `algorithms.md §0`（記号表）。本節はそれに厳密準拠し、本節固有の補助記号のみ追加定義する。
関連 Issue: ISSUE-024（観測者射影）/ ISSUE-025（C1–C4 接触）/ ISSUE-026（最大食）/ ISSUE-027（食分・食面積）/ ISSUE-028（高度方位・可視性）。供給源は ISSUE-037（直接評価, Standard 局地既定）/ ISSUE-022（多項式）。

> 状態: ドラフト（Milestone 0）。**局地時角の符号（H = μ − λ か μ + λ）は一次資料で式番号まで未確認のため「要確認」**として東経正の適応式を明示し、交差検証テストを指定する（後述 §A 注記・検証）。一次資料で式番号を確認できない箇所は「**要確認**」を残し、推測で式番号を書かない。

---

## 目的と入出力（単位・時刻系・座標系・フレーム）

- **目的**: 観測者（測地緯度 φ・東経 λ・楕円体高 h, ITRS 化は §10/ISSUE-011）を各時刻の月影軸を z 軸とするベッセル基本面へ射影し、局地接触・最大食・食分・食面積・高度方位・可視性を確定する（conventions §5/§7/§8/§9, architecture §6/§7）。
- **入力**: `Observer`（ITRS の Re 無次元化版, §10/ISSUE-011）、`BesselianSource`（任意 TT で x, y, d, μ, l1, l2, tan f1, tan f2 を供給。x, y, l1, l2 = Re, d, μ = rad）、太陽見かけ地心 α,δ（高度方位用, §3/ISSUE-015）、`TimeScales`（TT↔UTC）、`RefractionModel { None, Standard }`、`EngineConfig`（root_tolerance, k 値選択 conventions §9）。
- **出力**: `LocalCircumstances`（C1–C4 の `LocalContact`（TT+UTC, alt/az/PA/visible）、最大食 `LocalMaximum`（TT+UTC, m_min[Re], magnitude, obscuration）、`Visibility` 6 値）。
- **単位**: 角度 = rad 内部（公開は度も, conventions §2）。ξ, η, ζ, u, v, L1, L2, m = Re 無次元（conventions §1/§5）。微分の時間単位 = **SI 秒基準で固定**（Re/SI秒, conventions §1。Meeus の分単位は内部に持ち込まない, 境界で換算）。
- **時刻系**: 瞬時要素は **TT 基準**（conventions §6, accuracy §0(a)）。**ただし局地接触時刻は μ（= ERA(UT1)−α）経由で δUT1 に依存する**（accuracy §0(a) 脚注・§2.1L, ISSUE-021 I1）。各接触・最大食は **UTC と TT の両方**を返す（accuracy §0）。
- **フレーム**: FundamentalPlane（影軸 z、x̂=東、ŷ=天の北の射影, conventions §5, ISSUE-020）。高度方位は地平座標（観測者局所）。

---

## 記号（algorithms.md §0 を参照。本節固有の補助記号のみ定義）

algorithms.md §0「観測者（局地）」表（φ, λ, h, φ′, ρsinφ′/ρcosφ′, θ, ξ, η, ζ, u, v, L1, L2, m）と「食の量」表（magnitude, obscuration）を正本とする。本節で補助的に用いる記号:

| 記号 | 意味 | 出典/備考 |
|---|---|---|
| `θ` | 観測者における影軸の局地時角（= §0 の θ）。**θ = μ − λ（東経正, 要確認）** | 式 9.1・ISSUE-024 |
| `ξ', η'` | ξ, η の時間微分（Re/SI秒）。接触感度・継続時間に使用 | 式 9.5・ISSUE-024 |
| `u', v'` | u, v の時間微分（Re/SI秒）。u'=x'−ξ', v'=y'−η' | 式 9.9 |
| `n²` | 相対速度二乗 `u'² + v'²`（Re²/SI秒²）。接触補正・継続時間 | 式 9.10・Meeus Ch.54 |
| `g_C(t)` | 接触の求根対象 `m²(t) − L(t)²`（外接 L=L1 / 内接 L=L2） | 式 9.11 |
| `s_sun`, `s_moon` | 太陽/月 見かけ視半径（同単位, §0 天体位置表）。食面積に使用 | 式 9.16・ISSUE-027 |
| `Δ_sep` | 食面積用の太陽-月見かけ中心間距離（同単位） | 式 9.16 |
| `a_geom`, `a_app` | 幾何学的高度 / 大気差補正後高度（rad） | 式 9.19/9.21 |
| `A` | 太陽方位角（**北 0・東回り**, rad） | 式 9.20・conventions §7 |
| `PA` | 接触点の位置角（**天の北 0・東回り**, rad） | 式 9.22・conventions §7 |
| `ΔR` | 大気差量（rad） | 式 9.21・Meeus Ch.16 |

> ξ, η, ζ, u, v, L1, L2, m, magnitude, obscuration は §0 正本。θ は §0 にあるが**符号適応（μ−λ）が要確認**のため本節で明示。記号表への追記候補（θ の符号確定・ξ'/η'・n²）は末尾「報告」に記載。

---

## 数式（番号付き・各式に出典）

### A. 観測者の基本面射影 ξ, η, ζ（扁平補正込み）

観測者地心成分 `ρ cosφ′, ρ sinφ′`（扁平補正済み, §10/ISSUE-010・011）と影軸赤緯 d・局地時角 θ から基本面座標を作る。

**(9.1)** 局地時角（東経正の適応式, **要確認**）:
```
θ = μ − λ            （μ = 影軸の見かけグリニッジ時角 = ERA(UT1) − α_axis, §2/§0。λ = 観測者東経）
```

**(9.2)** 観測者の基本面座標（標準形）:
```
ξ = ρ cosφ′ · sin θ
η = ρ sinφ′ · cos d − ρ cosφ′ · sin d · cos θ
ζ = ρ sinφ′ · sin d + ρ cosφ′ · cos d · cos θ
```

- 出典: **Explanatory Supplement to the Astronomical Almanac (3rd ed.), §11「Eclipses」**（局地予報式: 観測者地心座標 ρcosφ′/ρsinφ′ と時角から基本面座標）/ **Meeus, *Astronomical Algorithms* (2nd ed.), Ch.54「Solar Eclipses」（式 54.1〜）** の u, v（=ξ, η）構成式。個別式番号は**要確認**（一次資料＝書籍で確認, 推測しない）。
- **扁平の効き**: ξ, η, ζ には地心緯度成分 `ρ cosφ′, ρ sinφ′` を通して **WGS84 扁平**が直接効く（測地緯度 φ → 地心緯度 φ′ の差・ρ の緯度依存, §10/ISSUE-010・011）。地心/測地の取り違えや球近似は ξ,η,ζ → 接触時刻・食分の系統誤差になる。`ρ cosφ′, ρ sinφ′` は §10 の WGS84 値を用い、扁平を含む点を実装コメントに明記（ISSUE-024）。
- **ζ の意味**: 観測者が基本面より月側（前方）か地球内側かの符号判定に使う（接触幾何・高度方位の前提, ISSUE-024）。
- **要確認（符号正本化, ISSUE-024/028）**: **局地時角を θ = μ − λ（東経正）と定める**（ξ = ρcosφ′·sin θ の sin の引数）。文献により H = μ + λ_east の符号定義が異なるため、**Explanatory Supplement §11 / Meeus Ch.54 のどちらを正本にするか決め、もう一方で交差検証**する（検証 §「時角符号×方位符号 交差検証」）。確証が取れるまで両流儀を併記し、東経正・conventions §3 に合わせて θ = μ − λ を採用、ISSUE-028 の方位符号と整合させる。実装コメントに採用式の章・式番号を転記。
- 正規化: θ は `[-π, π)` の signed 正規化（μ−λ の折返しで微分が壊れないよう連続化, conventions §2）。φ′, d 全域で sin/cos θ は連続。

### B. 相対位置 u, v と影半径 L1, L2

**(9.3)** 影軸と観測者の基本面相対位置（algorithms.md §0）:
```
u = x − ξ
v = y − η
```

**(9.4)** 中心間距離の二乗（numerical-policy §A5: 最小化対象は m²）:
```
m² = u² + v²            [Re²]
m  = √(u² + v²)         [Re]
```

**(9.5)** ζ 面での錐半径（観測者高さ＝基本面からの距離 ζ で補正, algorithms.md §0）:
```
L1 = l1 − ζ · tan f1            （半影縁・外接 C1/C4 用）
L2 = l2 − ζ · tan f2            （本影/反本影縁・内接 C2/C3 用）
```

- 出典: Explanatory Supplement §11（`(ξ−x)² + (η−y)² = L²`、`L = l − ζ·tan f` の高さ補正）/ Meeus Ch.54（u=ξ−x, v=η−y 構成。本実装は §0 の u=x−ξ に符号を揃え、m²=u²+v² は符号不変）。式番号は**要確認**。
- 注記: §0 は `u = x − ξ, v = y − η`。Meeus の `u = ξ − x` とは符号が逆だが、**m² = u²+v² は不変**で接触・最大食の求解結果は一致（実装は §0 規約に統一, ISSUE-025）。
- L2 符号: 金環で l2 < 0（algorithms.md §0 / ISSUE-021）。接触の内接条件は `m = |L2|`（conventions §8）。

### C. 局地接触 C1/C4（外接）と C2/C3（内接）

接触条件（conventions §8）: C1/C4 は外接 `m = L1`、C2/C3 は内接 `m = |L2|`。求根は二乗形で扱い acos 不要。

**(9.6)** 外接（C1: 部分食開始 / C4: 部分食終了）:
```
C1, C4 :  m² = L1²      ⇔   (x−ξ)² + (y−η)² = (l1 − ζ tan f1)²
```

**(9.7)** 内接（C2: 中心食開始 / C3: 中心食終了）:
```
C2, C3 :  m² = L2²      ⇔   (x−ξ)² + (y−η)² = (l2 − ζ tan f2)²
```

**(9.11)** 求根対象（差の連続関数, Brent）:
```
g_C(t) = m²(t) − L(t)²       （外接 L=L1 / 内接 L=L2）
```

- 出典: Explanatory Supplement §11（`(ξ−x)²+(η−y)² = L²`）/ Meeus Ch.54（u, v, n², 接触時刻補正の式 54.x）。conventions §8（外接/内接の定義）。式番号は**要確認**。
- **手法（粗走査→Brent, conventions §11 / numerical-policy §A5）**: 探索窓（全球接触 §8 ±マージン）を一定刻みで粗走査し `g_C(t)` を評価 → 符号変化区間をブラケット → **Brent 求根**（無条件 Newton 禁止）。Meeus の線形補正反復は**初期窓見積りのみ**に使用、確定は Brent。root_tolerance ≤ 0.01 s（接触 ±2 s の 1/10 以下, accuracy §2.1, numerical-policy §A5）。
- 接触種別の割当て: `g_C` の符号遷移方向（外側→内側 = C1/C2、内側→外側 = C3/C4）で機械的に。外接（L1）と内接（L2）は別 solver パス（L が異なる, ISSUE-025）。
- **中心食条件**: その地点で `m_min < |L2|`（本影が観測者を覆う）なら C2/C3 が存在、さもなくば **C2/C3 = None（部分食地点）**（api-draft §6 None 設計, ISSUE-025）。

### D. 最大食（dm/dt = 0 の求根）

最大食 = m が最小の瞬間。中心線尖点での √ 不連続を避けるため最小化対象は m²（D2, accuracy §2.1 / numerical-policy §A5）。

**(9.8)** 微分:
```
ξ', η' : ξ, η の時間微分（解析 or 中心差分＋Richardson, numerical-policy §A2(3)）。単位 Re/SI秒
```

**(9.9)** u, v の微分:
```
u' = x' − ξ'
v' = y' − η'
```

**(9.10)** 相対速度二乗（Meeus Ch.54）:
```
n² = u'² + v'²          [Re²/SI秒²]
```

**(9.12)** 最大食条件（dm/dt = 0 ⇔ d(m²)/dt = 0, Brent 求根）:
```
t_max : d(m²)/dt = 2(u·u' + v·v') = 0      ⇔   u·u' + v·v' = 0
m_min = m(t_max) = √(u(t_max)² + v(t_max)²)    [Re]
```

- 出典: Explanatory Supplement §11（投影面上の太陽中心-月中心距離最小）/ Meeus Ch.54（u=ξ−x, v=η−y, m=√(u²+v²) 最小、線形補正 `Δt = −(u·u'+v·v')/n²`）。式番号は**要確認**。
- **手法（D2 正式手法）**: `u·u'+v·v' = 0` を **Brent 求根**（ISSUE-008）。Meeus の線形補正 `Δt = −(u·u'+v·v')/n²` は**粗ブラケット（初期推定）のみ**に降格、距離最小化（黄金分割）も粗ブラケット併用（無条件 Newton 回避, conventions §11 / accuracy §2.1）。
- **皆既帯平底の退行扱い**: 皆既帯中心付近では `m < |L2|` の**平底**（dm/dt≈0 が区間で成立）があり得て**最小は一意でない**。平底区間の代表点（中央 or 全球最大食に最も近い点）を規約で定義し、解の非一意を許容（accuracy §2.1, numerical-policy §A5, ISSUE-026）。
- 最大食は**部分食地点でも必ず存在**（`maximum` は非 Option, api-draft §3.4）。可視域外でも幾何最大は計算し、可視性は §G が別判定。最大食は接触の内側（`c1 < t_max < c4`, ISSUE-026）。

### E. 食分 magnitude

**(9.13)** 食分（algorithms.md §0 / Explanatory Supplement §11 / Meeus Ch.54）:
```
magnitude = (L1 − m) / (L1 + L2)
```

- 出典: Explanatory Supplement §11 / Meeus Ch.54（外接縁から内接縁への食い込み割合）。**要式番号確認**（一次資料＝書籍で章・式番号を確定, 推測しない）。
- 範囲: 部分食 0 < magnitude < 1、皆既で > 1、金環で m≈|L2| 付近。`EclipseMagnitude` は 1 超を許容（皆既）。`m ≥ L1`（離隔）なら 0 にクランプ（食なし, ISSUE-027）。
- 最大食点では (9.13) を t_max の (m_min, L1, L2) で評価（§D 連携, ISSUE-026/027）。

### F. 食面積 obscuration（2 円交差面積）

obscuration = 月円が太陽円を覆う面積比。太陽見かけ半径 R = s_sun、月見かけ半径 r = s_moon、中心離隔 Δ_sep（同単位）。

**(9.16)** 入力（視半径平面・同単位, ISSUE-027）:
```
R = s_sun ,  r = s_moon ,  Δ_sep = 太陽-月 見かけ中心間距離
```

**(9.17)** 2 円の重なり面積（lens area, 部分重なり時のみ）:
```
A = R²·acos((Δ_sep² + R² − r²)/(2·Δ_sep·R))
  + r²·acos((Δ_sep² + r² − R²)/(2·Δ_sep·r))
  − ½·√( (−Δ_sep+r+R)(Δ_sep+r−R)(Δ_sep−r+R)(Δ_sep+r+R) )
```

**(9.18)** 食面積比:
```
obscuration = A_overlap / (π · R²)        （= A_overlap / 太陽円面積）
```

- 出典: 円-円交差面積の標準幾何式（**Weisstein, MathWorld "Circle-Circle Intersection"**。標準公式のため出典は公式名で明記, 章・式番号ではない）。日食文脈の食面積定義は Explanatory Supplement §11 / Meeus Ch.54。
- **acos クランプ必須**: 引数を `[-1, 1]` にクランプ（`x.clamp(-1.0, 1.0)`, 丸めで ±1 超え→NaN 回避, accuracy §2.2 / numerical-policy §A5）。判別式 `√(...)` の負値（丸め）も 0 クランプ。結果も `[0, 1]` クランプ。
- **5 境界の明示処理**（accuracy §2.2 / ISSUE-027, 0 除算 `2·Δ_sep·R`, `2·Δ_sep·r`, `Δ_sep=0` 回避）:
  1. **離隔** `Δ_sep ≥ R + r` → overlap = 0 → obscuration = 0。
  2. **内包** `Δ_sep ≤ |R − r|` → overlap = π·min(R,r)²。月が大（皆既近傍）→ obscuration = 1、太陽が大（金環）→ `r²/R²`。
  3. **外接** `Δ_sep = R + r` → overlap = 0（離隔境界）。
  4. **内接** `Δ_sep = |R − r|` → overlap = π·min(R,r)²（内包境界）。
  5. **同半径** `R = r` → lens 式の `2·Δ_sep·R = 2·Δ_sep·r`、判別式は対称形。0 除算は Δ_sep>0 で回避、Δ_sep=0 は内包分岐へ。
- 部分食地点で C2/C3 = None でも食分・食面積は定義される（最大時点 m で評価, ISSUE-025/026/027 整合）。
- 単位系の分離: **食分は基本面 Re（m, L1, L2）、食面積は視半径平面（s_sun, s_moon, Δ_sep）**。混在禁止（conventions §1, ISSUE-027）。

### G. 太陽高度・方位・可視性

赤道座標（太陽見かけ α, δ, §3/ISSUE-015）→ 地平座標。局地時角は §A と同一符号規約（CIO 系時角, 分点 GST 禁止, D4）。

**(9.19)** 幾何学的高度（球面天文標準, Meeus Ch.13 式 13.6）:
```
sin a_geom = sin φ · sin δ + cos φ · cos δ · cos H        （H = μ_sun − λ。要確認: §A と同符号）
```

**(9.20)** 方位角（**北 0・東回り**, conventions §7）:
```
tan A = sin H / (cos H · sin φ − tan δ · cos φ)        を atan2 で象限確定
A は南基準（Meeus Ch.13 式 13.5）→ 北 0・東回りへ 180° オフセット規約変換（要確認: 符号を実装で固定）
```

**(9.21)** 大気差補正後高度（`RefractionModel::Standard` のみ, Meeus Ch.16 式 16.3/16.4, Bennett 1982 / Sæmundsson）:
```
a_app = a_geom + ΔR        （ΔR = 標準大気差量。補正前後の両方を返す, conventions §7）
```

**(9.22)** 位置角（**天の北 0・東回り**, conventions §7）:
```
PA = （接触点の天の北からの角, 東回り）       出典: Explanatory Supplement §11 / Meeus Ch.54
```

- 出典: **Explanatory Supplement §7（座標系）/ Meeus AA 2nd ed. Ch.13「Transformation of Coordinates」式 13.5/13.6**（赤道→地平）。Meeus は**方位を南基準**で定義 → **北 0・東回りへ変換**（conventions §7）。大気差は **Meeus Ch.16 式 16.3/16.4**。SOFA 参照（移植しない）: `iauHd2ae`（hour angle/dec → az/el）。**時角・恒星時は ERA(iauEra00)経由 CIO ベース・分点 GST（`iauGst06a`）禁止**（D4, ISSUE-039 供給, ISSUE-028）。
- **既定は幾何学的高度**。大気差は `RefractionModel` で分離、補正前後を両方返す（conventions §7, accuracy §6 地平線付近は非保証）。
- **可視性 6 値**（`Visibility`, api-draft §3.4 / ISSUE-028 判定木）:
  - **NotVisible**: その地点で食域外（接触なし）。
  - **BelowHorizon**: 最大食時も太陽が地平下。
  - **SunriseEclipse**: 食の進行中に日の出（C1〜最大の一部が地平下で以降可視）。
  - **SunsetEclipse**: 食の進行中に日没。
  - **PartialVisible**: 一部の接触が地平下。
  - **FullyVisible**: C1〜C4 すべて地平上。
  - 日の出/日没境界の高度閾値（幾何 0 か 太陽縁＋大気差 −0.83°）を **1 箇所で固定**（既定は幾何学的高度, conventions §7, ISSUE-028）。

---

## 手順（実装順・数値注意=numerical-policy 参照）

1. **観測者射影 ξ, η, ζ**（最下層, ISSUE-024）: §10 から `ρ cosφ′, ρ sinφ′`（標高込み, WGS84 扁平）→ θ = μ − λ（要確認, 連続化）→ (9.2)。各時刻 t の純関数として solver が反復評価。
2. **u, v, m², L1, L2**: (9.3)–(9.5)。m² 基準で扱う（numerical-policy §A5）。
3. **局地接触 C1–C4**（ISSUE-025）: 探索窓を粗走査し `g_C = m² − L²` の符号変化をブラケット → Brent（外接 L1 / 内接 L2 を別パス）。符号遷移方向で C1/C2/C3/C4 割当て。`m_min < |L2|` で C2/C3 存在判定。
4. **最大食**（ISSUE-026）: `m²` 最小付近を粗ブラケット → `u·u'+v·v'=0` を Brent 求根（D2）。平底は代表点規約。m_min を取得。
5. **食分・食面積**（ISSUE-027）: 最大時点（or 各時点）の (m, L1, L2) で magnitude (9.13)、(s_sun, s_moon, Δ_sep) で obscuration (9.17/9.18)。acos/判別式クランプ・5 境界分岐。
6. **高度方位・可視性**（ISSUE-028）: 各接触/最大時点で (9.19)–(9.22) → 北 0 東回り変換 → 可視性 6 値判定木。
7. **時刻系**: 各接触・最大食を **TT と UTC の両方**で返す（accuracy §0）。**局地接触時刻は μ 経由で δUT1 依存**（accuracy §2.1L, ISSUE-021 I1）。将来日食は `delta_t_uncertainty_seconds` を metadata に。

**数値注意**（numerical-policy）:
- 求根は **無条件 Newton 禁止・Brent（ブラケット必須）**（conventions §11 / §A5）。粗走査刻みは掠め食（grazing, `m_min≈L1`）でも符号変化を逃さない細かさに（偽陰性ゼロ, architecture §3）。
- 最小化対象は m²（√ の中心線尖点の微分特異を回避, D2 / §A5）。
- ξ, η の ±π 折返し除去（θ = μ−λ の連続化）後に m², dm/dt を扱う（conventions §2）。
- 微分 ξ', η', x', y' は中心差分＋Richardson か解析微分（固定刻み禁止, §A2(3)）。**時間単位は SI 秒で 1 箇所に固定**（Re/SI秒, conventions §1）。
- acos/asin は `[-1, 1]` クランプ、判別式は 0 クランプ（§A5 / accuracy §2.2）。
- root_tolerance = 0.01 s 既定（接触 ±2 s の 1/10 以下, §A5）。
- magic number 禁止（k 値・大気差係数・閾値は出典付き定数, conventions §11）。西経正の内部持ち込み禁止（conventions §3）。地心/測地・km/Re 混在禁止。

---

## 境界・特異・異常系

- **部分食地点**: `m_min ≥ |L2|` → C2/C3 = None、C1/C4 は存在。食分・食面積は最大時点 m で定義（ISSUE-025/027）。
- **中心線上**: ξ=η≈0 → m_min≈0。皆既で magnitude ≥ 1、金環で m_min≈|L2|・magnitude≈1。順序 `c1 < c2 < t_max < c3 < c4`（プロパティ, ISSUE-025）。
- **掠め食（grazing, `m_min≈L1`）**: 粗走査刻みを十分細かく（見落とし防止, ISSUE-025）。
- **皆既帯平底**（`m < |L2|` 平坦, dm/dt≈0 区間）: 最大食時刻は代表点規約で決定的に（§D, accuracy §2.1, ISSUE-026）。
- **食面積 5 境界**: 離隔/内包/外接/内接/同半径で 0 除算（`2·Δ_sep·R`, `Δ_sep=0`）・NaN（acos, 判別式）を場合分けで回避（§F, ISSUE-027）。
- **可視性 6 値**: 最大時も地平下 → BelowHorizon、食域外 → NotVisible、日の出/日没中 → Sunrise/SunsetEclipse、一部地平下 → PartialVisible、全可視 → FullyVisible。日の出日没閾値を 1 箇所固定（ISSUE-028）。
- **極・日付変更線・天頂**: sin/cos(θ)・atan2 で連続。λ→λ+2π / μ→μ+2π で不変（プロパティ, ISSUE-024/028）。方位は atan2 で象限保持。
- **西経入力**: `EastLongitude::from_signed_degrees(−100)` が東経 260° 相当（conventions §3 西経吸収, ISSUE-024）。
- **標高差**: h=0 と h=4000 m で ζ が標高分変化（基本面前後位置, ISSUE-024）。
- **ブラケット不成立**: `RootNotBracketed` / 窓が取れない（食なし）→ 全接触 None で上位が可視性 NotVisible 判定（ISSUE-025）。

---

## 検証（受け入れ基準・基準値の出典）

> オラクル値は実装へコピーしない（conventions §11）。NASA/USNO 局地値・ゴールデン20・MockEphemeris・独立実装から算出。accuracy テストレベル **L6（局地条件）**（食面積は L1 純数学も）。

- **観測者射影（ISSUE-024）**: MockEphemeris で完全中心配置の影軸直下 → ξ=η≈0、横ずれ点で √(ξ²+η²) が独立計算と一致。既知地点（api-draft 例 岡山 34.507°N, 133.508°E, 10 m）で手計算 (ξ,η,ζ) と照合（オラクル = Explanatory Supplement §11 式の独立実装 or 公開ワークシート, 出典・取得日明記）。地点分類（中心線上/付近/北南限/部分食域/可視域外/標高差）。微分 ξ', η' を中心差分と一致。
- **接触 C1–C4（ISSUE-025）**: NASA 5千年カタログ / USNO の地点別接触時刻（第二義・整合, k・ΔT 慣習を Espenak へ揃える, accuracy §3.1）。部分食地点 c2==None && c3==None・c1/c4 存在。中心食地点 `c1 < c2 < max < c3 < c4`（順序 L8）。掠め食の刻み感度（偽陰性ガード）。`RootNotBracketed` 異常系。
- **最大食（ISSUE-026）**: 中心線上 m_min≈0・食分 ≥1（皆既）/≈1（金環）。**皆既帯平底 fixture**（dm/dt=0 求根が平底で破綻せず収束・代表点で決定的・食分安定）。最小性 `m(max) ≤ m(max±δ)`（平底では等号）。max は c1〜c4 の内側。
- **食分・食面積（ISSUE-027）**: 純幾何（L1, 解析解）— `Δ_sep=0,R=r`→overlap=πR²(obscuration=1)、`Δ_sep≥R+r`→0、`Δ_sep=|R−r|`→π·min²、acos 引数 1.0000001→クランプで NaN 出さない。5 境界（離隔/内包/外接/内接/同半径）で有限値。obscuration は Δ_sep 単調減少で単調増加・`[0,1]`。
- **高度方位・可視性（ISSUE-028）**: Meeus Ch.13 例題と alt/az 照合（北 0 東回り変換後）。**方位規約**: 真南正中 A=180°、真東 A=90°、真西 A=270°。**時角符号×方位符号 交差検証（必須）**: θ=μ−λ（東経正・正本）と方位（北 0 東回り）を一方を正本にもう一方で交差検証（正中 H=0 で A=180°、東地平 H<0 と西地平 H>0 で alt/az 符号整合, ISSUE-024 の時角符号と一致, CIO 系時角で固定）。大気差 Meeus Ch.16 例題（高度 0/5/45°）。可視性 6 値すべてを網羅する fixture。
- **ゴールデン20**（accuracy §3.4）の局地部分。**UTC/TT 両返し**・将来日食で `delta_t_uncertainty` が metadata（accuracy §0, **局地は μ 経由 δUT1 依存**, §2.1L）。

**許容誤差**（accuracy §2 / §4 層分解）:
- (ξ, η, ζ) 位置: sub-m 相当（Re 換算 ≤ 1.6e-7 Re, 目標 ≪ 0.1 m, accuracy §2.3）。既知点一致 ≤ 1e-9 Re（丸めのみ）。微分の数値一致 ≤ 1e-6（ISSUE-024）。
- C1–C4・最大食時刻（TT 基準・幾何相対）±2 s（接触）/ ±1〜2 s（最大食）。root_tolerance = 0.01 s（numerical-policy §A5 正本）。
- 食分 ±0.0005（0.001食分≈1.9″, accuracy §2.2）。食面積 ±0.0005 相当（要確認: accuracy は独立許容を明記せず食分基準に準拠と解釈, ISSUE-027）。純幾何（円交差）≤ 1e-10。
- 高度・方位（幾何学的）≤ 0.01°（表示・可視性用途, ISSUE-028）。可視性境界はフィクスチャで決定的に一致（pass/fail）。
- **UTC 絶対時刻は将来 ΔT/UT1 律速**（accuracy §0(b)/§2.3）。**局地接触の幾何分は δUT1（→μ）混入**（accuracy §2.1L, ISSUE-021 I1）。幾何（TT）と分離して報告。許容を通すための拡大禁止（conventions §11）。

---

## 出典

- **Explanatory Supplement to the Astronomical Almanac (3rd ed.), §11「Eclipses」**: 局地予報式（観測者基本面座標 ξ,η,ζ、接触 `(ξ−x)²+(η−y)²=L²`、`L=l−ζ·tan f`、最大食、magnitude、位置角）。**式番号は要確認**（一次資料＝書籍で確認, 推測しない）。**局地時角の符号（θ=μ−λ か μ+λ）も式番号まで要確認**。
- **Meeus, J. "Astronomical Algorithms", 2nd ed.**: Ch.54「Solar Eclipses」（u, v, n², 接触時刻補正, magnitude）、Ch.13「Transformation of Coordinates」式 13.5/13.6（赤道→地平, **方位は南基準** → 北 0 東回りへ変換）、Ch.16「Atmospheric Refraction」式 16.3/16.4（大気差, Bennett/Sæmundsson）。Ch.54 の個別式番号は**要確認**。
- **Weisstein, E. MathWorld "Circle-Circle Intersection"**: 2 円交差面積（lens area）の標準幾何式（食面積 obscuration）。標準公式のため出典は公式名で明記。
- **SOFA**（参照のみ・移植しない, data-sources §0）: `iauHd2ae` / `iauAe2hd`（hour angle/dec ↔ az/el）。**時角・恒星時は ERA(`iauEra00`)経由 CIO ベース・分点 GST（`iauGst06a`）禁止**（D4, §2）。
- **NASA / Espenak / USNO**（GSFC eclipse site, data-sources §4）: 地点別接触時刻・最大食時刻・食分・高度（第二義・整合チェック, k・ΔT 慣習を揃える, accuracy §3.1）。
- **WGS84**（conventions §4.1）: Re = a = 6 378 137.0 m（ξ,η,ζ,L1,L2,m の無次元化基準）。
- プロジェクト内: algorithms.md §0（記号表）, conventions §1/§2/§3/§5/§7/§8/§9/§11, accuracy §0/§2.1/§2.1L/§2.2/§2.3/§3/§4/§6, numerical-policy §A2/§A5, §8（08-global.md, 探索窓・全球接触）, §10（10-observer.md, ρcosφ′/ρsinφ′）, §2（02-frames.md, μ=ERA−α_axis・CIO 統一）, ISSUE-024/025/026/027/028, ISSUE-021/037/022/039/015/010/011。
