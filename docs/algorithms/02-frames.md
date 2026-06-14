# §2 フレーム変換（IAU2006 歳差 / IAU2000A 章動 / ERA・CIO）

> 正本: `algorithms.md §0`（記号表）・`conventions.md §5/§5.1/§5.2`。本セクションはそれらに**厳密準拠**する。
> 状態: ドラフト（Milestone 0）。一次資料で未確認の式番号・符号は「**要確認**」を残す。
> 関連 Issue: ISSUE-035（フレーム連鎖）, ISSUE-039（μ・GAST 素材供給）, ISSUE-040（IAU2000A 章動係数）。

---

## 目的と入出力（単位・時刻系・座標系・フレーム）

地球姿勢に基づくフレーム変換連鎖を**回転行列**として供給する。`conventions §5` の連鎖を実装する:

```
GCRS →[frame bias + 歳差(IAU2006) + 章動(IAU2000A)]→ CIRS →[ERA(UT1)]→ TIRS →[極運動 xp,yp; s′]→ ITRS
```

| 項目 | 内容 |
|---|---|
| 入力 | `time_tt: TtInstant`（歳差章動・s・s′・黄道傾斜）、`time_ut1: Ut1Instant`（ERA）、`xp, yp: Radians`（極運動, EOP）、影軸赤経 `α_axis: Radians`（CIRS 基準, §2 で μ 構成に使用） |
| 出力 | 回転 `Matrix3`（直交・右手系・無次元）。CIO ベース見かけ恒星時量 `GAST` と影軸時角 `μ`（`Radians`, `[0, 2π)`） |
| 単位 | 角度 = rad（conventions §2）、回転 = 無次元 Matrix3 |
| 時刻系 | 歳差章動・X,Y,s・s′・黄道傾斜 = **TT**、ERA = **UT1**（conventions §6）。EOP（xp,yp,δUT1）はデータ供給（accuracy §5, data-sources §3.1） |
| 座標系/フレーム | GCRS / CIRS / TIRS / ITRS（conventions §5 表）。**CIO ベースで統一**。分点ベース GST 経路は**作らない**（conventions §5.2 / D4） |

**確定方針（厳守）**:
- フレーム連鎖は CIO ベース（X, Y, s）で構成し、赤道分点ベースと混在させない（conventions §5.2 / D4, ISSUE-035, ISSUE-039）。
- 見かけ恒星時・μ・影軸赤経 `α_axis` はすべて **ERA 経由 CIO ベース**。**分点 GST（`iauGst06a` 等）は Standard で禁止**（conventions §5.2 / D4）。
- s′（TIO locator）は微小（数十 µas）だが省略しない（精度最優先, ISSUE-035 実装メモ）。

---

## 記号（algorithms.md §0 参照。固有補助記号のみ）

§0 で未定義の本セクション固有記号のみ列挙する（§0 へ追記候補は本ファイル末尾「記号表への追記提案」）。

| 記号 | 意味 | 単位 |
|---|---|---|
| `t` | J2000 からの TT ユリウス世紀 = `((JD_TT − 2451545.0))/36525`（§0, numerical-policy §A1） | 無次元 |
| `X`, `Y` | CIP（Celestial Intermediate Pole）の GCRS 単位方向座標 | 無次元（≈rad の方向余弦） |
| `s` | CIO locator（CIP 赤道上での CIO 位置を定める量） | rad |
| `s′` | TIO locator（TIRS の x 軸位置を定める量） | rad |
| `θ_ERA` | Earth Rotation Angle（= ERA, §0 の ERA と同義） | rad |
| `EO` | equation of the origins（CIO と分点の CIP 赤道上角距離） | rad |
| `Δψ`, `Δε` | 黄経章動・黄道傾斜章動（IAU2000A） | rad |
| `ε_A` | IAU2006 平均黄道傾斜（mean obliquity of date） | rad |
| `xp`, `yp` | 極運動（CIP の ITRS に対する位置, EOP） | rad |
| `Q(t)` | GCRS→CIRS の bias-precession-nutation（CIO ベース）回転行列 | 無次元 |
| `R(θ)` | ERA 回転（CIRS→TIRS） | 無次元 |
| `W(t)` | 極運動回転（TIRS→ITRS） | 無次元 |
| `α_axis` | 影軸赤経（CIRS 基準, §0 の α_axis） | rad |

行列記法: `R1(φ), R2(φ), R3(φ)` は x/y/z 軸まわりの右手系基本回転（IERS Conventions 2010, ch.5 の符号規約に従う）。

---

## 数式（番号付き・各式に出典）

### 2.1 GCRS → CIRS（frame bias + IAU2006 歳差 + IAU2000A 章動, CIO ベース）

CIP 座標 `X, Y`（GCRS におけるCIP方向）と CIO locator `s` から、CIO ベースの celestial-to-intermediate 行列を構成する。

**(F1) celestial-to-intermediate 行列（CIO ベース）**
出典: IERS Conventions 2010, ch.5, eq. (5.10)。SOFA `iauC2ixys(X, Y, s)`（= `iauC2i06a` 内部）。

```
              ⎡ 1 − a·X²      −a·X·Y       X ⎤
Q(t) = R3(s)· ⎢  −a·X·Y      1 − a·Y²      Y ⎥        (a = 1/(1 + Z),  Z = √(1 − X² − Y²))
              ⎣  −X            −Y           Z ⎦
```

ここで `R3(s)` は z 軸まわり角 `s` の回転。`a = 1/(1+Z)` は SOFA `iauC2ixys` の定義（IERS Conventions 2010, eq. (5.10) の `a = 1/(1+cos d) ≈ 1/2 + (X²+Y²)/8`、ただし**厳密形 `a=1/(1+Z)` を採用**, ISSUE-035）。
`r_CIRS = Q(t) · r_GCRS`。

**(F2) X, Y, s の供給**
出典: IERS Conventions 2010, ch.5。SOFA `iauXys06a(date1, date2) → X, Y, s`（内部で `iauPnm06a` の NPB 行列から `iaubpn2xy` で X,Y、`iauS06` で s）。

- NPB（bias-precession-nutation）行列: SOFA `iauPnm06a`（IAU2006 歳差 + IAU2000A 章動 + frame bias）。歳差: Capitaine, Wallace & Chapront (2003) A&A 412, 567（P03, IAU2006）、Fukushima-Williams 角は `iauPfw06`。章動: Mathews, Herring & Buffett (2002)（IAU2000A）、SOFA `iauNut00a`（lunisolar 678 + planetary 687 = 1365 項, ISSUE-040）。
- `X, Y`: NPB 行列第3行から `X = NPB[2][0]`, `Y = NPB[2][1]`（SOFA `iaubpn2xy`）。**要確認**: 行列要素の並び（行優先/列優先）を実装規約と一致させる。
- `s`: CIO locator 級数 SOFA `iauS06(date, X, Y)`（IERS Conventions 2010, table 5.2c の級数）。`s` は `XY/2` 項を含む（`s = s(t) − X·Y/2` の形, IERS eq. (5.11)）。**要確認**: 級数表の打切りと符号。

**(F3) 内部粗スキャン用（IAU2000B 章動・非公開）**
出典: SOFA `iauNut00b`（77 項）。**公開出力（Standard = IAU2006/2000A）では使用せず、内部粗スキャン（候補棄却専用・非公開）でのみ任意（非既定）に許容**する。`gcrs_to_cirs_matrix_2000b` で X,Y を近似（ISSUE-035, 内部粗スキャン用として存置）。追加誤差 ~1mas 級（accuracy §2.1 で 0.05″ 配分に余裕）。

### 2.2 CIRS → TIRS（Earth Rotation Angle）

**(F4) ERA**
出典: IAU2000 定義。IERS Conventions 2010, ch.5, eq. (5.15)。SOFA `iauEra00(UT1a, UT1b)`。

```
θ_ERA = 2π · ( 0.7790572732640 + 1.00273781191135448 · (JD_UT1 − 2451545.0) )      (rad, [0,2π) 正規化)
```

ここで `JD_UT1 − 2451545.0` は **JulianDate2 の整数部・小数部を分離して計算**（numerical-policy §A1: エポック減算を整数部側で厳密に。係数は magic number でなく IERS eq. (5.15) の定義値, ISSUE-039/041）。

**(F5) CIRS→TIRS 回転**
出典: IERS Conventions 2010, ch.5, eq. (5.5)。SOFA `iauC2tcio` の R3 部分。

```
R(θ_ERA) = R3(θ_ERA)
r_TIRS = R3(θ_ERA) · r_CIRS
```

### 2.3 TIRS → ITRS（極運動 + TIO locator s′）

**(F6) 極運動行列**
出典: IERS Conventions 2010, ch.5, eq. (5.3)。SOFA `iauPom00(xp, yp, sp)`。

```
W(t) = R3(−s′) · R2(xp) · R1(yp)
r_ITRS = W(t) · r_TIRS
```

**(F7) TIO locator s′**
出典: IERS Conventions 2010, ch.5, eq. (5.13)。SOFA `iauSp00(date1, date2)`。

```
s′ = −0.0015 · ( a_c²/1.2 + a_a² ) · t      [arcsec → rad に変換]
```

実装では SOFA `iauSp00` の主要項 `s′ ≈ −47 µas · t`（t = TT 世紀）を採用（IERS eq. (5.13) の数値係数, ISSUE-035 実装メモ。係数は出典付き定数, magic number 禁止）。**要確認**: 採用次数（主要項のみ vs 全項）。

### 2.4 連鎖合成 GCRS → ITRS

**(F8)** 出典: IERS Conventions 2010, ch.5, eq. (5.1)（CIO ベース）。SOFA `iauC2t06a` 相当。

```
[GCRS→ITRS](t) = W(t) · R3(θ_ERA) · Q(t)
r_ITRS = W(t) · R3(θ_ERA) · Q(t) · r_GCRS
```

逆変換は転置（直交行列）。

### 2.5 見かけ恒星時 GAST（ERA 経由 CIO ベース）

**分点 GST を経由せず**、ERA と equation of the origins `EO`（CIO locator・NPB 由来）から構成する（conventions §5.2 / D4, ISSUE-039）。

**(F9) GAST**
出典: Wallace & Capitaine (2006) A&A 459, 981, eq. (62)「GST = θ − EO」。SOFA `iauGst06`（内部で `iauEors` を使用）。

```
GAST = θ_ERA − EO            (rad, [0,2π) 正規化)
```

> 注意（CIO 統一の根拠）: SOFA `iauGst06a` は分点ベースの恒星時を返すため**使用しない**。本実装は ERA に CIO ベースの EO を差し引いて見かけ恒星時相当量を作り、**μ 構成と同一の CIO 基準**に揃える（conventions §5.2 / D4, ISSUE-039 受入テスト「分点経路を使っていないことの証明」）。

**(F10) equation of the origins EO**
出典: Wallace & Capitaine (2006) A&A 459, 981, eq. (66) および NPB+s からの構成（SOFA `iauEors(rnpb, s)`）。

```
EO = s − atan2( NPB[0][1], NPB[0][0] )      (SOFA iauEors の構成)
```

すなわち NPB 行列第1行（CIO ベース X 軸の GCRS 成分）と s から `EO` を得る。これにより `GAST = θ_ERA − EO` は CIO 基準で閉じる。**要確認**: `iauEors` の符号規約（`EO = s − E`, `E = atan2(NPB[0][1], NPB[0][0])`）を実装と一致させる。

### 2.6 影軸赤経 α_axis と影軸時角 μ（CIO 基準）

**(F11) μ（見かけグリニッジ時角）**
出典: conventions §5.2 / D4・algorithms.md §0。

```
μ = GAST − α_axis            (rad, [0,2π) 正規化)
```

- `α_axis` は **CIRS 基準の影軸赤経**（ISSUE-020 が供給、§6 で算出）。`GAST` も CIO ベース（F9）なので μ は CIO 一貫。
- μ は UT1（→δUT1）に依存（θ_ERA 経由）。この帰結は accuracy §0(a) 脚注・§2.1L 参照。**μ の UT1 混入を doc に明記**（ISSUE-039 実装メモ）。
- NASA ベッセル μ は分点基準・度/hour（D5 対応表）。本実装の CIO ラジアン量との**系統差は誤差化せず accuracy.md に記録**（ISSUE-039, I4/I5）。

---

## 手順（実装順・数値注意・適用順序）

1. **t 生成**（TT 世紀）と `JD_UT1`（UT1）を `JulianDate2` から生成。エポック 2451545.0 の減算は整数部側で厳密（numerical-policy §A1）。**生 f64 の JD を渡さない**（conventions §6, §11）。
2. **Q(t)**（GCRS→CIRS）: `iauPnm06a` 相当で NPB を構成 → X,Y を抽出 → `iauS06` で s → (F1) で `Q`。**公開出力（Standard）は IAU2006/2000A 固定、内部粗スキャン（非公開）のみ 2000B（F3）**。選択を `CalculationMetadata` に記録（ISSUE-035）。
3. **θ_ERA**（UT1）: (F4)。`[0,2π)` 正規化。
4. **W(t)**（TIRS→ITRS）: `iauSp00` で s′ → (F6)。xp,yp は EOP から（データ層, ISSUE-035 非目的）。
5. **連鎖合成** (F8)。逆は転置。
6. **GAST**: (F10) で EO → (F9) で `GAST = θ_ERA − EO`（**分点 GST を経由しない**）。
7. **μ**: `α_axis`（CIRS, §6）を受けて (F11) `μ = GAST − α_axis`。

**正規化の使い分け**（conventions §2）: 恒星時・赤経・μ は `[0,2π)`。求解側（合・接触）で連続性が要るときは呼び出し側で連続化（本層は素の値）。
**数値方針**: 角度係数・ERA 係数は出典付き定数（magic number 禁止, conventions §11 / ISSUE-039）。章動 1365 項は手書き直書き禁止、packed 係数を読む（ISSUE-040）。

---

## 境界・特異・異常系

- `Z = √(1 − X² − Y²)`: X²+Y² は ~10⁻⁶ 規模で `Z≈1`、`a=1/(1+Z)` は安定。クランプ不要だが `1 − X² − Y² < 0`（理論上起きない）は異常として検出（精度最優先, 安全弁）。
- `atan2` 由来（EO, μ）: `[0,2π)` 正規化を明示（conventions §2）。μ→μ+2π 不変をプロパティテスト（ISSUE-039 L8）。
- EOP 欠損（xp,yp,δUT1）: データ層の責務。欠損時は `OutOfSupportedRange` 等で上流に返す（本層は値を受けて回転を作るのみ, ISSUE-035 非目的）。
- s′ の符号・単位（arcsec→rad）の取り違えは数十 µas 系統誤差 → 単体テストでガード（ISSUE-035）。

---

## 検証（二段ゲート: M2 暫定 = Mock+SOFA+NASA / M10 = JPL DE 差分）

accuracy §3.1 / §3.3、ISSUE-035/039/040 受入テスト準拠。基準値は SOFA/IERS 公開ベクトルから取得し**実装へハードコードしない**（conventions §11）。

- **L3 / SOFA 突合（M2 暫定）**:
  - `Q(t)`（GCRS→CIRS）, `θ_ERA`, `W(t)`, 連鎖 `gcrs_to_itrs_matrix` の各要素を SOFA 参照値（`iauC2i06a`/`iauEra00`/`iauPom00`/`iauC2t06a`）と比較。残差 **~1mas 級**（許容 0.05″ に余裕, accuracy §2.1）。
  - 章動 Δψ, Δε（`iauNut00a`）を SOFA 参照値と **µas 級**で一致（ISSUE-040 ラウンドトリップ）。
  - **CIO 統一の証明**: `GAST`（F9）が分点 GST（`iauGst06a`）と「CIO−分点の既知量（EO の差・赤経起点差）」だけ差を持つことを確認 = 内部で分点経路を使っていない証明（ISSUE-039）。差は誤差でなく系統差として accuracy.md に記録（I4/I5）。
- **ラウンドトリップ**: GCRS→ITRS→GCRS が恒等（≲数 µas, ISSUE-035）。μ→μ+2π 不変（ISSUE-039）。
- **時刻系分離**: ERA は UT1、CIP/s は TT で評価していること（TT を ERA に渡すと既知量ずれる回帰, ISSUE-039）。
- **2000A vs 2000B 差分**が既知オーダー（~1mas 級）に収まる（内部粗スキャンで 2000B を許容する根拠。公開出力は 2000A, ISSUE-035）。
- **M10 DE 確定**: DE 差分パイプライン（accuracy §3.1）に組み込み、フレーム由来残差が層分解で **0.05″ 以下**に帰属（accuracy §3.3, §4）。

---

## 出典

- IERS Conventions (2010), ch.5（フレーム変換・歳差章動・ERA・極運動）: eq. (5.1) 連鎖, (5.3) 極運動, (5.5) ERA 回転, (5.10) C2I 行列, (5.11) s, (5.13) s′, (5.15) ERA。
- Capitaine, Wallace & Chapront (2003), A&A 412, 567（IAU2006/P03 歳差）。
- Mathews, Herring & Buffett (2002)（IAU2000A 章動）。
- Wallace & Capitaine (2006), A&A 459, 981, eq. (62)/(66)（GST = θ − EO, equation of the origins）。
- SOFA 関数: `iauPnm06a`, `iauPfw06`, `iauNut00a`, `iauNut00b`, `iauC2ixys`, `iauXys06a`, `iauS06`, `iauEra00`, `iauPom00`, `iauSp00`, `iauC2i06a`, `iauC2t06a`, `iauEors`, `iauObl06`（参照のみ・移植しない, data-sources §0 / conventions §11）。
- conventions §5/§5.1/§5.2, accuracy §2.1/§2.3/§3, numerical-policy §A1, data-sources §3.1。
- 関連 Issue: ISSUE-035, ISSUE-039, ISSUE-040。
