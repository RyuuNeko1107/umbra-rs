# §8 全球条件（Global circumstances, Step9）

本節は `umbra-rs` の**全球的な日食状況**（影軸の地心最小距離 gamma・種別判定・全球接触 P1/U1/最大/U4/P4・最大食地点・食帯幅・中心食継続時間）の数式・手順・符号規約を定める。地球を 1 天体として扱い、観測者は登場しない（観測者依存は §9 局地条件の責務）。
正本は `algorithms.md §0`（記号表）。本節はそれに厳密準拠し、本節固有の補助記号のみ追加定義する。
関連 Issue: ISSUE-023（全球分類・gamma・全球接触・最大食地点・帯幅・中心食継続）。供給源は ISSUE-037（直接評価, fit 誤差ゼロ・最大食精解の既定）/ ISSUE-022（多項式, 経路・帯幅多点）。

> 状態: ドラフト（Milestone 0）。**種別境界しきい値の数値は ISSUE-023 同様「要確認」**（l1/l2 依存。一次資料 Meeus Ch.54 / NASA で最終確認）。一次資料で式番号を確認できない箇所は「**要確認**」を残し、推測で式番号を書かない。

---

## 目的と入出力（単位・時刻系・座標系・フレーム）

- **目的**: ベッセル要素（瞬時 ISSUE-021 / 直接供給 ISSUE-037 / 多項式 ISSUE-022）から、日食の全球状況を確定する（architecture §3/§7, api-draft §3.4）。
  1. **gamma**: 影軸の地心最小距離（Re）。
  2. **種別**: Partial / Annular / Total / Hybrid / NonCentralAnnular / NonCentralTotal（`SolarEclipseKind`, api-draft §3.4）。
  3. **全球接触**: P1（部分食開始）/ U1（中心食開始）/ 最大食 / U4（中心食終了）/ P4（部分食終了）（conventions §8）。
  4. **最大食地点・食帯幅・中心食継続時間**（`GreatestEclipse`）。
- **入力**: `BesselianSource`（x, y, d, μ, l1, l2, tan f1, tan f2 を任意 TT で供給。x, y, l1, l2 は Re 無次元、d, μ は rad）、`TimeScales`（TT↔UTC, accuracy §0）、`EngineConfig`（k 値選択 conventions §9・`earth_model`）。
- **出力**: `GlobalCircumstances`（gamma[Re]、種別、各全球接触の `GeoPoint`+TT+UTC、`GreatestEclipse`）。
- **単位**: 角度 = rad（conventions §2）。ベッセル長 = Re（conventions §1/§4.1, Re = WGS84 a）。帯幅 = km、継続 = 秒、高度 = 度。
- **時刻系**: **gamma・最大食「時刻TT」は影軸幾何のみ＝純 TT**（δUT1 非依存, accuracy §0(a) 脚注 / §2.1L）。最大食「地点」（経度）は μ 経由で δUT1 依存。各接触・最大食は **UTC と TT の両方**を返す（conventions §6, accuracy §0）。
- **フレーム**: FundamentalPlane（影軸 z、x̂=東、ŷ=天の北の射影, conventions §5）で求解 → 地表貫通点を測地座標へ（geocentric→geodetic は §10/ISSUE-010・011, WGS84）。

---

## 記号（algorithms.md §0 を参照。本節固有の補助記号のみ定義）

algorithms.md §0「ベッセル要素」表（x, y, d, μ, l1, l2, f1, f2）と「食の量」表（gamma, magnitude, obscuration）を正本とする。本節で補助的に用いる記号:

| 記号 | 意味 | 出典/備考 |
|---|---|---|
| `g(t)` | 影軸の地心距離の二乗 `g(t) = x(t)² + y(t)²`（Re²） | 式 8.1。最小化対象（√ 回避, numerical-policy §A5） |
| `ρ_axis(t)` | 影軸の地心距離 `√(x²+y²)`（Re）。gamma = min ρ_axis | 式 8.2 |
| `t_max` | gamma を与える時刻（最大食 TT。`g'(t)=0` の根） | 式 8.3 |
| `ρ_g` | 地球縁の基本面射影半径（Re）。扁平を考慮した有効値 | 式 8.6・要確認 |
| `ρ_pen`, `ρ_umb` | 半影縁 l1 / 本影縁 l2 が地球縁に外接する時の影軸地心距離 | 式 8.7/8.8 |
| `c_pen` | 部分食限界しきい値（影軸が外れても半影が地球に触れる上限）。≈ 1 + l1（要確認） | 式 8.4・ISSUE-023 |
| `c_cen` | 中心食限界しきい値（影軸が地球に当たる上限）。≈ ρ_g（要確認） | 式 8.5・ISSUE-023 |
| `w_path` | 食帯幅（km）。最大食点での本影/反本影の地表幅 | 式 8.11 |
| `D_cen` | 中心食継続時間（秒）。最大食点を本影/反本影が通過する時間 | 式 8.12 |

> 本節は新規の独自記号を最小化する。gamma・x・y・l1・l2・d・μ は §0 正本。ρ_g, c_pen, c_cen はしきい値・地球縁射影の補助で、数値は要確認（ISSUE-023）。記号表への追記候補は末尾「報告」に記載。

---

## 数式（番号付き・各式に出典）

### A. gamma と最大食時刻（影軸の地心最小距離）

最大食 = 影軸が地心に最も近づく瞬間。距離の二乗を最小化（中心線尖点での √ 不連続を避けるため最小化対象は二乗, numerical-policy §A5 / accuracy §2.1 D2）。

**(8.1)** 影軸の地心距離の二乗:
```
g(t) = x(t)² + y(t)²            [Re²]
```

**(8.2)** 影軸の地心距離:
```
ρ_axis(t) = √(x(t)² + y(t)²)    [Re]
```

**(8.3)** 最大食時刻と gamma:
```
t_max : g'(t) = 0  ⇔  x·x' + y·y' = 0    （Brent 求根, 粗ブラケットに最小化併用）
gamma = ρ_axis(t_max) = √( x(t_max)² + y(t_max)² )    [Re]
```

- 出典: **Explanatory Supplement to the Astronomical Almanac (3rd ed.), Ch.11**（全球状況・gamma・central eclipse の条件）/ **Meeus, *Astronomical Algorithms* (2nd ed.), Ch.54「Eclipses」**（gamma 定義・u（=l2）・種別判定の実用式）/ NASA Espenak（gamma 定義, data-sources §4.1）。
- 数値手法（D2, accuracy §2.1）: 最大食は **`g'(t)=0`（= `x·x'+y·y'=0`）の Brent 求根を正式手法**とする。距離（最小化）は **粗ブラケット用に降格**（皆既帯の平底＝`g'≈0` 区間で最小化が劣化するため, numerical-policy §A5）。`x'`, `y'` は中心差分＋Richardson か解析微分（numerical-policy §A2(3)）。
- gamma の符号: 慣習上 gamma は符号付き（影軸が地心の北を通れば正・南を通れば負）で表記されることがある（NASA）。本実装は内部で `√(x²+y²) ≥ 0` を基本とし、**符号付き gamma を返す場合は最大食時の y の符号で与える**（要確認: NASA の符号定義 = 「最大食時に影軸が赤道面を横切る側」を一次資料で確認）。種別判定（§B）は `|gamma|` を用いる。

### B. 種別判定（皆既/金環/ハイブリッド/部分/非中心）

判定は **(i) `|gamma|` と限界しきい値の比較**（影軸が地球に当たるか）と、**(ii) 中心食区間中の l2 の符号**（皆既 vs 金環）の 2 段で行う。

**(8.4)** 部分食限界しきい値（半影が地球に届く上限）:
```
c_pen ≈ ρ_g + l1            [Re]   （要確認: 厳密形は l1 依存・地球扁平込み）
```

**(8.5)** 中心食限界しきい値（影軸が地球面に当たる上限）:
```
c_cen ≈ ρ_g                [Re]   （要確認: 非中心境界・地球扁平込み）
```

**(8.6)** 地球縁の基本面射影半径（扁平込み・要確認）:
```
ρ_g ≈ 1                     [Re]   （球近似。WGS84 扁平を考慮した有効半径は要確認, ISSUE-023）
```

**種別フロー**（Meeus Ch.54 種別フローチャート / NASA 種別定義）:

1. **`|gamma| > c_pen`** → 影軸も半影も地球を外れる（食なし）。
2. **`c_cen < |gamma| ≤ c_pen`** → 半影のみ地球に触れる（影軸は外れる）→ **Partial**（U1/U4 = None）。
3. **`|gamma| ≤ c_cen`** → 影軸が地球に当たる（中心食候補）。最大食時の l2 の符号と、中心食区間中の l2 符号変化で分岐:
   - 食の経過全体で **l2 > 0**（本影頂点が地表手前 = 反本影が地表）→ **Annular**（金環）。
   - 食の経過全体で **l2 < 0**（本影が地表に到達）→ **Total**（皆既）。
   - 食の経過中に **l2 の符号が反転**（金環⇄皆既）→ **Hybrid**（annular-total）。
   - 影軸は当たるが中心食条件を**地球縁で一部のみ満たす**（`|gamma|` が c_cen 近傍）→ **NonCentral**（NonCentralAnnular / NonCentralTotal, l2 符号で細分）。

- 出典: Meeus Ch.54 種別フローチャート / NASA 種別定義。**l2<0 ⇒ 皆既、l2>0 ⇒ 金環**（algorithms.md §0 / ISSUE-021 符号規約）。
- **要確認（ISSUE-023 同様）**: `c_pen`（≈1.55 が文献にあるが l1 依存・要確認）、`c_cen`（非中心境界）、`ρ_g`（扁平込み有効地球縁半径）の**数値しきい値は一次資料（Meeus Ch.54 / NASA）で最終確認**。閾値の正確な形（l1, l2, 地球扁平への依存）を実装コメントに転記し、推測値を確定しない。
- **k 値選択の系統差（conventions §9）**: l1/l2 は k（`IauMean`=0.2725076 / `EspenakUmbral`=0.272281 / `EspenakPenumbral`=0.2725076）に依存し、**皆既/金環/ハイブリッド境界が k で動く**。NASA 照合時は Espenak 2値慣習へ切替え、系統差を accuracy.md に記録（誤差を隠さない, accuracy §0/§3.1）。`metadata.lunar_radius_model` に必ず載せる。

### C. 全球接触 P1 / U1 / U4 / P4

全球接触 = 影錐の縁が地球縁に最初/最後に外接する時刻（conventions §8 外接の全球版）。

**(8.7)** 半影縁の地球縁外接（P1: 開始 / P4: 終了）:
```
P1, P4 :  ρ_axis(t) = ρ_g + l1(t)    ⇔   x² + y² = (ρ_g + l1)²
```

**(8.8)** 本影/反本影縁の地球縁外接（U1: 中心食開始 / U4: 中心食終了）:
```
U1, U4 :  ρ_axis(t) = ρ_g + |l2(t)|   ⇔   x² + y² = (ρ_g + |l2|)²
```

**(8.9)** 求根対象（差の連続関数, Brent）:
```
g_P(t) = (x² + y²) − (ρ_g + l1)²        （P1/P4 の零点）
g_U(t) = (x² + y²) − (ρ_g + |l2|)²      （U1/U4 の零点）
```

- 出典: Explanatory Supplement Ch.11 / Meeus Ch.54（全球接触）。conventions §8（外接=中心間距離 = 太陽視半径 + 月視半径 の全球版）。
- 手順: 探索窓（最大食 ±マージン or P1..P4 概算）で `g_P`, `g_U` を粗走査し符号変化区間を検出 → Brent（numerical-policy §A5, root_tolerance ≤ 0.01 s）。`g_P` は最大食を挟んで 2 根（P1 < t_max < P4）、`g_U` も同様（U1 < t_max < U4）。
- **U1/U4 は中心食でなければ None**（部分食では本影が地球縁に届かず `g_U` に符号変化なし, ISSUE-023）。Option を正しく埋める。
- 各接触点の地理座標: その時刻に影縁が地球縁に接する点を測地座標へ（geocentric→geodetic, §10/ISSUE-010・011）。
- **要確認**: `ρ_g`（地球縁の有効半径・扁平込み, 式 8.6）と l1/l2 への高さ補正の有無（全球版では観測者 ζ がないが、地球縁では扁平が効く）。一次資料 Explanatory Supplement Ch.11 で式番号と扁平の扱いを確認。

### D. 最大食地点・食帯幅・中心食継続時間

**(8.10)** 最大食地点（影軸の地表貫通点）:
```
最大食点 = （t_max における影軸が WGS84 地表を貫く点）を測地座標 (φ, λ_east) へ
```
- 影軸方向（d, μ）と影軸の基本面交点 (x, y) から地表貫通点を求め、geocentric→geodetic（§10）。経度は μ 経由（δUT1 依存）。出典: ISSUE-023 / Explanatory Supplement Ch.11。

**(8.11)** 食帯幅（km, 中心食のみ Some）:
```
w_path ≈ （最大食点での本影/反本影の地表での幅）
```
- 本影/反本影の基本面半径 |l2| を地表へ投影した幅。出典: Meeus Ch.54 / NASA（path width）。**要確認**: 地表傾斜（太陽高度）・扁平を含む厳密式は一次資料で確認。

**(8.12)** 中心食継続時間（秒, 中心食のみ Some）:
```
D_cen ≈ （本影/反本影が最大食点を通過する時間）
```
- 出典: Meeus Ch.54 / NASA（duration of totality/annularity）。**要確認**: 継続時間式（影の地表速度・観測点の自転速度の差）の式番号を一次資料で確認。
- **部分食・非中心では w_path = None, D_cen = None**（ISSUE-023, Option 設計）。

**(8.13)** 最大食点の magnitude / obscuration:
- 最大食点での食分 magnitude と食面積 obscuration は **§9 の式（9.x）を最大食点の (m, L1, L2)・視半径で評価**（acos クランプ・5 境界, §9 / ISSUE-027）。値の重複定義を避け §9 を参照する。

---

## 手順（実装順・数値注意=numerical-policy 参照）

1. **最大食時刻 t_max と gamma**: `g'(t)=x·x'+y·y'=0` を Brent 求根（D2, numerical-policy §A5）。粗ブラケットは `g(t)=x²+y²` の最小付近 3 点。`x', y'` は中心差分＋Richardson か解析微分（numerical-policy §A2(3)）。gamma = √(g(t_max))。**最大食供給源は ISSUE-037（直接, fit 誤差ゼロ）を既定**（精度）。
2. **種別判定**: `|gamma|` と c_pen/c_cen（式 8.4/8.5, **要確認**）で Partial / 中心食を分岐。中心食は中心食区間中の l2 符号で Annular / Total、符号反転で Hybrid、`|gamma|` が c_cen 近傍で NonCentral（§B フロー）。l2 符号境界（皆既↔金環）で種別が**連続に**切替わること。
3. **全球接触**: P1/P4（`g_P` 零点）, U1/U4（`g_U` 零点）を粗走査→Brent（式 8.9, numerical-policy §A5）。中心食でなければ U1/U4 = None。
4. **最大食地点**: 影軸地表貫通点を測地座標へ（式 8.10, §10/ISSUE-010・011）。太陽高度を算出（§9 高度方位 / ISSUE-028）。
5. **帯幅・中心食継続**: 中心食のみ式 8.11/8.12（要確認）。それ以外 None。
6. **magnitude / obscuration**: 最大食点で §9 の式を評価（acos クランプ, §9 / ISSUE-027）。
7. **時刻系**: 各接触・最大食を **TT と UTC の両方**で返す（accuracy §0）。将来日食は `delta_t_uncertainty_seconds` を metadata に（accuracy §2.3）。

**数値注意**（numerical-policy）:
- 求根・最小化は **無条件 Newton 禁止**（conventions §11）。Brent（ブラケット必須, §A5）。root_tolerance = 0.01 s 既定（solver 配分 0.05″ の 1/10, §A5）。
- acos/asin は `[-1, 1]` クランプ（§A5 / accuracy §2.2）。
- x, y の連続化（±π 折返し除去）後に g(t), g'(t) を扱う（conventions §2）。
- magic number 禁止: c_pen/c_cen/ρ_g は出典付き定数（要確認しきい値は出典コメントを残し推測確定しない, conventions §11）。

---

## 境界・特異・異常系

- **影が地球を外す（`|gamma| > c_pen`）**: 食なし。全接触 None・種別なし（上位で扱う）。
- **部分食（`c_cen < |gamma| ≤ c_pen`）**: U1/U4 = None、w_path/D_cen = None。P1/P4 のみ。
- **非中心（`|gamma| ≈ c_cen`）**: 影軸が地球縁すれすれ。中心食条件を地球縁で一部のみ満たす。NonCentral 系。v1.0 公開可否は api-draft §6 未確定（型は用意・分類ロジックは実装, ISSUE-023）。
- **l2 ≈ 0（皆既↔金環境界）**: Hybrid。中心食区間中の l2 符号反転を検出。境界で種別が連続に切替わること（ISSUE-023 種別境界テスト）。
- **皆既帯平底**（`g'≈0` が区間で成立）: 最大食時刻が一意でない。平底区間の代表点（中央 or 全球最大食に最も近い点）を規約で定義（accuracy §2.1 / numerical-policy §A5）。
- **gamma 符号境界**（最大食時 y≈0）: 符号付き gamma の符号が反転。`|gamma|` 判定は連続。
- **k 値境界**: `EspenakUmbral` vs `IauMean` で皆既/金環/ハイブリッド判定がずれる。系統差を accuracy.md に記録（誤差を隠さない, §B）。
- **求根ブラケット不成立**: `RootNotBracketed` 等で上流へ（粗走査刻みを掠め食でも符号変化を逃さない細かさに, 偽陰性ゼロ・architecture §3）。

---

## 検証（受け入れ基準・基準値の出典）

> オラクル値は実装へコピーしない（conventions §11）。NASA/USNO 公開値・MockEphemeris 人工配置・DE 差分から独立に算出する。accuracy テストレベル **L5（全球日食）**。

- **NASA 全球整合（第二義, accuracy §3.1 / data-sources §4.1）**: NASA 5千年カタログの既知日食（皆既/金環/ハイブリッド/部分/非中心）で **種別・gamma・最大食時刻・食分・帯幅**を比較。**k 慣習を Espenak（EspenakUmbral/Penumbral, conventions §9）に揃え、ΔT も合わせる**。系統差を accuracy.md に記録（絶対基準にしない, §3.1）。基準は fixtures。
- **DE 差分（第一義, accuracy §3.1）**: 解析暦 vs DE440 で同一全球パイプラインを通し gamma・最大食時刻の差を層分解（accuracy §4）。
- **MockEphemeris 人工ケース**: 完全中心皆既（gamma≈0, l2<0, Total）/ 明確な金環（l2>0, Annular）/ ハイブリッド（l2 符号反転, Hybrid）/ 部分（gamma 大, Partial, U1/U4=None）/ 非中心（gamma≈c_cen すれすれ, NonCentral）。各で種別・接触 Option を検証。
- **種別境界テスト（必須）**: l2≈0（皆既↔金環）と `|gamma|≈c_cen`（中心↔非中心）を掃引し、種別遷移が l2 符号・gamma しきい値で正しく切替わること。Hybrid の境界（区間中の l2 符号反転）を明示検証。
- **k 値系統差テスト（必須）**: `EspenakUmbral` vs `IauMean` で本影半径が変わり境界判定がずれることを定量化・記録（accuracy §0）。
- **ゴールデン20**（accuracy §3.4）の全球部分（種別・gamma・最大食・食分）。1900–2100 全日食で「有無・種別・最大食時刻・gamma・食分」を一括比較。
- **UTC/TT 両返し**・将来日食で `delta_t_uncertainty` が metadata に乗ること（accuracy §0）。

**許容誤差**（accuracy §2.1 / §4 層分解）:
- 最大食時刻（TT 基準）±1.5 s（gamma 最小化 solver 収束 0.05″, root_tolerance ≤ 目標の 1/10）。
- gamma: x, y 精度律速（影幾何＋暦, ≲0.49″ 合成 ≈ 1.0 s 相当）。
- 食分 ±0.0005（0.001食分 ≈ 1.9″, accuracy §2.2）。
- 中心線位置 sub-km（≲0.5 km, 最大食地点, accuracy §2.1）。
- UTC 絶対時刻は将来 ΔT/UT1 律速（accuracy §2.3）。幾何（TT）と分離して報告。許容を通すための拡大禁止（conventions §11）。

---

## 出典

- **Explanatory Supplement to the Astronomical Almanac (3rd ed.), Ch.11「Eclipses」**: 全球状況・gamma・central eclipse の条件・全球接触。該当式番号は**要確認**（一次資料＝書籍で確認, 推測しない）。
- **Meeus, J. "Astronomical Algorithms", 2nd ed., Ch.54「Eclipses」**: gamma 定義・u（=l2）・種別判定フローチャート・帯幅・継続時間の実用式。個別式番号（54.x）は**要確認**。
- **NASA / Espenak**（GSFC eclipse site / NASA TP-2006-214141, data-sources §4.1）: 種別境界・gamma 符号・path width・duration の慣習。**c_pen/c_cen/ρ_g の数値しきい値は要確認**。
- **WGS84**（conventions §4.1）: Re = a = 6 378 137.0 m（gamma・帯幅の無次元化/有次元化基準）。
- プロジェクト内: algorithms.md §0（記号表）, conventions §1/§4/§5/§8/§9/§11, accuracy §0/§2.1/§2.2/§2.3/§3/§4, numerical-policy §A2/§A5, §9（09-local.md, magnitude/obscuration・高度方位・geocentric→geodetic は §10）, ISSUE-021/023/037/022, ISSUE-010/011。
