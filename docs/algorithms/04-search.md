# §4 新月・合・候補探索（Step1-3）

> 正本: `algorithms.md §0`（記号表）・`conventions.md §2`（角度正規化・連続化）・`§6`（時刻系）・`numerical-policy.md §A5`（求根=Brent・連続化）・`accuracy.md §3.4`（偽陰性ゼロ・マージン D6）。本セクションはこれらに**厳密準拠**する。
> 状態: ドラフト（Milestone 0）。一次資料で未確認の式番号・符号は「**要確認**」を残し、推測で式番号を書かない。
> 関連 Issue: ISSUE-016（新月候補生成）, ISSUE-017（合 solver）, ISSUE-018（早期棄却フィルタ）。連携: ISSUE-029（全期間スキャンの実余裕統計）。

---

## 目的と入出力（単位・時刻系・座標系・フレーム）

期間 `UtcRange` を「新月（朔）」単位に分解し、各朔で**地心合**を精解し、日食が起こりうる候補のみを後段（§5–§9）へ渡す3段パイプライン。最重要要件は **偽陰性ゼロ（朔・日食を1件も落とさない）**（探索は偽陰性ゼロ＝accuracy.md §3.4 / §検索戦略, ISSUE-016/018）。

| Step | 内容 | 入力 | 出力 |
|---|---|---|---|
| Step1 | 朔の概算＋検索窓付与（ISSUE-016） | `TimeRange<UtcInstant>` | `Vec<NewMoonCandidate>`（`approx_tt`, `search_window: TimeInterval<TtInstant>`, `lunation_number`） |
| Step2 | 地心合の精解（ISSUE-017） | `NewMoonCandidate`, `ConjunctionKind` | `Conjunction`（`time_tt`, `separation`） |
| Step3 | 日食可能性の早期棄却（ISSUE-018） | `Conjunction`, `EngineConfig` | `EclipsePossibility`（`possible`, `approx_gamma`, `reason`） |

- **時刻系**: 概算・窓・合時刻すべて `TtInstant`（conventions §6）。期間入力は UTC、境界で TT へ変換（§1, ISSUE-006/007）。
- **角度**: ラジアン。月相 `Δλ`・赤経差 `Δα` は**連続化専用関数**で扱う（`[0,2π)`/signed と混在禁止, conventions §2）。
- **座標系/フレーム**: Step1 は黄道 of date 近似（概算で可）。Step2 黄経合は黄道 of date、赤経合は CIRS（§2/§3）。Step3 は見かけ地心（§3, ISSUE-015）。
- **暦精度**: Step1 は精度不問（簡易評価可, ISSUE-016）。Step2/Step3 は Standard 暦（§3）を使うが、Step3 の判定量は概算 gamma で粗く取り **保守側**（広め）に倒す。

---

## 記号（algorithms.md §0 参照。本節固有の補助記号のみ定義）

§0「時刻系」「天体位置」「食の量（gamma）」を正本とする。本節固有:

| 記号 | 意味 | 単位 |
|---|---|---|
| `k_syn` | lunation index（Meeus Ch.49 の `k`。新月で整数, 起点 k=0 = 2000-01-06 新月） | 無次元 |
| `T_syn` | 平均朔望月 = 29.530588 日（平均値, Meeus Ch.49） | 日 |
| `Δ_win` | 検索窓半幅（`search_window = [approx_tt − Δ_win, approx_tt + Δ_win]`） | 日 |
| `Δλ` | 月-太陽 黄経差 = `λ_moon − λ_sun`（合は 0 mod 2π） | rad |
| `Δα` | 月-太陽 赤経差 = `α_moon − α_sun`（赤経合は 0 mod 2π） | rad |
| `g(t)` | 連続化した合関数（`Δλ` または `Δα` を窓中心基準で unwrap したもの） | rad |
| `β` | 月の地心黄緯（食限判定） | rad |
| `dβ/dt` | 月黄緯の時間変化率（その上限を `(dβ/dt)_max`） | rad/s |
| `Δt_(合↔最大食)` | 合時刻と最大食時刻のずれ（D6 第2項） | s |
| `π_moon` | 月の地平視差 ≈ 0.95°（D6 第1項。地平視差 `arcsin(Re/r_moon)`） | rad |
| `s_sun`, `s_moon` | 太陽・月の見かけ視半径（§0, §3 E12/E13） | rad |
| `ε_eph` | 概算暦・概算 gamma の誤差上限（D6 第3項。角距離・gamma 各々） | rad / Re |
| `Δ_min` | 合付近の月-太陽 最小角距離（早期棄却の主量） | rad |

> 記号衝突注意: 本節の `β`（月黄緯）は §3 の VSOP/ELP の `B`（地球日心黄緯）・`U`（月黄緯）とは別ラベル（月地心黄緯を指す）。`g(t)` は ISSUE-017 の連続化関数。

---

## 数式（番号付き・各式に出典）

### 4.1 Step1: 朔の概算と検索窓（ISSUE-016）

**(S1) 平均新月時刻（Meeus Ch.49, 式 49.1）**
出典: Meeus, *Astronomical Algorithms* (2nd ed.), Ch.49「Phases of the Moon」式 49.1（JDE of mean phase）。

```
JDE_mean(k_syn) = 2451550.09766 + 29.530588861 · k_syn + 補正多項式(T)      (T = k_syn/1236.85)
```

- `k_syn` は新月で**整数**。`k_syn = 0` は 2000-01-06 の新月（Meeus Ch.49）。lunation index の起点を固定しコメント化（ISSUE-016 採番起点）。
- 周期補正項（太陽・月の平均近点角等, Meeus Ch.49）は**窓を狭めたい場合のみ**採用、必須でない（概算で足りる, ISSUE-016 非目的）。

**(S2) lunation index の初期推定（Meeus Ch.49）**
出典: Meeus Ch.49。
```
k_syn ≈ (year_decimal − 2000.0) · 12.3685
```
- 期間 start 直前の `k_syn` を逆算し、整数へ丸めて `k_0` とする。end を超えるまで `k_syn` を 1 ずつ増やす。

**(S3) 検索窓（ISSUE-016）**
出典: ISSUE-016（窓マージン= 平均朔の最大ずれ + 安全係数）。
```
approx_tt    = TT(JDE_mean(k_syn))                              # 概算朔時刻（TT）
search_window = [ approx_tt − Δ_win,  approx_tt + Δ_win ]       # 合 solver のブラケット窓
```
- **`Δ_win` の根拠**: 実際の朔は平均朔から **±約14時間（≈±0.6 日）ずれる**（Meeus Ch.49 解説）。これに安全係数を掛けた値（既定 **±1 日**）を採用し、magic number 化せず根拠コメントを付す（conventions §11, ISSUE-016）。**Δ_win が偽陰性ゼロを担保する唯一の砦**であり、テストで下限検証（下記）。
- 期間端で窓がはみ出す朔も取りこぼさないため、`k_0 − 1` 朔／`end + 1` 朔まで生成してから窓が範囲と交差するものだけを残す（ISSUE-016 端処理）。

### 4.2 Step2: 地心合の精解（ISSUE-017）

**(S4) 合関数（連続化前）**
出典: Meeus, Ch.54「Eclipses」冒頭（日食判定の伝統定義は**黄経合** `λ_moon = λ_sun`）/ Ch.49（朔）。
```
ConjunctionKind::EclipticLongitude :  f(t) = Δλ(t) = λ_moon(t) − λ_sun(t)
ConjunctionKind::RightAscension    :  f(t) = Δα(t) = α_moon(t) − α_sun(t)
```
- **黄経合を既定**、赤経合は照合用（NASA 表記突合, ISSUE-017）。`λ`/`α` は見かけ地心（§3, ISSUE-015）。

**(S5) 連続化（unwrap, conventions §2）**
出典: conventions §2「±π 折返しを除いた連続関数」/ ISSUE-017（NR 連続化前処理）。
```
g(t) = unwrap( f(t) ; 基準 = f(t_center) )       # 窓中心 t_center の値を基準に ±2π を引き去る
```
- 窓内では月の黄経角速度（≈13.2°/日）≫ 太陽（≈0.99°/日）のため `Δλ` は**単調増加**（Meeus Ch.47/25）。窓が小さい（±1 日, S3）ので折返しは高々1回。`g(t)` は合付近で連続・単調。

**(S6) ブラケット → Brent（numerical-policy §A5）**
出典: numerical-policy §A5（求根=Brent, Newton 単独禁止）/ ISSUE-008/017。
```
粗走査: 窓 [t0,t1] を等間隔サンプルし g(t_i)·g(t_{i+1}) < 0 を検出 → ブラケット [t_i, t_{i+1}]
精解  : Brent で g(t)=0 を root_tolerance まで（Newton 単独禁止, conventions §11）
```
- 粗走査刻みは月運動（13.2°/日）から「窓内で `g` が単調 1 通過」を保証する値とし根拠コメント化（刻み過大＝根の見落とし＝偽陰性, ISSUE-017）。
- `root_tolerance` は目標（最大食 ±1.5s）の 1/10 ＝ **≤0.15 s 相当**（accuracy.md §2.1, EngineConfig `root_tolerance_seconds`=0.01 s 既定が満たす, numerical-policy §A5）。
- 独立変数は窓内オフセット（日 or 秒, `JulianDate2` 差分で橋渡し, numerical-policy §A1/ISSUE-008）。

**(S7) 合時刻の角距離（早期棄却の入力）**
出典: ISSUE-017 / accuracy.md §2.2（acos クランプ）。
```
separation = acos( clamp( û_moon · û_sun , −1, 1 ) )        # 合時刻の月-太陽 角距離
```

### 4.3 Step3: 日食可能性の早期棄却（ISSUE-018）

**(S8) 視半径（§3 E12/E13, conventions §9）**
```
s_sun  = asin( R_sun_phys / R_sun )       # R_sun_phys = 696000 km（conventions §9）
s_moon = asin( k · Re / R_moon )          # k = 部分食=半影基準（大きい方）を採用（偽陰性回避, ISSUE-018）
```
- 棄却境界が皆既/金環の `k` 差（EspenakUmbral vs IauMean, conventions §9）に汚染されないよう、**部分食=半影基準の最大値**で判定（ISSUE-018）。

**(S9) 食限（ecliptic limits, Meeus Ch.54）**
出典: Meeus, Ch.54「Eclipses」食限。
- 月黄緯 `|β|` が「必ず起こらない限界」を超えれば即棄却 = `LatitudeTooHigh`。Meeus Ch.54 は太陽の交点離角に基づく近似限界（≈±18.5°=必ず起こる / ≈±15.4°=必ず起こらない）を与える。**要確認**: 採用する限界式の正確な係数・式番号は Meeus Ch.54 一次で確定（ここでは概念のみ。偽陰性回避のため「必ず起こらない限界」を保守側に取る）。

**(S10) 角距離棄却（Meeus Ch.54 食限の幾何等価）**
出典: Meeus Ch.54。
```
Δ_min > s_sun + s_moon + π_moon + margin   ⇒  SeparationTooLarge
```
- `Δ_min` は合付近の最小角距離（合時刻値で代用可だが最大食ずれを margin で吸収, ISSUE-018）。`margin` は D6（下記「許容誤差」）。

**(S11) 概算 gamma（影軸-地心距離, Re）**
出典: Espenak/NASA ベッセル定義（gamma = 影軸の地心最小距離, Re。§0, ISSUE-021 と共通）。
```
approx_gamma = (影軸の地心最小距離) / Re          # 概算（厳密 FundamentalPlane 基底=§6 は使わない, ISSUE-018）
|approx_gamma| > 1 + l1_概算 + margin   ⇒  ShadowAxisMissesEarth
```
- `approx_gamma` は §6 のフル gamma と一致不要。**フル gamma より必ず甘い（採用寄り）**ことをプロパティで保証（ISSUE-018）。
- いずれの棄却にも当たらなければ `possible = true`（`PossibleEclipse`）。グレーゾーンは必ず採用（偽陽性を許容, ISSUE-018）。

---

## 手順（実装順・数値注意）

1. **Step1（ISSUE-016）**: 期間 [start,end]（UTC）→ TT 境界変換（§1）。S2 で `k_0` 推定 → S1 で平均朔列生成（`k_0−1` … `end+1` 朔）。S3 で `approx_tt`・窓付与。窓が範囲と交差するものを時系列 `Vec` で返す。
2. **Step2（ISSUE-017）**: 各窓で S4→S5（連続化）→ S6（粗走査ブラケット → Brent）。S7 で `separation`。窓に符号変化が **0 回** なら `RootNotBracketed`（＝ Step1 のマージン不良が顕在化, ISSUE-017）、**複数回**なら粗走査刻みが粗すぎとして異常検出。
3. **Step3（ISSUE-018）**: 合時刻で見かけ位置・距離（§3）→ S8 視半径 → S9 食限 → S10 角距離 → S11 概算 gamma。全マージンは **偽陰性ゼロ側（広め）** に固定し根拠コメント化（magic number 禁止, conventions §11）。

**数値注意（横断, numerical-policy）**:
- 合関数は **連続化してから** Brent（±π 折返し除去, conventions §2 / §A5）。`[0,2π)`/signed 正規化と混在させない（用途別関数, conventions §2）。
- 求根は **Brent（要ブラケット）**。Newton 単独禁止（conventions §11 / §A5）。
- acos/asin 引数は `[-1,1]` クランプ（accuracy.md §2.2 / §A5）。
- 時刻独立変数は `JulianDate2` 差分（生 f64 JD 禁止, conventions §6/§11, numerical-policy §A1）。

---

## 境界・特異・異常系

- **窓に符号変化なし（Step2）**: `RootNotBracketed`（api-draft §3.5）。実運用で起きれば Step1 マージン不良のサイン → 偽陰性リスク直結なので回帰テストで監視（ISSUE-016/017）。
- **粗走査で符号変化が複数回**: 刻み過大。刻みを月運動から導出し直す（ISSUE-017）。
- **Brent 未収束**: `max_iterations` 不足 → `SolverDidNotConverge`（ISSUE-017）。
- **期間端・空範囲・1朔未満**: `k_0−1`/`end+1` 朔まで生成し交差判定（取りこぼし防止, ISSUE-016）。空範囲は空 `Vec`。
- **食限グレーゾーン（Step3）**: 必ず `possible = true`（偽陽性許容・偽陰性禁止, ISSUE-018）。
- **lunation_number**: 単調増加・一意（off-by-one を朔個数チェックで検出, ISSUE-016）。
- **±π 折返し（連続化）**: 窓が境界をまたぐケースで `g(t)` が連続・単調であること（ISSUE-017 L1）。

---

## 検証（基準値の出典。実装へ値コピー禁止 = conventions §11）

accuracy.md §3.1/§3.4、ISSUE-016/017/018 受入テスト準拠。基準値は fixtures/DE/NASA から動的取得（ハードコード禁止, conventions §11）。

- **偽陰性ゼロ網羅（最重要, L5 前段 / L8）**: NASA 5千年カタログ（data-sources §4.1, 第二義）の全日食朔（1900–2100）が、(i) Step1 のいずれかの窓に内包、(ii) Step3 で `possible=true` になること。1件でも外れれば fail。
- **朔個数チェック**: 100 年で `≈ 100·12.3685 ≈ 1237 件`（off-by-one 検出, ISSUE-016）。`lunation_number` 単調増加・一意。
- **窓マージン下限テスト**: 平均朔（S1）と精解朔（S6 or DE 由来の真の朔）の差が常に `Δ_win` 未満（Δ_win 不足＝偽陰性の直接検出, ISSUE-016）。
- **合 solver 検証**: DE440 で同一合関数を解いた朔時刻 vs 解析暦（Standard）の差を 1900–2100 で測定（第一義, 暦残差を §4 層分解, accuracy.md §3.1）。MockEphemeris（線形に動く月・太陽）の解析解と厳密比較（オラクル＝人工配置の数式, 実装非コピー, ISSUE-017）。連続化テスト（窓が ±π をまたぐ, L1）。収束 root_tolerance ≤0.15s（accuracy.md §2.1）。
- **実マージン余裕の統計出力（D6, ISSUE-029 連携）**: NASA 全朔（日食朔・非日食朔）で各朔の棄却境界までの余裕（角距離マージン残・gamma マージン残）を計算し、最小余裕・分布を統計出力。最小余裕が常に正（余裕付きで偽陰性ゼロ）であること、D6 各項（π_moon / 黄緯速度×ずれ / 概算暦誤差）が実データで妥当な余裕を持つことを確認（ISSUE-018/029）。
- **MockEphemeris 人工ケース（Step3）**: 影が地球を完全に外す → `ShadowAxisMissesEarth`。明確な部分/皆既/金環 → `possible=true`。境界（影縁が地球縁に接する）→ `possible=true`（保守, ISSUE-018）。
- **k 値非汚染（Step3）**: `EspenakUmbral`/`IauMean` で視半径が変わっても、部分食=半影基準（大きい方）判定なので棄却境界が皆既/金環の k 差に影響されないこと（ISSUE-018）。

---

## 許容誤差（偽陰性ゼロのマージン導出式 D6）

本節（§4）は精度バジェット（accuracy.md §2.1）に**直接寄与しない**（後段 §6–§9 が確定）。担保すべきは **網羅性＝偽陰性ゼロ**（ISSUE-016/018）。

**D6 マージン導出式（accuracy.md §3.4 / ISSUE-018 確定）** — Step3 の角距離・gamma マージンは下記 3 項の和で導出し、各項を根拠コメント化（偽陰性ゼロ側＝広めに固定, 誤差を隠さない accuracy.md §0）:

```
margin ≳ π_moon( ≈ 0.95° )
       + (dβ/dt)_max × Δt_(合↔最大食)
       + ε_eph
```

- **(1) 月地平視差 `π_moon ≈ 0.95°`**: 半影が地球に触れる条件は月の地平視差を含むため必須（落とすと偽陰性, ISSUE-018）。`π_moon = arcsin(Re / r_moon)`、近地点で最大 ≈0.95° → **最大値を採用**（保守側, conventions §1 Re 基準）。
- **(2) `(dβ/dt)_max × Δt_(合↔最大食)`**: 最大食は合（Step2）からずれる。合付近の最小角距離を「合時刻値」で代用する際の取りこぼしを、**月黄緯の最大変化率 × 合↔最大食の時間ずれ** で上乗せ。`(dβ/dt)_max` は月運動の上限（Meeus Ch.47, 要確認: 具体上限値は M2 実測で確定）、`Δt_(合↔最大食)` は数十分オーダー（ISSUE-018 実装メモ）。
- **(3) 概算暦誤差上限 `ε_eph`**: 本フィルタが使う概算暦・概算 gamma の誤差上限（角距離・gamma それぞれ）を見積もり必ず上乗せ（ISSUE-018, 要確認: 上限値は M2 実測）。

> Step1 の窓半幅 `Δ_win`（S3）も同種の保守設計（平均朔ずれ ±0.6 日 + 安全係数）。窓を不必要に広げると後段コストが増えるため、広すぎ/狭すぎは**性能指標**（精度ではない）として ISSUE-029 で実余裕を統計監視する（accuracy.md §3.4）。

許容は「通すためだけに拡大しない」（conventions §11）。

---

## 出典

- Meeus, *Astronomical Algorithms* (2nd ed.): Ch.49「Phases of the Moon」式 49.1（平均朔 JDE・k）, Ch.54「Eclipses」（合定義・食限）, Ch.47（月運動）, Ch.25（太陽運動）。
- Espenak/NASA（GSFC eclipse / EclipseWise）: gamma = 影軸の地心最小距離（Re）, ベッセル定義（§6 と共通, data-sources §4.1）。**第二義照合**（慣習差を accuracy.md に記録, accuracy.md §3.1）。
- conventions §2/§6/§9/§11, accuracy.md §0/§2.1/§2.2/§3.1/§3.4, numerical-policy §A1/§A5, algorithms.md §0。
- §3（見かけ地心位置 E12/E13 視半径）, §1（時刻系変換）, §6（フル gamma・ベッセル要素）。
- 関連 Issue: ISSUE-016, ISSUE-017, ISSUE-018（D6）, ISSUE-008（Brent）, ISSUE-029（実余裕統計）。
- **要確認**: Meeus Ch.54 食限の正確な係数・式番号（S9）。D6 (2)(3) の `(dβ/dt)_max`・`ε_eph` の数値上限（M2 実測, accuracy.md §3.4）。
