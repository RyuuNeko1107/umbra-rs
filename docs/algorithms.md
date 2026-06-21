# 数式仕様 (algorithms)

`umbra-rs` の計算パイプライン全体の**数式・手順・符号規約**を一箇所に集約する。レビューゲート「数式レビュー」の中核対象。
conventions.md（規約）・numerical-policy.md（数値方針）・accuracy.md（誤差バジェット）と対で読む。

> 状態: ドラフト（Milestone 0）。各式は出典（章・式番号）を必須とし、**一次資料で未確認の式番号・符号は「要確認」**を残す。
> セクション本体は `docs/algorithms/` に分割。本ファイルは記号表・規約・索引・レビュー観点を持つ正本。

---

## 0. 記法・記号表（全セクション共通コントラクト）

**この表が記号の単一の正本**。各セクションはこの記号・単位・符号に従う。逸脱が必要なら本表に追記してから使う。

### 一般記法
- ベクトルは太字小文字 `r`、単位ベクトルは `r̂`。スカラは細字。
- フレームは下付き/注記で示す（例: `r_GCRS`）。frame: ICRS/GCRS/CIRS/TIRS/ITRS/FundamentalPlane（conventions §5）。
- 角度は内部ラジアン（conventions §2）。距離は km、ただし**ベッセル長は地球赤道半径 Re で無次元化**（conventions §1/§4.1, Re = WGS84 a = 6 378 137.0 m）。
- 時刻は型付き（conventions §6）。`t` = J2000 からの **TT ユリウス世紀**（= (JD−2451545.0)/36525、JulianDate2 から生成: numerical-policy §A1）。
  - **時間基準の使い分け（重要・混同注意）**: `T` = **ユリウス千年**（=t/10、/365250）。**VSOP87 は T（千年）**、**ELP/MPP02 と歳差章動・一般式は t（世紀）** を引数に取る（暦の原著基準が異なる）。各セクションは引数が T か t かを式ごとに明記する。
  - `τ` は**光行時間**（§3, 秒）と**フィット正規化時刻**（§7, [−1,1] 無次元）の二用途。§7 では正規化時刻に限定し局所定義する（衝突回避）。

### 時刻系
| 記号 | 意味 |
|---|---|
| ΔT | TT − UT1（秒） |
| ΔAT | TAI − UTC（閏秒, 秒） |
| TT−TAI | 32.184 s（定数, conventions §4.1） |
| ERA | 地球回転角（UT1 由来, IAU2000, `iauEra00`） |
| GAST | 見かけ恒星時。**本実装では ERA 経由 CIO ベースで構成**（conventions §5.2, D4）。分点 GST は Standard で不使用 |

### 天体位置
| 記号 | 意味 |
|---|---|
| `r_sun`, `r_moon` | 見かけ地心位置（light-time→deflection→aberration→歳差章動。conventions §5.1） |
| `R_sun`, `R_moon` | 見かけ地心距離（km） |
| s_sun, s_moon | 見かけ視半径（rad）。太陽: R_sun_phys/距離（conventions §9, R_sun_phys=696000km）。月: k·Re/距離（k: conventions §9） |

### ベッセル要素（瞬時, §6）— 長さは Re 単位
| 記号 | 意味 | 符号/規約 |
|---|---|---|
| x, y | 影軸と基本面の交点座標（基本面基底 X̂,Ŷ 上） | X̂=東向き, Ŷ=天の北の基本面射影（conventions §5, ISSUE-020） |
| d | 影軸方向の赤緯（rad） | 影軸 = 太陽→月→地球向き。向き定義をセクションで明記 |
| μ | 影軸のグリニッジ時角（rad） | **μ = GAST − α_axis**、GAST = ERA − EO（CIO 経由, §2/D4）。**μ のみ UT1（δUT1）依存**（accuracy §0(a) 脚注/§2.1L）。NASA は度・hour（D5 対応表） |
| α_axis | 影軸方向の赤経（rad, CIRS 基準） | μ 構成に使用（§2/§6） |
| l1, l2 | 半影/本影錐の基本面での半径（Re） | **l2<0 ⇒ 皆既**（金環は l2>0）。符号定義をセクションで明記 |
| f1, f2 | 半影/本影錐の半頂角（rad） | tan f1, tan f2 を保持 |

### 観測者（局地, §9, §10）
| 記号 | 意味 |
|---|---|
| φ, λ, h | 測地緯度 / 東経（正） / 楕円体高（conventions §3） |
| φ′ | 地心緯度。ρ sinφ′, ρ cosφ′ = 扁平補正済み地心動径成分（ISSUE-010/011） |
| θ | 観測者における影軸の局地時角。**採用: θ = μ − λ（東経正）**。ただし一次資料（Explanatory Supplement §11）で式番号未確証のため **要確認**。両流儀の交差検証テスト（正中 θ=0→方位180°、東地平 θ<0 等）を §9 に置く |
| ξ, η, ζ | 観測者の基本面座標（Re） |
| u, v | u = x − ξ, v = y − η（影軸と観測者の基本面相対位置）。**Meeus Ch.54 は u=ξ−x（逆符号）だが m² は不変**。本書は x−ξ で統一 |
| u′, v′ | u, v の時間微分（Re/SI秒）。接触・継続・最大食 dm/dt に使用（§9） |
| L1, L2 | 観測者の ζ 面での錐半径。L1 = l1 − ζ tan f1, L2 = l2 − ζ tan f2 |
| m, n² | m² = u² + v²（numerical-policy §A5: 最小化対象は m²）。n² = u′² + v′²（相対速度二乗, §9 接触補正） |

### 食の量（§9）
| 記号 | 意味 |
|---|---|
| magnitude | 食分（EclipseMagnitude）。部分食 (L1 − m)/(L1 + L2)（要式番号確認） |
| obscuration | 食面積比（Obscuration, 0..1）。2円交差面積（§9, ISSUE-027） |
| gamma | 影軸の地心最小距離（Re, 全球, §8）。**内部は √≥0 の非負量**。NASA 慣習の符号付き（北正/南負）で返す場合は最大食時の y 符号で付与（**符号定義は要確認**, §8） |

### 追加記号（フレーム / 微分 / 楕円体）

| 記号 | 意味 | 参照 |
|---|---|---|
| X, Y, s | CIP の座標と CIO locator | §2 (IAU2006/2000A) |
| s′ | TIO locator | §2 (極運動) |
| EO | equation of the origins（GAST = ERA − EO） | §2 |
| Δψ, Δε, ε_A | 章動（黄経・傾斜）と IAU2006 平均黄道傾斜 | §2 |
| xp, yp | 極運動成分 | §2 (EOP) |
| L, B, R | 地球日心黄経・黄緯・動径（VSOP87, 引数 T=千年） | §3 |
| β_moon | 月の地心黄緯（早期棄却, §4。光行差 β=v/c とは別） | §4 |
| a, b, f, e², N, ρ | 楕円体長半径/短半径/扁平率/離心率²/卯酉線曲率半径/地心動径 | §10 (WGS84) |
| approx_gamma | 概算 gamma（早期棄却用・フル gamma より甘い） | §4 |
| ξ′, η′ | 観測者基本面座標の時間微分（Re/SI秒） | §9 |

### 記号衝突の回避規約（重要）

同一文字が文脈で別量を指す箇所がある。各セクションは下記の規約に従い、初出で明示する:

- **u, v**: 既定は基本面相対位置 `u=x−ξ, v=y−η`（§9）。§1 の `iauDtdb` 引数（観測者東距離・赤道面距離）は `u_obs, v_obs` と表記。§5 の月→太陽単位ベクトルは `û`（ハット）で区別。
- **g**: §4 の合連続化関数は `g(t)`（rad）、§6 の影軸地心位置中間量は `g_axis`（km ベクトル）と表記。
- **β**: §3 の光行差 `β=|v|/c`、§4 の月地心黄緯 `β_moon`。
- **T / t / τ**: T=千年（VSOP）、t=世紀（ELP・一般）、τ=光行時間（§3）／フィット正規化時刻（§7, 局所定義）。

---

## 1. 主要出典（一次資料）

各式は以下のいずれかに紐づけ、章・式番号を残す。未確認は「要確認」。

- **Explanatory Supplement to the Astronomical Almanac, 3rd ed.**（日食・ベッセル要素・局地条件の主典拠）
- **Meeus, "Astronomical Algorithms", 2nd ed.** Ch.54 Eclipses（局地条件の実装式）、Ch.7/22/25 等（時刻・章動・太陽）
- **IERS Conventions (2010)** ch.5（歳差章動・フレーム・ERA）、**IAU SOFA**（`iau*` 関数; 照合基準）
- **Bretagnon & Francou (1988)** A&A 202,309（VSOP87）、**Chapront & Francou (2002)**（ELP/MPP02）
- NASA/Espenak（ベッセル要素表記・μ単位・k 慣習。第二義照合, data-sources §4）

---

## 索引（セクション分割）

| § | 内容 | ファイル | 主 Issue |
|---|---|---|---|
| 1 | 時刻系変換（UTC/TAI/TT/UT1/TDB, ΔT, JD2） | algorithms/01-time-scales.md | 004-007, 042 |
| 2 | フレーム変換（IAU2006/2000A, ERA, CIO/GAST） | algorithms/02-frames.md | 035, 039, 040 |
| 3 | 天体暦・見かけ位置（VSOP87D/ELP-MPP02, 補正チェーン） | algorithms/03-ephemeris.md | 012-015, 033-034 |
| 4 | 新月・合・候補探索（Step1-3） | algorithms/04-search.md | 016-018 |
| 5 | 影円錐幾何（Step5） | algorithms/05-shadow-cone.md | 019 |
| 6 | 基本面・瞬時ベッセル要素（Step6-7） | algorithms/06-besselian.md | 020-021, 037 |
| 7 | ベッセル多項式（Step8） | algorithms/07-bessel-polynomial.md | 022 |
| 8 | 全球条件（Step9: P1/U1/最大/U4/P4, gamma, 種別） | algorithms/08-global.md | 023 |
| 9 | 局地条件（Step10: 射影/C1-C4/最大食/食分・食面積/高度方位/可視性） | algorithms/09-local.md | 024-028 |
| 10 | 観測者・楕円体（geodetic→ITRS, ρsinφ′/ρcosφ′） | algorithms/10-observer.md | 010-011 |
| 11 | 部分食域（partial_limit GeoPolygon・半影限界/rise-set 曲線・**設計ドラフト/未実装**） | algorithms/11-path-partial-domain.md | 045 |

---

## レビュー観点（数式レビューゲート）

- 各式に出典（章・式番号）。未確認は「要確認」を残し断定しない。
- **符号規約**: 時角 H/θ（東経正適応）、μ=GAST−α、l2<0=皆既、方位北0東回り、影軸向き — 一次資料と照合し本表と一致。
- **単位**: ベッセル長=Re、角度=rad、NASA μ=度・hour 対応表（D5）。
- **時刻系/フレーム**: TT基準のベッセル要素、μ のみ UT1 依存（accuracy §0(a) 脚注）。CIO 統一（D4）。
- **数値方針整合**: 級数和/微分/光行時間/フィット/求根は numerical-policy に従う。
- 誤差の層分解（accuracy §4）が式の段階に対応していること。
