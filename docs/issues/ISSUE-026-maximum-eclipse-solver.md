# ISSUE-026: Maximum eclipse solver（距離最小化・最大食時刻・食分）

- crate: umbra-eclipse
- 依存: ISSUE-008（Brent 求根: D2 正式手法 dm/dt=0 の求根に使用）, ISSUE-009（1次元最小化: D2 で粗ブラケット用に降格・併用）, ISSUE-024（観測者→基本面射影）, ISSUE-021（瞬時ベッセル要素）, ISSUE-037（`BesselianSource`）, ISSUE-027（食分・食面積／最大時点の値算出）, ISSUE-007（UTC↔TT）, ISSUE-001
- モード(tdd-workflow): strict（最大食時刻 ±1〜2s・食分 ±0.0005 を律速。最小化の数値安定性が精度を決めるため strict）

## 目的
観測地点で、太陽-月の **基本面上中心間距離 `m(t)=√((ξ−x)²+(η−y)²)`** が最小となる最大食の時刻・その時点の食分を求める（conventions §8 最大食、architecture §7）。
- **手法（D2 確定）**: 正式手法は **dm/dt=0 の求根（Brent, ISSUE-008）**。最小化対象は **m²(=u²+v², u=ξ−x, v=η−y) に統一**し（ISSUE-025 の g=m²−L² 系と揃える。中心線尖点で √ を直接最小化すると微分特異になる回避）、d(m²)/dt=0 を Brent で求根する。距離最小化（黄金分割・ISSUE-009）は **粗ブラケット（最小付近の括り出し）用に降格／併用**し、確定は dm/dt=0 求根で行う（架構 §12）。
- 最大食は部分食地点でも必ず存在（`LocalCircumstances.maximum: LocalContact` は `Option` でない・api-draft §3.4）。
- 時刻は UTC+TT 両方（accuracy.md §0）。食分・食面積は ISSUE-027 を最大時点で評価。

## 非目的
- 接触 C1〜C4（ISSUE-025）。
- 食分・食面積の式そのもの（ISSUE-027）。本 issue は最大時刻を出し、その時刻で ISSUE-027 を呼ぶ。
- 全球の最大食（`GreatestEclipse`・architecture §7）。本 issue は**観測地点局地**の最大食。
- 高度・方位・可視性（ISSUE-028。最大時点の alt/az は呼んで充足）。

## 公開インターフェース
api-draft §3.4 `LocalCircumstances.maximum` / `maximum_altitude`、architecture §7 に整合:
```rust
pub(crate) struct LocalMaximum {
    pub time_tt: TtInstant, pub time_utc: UtcInstant,
    pub min_separation: f64,        // m_min（Re）
    pub magnitude: EclipseMagnitude,
    pub obscuration: Obscuration,
    pub contact: LocalContact,      // 最大時点の alt/az/PA/visible（ISSUE-028）
}

pub(crate) fn solve_local_maximum<B: BesselianSource>(
    source: &B,
    observer_itrs_re: Position<Itrs>,
    search: TimeInterval<TtInstant>,   // 全球最大 ±マージン
    config: &EngineConfig,
) -> Result<LocalMaximum, EclipseError>;
```
- `EclipseMagnitude`（皆既で 1 超を許容）, `Obscuration`（0..1）は api-draft §3.4 / ISSUE-027。
- `min_separation` は食分/食面積算出（ISSUE-027）の入力。

## 数式・アルゴリズムの出典
- **Explanatory Supplement to the Astronomical Almanac, §11 "Eclipses"**: 最大食 = 投影面上の太陽中心-月中心距離最小。Meeus の `m` 最小化と等価。
- **Meeus, Astronomical Algorithms 2nd ed., Ch.54**: `u=ξ−x, v=η−y`、最大食は `m=√(u²+v²)` 最小。**D2: 最小化対象は m²=u²+v² に統一**。極小条件は `d(m²)/dt = 2(u·u'+v·v') = 0`、すなわち **dm/dt=0 ⇔ (u·u'+v·v')=0**。本 issue は **この (u·u'+v·v')=0 を Brent（ISSUE-008）で求根するのを正式手法**とする。Meeus の `t = −(u·u'+v·v')/(u'²+v'²)` 線形補正反復（式 54.x）は**初期推定（粗ブラケット）にのみ使用**し、確定は dm/dt=0 求根（無条件 Newton 回避・conventions §11、堅牢性）。黄金分割（ISSUE-009）は粗ブラケット併用に降格。
- 食分（最大時点で ISSUE-027 へ）: `magnitude = (L1 − m)/(L1 + L2)`（外接縁から内接縁への食い込み割合。Explanatory Supplement / Meeus Ch.54）。皆既/金環で 1 を超える/満たす扱いは ISSUE-027。出典式は実装コメントに章・式番号。

## 単位 / 時刻系 / 座標系
- 入力: 観測者 ITRS（Re）、`BesselianSource`（x,y,l1,l2=Re, d,μ=`Radians`）。
- 最小化独立変数: TT 秒（or 日）。`JulianDate2` 差分を f64 へ橋渡し（conventions §6）。
- 出力: `time_tt`/`time_utc`（TT→UTC は ISSUE-007）、`min_separation`（Re）、食分・食面積（ISSUE-027）。
- 座標系: FundamentalPlane（conventions §5）。

## アルゴリズム概要
1. 探索窓を粗走査し `m²(t)` の最小付近の 3 点ブラケット（左>中<右）を検出（粗ブラケット用に ISSUE-009 黄金分割を併用可）。
2. **（D2 正式手法）** `d(m²)/dt = 2(u·u'+v·v') = 0` を **Brent 求根（ISSUE-008）** で確定＝最大食時刻。距離最小化は粗ブラケットのみに使用。
3. 最大時点の瞬時要素で食分・食面積を ISSUE-027 で算出。
4. 最大時点の alt/az/PA/visible を ISSUE-028 で算出。
5. TT→UTC 変換、`LocalMaximum` 返却。
- 数値安定性（D2）: `m²` は ξ,η の連続化（±π 折返し除去・conventions §2）後に扱う。正式手法 dm/dt=0 求根（Brent, ISSUE-008）は m² の極小を導関数ゼロ点として捉え、√ の中心線尖点の微分特異を回避（m² は滑らか）。粗ブラケットの黄金分割（ISSUE-009）は導関数不要・無条件収束。Meeus 線形補正は初期推定のみ。**皆既帯では `m` が |L2| を下回る平底（mがほぼ一定の平坦域）があり得る**: この場合 dm/dt≈0 が区間で成り立ち**最小は一意とは限らない（平底）**。平底ケースでは「最大食時刻」を平底区間の代表点（例: 中央 or 全球最大食に最も近い点）として規約で定義し、解の非一意を許容する（旧「最小は一意」断定は平底前提に修正）。掠め食で最小が浅い場合も粗走査刻みで括れること。
- 部分食地点でも最大食は存在（`maximum` は非 Option）。可視域外でも幾何最大は計算し、可視性は ISSUE-028 が別判定。

## 受け入れテスト
accuracy.md テストレベル **L6（局地条件）**。基準値は実装へコピー禁止。
- 地点分類（accuracy.md L6）: 中心線上（m_min≈0）/ 付近 / 北南限 / 部分食域 / 限界 / 可視域外 / 標高差。比較項目: **最大食時刻、食分**（食面積は ISSUE-027 連携）。
- オラクル（第二義・整合・accuracy.md §3.1）: NASA 5千年カタログ / USNO の地点別最大食時刻・食分（data-sources §4）。ΔT・k 慣習を揃える（conventions §9）。fixtures 転記・出典/取得日明記（ISSUE-029）。
- 中心線上: `m_min ≈ 0`、食分 ≥ 1（皆既）or ≈1（金環は m_min≈|L2|）。
- **皆既帯の平底 fixture（D2・受入テスト）**: `m` が |L2| を下回り平坦（dm/dt≈0 が区間で成立）になる皆既帯中心付近の配置を MockEphemeris で構成し、(a) dm/dt=0 求根が平底区間で破綻せず収束する、(b) 「最大食時刻」が規約（平底区間の代表点）で決定的に定まる、(c) 食分が皆既値で安定、を確認。平底では最小が一意でない前提を明示テスト化。
- 最小性プロパティ（L8）: `m(max) ≤ m(max±δ)`（両側で増加 **または平底では等しい**）。平底以外では狭義に両側増加、平底では規約代表点が区間内に入ること。最大食は C1〜C4（ISSUE-025）の間に入る（`c1 < max < c4`）。
- 異常系: 探索窓に極小がない（その地点で食なし）→ `DegenerateGeometry` or 上位で可視性 NotVisible。
- 時刻は TT/UTC 両方（accuracy.md §0）。

## 許容誤差
accuracy.md §2.1「最大食時刻 ±1〜2 秒（TT基準）」、§2.2「食分 ±0.0005」、§0(b) UTC 律速:
- **最大食時刻（TT 基準）: ±1〜2 s**。最小化収束は目標の 1/10（時刻 ≤0.15〜0.2s、§2.1 root_tolerance 規約）。
- **食分: ±0.0005**（ISSUE-027 の精度。相対位置 ≲1″、§2.2）。
- **UTC 絶対時刻**: ΔT/UT1 律速（§0(b)/§2.3）。不確実性帯を metadata に。TT 基準を一級保持。
- 根拠: 最大食時刻は相対形状（計算律速・§0(a)）。最小化器の収束ガードで余裕を確保し、許容を通すためだけに拡大しない（conventions §11）。

## 実装メモ
- 最小化の独立変数連続化（ξ,η の折返し除去）を ISSUE-024 の出力で担保。求根／最小化器へ渡す `m²(t)`・`dm/dt` は連続前提。
- **D2: 確定は dm/dt=0 の Brent 求根（ISSUE-008）**。Meeus 線形補正反復・黄金分割（ISSUE-009）は粗ブラケット（初期推定）限定に降格。最小化対象は m²=u²+v² に統一（ISSUE-025 と整合, architecture §12）。
- 食分・食面積は最大時点の瞬時要素を ISSUE-027 へ渡して算出（重複評価を避け要素をキャッシュ）。
- 中心食帯では magnitude の 1 超（皆既）/ ≈1（金環）を ISSUE-027 が境界処理。
- レビュー重点: **D2（dm/dt=0 求根が正式手法、m² 統一、距離最小化は粗ブラケットに降格、皆既帯平底の非一意扱い）**、粗ブラケットの括り（3点）、連続化、Newton 単独不使用、UTC+TT 両方、最大食が接触の内側に入る検証。
