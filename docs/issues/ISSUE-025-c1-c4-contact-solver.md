# ISSUE-025: C1/C4 contact solver（外接・粗走査→接触候補区間→Brent）

- crate: umbra-eclipse
- 依存: ISSUE-008（Brent 求根）, ISSUE-024（観測者→基本面射影）, ISSUE-021（瞬時ベッセル要素）, ISSUE-037（`BesselianSource`／直接評価器・Standard 局地既定）, ISSUE-023（全球分類＝探索窓の元）, ISSUE-007（UT1/ΔT・UTC↔TT）, ISSUE-001（規約）
- モード(tdd-workflow): strict（局地接触時刻の正本。conventions §11 の無条件 Newton 禁止・ブラケット必須、accuracy.md §2 接触 ±2s を律速。安全性=収束保証に関わるため strict）

## 目的
観測地点の **部分食 開始 C1 / 終了 C4（外接）**、および中心食地点の **C2/C3（内接）** を、ベッセル基本面上で求解する（conventions §8, architecture §7）。
- 幾何条件（conventions §8）: C1/C4 は中心間距離 `m = √((ξ−x)² + (η−y)²) = L1`（外接、半影縁=部分食限界）。C2/C3 は `m = |L2|`（内接、本影/反本影縁）。
- 手順: **粗走査で `g(t)=m(t)−L(t)` の符号変化区間を検出 → ブラケット → Brent（ISSUE-008）**。無条件 Newton 禁止（conventions §11）。
- 接触時刻は **UTC と TT の両方**を返す（accuracy.md §0, conventions §6）。

## 非目的
- 最大食時刻・食分（ISSUE-026/027）。
- 高度・方位・可視性（ISSUE-028）。ただし各接触の `LocalContact` は alt/az/position_angle/visible フィールドを持つため、それらの計算は ISSUE-028 を呼んで埋める（本 issue は時刻求解が主、フィールド充足は連携）。
- 全球 P1/U1/U4/P4（architecture §7、別 issue）。本 issue は局地接触。
- 瞬時ベッセル要素の生成（ISSUE-021）と供給源（ISSUE-037）。本 issue は `BesselianSource::at(t)` を時刻関数として呼ぶ。

## 公開インターフェース
api-draft §3.4 `LocalContact` / `LocalCircumstances`、architecture §7 に整合:
```rust
/// 局地接触（api-draft §3.4）。UTC+TT 両方（accuracy.md §0）。
#[derive(Clone, Copy, Debug)]
pub struct LocalContact {
    pub time_utc: UtcInstant, pub time_tt: TtInstant,
    pub sun_altitude: Degrees, pub sun_azimuth: Degrees,
    pub position_angle: Degrees, pub visible: bool,
}

/// 局地接触求解の結果（部分食地点では c2/c3 = None）。
pub(crate) fn solve_local_contacts<B: BesselianSource>(
    source: &B,
    observer_itrs_re: Position<Itrs>,
    search: TimeInterval<TtInstant>,    // 全球接触から絞った探索窓
    config: &EngineConfig,              // root_tolerance_seconds（目標の1/10、accuracy §2.1）
) -> Result<LocalContactSet, EclipseError>;

pub struct LocalContactSet {
    pub c1: Option<LocalContact>, pub c2: Option<LocalContact>,
    pub c3: Option<LocalContact>, pub c4: Option<LocalContact>,
}
```
- `c2/c3` は中心食域外（部分食地点）で `None`（api-draft §6 未確定事項「部分食地点で c2/c3 が None」を型で表現）。
- `BesselianSource`（api-draft §3.3）。既定 Standard は `InstantaneousEvaluator`（ISSUE-037、直接評価・fit 誤差ゼロ、architecture §6.1）。
- `EclipseError`（api-draft §3.5）: `RootNotBracketed`, `SolverDidNotConverge`, `DegenerateGeometry` 等。

## 数式・アルゴリズムの出典
- **Explanatory Supplement to the Astronomical Almanac, §11 "Eclipses"** 局地予報: 接触は `(ξ−x)² + (η−y)² = L²`（L=L1 で外接 C1/C4、L=|L2| で内接 C2/C3）。L1/L2 は基本面上の影半径（観測者 ζ・tan f1/f2 で高さ補正: `L = l − ζ·tan f`）。
- **Meeus, Astronomical Algorithms 2nd ed., Ch.54** 局地状況の接触時刻計算（式 54.x: `u=ξ−x`, `v=η−y`, `n²=u'²+v'²`, 接触の時刻補正）。Meeus は反復補正法だが、本 issue は **粗走査+Brent で堅牢化**（conventions §11、Newton 単独禁止）。Meeus 式は初期ブラケット範囲の見積りに使用。
- 求根対象関数: `g(t) = ((ξ(t)−x(t))² + (η(t)−y(t))²) − L(t)²` = `m²(t) − L(t)²`（差の連続関数）。出典式は実装コメントに章・式番号転記。
- **D2（最小化対象の統一）**: 本 Issue の `g = m² − L²` 系と ISSUE-026 の最大食解を **m²(=u²+v², u=ξ−x, v=η−y) 基準に統一**する（中心線尖点で m=√(…) を直接最小化すると微分が特異になる回避）。接触は m²=L² の求根、最大食は m² の極小（ISSUE-026, dm/dt=0 ⇔ d(m²)/dt=0 の求根）として整合させる。
- ζ による影半径補正 `L1 = l1 − ζ·tan f1`, `L2 = l2 − ζ·tan f2`（Explanatory Supplement、観測者高さ＝基本面からの距離）。tan f1/f2 は瞬時要素（api-draft §3.3）。

## 単位 / 時刻系 / 座標系
- 入力: 観測者 ITRS（Re、ISSUE-011/024）、`BesselianSource`（要素 x,y,l1,l2 は Re、d,μ は `Radians`）。
- 求根独立変数: TT 秒（or 日）。`JulianDate2` 差分を f64 へ橋渡し（conventions §6、ISSUE-008 の境界）。
- 出力: `LocalContact`（`time_utc`/`time_tt` 併記。TT→UTC は ISSUE-007 の `TimeScales::tt_to_utc`）。
- 座標系: FundamentalPlane（conventions §5）。

## アルゴリズム概要
1. 探索窓（全球接触 ±マージン）を一定刻みで粗走査し `g(t)` を評価。
2. `g` の符号変化区間を検出 → 各区間を `Bracket` 化（ISSUE-008、`f(a)·f(b)<0` 保証）。
3. 各ブラケットで Brent 求根（root_tolerance は目標の 1/10 以下＝接触 ±2s に対し ≤0.2s、accuracy §2.1）。
4. 求まった時刻の前後で `g` の符号（外側→内側＝C1/C2、内側→外側＝C3/C4）で接触種別を割当て。
5. 中心食条件: その地点で `m_min < |L2|`（本影が観測者を覆う）なら C2/C3 が存在、さもなくば `None`（部分食地点）。
6. 各接触時刻で alt/az/position_angle/visible（ISSUE-028）を計算し充足。TT→UTC 変換。
- 数値安定性: `g` は ξ,η を ±π 折返し除去後の連続関数で（conventions §2, §24）。L²（2乗形）で扱い acos 不要。`m≈L` の接触付近は勾配が小さいが Brent の二分法フォールバックで堅牢（無条件 Newton 禁止・conventions §11）。
- 部分食地点で c2/c3=None: m_min と |L2| の比較で明示的に分岐（api-draft §6 の None 設計）。掠める接触（grazing, m_min≈L1）は粗走査刻みを十分細かく取り見落とし防止（偽陰性不可・architecture §3）。

## 受け入れテスト
accuracy.md テストレベル **L6（局地条件）**。基準値は実装へコピー禁止（出典明記）。
- 地点分類（accuracy.md L6）: 中心線上 / 中心線付近 / 北限・南限 / 部分食域 / 限界（grazing）/ 可視域外 / 日の出中 / 日没中 / 標高差。比較項目: **C1, C2, C3, C4**。
- オラクル（第二義・整合チェック、絶対基準にしない・accuracy.md §3.1）: **NASA 5千年カタログ / USNO の地点別接触時刻**（data-sources §4）。ΔT・k 慣習を揃えた上で比較（conventions §9、Espenak 慣習へ切替）。fixtures に数値転記・出典/取得日明記（ISSUE-029）。
- 部分食地点: `c2 == None && c3 == None`、`c1`/`c4` は存在。
- 中心食地点: `c2`/`c3` が存在し `c1 < c2 < max < c3 < c4`（順序プロパティ、L8）。
- 異常系: ブラケット不成立 → `RootNotBracketed`。掠め食で粗走査刻みを粗くすると見落とす → 刻み感度テストで偽陰性ガード。
- 接触時刻は TT と UTC の両方を検証（accuracy.md §0）。UTC は ΔT 経由（ISSUE-007）。

## 許容誤差
accuracy.md §2「局地接触時刻 ±2 秒（幾何分・TT基準）」、§0(b)「UTC 絶対は ΔT/UT1 律速」から:
- **C1〜C4 時刻（TT 基準・幾何相対）: ±2 s**。root_tolerance はその 1/10 以下（≤0.2 s、§2.1）。角度感度 0.5″/s → 0.2s≈0.1″ で solver バジェット 0.05″ に収める設定を推奨。
- **UTC 絶対時刻**: ΔT/UT1 予測律速（accuracy.md §0(b)/§2.3）。過去・近傍は EOP 実測で <0.1s、将来（〜2100）は数秒に増大。`CalculationMetadata.delta_t_uncertainty_seconds` を併記。**TT 基準の値を一級保持**（conventions §6）。
- 根拠: 接触は「計算律速の幾何分」と「予測律速の UTC 分」を分離（accuracy.md §0）。許容を通すためだけに拡大しない（conventions §11）。

## 実装メモ
- 全球接触（P1/U1/U4/P4）で探索窓を絞り、無駄な粗走査を避ける。窓が取れない（その地点で日食なし）なら全接触 None を返し、上位で可視性判定（ISSUE-028）。
- 接触種別の割当ては `g` の符号遷移方向で機械的に。外接（L1）と内接（L2）を別 solver パスで（L が異なる）。
- TT→UTC は ISSUE-007 の `TimeScales`。将来日食は不確実性帯を必ず metadata に。
- Newton 単独禁止（conventions §11）。Meeus の反復補正は初期窓見積りのみに使い、確定は Brent。
- レビュー重点: ブラケット維持、粗走査刻みの偽陰性ガード（grazing）、c2/c3=None 分岐、L=l−ζ·tan f の高さ補正、UTC+TT 両方、時角符号（ISSUE-024 と整合）。
