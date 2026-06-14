# ISSUE-037: Instantaneous Besselian evaluator（直接: 各反復で暦再評価・BesselianSource 実装・Standard 局地既定）

- crate: umbra-eclipse
- 依存: ISSUE-021（瞬時要素カーネル `besselian_elements_at`・定義）, ISSUE-019（影円錐）, ISSUE-020（基本面基底）, ISSUE-015（見かけ地心位置）, ISSUE-012（Ephemeris）, ISSUE-007（恒星時/UT1）, umbra-core（TtInstant, TimeInterval）
- モード(tdd-workflow): strict（`BesselianSource` は局地 solver の値供給契約（公開 trait, api-draft §3.3）。Standard 局地の**既定供給源**であり fit 誤差ゼロを担保する層。trait 契約は SemVer 境界。strict）

## 目的
任意の TT 時刻について、**各時刻で暦を再評価して**瞬時ベッセル要素を供給する直接評価層 `InstantaneousEvaluator` を実装し、`BesselianSource` trait を満たす（api-draft §3.3, architecture §6.1）。
- 局地 solver の各反復で `at(time_tt)` を呼ぶと、その都度 ISSUE-015→019→020→021 のパイプラインを回して瞬時要素を返す（**多項式 fit 誤差ゼロ**, architecture §6.1）。
- **Standard プロファイルの局地計算の既定供給源**（architecture §6.1, api-draft §3.2）。
- **責務**: ISSUE-021（1時刻の値計算カーネル＋定義）を「任意 TT で連続供給する層」。ISSUE-022（多項式）と対をなす供給源の片割れ。

## 非目的
- 瞬時要素の定義・1時刻の値計算式（= ISSUE-021。本 Issue は ISSUE-021 を任意 TT で繰り返し呼ぶ層）。
- 多項式近似・係数 fit（= ISSUE-022。経路/エクスポート用の別供給源）。
- 接触/最大食の求解そのもの（局地 solver。本 Issue は値供給に徹する）。
- 暦・見かけ位置の実装（ISSUE-012/015 を利用）。

## 公開インターフェース
api-draft §3.3 / architecture §6.1 に準拠（公開 trait ＋実装型）。

```rust
pub trait BesselianSource {
    /// 任意 TT における瞬時ベッセル要素を返す。
    fn at(&self, time: TtInstant) -> Result<InstantaneousBesselianElements, EclipseError>;
    fn fit_interval(&self) -> TimeInterval<TtInstant>;
}

/// 各時刻で暦を直接再評価（fit 誤差ゼロ）。Standard 局地の既定供給源。
pub struct InstantaneousEvaluator<'e, E: Ephemeris> {
    /* eph: &'e E, astrometry: AstrometryOptions(Standard), config 参照,
       time_scales(恒星時/UT1), 有効区間 */
}

impl<'e, E: Ephemeris> InstantaneousEvaluator<'e, E> {
    pub fn new(eph: &'e E, time_scales: &'e TimeScales, config: &'e EngineConfig,
               interval: TimeInterval<TtInstant>) -> Self;
}

impl<'e, E: Ephemeris> BesselianSource for InstantaneousEvaluator<'e, E> {
    fn at(&self, time: TtInstant) -> Result<InstantaneousBesselianElements, EclipseError>;
    fn fit_interval(&self) -> TimeInterval<TtInstant>;  // = 構築時 interval（近傍に制限）
}
```

- `at` は内部で ISSUE-015（見かけ位置）→ ISSUE-019（影円錐）→ ISSUE-020（基底）→ ISSUE-021（`besselian_elements_at`）を順に呼ぶ。
- `fit_interval` は「この供給源が妥当な区間」（多項式と API を揃えるためのもの。直接評価では暦対応範囲・日食近傍を返す）。
- Standard 局地計算（`EclipseEngine::local_circumstances`, api-draft §3.2）は本型を既定で使用。

## 数式・アルゴリズムの出典
- 計算式は ISSUE-021（Explanatory Supplement Ch.11 / Meeus Ch.54 / NASA Espenak）に全面依拠。本 Issue は**式を持たず**、ISSUE-021 カーネルを任意 TT で評価するオーケストレーション層。
- 見かけ位置補正は ISSUE-015（Standard: light_time/aberration/precession_nutation 必須 ON, accuracy.md §1）。
- 恒星時（μ 用 GAST）は ISSUE-007/035（UT1, IAU2006/2000A）。
- 出典の章・式番号は ISSUE-021 を参照（重複記述しないが、実装コメントは ISSUE-021 へリンク, conventions §10）。

## 単位 / 時刻系 / 座標系
- ISSUE-021 と完全一致: 要素は **TT 基準**（conventions §6）、**FundamentalPlane**（conventions §5）、x,y,l1,l2 は **Re 無次元**、d,μ ラジアン。
- `at(time)` の入力は TtInstant。内部で暦評価は TDB（TT≈TDB 許容, metadata 記録）、μ 恒星時は UT1。
- 供給源差し替えで単位/座標系が変わってはならない（多項式 ISSUE-022 と同一契約）。

## アルゴリズム概要
1. `at(time_tt)`: 入力時刻が `fit_interval` 外なら範囲チェック（呼出側の責務だが防御的に検証）。
2. ISSUE-015 で太陽・月の見かけ地心位置・距離を取得（Standard 補正）。
3. ISSUE-019 で影円錐（config の k/半径モデル, conventions §9）。
4. ISSUE-020 で基本面基底（d, α_z）。
5. ISSUE-007/035 で GAST（μ 用）。
6. ISSUE-021 `besselian_elements_at` を呼び `InstantaneousBesselianElements` を構成。
7. 各呼出は独立・副作用なし（純関数的）。**fit 誤差ゼロ**（毎回暦を引く）。
- 数値安定性: ISSUE-015/019/020/021 の安定性に従う。失敗は `EclipseError`（EphemerisUnavailable/DegenerateGeometry 等）を伝播。

## 受け入れテスト
accuracy.md テストレベル **L4（ベッセル）** ＋ **L7 サブテスト（直接 vs 多項式残差, accuracy.md §3.2）**。
- **ISSUE-021 一致テスト**: `InstantaneousEvaluator::at(t)` の出力が、同一入力で ISSUE-021 `besselian_elements_at` を直接組んだ値と**完全一致**（オーケストレーションが式を変えない保証）。
- **直接 vs 多項式 残差（L7, accuracy.md §3.2, architecture §6.1）**: 同一日食で本 Issue（直接）と ISSUE-022（多項式）の x,y,l1,l2 を fit 区間で比較し残差を実測 → profile 毎の採用根拠に。直接側を基準（fit 誤差ゼロ）とする。
- **NASA ベッセル値比較（第二義, data-sources §4.1）**: 任意 TT（最大食・接触付近）の評価値を NASA 公開値と整合（k/ΔT 慣習を揃える, accuracy.md §3.1）。
- **MockEphemeris（accuracy.md §3.1）**: 人工配置で `at` を複数時刻評価し、要素の時間変化（x,y が直線的に最接近へ）が物理的に妥当。
- `fit_interval` 外の時刻、暦範囲外 → 適切な `EclipseError`。
- trait オブジェクト互換: `&dyn BesselianSource` 経由でも ISSUE-022 と差し替え可能（局地 solver が供給源にジェネリック, architecture §6.1）。

## 許容誤差
- accuracy.md §2.1: **多項式 fit 誤差 0**（直接評価のため）。よって本供給源の許容は ISSUE-021（影幾何誤差）＋ ISSUE-015（暦・補正誤差）に帰着し、本 Issue 自体は追加誤差を入れない。
- L7 サブテスト（accuracy.md §3.2）: 「直接 vs 多項式」残差を本 Issue を基準に測定。多項式側 <0.10″ 相当（accuracy.md §2.1 多項式 fit 行）を満たすか判定。
- Standard 局地の既定が本 Issue（直接）である根拠＝fit 誤差ゼロ（精度最優先, architecture §6.1）。性能が問題な場合のみ多項式へ（経路用途）。

## 実装メモ
- **責務分担（品質基準）**: ISSUE-021=定義＋1時刻カーネル / 本 Issue(037)=任意 TT 直接供給層 / ISSUE-022=多項式供給層。Standard 局地既定=本 Issue、経路=ISSUE-022（architecture §6.1, api-draft §3.2/§3.3）。
- 局地 solver は `impl BesselianSource` にジェネリック（or `&dyn`）。本 Issue と ISSUE-022 を無改修で差し替え可能にする（architecture §6.1）。
- 各 `at` 呼出で暦をフル評価＝計算コスト大。Standard 局地は接触5点＋反復のみなので許容（経路の多点は多項式へ）。キャッシュは入れない（純関数性・fit誤差ゼロ優先, 入れるなら要レビュー）。
- ライフタイム/所有権（`&'e E` 参照保持 vs `Arc`）は `EclipseEngine` のジェネリック構成（api-draft §6 未確定事項）と整合させる。レビューで確定。
- metadata に「直接評価（fit誤差0）」を記録（CalculationMetadata, accuracy.md §4 層分解）。
- レビュー重点: ISSUE-021 との値一致、trait 契約（多項式と差し替え可）、fit誤差ゼロの担保、単位/座標系が供給源で不変。
