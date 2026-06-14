# ISSUE-035: IAU2006歳差 + IAU2000A章動 + フレーム変換（GCRS→CIRS→TIRS→ITRS）

- crate: umbra-ephemeris（フレーム型は umbra-core 由来）
- 依存: ISSUE-012（基盤）, umbra-core（Matrix3/Vector3/Position<F>/Radians/Ut1Instant/TtInstant）, EarthOrientation（極運動・UT1-UTC。api-draft §2）
- モード(tdd-workflow): **strict**（公開仕様: フレーム連鎖の定義・Standard 必須ON・採用モデル（公開出力 = 2006/2000A。2000B は内部粗スキャン用として存置・非公開）。conventions §5 / accuracy.md §1。strict）

## 目的
地球姿勢に基づくフレーム変換連鎖を実装する: `GCRS →[frame bias + 歳差(IAU2006) + 章動(IAU2000A)]→ CIRS →[ERA(UT1)]→ TIRS →[極運動]→ ITRS`（conventions §5）。見かけ位置（ISSUE-015）とベッセル要素生成（umbra-eclipse）が利用する。**公開出力（Standard）では IAU2006/2000A を必須**、**内部粗スキャン（非公開・候補棄却専用）のみ IAU2000B を許容**（architecture §5）。

## 非目的
- 暦・見かけ位置補正の適用順序管理（ISSUE-015 が本 Issue の回転を呼ぶ）。
- EOP/ΔT データ読込（IersEopData / DeltaTModel。別担当。本 Issue は値を受け取り回転を作る）。
- 黄道座標 of date ↔ ICRS（ISSUE-013/014 が利用するが、黄道傾斜・frame bias の提供は本 Issue が担う）。

## 公開インターフェース（※署名はレビュー確定）
- `pub fn gcrs_to_cirs_matrix(time_tt: TtInstant) -> Matrix3`（frame bias + IAU2006 歳差 + IAU2000A 章動。CIO ベース）。
- `pub fn era_rotation(time_ut1: Ut1Instant) -> Matrix3`（CIRS→TIRS。Earth Rotation Angle）。
- `pub fn polar_motion_matrix(xp: Radians, yp: Radians, time_tt: TtInstant) -> Matrix3`（TIRS→ITRS。TIO locator s′ 含む）。
- `pub fn gcrs_to_itrs_matrix(time_tt, time_ut1, xp, yp) -> Matrix3`（連鎖合成）。
- 内部粗スキャン用（非公開・任意/非既定）: `gcrs_to_cirs_matrix_2000b(time_tt) -> Matrix3`（IAU2000B 章動。内部粗スキャン用として存置。公開出力では使わない）。
- 黄道傾斜 / 平均黄道 of date 提供: `pub fn mean_obliquity_iau2006(time_tt) -> Radians`（ISSUE-013/014 の黄道→ICRS 用）。
- `Position<F>` 変換ヘルパ（型レベルでフレーム遷移を表現、conventions §5）。

## 数式・アルゴリズムの出典（SOFA 関数名まで特定）
- **IAU2006 歳差**: Capitaine, Wallace & Chapront (2003) A&A 412, 567 / IAU2006 (P03)。SOFA: `iauP06e`（歳差角）, `iauPmat06`（歳差行列）, `iauPfw06`（Fukushima-Williams 角）。IERS Conventions 2010 ch.5。
- **IAU2000A 章動**: Mathews, Herring, Buffett (2002) / IAU2000A。SOFA: `iauNut00a`（Δψ, Δε, 1365項級数）。内部粗スキャン用 IAU2000B: `iauNut00b`（77項・非公開・存置）。
- **frame bias + 歳差 + 章動（NPB）合成**: SOFA `iauPnm06a`（IAU2006/2000A の bias-precession-nutation 行列）。CIO ベース: `iauC2i06a`（celestial-to-intermediate 行列, GCRS→CIRS）。CIP 座標 X,Y + CIO locator s: `iauXys06a`。
- **ERA（Earth Rotation Angle）**: IAU2000 定義。SOFA: `iauEra00(UT1)`。ERA = 2π(0.7790572732640 + 1.00273781191135448·(JD_UT1 − 2451545.0))。
- **極運動 + TIO locator s′**: SOFA `iauPom00(xp, yp, sp)`、`iauSp00(time)` で s′。極運動 (xp,yp) は EOP（accuracy.md §2.3, data-sources §3.1）。
- **平均黄道傾斜 IAU2006**: SOFA `iauObl06`。
- 時間引数: 歳差章動は TT（ユリウス世紀 from J2000）、ERA は UT1。conventions §6。

## 単位 / 時刻系 / 座標系
- 単位: 角度 rad、回転は無次元 Matrix3（直交・右手系, conventions §5）。
- 時刻系: 歳差章動・s/s′・黄道傾斜 = TT、ERA = UT1（ΔT/UT1-UTC は EOP/DeltaTModel から、accuracy.md §5）。
- 座標系: GCRS / CIRS / TIRS / ITRS（conventions §5 表）。CIO ベース変換で統一（赤道分点ベースの恒星時は使わない。一貫性のため ※要確認: 恒星時版を内部粗スキャンに許すか）。

## アルゴリズム概要
1. GCRS→CIRS: `iauC2i06a` 相当（X,Y,s から CIP/CIO 行列）。内部粗スキャン（非公開）は 2000B 章動で X,Y を近似。
2. CIRS→TIRS: ERA(UT1) 回りの R3 回転。
3. TIRS→ITRS: 極運動行列（xp,yp,s′）。
4. 連鎖合成 `gcrs_to_itrs_matrix` を提供。逆行列（転置）も提供。
5. 黄道傾斜・frame bias を ISSUE-013/014 の黄道↔ICRS に供給。

## 受け入れテスト
- **L3 / SOFA 値突合（要確認: SOFA をテスト依存に入れてよいか。入れない場合は IERS/SOFA 公開テストベクトルを fixtures へ転記）**: 既知 TT/UT1/xp/yp に対し `gcrs_to_itrs_matrix` 要素を SOFA 参照値と比較。**残差 ~1mas 級**（accuracy.md §2.1 で 0.05″ 配分に対し実力 ~1mas で余裕）。
- 各段単体: `iauPnm06a` / `iauEra00` / `iauPom00` 相当をそれぞれ SOFA テストベクトルで検証。
- IAU2000A vs 2000B 差分が既知オーダー（~1mas 級）に収まることを確認（内部粗スキャンで 2000B を許容する根拠。公開出力は 2000A）。
- ラウンドトリップ: GCRS→ITRS→GCRS が恒等（≲数 µas）。
- DE 差分パイプライン（accuracy.md §3.1）に組み込んだとき、フレーム由来残差が層分解で 0.05″ 以下に帰属することを確認（§3.3, §4）。
- 基準値は SOFA/IERS 公開ベクトルから（実装へハードコードのコピーをしない。conventions §11）。

## 許容誤差
- 歳差章動 + フレーム = **0.05″**（accuracy.md §2.1 RSS 配分）。IAU2006/2000A 実力 ~1mas のため大幅な余裕。
- 内部粗スキャン（非公開・候補棄却専用）の IAU2000B は ~1mas 級の追加誤差だが、内部粗スキャンは候補棄却用途のため無視可（出力精度目標ではない。報告値は Standard 再計算）。
- 数値根拠: §2.1 RSS（月0.40/太陽0.20/光行差0.10/フレーム0.05/solver0.05）→ 合成 ~0.49″ ≈ 1.0s、目標 ±1.5s に余裕。

## 実装メモ
- 純Rust で SOFA を移植せず**式から実装**（SOFA は MIT 類似だが C。係数級数 nut00a は 1365項 → 生成 or 定数表。表は ISSUE-033/034 同様に provenance 付きで管理推奨。要確認: nut00a 係数の出典 IERS Conventions 2010 表をデータ化）。
- CIO ベース（X,Y,s）で統一し、赤道分点/GST 経路と混在させない（conventions §2 用途別正規化と同精神）。
- 公開出力（Standard）は IAU2006/2000A を EngineConfig で固定。内部粗スキャン（非公開・候補棄却専用）のみ 2000B 切替を許可（私的設定、公開 API には出さない）し、選択を `CalculationMetadata`（ephemeris/フレームモデル欄）へ記録。
- 相対論偏向はここではなく ISSUE-015。本 Issue は純粋に幾何回転。
- s′（TIO locator）は微小（数十 µas）だが省略せず実装（精度最優先方針）。省略するなら metadata 注記。
