# ISSUE-014: Analytical Moon ephemeris（ELP/MPP02 DE-fit 版）

- crate: umbra-ephemeris
- 依存: ISSUE-012（Ephemeris trait）, ISSUE-034（ELP/MPP02 係数生成パイプライン）, ISSUE-036（JPL DE テストオラクル: **最小 SPK reader=M2 必須で打切り次数決定**, フル版 DE 差分=M10 最終ゲート）, ISSUE-047（暫定オラクル戦略: M2 暫定ゲート = Mock+SOFA+NASA。**幾何層は §3.1-2 独立オラクルで別途検証**）, umbra-core（Vector3 / Radians / Kilometers / TdbInstant）
- モード(tdd-workflow): standard（級数評価ロジックをテスト駆動。係数生成は ISSUE-034 が strict。ライセンス由来の混入禁止は ISSUE-034 側で担保）

## 目的
ELP/MPP02（DE fit 版）月暦を実装し、`Ephemeris` trait の `Body::Moon, Origin::Geocenter` を `AnalyticalEphemeris` の一部として担当する。月の地心位置（経度・緯度・距離）を主問題級数＋摂動級数で評価する。

- ELP/MPP02 主問題（main problem）＋惑星/地球摂動級数の評価。
- 経度 V・緯度 U・距離 r の合成と地心黄道（of date）直交への変換。
- 速度供給（解析微分 or 対称差分）。

## 非目的
- 太陽暦（ISSUE-013）。
- 見かけ位置補正（光行時間・歳差章動・光行差。ISSUE-015。本 Issue は補正前の幾何地心位置まで）。
- ELP 原データ取込・打切り・packed 化・checksum・**GPL 非混入の保証**（ISSUE-034）。
- 月縁地形（ベイリービーズ等。accuracy.md §6 非保証）。

## 公開インターフェース
- 公開境界は ISSUE-012 の `Ephemeris` trait 経由（`AnalyticalEphemeris::state(Body::Moon, t, Origin::Geocenter, frame)`）。
- 月専用 pub 型は作らない（暦統合）。内部（非公開）:
  - `fn moon_geocentric_ecliptic_of_date(t: TdbInstant) -> (Radians /*V経度*/, Radians /*U緯度*/, Kilometers /*r*/)`
  - `fn moon_geocentric(t, frame) -> StateVector`（フレーム変換含む）。

## 数式・アルゴリズムの出典
- **ELP/MPP02**: Chapront & Francou (2002),「The lunar theory ELP revisited. Improvement of the main problem」, *A&A* 387, 700（および IMCCE 配布のMPP02説明）。data-sources §2.2/§0。
- 構造: 主問題級数（Delaunay 引数 D, l, l′, F の整数結合に対する正弦/余弦項）＋ 摂動級数（惑星摂動・地球扁平・潮汐・相対論・固有摂動）。経度 V、緯度 U、距離 r の3系列。
- 採用版: **MPP02 の DE fit 版**（LLR 版ではなく、Reference=JPL DE と整合。data-sources §2.2）。fit パラメータ系列の選択を metadata に固定。
- Delaunay 基本引数の時間多項式（D, l, l′, F, ζ/平均経度）の出典は Chapront-Touzé & Chapront / ELP2000-82 系の引数式に準拠（MPP02 が採用する版を ISSUE-034 が原データから抽出。係数を事実として取込）。
- **重要（ライセンス）**: 係数・引数式は **IMCCE 原データから ISSUE-034 が自前生成**。Yuk Tung Liu（ytliu0/ElpMpp02）の **GPL-3.0** 実装/係数は参照（理解・検証）のみで移植・取込禁止（data-sources §0/§2.2）。
- SOFA に月解析暦の直接対応関数はない（DE 系 `iauMoon98` 等は別系列）。検証は DE 差分で行い、SOFA 月関数は使わない ※要確認。

## 単位 / 時刻系 / 座標系
- 単位: V, U = ラジアン、r = km（ELP は km または地球半径単位 → 境界で km に統一。conventions §1）。
- 時刻系: TDB（ELP/MPP02 基準。TT 差 ≲2ms → 月位置 ≲0.03″。許容しつつ metadata に記録）。
- 座標系: native = 平均黄道・平均分点 of date（J2000 黄道版か of date 版かを ISSUE-034 の採用系列に合わせ固定。※要確認）、要求に応じ ICRS（ISSUE-035 経由）。原点 Geocenter。

## アルゴリズム概要
1. T（ユリウス世紀 from J2000 TDB。ELP は世紀単位が一般的。係数の時間単位を ISSUE-034 packed のメタに合わせる）を算出。
2. Delaunay 引数等の基本引数を多項式評価。
3. 主問題級数 + 摂動級数を V, U, r それぞれ累積（各項 = 振幅×sin/cos(引数の整数結合)）。
4. 平均経度を加算して V を完成 → 球面 (V,U,r) → 黄道 of date 直交。
5. frame=Icrs なら ISSUE-035 経由で ICRS。
6. velocity: 既定 対称差分（Δ は accuracy.md §3.3 でテスト決定）。供給方式名を metadata に。

## 受け入れテスト
- **二段ゲート（D / ISSUE-047, accuracy.md §3.3）**: 暦確定オラクルを二段構成にする。M10 まで盲目にはしない。
  - **M2 ゲート**: **最小 DE reader（ISSUE-036, 太陽/月のみ・M2 必須）で打切り次数を DE 実値により決定**する（暦層の第一義オラクルを M2 から使用）。併せて MockEphemeris（幾何足場）＋ SOFA 参照値（検証参照のみ・移植禁止）＋ NASA 公開値、**および幾何層の独立オラクル（accuracy.md §3.1-2, NASA ベッセル要素表 / 独立第二実装）**で検証。保証値化はしない。
  - **M10 最終ゲート（ISSUE-036 フル版）**: 下記フル版 DE 差分で M2 項数を再検証し最終確定。
- **L3 / DE 差分（第一義=M10 最終ゲート, accuracy.md §3.1）**: `JplEphemeris`（DE440, ISSUE-036）の月地心方向を基準に 1900–2100 多数サンプルで角度差測定。**残差 0.1″ 級**を確認。§3.3 手順で項を寄与順に並べ累積残差 0.1″ を切る最小項数を採用（ISSUE-034 と往復で打切り次数確定）。
- 距離残差: r の DE 差分が距離許容内（視半径経由で食分に効くため、相対 ≲ 数 ppm を別途ガード）。
- **距離残差→地平視差→局地接触時刻 感度（D1, 局地 topocentric バジェット）**: 距離 r を経度・緯度の角度級数と同基準で打切ると、距離残差 δr/r が **地平視差 π ≈ 3422″·(δr/r)** を介して局地接触時刻にまで効く（地心 gamma・最大食時刻 TT には現れない局地固有経路）。受入テストで δr/r を掃引し、π の変化 → 局地接触時刻の変化を実測して局地バジェット（別建て RSS）の許容内に収まることを確認。距離は角度級数と同精度では足りない場合があるため、距離側の打切り残差許容を視差感度から逆算して設定する。
- L1: J2000.0 等で V, U, r が暦表値域に入るサニティ。
- フレーム整合: EclipticOfDate↔Icrs ラウンドトリップ ≲1mas（ISSUE-035 連携）。
- 速度: 対称差分 vs 解析微分（実装すれば）一致、または DE 速度比較。
- 基準値は DE オラクルから動的取得。ハードコード禁止（conventions §11）。

## 許容誤差
- 月 地心見かけ位置バジェット = **0.40″**（accuracy.md §2.1 RSS 最大配分）。本 Issue（補正前幾何）は打切り残差 **0.1″ 級**（§2.4）を担保、残りは ISSUE-015/035 と RSS。
- 数値根拠: 月は相対角速度の主因（1″≈2s）。0.40″ 配分は最大食 ±1.5s 目標の支配項。打切り 0.1″ は DE 差分実測下限狙い（過剰打切りを避ける、§2.4 方針）。
- **距離 r の許容（D1, 局地 topocentric 経由）**: 角度（V, U）残差許容とは別に、距離残差 δr/r に独立の許容を置く。δr/r は地平視差 π ≈ 3422″·(δr/r) を介して局地接触時刻に効くため、距離を角度級数と同基準で切ると視差経由で局地時刻に効く点を明示し、許容は視差感度（局地バジェット別建て RSS, accuracy.md D1）から逆算する。

## 実装メモ
- 係数規模が大（主問題＋摂動で数千～数万項）。ISSUE-034 の packed 形式（整数引数係数テーブル＋実係数）を組込み、評価ループはキャッシュ効率を意識。
- sin/cos 評価が律速。引数の整数結合は和角漸化で高速化可だが精度劣化に注意（まず libm、ベンチ後に最適化。要確認）。
- 時間単位（世紀 vs 千年）・黄道版（J2000 vs of date）は ISSUE-034 の生成メタと厳密に一致させる。不一致は系統誤差になるため、ラウンドトリップ＋DE 差分で検出。
- 打切り次数・達成残差・採用 fit 系列を `EphemerisMetadata.version` に記録（accuracy.md §2.4）。
- GPL 実装の数値を「期待値」としてテストに貼らない（参照のみ。基準は DE）。
