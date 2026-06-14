# ISSUE-013: Analytical Sun ephemeris（VSOP87D、地球日心→太陽地心の反転）

- crate: umbra-ephemeris
- 依存: ISSUE-012（Ephemeris trait）, ISSUE-033（VSOP87 係数生成パイプライン）, ISSUE-036（JPL DE テストオラクル: **最小 SPK reader=M2 必須で打切り次数決定**, フル版 DE 差分=M10 最終ゲート）, ISSUE-047（暫定オラクル戦略: M2 暫定ゲート = Mock+SOFA+NASA。**幾何層は §3.1-2 独立オラクルで別途検証**）, umbra-core（Vector3 / Radians / AstronomicalUnits / TdbInstant）
- モード(tdd-workflow): standard（数式実装。係数は ISSUE-033 が strict で供給。本 Issue は級数評価ロジックの正しさをテスト駆動で詰める standard）

## 目的
VSOP87D 太陽暦を実装する。VSOP87D は地球の日心黄道座標（L, B, R, of date）を与えるため、これを反転して太陽の地心位置を得る。`Ephemeris` trait の `Body::Sun` 系を `AnalyticalEphemeris` の一部として担当する。

- VSOP87D の L（黄経）・B（黄緯）・R（動径 AU）級数評価（Horner で T のべき、各項は振幅×cos(位相+振動数×T)）。
- 地球日心 → 太陽地心の反転（ベクトル符号反転）。
- ICRS / EclipticOfDate 両 `EphemerisFrame` への対応（VSOP87D は黄道 of date が native）。
- 速度供給（解析微分 or 対称差分。architecture §4）。

## 非目的
- 月暦（ISSUE-014）。
- 光行時間・歳差章動・光行差（ISSUE-015。本 Issue は「補正前の幾何的地心位置」まで）。
- VSOP87 原データ取込・打切り・packed 化・checksum（ISSUE-033）。
- 惑星位置（v1.0 不要。Earth のみ使用）。

## 公開インターフェース
- 公開境界は ISSUE-012 の `Ephemeris` trait 経由（`AnalyticalEphemeris::state(Body::Sun, t, Origin::Geocenter, frame)`）。
- 太陽専用の追加 pub 型は作らない（暦は統合）。内部関数（非公開）:
  - `fn earth_heliocentric_ecliptic_of_date(t: TdbInstant) -> (Radians /*L*/, Radians /*B*/, AstronomicalUnits /*R*/)`
  - `fn sun_geocentric(t, frame) -> StateVector`（反転＋フレーム変換）。

## 数式・アルゴリズムの出典
- **VSOP87**: Bretagnon & Francou (1988), *A&A* 202, 309。級数形式: 各座標変数 `s ∈ {L, B, R}` を `s = Σ_{α=0..5} T^α Σ_k A_{k} cos(B_{k} + C_{k} T)`（T = ユリウス千年 from J2000 TDB）。VSOP87**D** = 黄道・平均分点 of date・球面（heliocentric, spherical）。data-sources §2.1/§5。
- **前提（版D・地球系列。B4(c)）**: 本実装が消費するのは **VSOP87 版D の地球系列（地球の日心黄道座標。地球中心＝Earth であり EMB ではない）**である。**EMB（地球–月重心）主系列を取り違えると 6.4″/月オーダーの系統誤差**になる（地球–EMB 差は月軌道による月次振動）。版/body の取り違えは**太陽地心＝地球日心の符号反転テストでは検出できない**（反転は絶対値の取り違えを打ち消す）。よって ISSUE-033 が生成時に版D・地球系列を明示チェックし、本実装はロード時にその版識別子・body をメタから再検査する（不一致は失敗）。
- 変数定義・係数列（A=振幅, B=位相 rad, C=振動数 rad/千年）は ISSUE-033 が生成する packed テーブルから読む。係数の出典・checksum は generated 側 NOTICE に従う。
- 時間引数 T: `T = (JD_TDB − 2451545.0) / 365250`（ユリウス千年）。出典 Bretagnon & Francou 同上。
- 太陽地心位置の反転: 太陽の地心方向ベクトル = −（地球の日心位置ベクトル）。SOFA 参照: `iauEpv00`（地球位置）に相当する概念だが本実装は VSOP87D 系列を直接評価（SOFA は使わず純Rust。SOFA は検証参照のみ ※要確認）。
- 黄道 of date → ICRS への変換は ISSUE-035 のフレーム変換（黄道傾斜・frame bias）を利用。本 Issue では `EclipticOfDate` を native とし、`Icrs` 要求時は ISSUE-035 へ委譲。

## 単位 / 時刻系 / 座標系
- 単位: L, B = ラジアン（内部）、R = AU（内部）→ trait 境界で km（× 149597870.7、conventions §4 境界変換）。
- 時刻系: TDB（VSOP87 は TDB 基準。TT との差は ≲2ms で太陽位置 ≲0.001″、許容。差異の扱いを metadata に記録）。
- 座標系: native = 黄道・平均分点 of date（EclipticOfDate）、要求に応じ ICRS（ISSUE-035 経由）。原点 = Geocenter（反転後）/ Heliocentric（内部中間）。

## アルゴリズム概要
1. T（ユリウス千年, TDB）を算出。
2. 地球の L, B, R を各べき α=0..5 の係数群で評価（Horner で T、各項は cos 加算）。正規化は L→[0,2π)（conventions §2 用途別）。
3. 球面（L,B,R）→ 日心黄道直交（of date）へ。
4. 太陽地心 = 符号反転。
5. frame=Icrs なら ISSUE-035 のフレーム変換で ICRS へ。
6. velocity: 既定は対称差分（Δ は accuracy.md §3.3 手順でテスト決定）。解析微分実装時は係数の d/dT を別経路で。供給方式名を metadata に。

## 受け入れテスト
- **二段ゲート（D / ISSUE-047, accuracy.md §3.3）**: 暦確定オラクルを二段構成にする。M10 まで盲目にはしない。
  - **M2 ゲート**: **最小 DE reader（ISSUE-036, 太陽/月のみ・M2 必須）で打切り次数を DE 実値により決定**する（暦層の第一義オラクルを M2 から使用）。併せて MockEphemeris（幾何足場）＋ SOFA 参照値（検証参照のみ・移植禁止）＋ NASA 公開値、**および幾何層の独立オラクル（accuracy.md §3.1-2, NASA ベッセル要素表 / 独立第二実装）**で検証。保証値化はしない（実測前は保証しない方針）。
  - **M10 最終ゲート（ISSUE-036 フル版）**: 下記フル版 DE 差分で M2 項数を再検証し最終確定。
- **L3 天体位置 / DE 差分（第一義オラクル=M10 最終ゲート, accuracy.md §3.1）**: feature `jpl`（ISSUE-036）の `JplEphemeris`（DE440）を基準に、1900–2100 を多数サンプリングし太陽地心方向の角度差を測定。**残差 0.05″ 級**を確認（§3.3 手順: 項を寄与順に並べ累積残差 0.05″ を切る最小項数を採用、ISSUE-033 と往復）。
- L1: 既知の単発時刻（J2000.0 等）で R≈1 AU・L が暦表値域に入ることのサニティ。
- 反転テスト: |地球日心方向 + 太陽地心方向| ≈ 0（同一直線・逆符号）を厳密に確認。
- フレーム整合: EclipticOfDate→Icrs→EclipticOfDate のラウンドトリップが ≲1mas（ISSUE-035 のテストと連携）。
- 速度テスト: 対称差分速度と解析微分（実装すれば）の一致、または DE 速度との比較を ≲（角速度許容）で確認。
- 基準値はすべて DE オラクルから動的取得。実装側へハードコードしない（conventions §11）。

## 許容誤差
- 太陽 地心見かけ位置バジェット = **0.20″**（accuracy.md §2.1 RSS 配分のうち太陽分）。本 Issue（補正前幾何位置）はそのうち **打切り残差 0.05″ 級**（§2.4）を担保し、残りは ISSUE-015（光行差等 0.10″ 共有）・ISSUE-035（0.05″）へ。
- 数値根拠: 感度 1″≈2s（相対角速度 0.5″/s）。太陽 0.20″ は最大食 ±1.5s 目標への RSS 寄与。打切り 0.05″ は DE 差分の実測下限狙い。

## 実装メモ
- 係数テーブルは ISSUE-033 の packed 形式（f64 配列＋オフセット表）を `include_bytes!` 等で組込み。手書き配列禁止（architecture §11）。
- 多数項の cos 評価が律速。Horner（T のべき）＋必要なら sin/cos の和角漸化（要確認: 精度劣化に注意、まずは素直な libm cos）。
- L の正規化と連続性: 暦評価では [0,2π) でよいが、合 solver へ渡る経路（別 crate）では連続関数化が必要（conventions §2）。本層は素の値を返し、連続化は呼び出し側。
- TDB/TT 同一視の判断と max_residual_arcsec を metadata に記録（accuracy.md §2.4）。
