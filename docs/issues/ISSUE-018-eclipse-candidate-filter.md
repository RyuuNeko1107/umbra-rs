# ISSUE-018: Eclipse candidate filter（早期棄却・偽陽性可/偽陰性不可）

- crate: umbra-eclipse
- 依存: ISSUE-017（合の精解時刻と角距離）, ISSUE-015（見かけ地心位置・距離）, ISSUE-012（Ephemeris）, ISSUE-002（Radians, acos クランプ前提）, umbra-core（TtInstant）
- モード(tdd-workflow): standard（粗フィルタ層。閾値を保守側に取れば最終結果は変わらない。ただし「偽陰性ゼロ」は search の正しさを担保するため網羅性テストは厳密。standard）

## 目的
合の各候補（ISSUE-017）について、日食が**起こりうるか**を安価な幾何量で早期判定し、明らかに非該当の朔を棄却する（architecture §3「日食可能性 早期棄却（偽陽性可・偽陰性不可）」）。
- 判定量: 月の地心黄緯（β）、月-太陽角距離、太陽・月の視半径、影軸と地球中心の距離（gamma 概算）。
- **偽陽性は可**（後段ベッセル/全球判定で確定棄却, ISSUE-021/023）、**偽陰性は不可**（実際の日食を落とさない。探索は偽陰性ゼロ＝accuracy.md §3.4）。

## 非目的
- 日食種別の確定（皆既/金環/部分/ハイブリッド）= ISSUE-023。本 Issue は「日食の可能性あり/なし」の2値。
- ベッセル要素の生成（ISSUE-021）。判定は概算幾何量で行い、フル精度は不要。
- gamma の精密値（ISSUE-023 の全球判定）。本 Issue は概算 gamma で粗判定。

## 公開インターフェース
crate 内部 API（`pub(crate)`）。判定根拠を保持して偽陽性/偽陰性のデバッグを可能にする。

```rust
pub(crate) struct EclipsePossibility {
    pub possible: bool,
    pub conjunction: Conjunction,        // ISSUE-017
    pub min_separation: Radians,         // 合付近の最小角距離（概算）
    pub sum_semidiameters: Radians,      // 太陽視半径 + 月視半径（概算）
    pub approx_gamma: f64,               // 影軸-地球中心 概算（単位 Re, ISSUE-021 と整合）
    pub reason: PossibilityReason,       // 棄却/採用の理由（デバッグ・テスト用）
}

pub(crate) enum PossibilityReason {
    PossibleEclipse,
    LatitudeTooHigh,         // 月黄緯が大きく食帯が地球を外す
    SeparationTooLarge,      // 角距離 > 視半径和 + マージン
    ShadowAxisMissesEarth,   // |gamma| > 閾値（地球半径 + 半影マージン）
}

pub(crate) fn assess_eclipse_possibility(
    eph: &impl Ephemeris,
    conjunction: &Conjunction,
    config: &EngineConfig,   // lunar/solar radius model（k 値, conventions §9）を参照
) -> Result<EclipsePossibility, EclipseError>;
```

- `config` から `LunarRadiusModel`/`SolarRadiusModel`（conventions §9）を取り、視半径を概算。**偽陰性回避のためマージンは保守側（部分食=半影基準の最大値）に取る**。

## 数式・アルゴリズムの出典
- **日食可能条件（ecliptic limits）**: Meeus, *Astronomical Algorithms* (2nd ed.), **Ch.54「Eclipses」** の食限（eclipse limits）。月の黄緯 β と昇交点離角に基づく「日食が起こりうる β の上限」。Meeus Ch.54 では太陽が交点付近（≈±18.5°: 必ず起こる限界 / ±15.4°: 必ず起こらない限界）にあるかで判定する近似式が与えられる。
- **角距離判定**: 合付近の月-太陽最小角距離 `Δ_min` < 視半径和 `(s_sun + s_moon)` + 地球半径分の視差マージンなら、地球上のどこかで部分食。Meeus Ch.54 の食限の幾何的等価。
- **概算 gamma（影軸-地球中心距離, 単位 Re）**: 月→太陽方向の影軸が地球中心から外れる距離。`|gamma| ≲ 1 + l1`（半影外半径 l1, ISSUE-021）なら半影が地球に触れる＝部分食以上。出典: Espenak/NASA ベッセル要素定義（gamma = 最小影軸距離, 単位 Re。ISSUE-021 §出典と共通）。本 Issue は概算 gamma で `|gamma| < 閾値` を粗判定。
- **視半径**: 太陽 `s_sun = asin(R_sun / d_sun)`（R_sun=696000km, conventions §9）、月 `s_moon = asin(k·Re / d_moon)`（k は conventions §9 の選択値）。距離は ISSUE-015 の地心距離。

## 単位 / 時刻系 / 座標系
- 時刻系: 合時刻 TtInstant（ISSUE-017）。
- 角度: ラジアン。角距離は `acos(clamp(û_m·û_s, -1, 1))`（accuracy.md §2.2 クランプ必須）。
- 単位: gamma は **Re 無次元**（conventions §1, ISSUE-021 と統一）。視半径はラジアン。
- 座標系: 見かけ地心（ISSUE-015）。gamma 概算の影軸は地心ベクトルから構成（厳密な FundamentalPlane 基底は ISSUE-020、本 Issue は概算で可）。

## アルゴリズム概要
1. 合時刻で月・太陽の見かけ地心方向・距離を取得（ISSUE-015）。
2. 視半径 `s_sun, s_moon` を config の半径モデルで算出。
3. 月黄緯 β を評価し、Meeus Ch.54 食限と比較（β が「必ず起こらない限界」超なら即棄却 = `LatitudeTooHigh`）。
4. 角距離 `Δ`（合付近の最小, 合時刻で代表 or 微小区間で最小化）を計算。`Δ > s_sun + s_moon + 地球視差マージン + 安全マージン` なら `SeparationTooLarge`。
5. 概算 gamma を計算し `|gamma| > 1 + l1_概算 + 安全マージン` なら `ShadowAxisMissesEarth`。
6. いずれの棄却にも当たらなければ `possible=true`（`PossibleEclipse`）。
- 数値安定性: acos クランプ。全マージンは**偽陰性ゼロを保証する側（広め）**に固定し根拠コメント化（magic number 禁止, conventions §11）。グレーゾーン（食限ぎりぎり）は必ず `possible=true`（偽陽性を許容）。

## 受け入れテスト
accuracy.md テストレベル **L5 前段（網羅性）** ＋ L8（プロパティ: 偽陰性ゼロ）。
- **偽陰性ゼロ網羅テスト（最重要）**: NASA 5千年カタログ（data-sources §4.1）の全日食朔（1900–2100）が `possible=true` になること。1件でも `false` なら fail。基準は fixtures から取得（conventions §11）。
- **実マージン余裕の統計出力テスト（D6, ISSUE-029 連携）**: ISSUE-029 の **NASA 全朔網羅**（日食朔・非日食朔とも）に対し、各朔で「棄却境界までの余裕」（角距離マージン残・gamma マージン残）を計算し、最小余裕・分布（ヒストグラム/分位）を統計出力する。最小余裕が常に正（偽陰性ゼロが余裕付きで成立）であること、および D6 導出式の各項（視差 0.95° / 黄緯速度×ずれ / 概算暦誤差）が実データで妥当な余裕を持つことを確認する。
- **偽陽性の許容上限テスト**: 日食でない朔のうち `possible=true` になる割合（無駄な後段計算量）を測定。許容は性能指標（精度ではない）。0% は要求しない。
- **MockEphemeris 人工ケース（accuracy.md §3.1）**: 影が地球を完全に外す配置 → `ShadowAxisMissesEarth`。明確な部分食/皆既/金環配置 → `possible=true`。境界（影縁が地球縁にちょうど接する）→ `possible=true`（保守判定）。
- 食限境界テスト（L8）: β を食限近傍で掃引し、棄却→採用の遷移が「真の食限より外側で起きない」（偽陰性なし）こと。
- k 値選択（conventions §9）の影響: `EspenakUmbral`/`IauMean` で視半径が変わるが、本 Issue は部分食=半影基準（大きい方）で判定するため棄却境界が皆既/金環の k 差に影響されないこと。

## 許容誤差
- 本 Issue は精度バジェット（accuracy.md §2.1）に**寄与しない**（後段が確定）。担保は**偽陰性ゼロ**。
- **マージン設計（D6 確定: 偽陰性ゼロのマージン導出式）**: マージンは下記 3 項の和で導出し、各項を根拠コメント化する（偽陰性ゼロを保証する側＝広めに固定）。
  - **(1) 月地平視差 ≈ 0.95°**: 半影が地球に触れる条件は月の地平視差を含むため必須（落とすと偽陰性）。
  - **(2) 月最大黄緯速度 × 合↔最大食ずれ時間**: 最大食は合（ISSUE-017）からずれるため、「合付近の最小角距離」を合時刻値で代用する際の取りこぼしを、月の最大黄緯速度（dβ/dt 上限）× 合と最大食の時間ずれ で上乗せする。
  - **(3) 概算暦誤差上限**: 本フィルタが使う概算暦・概算 gamma の誤差上限（角距離・gamma それぞれ）を見積もり、その分を必ず上乗せ。
  - すなわち `margin ≳ π_moon(≈0.95°) + (dβ/dt)_max × Δt_(合↔最大食) + ε_ephemeris_approx`。誤差を隠さず明示（accuracy.md §0）。
- 偽陽性率の目標は性能上の緩い上限（例: 後段コストが許容範囲）であり、精度ゲートではない。

## 実装メモ
- 「合付近の最小角距離」は合時刻の値で代用可だが、最大食は合からずれるため**微小区間（±数十分）で最小化** or 保守マージンで吸収する。偽陰性回避のため後者（広いマージン）が安全。
- 概算 gamma は ISSUE-021 のフル gamma と一致する必要はない。**フル gamma より必ず甘い（採用寄り）**ことをプロパティで保証。
- マージンに使う「地球視差」は月の地平視差（≈0.95°）相当。半影が地球に触れる条件は月の視差を含むため、ここを落とすと偽陰性になる（重要）。
- レビュー重点: 偽陰性ゼロ（マージンが全概算誤差＋視差を包含）、グレーゾーンは必ず採用、k 差が部分食判定に漏れ込まないこと。
