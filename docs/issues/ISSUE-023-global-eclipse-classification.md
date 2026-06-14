# ISSUE-023: Global eclipse classification（種別/gamma/全球 P1U1U4P4/最大食地点/食帯幅/中心食継続）

- crate: umbra-eclipse
- 依存: ISSUE-022（BesselianPolynomial, 経路/全球用供給源）, ISSUE-037（直接供給, 最大食精解）, ISSUE-021（瞬時要素・l2 符号）, ISSUE-008（Brent）, ISSUE-009（1次元最小化, 最大食）, ISSUE-010/011（WGS84・地上射影）, ISSUE-007（UT1/恒星時, 地理座標）, umbra-core（UtcInstant/TtInstant）
- モード(tdd-workflow): strict（種別・gamma・全球接触・最大食は `GlobalCircumstances`（公開型, api-draft §3.4）の中核。種別境界（皆既/金環/ハイブリッド/非中心）と k 値選択の系統差が結果を左右。strict）

## 目的
ベッセル要素（直接 ISSUE-037 / 多項式 ISSUE-022）から日食の**全球的状況**を確定する（architecture §3/§7, api-draft §3.4）。
- **種別判定**: Partial / Annular / Total / Hybrid / NonCentralAnnular / NonCentralTotal（api-draft §3.4 `SolarEclipseKind`）。
- **gamma**（影軸の地球中心最小距離, Re）。
- **全球接触 P1/U1/U4/P4**（部分食開始/中心食開始/中心食終了/部分食終了, conventions §8）。
- **最大食地点・食帯幅・中心食継続時間**（`GreatestEclipse`, api-draft §3.4）。

## 非目的
- 局地接触 C1–C4・観測者条件（別 Issue, conventions §8）。本 Issue は全球（地球を1天体として扱う）。
- 中心線/限界線の経路サンプリング（umbra-geo。本 Issue は最大食点・帯幅・全球接触まで）。
- ベッセル要素の生成（ISSUE-021/037/022 を利用）。

## 公開インターフェース
api-draft §3.4 に準拠（公開型）。

```rust
#[non_exhaustive] #[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SolarEclipseKind { Partial, Annular, Total, Hybrid, NonCentralAnnular, NonCentralTotal }

#[derive(Clone, Debug)]
pub struct GlobalCircumstances {
    pub kind: SolarEclipseKind,
    pub partial_begin: Option<GlobalContact>, pub central_begin: Option<GlobalContact>, // P1 / U1
    pub greatest: GreatestEclipse,
    pub central_end: Option<GlobalContact>,   pub partial_end: Option<GlobalContact>,   // U4 / P4
    pub gamma: f64,    // 単位 Re
}
#[derive(Clone, Copy, Debug)] pub struct GlobalContact { pub time_utc: UtcInstant, pub time_tt: TtInstant, pub position: GeoPoint }
#[derive(Clone, Copy, Debug)]
pub struct GreatestEclipse {
    pub time_utc: UtcInstant, pub time_tt: TtInstant,
    pub position: GeoPoint, pub magnitude: EclipseMagnitude, pub obscuration: Obscuration,
    pub path_width: Option<Kilometers>, pub central_duration: Option<f64 /* s */>,
    pub sun_altitude: Degrees,
}

pub(crate) fn classify_global(
    source: &impl BesselianSource,        // ISSUE-037 直接（最大食精解）/ ISSUE-022
    time_scales: &TimeScales,             // TT↔UTC（両方返す, accuracy.md §0）
    config: &EngineConfig,                // k 値選択（conventions §9）, earth_model
) -> Result<GlobalCircumstances, EclipseError>;
```

- 接触時刻は **UTC と TT の両方**（conventions §6, accuracy.md §0）。
- gamma・帯幅は Re/km（conventions §1/§4）。

## 数式・アルゴリズムの出典
- **第一義: Explanatory Supplement to the Astronomical Almanac (3rd ed.), Ch.11**（全球状況・gamma・central eclipse の条件）。**補助: Meeus, *Astronomical Algorithms* (2nd ed.), Ch.54「Eclipses」**（gamma, u（=l2）, 種別判定の実用式）。**NASA Espenak**（種別境界・帯幅・継続時間, data-sources §4.1）。
- **gamma**: 最大食時刻における影軸の地球中心最小距離（基本面上）`gamma = min sqrt(x²+y²)`（時間最小化）。単位 Re。出典: Meeus Ch.54 / NASA（gamma 定義）。
- **種別判定（Meeus Ch.54 / NASA 境界）**:
  - `|gamma| > ~1.55`（要確認, l1 依存）→ 影軸が地球を外れる → 半影のみ → **Partial**。
  - 影軸が地球に当たる（中心食）かつ最大食時の **l2（本影縁半径, ISSUE-021, 正本 B1: 皆既で負・金環で正）** の符号:
    - l2 > 0（本影頂点が地球側手前 = 反本影が地表）→ **Annular**。
    - l2 < 0（本影が地表に到達）→ **Total**。
    - 食の経過中に l2 の符号が反転（金環⇄皆既）→ **Hybrid**（annular-total）。
  - 中心線が地球縁をかすめる（影軸は当たるが中心食条件を一部のみ満たす）→ **NonCentral**（NonCentralAnnular/NonCentralTotal）。
  - 出典: Meeus Ch.54 種別フローチャート / NASA 種別定義。**境界しきい値と l2 符号規約（ISSUE-021）を実装コメントに明記**。
- **全球接触 P1/U1/U4/P4**: 半影縁 l1（P1/P4: 半影が地球縁に外接）・本影縁 l2（U1/U4: 本影が地球縁に外接）が地球面に最初/最後に触れる時刻。`(x²+y²)` と地球縁の交差を Brent で求解（conventions §8 外接/内接の全球版）。出典: Explanatory Supplement Ch.11 / Meeus Ch.54。
- **最大食地点**: 最大食時刻の影軸が地表を貫く点（地理緯度経度）。地心→測地変換（ISSUE-010/011, WGS84）。
- **食帯幅・中心食継続**: 最大食点での本影/反本影の地表での幅、影が通過する継続時間。出典: Meeus Ch.54 / NASA（path width, duration of totality/annularity）。
- **magnitude/obscuration**: 最大食点での食分・面積比。境界条件明示（accuracy.md §2.2, acos クランプ）。

## 単位 / 時刻系 / 座標系
- 時刻系: ベッセル要素・最大食求解は **TT 基準**（accuracy.md §0(a) 幾何精度）。接触/最大食は **UTC も併記**（ΔT/UT1 経由, accuracy.md §0(b)/§2.3 将来は予測律速）。
- 座標系: 基本面（x,y,l1,l2, Re）＋地理座標（GeoPoint, 測地, conventions §3）。地心→測地は ISSUE-010/011。
- 単位: gamma Re、帯幅 km（conventions §4）、継続 秒、高度 度。
- 中心食継続・帯幅の地上計算は WGS84（earth_model, conventions §4）。

## アルゴリズム概要
1. 最大食時刻を求解: `g(t) = x(t)²+y(t)²` を ISSUE-009（1次元最小化, Brent/黄金分割）で最小化（直接供給 ISSUE-037 推奨, 精度）。→ gamma, greatest time。
2. 種別判定: gamma と最大食時の l1, l2（ISSUE-021 符号）で Partial/中心食を分岐。中心食は l2 符号で Annular/Total、経過中の l2 符号反転で Hybrid。中心線が地球縁すれすれなら NonCentral。
3. 全球接触: P1/P4（半影 l1 が地球縁外接）、U1/U4（本影 l2 が地球縁外接）を `(x²+y²)` と地球縁の交差で Brent 求解（conventions §8）。中心食でなければ U1/U4 は None。
4. 最大食地点: 影軸の地表貫通点を測地座標へ（ISSUE-010/011）。太陽高度算出。
5. 帯幅・中心食継続: 最大食点での本影/反本影地表幅・通過継続を算出（中心食のみ Some, それ以外 None）。
6. magnitude/obscuration を最大食点で算出（acos クランプ, accuracy.md §2.2）。
7. 各時刻を TT と UTC で返す（accuracy.md §0）。
- 数値安定性: gamma 最小化は Newton 単独禁止（ISSUE-009, conventions §11）。acos/asin クランプ（accuracy.md §2.2）。l2 符号境界（皆既↔金環）で種別が連続に切替わること。地球縁交差が無い（部分食で U1/U4 None）を正しく扱う。

## 受け入れテスト
accuracy.md テストレベル **L5（全球日食）**。**NASA 公開ベッセル値・全球状況との比較**（品質基準）。
- **NASA 全球整合（第二義, data-sources §4.1）**: NASA 5千年カタログの既知日食（皆既/金環/ハイブリッド/部分/非中心）で、種別・gamma・最大食時刻・食分・帯幅を比較。**k 慣習を Espenak（conventions §9: 本影 EspenakUmbral=0.272281 / 半影 EspenakPenumbral）に揃え、ΔT も合わせる**。系統差を accuracy.md に記録（絶対基準にしない, accuracy.md §3.1）。基準は fixtures（conventions §11）。
- **DE 差分（第一義, accuracy.md §3.1）**: 解析暦 vs DE440 で同一全球パイプラインを通し gamma・最大食時刻の差を層分解（§4）。
- **MockEphemeris 人工ケース（accuracy.md §3.1）**: 完全中心皆既（gamma≈0, l2<0, Total）/ 明確な金環（l2>0, Annular）/ ハイブリッド（l2 符号反転, Hybrid）/ 部分（gamma 大, Partial, U1/U4=None）/ 非中心（gamma≈1.0 すれすれ, NonCentral）。各で種別・接触の Option を検証。
- **種別境界テスト（必須, 品質基準）**: l2≈0（皆既↔金環境界）と gamma≈中心食限界（中心↔非中心境界）を掃引し、種別遷移が l2 符号・gamma しきい値で正しく切替わること。**Hybrid の境界条件（経過中の l2 符号反転）**を明示検証。
- **k 値系統差テスト（必須, 品質基準）**: `EspenakUmbral` vs `IauMean`（conventions §9）で本影半径が変わり、皆既/金環/ハイブリッド境界の判定がずれることを定量化し、系統差を accuracy.md へ記録（誤差を隠さない, accuracy.md §0）。
- ゴールデン20（accuracy.md §3.4）の全球部分（種別・gamma・最大食・食分）。
- UTC/TT 両方が返ること、将来日食で delta_t_uncertainty が metadata に乗ること（accuracy.md §0）。

## 許容誤差
- accuracy.md §2.1 幾何バジェット（§4 層分解）:
  - 最大食時刻 ±1.5s（gamma 最小化 solver 収束 0.05″, root_tolerance 目標の 1/10）。
  - gamma: x,y 精度律速（影幾何＋暦, ≲0.49″ 合成 ≈ 1.0s 相当）。
  - 食分 ±0.0005（0.001食分≈1.9″, accuracy.md §2.2）。
  - 中心線位置 sub-km（≲0.5km, 最大食地点に適用, accuracy.md §2.1）。
  - 帯幅・継続は l2/l1 と地上射影律速。
- **k 値選択による系統差**（conventions §9）は皆既/金環/ハイブリッド境界に出る。許容で吸収せず、慣習を揃えて比較し系統差を記録（accuracy.md §0/§3.1）。
- UTC 絶対時刻は将来 ΔT/UT1 律速（accuracy.md §2.3）。幾何（TT）精度と分離して報告。許容を通すための拡大禁止（conventions §11）。

## 実装メモ
- **最大食の供給源は ISSUE-037（直接, fit誤差ゼロ）を既定**（精度）。経路/帯幅の多点は ISSUE-022（多項式）でも可（L7 残差で確認）。
- **種別境界条件を明示（品質基準）**: Total/Annular = 最大食時 l2 符号、Hybrid = 経過中の l2 符号反転、NonCentral = 影軸は当たるが中心食条件を地球縁で一部のみ満たす（gamma が中心食限界近傍）。各しきい値の出典（Meeus Ch.54 / NASA）と l2 符号規約（ISSUE-021）をコメント化。
- **k 値選択の系統差を明記（品質基準）**: 既定 IauMean（単一 k）と Espenak 2値（本影/半影で別 k）で皆既/金環判定が変わりうる。NASA 照合は Espenak へ切替、系統差を accuracy.md に記録。metadata に lunar_radius_model を必ず載せる（architecture §9）。
- U1/U4（中心食接触）は部分食では None。Option を正しく埋める（api-draft §3.4, §6 未確定: Option 設計）。
- NonCentral 系を v1.0 で公開するかは api-draft §6 未確定（v0.1 は Partial 中心）。型は用意し、分類ロジックは実装。
- gamma 最小化・接触求解は Newton 単独禁止（ISSUE-008/009, conventions §11）。acos クランプ必須。
- レビュー重点: 種別境界（l2 符号・gamma しきい値・Hybrid）、k 系統差の記録、UTC/TT 両返し、U1/U4 の Option、最大食供給源（直接既定）。
