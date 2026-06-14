# ISSUE-041: 物理定数・規約定数の集約（出典付き定数モジュール）

- crate: umbra-core
- 依存: ISSUE-001（規約・単位）, ISSUE-002（`Radians`/`Degrees` newtype）, ISSUE-010（WGS84 楕円体 a・1/f）, conventions §4/§9/§11
- milestone: M1（基盤層。多数の Issue が参照するため早期に固定。umbra-core 基礎）
- モード(tdd-workflow): **strict**（定数の値・出典・単位は全計算の基準で SemVer 公開境界。magic number 禁止（conventions §11）の全体担保なので strict）

## 目的
`umbra-rs` 全体で使う物理定数・規約定数を、**出典付きの単一定数モジュール**へ集約する。conventions §11「出典不明の数式 / magic number 禁止」を**全クレートで担保**する基盤（architecture §11, conventions §10/§11）。
- 主要定数（出典付き）:
  - 光速 `c = 299792.458 km/s`（CODATA/SI 定義値）。
  - 天文単位 `AU = 149597870.7 km`（IAU2012 定義値）。
  - `TT − TAI = 32.184 s`（IAU 定義、ISSUE-006）。
  - 地球赤道半径 `Re = WGS84 a = 6378137.0 m`（ベッセル無次元化基準、conventions §4。ISSUE-010 と単一ソース）。
  - 月半径係数 k（`IauMean=0.2725076` / `EspenakUmbral=0.272281` / `EspenakPenumbral=0.2725076`、conventions §9）。
  - 太陽半径 `R_sun = 696000 km`（IAU2015 公称、conventions §9）。
  - WGS84 扁平率 `1/f = 298.257223563`（conventions §4）。
  - 角度変換・J2000 基準エポック・ユリウス世紀/千年の日数 等の規約定数。
- 各定数に**出典（文献・規約・式番号）・単位・採用根拠**を doc コメントで付す（conventions §10）。

## 非目的
- 半径**モデル**の選択ロジック（`LunarRadiusModel`/`SolarRadiusModel` enum と既定選択。conventions §9・umbra-eclipse の設定。本 Issue は値の保持元）。
- 章動/歳差の級数係数（ISSUE-040/033/034。大量係数はデータパイプライン管理、本 Issue は**スカラ規約定数**のみ）。
- ΔT/EOP の時変値（ISSUE-007。本 Issue は `TT−TAI` 等の固定オフセットのみ）。
- WGS84 楕円体の幾何計算（ISSUE-010。本 Issue は a・1/f の**値の単一ソース**を提供、ISSUE-010 と重複定義しない）。

## 公開インターフェース（※署名はレビュー確定）
api-draft に明示の章はないが conventions §4/§9/§11 を型化。newtype 付きで提供（生 f64 を単位付き量として配らない、conventions §1）:
```rust
pub mod constants {
    use crate::{Kilometers, Meters, Radians};

    // --- 物理定数（出典コメント必須）---
    pub const SPEED_OF_LIGHT_KM_S: f64 = 299_792.458;          // SI 定義値（CODATA）
    pub const ASTRONOMICAL_UNIT_KM: f64 = 149_597_870.7;       // IAU2012 定義値
    pub const TT_MINUS_TAI_SECONDS: f64 = 32.184;              // IAU 定義（ISSUE-006）

    // --- 地球モデル（conventions §4・ISSUE-010 と単一ソース）---
    pub const WGS84_SEMI_MAJOR_AXIS_M: f64 = 6_378_137.0;      // a
    pub const WGS84_INVERSE_FLATTENING: f64 = 298.257_223_563; // 1/f
    pub fn earth_equatorial_radius() -> Kilometers;            // Re（ベッセル無次元化基準）

    // --- 半径モデル係数（conventions §9）---
    pub const LUNAR_K_IAU_MEAN: f64 = 0.272_507_6;
    pub const LUNAR_K_ESPENAK_UMBRAL: f64 = 0.272_281;
    pub const LUNAR_K_ESPENAK_PENUMBRAL: f64 = 0.272_507_6;
    pub const SOLAR_RADIUS_KM: f64 = 696_000.0;                // IAU2015 公称

    // --- 時刻・エポック規約 ---
    pub const J2000_JD_TT: f64 = 2_451_545.0;                  // J2000.0
    pub const DAYS_PER_JULIAN_CENTURY: f64 = 36_525.0;
    pub const DAYS_PER_JULIAN_MILLENNIUM: f64 = 365_250.0;
}
```
- 値は `f64` 定数だが、外部公開境界では可能な限り newtype（`Kilometers` 等）でラップして渡す（conventions §1）。各定数の doc に出典・単位・採用根拠。

## 数式・アルゴリズムの出典（定数の一次出典）
- **光速 c**: SI 定義値 299792458 m/s = 299792.458 km/s（CODATA / SI）。light-time で使用（ISSUE-015）。
- **AU**: IAU2012 Resolution B2 = 149597870700 m = 149597870.7 km。暦境界の AU→km 変換（conventions §4, ISSUE-013）。
- **TT − TAI = 32.184 s**: IAU 定義（ISSUE-006, conventions §6）。
- **Re = WGS84 a = 6378137.0 m, 1/f = 298.257223563**: WGS84（conventions §4）。NASA 慣習（赤道半径 6378.137 km 基準）との差異は accuracy.md/conventions §4 に記録（ISSUE-021 §63）。
- **月半径係数 k**: conventions §9（IauMean 0.2725076 / Espenak 2 値）。出典 = Espenak/NASA TP-2006-214141・IAU 平均半径（data-sources §4.1）。
- **太陽半径 696000 km**: IAU2015 Resolution B3 公称太陽半径（conventions §9）。
- **J2000・ユリウス世紀/千年**: IAU 慣習（VSOP87/ELP の T 単位、ISSUE-013/014/035）。
- 要確認: 各定数の最終桁・有効数字を一次規約文書で再確認（特に k の出典桁、NASA Re 慣習差の記録方法）。

## 単位 / 時刻系 / 座標系
- 単位は定数名に明示（`_KM_S` / `_KM` / `_M` / `_SECONDS`）。conventions §11「単位はフィールド名/定数名で明示」。
- 角度規約定数を置く場合はラジアン基準（conventions §1/§2）。
- 時刻系: `TT−TAI`・J2000 は TT 基準（conventions §6）。本 Issue は固定オフセット/エポックのみ（時変は ISSUE-007）。
- 座標系: 定数自体はフレーム非依存。Re はベッセル無次元化（FundamentalPlane, conventions §5）の基準。

## アルゴリズム概要
1. 定数を `umbra-core::constants` に集約し、各々に出典・単位・採用根拠の doc を付す。
2. WGS84 a・1/f は ISSUE-010 と**単一ソース**（重複定義禁止。ISSUE-010 が本モジュールを参照する構成 or 本モジュールが正本）。レビューで正本側を確定。
3. 半径係数 k・太陽半径は umbra-eclipse の `LunarRadiusModel`/`SolarRadiusModel`（conventions §9）が値ソースとして参照。
4. AU/c は暦・light-time（ISSUE-013/015）が参照。
- 数値安定性: 定数は SI/IAU 定義値の最大有効桁で記述。`as` 暗黙変換に頼らず単位明示。禁止: 同一定数の二重定義、出典なし定数（conventions §11）。

## 受け入れテスト
accuracy.md テストレベル **L1（純数学・定数）**。基準は一次規約文書（実装からのコピーでなく規約値）。
- 各定数が一次規約値と一致（c, AU, TT−TAI, Re, 1/f, k 3 値, R_sun, J2000, 日数）。
- **単一ソース検査**: WGS84 a が ISSUE-010 と本モジュールで同一値（二重定義による不一致がコンパイル/テストで検出される構成）。
- 単位整合: `earth_equatorial_radius()` が `Kilometers(6378.137)`（m→km 変換）を返す。
- doc 出典の存在: 各定数 doc に出典文字列が含まれること（doc テスト or lint。conventions §10）。
- magic number 回帰ガード: 下流クレートで同値の生リテラル直書きが無いことをレビュー/grep ガードで担保（conventions §11、CI 補助）。

## 許容誤差
- 定数は**定義値・公称値**のため誤差なし（厳密一致）。丸めは f64 表現限界のみ。
- 半径モデル k の選択差（IauMean vs Espenak）は系統差として accuracy.md §2.2 に記録（誤差化しない、conventions §9）。
- NASA Re 慣習（6378.137 km）との差は accuracy.md/conventions §4 に記録（誤差に誤帰属しない、ISSUE-021）。

## 実装メモ
- **magic number 禁止の全体担保**（conventions §11）が本 Issue の主眼。全クレートの数式実装は本モジュールを参照し、生リテラルを散らさない。
- WGS84 a の正本側（本モジュール or ISSUE-010）をレビューで一意に決め、もう一方は参照のみ（二重定義禁止）。
- 半径係数 k・太陽半径は値のみ保持。モデル選択 enum（conventions §9）は umbra-eclipse 側（責務分離）。
- 定数の出典・単位・採用根拠を doc に必ず残す（conventions §10）。NASA/Espenak 慣習差は記録（隠さない、accuracy.md §0）。
- レビュー重点: 出典の一次性・有効桁、単一ソース（Re）、単位明示、newtype ラップ、下流 magic number 排除の担保手段。
