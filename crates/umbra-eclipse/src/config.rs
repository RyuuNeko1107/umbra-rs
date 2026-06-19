//! 日食エンジン設定型（`docs/api-draft.md` §3.1・`docs/accuracy.md` §1, ISSUE-043 S4a）。
//!
//! 公開 2 層の精度プロファイル（`AccuracyProfile`）、月・太陽半径モデル、地球モデル、
//! 大気差モデル、根探索許容・経路サンプル間隔をまとめた [`EngineConfig`] を提供する。
//! 既定は本番標準（`EngineConfig::standard()`）。
//!
//! 本体型はメイン実装が追加する（ISSUE-043 S4a）。

use umbra_core::constants::SOLAR_RADIUS_KM;
use umbra_core::EarthModel;

use crate::horizontal::RefractionModel;

/// 精度プロファイル（公開 2 層, accuracy.md §1）。探索段の高速化は非公開の内部粗スキャンで扱う。
#[non_exhaustive]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize))]
#[cfg_attr(feature = "serde", serde(tag = "type"))]
pub enum AccuracyProfile {
    /// 本番標準（VSOP87D+ELP/MPP02 フル・IAU2006/2000A・見かけ補正 ON）。
    Standard,
    /// 高精度参照（回帰・ベンチ・差分テストの第一義オラクル）。
    Reference,
}

/// 月半径モデル（k = R_moon / R_earth, conventions §8 / physical-models）。
#[non_exhaustive]
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum LunarRadiusModel {
    /// IAU 平均（k = 0.2725076, 既定）。
    IauMean,
    /// Espenak 本影用（k = 0.272281）。
    EspenakUmbral,
    /// Espenak 半影用（k = 0.2725076）。
    EspenakPenumbral,
}

impl LunarRadiusModel {
    /// 月半径係数 k（= R_moon / R_earth）。
    ///
    /// 出典（conventions §8/§9・physical-models §C）: IauMean は IAU 採択 k = 0.2725076、
    /// EspenakUmbral は Espenak/NASA の本影用 k = 0.272281（縮小補正でやや小）、
    /// EspenakPenumbral は半影用で IauMean と同値 0.2725076。
    pub fn k(&self) -> f64 {
        match self {
            LunarRadiusModel::IauMean | LunarRadiusModel::EspenakPenumbral => 0.272_507_6,
            LunarRadiusModel::EspenakUmbral => 0.272_281,
        }
    }

    /// モデル名（`CalculationMetadata.lunar_radius_model` のレシピ識別子）。
    pub fn name(&self) -> &'static str {
        match self {
            LunarRadiusModel::IauMean => "IauMean",
            LunarRadiusModel::EspenakUmbral => "EspenakUmbral",
            LunarRadiusModel::EspenakPenumbral => "EspenakPenumbral",
        }
    }
}

/// 太陽半径モデル。
#[non_exhaustive]
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum SolarRadiusModel {
    /// IAU 2015 公称太陽半径（696000 km）。
    Iau2015,
}

impl SolarRadiusModel {
    /// 太陽物理半径 \[km\]。
    pub fn radius_km(&self) -> f64 {
        match self {
            SolarRadiusModel::Iau2015 => SOLAR_RADIUS_KM,
        }
    }
}

/// 日食エンジン設定。`standard()`（既定）／`reference()` のショートカットを持つ。
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EngineConfig {
    /// 精度プロファイル。
    pub accuracy: AccuracyProfile,
    /// 地球モデル（観測者座標・楕円体）。
    pub earth_model: EarthModel,
    /// 月半径モデル（k）。
    pub lunar_radius_model: LunarRadiusModel,
    /// 太陽半径モデル。
    pub solar_radius_model: SolarRadiusModel,
    /// 大気差モデル（局地高度・可視性）。
    pub refraction: RefractionModel,
    /// 根探索許容（秒）。目標の 1/10 以下（accuracy.md §2.1）。
    pub root_tolerance_seconds: f64,
    /// 経路サンプル間隔（秒, path 用）。
    pub path_sample_interval_seconds: f64,
}

impl EngineConfig {
    /// 本番標準設定。
    pub fn standard() -> Self {
        Self {
            accuracy: AccuracyProfile::Standard,
            earth_model: EarthModel::Wgs84,
            lunar_radius_model: LunarRadiusModel::IauMean,
            solar_radius_model: SolarRadiusModel::Iau2015,
            refraction: RefractionModel::Standard,
            root_tolerance_seconds: 0.1,
            path_sample_interval_seconds: 60.0,
        }
    }

    /// 高精度参照設定（根探索許容・サンプル間隔を標準より厳しく）。
    pub fn reference() -> Self {
        Self {
            accuracy: AccuracyProfile::Reference,
            earth_model: EarthModel::Wgs84,
            lunar_radius_model: LunarRadiusModel::IauMean,
            solar_radius_model: SolarRadiusModel::Iau2015,
            refraction: RefractionModel::Standard,
            root_tolerance_seconds: 0.01,
            path_sample_interval_seconds: 10.0,
        }
    }
}

impl Default for EngineConfig {
    /// 既定は `standard()`。
    fn default() -> Self {
        Self::standard()
    }
}

#[cfg(test)]
mod tests {
    //! ISSUE-043 S4a 受け入れテスト（strict・エンジン設定型）。
    //!
    //! ## オラクル戦略
    //! - `LunarRadiusModel::k()` / `SolarRadiusModel::radius_km()` は **api-draft §3.1 の
    //!   既知 k 値・SOLAR_RADIUS_KM（umbra_core::constants）** を独立オラクルとする。
    //! - `EngineConfig::standard()`/`reference()` は確定仕様の全フィールド値を exact 比較する
    //!   （実装方針には立ち入らず、公開 IF が約束する値のみを縛る）。
    //!
    //! 本体型未実装のため、これらは red（未存在シンボルでコンパイルエラー）想定。

    use crate::config::{AccuracyProfile, EngineConfig, LunarRadiusModel, SolarRadiusModel};
    use crate::horizontal::RefractionModel;
    use umbra_core::constants::SOLAR_RADIUS_KM;
    use umbra_core::EarthModel;

    // ============================================================
    // LunarRadiusModel::k()（月半径モデルの k = R_moon / R_earth）
    // ============================================================

    /// `LunarRadiusModel::k()` の 3 値 exact（api-draft §3.1）:
    /// IauMean=0.2725076, EspenakUmbral=0.272281, EspenakPenumbral=0.2725076。
    /// 各バリアントの k を個別に exact 比較し、k 値の取り違え（バリアント間の値入れ替え・
    /// 定数置換）変異を撃破する。
    #[test]
    fn lunar_radius_model_k_values_are_exact() {
        assert_eq!(LunarRadiusModel::IauMean.k(), 0.2725076);
        assert_eq!(LunarRadiusModel::EspenakUmbral.k(), 0.272281);
        assert_eq!(LunarRadiusModel::EspenakPenumbral.k(), 0.2725076);
    }

    /// IauMean と EspenakPenumbral は同値（0.2725076）、EspenakUmbral のみ別値（0.272281）。
    /// 「IauMean==EspenakPenumbral かつ ≠ EspenakUmbral」という関係を縛り、
    /// 3 値のうち 2 値が一致する構造の取り違え（例: Umbral と Penumbral の値入れ替え）を撃破する。
    #[test]
    fn lunar_radius_model_iaumean_equals_penumbral_and_differs_from_umbral() {
        assert_eq!(
            LunarRadiusModel::IauMean.k(),
            LunarRadiusModel::EspenakPenumbral.k(),
            "IauMean と EspenakPenumbral は同じ k (0.2725076)"
        );
        assert_ne!(
            LunarRadiusModel::EspenakUmbral.k(),
            LunarRadiusModel::IauMean.k(),
            "EspenakUmbral は umbral 補正で別値 (0.272281)"
        );
        // umbral は mean より小さい（影縮小補正）。値の大小関係も固定。
        assert!(
            LunarRadiusModel::EspenakUmbral.k() < LunarRadiusModel::IauMean.k(),
            "EspenakUmbral の k は IauMean より小さい"
        );
    }

    /// `LunarRadiusModel::name()` の 3 値 exact（レシピ識別子）:
    /// IauMean="IauMean", EspenakUmbral="EspenakUmbral", EspenakPenumbral="EspenakPenumbral"。
    /// 各バリアントの name を個別に exact 比較し、name アームの取り違え（文字列入れ替え・
    /// 定数置換）変異を撃破する。
    #[test]
    fn lunar_radius_model_name_values_are_exact() {
        assert_eq!(LunarRadiusModel::IauMean.name(), "IauMean");
        assert_eq!(LunarRadiusModel::EspenakUmbral.name(), "EspenakUmbral");
        assert_eq!(
            LunarRadiusModel::EspenakPenumbral.name(),
            "EspenakPenumbral"
        );
    }

    // ============================================================
    // SolarRadiusModel::radius_km()（太陽物理半径 [km]）
    // ============================================================

    /// `SolarRadiusModel::radius_km()`: Iau2015=696000.0、かつ
    /// umbra_core::constants::SOLAR_RADIUS_KM と一致（独立定数オラクル）。
    /// ハードコード値の取り違え・別定数参照を撃破する。
    #[test]
    fn solar_radius_model_iau2015_matches_constant() {
        assert_eq!(SolarRadiusModel::Iau2015.radius_km(), 696_000.0);
        assert_eq!(
            SolarRadiusModel::Iau2015.radius_km(),
            SOLAR_RADIUS_KM,
            "Iau2015 は umbra_core::constants::SOLAR_RADIUS_KM と一致"
        );
    }

    // ============================================================
    // EngineConfig::standard() / reference() / Default の全フィールド exact
    // ============================================================

    /// `EngineConfig::standard()` の全フィールド exact（確定仕様）:
    /// {Standard, Wgs84, IauMean, Iau2015, Standard(refraction),
    ///  root_tolerance_seconds: 0.1, path_sample_interval_seconds: 60.0}。
    /// 各既定値の取り違え（profile・モデル・許容・間隔の差し替え）を個別フィールドで撃破する。
    #[test]
    fn engine_config_standard_has_exact_fields() {
        let c = EngineConfig::standard();
        assert_eq!(c.accuracy, AccuracyProfile::Standard, "accuracy");
        assert_eq!(c.earth_model, EarthModel::Wgs84, "earth_model");
        assert_eq!(c.lunar_radius_model, LunarRadiusModel::IauMean, "lunar");
        assert_eq!(c.solar_radius_model, SolarRadiusModel::Iau2015, "solar");
        assert_eq!(c.refraction, RefractionModel::Standard, "refraction");
        assert_eq!(c.root_tolerance_seconds, 0.1, "root_tolerance_seconds");
        assert_eq!(
            c.path_sample_interval_seconds, 60.0,
            "path_sample_interval_seconds"
        );
    }

    /// `EngineConfig::reference()` の全フィールド exact（確定仕様）:
    /// {Reference, Wgs84, IauMean, Iau2015, Standard(refraction),
    ///  root_tolerance_seconds: 0.01, path_sample_interval_seconds: 10.0}。
    /// reference 固有の厳しい許容(0.01)・細かい間隔(10.0)・Reference profile を撃破対象に。
    #[test]
    fn engine_config_reference_has_exact_fields() {
        let c = EngineConfig::reference();
        assert_eq!(c.accuracy, AccuracyProfile::Reference, "accuracy");
        assert_eq!(c.earth_model, EarthModel::Wgs84, "earth_model");
        assert_eq!(c.lunar_radius_model, LunarRadiusModel::IauMean, "lunar");
        assert_eq!(c.solar_radius_model, SolarRadiusModel::Iau2015, "solar");
        assert_eq!(c.refraction, RefractionModel::Standard, "refraction");
        assert_eq!(c.root_tolerance_seconds, 0.01, "root_tolerance_seconds");
        assert_eq!(
            c.path_sample_interval_seconds, 10.0,
            "path_sample_interval_seconds"
        );
    }

    /// `Default::default() == EngineConfig::standard()`（既定は本番標準）。
    /// Default 実装が reference() 等を返す取り違えを撃破する。
    #[test]
    fn engine_config_default_equals_standard() {
        assert_eq!(EngineConfig::default(), EngineConfig::standard());
    }

    /// standard() ≠ reference()（accuracy/tolerance/interval が異なる別プロファイル）。
    /// 2 メソッドが同一値を返す退化（両者とも standard を返す等）を撃破する。
    /// 差は accuracy（Standard vs Reference）・許容（0.1 vs 0.01）・間隔（60 vs 10）に表れる。
    #[test]
    fn engine_config_standard_differs_from_reference() {
        let s = EngineConfig::standard();
        let r = EngineConfig::reference();
        assert_ne!(s, r, "standard と reference は別プロファイル");
        assert_ne!(s.accuracy, r.accuracy, "accuracy が異なる");
        assert_ne!(
            s.root_tolerance_seconds, r.root_tolerance_seconds,
            "root_tolerance_seconds が異なる"
        );
        assert_ne!(
            s.path_sample_interval_seconds, r.path_sample_interval_seconds,
            "path_sample_interval_seconds が異なる"
        );
    }

    /// reference の根探索許容は standard より厳しい（小さい）、サンプル間隔も細かい（小さい）。
    /// 「Reference = より高精度」の方向性（accuracy.md §1）を大小関係で固定し、
    /// 値の入れ替え（standard と reference の tolerance/interval 取り違え）を撃破する。
    #[test]
    fn reference_is_tighter_than_standard() {
        let s = EngineConfig::standard();
        let r = EngineConfig::reference();
        assert!(
            r.root_tolerance_seconds < s.root_tolerance_seconds,
            "reference の root_tolerance は standard より小さい"
        );
        assert!(
            r.path_sample_interval_seconds < s.path_sample_interval_seconds,
            "reference のサンプル間隔は standard より小さい"
        );
    }

    /// `EngineConfig` が `Copy`（フィールド全て Copy 由来）。設定を値渡しで複製できることを
    /// コンパイル時に縛り、`#[derive(Copy)]` 脱落を撃破する。
    #[test]
    fn engine_config_is_copy() {
        fn assert_copy<T: Copy>(_: T) {}
        let c = EngineConfig::standard();
        let d = c; // Copy
        assert_copy(d);
        assert_eq!(c, d); // c はムーブされず有効
    }

    // ============================================================
    // AccuracyProfile（精度プロファイル列挙）
    // ============================================================

    /// `AccuracyProfile` の `PartialEq`/`Copy`: Standard と Reference は異なり、
    /// 同一バリアントは等しい。Copy 境界へ通して導出脱落を撃破する。
    #[test]
    fn accuracy_profile_eq_and_copy() {
        fn assert_copy<T: Copy>(_: T) {}
        let p = AccuracyProfile::Standard;
        assert_copy(p);
        assert_eq!(AccuracyProfile::Standard, AccuracyProfile::Standard);
        assert_eq!(AccuracyProfile::Reference, AccuracyProfile::Reference);
        assert_ne!(
            AccuracyProfile::Standard,
            AccuracyProfile::Reference,
            "2 つの精度プロファイルは区別される"
        );
        // p はムーブされず有効（Copy）。
        assert_eq!(p, AccuracyProfile::Standard);
    }
}
