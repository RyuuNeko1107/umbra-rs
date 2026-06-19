//! 計算メタデータ型（結果型に付随, `docs/accuracy.md` §0, ISSUE-043 S4a）。
//!
//! [`CalculationMetadata`] は計算結果に付随し、暦・ΔT・地球/月モデル・精度プロファイル・
//! ライブラリ版などの「レシピ」と、生成時刻印・ΔT 不確かさを保持する。`fingerprint()` は
//! DB 差分再生成用に **レシピ部分のみ** から決定的なフィンガープリント（16 進小文字）を返し、
//! 時刻印（`generated_at`）と瞬時出力（`delta_t_uncertainty_seconds`）は除外する。
//!
//! 本体型はメイン実装が追加する（ISSUE-043 S4a）。

use umbra_core::UtcInstant;

use crate::config::AccuracyProfile;

/// 計算結果に付随するメタデータ（accuracy.md §0）。暦・ΔT・地球/月モデル・精度プロファイル・
/// ライブラリ版の「レシピ」と、生成時刻印・ΔT 不確かさ（将来 UTC 律速）を保持する。
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize))]
pub struct CalculationMetadata {
    /// ライブラリ版（再現性キー）。
    pub library_version: String,
    /// 暦モデル名（例 `"VSOP87D + ELP2000-82B"`）。
    pub ephemeris_model: String,
    /// 暦版（採用打切り次数・達成残差を含む識別子）。
    pub ephemeris_version: String,
    /// ΔT モデル名。
    pub delta_t_model: String,
    /// ΔT 不確実性帯（秒, accuracy.md §0。将来 UTC 律速。瞬時ごとの出力）。
    pub delta_t_uncertainty_seconds: f64,
    /// 地球モデル名。
    pub earth_model: String,
    /// 月半径モデル名。
    pub lunar_radius_model: String,
    /// 精度プロファイル。
    pub accuracy_profile: AccuracyProfile,
    /// 生成時刻（UTC）。
    pub generated_at: UtcInstant,
}

impl CalculationMetadata {
    /// DB 差分再生成用の決定的フィンガープリント（16 進小文字）。
    ///
    /// 計算の**レシピ**（`library_version` / `ephemeris_model` / `ephemeris_version` /
    /// `delta_t_model` / `earth_model` / `lunar_radius_model` / `accuracy_profile`）のみを識別する。
    /// **`generated_at`（時刻印）と `delta_t_uncertainty_seconds`（瞬時ごとの出力）は除外**するため、
    /// 同一レシピなら計算時刻・対象瞬時に依らず同一値を返す（plan §22）。
    ///
    /// FNV-1a 64-bit（決定的・依存なし。暗号用途ではなく差分検出用）。フィールド境界は
    /// ユニットセパレータ `U+001F` で区切り、`"ab"+"c"` と `"a"+"bc"` の衝突を避ける。
    pub fn fingerprint(&self) -> String {
        let profile = match self.accuracy_profile {
            AccuracyProfile::Standard => "Standard",
            AccuracyProfile::Reference => "Reference",
        };
        let recipe = [
            self.library_version.as_str(),
            self.ephemeris_model.as_str(),
            self.ephemeris_version.as_str(),
            self.delta_t_model.as_str(),
            self.earth_model.as_str(),
            self.lunar_radius_model.as_str(),
            profile,
        ]
        .join("\u{1f}");
        fnv1a_hex(recipe.as_bytes())
    }
}

/// FNV-1a 64-bit を 16 進小文字（16 桁ゼロ詰め）で返す。決定的・依存なし。
///
/// 定数出典: Fowler–Noll–Vo hash 64-bit 変種（IETF `draft-eastlake-fnv` §2.1 / FNV 公式）。
/// offset_basis = 14695981039346656037 = `0xcbf29ce484222325`、
/// FNV_prime = 1099511628211 = `0x100000001b3`。暗号用途ではなく DB 差分検出用。
fn fnv1a_hex(bytes: &[u8]) -> String {
    // FNV-1a 64-bit offset basis（IETF draft-eastlake-fnv §2.1）。
    let mut hash: u64 = 0xcbf2_9ce4_8422_2325;
    for &b in bytes {
        hash ^= u64::from(b);
        // FNV prime 64-bit（IETF draft-eastlake-fnv §2.1）。
        hash = hash.wrapping_mul(0x0000_0100_0000_01b3);
    }
    format!("{hash:016x}")
}

#[cfg(test)]
mod tests {
    //! ISSUE-043 S4a 受け入れテスト（strict・計算メタデータ / fingerprint）。
    //!
    //! ## オラクル戦略
    //! `fingerprint()` は実装依存のハッシュであり**厳密値は pin しない**。代わりに公開仕様が
    //! 約束する **性質**で縛る:
    //!   1. 決定性（同一フィールド → 同一文字列）
    //!   2. 形式（非空・16 進小文字）
    //!   3. 除外（generated_at / delta_t_uncertainty_seconds だけが違う → 同一）
    //!   4. レシピ感度（library_version / ephemeris_model / ephemeris_version / delta_t_model /
    //!      earth_model / lunar_radius_model / accuracy_profile のいずれかが違う → 異なる）
    //!
    //! 本体型未実装のため red（未存在シンボルでコンパイルエラー）想定。

    use crate::calc_metadata::CalculationMetadata;
    use crate::config::AccuracyProfile;
    use umbra_core::UtcInstant;

    /// 代表的な完全メタデータ（レシピ全フィールド非空）。generated_at は固定 UTC。
    fn base() -> CalculationMetadata {
        CalculationMetadata {
            library_version: "0.1.0".to_string(),
            ephemeris_model: "ELP/MPP02+VSOP87D".to_string(),
            ephemeris_version: "2024a".to_string(),
            delta_t_model: "EspenakMeeus".to_string(),
            delta_t_uncertainty_seconds: 0.5,
            earth_model: "WGS84".to_string(),
            lunar_radius_model: "IauMean".to_string(),
            accuracy_profile: AccuracyProfile::Standard,
            generated_at: UtcInstant::from_gregorian(2026, 6, 18, 0, 0, 0.0).unwrap(),
        }
    }

    /// 16 進小文字（0-9a-f）のみで構成され、かつ非空かを判定するヘルパ。
    fn is_nonempty_lower_hex(s: &str) -> bool {
        !s.is_empty()
            && s.chars()
                .all(|c| c.is_ascii_digit() || ('a'..='f').contains(&c))
    }

    // ============================================================
    // 形式・決定性
    // ============================================================

    /// `fingerprint()` は非空・16 進小文字（大文字・非 hex 文字を含まない）。
    /// 形式契約（16 進小文字）を縛り、大文字化や空文字列を返す変異を撃破する。
    #[test]
    fn fingerprint_is_nonempty_lowercase_hex() {
        let fp = base().fingerprint();
        assert!(
            is_nonempty_lower_hex(&fp),
            "fingerprint は非空・16進小文字であること: {fp:?}"
        );
    }

    /// 決定性: 同一レシピフィールドの 2 インスタンスは同一 fingerprint。
    /// 同一 base を 2 回構築して一致を縛り、非決定（乱数・時刻混入・アドレス依存）変異を撃破する。
    #[test]
    fn fingerprint_is_deterministic_for_identical_recipe() {
        assert_eq!(
            base().fingerprint(),
            base().fingerprint(),
            "同一フィールドなら fingerprint は常に同一"
        );
    }

    /// 既知レシピの fingerprint を 1 点 pin する（FNV-1a 64-bit は安定アルゴリズムなので値固定可）。
    /// 期待値は base() のレシピ文字列 "0.1.0\u{1f}ELP/MPP02+VSOP87D\u{1f}2024a\u{1f}EspenakMeeus\
    /// \u{1f}WGS84\u{1f}IauMean\u{1f}Standard" の UTF-8 を FNV-1a 64-bit した独立計算値。
    /// 演算子取り違え（`^=`→`|=` 等、決定性/感度テストでは生存する FNV 内部演算の変異）を撃破する。
    #[test]
    fn fingerprint_pins_known_recipe_value() {
        assert_eq!(
            base().fingerprint(),
            "bc28de7fd2b5570c",
            "FNV-1a 64-bit of the base() recipe (独立計算オラクル)"
        );
    }

    // ============================================================
    // 除外フィールド（generated_at / delta_t_uncertainty_seconds）
    // ============================================================

    /// generated_at だけが異なる 2 インスタンスは **同一** fingerprint（時刻印はレシピ外）。
    /// fingerprint に generated_at を誤って含める変異を撃破する。
    #[test]
    fn fingerprint_ignores_generated_at() {
        let a = base();
        let mut b = base();
        b.generated_at = UtcInstant::from_gregorian(1999, 1, 1, 12, 30, 15.0).unwrap();
        // 前提: generated_at は実際に異なる（除外検証の有効性確認）。
        assert_ne!(
            a.generated_at, b.generated_at,
            "テスト前提: generated_at は異なる値"
        );
        assert_eq!(
            a.fingerprint(),
            b.fingerprint(),
            "generated_at の違いは fingerprint に影響しない"
        );
    }

    /// delta_t_uncertainty_seconds だけが異なる 2 インスタンスは **同一** fingerprint
    /// （瞬時ごとの出力はレシピ外）。fingerprint に不確かさを誤って含める変異を撃破する。
    #[test]
    fn fingerprint_ignores_delta_t_uncertainty() {
        let a = base();
        let mut b = base();
        b.delta_t_uncertainty_seconds = a.delta_t_uncertainty_seconds + 1.0;
        assert_ne!(
            a.delta_t_uncertainty_seconds, b.delta_t_uncertainty_seconds,
            "テスト前提: 不確かさは異なる値"
        );
        assert_eq!(
            a.fingerprint(),
            b.fingerprint(),
            "delta_t_uncertainty_seconds の違いは fingerprint に影響しない"
        );
    }

    /// generated_at と delta_t_uncertainty の **両方**を同時に変えても fingerprint 不変
    /// （除外フィールドの組合せ）。片方しか除外しない変異を撃破する。
    #[test]
    fn fingerprint_ignores_both_excluded_fields_together() {
        let a = base();
        let mut b = base();
        b.generated_at = UtcInstant::from_gregorian(1980, 12, 31, 23, 59, 59.0).unwrap();
        b.delta_t_uncertainty_seconds = a.delta_t_uncertainty_seconds + 42.0;
        assert_eq!(
            a.fingerprint(),
            b.fingerprint(),
            "除外フィールドを同時に変えても fingerprint は不変"
        );
    }

    // ============================================================
    // レシピ感度（各フィールドを 1 つずつ変えて fingerprint が変わる）
    // ============================================================

    /// レシピ依存: library_version / ephemeris_model / ephemeris_version / delta_t_model /
    /// earth_model / lunar_radius_model / accuracy_profile を **1 つずつ** 変えると
    /// fingerprint が変わる。各レシピフィールドを fingerprint から取りこぼす変異を、
    /// フィールド名付きで個別に撃破する（どのフィールドが効いていないか特定可能）。
    #[test]
    fn fingerprint_is_sensitive_to_each_recipe_field() {
        let baseline = base().fingerprint();

        // (フィールド名, 当該フィールドのみを変える変異関数) の列。accuracy_profile 以外は
        // String フィールドなので別文字列へ、accuracy_profile は別バリアントへ変える。
        type Mutator = fn(&mut CalculationMetadata);
        let mutators: &[(&str, Mutator)] = &[
            ("library_version", |m| {
                m.library_version = "9.9.9".to_string()
            }),
            ("ephemeris_model", |m| {
                m.ephemeris_model = "DE440".to_string()
            }),
            ("ephemeris_version", |m| {
                m.ephemeris_version = "9999z".to_string()
            }),
            ("delta_t_model", |m| {
                m.delta_t_model = "MorrisonStephenson".to_string()
            }),
            ("earth_model", |m| m.earth_model = "GRS80".to_string()),
            ("lunar_radius_model", |m| {
                m.lunar_radius_model = "EspenakUmbral".to_string()
            }),
            ("accuracy_profile", |m| {
                m.accuracy_profile = AccuracyProfile::Reference
            }),
        ];

        for (field, mutate) in mutators {
            let mut m = base();
            mutate(&mut m);
            assert_ne!(
                m.fingerprint(),
                baseline,
                "レシピフィールド `{field}` の変化は fingerprint を変えるべき"
            );
            // 変えた値でも形式（16進小文字・非空）は保つ。
            assert!(
                is_nonempty_lower_hex(&m.fingerprint()),
                "`{field}` 変更後も fingerprint は 16進小文字・非空であること"
            );
        }
    }

    /// 異なるレシピ同士は互いに異なる fingerprint（衝突しない）の補強: accuracy_profile を
    /// Standard と Reference にした 2 インスタンスは異なる。除外フィールド検証とは独立に
    /// 「レシピが違えば違う」方向を念押しする。
    #[test]
    fn fingerprint_distinguishes_accuracy_profiles() {
        let mut std_m = base();
        std_m.accuracy_profile = AccuracyProfile::Standard;
        let mut ref_m = base();
        ref_m.accuracy_profile = AccuracyProfile::Reference;
        assert_ne!(
            std_m.fingerprint(),
            ref_m.fingerprint(),
            "Standard と Reference の accuracy_profile は fingerprint を区別する"
        );
    }
}
