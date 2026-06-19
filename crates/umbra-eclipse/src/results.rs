//! 全球日食結果型（`search()` が返すデータコンテナ群, `docs/api-draft.md` §3.4, ISSUE-043 S4d）。
//!
//! 本モジュールは `search()` の返却要素 [`SolarEclipse`] と、その内訳である全球条件
//! [`GlobalCircumstances`]・最大食 [`GreatestEclipse`]・接触点 [`GlobalContact`] を提供する。
//! いずれも **pub フィールドのデータコンテナ**（コンストラクタを持たず、エンジンが struct
//! リテラルで構築し、各フィールドを読み出す）。
//!
//! - [`GlobalContact`] / [`GreatestEclipse`] / [`GlobalCircumstances`] は `Copy`・`PartialEq`
//!   （値比較）。
//! - [`SolarEclipse`] は [`BesselianPolynomial`] が `PartialEq` 非実装のため `Clone`・`Debug`
//!   のみ。
//!
//! 注（derive）: `PartialEq` と [`GlobalCircumstances`] の `Copy` は api-draft §3.4 の
//! `Clone`/`Debug` を **DB 差分比較・テスト・人間工学のため拡張**したもの。全フィールドが
//! `Copy`/`PartialEq` を満たすため無害（[`SolarEclipse`] のみ `BesselianPolynomial` 制約で拡張不可）。
//!
//! 全球の地表点は [`umbra_geo::GeoPoint`]（eclipse→geo 依存・循環なし）。

use umbra_core::{Degrees, Kilometers, TtInstant, UtcInstant};
use umbra_geo::GeoPoint;

use crate::bessel_poly::BesselianPolynomial;
use crate::calc_metadata::CalculationMetadata;
use crate::global::SolarEclipseKind;
use crate::horizontal::Visibility;
use crate::magnitude::{EclipseMagnitude, Obscuration};

/// 全球接触点（時刻 TT/UTC ＋ 地表点）。
#[derive(Clone, Copy, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize))]
pub struct GlobalContact {
    /// 接触の UTC 時刻。
    pub time_utc: UtcInstant,
    /// 接触の TT 時刻（幾何相対の一級値, conventions §6）。
    pub time_tt: TtInstant,
    /// 接触の地表点（中心線が地表に達する点）。
    pub position: GeoPoint,
}

/// 最大食の全球条件。
#[derive(Clone, Copy, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize))]
pub struct GreatestEclipse {
    /// 最大食の UTC 時刻。
    pub time_utc: UtcInstant,
    /// 最大食の TT 時刻。
    pub time_tt: TtInstant,
    /// 最大食地点（地表）。
    pub position: GeoPoint,
    /// 食分（皆既で 1 超可）。
    pub magnitude: EclipseMagnitude,
    /// 食面積（0..1）。
    pub obscuration: Obscuration,
    /// 帯幅 \[km\]（中心食のみ Some）。JSON は単位明示で `path_width_km`（A7）。
    #[cfg_attr(feature = "serde", serde(rename = "path_width_km"))]
    pub path_width: Option<Kilometers>,
    /// 中心食継続時間 \[s\]（中心食のみ Some）。JSON は単位明示で `central_duration_seconds`（A7）。
    #[cfg_attr(feature = "serde", serde(rename = "central_duration_seconds"))]
    pub central_duration: Option<f64>,
    /// 最大食地点での太陽高度。JSON は単位明示で `sun_altitude_deg`（A7）。
    #[cfg_attr(feature = "serde", serde(rename = "sun_altitude_deg"))]
    pub sun_altitude: Degrees,
}

/// 全球条件（種別・P1/U1/最大食/U4/P4・gamma）。
#[derive(Clone, Copy, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize))]
pub struct GlobalCircumstances {
    /// 日食種別。
    pub kind: SolarEclipseKind,
    /// 部分食開始 P1（外接, 常に存在しうる）。
    pub partial_begin: Option<GlobalContact>,
    /// 中心食開始 U1（中心食のみ）。
    pub central_begin: Option<GlobalContact>,
    /// 最大食。
    pub greatest: GreatestEclipse,
    /// 中心食終了 U4（中心食のみ）。
    pub central_end: Option<GlobalContact>,
    /// 部分食終了 P4（外接）。
    pub partial_end: Option<GlobalContact>,
    /// 影軸の地心最小距離 gamma（Re）。
    pub gamma: f64,
}

/// `search()` の各要素（全球日食）。
#[derive(Clone, Debug)]
#[cfg_attr(feature = "serde", derive(serde::Serialize))]
pub struct SolarEclipse {
    /// 安定キー（DB 用, A4: 最大食 UTC 日付 + lunation 番号）。
    pub event_key: String,
    /// 日食種別。
    pub kind: SolarEclipseKind,
    /// 全球条件。
    pub global: GlobalCircumstances,
    /// ベッセル多項式（経路・エクスポート用, ISSUE-022）。
    pub bessel: BesselianPolynomial,
    /// 計算メタデータ（accuracy.md §0）。
    pub metadata: CalculationMetadata,
}

/// 局地接触（時刻 TT/UTC ＋ 観測フィールド, api-draft §3.4）。
///
/// `local_contacts::ContactInstant`（時刻のみの幾何ソルバ出力）に、ISSUE-028 の太陽地平座標・
/// position angle・可視を付与した公開 result 型。EclipseEngine（S7）が組み立てる。
#[derive(Clone, Copy, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize))]
pub struct LocalContact {
    /// 接触の UTC 時刻。
    pub time_utc: UtcInstant,
    /// 接触の TT 時刻（幾何相対の一級値, conventions §6）。
    pub time_tt: TtInstant,
    /// 接触時の太陽高度（conventions §7）。JSON は単位明示で `sun_altitude_deg`（A7）。
    #[cfg_attr(feature = "serde", serde(rename = "sun_altitude_deg"))]
    pub sun_altitude: Degrees,
    /// 接触時の太陽方位（北 0°・東回り）。JSON は単位明示で `sun_azimuth_deg`（A7）。
    #[cfg_attr(feature = "serde", serde(rename = "sun_azimuth_deg"))]
    pub sun_azimuth: Degrees,
    /// 接触の位置角（太陽周縁上, 北 0°・東回り）。JSON は単位明示で `position_angle_deg`（A7）。
    #[cfg_attr(feature = "serde", serde(rename = "position_angle_deg"))]
    pub position_angle: Degrees,
    /// 接触時に太陽が地平上か（可視）。
    pub visible: bool,
}

/// 局地接触集合（A3: c1..c4 は部分食地点で `None`、`maximum` は常に存在＝非 Option）。
/// フィールドは時系列順（c1, c2, maximum, c3, c4）。
#[derive(Clone, Copy, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize))]
pub struct LocalContactSet {
    /// 第1接触 C1（部分食開始・外接）。
    pub c1: Option<LocalContact>,
    /// 第2接触 C2（皆既/金環開始・内接, 中心食地点のみ）。
    pub c2: Option<LocalContact>,
    /// 最大食（どの地点でも定義される＝非 Option）。
    pub maximum: LocalContact,
    /// 第3接触 C3（皆既/金環終了・内接, 中心食地点のみ）。
    pub c3: Option<LocalContact>,
    /// 第4接触 C4（部分食終了・外接）。
    pub c4: Option<LocalContact>,
}

/// 観測地点の局地条件（`local_circumstances()` の結果）。
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize))]
pub struct LocalCircumstances {
    /// 接触集合（C1〜C4 ＋ 最大食）。
    pub contacts: LocalContactSet,
    /// 最大食での食分（皆既で 1 超可）。
    pub magnitude: EclipseMagnitude,
    /// 最大食での食面積（0..1）。
    pub obscuration: Obscuration,
    /// 最大食での太陽高度。JSON は単位明示で `maximum_altitude_deg`（A7）。
    #[cfg_attr(feature = "serde", serde(rename = "maximum_altitude_deg"))]
    pub maximum_altitude: Degrees,
    /// 可視性（6 値）。
    pub visibility: Visibility,
    /// 計算メタデータ（accuracy.md §0）。
    pub metadata: CalculationMetadata,
}

/// 可視日食（`next_visible_eclipse()` の結果）。日食とその観測地点の局地条件。
#[derive(Clone, Debug)]
#[cfg_attr(feature = "serde", derive(serde::Serialize))]
pub struct VisibleSolarEclipse {
    /// 全球日食。
    pub eclipse: SolarEclipse,
    /// 観測地点の局地条件。
    pub local: LocalCircumstances,
}

#[cfg(test)]
mod tests {
    //! ISSUE-043 S4d 受け入れテスト（strict・全球 result 型）。
    //!
    //! ## オラクル戦略
    //! 本群はすべて **pub フィールドのデータコンテナ**であり、振る舞いは「struct リテラルで
    //! 構築した各フィールドを取り違えずに保持・読み出す」「derive（Copy/PartialEq/Clone）が
    //! 仕様どおり付く」ことに尽きる。よって:
    //!   1. 非対称な値で各フィールドを構築し、フィールドごとに read-back（取り違え変異を撃破）。
    //!   2. `Option` フィールドは Some/None 両方を構築し保持を確認（Option 取り違えを撃破）。
    //!   3. 異値で `ne`・同値で `eq`、`Copy`/`Clone` をコンパイル時境界で縛る（derive 脱落を撃破）。
    //!
    //! 期待値はすべてテスト側で組み立てた構築値そのもの（外部マジック値なし）。
    //! 本体型未実装のため red（未存在シンボルでコンパイルエラー）想定。テストは
    //! `crate::results::{...}` を参照する。

    #![allow(clippy::excessive_precision)]

    use crate::bessel_poly::{BesselFitError, BesselianPolynomial};
    use crate::calc_metadata::CalculationMetadata;
    use crate::config::AccuracyProfile;
    use crate::global::SolarEclipseKind;
    use crate::horizontal::Visibility;
    use crate::magnitude::{EclipseMagnitude, Obscuration};
    use crate::polynomial::Polynomial;
    use crate::results::{
        GlobalCircumstances, GlobalContact, GreatestEclipse, LocalCircumstances, LocalContact,
        LocalContactSet, SolarEclipse, VisibleSolarEclipse,
    };

    use umbra_core::{Degrees, JulianDate2, Kilometers, TimeInterval, TtInstant, UtcInstant};
    use umbra_geo::GeoPoint;

    // ------------------------------------------------------------------
    // 構築ヘルパ（非対称な既知値で各部品を作る）
    // ------------------------------------------------------------------

    /// UTC 瞬時を整数引数で組む（取り違え判別のため互いに異なる値を選ぶ）。
    fn utc(year: i32, month: u8, day: u8, hour: u8, minute: u8, second: f64) -> UtcInstant {
        UtcInstant::from_gregorian(year, month, day, hour, minute, second).expect("有効な UTC 日時")
    }

    /// TT 瞬時を 2 要素 JD で組む（UTC と区別できる別スケール値）。
    fn tt(jd1: f64, jd2: f64) -> TtInstant {
        TtInstant::from_jd2(JulianDate2::new(jd1, jd2))
    }

    /// 地表点（lat, lon）を度から組む（S4c。緯度・経度の取り違えを判別するため lat≠lon）。
    fn geo(lat: f64, lon: f64) -> GeoPoint {
        GeoPoint::from_degrees(lat, lon).expect("有効な地表点")
    }

    /// 最小の [`BesselianPolynomial`] を 1 つ作るヘルパ（実コードの pub フィールドで直接構築）。
    /// `Polynomial::new` は存在しないため pub フィールド `coefficients` で構築する（polynomial.rs）。
    /// 各成分は定数多項式、tan f1/f2 は定数、fit_interval/fit_error も pub フィールド構築。
    fn minimal_bessel() -> BesselianPolynomial {
        let c = |v: f64| Polynomial {
            coefficients: vec![v],
        };
        BesselianPolynomial {
            epoch_tt: tt(2_451_545.0, 0.0),
            x: c(0.20),
            y: c(-0.30),
            d: c(0.2070),
            mu: c(1.2),
            l1: c(0.5400),
            l2: c(-0.0090),
            tan_f1: 0.004_65,
            tan_f2: 0.004_63,
            fit_interval: TimeInterval {
                start: tt(2_451_544.9, 0.0),
                end: tt(2_451_545.1, 0.0),
            },
            fit_error: BesselFitError {
                max_x: 1.0e-7,
                max_y: 2.0e-7,
                max_l1: 3.0e-7,
                max_l2: 4.0e-7,
            },
        }
    }

    /// 代表的なメタデータ（S4a。レシピ全フィールド非空, generated_at 固定）。
    fn metadata() -> CalculationMetadata {
        CalculationMetadata {
            library_version: "0.1.0".to_string(),
            ephemeris_model: "ELP/MPP02+VSOP87D".to_string(),
            ephemeris_version: "2024a".to_string(),
            delta_t_model: "EspenakMeeus".to_string(),
            delta_t_uncertainty_seconds: 0.5,
            earth_model: "WGS84".to_string(),
            lunar_radius_model: "IauMean".to_string(),
            accuracy_profile: AccuracyProfile::Standard,
            generated_at: utc(2026, 6, 18, 0, 0, 0.0),
        }
    }

    /// 中心食の最大食（path_width / central_duration が Some）を非対称値で組む。
    fn greatest_central() -> GreatestEclipse {
        GreatestEclipse {
            time_utc: utc(2024, 4, 8, 18, 17, 0.0),
            time_tt: tt(2_460_409.0, 0.123),
            position: geo(25.0, -104.0),
            magnitude: EclipseMagnitude(1.0566),
            obscuration: Obscuration(1.0),
            path_width: Some(Kilometers(197.0)),
            central_duration: Some(268.0),
            sun_altitude: Degrees(70.3),
        }
    }

    /// 部分食の最大食（path_width / central_duration が None）を非対称値で組む。
    fn greatest_partial() -> GreatestEclipse {
        GreatestEclipse {
            time_utc: utc(2025, 3, 29, 10, 47, 0.0),
            time_tt: tt(2_460_763.0, 0.456),
            position: geo(61.0, -77.0),
            magnitude: EclipseMagnitude(0.938),
            obscuration: Obscuration(0.91),
            path_width: None,
            central_duration: None,
            sun_altitude: Degrees(12.0),
        }
    }

    /// 4 接触点（P1/U1/U4/P4）を **互いに異なる time** で組む（取り違え判別用）。
    fn contact(h: u8, lat: f64) -> GlobalContact {
        GlobalContact {
            time_utc: utc(2024, 4, 8, h, 0, 0.0),
            time_tt: tt(2_460_409.0, f64::from(h) * 0.01),
            position: geo(lat, -100.0),
        }
    }

    // ==================================================================
    // GlobalContact: 3 フィールド保持・取り違えない
    // ==================================================================

    /// `GlobalContact` の 3 フィールド（time_utc / time_tt / position）が構築値どおりに保持され、
    /// 互いに取り違えられない。time_utc≠time_tt（別スケール）・position は非自明値で縛る。
    /// 殺す変異: フィールド read-back の入れ替え（time_utc↔time_tt 等）、position の取り違え。
    #[test]
    fn global_contact_holds_each_field() {
        let t_utc = utc(2017, 8, 21, 18, 25, 30.0);
        let t_tt = tt(2_457_987.0, 0.2685);
        let pos = geo(36.0, -80.0);
        let c = GlobalContact {
            time_utc: t_utc,
            time_tt: t_tt,
            position: pos,
        };
        assert_eq!(c.time_utc, t_utc, "time_utc を保持");
        assert_eq!(c.time_tt, t_tt, "time_tt を保持");
        assert_eq!(c.position, pos, "position を保持");
    }

    /// `GlobalContact` は `Copy`・`PartialEq`。異なる接触点は `ne`、同値は `eq`。
    /// 殺す変異: `#[derive(Copy)]`/`#[derive(PartialEq)]` の脱落。
    #[test]
    fn global_contact_is_copy_and_partial_eq() {
        fn assert_copy<T: Copy>(_: T) {}
        let a = contact(15, 30.0);
        let b = a; // Copy（move されない）
        assert_copy(b);
        assert_eq!(a, b, "Copy 後も等しい（a は有効）");
        // 同値で eq・異値で ne。
        assert_eq!(contact(15, 30.0), contact(15, 30.0), "同一構築値は eq");
        assert_ne!(contact(15, 30.0), contact(16, 30.0), "time 違いは ne");
        assert_ne!(contact(15, 30.0), contact(15, 31.0), "position 違いは ne");
    }

    // ==================================================================
    // GreatestEclipse: 8 フィールド保持 + Option 両方 + magnitude≠obscuration
    // ==================================================================

    /// `GreatestEclipse` の 8 フィールド（time_utc / time_tt / position / magnitude /
    /// obscuration / path_width / central_duration / sun_altitude）が構築値どおり保持される。
    /// magnitude(1.0566) と obscuration(1.0) を **異なる値**で与え、両者の取り違えを殺す。
    /// path_width=Some / central_duration=Some（中心食）を確認。
    #[test]
    fn greatest_eclipse_holds_all_fields_central() {
        let g = greatest_central();
        assert_eq!(g.time_utc, utc(2024, 4, 8, 18, 17, 0.0), "time_utc");
        assert_eq!(g.time_tt, tt(2_460_409.0, 0.123), "time_tt");
        assert_eq!(g.position, geo(25.0, -104.0), "position");
        assert_eq!(g.magnitude, EclipseMagnitude(1.0566), "magnitude");
        assert_eq!(g.obscuration, Obscuration(1.0), "obscuration");
        // magnitude と obscuration は異なる値（取り違え変異を撃破）。
        assert!(
            (g.magnitude.0 - g.obscuration.0).abs() > 1e-9,
            "magnitude と obscuration は別フィールド・別値であること"
        );
        assert_eq!(g.path_width, Some(Kilometers(197.0)), "path_width=Some");
        assert_eq!(g.central_duration, Some(268.0), "central_duration=Some");
        assert_eq!(g.sun_altitude, Degrees(70.3), "sun_altitude");
    }

    /// 部分食の `GreatestEclipse` は path_width=None・central_duration=None を保持する
    /// （中心食のみ Some の仕様）。他フィールドは保持。
    /// 殺す変異: None を Some に取り違える / 中心食用フィールドを常に Some にする。
    #[test]
    fn greatest_eclipse_holds_none_for_partial() {
        let g = greatest_partial();
        assert_eq!(g.path_width, None, "部分食 path_width は None");
        assert_eq!(g.central_duration, None, "部分食 central_duration は None");
        // 他フィールドは構築値どおり。
        assert_eq!(g.magnitude, EclipseMagnitude(0.938), "magnitude");
        assert_eq!(g.obscuration, Obscuration(0.91), "obscuration");
        assert_eq!(g.sun_altitude, Degrees(12.0), "sun_altitude");
        assert_eq!(g.position, geo(61.0, -77.0), "position");
    }

    /// `GreatestEclipse` は `Copy`・`PartialEq`。異なる最大食（中心食 vs 部分食）は `ne`、
    /// 同値は `eq`。Option の Some/None 差も `ne` に反映される。
    /// 殺す変異: `#[derive(Copy)]`/`#[derive(PartialEq)]` の脱落。
    #[test]
    fn greatest_eclipse_is_copy_and_partial_eq() {
        fn assert_copy<T: Copy>(_: T) {}
        let a = greatest_central();
        let b = a; // Copy
        assert_copy(b);
        assert_eq!(a, b, "Copy 後も等しい（a は有効）");
        assert_eq!(greatest_central(), greatest_central(), "同一構築値は eq");
        assert_ne!(
            greatest_central(),
            greatest_partial(),
            "中心食(Some) と部分食(None) は ne"
        );
    }

    // ==================================================================
    // GlobalCircumstances: kind/gamma/greatest + 5 接触の Some/None + 取り違え
    // ==================================================================

    /// 中心食の `GlobalCircumstances`: 全接触（P1/U1/最大食/U4/P4）が埋まり、
    /// kind=Total・gamma 保持・greatest 保持。**5 接触フィールドを互いに異なる time** で構築し、
    /// 各 Option フィールドから正しい接触点を read-back する（取り違えを殺す）。
    #[test]
    fn global_circumstances_central_holds_all_contacts() {
        let p1 = contact(16, 10.0);
        let u1 = contact(17, 20.0);
        let u4 = contact(19, 40.0);
        let p4 = contact(20, 50.0);
        let great = greatest_central();
        let gc = GlobalCircumstances {
            kind: SolarEclipseKind::Total,
            partial_begin: Some(p1),
            central_begin: Some(u1),
            greatest: great,
            central_end: Some(u4),
            partial_end: Some(p4),
            gamma: 0.3431,
        };
        assert_eq!(gc.kind, SolarEclipseKind::Total, "kind");
        assert_eq!(gc.gamma, 0.3431, "gamma");
        assert_eq!(gc.greatest, great, "greatest");
        // 5 接触の read-back（P1/U1/U4/P4 を互いに区別。取り違えを殺す）。
        assert_eq!(gc.partial_begin, Some(p1), "P1=partial_begin");
        assert_eq!(gc.central_begin, Some(u1), "U1=central_begin");
        assert_eq!(gc.central_end, Some(u4), "U4=central_end");
        assert_eq!(gc.partial_end, Some(p4), "P4=partial_end");
        // 取り違え検出: 異なる time の接触なので相互に等しくない。
        assert_ne!(gc.partial_begin, gc.central_begin, "P1 と U1 は別接触");
        assert_ne!(gc.central_end, gc.partial_end, "U4 と P4 は別接触");
    }

    /// 部分食の `GlobalCircumstances`: central 系（central_begin/central_end）が None、
    /// partial 系（partial_begin/partial_end）は Some。kind=Partial・gamma 保持。
    /// 殺す変異: 中心食用 Option（U1/U4）を常に Some にする / partial と central の取り違え。
    #[test]
    fn global_circumstances_partial_has_none_central() {
        let p1 = contact(10, 60.0);
        let p4 = contact(12, 62.0);
        let gc = GlobalCircumstances {
            kind: SolarEclipseKind::Partial,
            partial_begin: Some(p1),
            central_begin: None,
            greatest: greatest_partial(),
            central_end: None,
            partial_end: Some(p4),
            gamma: 1.21,
        };
        assert_eq!(gc.kind, SolarEclipseKind::Partial, "kind");
        assert_eq!(gc.gamma, 1.21, "gamma");
        assert_eq!(gc.partial_begin, Some(p1), "P1=Some");
        assert_eq!(gc.partial_end, Some(p4), "P4=Some");
        assert_eq!(gc.central_begin, None, "U1=None（部分食）");
        assert_eq!(gc.central_end, None, "U4=None（部分食）");
    }

    /// `GlobalCircumstances` は `Copy`・`PartialEq`。中心食と部分食の条件は `ne`、同値は `eq`。
    /// 接触 Option の Some/None 差・gamma 差・kind 差が `ne` に表れる。
    /// 殺す変異: `#[derive(Copy)]`/`#[derive(PartialEq)]` の脱落。
    #[test]
    fn global_circumstances_is_copy_and_partial_eq() {
        fn assert_copy<T: Copy>(_: T) {}
        let central = GlobalCircumstances {
            kind: SolarEclipseKind::Total,
            partial_begin: Some(contact(16, 10.0)),
            central_begin: Some(contact(17, 20.0)),
            greatest: greatest_central(),
            central_end: Some(contact(19, 40.0)),
            partial_end: Some(contact(20, 50.0)),
            gamma: 0.3431,
        };
        let copied = central; // Copy
        assert_copy(copied);
        assert_eq!(central, copied, "Copy 後も等しい（central は有効）");
        assert_eq!(central, copied, "同値は eq");

        let partial = GlobalCircumstances {
            kind: SolarEclipseKind::Partial,
            partial_begin: Some(contact(10, 60.0)),
            central_begin: None,
            greatest: greatest_partial(),
            central_end: None,
            partial_end: Some(contact(12, 62.0)),
            gamma: 1.21,
        };
        assert_ne!(central, partial, "中心食条件 と 部分食条件 は ne");

        // gamma だけ違えても ne（gamma がフィールドとして比較に効く）。
        let mut gamma_diff = central;
        gamma_diff.gamma = 0.9999;
        assert_ne!(central, gamma_diff, "gamma 違いは ne");
    }

    // ==================================================================
    // SolarEclipse: 5 フィールド保持 + Clone 後一致（PartialEq なし）
    // ==================================================================

    /// `SolarEclipse` の 5 フィールド（event_key / kind / global / bessel / metadata）を
    /// 構築値どおり保持する。`SolarEclipse` は PartialEq 非実装のため、比較可能なフィールド
    /// （event_key:String, kind:PartialEq, global:PartialEq, metadata:PartialEq）を個別に検証し、
    /// bessel は pub フィールド経由（fit_error 等）で値を読み出して確認する。
    /// 殺す変異: event_key/kind/global/metadata の取り違え・bessel フィールドの差し替え。
    #[test]
    fn solar_eclipse_holds_each_field() {
        let key = "2024-04-08#1252".to_string();
        let global = GlobalCircumstances {
            kind: SolarEclipseKind::Total,
            partial_begin: Some(contact(16, 10.0)),
            central_begin: Some(contact(17, 20.0)),
            greatest: greatest_central(),
            central_end: Some(contact(19, 40.0)),
            partial_end: Some(contact(20, 50.0)),
            gamma: 0.3431,
        };
        let bessel = minimal_bessel();
        let meta = metadata();
        let se = SolarEclipse {
            event_key: key.clone(),
            kind: SolarEclipseKind::Total,
            global,
            bessel: bessel.clone(),
            metadata: meta.clone(),
        };
        assert_eq!(se.event_key, key, "event_key（安定キー文字列）を保持");
        assert_eq!(se.kind, SolarEclipseKind::Total, "kind を保持");
        assert_eq!(se.global, global, "global を保持");
        assert_eq!(se.metadata, meta, "metadata を保持");
        // bessel は PartialEq 非実装。pub フィールド経由で構築値が読めることを縛る。
        assert_eq!(
            se.bessel.epoch_tt, bessel.epoch_tt,
            "bessel.epoch_tt を保持"
        );
        assert_eq!(se.bessel.tan_f1, bessel.tan_f1, "bessel.tan_f1 を保持");
        assert_eq!(
            se.bessel.fit_error, bessel.fit_error,
            "bessel.fit_error を保持"
        );
        assert_eq!(
            se.bessel.x.coefficients, bessel.x.coefficients,
            "bessel.x 多項式係数を保持"
        );
    }

    /// `SolarEclipse` は `Clone`。Clone 後も全フィールドが一致する（event_key 文字列・metadata
    /// （String フィールド含む）・global・bessel の pub フィールド）。
    /// 殺す変異: `#[derive(Clone)]` の脱落（コンパイルで露見）、Clone がフィールドを取りこぼす。
    #[test]
    fn solar_eclipse_clone_preserves_all_fields() {
        let global = GlobalCircumstances {
            kind: SolarEclipseKind::Annular,
            partial_begin: Some(contact(8, 5.0)),
            central_begin: Some(contact(9, 6.0)),
            greatest: greatest_central(),
            central_end: Some(contact(11, 7.0)),
            partial_end: Some(contact(12, 8.0)),
            gamma: -0.4,
        };
        let se = SolarEclipse {
            event_key: "2023-10-14#1248".to_string(),
            kind: SolarEclipseKind::Annular,
            global,
            bessel: minimal_bessel(),
            metadata: metadata(),
        };
        let cloned = se.clone();
        // 文字列キー（Clone で複製される値）。
        assert_eq!(cloned.event_key, se.event_key, "event_key が Clone で一致");
        assert_eq!(cloned.kind, se.kind, "kind が Clone で一致");
        assert_eq!(cloned.global, se.global, "global が Clone で一致");
        // metadata は String フィールド込みで一致。
        assert_eq!(cloned.metadata, se.metadata, "metadata が Clone で一致");
        // bessel は pub フィールド経由（PartialEq なし）。
        assert_eq!(
            cloned.bessel.fit_interval, se.bessel.fit_interval,
            "bessel.fit_interval が Clone で一致"
        );
        assert_eq!(
            cloned.bessel.l2.coefficients, se.bessel.l2.coefficients,
            "bessel.l2 係数が Clone で一致"
        );
        // 原本も有効（Clone であって move でない）。
        assert_eq!(se.event_key, "2023-10-14#1248", "原本 se は Clone 後も有効");
    }

    // ==================================================================
    // 局地 result 型（ISSUE-043 S4e, api-draft §3.4 rich 版）
    // ==================================================================
    //
    // ## オラクル戦略（S4d と同方針）
    // 局地型もすべて **pub フィールドのデータコンテナ**。振る舞いは「struct リテラルで構築した
    // 各フィールドを取り違えず保持・読み出す」「derive が仕様どおり付く」「Option/非 Option の
    // 区別」に尽きる。よって:
    //   1. 各フィールドを **互いに異なる非対称値**で構築し read-back（フィールド取り違え変異を撃破）。
    //   2. `maximum` は **非 Option**（型レベルで `LocalContact`）、c1..c4 は `Option`。部分食地点で
    //      c2/c3=None・c1/c4=Some、中心食地点で全 Some の 2 パターンを構築（Option 取り違え・
    //      maximum の Option 化変異を撃破）。
    //   3. 異値で `ne`・同値で `eq`、`Copy`/`Clone` をコンパイル時境界で縛る（derive 脱落を撃破）。
    //   4. `VisibleSolarEclipse` は PartialEq 非実装（SolarEclipse が PartialEq なし）。Clone 後に
    //      比較可能フィールドを個別照合する。

    // ------------------------------------------------------------------
    // 局地用 構築ヘルパ（非対称な既知値で各部品を作る）
    // ------------------------------------------------------------------

    /// `LocalContact` を **6 フィールドすべて異なる値**で組む（取り違え判別用）。
    /// `seed` で各フィールドが互いにずれるようにし、`vis` で visible を真偽切替する。
    /// sun_altitude / sun_azimuth / position_angle は **互いに異なる値域**を選び、
    /// time_utc と time_tt は別スケール（UTC は暦・TT は JD）で区別する。
    fn local_contact(seed: f64, vis: bool) -> LocalContact {
        LocalContact {
            // 暦由来（秒だけ seed に連動）。
            time_utc: utc(2024, 4, 8, 18, 17, seed),
            // JD 由来・別スケール（time_utc とは取り違え不能な値）。
            time_tt: tt(2_460_409.0, 0.1 + seed * 0.001),
            // 高度・方位・位置角は互いに大きく異なる値域（取り違えを撃破）。
            sun_altitude: Degrees(10.0 + seed),
            sun_azimuth: Degrees(200.0 + seed),
            position_angle: Degrees(300.0 + seed),
            visible: vis,
        }
    }

    /// 中心食地点の `LocalContactSet`（c1..c4 すべて Some・maximum は値）。
    /// 5 接触を **互いに異なる time（seed）** で組み、時系列の取り違え（c1↔c4, c2↔c3）を撃破する。
    fn local_set_central() -> LocalContactSet {
        LocalContactSet {
            c1: Some(local_contact(1.0, true)),
            c2: Some(local_contact(2.0, true)),
            maximum: local_contact(3.0, true),
            c3: Some(local_contact(4.0, true)),
            c4: Some(local_contact(5.0, true)),
        }
    }

    /// 部分食地点の `LocalContactSet`（c2/c3 は None・c1/c4 は Some・maximum は値）。
    /// 部分食では内接（C2/C3）が無いため None、外接（C1/C4）と最大食のみ。
    fn local_set_partial() -> LocalContactSet {
        LocalContactSet {
            c1: Some(local_contact(11.0, true)),
            c2: None,
            maximum: local_contact(13.0, false),
            c3: None,
            c4: Some(local_contact(15.0, true)),
        }
    }

    /// 代表的な `LocalCircumstances`（6 フィールド非対称・magnitude≠obscuration）。
    fn local_circumstances() -> LocalCircumstances {
        LocalCircumstances {
            contacts: local_set_central(),
            magnitude: EclipseMagnitude(1.0234),
            obscuration: Obscuration(0.987),
            maximum_altitude: Degrees(55.5),
            visibility: Visibility::FullyVisible,
            metadata: metadata(),
        }
    }

    // ==================================================================
    // LocalContact: 6 フィールド保持・取り違えない / Copy・PartialEq
    // ==================================================================

    /// `LocalContact` の 6 フィールド（time_utc / time_tt / sun_altitude / sun_azimuth /
    /// position_angle / visible）が構築値どおり保持され、互いに取り違えられない。
    /// 高度・方位・位置角は **異なる値域**（10°台 / 200°台 / 300°台）で与え、3 つの Degrees
    /// フィールドの相互取り違えを撃破。time_utc≠time_tt（別スケール）、visible 真偽も縛る。
    /// 殺す変異: sun_altitude↔sun_azimuth↔position_angle の入れ替え、time_utc↔time_tt の入れ替え、
    /// visible の固定化（常に true/false）。
    #[test]
    fn local_contact_holds_each_field() {
        let t_utc = utc(2024, 4, 8, 18, 17, 6.5);
        let t_tt = tt(2_460_409.0, 0.234);
        let c = LocalContact {
            time_utc: t_utc,
            time_tt: t_tt,
            sun_altitude: Degrees(42.0),
            sun_azimuth: Degrees(210.0),
            position_angle: Degrees(305.0),
            visible: true,
        };
        assert_eq!(c.time_utc, t_utc, "time_utc を保持");
        assert_eq!(c.time_tt, t_tt, "time_tt を保持");
        assert_eq!(c.sun_altitude, Degrees(42.0), "sun_altitude を保持");
        assert_eq!(c.sun_azimuth, Degrees(210.0), "sun_azimuth を保持");
        assert_eq!(c.position_angle, Degrees(305.0), "position_angle を保持");
        assert!(c.visible, "visible=true を保持");
        // 3 つの角度フィールドは互いに別値（取り違え変異を撃破）。
        assert_ne!(
            c.sun_altitude, c.sun_azimuth,
            "sun_altitude と sun_azimuth は別フィールド・別値"
        );
        assert_ne!(
            c.sun_azimuth, c.position_angle,
            "sun_azimuth と position_angle は別フィールド・別値"
        );
        assert_ne!(
            c.sun_altitude, c.position_angle,
            "sun_altitude と position_angle は別フィールド・別値"
        );

        // visible=false 側も保持する（真偽の固定化変異を撃破）。
        let invisible = LocalContact {
            visible: false,
            ..c
        };
        assert!(!invisible.visible, "visible=false を保持");
    }

    /// `LocalContact` は `Copy`・`PartialEq`。同値は `eq`、各フィールド違いは `ne`。
    /// 殺す変異: `#[derive(Copy)]`/`#[derive(PartialEq)]` の脱落、visible を比較から外す。
    #[test]
    fn local_contact_is_copy_and_partial_eq() {
        fn assert_copy<T: Copy>(_: T) {}
        let a = local_contact(1.0, true);
        let b = a; // Copy（move されない）
        assert_copy(b);
        assert_eq!(a, b, "Copy 後も等しい（a は有効）");
        // 同一構築値は eq。
        assert_eq!(
            local_contact(1.0, true),
            local_contact(1.0, true),
            "同一構築値は eq"
        );
        // seed 違い（time/角度すべてずれる）は ne。
        assert_ne!(
            local_contact(1.0, true),
            local_contact(2.0, true),
            "seed 違いは ne"
        );
        // visible だけ違っても ne（visible が比較に効く）。
        assert_ne!(
            local_contact(1.0, true),
            local_contact(1.0, false),
            "visible 違いは ne"
        );
    }

    // ==================================================================
    // LocalContactSet: c1/c2/maximum/c3/c4 保持・maximum は非 Option
    // ==================================================================

    /// 中心食地点の `LocalContactSet`: c1..c4 がすべて Some、maximum は値。**5 接触を互いに
    /// 異なる time** で構築し、各フィールドから正しい接触を read-back（時系列の取り違えを撃破）。
    /// 殺す変異: c1↔c4 / c2↔c3 の入れ替え、maximum と他接触の取り違え。
    #[test]
    fn local_contact_set_central_holds_all_contacts() {
        let c1 = local_contact(1.0, true);
        let c2 = local_contact(2.0, true);
        let mx = local_contact(3.0, true);
        let c3 = local_contact(4.0, true);
        let c4 = local_contact(5.0, true);
        let set = LocalContactSet {
            c1: Some(c1),
            c2: Some(c2),
            maximum: mx,
            c3: Some(c3),
            c4: Some(c4),
        };
        assert_eq!(set.c1, Some(c1), "c1 を保持");
        assert_eq!(set.c2, Some(c2), "c2 を保持");
        assert_eq!(set.maximum, mx, "maximum を保持（非 Option で値）");
        assert_eq!(set.c3, Some(c3), "c3 を保持");
        assert_eq!(set.c4, Some(c4), "c4 を保持");
        // 時系列の取り違え検出（互いに異なる time なので相互に ne）。
        assert_ne!(set.c1, set.c4, "c1 と c4 は別接触（時系列取り違え撃破）");
        assert_ne!(set.c2, set.c3, "c2 と c3 は別接触（時系列取り違え撃破）");
        assert_ne!(set.maximum, c1, "maximum と c1 は別接触");
    }

    /// 部分食地点の `LocalContactSet`: c2/c3 が None（内接なし）、c1/c4 は Some、maximum は値。
    /// **maximum はどの地点でも非 Option で存在する**こと（部分食でも値を持つ）を縛る。
    /// 殺す変異: 部分食で c2/c3 を常に Some にする、maximum を Option 化して部分食で None にする、
    /// c1/c4 を None と取り違える。
    #[test]
    fn local_contact_set_partial_has_none_inner_contacts() {
        let c1 = local_contact(11.0, true);
        let mx = local_contact(13.0, false);
        let c4 = local_contact(15.0, true);
        let set = LocalContactSet {
            c1: Some(c1),
            c2: None,
            maximum: mx,
            c3: None,
            c4: Some(c4),
        };
        assert_eq!(set.c1, Some(c1), "部分食 c1=Some");
        assert_eq!(set.c2, None, "部分食 c2=None（内接なし）");
        assert_eq!(set.c3, None, "部分食 c3=None（内接なし）");
        assert_eq!(set.c4, Some(c4), "部分食 c4=Some");
        // maximum は部分食でも非 Option の値として存在（Option 化変異を撃破）。
        assert_eq!(
            set.maximum, mx,
            "部分食でも maximum は値（非 Option・常に存在）"
        );
    }

    /// `LocalContactSet` は `Copy`・`PartialEq`。中心食集合と部分食集合（c2/c3 の Some/None 差）は
    /// `ne`、同値は `eq`。
    /// 殺す変異: `#[derive(Copy)]`/`#[derive(PartialEq)]` の脱落、Option フィールドを比較から外す。
    #[test]
    fn local_contact_set_is_copy_and_partial_eq() {
        fn assert_copy<T: Copy>(_: T) {}
        let central = local_set_central();
        let copied = central; // Copy
        assert_copy(copied);
        assert_eq!(central, copied, "Copy 後も等しい（central は有効）");
        assert_eq!(local_set_central(), local_set_central(), "同一構築値は eq");
        // 中心食（c2/c3=Some）と部分食（c2/c3=None）は ne。
        assert_ne!(
            local_set_central(),
            local_set_partial(),
            "中心食集合 と 部分食集合 は ne（Some/None 差）"
        );

        // maximum だけ違っても ne（maximum が比較に効く）。seed は秒域 [0,60) 内の別値。
        let mut max_diff = central;
        max_diff.maximum = local_contact(42.0, true);
        assert_ne!(central, max_diff, "maximum 違いは ne");
    }

    // ==================================================================
    // LocalCircumstances: 6 フィールド保持 / Clone・PartialEq（Copy 不可）
    // ==================================================================

    /// `LocalCircumstances` の 6 フィールド（contacts / magnitude / obscuration /
    /// maximum_altitude / visibility / metadata）が構築値どおり保持される。
    /// magnitude(1.0234) と obscuration(0.987) を **異なる値**で与え、両者の取り違えを撃破。
    /// 殺す変異: magnitude↔obscuration の取り違え、各フィールドの read-back 入れ替え。
    #[test]
    fn local_circumstances_holds_each_field() {
        let contacts = local_set_central();
        let meta = metadata();
        let lc = LocalCircumstances {
            contacts,
            magnitude: EclipseMagnitude(1.0234),
            obscuration: Obscuration(0.987),
            maximum_altitude: Degrees(55.5),
            visibility: Visibility::FullyVisible,
            metadata: meta.clone(),
        };
        assert_eq!(lc.contacts, contacts, "contacts を保持");
        assert_eq!(lc.magnitude, EclipseMagnitude(1.0234), "magnitude を保持");
        assert_eq!(lc.obscuration, Obscuration(0.987), "obscuration を保持");
        // magnitude と obscuration は別フィールド・別値（取り違え撃破）。
        assert!(
            (lc.magnitude.0 - lc.obscuration.0).abs() > 1e-9,
            "magnitude と obscuration は別フィールド・別値であること"
        );
        assert_eq!(
            lc.maximum_altitude,
            Degrees(55.5),
            "maximum_altitude を保持"
        );
        assert_eq!(lc.visibility, Visibility::FullyVisible, "visibility を保持");
        assert_eq!(lc.metadata, meta, "metadata を保持");
    }

    /// `LocalCircumstances` は `Clone`・`PartialEq`（CalculationMetadata の String ゆえ Copy 不可）。
    /// Clone 後も全フィールドが一致し、異なる circumstances は `ne`、同値は `eq`。
    /// 殺す変異: `#[derive(Clone)]`/`#[derive(PartialEq)]` の脱落、フィールドの取りこぼし。
    #[test]
    fn local_circumstances_is_clone_and_partial_eq() {
        let a = local_circumstances();
        let cloned = a.clone();
        assert_eq!(a, cloned, "Clone 後も等しい");
        assert_eq!(a.metadata, cloned.metadata, "metadata（String 込み）が一致");
        // 原本も有効（Clone であって move でない）。
        assert_eq!(
            a.maximum_altitude,
            Degrees(55.5),
            "原本 a は Clone 後も有効"
        );
        // 同一構築値は eq。
        assert_eq!(local_circumstances(), local_circumstances(), "同値は eq");

        // visibility だけ違えば ne（visibility が比較に効く）。
        let mut vis_diff = a.clone();
        vis_diff.visibility = Visibility::PartialVisible;
        assert_ne!(a, vis_diff, "visibility 違いは ne");

        // contacts（中心食→部分食）が違えば ne。
        let mut contacts_diff = a.clone();
        contacts_diff.contacts = local_set_partial();
        assert_ne!(a, contacts_diff, "contacts 違いは ne");
    }

    // ==================================================================
    // VisibleSolarEclipse: eclipse/local 保持 + Clone 後一致（PartialEq なし）
    // ==================================================================

    /// `VisibleSolarEclipse` の 2 フィールド（eclipse / local）を構築値どおり保持する。
    /// `VisibleSolarEclipse` は PartialEq 非実装（SolarEclipse が PartialEq なし）のため、
    /// 比較可能な内部フィールドを個別照合する（eclipse の event_key・local の各フィールド）。
    /// 殺す変異: eclipse↔local のフィールド取り違え、フィールドの差し替え。
    #[test]
    fn visible_solar_eclipse_holds_each_field() {
        let global = GlobalCircumstances {
            kind: SolarEclipseKind::Total,
            partial_begin: Some(contact(16, 10.0)),
            central_begin: Some(contact(17, 20.0)),
            greatest: greatest_central(),
            central_end: Some(contact(19, 40.0)),
            partial_end: Some(contact(20, 50.0)),
            gamma: 0.3431,
        };
        let eclipse = SolarEclipse {
            event_key: "2024-04-08#1252".to_string(),
            kind: SolarEclipseKind::Total,
            global,
            bessel: minimal_bessel(),
            metadata: metadata(),
        };
        let local = local_circumstances();
        let vse = VisibleSolarEclipse {
            eclipse: eclipse.clone(),
            local: local.clone(),
        };
        // eclipse フィールド（PartialEq なし）の比較可能な部分を照合。
        assert_eq!(
            vse.eclipse.event_key, "2024-04-08#1252",
            "eclipse.event_key を保持"
        );
        assert_eq!(
            vse.eclipse.kind,
            SolarEclipseKind::Total,
            "eclipse.kind を保持"
        );
        assert_eq!(vse.eclipse.global, global, "eclipse.global を保持");
        // local フィールド（LocalCircumstances は PartialEq）を丸ごと照合。
        assert_eq!(vse.local, local, "local を保持");
    }

    /// `VisibleSolarEclipse` は `Clone`。Clone 後も eclipse・local が一致する。
    /// PartialEq は無い前提（SolarEclipse が PartialEq 非実装）でフィールド個別比較する。
    /// 殺す変異: `#[derive(Clone)]` の脱落（コンパイルで露見）、Clone がフィールドを取りこぼす。
    #[test]
    fn visible_solar_eclipse_clone_preserves_fields() {
        let global = GlobalCircumstances {
            kind: SolarEclipseKind::Annular,
            partial_begin: Some(contact(8, 5.0)),
            central_begin: Some(contact(9, 6.0)),
            greatest: greatest_central(),
            central_end: Some(contact(11, 7.0)),
            partial_end: Some(contact(12, 8.0)),
            gamma: -0.4,
        };
        let vse = VisibleSolarEclipse {
            eclipse: SolarEclipse {
                event_key: "2023-10-14#1248".to_string(),
                kind: SolarEclipseKind::Annular,
                global,
                bessel: minimal_bessel(),
                metadata: metadata(),
            },
            local: local_circumstances(),
        };
        let cloned = vse.clone();
        // eclipse（PartialEq なし）は比較可能フィールドで照合。
        assert_eq!(
            cloned.eclipse.event_key, vse.eclipse.event_key,
            "eclipse.event_key が Clone で一致"
        );
        assert_eq!(
            cloned.eclipse.kind, vse.eclipse.kind,
            "eclipse.kind が Clone で一致"
        );
        assert_eq!(
            cloned.eclipse.global, vse.eclipse.global,
            "eclipse.global が Clone で一致"
        );
        // local（PartialEq あり）は丸ごと照合。
        assert_eq!(cloned.local, vse.local, "local が Clone で一致");
        // 原本も有効（Clone であって move でない）。
        assert_eq!(
            vse.eclipse.event_key, "2023-10-14#1248",
            "原本 vse は Clone 後も有効"
        );
    }
}
