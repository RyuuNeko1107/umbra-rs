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
use crate::magnitude::{EclipseMagnitude, Obscuration};

/// 全球接触点（時刻 TT/UTC ＋ 地表点）。
#[derive(Clone, Copy, Debug, PartialEq)]
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
    /// 帯幅 \[km\]（中心食のみ Some）。
    pub path_width: Option<Kilometers>,
    /// 中心食継続時間 \[s\]（中心食のみ Some）。
    pub central_duration: Option<f64>,
    /// 最大食地点での太陽高度。
    pub sun_altitude: Degrees,
}

/// 全球条件（種別・P1/U1/最大食/U4/P4・gamma）。
#[derive(Clone, Copy, Debug, PartialEq)]
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
    use crate::magnitude::{EclipseMagnitude, Obscuration};
    use crate::polynomial::Polynomial;
    use crate::results::{GlobalCircumstances, GlobalContact, GreatestEclipse, SolarEclipse};

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
}
