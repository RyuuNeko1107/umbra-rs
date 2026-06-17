//! 解析暦バックエンド（VSOP87D 太陽 + ELP2000-82B 月）の [`Ephemeris`] 実装（ISSUE-043 S1）。
//!
//! [`crate::apparent`] / [`crate::sun`] の既存幾何位置関数を [`Ephemeris`] 形へ**薄くラップ**する
//! だけで、新しい天文計算は行わない。TT≈TDB として `TtInstant::from_jd2(time.jd2())` で既存関数を
//! 呼ぶ。Origin は Geocenter のみサポート（SolarSystemBarycenter は [`EphemerisError::DataUnavailable`]）。
//!
//! [`Ephemeris`]: crate::ephemeris::Ephemeris
//! [`EphemerisError::DataUnavailable`]: crate::ephemeris::EphemerisError::DataUnavailable

use umbra_core::constants::ASTRONOMICAL_UNIT_KM;
use umbra_core::{JulianDate2, TdbInstant, TimeRange, TtInstant, Vector3};

use crate::apparent::{moon_geocentric_gcrs, sun_geocentric_gcrs};
use crate::ephemeris::{
    Body, Ephemeris, EphemerisError, EphemerisFrame, EphemerisMetadata, Origin, StateVector,
};
use crate::frames::ecliptic_to_gcrs_matrix;
use crate::sun::{earth_heliocentric_velocity_ecliptic_of_date, sun_geocentric_ecliptic_of_date};

/// 対応 TDB 範囲の下端＝1900-01-01 0h（JD）。ELP2000-82B の実用域に合わせる。
const SUPPORTED_START_JD: f64 = 2_415_020.5;
/// 対応 TDB 範囲の上端＝2100-01-01 0h（JD）。
const SUPPORTED_END_JD: f64 = 2_488_069.5;

/// VSOP87D（太陽）+ ELP2000-82B（月）の解析暦バックエンド（幾何位置）。
///
/// [`crate::apparent`] / [`crate::sun`] の既存幾何位置関数を [`Ephemeris`] 形に薄くラップする
/// 純粋型（保持状態なし）。見かけ補正（光行時間・光行差・歳差章動）は本型では行わず、上位の
/// apparent 層（ISSUE-043 S2 でジェネリック化）が `state()` の幾何位置に適用する。
#[derive(Debug, Clone, Copy, Default)]
pub struct AnalyticalEphemeris;

impl AnalyticalEphemeris {
    /// 解析暦バックエンドを構築する（保持状態なし。[`Default`] と同一）。
    pub fn new() -> Self {
        Self
    }
}

impl Ephemeris for AnalyticalEphemeris {
    /// `body` の `time`（TT≈TDB）における幾何状態を返す。**地心のみ**供給（SSB は
    /// [`EphemerisError::DataUnavailable`]）。月の黄道 of date・EMB は v0.1 未提供。
    fn state(
        &self,
        body: Body,
        time: TdbInstant,
        origin: Origin,
        frame: EphemerisFrame,
    ) -> Result<StateVector, EphemerisError> {
        // 解析暦は地心位置のみ供給する（太陽系重心基準は JPL DE バックエンドの責務）。
        if origin != Origin::Geocenter {
            return Err(EphemerisError::DataUnavailable);
        }
        // 既存の幾何位置関数は TtInstant 入力。TT≈TDB として同一 JD で評価する。
        let tt = TtInstant::from_jd2(time.jd2());
        match (body, frame) {
            // 地心基準の地球は原点（速度ゼロ）。frame に依らない。
            (Body::Earth, _) => Ok(StateVector {
                position: Vector3::ZERO,
                velocity: Some(Vector3::ZERO),
            }),
            // 太陽 GCRS（ICRS）幾何位置。地心太陽速度 = −(地球日心速度) を GCRS へ回す。
            // 回転は `tt`（TT）、速度評価は `time`（TDB）だが TT≈TDB かつ同一 JD のため整合する
            // （位置側 `sun_geocentric_gcrs` の内部評価と同じ JD）。
            (Body::Sun, EphemerisFrame::Icrs) => {
                let earth_v_gcrs = ecliptic_to_gcrs_matrix(tt)
                    .mul_vec(earth_heliocentric_velocity_ecliptic_of_date(time));
                Ok(StateVector {
                    position: sun_geocentric_gcrs(tt),
                    velocity: Some(earth_v_gcrs.scale(-1.0)),
                })
            }
            // 太陽 黄道 of date 幾何位置（AU→km）。速度も同フレームの −(地球日心速度)。
            (Body::Sun, EphemerisFrame::EclipticOfDate) => Ok(StateVector {
                position: sun_geocentric_ecliptic_of_date(time).scale(ASTRONOMICAL_UNIT_KM),
                velocity: Some(earth_heliocentric_velocity_ecliptic_of_date(time).scale(-1.0)),
            }),
            // 月 GCRS（ICRS）幾何位置。月速度は未実装（None）。
            (Body::Moon, EphemerisFrame::Icrs) => Ok(StateVector {
                position: moon_geocentric_gcrs(tt),
                velocity: None,
            }),
            // 月の黄道 of date は J2000→of date 変換が未実装（v0.1 非対応）。
            (Body::Moon, EphemerisFrame::EclipticOfDate) => Err(EphemerisError::DataUnavailable),
            // EMB は解析暦 v0.1 では未提供。
            (Body::EarthMoonBarycenter, _) => Err(EphemerisError::DataUnavailable),
        }
    }

    /// 対応 TDB 範囲 [1900-01-01, 2100-01-01]（ELP2000-82B の実用域が律速）。
    fn supported_range(&self) -> TimeRange<TdbInstant> {
        TimeRange {
            start: TdbInstant::from_jd2(JulianDate2::from_jd(SUPPORTED_START_JD)),
            end: TdbInstant::from_jd2(JulianDate2::from_jd(SUPPORTED_END_JD)),
        }
    }

    /// バックエンドのメタデータ（`CalculationMetadata` へ転記）。最大残差は M10（DE 差分）で確定。
    fn metadata(&self) -> EphemerisMetadata {
        EphemerisMetadata {
            model: "VSOP87D + ELP2000-82B".to_string(),
            version: "VSOP87D earth (full) + ELP2000-82B (full); see generated artifacts".to_string(),
            source: "VSOP87 (IMCCE / Bretagnon & Francou) + ELP2000-82B (IMCCE / Chapront-Touzé & Chapront)"
                .to_string(),
            license: "IMCCE published scientific data (attribution); not a GPL derivative"
                .to_string(),
            supported: self.supported_range(),
            // 達成残差は未測定（M10 の DE 差分で確定）。
            max_residual_arcsec: f64::NAN,
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::analytical::AnalyticalEphemeris;
    use crate::apparent::{moon_geocentric_gcrs, sun_geocentric_gcrs};
    use crate::ephemeris::{Body, Ephemeris, EphemerisError, EphemerisFrame, Origin, StateVector};
    use crate::frames::ecliptic_to_gcrs_matrix;
    use crate::sun::{
        earth_heliocentric_velocity_ecliptic_of_date, sun_geocentric_ecliptic_of_date,
    };
    use umbra_core::constants::ASTRONOMICAL_UNIT_KM;
    use umbra_core::{JulianDate2, TdbInstant, TtInstant, Vector3};

    // ------------------------------------------------------------------
    // 検証エポックとヘルパ。
    //   J2000 (2451545.0) と 2017-08-21 皆既日食付近 (≈2457987.0) の 2 点で、
    //   位置・速度の組立同一性を確認する（時刻依存・回転日の取り違えを励起）。
    // ------------------------------------------------------------------

    const JD_J2000: f64 = 2_451_545.0;
    /// 2017-08-21（北米皆既日食）付近。
    const JD_2017_ECLIPSE: f64 = 2_457_987.0;
    /// 検証エポック一覧。
    const EPOCHS: [f64; 2] = [JD_J2000, JD_2017_ECLIPSE];

    fn tdb(jd: f64) -> TdbInstant {
        TdbInstant::from_jd2(JulianDate2::new(jd, 0.0))
    }

    /// TT≈TDB として同 JD の TtInstant を構築（既存関数呼び出し用）。
    fn tt(jd: f64) -> TtInstant {
        TtInstant::from_jd2(JulianDate2::new(jd, 0.0))
    }

    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    fn vec_close(a: Vector3, b: Vector3, tol: f64) -> bool {
        close(a.x, b.x, tol) && close(a.y, b.y, tol) && close(a.z, b.z, tol)
    }

    /// state を取り、Ok の StateVector を返す（unwrap）。
    fn state(
        eph: &AnalyticalEphemeris,
        body: Body,
        jd: f64,
        origin: Origin,
        frame: EphemerisFrame,
    ) -> StateVector {
        eph.state(body, tdb(jd), origin, frame).unwrap()
    }

    // ==================================================================
    // 1. サポート組: 位置の組立同一性（既存検証済み関数がオラクル）
    // ==================================================================

    /// Sun/Geocenter/Icrs の position == apparent::sun_geocentric_gcrs(tt)（成分一致）。
    /// 殺す変異: 別関数を呼ぶ・frame 取り違え（Ecliptic を返す）・body 取り違え（Moon 値）。
    #[test]
    fn sun_geocenter_icrs_position_equals_sun_geocentric_gcrs() {
        let eph = AnalyticalEphemeris::new();
        for &jd in &EPOCHS {
            let got = state(&eph, Body::Sun, jd, Origin::Geocenter, EphemerisFrame::Icrs).position;
            let expected = sun_geocentric_gcrs(tt(jd));
            assert!(
                vec_close(got, expected, 1e-9),
                "sun icrs pos(jd={jd}) = {got:?}, expected {expected:?}"
            );
        }
    }

    /// Moon/Geocenter/Icrs の position == apparent::moon_geocentric_gcrs(tt)（成分一致）。
    /// 殺す変異: body 取り違え（Sun 値を返す）・別関数。
    #[test]
    fn moon_geocenter_icrs_position_equals_moon_geocentric_gcrs() {
        let eph = AnalyticalEphemeris::new();
        for &jd in &EPOCHS {
            let got = state(
                &eph,
                Body::Moon,
                jd,
                Origin::Geocenter,
                EphemerisFrame::Icrs,
            )
            .position;
            let expected = moon_geocentric_gcrs(tt(jd));
            assert!(
                vec_close(got, expected, 1e-9),
                "moon icrs pos(jd={jd}) = {got:?}, expected {expected:?}"
            );
        }
    }

    /// Earth/Geocenter は position=ZERO, velocity=Some(ZERO)（任意 frame で）。
    /// 殺す変異: 非ゼロ返却・velocity を None にする・frame 依存にする。
    #[test]
    fn earth_geocenter_is_zero_with_zero_velocity_any_frame() {
        let eph = AnalyticalEphemeris::new();
        for &frame in &[EphemerisFrame::Icrs, EphemerisFrame::EclipticOfDate] {
            for &jd in &EPOCHS {
                let s = state(&eph, Body::Earth, jd, Origin::Geocenter, frame);
                assert_eq!(s.position, Vector3::ZERO, "earth pos(jd={jd}, {frame:?})");
                assert_eq!(
                    s.velocity,
                    Some(Vector3::ZERO),
                    "earth vel(jd={jd}, {frame:?})"
                );
            }
        }
    }

    /// Sun/Geocenter/EclipticOfDate の position == sun_geocentric_ecliptic_of_date(tdb) * AU_KM。
    /// 殺す変異: AU→km 換算漏れ（AU のまま返す）・frame 取り違え（GCRS を返す）。
    #[test]
    fn sun_geocenter_ecliptic_position_equals_ecliptic_times_au() {
        let eph = AnalyticalEphemeris::new();
        for &jd in &EPOCHS {
            let got = state(
                &eph,
                Body::Sun,
                jd,
                Origin::Geocenter,
                EphemerisFrame::EclipticOfDate,
            )
            .position;
            let expected = sun_geocentric_ecliptic_of_date(tdb(jd)).scale(ASTRONOMICAL_UNIT_KM);
            assert!(
                vec_close(got, expected, 1e-6),
                "sun ecliptic pos(jd={jd}) = {got:?}, expected {expected:?}"
            );
        }
    }

    // ==================================================================
    // 2. サポート組: 速度の組立同一性
    // ==================================================================

    /// Sun/Geocenter/Icrs の velocity = -(ecliptic_to_gcrs_matrix(tt) *
    ///   earth_heliocentric_velocity_ecliptic_of_date(tdb))（既に km/s）。独立に組立てて一致。
    /// 殺す変異: 符号反転漏れ（地球日心速度のまま）・回転行列の日取り違え・None にする・km/s 換算誤り。
    #[test]
    fn sun_geocenter_icrs_velocity_equals_negated_earth_velocity_gcrs() {
        let eph = AnalyticalEphemeris::new();
        for &jd in &EPOCHS {
            let s = state(&eph, Body::Sun, jd, Origin::Geocenter, EphemerisFrame::Icrs);
            let earth_v_gcrs = ecliptic_to_gcrs_matrix(tt(jd))
                .mul_vec(earth_heliocentric_velocity_ecliptic_of_date(tdb(jd)));
            let expected = earth_v_gcrs.scale(-1.0);
            let got = s
                .velocity
                .unwrap_or_else(|| panic!("sun icrs velocity is Some (jd={jd})"));
            assert!(
                vec_close(got, expected, 1e-9),
                "sun icrs vel(jd={jd}) = {got:?}, expected {expected:?}"
            );
        }
    }

    /// Sun Icrs velocity の符号: 太陽地心速度は地球日心速度の符号反転（cos < 0）かつノルム ≈ 29-30 km/s。
    /// 殺す変異: 符号反転漏れ（同方向 cos>0 になる）・ゼロ速度。
    #[test]
    fn sun_geocenter_icrs_velocity_sign_and_magnitude() {
        let eph = AnalyticalEphemeris::new();
        for &jd in &EPOCHS {
            let got = state(&eph, Body::Sun, jd, Origin::Geocenter, EphemerisFrame::Icrs)
                .velocity
                .unwrap();
            let earth_v_gcrs = ecliptic_to_gcrs_matrix(tt(jd))
                .mul_vec(earth_heliocentric_velocity_ecliptic_of_date(tdb(jd)));
            // 太陽地心速度は地球日心速度と逆向き。
            let cos = got.dot(earth_v_gcrs) / (got.norm() * earth_v_gcrs.norm());
            assert!(
                cos < -0.99,
                "sun vel direction(jd={jd}) cos = {cos}, want ≈ -1"
            );
            // 地球公転速度オーダー 29-30 km/s。
            assert!(
                (29.0..30.5).contains(&got.norm()),
                "sun geocentric speed(jd={jd}) = {} km/s out of [29,30.5]",
                got.norm()
            );
        }
    }

    /// Moon/Geocenter/Icrs の velocity == None（月速度は未実装）。
    /// 殺す変異: Some を返す（Sun と同様に velocity を埋める）。
    #[test]
    fn moon_geocenter_icrs_velocity_is_none() {
        let eph = AnalyticalEphemeris::new();
        for &jd in &EPOCHS {
            let s = state(
                &eph,
                Body::Moon,
                jd,
                Origin::Geocenter,
                EphemerisFrame::Icrs,
            );
            assert_eq!(s.velocity, None, "moon icrs velocity(jd={jd}) must be None");
        }
    }

    /// Sun/Geocenter/EclipticOfDate の velocity == -(earth_heliocentric_velocity_ecliptic_of_date)
    /// （同フレームで符号反転のみ・回転なし）。
    /// 殺す変異: 符号反転漏れ（地球日心速度のまま正方向）・None 化・回転を誤って掛ける。
    #[test]
    fn sun_geocenter_ecliptic_velocity_equals_negated_earth_velocity() {
        let eph = AnalyticalEphemeris::new();
        for &jd in &EPOCHS {
            let got = state(
                &eph,
                Body::Sun,
                jd,
                Origin::Geocenter,
                EphemerisFrame::EclipticOfDate,
            )
            .velocity
            .unwrap_or_else(|| panic!("sun ecliptic velocity is Some (jd={jd})"));
            let expected = earth_heliocentric_velocity_ecliptic_of_date(tdb(jd)).scale(-1.0);
            assert!(
                vec_close(got, expected, 1e-9),
                "sun ecliptic vel(jd={jd}) = {got:?}, expected {expected:?}"
            );
            // 符号反転されている（地球日心速度そのものではない）ことも固定する。
            let earth_v = earth_heliocentric_velocity_ecliptic_of_date(tdb(jd));
            assert!(
                got.dot(earth_v) < 0.0,
                "sun ecliptic vel(jd={jd}) must oppose earth heliocentric velocity"
            );
        }
    }

    // ==================================================================
    // 3. body/frame/origin 取り違えを殺す（識別）
    // ==================================================================

    /// Sun Icrs と Moon Icrs の position は明確に異なる（body 取り違えの相互混同を殺す）。
    #[test]
    fn sun_and_moon_icrs_positions_are_distinct() {
        let eph = AnalyticalEphemeris::new();
        for &jd in &EPOCHS {
            let sun = state(&eph, Body::Sun, jd, Origin::Geocenter, EphemerisFrame::Icrs).position;
            let moon = state(
                &eph,
                Body::Moon,
                jd,
                Origin::Geocenter,
                EphemerisFrame::Icrs,
            )
            .position;
            // 太陽 ≈ 1.5e8 km、月 ≈ 3.8e5 km。ノルム差で明確に分離。
            assert!(
                (sun - moon).norm() > 1e8,
                "sun/moon icrs positions too close(jd={jd}): sun={sun:?}, moon={moon:?}"
            );
        }
    }

    /// Sun Icrs と Sun EclipticOfDate の position は異なる（frame 取り違えを殺す。
    ///   GCRS と黄道 of date は回転で別ベクトル）。
    #[test]
    fn sun_icrs_and_ecliptic_positions_differ() {
        let eph = AnalyticalEphemeris::new();
        for &jd in &EPOCHS {
            let icrs = state(&eph, Body::Sun, jd, Origin::Geocenter, EphemerisFrame::Icrs).position;
            let ecl = state(
                &eph,
                Body::Sun,
                jd,
                Origin::Geocenter,
                EphemerisFrame::EclipticOfDate,
            )
            .position;
            // ノルムはほぼ同じ（回転）が、ベクトルは一致しない。
            assert!(
                !vec_close(icrs, ecl, 1.0),
                "sun icrs == ecliptic(jd={jd}) (frame ignored?): icrs={icrs:?}, ecl={ecl:?}"
            );
        }
    }

    // ==================================================================
    // 4. サニティ: 距離オーダー
    // ==================================================================

    /// Sun Icrs position のノルム ≈ 1 AU（1.4e8..1.6e8 km）。
    /// 殺す変異: スケール暴走・AU/km 取り違え・ゼロ返却。
    #[test]
    fn sun_icrs_distance_is_about_one_au() {
        let eph = AnalyticalEphemeris::new();
        for &jd in &EPOCHS {
            let r = state(&eph, Body::Sun, jd, Origin::Geocenter, EphemerisFrame::Icrs)
                .position
                .norm();
            assert!(
                (1.4e8..1.6e8).contains(&r),
                "sun icrs distance(jd={jd}) = {r} km out of [1.4e8,1.6e8]"
            );
        }
    }

    /// Moon Icrs position のノルム ≈ 月距離（356000..407000 km）。
    #[test]
    fn moon_icrs_distance_is_about_lunar_distance() {
        let eph = AnalyticalEphemeris::new();
        for &jd in &EPOCHS {
            let r = state(
                &eph,
                Body::Moon,
                jd,
                Origin::Geocenter,
                EphemerisFrame::Icrs,
            )
            .position
            .norm();
            assert!(
                (356_000.0..407_000.0).contains(&r),
                "moon icrs distance(jd={jd}) = {r} km out of [356000,407000]"
            );
        }
    }

    // ==================================================================
    // 5. 非サポート組: DataUnavailable
    // ==================================================================

    /// origin == SolarSystemBarycenter は全 body で DataUnavailable。
    /// 殺す変異: SSB を Geocenter と同一視して値を返す・OutOfSupportedRange を返す。
    #[test]
    fn solar_system_barycenter_is_unavailable_for_all_bodies() {
        let eph = AnalyticalEphemeris::new();
        let bodies = [
            Body::Sun,
            Body::Moon,
            Body::Earth,
            Body::EarthMoonBarycenter,
        ];
        for &frame in &[EphemerisFrame::Icrs, EphemerisFrame::EclipticOfDate] {
            for &body in &bodies {
                let r = eph.state(body, tdb(JD_J2000), Origin::SolarSystemBarycenter, frame);
                assert_eq!(
                    r,
                    Err(EphemerisError::DataUnavailable),
                    "SSB {body:?} {frame:?} must be DataUnavailable, got {r:?}"
                );
            }
        }
    }

    /// Moon/Geocenter/EclipticOfDate は DataUnavailable（of-date 変換未実装）。
    /// 殺す変異: Icrs と同じ値を返す・誤って Ok にする。
    #[test]
    fn moon_geocenter_ecliptic_is_unavailable() {
        let eph = AnalyticalEphemeris::new();
        for &jd in &EPOCHS {
            let r = eph.state(
                Body::Moon,
                tdb(jd),
                Origin::Geocenter,
                EphemerisFrame::EclipticOfDate,
            );
            assert_eq!(
                r,
                Err(EphemerisError::DataUnavailable),
                "moon geocenter ecliptic(jd={jd}) must be DataUnavailable, got {r:?}"
            );
        }
    }

    /// EarthMoonBarycenter は（Geocenter でも）全 frame で DataUnavailable（v0.1 EMB 未提供）。
    /// 殺す変異: 月位置や中点を返す・Ok にする。
    #[test]
    fn earth_moon_barycenter_is_unavailable() {
        let eph = AnalyticalEphemeris::new();
        for &frame in &[EphemerisFrame::Icrs, EphemerisFrame::EclipticOfDate] {
            for &jd in &EPOCHS {
                let r = eph.state(Body::EarthMoonBarycenter, tdb(jd), Origin::Geocenter, frame);
                assert_eq!(
                    r,
                    Err(EphemerisError::DataUnavailable),
                    "EMB(jd={jd}, {frame:?}) must be DataUnavailable, got {r:?}"
                );
            }
        }
    }

    // ==================================================================
    // 6. supported_range
    // ==================================================================

    /// supported_range は正順序（start < end）かつ 2000-01-01 と 2017 検証エポックを含む。
    /// 殺す変異: start/end の入れ替え・空範囲・検証エポックを外す範囲。
    #[test]
    fn supported_range_is_ordered_and_contains_validation_epochs() {
        let eph = AnalyticalEphemeris::new();
        let r = eph.supported_range();
        let start = r.start.jd2().jd();
        let end = r.end.jd2().jd();
        assert!(
            start < end,
            "supported_range not ordered: start={start}, end={end}"
        );
        // 2000-01-01（J2000）を含む。
        assert!(
            start <= JD_J2000 && JD_J2000 <= end,
            "supported_range must contain J2000 ({JD_J2000}): [{start}, {end}]"
        );
        // 2017 皆既日食エポックを含む。
        assert!(
            start <= JD_2017_ECLIPSE && JD_2017_ECLIPSE <= end,
            "supported_range must contain 2017 eclipse ({JD_2017_ECLIPSE}): [{start}, {end}]"
        );
    }

    // ==================================================================
    // 7. metadata
    // ==================================================================

    /// metadata: model に "VSOP87D" と "ELP" を含む、provenance（source/license/version）非空、
    ///   supported == supported_range。max_residual_arcsec は NaN 可（未測定）。
    /// 殺す変異: model 文字列の取り違え・provenance 空文字・supported を別範囲にする。
    #[test]
    fn metadata_has_model_provenance_and_matching_range() {
        let eph = AnalyticalEphemeris::new();
        let meta = eph.metadata();
        assert!(
            meta.model.contains("VSOP87D"),
            "metadata.model must mention VSOP87D: {:?}",
            meta.model
        );
        assert!(
            meta.model.contains("ELP"),
            "metadata.model must mention ELP: {:?}",
            meta.model
        );
        assert!(!meta.source.is_empty(), "metadata.source must be non-empty");
        assert!(
            !meta.license.is_empty(),
            "metadata.license must be non-empty"
        );
        assert!(
            !meta.version.is_empty(),
            "metadata.version must be non-empty"
        );
        assert_eq!(
            meta.supported,
            eph.supported_range(),
            "metadata.supported must equal supported_range()"
        );
    }

    // ==================================================================
    // 8. Send + Sync コンパイル時アサーション
    // ==================================================================

    /// AnalyticalEphemeris: Send + Sync（Ephemeris trait の境界）。
    #[test]
    fn is_send_and_sync() {
        fn assert_send_sync<T: Send + Sync>() {}
        assert_send_sync::<AnalyticalEphemeris>();
    }
}
