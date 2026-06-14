//! 天体暦バックエンドの抽象（`docs/issues/ISSUE-012`、`docs/architecture.md` §4）。
//!
//! 太陽・月を別 trait にせず統合した [`Ephemeris`] で扱う。実装は `AnalyticalEphemeris`
//! （VSOP87D+ELP/MPP02, Milestone 2）, `JplEphemeris`（Reference, feature `jpl`）,
//! [`crate::mock::MockEphemeris`]（テスト）。

use thiserror::Error;
use umbra_core::{TdbInstant, TimeRange, Vector3};

/// 天体。
#[non_exhaustive]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Body {
    /// 太陽。
    Sun,
    /// 地球。
    Earth,
    /// 月。
    Moon,
    /// 地球・月の重心。
    EarthMoonBarycenter,
}

/// 位置の原点。
#[non_exhaustive]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Origin {
    /// 太陽系重心。
    SolarSystemBarycenter,
    /// 地心。
    Geocenter,
}

/// 出力フレーム。
#[non_exhaustive]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum EphemerisFrame {
    /// ICRS（慣性赤道座標）。
    Icrs,
    /// その時点の黄道座標。
    EclipticOfDate,
}

/// 状態ベクトル（位置と任意の速度）。位置・速度は km, km/s。
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StateVector {
    /// 位置 \[km\]。
    pub position: Vector3,
    /// 速度 \[km/s\]（解析微分などを提供できないバックエンドでは `None`）。
    pub velocity: Option<Vector3>,
}

/// バックエンドのメタデータ（結果の `CalculationMetadata` に転記。accuracy.md §2.4）。
#[derive(Clone, Debug, PartialEq)]
pub struct EphemerisMetadata {
    /// モデル名（例 `"VSOP87D+ELP/MPP02"`）。
    pub model: String,
    /// 版（採用打切り次数・達成残差を含む識別子）。
    pub version: String,
    /// 出典。
    pub source: String,
    /// ライセンス。
    pub license: String,
    /// 対応する TDB 範囲。
    pub supported: TimeRange<TdbInstant>,
    /// 達成した最大残差（秒角, 実測値。テスト用は NaN/未申告）。
    pub max_residual_arcsec: f64,
}

/// 天体暦バックエンドのエラー。
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq, Eq, Error)]
pub enum EphemerisError {
    /// 要求時刻が対応範囲外。
    #[error("requested time is outside the ephemeris supported range")]
    OutOfSupportedRange,
    /// 暦データが利用不可。
    #[error("ephemeris data unavailable")]
    DataUnavailable,
}

/// 天体暦バックエンド。
pub trait Ephemeris: Send + Sync {
    /// `body` の `time`（TDB）における状態ベクトルを `origin` 基準・`frame` で返す。
    fn state(
        &self,
        body: Body,
        time: TdbInstant,
        origin: Origin,
        frame: EphemerisFrame,
    ) -> Result<StateVector, EphemerisError>;

    /// 対応する TDB 時刻範囲。
    fn supported_range(&self) -> TimeRange<TdbInstant>;

    /// メタデータ。
    fn metadata(&self) -> EphemerisMetadata;
}
