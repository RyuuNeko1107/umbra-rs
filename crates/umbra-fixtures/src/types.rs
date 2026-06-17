//! ゴールデンフィクスチャの公開型（`docs/issues/ISSUE-029`）。
//!
//! 外部オラクル（NASA 5MCSE / USNO Solar Eclipse Calculator）の**数値事実のみ転記**した
//! 固定回帰データの型定義。計算ロジックは持たない（実装出力との一致比較は ISSUE-030）。
//!
//! 単位・座標系（ISSUE-029 §単位 / `docs/conventions.md`）:
//! - 角度は**度**で保持（転記しやすさ。内部ラジアン規約 §1 とは別に、フィクスチャは度で境界保持）。
//! - 方位は北 0°・東回り `[0,360)`（§7）。緯度は北正、経度は**東経正** `(-180,180]`（§3）。
//! - 時刻は UTC（オラクルが TT/TD を持つ場合のみ併記。持たなければ `None`・捏造しない §0）。

use umbra_core::{TtInstant, UtcInstant};
use umbra_eclipse::{SolarEclipseKind, Visibility};

/// 出典付きのゴールデン日食フィクスチャ（全球パラメータ＋地点別局地状況）。
///
/// 全球パラメータ（種別・最大食時刻・gamma・食分）は NASA 5MCSE から、地点別の局地状況は
/// USNO Solar Eclipse Calculator から転記する（`source` に両者を明記）。
#[derive(Clone, Debug, PartialEq)]
pub struct GoldenEclipse {
    /// 一意キー（例 `"2017-08-21-total"`）。
    pub event_key: String,
    /// 全球の日食種別（NASA 種別）。
    pub kind_expected: SolarEclipseKind,
    /// 最大食（greatest eclipse）の UTC 時刻。
    pub greatest_time_utc: UtcInstant,
    /// 最大食の TT(=TD) 時刻（オラクルが TD を与える場合のみ）。
    pub greatest_time_tt: Option<TtInstant>,
    /// 影軸の地心最小距離 γ（Re, 符号付き）。
    pub gamma: f64,
    /// 全球最大食での食分（無次元）。
    pub magnitude: f64,
    /// オラクルが採用した ΔT（秒。NASA 5MCSE 値）。
    pub delta_t_seconds: Option<f64>,
    /// 地点別の局地状況（5〜10 地点目標、seed では 3〜）。
    pub locations: Vec<GoldenLocation>,
    /// 出典・取得日・k/ΔT 慣習・ライセンス注記。
    pub source: OracleSource,
}

/// 1 地点のゴールデン局地状況（USNO Solar Eclipse Calculator 転記）。
#[derive(Clone, Debug, PartialEq)]
pub struct GoldenLocation {
    /// 地点名（都市名・国）。
    pub name: String,
    /// 測地緯度（度・北正）。
    pub latitude_deg: f64,
    /// 東経（度・東経正・西経は東経正へ変換し記録, §3）。範囲 `(-180,180]`。
    pub east_longitude_deg: f64,
    /// 楕円体高（m）。
    pub elevation_m: f64,
    /// 地点条件の分類（被覆メタテスト用, accuracy §3.4 / L6）。
    pub location_class: LocationClass,
    /// 第 1 接触 C1（部分食開始）。地平下・未到来なら `None`。
    pub c1: Option<GoldenContact>,
    /// 第 2 接触 C2（皆既/金環開始＝内接）。中心相が無ければ `None`。
    pub c2: Option<GoldenContact>,
    /// 最大食（常に存在）。
    pub maximum: GoldenContact,
    /// 第 3 接触 C3（皆既/金環終了）。中心相が無ければ `None`。
    pub c3: Option<GoldenContact>,
    /// 第 4 接触 C4（部分食終了）。地平下・未到来なら `None`。
    pub c4: Option<GoldenContact>,
    /// この地点の食分（無次元）。
    pub magnitude: f64,
    /// この地点の食面積比（obscuration, `0..=1`）。
    pub obscuration: f64,
    /// 最大食時の太陽高度（度。大気差込みの USNO 値）。
    pub max_altitude_deg: f64,
    /// 最大食時の太陽方位（度・北 0 東回り `[0,360)`, §7）。
    pub max_azimuth_deg: f64,
    /// 期待される可視性（`umbra_eclipse::Visibility`）。
    pub visibility_expected: Visibility,
}

/// 接触/最大食の 1 時刻（UTC、オラクルが TT を持てば併記）＋太陽高度。
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GoldenContact {
    /// 接触時刻（UTC）。
    pub time_utc: UtcInstant,
    /// 接触時刻（TT(=TD)）。USNO は UT のみのため通常 `None`。
    pub time_tt: Option<TtInstant>,
    /// この接触時の太陽高度（度）。
    pub altitude_deg: f64,
}

/// 地点条件の分類（accuracy.md §3.4 / L6・被覆メタテスト）。
///
/// 「可視域外（食域外）」は転記すべき数値事実を持たないため独立分類を置かない。観測不能側は
/// `BelowHorizon`（最大食が地平下＝負の高度）で表現する（`maximum` を捏造しない）。
#[non_exhaustive]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum LocationClass {
    /// 中心線上（または極近傍。皆既/金環の最長級）。
    Centerline,
    /// 限界線近傍（中心相が極端に短い）。
    NearLimit,
    /// 部分食域（中心相なし）。
    PartialZone,
    /// 日の出中に食が進行。
    Sunrise,
    /// 日没中に食が進行。
    Sunset,
    /// 最大食時に太陽が地平下（観測不能側の表現）。
    BelowHorizon,
    /// 高標高地点（標高補正の確認）。
    HighElevation,
}

/// オラクル出典・取得日・慣習・ライセンス注記（ハードコード防止・data-sources §4）。
#[derive(Clone, Debug, PartialEq)]
pub struct OracleSource {
    /// 出典名（例 `"NASA 5MCSE (global) + USNO Solar Eclipse Calculator (local)"`）。
    pub name: String,
    /// 取得元 URL。
    pub url: String,
    /// 取得日（`YYYY-MM-DD`）。
    pub retrieved: String,
    /// ΔT 慣習（accuracy.md §3.1。全球=NASA / 局地=USNO の差異を明記）。
    pub delta_t_convention: String,
    /// k 慣習（Espenak 2 値 等, conventions §9）。
    pub k_convention: String,
    /// ライセンス注記（data-sources §0/§4。数値事実のみ転記の旨）。
    pub license_note: String,
}
