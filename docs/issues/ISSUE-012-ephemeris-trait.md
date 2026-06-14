# ISSUE-012: Ephemeris trait（太陽・月統合の暦バックエンド抽象）

- crate: umbra-ephemeris
- 依存: ISSUE-001相当（umbra-core 基盤型: Vector3 / TdbInstant / TimeRange / Radians）※要確認（core側の起票番号は別担当）
- モード(tdd-workflow): strict（公開仕様＝trait署名・StateVector・metadata 契約はSemVer境界。api-draft §2 を一字一句固定するため strict）

## 目的
太陽・月（および Earth / EarthMoonBarycenter）の状態ベクトルを単一 trait で供給する暦バックエンド抽象を定義する。解析暦（VSOP87D+ELP/MPP02）・JPL DE・Mock を同一インターフェースで差し替え可能にし、暦由来誤差と幾何由来誤差を分離する基盤（accuracy.md §4 層分解）を提供する。

- `Ephemeris` trait（`state` / `supported_range` / `metadata`）の確定。
- `Body` / `Origin` / `EphemerisFrame` / `StateVector` / `EphemerisMetadata` / `EphemerisError` の型定義。
- 速度供給方針（解析微分 / 対称差分 / 補間微分のいずれかをバックエンドが選択）の契約化（StateVector.velocity = Option）。

## 非目的
- VSOP87D / ELP/MPP02 の数式実装（ISSUE-013/014）。
- 見かけ位置補正（光行時間・歳差章動・光行差。ISSUE-015）。
- 係数生成パイプライン（ISSUE-033/034）。
- JPL バックエンド実体・SPK reader（ISSUE-036、Milestone 10）。
- 時刻系変換（TimeScales。別担当）。本 trait は `TdbInstant` を受けるのみ。

## 公開インターフェース
api-draft §2 に準拠（一字一句これを正本とする）:

```rust
#[non_exhaustive] #[derive(Clone, Copy, Debug, PartialEq, Eq)] pub enum Body { Sun, Earth, Moon, EarthMoonBarycenter }
#[non_exhaustive] #[derive(Clone, Copy, Debug)] pub enum Origin { SolarSystemBarycenter, Geocenter }
#[non_exhaustive] #[derive(Clone, Copy, Debug)] pub enum EphemerisFrame { Icrs, EclipticOfDate }

#[derive(Clone, Copy, Debug)] pub struct StateVector { pub position: Vector3, pub velocity: Option<Vector3> }

#[derive(Clone, Debug)] pub struct EphemerisMetadata {
    pub model: String,        // 例 "VSOP87D+ELP/MPP02"
    pub version: String,      // 採用打切り次数・達成残差を含む識別子
    pub source: String, pub license: String,
    pub supported: TimeRange<TdbInstant>,
    pub max_residual_arcsec: f64,   // accuracy.md §2.4 実測値
}

pub trait Ephemeris: Send + Sync {
    fn state(&self, body: Body, time: TdbInstant, origin: Origin, frame: EphemerisFrame)
        -> Result<StateVector, EphemerisError>;
    fn supported_range(&self) -> TimeRange<TdbInstant>;
    fn metadata(&self) -> EphemerisMetadata;
}

#[non_exhaustive] #[derive(Debug)] pub enum EphemerisError {
    OutOfSupportedRange, DataUnavailable, Io(/* ... */),
}
```

- `Position<F>` ではなく素の `Vector3` を返す点に注意（フレームは `EphemerisFrame` 引数で表現。型レベル区別は見かけ位置層 ISSUE-015 で `Position<Gcrs>` 等へ持ち上げる）。
- 単位は km（conventions §1: 天体暦内部 AU 許容だが trait 境界で km へ変換して返す）。これを doc コメントとテストで固定する。

## 数式・アルゴリズムの出典
- 本 Issue は抽象定義のため数式なし。ただし契約（座標系・単位・原点の意味）の出典:
  - ICRS / Geocenter / Barycenter の定義: IERS Conventions 2010, ch.2 / IAU 1997, 2006 resolutions。
  - VSOP87D の「黄道・平均分点 of date・球面」: Bretagnon & Francou (1988) A&A 202, 309（`EphemerisFrame::EclipticOfDate` の意味づけ。data-sources §2.1）。
  - 太陽地心位置＝地球日心位置の反転: ISSUE-013 で実装（本 trait は `Body::Sun, Origin::Geocenter` の呼び出しでそれを保証する契約のみ）。

## 単位 / 時刻系 / 座標系
- 単位: position = km, velocity = km/s（Option）。AU は内部のみ、trait 境界で km。
- 時刻系: 入力 `TdbInstant`（conventions §6: Reference/暦評価は TDB。VSOP87D/ELP は実用上 TT≈TDB の差〈≲2ms〉を許容し、metadata に「TDB 引数を TT 同一視するか」を記録 ※要確認）。
- 座標系: `EphemerisFrame::Icrs`（ICRS 軸）/ `EclipticOfDate`（黄道・平均分点 of date）。`Origin` は SSB または Geocenter。

## アルゴリズム概要
1. `state()` は (body, time, origin, frame) を受け、バックエンド固有の評価へディスパッチ。
2. 範囲外時刻は `OutOfSupportedRange`。`supported_range()` と整合させる。
3. velocity は供給方式をバックエンドが選択（解析モデルは解析微分 or 対称差分、DE は SPK の速度成分／補間微分）。供給しない場合は `None`。供給方式は metadata.version 文字列に含める。
4. metadata は採用打切り次数・達成残差（max_residual_arcsec, accuracy.md §2.4）・source・license を必ず埋める。

## 受け入れテスト
- L1/契約テスト: trait オブジェクト安全性（`dyn Ephemeris` 可）・`Send + Sync` 境界をコンパイルテストで保証。
- L3（天体位置・スモーク）: `MockEphemeris`（ISSUE 別）実装を1つ用意し、`state()` が単位 km の `Vector3` を返すこと、範囲外で `OutOfSupportedRange` を返すことを確認。
- 単位境界テスト: AU で内部保持するダミー実装が境界で km へ正しく変換することを 1 桁の既知値で検証。
- metadata 必須項目が空でないこと（model/version/source/license/supported/max_residual_arcsec）を property 的にガード。
- 速度 Option 契約: velocity = None と Some(対称差分) の両ケースが呼び出し側で扱えること。

## 許容誤差
- 抽象層のため数値許容誤差なし。ただし「単位は km、誤差はバックエンド側が metadata.max_residual_arcsec で申告する」契約を固定（accuracy.md §2.1 月0.40″/太陽0.20″配分は ISSUE-013/014 で担保、本層は申告経路を保証するのみ）。

## 実装メモ
- `EphemerisError::Io` の中身は `std::io::Error` を直接持たず軽量化する（feature `jpl` 無効時に std::io 依存を引かない設計を検討。※要確認）。
- `non_exhaustive` を Body/Origin/EphemerisFrame に付与（api-draft §0 前方互換）。Origin は将来 `TopoCenter` 追加余地を残す。
- velocity の供給方式名（"analytic-derivative" / "central-difference" / "spk-velocity"）は文字列定数化して metadata.version に埋め込む。差分幅は固定せずテストで決める（architecture §4/§12）。
- フレーム型（`Position<F>`）への持ち上げは ISSUE-015 で行う。本層で PhantomData を持ち込まない（trait を薄く保つ）。
