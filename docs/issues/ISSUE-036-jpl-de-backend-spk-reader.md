# ISSUE-036: JPL DE バックエンド / SPK reader（feature jpl・Reference オラクル）

- crate: umbra-ephemeris（feature `jpl`）
- 依存: ISSUE-012（Ephemeris trait）, umbra-core（Vector3/TdbInstant/TimeRange）
- モード(tdd-workflow): **strict**（永続形式の読取り＝SPK バイナリ仕様の解釈は厳密性必須。差分テストの第一義オラクルとして正しさが全層の検証根拠になる。strict）
- マイルストン（二段構成）:
  - **最小 SPK reader（太陽/月のみ）＝Milestone 2 必須**（accuracy.md §3.3 二段ゲート）。DE440 の太陽(10)・月(301)・地球(399)・EMB(3) に限定し、Chebyshev type2/3 評価のみを実装。**暦層の第一義オラクルとして ISSUE-013/014 の打切り次数を M2 から DE 実値で決める**ために必須。reader 選定（自前 vs 既存 crate）という重い問題は**フル版へ切り離す**——最小版は自前 DAF/SPK パース（太陽/月セグメントのみ）で十分。
  - **フル版 / Reference 用途＝Milestone 10 据置**（data-sources §2.3, accuracy.md §1 Reference）。全 body 対応・reader 実装選定（自前 vs 既存純Rust reader）・Reference バックエンド完成・最終ゲート（accuracy.md §3.3 M10）は M10。
- 注: 起票時の「実装は M10」方針を、**最小版 M2 必須 / フル版 M10 据置**へ更新（B4(a)）。

## 実装状況（完了・2026-06-20）

**フル版 Reference バックエンド完了**（feature `jpl`・既定 off・DE データ非同梱）。自前 DAF/SPK パーサ＋
type2 Chebyshev 評価＋`JplEphemeris`（`Ephemeris` 実装）を **3 スライス strict TDD**（テスト設計/
実装/レビューをサブエージェント分離・各スライス cargo mutants）で実装。`xtask fetch-de440s`/
`verify-de440s` で取得・SHA-256 検証（cc34384）。

- **S1** `parse_spk_segments`（DAF 構造解析）— commit `3f73bbc`
- **S2** `eval_type2`（SPK type2 Chebyshev 位置/速度）— commit `d1f91f2`
- **S3** `JplEphemeris`＋`Ephemeris`（body 差・原点・metadata）— commit `6d9b202`

**正当性ゲート（accuracy.md §3.3 M10）**: 自前 reader を **SPICE（spiceypy 8.1.2/CSPICE N0067）の
`spkgeo`** と突合し、合成 Chebyshev（厳密値）＋実 DE440s で **位置 < 10 m・速度 < 1 mm/s 一致**を確認
（Sun/EMB/Earth wrt SSB、Sun/Moon wrt 地心）。オラクル誤差 ≪ §2.1 各層配分目標を満たす。

`from_spk_path` で利用者が任意取得の `.bsp` を読む（実行時ネットワーク禁止・非同梱方針 §6 維持）。
`EphemerisFrame::EclipticOfDate` は ISSUE-035 変換経由のため本バックエンド未提供（ICRS 直接利用が主）。

## 目的
JPL DE（既定 DE440）を読む `JplEphemeris` を feature `jpl` で実装し、`Ephemeris` trait の Reference バックエンド＝**差分テストの第一義オラクル**（accuracy.md §3.1）を提供する。解析暦（VSOP87D/ELP/MPP02）と同一のベッセル/接触パイプラインを通して差分を取り、暦由来誤差と幾何由来誤差を分離する（§4 層分解）。

## 非目的
- crate へのDEデータ同梱（**しない**。利用者が任意 DL。data-sources §2.3/§6）。
- 解析暦の実装（ISSUE-013/014）。
- Standard 経路（公開）への組込み（Reference 専用。本番は AnalyticalEphemeris）。
- 全惑星対応（太陽・月・地球・EMB に限定して十分。日食用途）。

## 公開インターフェース
api-draft §2 準拠:

```rust
#[cfg(feature = "jpl")] pub struct JplEphemeris { /* SPK reader */ }
#[cfg(feature = "jpl")] impl JplEphemeris {
    pub fn from_spk_path(path: &std::path::Path) -> Result<Self, EphemerisError>;
}
#[cfg(feature = "jpl")] impl Ephemeris for JplEphemeris { /* state/supported_range/metadata */ }
```

- feature ゲートで `jpl` 無効時はコンパイル対象外（標準ビルドに DE 依存・巨大データを引かない）。
- `state()` は ISSUE-012 契約通り km / TdbInstant / Origin / EphemerisFrame。

## 数式・アルゴリズムの出典
- **JPL DE440/441**: Park, Folkner, Williams, Boggs (2021), *AJ* 161:105, DOI:10.3847/1538-3881/abd414。data-sources §2.3。
- **SPK / SPICE カーネル形式**: NAIF SPK 仕様（DAF/SPK, SPK type 2 = Chebyshev position-only / type 3 = position+velocity）。配布: `naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/`。
  - 評価 = 区間ごとの Chebyshev 多項式（位置: type 2、位置+速度: type 3）。各セグメントの中点・半区間・係数から Clenshaw/Horner で評価。
- 時刻: DE は **TDB** 基準（barycentric dynamical time）。conventions §6（Reference は TDB）。
- 月地心位置 = 月（DE では EMB 相対 or geocentric）の取り出し方は DE の body ID（301=Moon, 399=Earth, 10=Sun, 3=EMB）に従い、SSB 基準ベクトルの差で Geocenter 原点へ。
- **最小版（M2）**: 太陽/月セグメントのみを読む自前 DAF/SPK パース + Chebyshev type2/3 評価。reader 選定問題は持ち込まない（最小自前で完結）。
- **reader 実装の選定はフル版＝Milestone 10**: 自前 DAF/SPK パーサ vs 既存純Rust reader（例: spice 系 crate の評価）。ライセンス・純Rust 方針との整合を確認（data-sources §2.3 では「reader 実装の選定は Milestone 10」）。

## 単位 / 時刻系 / 座標系
- 単位: SPK は km / km/s（既に km。AU 変換不要）。
- 時刻系: 入力 TdbInstant。SPK の ET(=TDB) と直接対応。
- 座標系: SPK は ICRF（≈ICRS）/ SSB 原点。`EphemerisFrame::Icrs` を native とし、Geocenter 原点は body 差で算出。`EclipticOfDate` 要求は ISSUE-035 経由で変換（ただし Reference は ICRS 直接利用が主）。

## アルゴリズム概要
1. `from_spk_path` で DAF/SPK ファイルを mmap or 読込、セグメント・サマリレコードを解析。
2. body ID と TDB 時刻からセグメント特定、Chebyshev 係数取得。
3. Clenshaw 評価で位置（type3 は速度も）。type2 のとき速度は係数微分 or 対称差分。
4. Origin/Frame に応じ body 差・フレーム変換。
5. supported_range は DE440 = 1550–2650（v1.0 1900–2100 を包含）。範囲外で `OutOfSupportedRange`。
6. metadata に DE 版・source・license（米政府/Caltech-JPL・同梱しない方針）・supported を記録。

## 受け入れテスト
- **第一義オラクルとしての自己検証**: NAIF/JPL 公開のテスト状態ベクトル（既知 ET の太陽・月位置）と `state()` を突合（µas/m 級一致 = reader 正当性）。基準値は JPL 公開ベクトルから（ハードコードのコピー流用でなく出典明記、conventions §11）。
- **差分テスト基盤（accuracy.md §3.1/§3.3 を駆動）**: ISSUE-013/014/015 の DE 差分テスト（月0.1″/太陽0.05″ 残差確認）の**基準提供側**として機能。同一ベッセル/接触パイプライン経由の差分が層分解できることを確認（§4）。**最小版は M2 でこの基準提供側（暦層の打切り次数決定）を担う**（accuracy.md §3.3 M2）。
- **最小版受入（M2）**: 太陽/月のみの `state()` が NAIF/JPL 公開テストベクトルと µas/m 級一致すること（reader 正当性）。type2/type3 評価とセグメント境界連続性を太陽/月セグメントで確認。これを ISSUE-013/014 の M2 打切り次数決定の基準に用いる。
- **注（幾何層は別オラクル）**: 本 reader（暦層）は accuracy.md §3.1-1 の限界どおり**暦層しか検証できない**。幾何共通バグは accuracy.md §3.1-2 の幾何独立オラクル（NASA ベッセル要素表 / 独立第二実装）で別途検証する（本 Issue の範囲外）。
- type2/type3 両対応、セグメント境界跨ぎの連続性、範囲外エラーのテスト。
- feature `jpl` 無効ビルドで本コードが除外されることのコンパイルテスト。
- CI 必須にしない（巨大データ）。nightly / 手動ジョブ（accuracy.md §3.1）。

## 許容誤差
- DE 読取り自体は**オラクル**のため「許容誤差」ではなく**正当性**を要求: 公開テストベクトルと µas / sub-m 級一致（Chebyshev 評価の数値誤差のみ許容）。
- これが accuracy.md §2.1 の各層配分（月0.40/太陽0.20/フレーム0.05/光行差0.10）の**測定基準**になるため、オラクル誤差は目標より十分小さい（< 0.005″ 目安）必要がある。§3.1「目標より良いオラクルを使う」。

## 実装メモ（ライセンス・同梱方針 — data-sources §2.3/§6）
- **DE データは crate に同梱しない**。利用者が任意取得（`.bsp`）、`from_spk_path` で読む。これによりライセンスリスク（米政府/Caltech-JPL の OSS 同梱可否未確定 §6）を回避。
- 引用: Park et al. (2021), AJ 161:105。metadata.license に「JPL/Caltech・非同梱・任意DL」明記。
- reader は純Rust 方針。既存 crate 採用時は**そのライセンスを `cargo-deny` allow-list で検査**（GPL 不可、data-sources §0）。reader 選定（自前 vs 既存）は Milestone 10 で決定。
- DE441（-13200–+17191）は将来の長期拡張用。v1.0 は DE440 既定。
- 巨大データ前提のため差分テストは nightly/手動（CI 必須にしない）。
- 実行時ネットワーク禁止（accuracy.md §5）。DL は利用者手順 / xtask の明示手順に隔離。
