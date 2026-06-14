# 性能・運用設計 (operations)

精度設計（accuracy/numerical-policy/algorithms/physical-models）と別に、**現実に運用できるか**を決める。
対象: (D1) 性能設計、(D2) データ版管理と再現性／fingerprint、(D4) tz/CLI、(D5) MSRV/no_std/feature。
（D3「事前計算→DB連携」は**利用側アプリの責務**でライブラリ範囲外のため削除。ライブラリが負う要件は D2 に集約: 決定性の結果＋fingerprint＋serde を提供するのみ。事前計算・キャッシュ・DB・通知スケジューラは利用側で実装する。番号 D4/D5 は既存参照のため据置。）

> 状態: ドラフト（Milestone 0）。性能目標は plan §19 の暫定値。**性能より正確性を優先**（plan §19）。
> **確定（B2）**: Standard 局地の既定は **直接評価器 InstantaneousEvaluator（037, fit 誤差ゼロ）**。多項式（022）は経路/GeoJSON/NASA エクスポート用＋「M2 で fit 残差 <0.01″ を実測実証後のバッチ局地の任意最適化（未達はフォールバック）」。本セッションで一時的に置いた「局地＝多項式既定」（旧 D1）は**撤回**し、architecture.md §6.1 の直接既定へ整合させる。

---

## D1. 性能設計

### 性能目標（plan §19 暫定）
- 局地条件 1 地点: Standard で **<1 ms** 目標
- 100 年間の日食検索: **数秒以内**
- 1 万地点の局地条件、全球経路 1 件、太陽/月位置 1 回、ベッセル生成

### コストの所在
支配項は **天体暦評価**（VSOP87D 数千項 + ELP/MPP02 フル数万項の Σ A·cos）。1 回 ~10–100 µs。

- **直接評価器（037）のコスト**: 局地は C1–C4＋最大の Brent 求根（各 ~10–30 反復）で、反復ごとに瞬時ベッセル要素＝太陽・月の見かけ位置を再評価。1 地点 ≈ 5 接触 × 20 反復 × ~50 µs ≈ **5 ms**。Standard 局地は接触 5 点＋反復のみのため、単一〜少数地点では許容（精度最優先＝fit 誤差ゼロ, B2）。多地点バッチで律速になる場合に限り後述の最適化を検討する。

### 設計判断（性能と精度の両立）

1. **Standard 局地の既定 = 直接評価器（037, fit 誤差ゼロ, B2）**。精度最優先のため各地点で暦を直接再評価する（architecture §6.1）。
   - **多項式（022）の局地転用は「実証ゲート後の任意最適化」に限定**: 日食ごとに瞬時要素（021）を直接サンプリングし Chebyshev 最小二乗で多項式化（022 / numerical-policy §A4）。**M2 で fit 残差 <0.01″（目標の 1/10）を実測実証してはじめて**、バッチ局地（多地点）の任意最適化として採用可。**未達は直接（037）へフォールバック**。実証前は局地で多項式を既定にしない。
   - 多項式は経路/GeoJSON/NASA エクスポートでは無条件に本番経路（多点が µs オーダーの Horner 評価で済む）。
   - L7 の「直接 vs 多項式」残差比較で多項式の妥当性を裏取りし、局地転用の採否を最終決定する。
2. **探索性能は Fast プロファイル棄却＋並列で追求**（局地多項式化に依存しない）。
   - 100 年 ≈ 1240 朔。各朔は Fast 暦で合 solver＋早期棄却（偽陰性ゼロのマージン: physical/accuracy）。実日食（~240 件）のみフル Standard ベッセル。100 年検索が数秒に収まる設計。
3. **項数階層**: `AnalyticalEphemeris` は AccuracyProfile で級数打切りを切替（Fast=少項 / Standard=フル / Reference=DE）。
4. **並列化**: 朔の探索・地点ループは独立 → **rayon を `parallel` feature で**。`Ephemeris: Send + Sync`（既定 trait 制約）で並列安全。core は依存を増やさないため parallel は上位 crate のオプション。多地点バッチの直接評価コストは rayon 並列で吸収するのを第一手段とする（局地多項式化は実証ゲート後の任意最適化）。
5. **キャッシュ**: 同一日食内で太陽・月の見かけ位置をサンプル時刻でメモ化（多項式フィット用サンプルの再利用）。実行時グローバルキャッシュは持たない（決定性・スレッド安全のため）。

### ベンチマーク（Criterion, ISSUE 追加候補）
plan §19 の各項を Criterion で計測し、回帰検出を CI（または nightly）に置く。目標は暫定、正確性を優先。

---

## D2. データ版管理と再現性 / fingerprint

結果は **ライブラリ版** と **データ版**（係数・章動・EOP・閏秒・ΔT）の両方に依存する。これを明示管理する。

### データの2分類
| 区分 | 例 | 更新頻度 | 管理 |
|---|---|---|---|
| 係数（不変寄り） | VSOP87/ELP/MPP02/章動係数 | ライブラリ改訂時のみ | crate 同梱、ephemeris_version に紐付け |
| 時変データ | EOP（UT1/極運動）、閏秒、ΔT 予測 | 週〜半年・随時 | **versioned スナップショット同梱 ＋ from_path 上書き** |

### 再現性の保証
- **同一（ライブラリ版 + データ版 + config）→ ビット同一結果**（決定性: 実行時ネットワーク禁止・Date::now を計算に使わない・スレッド順序非依存）。回帰試験と DB の前提。
- `CalculationMetadata.fingerprint()` = ハッシュ(library_version, ephemeris_model+version, ΔT モデル+データ版, EOP データ版, earth_model, lunar/solar_radius_model, 主要 config)。

### 更新と SemVer
- **係数/打切りの変更** → 新 ephemeris_version → fingerprint 変化 → 出力が動く可能性。**MINOR 以上**で CHANGELOG に期待差分を明記。許容超で動くなら移行ノート。
- **EOP/閏秒スナップショット更新** → 出力が（特に近傍/将来日付で）動く。**MINOR 扱い**＋期待差分の桁を明記。旧データに固定したい利用者向けに `from_path` で版ピン可能。
- **将来日付の結果は暫定**: ΔT/EOP 予測更新で将来日食 UTC は変わる（accuracy §0(b)）。fingerprint に ΔT/EOP データ版を含めるので**陳腐化を検出可能**。「将来日付は不確実性帯つき・EOP 更新時に再計算推奨」を README に明記。

### データ期限
- 閏秒・EOP は `valid_to` を持つ。超過時は `Missing*Data` を返す（api-draft A2）。同梱データの鮮度を CI で警告（ISSUE-046）。

---

## D4. タイムゾーン / CLI

- **core ライブラリは tz 非依存**。全ての時刻 I/O は UTC/TT。tz 変換は持たない。
- **CLI（umbra-cli）のみ tz 依存**: `umbra local --timezone Asia/Tokyo` 等の表示・入力変換に **chrono-tz（または jiff）** を使用。
  - ライセンス: chrono-tz = MIT/Apache-2.0（許諾型・OK）。同梱 IANA tz データはパブリックドメイン。jiff も Apache-2.0/MIT。**cargo-deny allow-list に追加**（ISSUE-046 / data-sources §6）。
  - 西経・各種入力の吸収は CLI 層（conventions §3、ISSUE-032）。

---

## D5. MSRV / no_std / feature 構成

### crate 別 no_std 方針
| crate | 方針 |
|---|---|
| umbra-core | **no_std 互換・純粋**（trig は `libm`）。純数学・時刻演算・solver。**型 `TimeData`（データ束ね）/`TimeScales`（変換 facade, TimeData から構築, 変換は Result）/`DataSetMetadata` は umbra-core の純粋型**（B3）。**同梱バイトは core に置かない**。 |
| umbra-ephemeris | 既定 std。no_std は将来目標。**`TimeData::bundled()`（同梱バイト）は本 crate に置き re-export**（`bundled-data` feature ゲート, B3）。同梱データは include_bytes で可能だが v1.0 で no_std 保証はしない。 |
| umbra-eclipse | 既定 std。 |
| umbra-geo / cli | std 前提。 |

> **型と同梱データの crate 帰属（確定 B3）**: `TimeData`/`TimeScales`/`DataSetMetadata` 型定義は umbra-core（純粋・no_std 互換, データを持たない）。`bundled()`（同梱バイトを返すコンストラクタ）は **umbra-ephemeris に置き `bundled-data` feature でゲート**し、上位 crate へ re-export する（core にデータを置かない原則を守る）。off 時の挙動は下記。

### feature マトリクス
| feature | 既定 | 効果 |
|---|---|---|
| `std` | on | 標準ライブラリ。off で no_std（core） |
| `libm` | (no_std時) | no_std の三角関数 |
| `serde` | off | 直列化（tag="type"、単位はフィールド名、api-draft A7） |
| `jpl` | off | Reference DE バックエンド（同梱しない、data-sources §2.3） |
| `geojson` | off | umbra-geo 出力（M9） |
| `parallel` | off | rayon 並列（探索・地点ループ） |
| `bundled-data` | on | EOP/閏秒/係数スナップショット同梱（umbra-ephemeris, B3）。**off 時は `bundled()` を `#[cfg(feature = "bundled-data")]` で消す**（API から不在＝コンパイル時に存在しない）→ 利用者は `TimeData::from_path()` 必須。on/off で同一データなら `bundled()` と `from_path()` は同一結果（ISSUE-042 同等性テスト） |

> **bundled-data ゲート挙動（B3, 要追記事項）**:
> - **on**: `umbra-ephemeris` が `include_bytes!` で EOP C04 / 閏秒 / ΔT（および係数）スナップショットを埋込み、`TimeData::bundled()` を提供・re-export。
> - **off**: `bundled()` は `#[cfg(feature = "bundled-data")]` で**シンボルごと消える**（実行時エラーではなくコンパイル時に不在）。利用者は `TimeData::from_path(dir)` で外部ファイルを供給する。
> - **同梱データサイズ概算（要追記・要実測で確定）**: EOP C04（1962–現在, 日次 ~2.3 万行 × 数値数列）≈ 数 MB オーダー（packed で圧縮余地あり）、閏秒テーブル ≈ 数 KB、ΔT 多項式係数 ≈ 数 KB。係数スナップショット（VSOP87D/ELP-MPP02/章動）は別ゲート（M2 で項数確定後にサイズ確定）。**実バイト数は M2 の packed 生成後に DataSetMetadata.checksum とともに確定（要確認）**。バイナリ肥大を避けたい利用者向けに off + from_path を用意する理由。

### MSRV・依存方針
- **MSRV**: 保守的に固定（例: Rust 1.75、CI で検証。確定は要決定）。引上げは MINOR ＋ CHANGELOG 明記。
- **依存最小化**: core は `libm`（任意）・`thiserror` 程度。重い依存を避ける。**cargo-deny でライセンス allow-list（MIT/Apache/PD）と cargo-audit を CI 化**（ISSUE-046）。

---

## 反映先（本書が改訂/制約するもの）

- **architecture.md §6.1 / api-draft.md**: 「**Standard 局地＝直接(037, fit 誤差ゼロ)**」を**維持**（B2）。旧 D1 の「局地＝高精度フィット多項式既定」改訂は**撤回**。多項式(022)はエクスポート＋「M2 で fit 残差 <0.01″ を実証後のバッチ局地の任意最適化（未達はフォールバック）」と位置づける。`EngineConfig` に残差ゲート閾値の場は実証後に検討。
- **ISSUE-022**: エクスポート（経路/GeoJSON/NASA）用＋M2 実証後の局地任意最適化（残差 <0.01″ ゲート, 未達はフォールバック）。**ISSUE-037**: Standard 局地の既定供給源（直接, fit 誤差ゼロ）／多項式ゲート未達時のフォールバック。
- **ISSUE-031/032（CLI）**: chrono-tz/jiff 依存（D4）。**ISSUE-036**: jpl feature（D5）。**ISSUE-046**: cargo-deny allow-list に chrono-tz 等、bundled-data 鮮度 CI。
- **ISSUE-042（B3）**: `TimeData`/`TimeScales`/`DataSetMetadata` は umbra-core 純粋型、`bundled()` は umbra-ephemeris で `bundled-data` ゲート＋re-export（core にデータを置かない）。
- **ISSUE 追加候補**: Criterion ベンチ（D1）、`parallel`/`bundled-data` feature 整備、MSRV/CI マトリクス。
- **利用側連携（ライブラリ範囲外）**: 事前計算・キャッシュ・DB・通知は利用側アプリの責務。core は「決定性結果＋fingerprint＋serde」のみ提供（D2）。
