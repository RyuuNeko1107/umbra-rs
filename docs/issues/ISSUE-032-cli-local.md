# ISSUE-032: CLI local（umbra local --date --lat --lon --elevation --timezone・西経入力吸収）

- crate: umbra-cli
- 依存: ISSUE-024〜028（局地条件）, ISSUE-031（CLI 基盤・引数共通化）, ISSUE-002（角度/緯度経度 newtype・西経吸収）, ISSUE-007, ISSUE-001
- モード(tdd-workflow): standard（CLI ラッパ。局地計算はコアが担保。西経/タイムゾーン/標高の入力解釈と UTC+TT/可視性表示の正しさが要点。数式は持たないため standard）

## 目的
指定地点・指定日の局地条件を表示する CLI サブコマンド `umbra local` を実装する（architecture §7, api-draft §3.2 `local_circumstances`/`next_visible_eclipse`）。
- 引数: `--date <date> --lat <deg> --lon <deg> --elevation <m> --timezone <tz>`。
- **西経入力を吸収**（`--lon` 負＝西経 → 東経正へ・conventions §3, `EastLongitude::from_signed_degrees`・api-draft §1.2）。
- 出力: C1/C2/最大/C3/C4（UTC+TT、ローカル時刻＝`--timezone`）・食分・食面積・最大高度方位・**可視性 Visibility**（api-draft §3.4）。

## 非目的
- search の列挙（ISSUE-031）。`--date` 周辺の該当日食を `next_visible_eclipse`/`local_circumstances` で扱う。
- 局地計算本体（ISSUE-024〜028）。本 issue は CLI 境界＋入力正規化＋出力整形。
- 経路/GeoJSON（別 issue）。

## 公開インターフェース
api-draft §1.2/§3.2/§3.4、conventions §3/§7、CLI 仕様:
```
umbra local --date <DATE> --lat <DEG> --lon <DEG> [--elevation <M>]
            [--timezone <TZ>] [--accuracy <...>] [--format <text|json>]
            [--refraction <none|standard>]
```
- `--lat`: 測地緯度（`GeodeticLatitude::from_degrees`・api-draft §1.2）。
- `--lon`: 東経正、**負＝西経も受理**（`EastLongitude::from_signed_degrees`・conventions §3, api-draft §1.2）。
- `--elevation`: 楕円体高 m（既定 0、`Meters`・conventions §4）。
- `--timezone`: ローカル表示用 TZ（IANA 名 or UTC オフセット。出力のローカル時刻表示用。内部計算は UTC/TT・conventions §6）。
- `--refraction`: `RefractionModel`（既定 standard? 要確認。conventions §7 既定は幾何高度だが通知用途で補正併記）。
```rust
pub struct LocalArgs { pub date: String, pub lat: f64, pub lon: f64, pub elevation: f64, pub timezone: Option<String>, pub accuracy: AccuracyArg, pub format: FormatArg, pub refraction: RefractionArg }
pub fn run_local(args: LocalArgs) -> Result<(), CliError>;
```
- 出力: `LocalCircumstances`（api-draft §3.4）＋ `CalculationMetadata`。可視性 6 値（api-draft §3.4）。

## 数式・アルゴリズムの出典
- 本 issue は計算を持たない（CLI 境界）。局地条件は `local_circumstances`（ISSUE-024〜028）。
- 西経吸収: conventions §3（東経正・西経は境界で吸収）。`from_signed_degrees`（api-draft §1.2）。
- タイムゾーン変換: ローカル時刻表示は IANA tz / オフセット（内部は UTC・conventions §6）。要確認: tz データベース依存（chrono-tz 等）の採否とライセンス（data-sources §6）。UTC オフセットのみ対応にする簡素案も可。
- 可視性・高度方位・接触の出典はコア側（ISSUE-025/026/027/028、Explanatory Supplement §11 / Meeus Ch.54 / 球面天文標準）。

## 単位 / 時刻系 / 座標系
- 入力: 緯度=度（測地）、経度=度（東経正・西経吸収）、標高=m（楕円体高）。conventions §3/§4。
- 出力時刻: 接触・最大は **UTC+TT 両方**（accuracy.md §0）＋ローカル時刻（`--timezone`）。将来日食は ΔT 不確実性帯併記（accuracy.md §0/§2.3）。
- 角度: 高度・方位（**北0東回り**・conventions §7）、位置角（天の北0東回り）、食分/食面積 無次元。
- 座標系: 観測者 ITRS 経由（conventions §5、ISSUE-011）。

## アルゴリズム概要
1. 引数パース。`--lon` を `from_signed_degrees` で東経正へ（西経吸収・conventions §3）。`Observer` 構築（測地緯度・東経・楕円体高）。
2. `--date` 周辺の該当日食を取得（`next_visible_eclipse` or 指定日の search→`local_circumstances`）。
3. 局地条件計算（ISSUE-024〜028、`--accuracy`/`--refraction` 反映）。
4. 出力: C1〜C4（UTC+TT+ローカル）・食分・食面積・最大高度方位（北0東回り）・可視性 6 値・metadata。部分食地点は c2/c3 を "—"（None）表示。
5. 該当日食なし → 「指定日に当地点で食なし」を明示（エラーにしない・api-draft §3.2）。
- 境界/堅牢性: 緯度 ±90°/経度 ±180° 範囲外 → `DomainError`。西経/東経の取り違え防止（負値=西経を明示ヘルプ）。可視域外 → `Visibility::NotVisible`/`BelowHorizon` を表示。日の出日没中 → `SunriseEclipse`/`SunsetEclipse`（ISSUE-028）。

## 受け入れテスト
テストレベル **L6（局地）＋ CLI 統合**。基準値は実装へコピー禁止（ISSUE-029 フィクスチャ・出典明記）:
- 既知地点・既知日食（ゴールデン20・ISSUE-029）: C1〜C4・最大・食分・最大高度方位・可視性が NASA/USNO と整合（慣習を揃えて・accuracy.md §3.1）。
- **西経吸収**: `--lon -100`（西経100°）が東経260°相当として計算され、`--lon 260` と一致（conventions §3）。
- タイムゾーン: `--timezone` でローカル時刻が UTC から正しくオフセット（UTC/TT は不変）。
- 標高差: `--elevation 0` と `4000` で接触時刻が標高分わずかに変化（ISSUE-024 の ζ 補正反映）。
- 可視性 6 値: 中心線上（FullyVisible）/日の出中（SunriseEclipse）/日没中（SunsetEclipse）/一部地平下（PartialVisible）/常時地平下（BelowHorizon）/食域外（NotVisible）を各フィクスチャで（ISSUE-028 と整合）。
- 部分食地点: c2/c3 が "—"（None）表示。
- 出力に UTC+TT 両方・ΔT 不確実性帯・metadata（accuracy.md §0）。
- 異常系: 緯度/経度範囲外、不正日付、対応年代外 → 非0終了＋明確メッセージ。

## 許容誤差
CLI は計算許容を新設しない。値の許容はコア＝accuracy.md §2:
- 接触 **±2s**（TT 基準・幾何）、最大食 **±1〜2s**、食分 **±0.0005**、高度方位は表示精度（ISSUE-028）。
- **UTC 絶対は ΔT/UT1 律速**（accuracy.md §0/§2.3）。将来日食は不確実性帯を出力。
- CLI 要件: 西経↔東経・標高・タイムゾーン変換で**値を歪めない**（境界変換の正しさ）。有効桁を metadata 精度に合わせ捏造しない。
- 根拠: 表示層としてコア精度を素通しし、UTC/TT 分離・可視性・西経吸収を正しく提示（conventions §3/§6/§7、accuracy.md §0）。

## 実装メモ
- 西経吸収は `EastLongitude::from_signed_degrees`（api-draft §1.2、conventions §3）。ヘルプに「負=西経」を明記。
- タイムゾーンは表示用のみ（内部 UTC/TT 不変・conventions §6）。tz DB 依存（chrono-tz 等）の採否・ライセンスは要確認（data-sources §6）。最小実装は UTC オフセット指定のみでも可。
- `--refraction` 既定の選択（幾何高度 conventions §7 既定 vs 通知用途 standard 併記）を要確認・固定。補正前後を出せると親切（conventions §7）。
- 部分食地点 c2/c3=None を "—" 表示（api-draft §6 None 設計）。可視性 6 値を必ず提示。
- 実行時ネットワーク禁止（data-sources §0）。同梱データ前提（ISSUE-031 と共通化）。
- レビュー重点: 西経吸収、タイムゾーン表示の内部不変、標高反映、可視性 6 値、UTC+TT+ローカル、c2/c3=None 表示、範囲外エラー。
