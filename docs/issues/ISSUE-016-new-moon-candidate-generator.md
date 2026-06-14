# ISSUE-016: New moon candidate generator（朔の概算→検索窓付与・期間→新月単位分解）

- crate: umbra-eclipse
- 依存: ISSUE-006（UTC/TAI/TT）, ISSUE-007（ΔT/UT1/EOP, 時刻変換）, ISSUE-012（Ephemeris trait, 概算位置に Mock/Analytical 利用可）, umbra-core（TtInstant/UtcInstant, TimeRange/TimeInterval）
- モード(tdd-workflow): standard（内部の候補生成層。出力は精解 solver（ISSUE-017）への入力窓であり、最終精度を直接律速しない。ただし「偽陰性なし＝朔を1件も取りこぼさない」ことが search 全体の正しさを担保するため、網羅性テストは厳密に行う。standard）

## 目的
期間 `UtcRange` を「新月（朔）」単位に分解し、各朔に対して合 solver（ISSUE-017）が確実にブラケットを取れる**検索窓 `TimeInterval<TtInstant>`** を付与した候補列を生成する（architecture §3 データフロー先頭）。
- 朔の概算時刻を朔望月周期と簡易月相から算出し、その近傍に十分なマージンの窓を張る。
- 偽陽性は可（後段で棄却）、**偽陰性は不可**（朔を1件も落とさない。accuracy.md Fast「日食の見落とし: なし」/ architecture §3）。

## 非目的
- 朔の精密求解（地心合のゼロ点）= ISSUE-017。本 Issue は概算と窓付与のみ。
- 日食可能性の早期棄却（月黄緯・角距離フィルタ）= ISSUE-018。
- 高精度な月相計算。概算は VSOP87D/ELP の簡易評価 or 平均朔望式で足りる（窓マージンで吸収）。
- 時刻系変換の実装（ISSUE-006/007 を利用するクライアント）。

## 公開インターフェース
本 Issue は主に crate 内部 API（`pub(crate)`）。検証のため一部 `pub` 化候補。api-draft に直接の型はないが `TimeInterval<TtInstant>`（api-draft §1.3）を窓に使う。

```rust
/// 1 朔の候補。合 solver への入力窓を含む。
pub(crate) struct NewMoonCandidate {
    pub approx_tt: TtInstant,                 // 朔の概算時刻（地心合の概算）
    pub search_window: TimeInterval<TtInstant>, // ISSUE-017 がブラケットを取る窓
    pub lunation_number: i64,                 // Brown lunation 等の安定番号（event_key 由来候補）
}

/// 期間内の全朔候補を時系列順で返す（偽陰性なし）。
pub(crate) fn new_moon_candidates(
    eph: &impl Ephemeris,
    time_scales: &TimeScales,
    range: TimeRange<UtcInstant>,
) -> Result<Vec<NewMoonCandidate>, EclipseError>;
```

- 窓は TT 基準（conventions §6, ベッセル/位置計算は TT）。期間境界は UTC 入力→TT へ境界変換。
- `lunation_number` は `SolarEclipse.event_key`（api-draft §3.4）生成の素材候補（採番規則はレビューで確定）。

## 数式・アルゴリズムの出典
- **平均朔望月（synodic month）= 29.530588 日**（平均値）。朔の概算採番に使用。出典: Meeus, *Astronomical Algorithms* (2nd ed.), Ch.49「Phases of the Moon」式 49.1（JDE of mean phase, k からの多項式）。
- **平均新月時刻**: Meeus Ch.49 の `k`（lunation index, 新月で整数）→ `JDE = 2451550.09766 + 29.530588861·k + 補正多項式`（式 49.1）。初期窓中心はこの平均式で十分（精解は ISSUE-017）。Meeus Ch.49 の周期補正項（太陽・月の平均近点角等）は窓を狭めたい場合のみ採用、必須でない。
- 代替（暦直接）: 太陽黄経 λ_sun と月黄経 λ_moon の差 `Δλ = λ_moon − λ_sun` を概算評価し、`Δλ` が 0（mod 2π）となる近傍を平均周期で刻む。どちらでも可。**採用方式と窓マージンの根拠を実装コメントに残す**（conventions §10）。
- 朔望月の変動（実際の朔は平均から ±約14時間ずれる, Meeus Ch.49 解説）を窓マージンの下限根拠とする。

## 単位 / 時刻系 / 座標系
- 時刻系: 窓・概算時刻ともに **TtInstant**（conventions §6）。期間入力は UTC、境界で TT へ変換（ISSUE-006/007）。
- 角度: 月相 `Δλ` はラジアン、循環量として `[0, 2π)`/連続化は ISSUE-017 側の責務（本 Issue は概算のみ）。
- 座標系: 黄道座標（of date 近似で可。概算のため frame 厳密性不要）。

## アルゴリズム概要
1. 期間 [start, end]（UTC）を TT へ境界変換。
2. start 直前の lunation index `k_0` を Meeus 式 49.1 から逆算（`k ≈ (year − 2000)·12.3685`、Ch.49）。end を超えるまで `k` を 1 ずつ増やし平均朔 JDE を生成。
3. 各平均朔に対し概算時刻 `approx_tt` と窓 `[approx_tt − Δ, approx_tt + Δ]` を付与。**Δ は平均朔の最大ずれ（≈±0.6 日）に安全係数を掛けた値**（例 ±1 日）を既定とし、定数の根拠をコメント化（magic number 禁止, conventions §11）。
4. 期間端の朔（窓が範囲外へはみ出す朔）も**取りこぼさない**ため、start−1朔／end+1朔まで生成してから窓が範囲と交差するものを残す。
5. 時系列順 `Vec` で返す。
- 数値安定性: 平均式のみで折返し問題は発生しない。窓マージンが「偽陰性なし」を担保する唯一の砦なので、マージンはテストで下限検証（下記）。

## 受け入れテスト
accuracy.md テストレベル **L5（全球日食）前段の網羅性**／補助的に L2（時刻）。
- **偽陰性なしの網羅テスト（最重要）**: 既知の日食日（NASA 5千年カタログの朔日, data-sources §4.1, 第二義オラクル）が、生成された候補窓のいずれかに**必ず内包**されることを 1900–2100 で確認。1件でも窓外なら fail。基準値は fixtures から取得（実装にハードコードしない, conventions §11）。
- 朔の総数チェック: 100 年間の朔の個数が `≈ 100·12.3685 ≈ 1237 件`（朔望月数）と一致することを確認（off-by-one を検出）。
- 窓マージン下限テスト: 平均朔と精解朔（ISSUE-017 or DE 由来の真の朔）の差が常に窓半幅 Δ 未満であることを多数サンプルで確認（Δ が不足＝偽陰性リスクの直接検出）。
- 期間端テスト: 範囲開始/終了直後の朔、範囲が1朔未満、空範囲、を MockEphemeris/平均式で検証。
- `lunation_number` の単調増加・一意性。

## 許容誤差
- 本 Issue は精度バジェット（accuracy.md §2.1）に**直接寄与しない**（後段 solver が精度を担保）。担保すべきは**網羅性**: 窓が真の朔を100%内包すること。
- 概算時刻の許容: 真の朔に対し窓半幅 Δ 未満であれば良い（精度要求なし）。Δ の既定値は「平均朔の最大ずれ + 安全マージン」で設定し、テスト（窓マージン下限テスト）で 0 件の窓外を保証。
- 窓を不必要に広げると後段 solver の粗走査コストが増える。広すぎ/狭すぎのトレードオフは性能テストで調整（精度ではなくコスト指標）。

## 実装メモ
- 窓半幅 Δ は `EngineConfig` 由来にはせず本層の定数（根拠コメント付き）。将来 Fast/Reference で刻みを変えるなら設定化（要レビュー）。
- 概算に使う暦は精度不問のため、Fast 相当の簡易評価で良い（Standard 暦をフルに回すのは無駄）。ただしバックエンド抽象（`impl Ephemeris`）は維持し、MockEphemeris でもテスト可能にする。
- `lunation_number` は event_key（api-draft §3.4, plan §22 DB キー）の安定素材。採番起点（Meeus k の 0 = 2000-01-06 新月）を固定しコメント化。
- レビュー重点: 偽陰性ゼロ（窓マージンの数値根拠）、期間端の取りこぼし、朔個数の off-by-one。
