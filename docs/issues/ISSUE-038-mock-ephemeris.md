# ISSUE-038: MockEphemeris 実装（人工配置で幾何ロジックを検証する足場）

- crate: umbra-ephemeris
- 依存: ISSUE-012（`Ephemeris` trait / `Body` / `Origin` / `EphemerisFrame` / `StateVector` / `EphemerisMetadata` / `EphemerisError`）, umbra-core（`Vector3` / `TdbInstant` / `TimeRange` / `Radians`）
- milestone: M2（暦バックエンド層。DE 無しで幾何検証を回す足場のため、解析暦実体 013/014 と並走で最優先に近い位置づけ）
- モード(tdd-workflow): standard（公開型は `MockEphemeris` の4コンストラクタのみで `Ephemeris` 契約は ISSUE-012 が strict 確定済み。本 Issue は人工配置の値設計が中心でテスト足場のため standard）

## 目的
`Ephemeris` trait（ISSUE-012）の人工配置実装 `MockEphemeris` を提供し、**DE 無し・純解析オラクルで日食幾何ロジック（影円錐・ベッセル要素・分類・局地）を単体検証する足場**を作る（architecture §4, accuracy.md §3.1）。
- 4 つの人工配置コンストラクタ:
  - `central_total()` … 完全皆既（影軸が地心を貫く中心配置、l2<0、gamma≈0）。
  - `clear_annular()` … 明確な金環（本影頂点が基本面の地球側、l2>0）。
  - `clear_partial()` … 明確な部分食（影軸が地球をかすめる、|gamma| が 1 近傍だが食は起きる）。
  - `shadow_misses_earth()` … 影が地球を完全に外す（|gamma| > 1 + l1、日食なし）。
- 各配置は**解析的に既知の太陽・月の位置**を返し、下流の幾何結果（x,y,l1,l2,gamma,種別）を解析オラクルと照合できるようにする。

## 非目的
- 実天体暦（VSOP87D/ELP/MPP02。ISSUE-013/014）や DE バックエンド（ISSUE-036）の代替ではない。**精度検証ではなく幾何ロジック検証**用。
- 見かけ位置補正（light-time/aberration/歳差章動。ISSUE-015/035）。Mock は補正前の素の位置を返す契約（補正側のテストは別途、補正を識別恒等化できる単純配置で行う）。
- 時刻系変換（TimeScales。別担当）。trait 契約どおり `TdbInstant` を受けるのみ。
- 現実的な軌道運動の再現。各コンストラクタは**最大食付近の静的配置＋一定速度**で足りる（速度は対称差分でも可、velocity=Some/None 両対応）。

## 公開インターフェース
api-draft §2 を転記・具体化（`Ephemeris` 契約は ISSUE-012 正本）:
```rust
pub struct MockEphemeris { /* 人工配置: 各 body の解析的状態関数 */ }

impl MockEphemeris {
    pub fn central_total() -> Self;       // 皆既・中心（l2<0, gamma≈0）
    pub fn clear_annular() -> Self;       // 金環（l2>0）
    pub fn clear_partial() -> Self;       // 部分（食あり・中心食なし）
    pub fn shadow_misses_earth() -> Self; // 日食なし（|gamma|>1+l1）
}

impl Ephemeris for MockEphemeris {
    fn state(&self, body: Body, time: TdbInstant, origin: Origin, frame: EphemerisFrame)
        -> Result<StateVector, EphemerisError>;
    fn supported_range(&self) -> TimeRange<TdbInstant>;
    fn metadata(&self) -> EphemerisMetadata;
}
```
- 単位は km（ISSUE-012 の trait 境界契約）。`EphemerisFrame::Icrs`（ICRS 軸）を既定で返す（黄道版は v0.1 で不要なら未対応＋metadata 注記。要確認）。
- `metadata()` は model="MockEphemeris"、source/license に「テスト専用・人工配置」、`max_residual_arcsec` は N/A 扱い（幾何検証用で精度申告対象外である旨を version 文字列に明記）。

## 数式・アルゴリズムの出典
- 数式実装ではなく**解析配置の設計**。出典は幾何条件の定義のみ:
  - gamma・l1・l2 と影軸の幾何関係: Explanatory Supplement to the Astronomical Almanac (3rd ed.) Ch.11「Eclipses」, Meeus *Astronomical Algorithms* (2nd ed.) Ch.54（ISSUE-021/023 と同一正本）。
  - 各配置の判定境界: `central` ⇔ |gamma| < 1−|l2|、`annular` ⇔ l2>0（**正本 B1: l2<0=皆既 / l2>0=金環**）、`partial` ⇔ 1−l1 < |gamma| < 1+l1 かつ中心食条件不成立、`no eclipse` ⇔ |gamma| > 1+l1（ISSUE-023 の分類境界と整合させる）。
- 配置値（太陽距離≈1 AU、月距離≈平均地心距離、半径 k=conventions §9）は**マジックナンバー禁止に従い ISSUE-041 の定数モジュール由来**で構成し、配置のオフセット量のみ Mock 固有定数として出典コメント付きで定義（conventions §11）。

## 単位 / 時刻系 / 座標系
- 単位: position = km、velocity = km/s（Option。解析微分 or 対称差分。ISSUE-012 契約）。
- 時刻系: 入力 `TdbInstant`（ISSUE-012 契約。Mock は静的配置中心なので時刻依存は最小、最大食時刻付近で連続な軌跡を返す）。
- 座標系: `EphemerisFrame::Icrs`（ICRS 軸、Geocenter 原点で太陽・月の地心位置を構成）。

## アルゴリズム概要
1. 各コンストラクタで、目標とする (gamma, l2 符号) を満たす太陽・月の地心位置（と速度）を解析的に逆算して内部に保持。
2. `state()`: body=Sun/Moon/Earth に対し、時刻 t での位置を解析関数（静的位置＋一定相対速度の線形 or 単純円運動）で返す。範囲外は `OutOfSupportedRange`。
3. `central_total`: 影軸が地心を通る（gamma≈0）配置。本影が地表へ届く距離関係で l2<0。
4. `clear_annular`: 月をやや遠方に置き本影頂点を基本面の地球側へ（l2>0）。
5. `clear_partial`: 影軸を地球縁近傍へオフセット（1−l1<|gamma|<1+l1、中心食不成立）。
6. `shadow_misses_earth`: 影軸を地球半径＋半影外へオフセット（|gamma|>1+l1）。
- 数値安定性: 各配置の gamma/l2 を構築時に一度計算し、解析オラクル値として**テスト側が独立再計算**できるよう設計値を Issue/コメントに記録（実装からのコピー禁止、conventions §11）。

## 受け入れテスト
accuracy.md テストレベル **L3（天体位置・契約スモーク）＋ L4/L5 の足場**。MockEphemeris 自体は精度オラクルではなく幾何ロジックの解析オラクル供給源（accuracy.md §3.1）。
- 契約: `state()` が km 単位 `Vector3` を返し、範囲外で `OutOfSupportedRange`。`dyn Ephemeris` 可・`Send+Sync`（ISSUE-012 と同等）。
- 各配置の幾何整合（下流 Issue が利用、本 Issue ではスモーク）:
  - `central_total`: 影軸を ISSUE-019/020/021 へ通すと x≈y≈0、l2<0、gamma≈0。
  - `clear_annular`: l2>0（金環＝正符号、ISSUE-021 の正本「l2<0=皆既 / l2>0=金環」と一致）。
  - `clear_partial`: 中心食条件不成立かつ部分食成立（ISSUE-023 分類が Partial）。
  - `shadow_misses_earth`: |gamma|>1+l1 で「日食なし」（ISSUE-018/023 が棄却/非該当）。
- velocity: `None` と `Some(対称差分/解析微分)` の両方を呼び出し側が扱えること。
- 設計値の独立照合: 各配置の目標 gamma/l2 を Issue 記載の解析式でテストが再計算し、Mock 状態から組んだ要素と一致（実装値コピー禁止）。

## 許容誤差
- 本 Issue は**精度検証対象外**（幾何ロジック検証用、accuracy.md §3.1 MockEphemeris）。
- 解析配置の内部整合（構築時 gamma/l2 と state から再構成した値）は**純幾何の丸め誤差のみ**（≤ 1e-9 相当）を目標。
- accuracy.md §2 の角度バジェット（最大食 ±1.5s 等）は実暦（013/014/036）で担保し、本 Issue では適用しない（区別を実装メモに明記）。

## 実装メモ
- **DE 無しで回せること**が要件（accuracy.md §3.1, architecture §4）。CI 必須経路に組み込み、暦データ非依存の幾何回帰を成立させる（DE 差分は nightly/手動、ISSUE-029/036）。
- 各配置の (gamma, l2 符号, 種別) を**設計意図として Issue とコード doc に明記**し、テストは設計式から独立再計算（conventions §11 オラクル非コピー）。
- マジックナンバー禁止（conventions §11）: 半径・距離の基準値は ISSUE-041 定数由来、配置オフセットのみ Mock 固有定数として出典/意図コメント付き。
- フレーム/補正の扱い（Mock は補正前素位置）を doc に固定し、ISSUE-015/035 のテストが Mock を識別恒等的に使えるよう注記。
- レビュー重点: 4 配置の gamma/l2/種別が分類境界（ISSUE-023）と矛盾しないか、km 単位境界、`Send+Sync`、設計値のオラクル非コピー。
