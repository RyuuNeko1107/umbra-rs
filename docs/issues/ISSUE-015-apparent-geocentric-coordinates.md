# ISSUE-015: Apparent geocentric coordinates（光行時間→歳差章動→光行差）

- crate: umbra-ephemeris
- 依存: ISSUE-012（Ephemeris trait）, ISSUE-013（太陽）, ISSUE-014（月）, ISSUE-035（歳差章動・フレーム変換）, umbra-core（Position<F>, Vector3, TtInstant/TdbInstant）
- モード(tdd-workflow): strict（公開仕様 `AstrometryOptions` と Standard 必須ON契約は SemVer 境界。補正の有無は結果値を左右し metadata に残す＝公開仕様。strict）

## 目的
暦の幾何地心位置（ISSUE-013/014）に、見かけ位置補正（光行時間→歳差章動→光行差→任意で相対論偏向）を適用し、観測時刻における太陽・月の**見かけ地心位置**を返す。Standard プロファイルでは light_time / aberration / precession_nutation を**必須 ON**で固定する（accuracy.md §1, architecture §5）。

## 非目的
- 暦級数の評価（ISSUE-013/014）。
- 歳差章動行列そのものの実装（ISSUE-035。本 Issue は ISSUE-035 を利用するクライアント）。
- 地心→地上（topocentric）視差・観測者補正（umbra-eclipse 側。本 Issue は geocentric まで）。
- 大気差（RefractionModel。umbra-eclipse / conventions §7）。

## 公開インターフェース
api-draft §2 に準拠:

```rust
#[derive(Clone, Copy, Debug)]
pub struct AstrometryOptions {
    pub light_time: bool, pub aberration: bool,
    pub precession_nutation: bool, pub relativistic_deflection: bool,
}
impl AstrometryOptions {
    pub fn standard() -> Self;  // light_time/aberration/precession_nutation = true, deflection = false
    // fast() は公開 API から削除。低次・一部 OFF の簡易評価は内部粗スキャン（非公開・候補棄却専用）が
    // 私的設定（個別 bool）で構成する（AstrometryOptions 型自体は残す）。
}
```

- 見かけ位置を返す関数（公開・※署名は要レビュー）:
  - `pub fn apparent_geocentric(eph: &impl Ephemeris, body: Body, time_tt: TtInstant, opts: AstrometryOptions) -> Result<Position<Gcrs>, EphemerisError>`
  - もしくは出力フレームを CIRS/of date まで進める版（ベッセル要素生成側の要求に合わせる。ISSUE-035 連鎖と整合）。返却フレーム型は load-bearing なので `Position<F>` で固定。

## 数式・アルゴリズムの出典
- **光行時間（light-time）**: 観測時刻 t に対し、天体放射時刻 t − τ を反復で解く。τ = |r(t−τ)| / c。月は τ≈1.3s、太陽 τ≈499s。SOFA 参照: `iauLtpequ` ではなく反復は手実装（SOFA の天体測位例 `iauApcg`/`iauAtciq` 群が aberration/light-time を内包。本実装は分解して適用）。c = 299792.458 km/s（IAU 公称, conventions に追記 ※要確認）。
- **歳差（IAU2006）+ 章動（IAU2000A）**: ISSUE-035 が供給する frame bias + 歳差 + 章動の回転（GCRS→CIRS）。出典: Capitaine et al. (2003/2006) IAU2006 precession; IAU2000A nutation; IERS Conventions 2010 ch.5。SOFA: `iauPnm06a`（NPB 行列）/ `iauC2i06a`（celestial→intermediate）相当。
- **光行差（aberration）**: 観測者（地心）の速度による恒星光行差。一次（古典）+ 相対論補正項。観測者速度 = 地球の重心速度（VSOP87D 地球速度 or DE）。SOFA: `iauAb`（相対論的光行差）の式に準拠。式: u' = (u/β + (1+(u·v)/(1+1/γ))·v/c) / (1 + u·v/c)（SOFA `iauAb` 形）。
- **相対論偏向（deflection）**: 太陽重力による光の曲がり。日食では太陽近傍を扱うため寄与は要評価だが、初期は省略可（architecture §5）。SOFA: `iauLd`。省略時は metadata に記録（必須）。
- **適用順序（D3 確定: SOFA `iauAtciq` に固定）**: light-time（放射時刻で暦再評価）→ **GCRS 内で relativistic_deflection（`iauLd` 相当, M1 扱いで既定 OFF）→ aberration（`iauAb` 相当）** → その後 **frame bias + IAU2006 歳差 + IAU2000A 章動（GCRS→CIRS, ISSUE-035）**。すなわち偏向・光行差は GCRS（J2000 軸）座標で先に適用し、bias/歳差/章動は最後にまとめて回す。これは SOFA `iauAtciq`（および `iauAtci13`）の標準順序であり、本実装はこれを確定仕様とする（順序ミスは数″誤差）。

## 単位 / 時刻系 / 座標系
- 単位: 位置 km、速度 km/s、角度ラジアン。c = km/s。
- 時刻系: 入力 `TtInstant`（位置計算標準, conventions §6）。暦評価は TDB へ（TT≈TDB 許容、metadata 記録）。光行時間補正後の放射時刻も TDB 系で暦再評価。
- 座標系: 入力 GCRS（ICRS 軸・地心）→ 出力 GCRS の見かけ位置、または CIRS（of date）まで（ISSUE-035 連鎖、conventions §5）。返却フレームを型で明示。

## アルゴリズム概要
1. opts に従い分岐（Standard は前3つ強制 ON。内部粗スキャン（非公開・候補棄却専用）は緩和、Reference は DE + 全補正）。
2. light_time: τ を反復（初期 τ0 = |r(t)|/c、2–3反復で収束）。放射時刻で暦再評価し方向ベクトル確定。
3. **（D3 順序確定: SOFA `iauAtciq`）** GCRS 内で relativistic_deflection: opts で ON のとき `iauLd`、OFF なら metadata に省略記録（既定 OFF, M1 扱い）。
4. **（同・GCRS 内）** aberration: 地球重心速度で SOFA `iauAb` 式を適用。
5. precession_nutation: ISSUE-035 の GCRS→CIRS 回転を**最後に**適用（frame bias + IAU2006 歳差 + IAU2000A 章動。内部粗スキャン（非公開）は IAU2000B 許容。公開出力は 2000A）。
6. 出力 `Position<F>` と適用補正の記録を返す。

## 受け入れテスト
- **L3 / DE 差分（第一義, accuracy.md §3.1, §3.3 手順）**: DE440 + 同等補正パイプラインを基準に、Standard 補正後の月/太陽見かけ地心方向の角度差を 1900–2100 で測定。**月 ≲0.1″ / 太陽 ≲0.05″ 級**の残差確認（暦残差と補正残差を層分解, §4）。
- 補正分解テスト: light_time のみ / +歳差章動 / +aberration を段階適用し、各段の寄与量が既知オーダー（月 light-time ≈数″移動、aberration ≈20″ 級, 太陽 aberration ≈20.5″）に一致することを確認。
- standard() が light_time/aberration/precession_nutation = true を返すこと、Standard プロファイルでこれらを OFF にできない（型/エンジン側で固定）ことのテスト。
- 省略補正（deflection OFF）が metadata に記録されることを確認。
- **適用順序固定の回帰テスト（D3）**: 補正を SOFA `iauAtciq` 順（GCRS 内 deflection→aberration → その後 bias+IAU2006 歳差+IAU2000A 章動）で適用した結果を固定し、順序入替（例: 章動後に aberration）との差分が数″オーダーで出ること、および正順が DE 同等パイプラインと一致することを回帰で固定。順序が将来変わったら fail する。
- **相対論偏向の省略上限 実測テスト（M1, reviews M1）**: 日食＝太陽縁近傍配置で deflection ON/OFF の差を一度実測し、省略上限（≪0.05″ 想定）を確認して metadata 注記に反映することをテスト化（無視可の根拠を実測値で残す）。
- 基準値は DE オラクルから動的取得（conventions §11）。

## 許容誤差
- accuracy.md §2.1 RSS 配分のうち本 Issue 担当:
  - 光行時間 + 光行差 = **0.10″**（Standard 必須）。
  - （歳差章動 + フレーム = 0.05″ は ISSUE-035 担当だが本 Issue が連鎖適用）。
- 数値根拠: aberration は ~20″ の系統項で、これを ON にしないと 40s 級（20″×2s/″）の誤差。必須化で残差を 0.10″ 以下（実装精度・速度入力精度律速）に抑える。月 0.1″/太陽 0.05″ の最終残差目標へ RSS 寄与。

## 実装メモ
- light-time 反復は固定回数（2–3）でなく相対収束判定（角度変化 < 許容/10）にする。月は1反復でほぼ収束、太陽も少回数。
- aberration の観測者速度は暦と同一バックエンドの地球速度（StateVector.velocity）を使う。速度が None のバックエンドでは対称差分速度を要求（ISSUE-013/014 が供給。差分幅はテスト決定 accuracy.md §3.3）。
- 補正の適用順序は **D3 で確定済み**: SOFA `iauAtciq` 順（GCRS 内 deflection→aberration → その後 bias+IAU2006 歳差+IAU2000A 章動）。実装はこの順に固定し、上記回帰テストで保護する（順序ミスは数″誤差）。SOFA 関数の正確な式番号・引数規約は一次資料（SOFA Cookbook）で最終確認（要確認）。
- 相対論偏向（deflection）は **M1 扱い**（reviews M1, 既定 OFF）。accuracy.md §2.1 では明示配分なし（≪0.05″ 想定）。日食＝太陽縁近傍のため寄与上限を一度実測し（上記 M1 受入テスト）、無視可を metadata 注記付きで確定する。
- Standard 強制 ON は EngineConfig::standard()（umbra-eclipse, ISSUE 別）側でも二重に固定。
