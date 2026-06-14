# §7 ベッセル多項式（Step8）

> 正本: `algorithms.md §0`（記号表 = x,y,d,μ,l1,l2,tan f1,tan f2）・`numerical-policy.md §A4`（**Chebyshev 最小二乗フィット**・ノード・次数選択・残差ゲート）・`conventions.md §2`（μ 連続化）・`§6`（時刻系）・`accuracy.md §2.1`（多項式 fit <0.10″）/`§3.2`（L7 直接 vs 多項式）。本セクションはこれらに**厳密準拠**する。
> 状態: ドラフト（Milestone 0）。一次資料で未確認の式番号・符号は「**要確認**」を残し、推測で式番号を書かない。
> 関連 Issue: ISSUE-022（ベッセル多項式 fit）。サンプリング元: ISSUE-037（直接評価器）/ ISSUE-021（瞬時要素定義）。

---

## 目的と入出力（単位・時刻系・座標系・フレーム）

最大食付近の時系列で瞬時ベッセル要素（§6/ISSUE-021/037）をサンプリングし、各成分（x, y, d, μ, l1, l2）を時間多項式へ**Chebyshev 最小二乗フィット**して `BesselianPolynomial` を生成する（architecture §6.1, api-draft §3.3, ISSUE-022）。経路（中心線・限界線）/GeoJSON/NASA 形式エクスポート用の高速供給源。局地の既定は直接評価（§6/ISSUE-037, fit 誤差ゼロ）であり、本層は経路/エクスポート用途（numerical-policy §A4）。

| 項目 | 内容 |
|---|---|
| 入力 | `source: impl BesselianSource`（**ISSUE-037 直接が基準**, fit 誤差ゼロ）, `epoch_tt`(t0), `fit_interval`, `degree`（NASA 低次から）, `tolerance: BesselFitError` |
| 出力 | `BesselianPolynomial`: x,y,d,μ,l1,l2 の `Polynomial`, tan_f1/tan_f2（定数, NASA 慣習）, `fit_interval`, **`fit_error: BesselFitError`（残差を必ず保持）** |
| 単位 | x,y,l1,l2 = **Re 無次元**, d,μ = rad（§6 と同一）。多項式変数 = epoch_tt からの経過時間（NASA 慣習＝**時 hour** 候補, 要確認・実装で固定しコメント） |
| 時刻系 | **TT 基準**（conventions §6, epoch_tt 起点）。μ の恒星時部分のみ UT1 由来（§6 B5 と同断り） |
| 座標系/フレーム | **FundamentalPlane**（conventions §5, §6 と同一）。供給源差し替え（ISSUE-037 と）で単位/座標系不変（同一 `BesselianSource` 契約） |

**確定方針（厳守, numerical-policy §A4）**: 正規化時刻 `τ ∈ [−1,1]` 上の **Chebyshev 最小二乗フィット**（単項式 Vandermonde 直接 fit は悪条件ゆえ禁止）。出力時に NASA 互換単項式へ変換。残差 `BesselFitError` を必ず保持し、許容超は `BesselFitExceededTolerance`（誤差を隠さない, conventions §11）。

---

## 記号（algorithms.md §0 参照。本節固有の補助記号のみ定義）

§0「ベッセル要素」表を正本とする。本節固有:

| 記号 | 意味 | 単位 |
|---|---|---|
| `t0` | 多項式の基準時刻 = `epoch_tt`（NASA 形式, 最大食付近） | TT |
| `t` | t0 からの経過時間（NASA 単項式変数。単位 hour 候補, 要確認） | hour（or s/日, 実装固定） |
| `τ` | 正規化時刻 = `(t_inst − t0)/Δ ∈ [−1,1]`（Chebyshev 基底変数, numerical-policy §A4） | 無次元 |
| `Δ` | fit 区間の半幅（τ スケーリング分母。`fit_interval` から） | 同 t |
| `T_n(τ)` | n 次 Chebyshev 多項式（第一種） | 無次元 |
| `c_n` | Chebyshev 係数（成分ごと） | 成分単位 |
| `N_deg` | 採用多項式次数（NASA 低次=3 から開始, 上限あり） | 無次元 |
| `BesselFitError` | 残差（max_x, max_y, max_l1, max_l2） | Re 無次元 |
| `s ∈ {x,y,d,μ,l1,l2}` | fit 対象成分（tan f1/f2 は定数扱い） | 各成分単位 |

---

## 数式（番号付き・各式に出典）

### 7.1 NASA ベッセル多項式形式（出力形）

**(P1) 単項式多項式（NASA 形式）**
出典: Espenak の Besselian Elements（NASA GSFC / NASA TP-2006-214141, data-sources §4.1）。要素を t0 からの時間 t の多項式で表す（通常 **3 次**）:
```
s(t) = s0 + s1·t + s2·t² + s3·t³          ( s ∈ {x, y, d, μ, l1, l2} )
tan f1, tan f2 = 定数                       ( NASA 慣習 )
μ(t) は t の 1 次が主（地球自転 ≈15°/hour, §6 D5 μ′）
```
- 時間単位は NASA 慣習＝**時 hour**（**要確認**, 一次で確定。実装で固定しコメント, conventions §10）。t0 起点も固定。

### 7.2 Chebyshev 最小二乗フィット（内部基底, numerical-policy §A4）

**(P2) 正規化時刻**
出典: numerical-policy §A4（τ ∈ [−1,1] Chebyshev）。
```
τ = ( t_inst − t0 ) / Δ ∈ [−1, 1]          # 区間端を ±1 に写す
```
- 単項式（Vandermonde）直接 fit は**悪条件**ゆえ禁止。Chebyshev で条件数を抑える（numerical-policy §A4）。

**(P3) Chebyshev 展開と最小二乗**
出典: numerical-policy §A4 / 標準 Chebyshev LS（Numerical Recipes 多項式 fit）。
```
s(τ) ≈ Σ_{n=0..N_deg} c_n · T_n(τ)
c_n = 最小二乗解（normal equations or QR。Chebyshev ノードなら直交性で安定）
```

**(P4) Chebyshev ノード**
出典: numerical-policy §A4（Chebyshev ノードで Runge 現象回避）。瞬時要素（§6/ISSUE-037）を Chebyshev ノードで評価（または密な一様＋LS）:
```
τ_j = cos( π (j + 1/2) / M ),  j = 0..M−1       # M 個のノード（M ≥ N_deg+1, 余裕を持つ）
t_inst,j = t0 + Δ · τ_j                          # 対応 TT で source.at() 評価
```

**(P5) NASA 単項式への変換**
出典: numerical-policy §A4（出力時に NASA 互換単項式へ変換, D5）。`τ = (t − t0)/Δ` を代入し `T_n(τ)` を `t` の冪に展開して `s0..s3`（P1）を得る:
```
Σ c_n T_n( (t−t0)/Δ )  ⇒  s0 + s1 t + s2 t² + s3 t³      ( Clenshaw/基底変換 )
```
- μ は度・hour 系（NASA, §6 D5）。変換時に単位を NASA 表記へ揃える（系統差を記録, accuracy.md §0）。

### 7.3 μ の連続化（fit 前処理）

**(P6) μ unwrap**
出典: conventions §2（連続化）/ ISSUE-022。
```
μ_continuous(τ) = unwrap( μ(τ) )          # [0,2π) 折返しを fit 区間内で除去
```
- μ は赤経基準で `[0,2π)` 折返しがあるため、**fit 前に区間内で unwrap**（連続化しないと多項式が破綻, conventions §2, ISSUE-022）。

### 7.4 残差ゲート（BesselFitError）

**(P7) 残差評価**
出典: accuracy.md §3 実測ガード / numerical-policy §A4 / ISSUE-022。fit 区間を密にサンプルし、直接値（§6/ISSUE-037, fit 誤差ゼロ）と多項式値の**最大絶対差**を保持:
```
max_x  = max_{τ∈fit} | x_poly(τ)  − x_direct(τ)  |
max_y  = max_{τ∈fit} | y_poly(τ)  − y_direct(τ)  |
max_l1 = max_{τ∈fit} | l1_poly(τ) − l1_direct(τ) |
max_l2 = max_{τ∈fit} | l2_poly(τ) − l2_direct(τ) |
```
- **次数選択**: NASA 低次（3 次）から開始し、`BesselFitError` が **fit 配分（x,y で <0.10″ 相当）** を満たす最小 `N_deg` を採用（numerical-policy §A4, accuracy.md §2.1）。満たさなければ次数を上げて再 fit、**最大次数を上限**（過剰適合/Runge 防止, 要確認: 上限値は M2 実測, numerical-policy §A4）。
- それでも超過なら `BesselFitExceededTolerance`（api-draft §3.5, conventions「誤差を隠さない」, ISSUE-022）。
- **0.10″ → Re 換算**: x,y の 0.10″ 相当を Re へ換算した値を `tolerance` 既定に（要 M2 実測, accuracy.md §2.1 注）。

---

## 手順（実装順・数値注意）

1. **サンプリング (P4)**: `fit_interval` を Chebyshev ノードで `source.at(t_inst,j)`（**ISSUE-037 直接, fit 誤差ゼロ基準**）評価。最大食を中心に配置。**暦を直接サンプリングしない**（層分解のため必ず直接供給源を経由, accuracy.md §4, ISSUE-022）。
2. **μ 連続化 (P6)**: μ を区間内 unwrap（必須）。
3. **正規化 (P2)**: `τ = (t_inst − t0)/Δ ∈ [−1,1]`。
4. **Chebyshev LS (P3)**: 各成分 `s ∈ {x,y,d,μ,l1,l2}` を N_deg=3 から fit。tan f1/f2 は区間平均（定数, NASA 慣習。低次にするかは残差で決定, 要レビュー ISSUE-022）。
5. **単項式変換 (P5)**: NASA 形式 `s0..s3` へ。μ は度・hour へ（D5）。
6. **残差ゲート (P7)**: 密サンプルで `BesselFitError` 測定 → 許容超なら次数+1 再 fit（上限まで）→ なお超なら `BesselFitExceededTolerance`。
7. **構成**: `BesselianPolynomial`（`fit_error` を必ず同梱）。`at()` は **Horner 評価**（architecture §12）で `BesselianSource` を満たし、§6/ISSUE-037 と差し替え可能。

**fit 区間 Δ（numerical-policy §A4）**: 中心食前後（おおむね P1..P4 か少し広め, §8/ISSUE-023）。区間端での発散を残差ゲート（P7）で検出。

**数値注意（横断, numerical-policy §A4）**:
- **Chebyshev 必須**（単項式 Vandermonde 直接 fit は悪条件ゆえ禁止）。τ ∈ [−1,1] スケーリングで条件数を抑える。
- **μ unwrap 必須**（折返しで fit 破綻, conventions §2）。
- 残差ゲートを「通すため」に緩めない（`fit_error` で必ずガード, conventions §11）。`fit_error` は常に真の残差を報告（誤差を隠さない）。
- Horner 評価（architecture §12）。

---

## 境界・特異・異常系

- **fit 残差許容超**: `BesselFitExceededTolerance`（fit 区間過大・次数不足・端発散, ISSUE-022）。
- **μ の ±π/2π 折返し**: unwrap が効かないと fit 破綻 → 区間内 unwrap で連続化（L1 テスト, ISSUE-022）。
- **Runge 現象/過剰適合**: 次数上限と Chebyshev ノードで防止（numerical-policy §A4）。一様サンプルの高次は Runge → ノード使用。
- **fit 区間外の `at()`**: `fit_interval()` 外は範囲チェック（多項式は区間内のみ妥当, ISSUE-022/037）。
- **`fit_error` 未充填**: 常に非ゼロで埋め結果に同梱（誤差を隠さない, conventions §11, ISSUE-022）。
- **tan f1/f2 を定数にする近似**: 区間で変動が大きい場合は残差に現れる → 低次多項式化を検討（NASA 慣習＋残差で決定, ISSUE-022）。

---

## 検証（基準値の出典。実装へ値コピー禁止 = conventions §11）

accuracy.md §3.2（L7）、ISSUE-022 受入テスト準拠。基準値は直接供給源（§6/ISSUE-037）/NASA から取得し**実装へハードコードしない**（conventions §11）。

- **直接 vs 多項式 残差（最重要, L7, accuracy.md §3.2, architecture §6.1）**: ISSUE-037（直接, fit 誤差ゼロ）を基準に `BesselianPolynomial.at(t)` の x,y,l1,l2 残差を fit 区間で実測。`fit_error` が実残差を正しく報告し、許容超で `BesselFitExceededTolerance` を返すこと。profile 毎の採用（経路=多項式可, 局地=直接既定）を残差で最終決定。
- **NASA 係数比較（第二義, data-sources §4.1, 品質基準「係数比較＋瞬時値比較の両方」）**: 既知日食で生成多項式係数（x0,x1,x2,x3,...）を NASA 公開ベッセル多項式と比較（**k/ΔT 慣習を揃える**, accuracy.md §3.1）。係数比較と評価瞬時値比較の両方。
- **MockEphemeris（accuracy.md §3.1）**: 人工配置でサンプリング→fit→評価が直接値に一致（fit 誤差が許容内）。
- **fit_error 保持テスト**: `BesselFitError` が必ず非ゼロで埋まり結果に同梱（誤差を隠さない, conventions §11）。
- **許容超ケース**: fit 区間を過大にして残差を悪化させ `BesselFitExceededTolerance` を確認。
- **μ 連続化テスト（L1）**: μ が ±π/2π 境界をまたぐ区間で unwrap が効き fit が破綻しないこと。
- **基準値コピー禁止**: GPL 実装/NASA 値を期待値に貼らない（基準は直接供給源 ISSUE-037 = 数式オラクル系, conventions §11）。

---

## 許容誤差

- accuracy.md §2.1「多項式 fit（使う場合）**<0.10″**」。fit 区間で残差を実測ガード（§3）。x,y の 0.10″ 相当を Re 換算した値を `tolerance` 既定に（要 M2 実測, accuracy.md §2.1 注）。
- **直接（ISSUE-037, fit 誤差0）vs 多項式の残差 = L7 サブテスト**（accuracy.md §3.2）。profile 毎の採用を残差で決定（architecture §6.1）。
- 食分 0.001 ≈ 1.9″（accuracy.md §2.2）に対し l1,l2 の fit 残差は十分小さく（<0.10″ 相当）。
- 許容を通すための拡大禁止（conventions §11）。`fit_error` は常に真の残差を報告。
- 本層は**経路/エクスポート用途**であり、局地接触の精度は直接評価（§6/ISSUE-037, fit 誤差ゼロ）が担保（精度最優先, architecture §6.1）。

---

## 出典

- Espenak/NASA: Besselian Elements 多項式形式（GSFC eclipse / NASA TP-2006-214141, data-sources §4.1, 第二義照合）。**要確認**: 多項式時間単位（hour）・3 次慣習・tan f 定数扱いの正確な記述。
- Numerical Recipes（多項式最小二乗 / Chebyshev fit / Clenshaw 評価）。numerical-policy §A4（Chebyshev・ノード・次数選択・残差ゲート・NASA 単項式変換）。
- conventions §2/§5/§6/§10/§11, accuracy.md §0/§2.1/§2.2/§3.1/§3.2/§4, algorithms.md §0。
- §6（瞬時ベッセル要素 ISSUE-021・直接評価器 ISSUE-037・NASA 対応表 D5）, §8（fit 区間 P1..P4 ISSUE-023）, architecture §6.1/§12（Horner）。
- 関連 Issue: ISSUE-022, ISSUE-037（サンプリング元）, ISSUE-021（定義）, ISSUE-008（数値基盤）。
- **要確認**: NASA 多項式の時間単位（hour vs 日/秒）と t0 起点（P1）。多項式最大次数の上限（P7, M2 実測, numerical-policy §A4）。tan f1/f2 を定数にするか低次にするか（残差で決定）。
