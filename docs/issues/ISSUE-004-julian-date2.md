# ISSUE-004: JulianDate2（2要素表現）

- crate: umbra-core
- 依存: ISSUE-001
- モード(tdd-workflow): strict（永続形式・公開仕様。全時刻系の内部表現の土台で、桁落ちが全時刻変換精度を律速するため strict）

## 目的
ユリウス日を 2 要素 `JulianDate2 { part1, part2 }` で表現し、巨大 JD と微小差を 1 つの `f64` に押し込まない（conventions §6 / architecture §2）。
- 正規化（part1 を整数寄り or 既定基準に、part2 を `[0,1)` or `[-0.5,0.5)` に再配分）。
- 2 要素同士の加減算（秒/日の微小増分を精度を保ったまま加える）。
- 2 要素を保ったままの比較・差分（`(a - b)` を高精度日数で返す）。

## 非目的
- 暦⇔JD 変換（ISSUE-005）。
- 時刻系（UTC/TAI/TT/UT1/TDB）の意味付け（ISSUE-006/007 と上位型）。`JulianDate2` は時刻系を持たない純数値表現。
- 閏秒・ΔT の適用。

## 公開インターフェース
api-draft §1.3 を転記・具体化:
```rust
#[derive(Clone, Copy, Debug, PartialEq, PartialOrd)]
pub struct JulianDate2 { pub part1: f64, pub part2: f64 }

impl JulianDate2 {
    pub fn new(part1: f64, part2: f64) -> Self;
    /// 高精度な再配分（part1 を基準へ、part2 を小数部へ寄せる）
    pub fn normalized(self) -> Self;
    /// 日数（SI 日）を高精度に加える
    pub fn add_days(self, days: f64) -> Self;
    /// 秒を高精度に加える（内部で /86400、桁分離）
    pub fn add_seconds(self, seconds: f64) -> Self;
    /// 高精度差分（日）。a.diff_days(b) = (a - b) [日]
    pub fn diff_days(self, other: Self) -> f64;
    /// 合算（互換目的。精度が必要な所では part を保つ）
    pub fn to_f64(self) -> f64;  // part1 + part2（桁落ちあり・表示/概算専用）
}
```

## 数式・アルゴリズムの出典
- 2 要素 JD 設計と分割規約: **IAU SOFA** の `(date1, date2)` 二引数規約に準拠（SOFA `iauJd2cal` / `iauCal2jd` 等が採用。`date1+date2 = JD`、可搬な分割を許容）。
- 桁分離加算: SOFA 内部および一般的な **two-sum / compensated summation**（Dekker 1971, Knuth）の考え方。part2 への微小量加算で part1 の有効桁を消費しない。
- 要確認: part1/part2 の既定基準（SOFA 慣習は part1 = 2400000.5 等の固定 or JD 整数部）。本プロジェクトの「正規化」既定（part1 整数寄り / part2 小数部）をコメントで一意に固定。

## 単位 / 時刻系 / 座標系
- 入力/出力単位: 日（part1, part2）。`add_seconds` は秒入力（SI 秒）。
- 時刻系: なし（時刻系非依存の数値表現）。上位の `*Instant` 型がスケールを付与。
- 座標系: 無関係。

## アルゴリズム概要
1. `new`: そのまま保持。`normalized` で再配分。
2. `add_seconds(s)`: `days = s / 86400.0`（86400 は SI 日定数＝magic number ではなく定義値、定数化しコメント）。part2 に加え、必要なら再正規化。**part1 を直接いじらない**。
3. `diff_days`: `(self.part1 - other.part1) + (self.part2 - other.part2)`。大きい part1 同士を先に引いて桁落ちを抑える順序を固定。
4. `to_f64`: `part1 + part2`。**桁落ちが起きうる**旨をドキュメントし、計算経路では使わせない（表示/概算限定）。
- 数値安定性: 微小増分は必ず part2 側へ。連続加算で part2 が肥大したら `normalized` で part1 へ繰り上げ。禁止: 1 個の f64 JD への退避、生 f64 JD の関数間受け渡し（conventions §6 / §11）。

## 受け入れテスト
accuracy.md テストレベル **L2（時刻）** の基盤。
- 精度テスト（オラクル＝高精度参照: 同一値を `f128`/有理数 or SOFA 同等手順で算出、実装からコピーしない）:
  - 現代 JD（≈ 2.46e6）に 1 マイクロ秒を `add_seconds` し、`diff_days` で 1µs（= 1e-6/86400 日）が**復元できる**こと。単一 f64 JD ではこの分解能が出ない（対比テスト）。
  - `add_seconds(86400)` ≡ `add_days(1)` ≡ part 繰上げ後に元 +1 日。
  - `normalized` 後も `to_f64` が（許容内で）保存。
- 境界値: 非常に大きい part1（DE441 範囲 -13200〜+17191 年相当の JD）、負の JD、part2 が `[0,1)` 外の入力 → 正規化で吸収。
- 異常系: `NaN`/`inf` の伝播仕様を固定。
- プロパティ（L8）: `(jd.add_seconds(s)).add_seconds(-s) ≈ jd`（往復、許容内）。

## 許容誤差
accuracy.md §2.1「solver 収束 0.05″」「最大食 ±1.5s」、conventions/§6 の時刻方針から逆算:
- 時刻変換の合成目標は後続 issue で「µs オーダ」を狙う（角速度 0.5″/s より 1µs ≈ 5e-7″ で完全に無視可）。よって `JulianDate2` の加減算・差分は **相対誤差 ≤ 数 ULP**、現代 JD で **絶対 ≤ 1e-11 日（≈ 1µs 未満）** を保証。
- 根拠: 1µs = 1.157e-11 日。f64 の現代 JD 直接表現の分解能は ~約 10µs（仮数 52bit, 2.46e6 日 → ULP ≈ 5e-10 日 ≈ 43µs）。2 要素化で part2 が µs を保持できることをこの誤差で担保する（= 2 要素表現の存在意義）。

## 実装メモ
- `to_f64` は「使ってはいけない場所」を doc で明示（計算経路禁止、conventions §11 の生 f64 時刻混在禁止に連なる）。
- `PartialOrd` は 2 要素合算順序での比較になるため、`normalized` 済み前提か、比較時に内部で差分を取る実装かを固定（推奨: `diff_days` の符号で比較し、桁落ちを避ける）。
- レビュー重点: add 系で part1 の有効桁が削られないか（two-sum 的補償の有無）。SOFA 二引数規約との整合コメント。
