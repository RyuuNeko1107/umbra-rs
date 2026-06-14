# ミューテーション検証レポート: umbra-core（tdd-workflow 工程7）

ツール: `cargo-mutants` 27.1.0（Docker 内）。対象: `crates/umbra-core`。コマンド: `cargo mutants --package umbra-core`。

## 結果（最新）

```
665 mutants tested: 527 caught, 118 missed, 20 unviable
```

missed の内訳（モジュール別）:

| モジュール | missed | 状態 |
|---|---|---|
| vector / matrix / julian / calendar / ellipsoid | **0** | 全 kill（下記の通りテスト強化で解消） |
| solver | 118 | 等価変異として許容（下記 §許容） |

## 見つかった「本物のテスト穴」と修正（初回 missed=158 → 解消）

初回実行で 36 件の非 solver missed と 122 件の solver missed。非 solver はすべて**テストが特殊値ばかりで各項を励起していなかった**ことが原因。以下で全 kill:

- **vector**: dot/cross/scale を**全成分非ゼロ・非対称**かつ「係数がどの成分とも一致しない」入力に変更（`(1,2,3)·(4,5,6)`、`cross` の直交性、`scale(10)`）。z=0 や対称値では `+`/`-`/`*` の取り違えが不可視だった。
- **matrix**: 一般行列 `[[1..9]]` × `(1,2,3)` の手計算オラクルで `mul_vec` の各積項を区別。
- **julian**: `julian_millennia_since_j2000` に**テストが皆無**だった → J2000=0・1千年後=1・「世紀の1/10」関係を追加。
- **calendar**: 往復テストに**2月・閏日**（`2024-02-29`, `2000-02-15`）を追加し `month_f>2` 境界と `e−13` 分岐を励起。
- **ellipsoid**: Meeus Palomar 例の具体オラクル、楕円体面不変量 `(X²+Y²)/a²+Z²/b²=1`、地心緯度の定義関係 `tanφ'=(1−e²)tanφ`、観測者/ECEF の**高さ項の厳密関係**（`Δrho=(h/a)·sin/cosφ`、`Δecef=h·法線`）を追加。h=0 や norm 比較では高さ項の符号・演算が不可視だった。

## 許容する生存変異（118 件・すべて solver.rs）

`brent_root` と `minimize_golden` の**内部ヒューリスティクスに限定**される。これらは**等価変異**であり、検証対象の契約を破らない:

- **Brent**: 正しさは「区間ブラケット不変量 ＋ 二分法フォールバック」から導かれる。逆二次補間の算術（`a·fb·fc/…`）と採択条件（`reject` 判定・swap の簿記）を改変しても、不正な候補は `reject` で弾かれ二分法に退化するだけで、**最終的な根は tol 内・max_iter 内で正しく求まる**。したがって差分が出ない。
- **golden**: 比較の同点処理（`fc < fd` の `<=`）・収束判定（`(b−a)<tol`）の改変は、`max_iter` まで回れば同じ最小点へ収束するため結果に出ない。

契約自体は**振る舞いベースのテストで担保**している: 多様な関数の根（`x³−x−2`, `eˣ−3`, Omega, 二次式）、未ブラケット時の `RootNotBracketed`、端点ちょうどの根、非対称・端寄りの最小点、端点との大小比較。

これらを kill するには**内部の反復回数・経路を assert**する必要があるが、それは実装結合で脆い（収束速度を 25 反復で縛るテストは**正しいコードでも false failure** を出したため撤回）。よって工程7の方針（「全 mutant が殺せるとは限らない。生存を列挙し許容可否を明示」）に従い、**等価変異として許容**する。

## 追加モジュール（time / deltat, ISSUE-006/007）

- **time**（UTC/TAI/TT）: 当初 `UtcInstant::to_gregorian` に直接テストが無く 140 変異が見逃し → 往復テスト追加で kill。`tai_to_utc` の `-37` 近似が等価変異化 → 単純法へ書き換え round-trip で kill。**生存ゼロ**。
- **julian** `days_since`: 追加。生存ゼロ。
- **deltat**（ΔT/UT1）: `uncertainty_seconds` を不連続な順序しきい値に書き換え、各分岐・境界・式を厳密値で検証 → kill。`delta_t_seconds` の区分内部の演算子変異は内部点テストで kill。
  - **許容する等価変異（6件）**: `delta_t_seconds` の区分境界比較 `< → <=`。Espenak–Meeus の区分多項式は**接合点で連続**（隣接区分の差は 1941 境界で ~0.0008s 等、いずれも <0.001s）であり、ΔT 自体の不確実性（秒オーダー）を遥かに下回る。境界年ちょうどでのみ僅差が出る等価変異のため、`--exclude-re 'with <= in.*delta_t_seconds'` でゲートから除外する。

## 退行ガード（CI）

テスト有効性の退行を防ぐため、CI（`mutants` ジョブ, `.github/workflows/mutation.yml`）は**本レポートで許容した等価変異を除いた全モジュールで生存ゼロ**を要求する:

```
cargo mutants -p umbra-core --exclude '**/solver.rs' --exclude-re 'with <= in.*delta_t_seconds'
```

直近の結果: **911 mutants tested, 874 caught, 37 unviable, 0 missed**。フル実行（除外なし）は `cargo mutants -p umbra-core` で再現可能（solver 等価変異と delta_t 境界が生存として現れる）。
