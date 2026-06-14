# ISSUE-017: Conjunction solver（地心合の精解・連続関数化＋粗走査→Brent）

- crate: umbra-eclipse
- 依存: ISSUE-016（新月候補窓）, ISSUE-008（Brent root solver）, ISSUE-015（見かけ地心位置）, ISSUE-012（Ephemeris）, ISSUE-002（Radians 正規化）, umbra-core（TtInstant, 角度連続化）
- モード(tdd-workflow): strict（朔の精解時刻はベッセル要素の基準時刻・最大食概算の起点。求解の連続性・収束（無条件 Newton 禁止, conventions §11）が後段精度を律速。strict）

## 目的
ISSUE-016 の各候補窓内で、太陽と月の**地心合**（黄経差または赤経差のゼロ点）を Brent 法で精密に解き、合の TT 時刻を返す（architecture §3, §12）。
- 黄経差／赤経差を **±π 折返しを除いた連続関数**へ変換してから解く（conventions §2, architecture §12）。
- 粗走査で符号変化区間を検出 → Brent でゼロ点（ISSUE-008）。**Newton 単独禁止**（conventions §11）。

## 非目的
- Brent 求根器そのものの実装（ISSUE-008 を利用）。
- 見かけ位置補正の実装（ISSUE-015 を利用）。
- 日食可能性判定・早期棄却（ISSUE-018）。合は日食でなくても起きる（毎朔）。
- 最大食時刻（ベッセル面上の最接近）の精解 = ISSUE-023/局地。本 Issue の「合」は最大食の良い初期値だが同一ではない。

## 公開インターフェース
crate 内部 API 中心（`pub(crate)`）。検証用に `pub` 化候補。

```rust
/// 合の種別（どの座標で合を定義するか）。NASA/Meeus 慣習に合わせ選択。
pub(crate) enum ConjunctionKind {
    EclipticLongitude,  // Δλ = λ_moon − λ_sun = 0（黄経合, Meeus Ch.54 の伝統的定義）
    RightAscension,     // Δα = α_moon − α_sun = 0（赤経合）
}

pub(crate) struct Conjunction {
    pub time_tt: TtInstant,
    pub kind: ConjunctionKind,
    pub separation: Radians,   // 合時刻の月-太陽 角距離（早期棄却 ISSUE-018 の入力）
}

pub(crate) fn solve_conjunction(
    eph: &impl Ephemeris,
    candidate: &NewMoonCandidate,   // ISSUE-016 の窓
    kind: ConjunctionKind,
    config: RootConfig,             // ISSUE-008（x_tolerance, max_iterations）
) -> Result<Conjunction, EclipseError>;  // RootNotBracketed / SolverDidNotConverge を伝播
```

- `RootConfig`（ISSUE-008）の `x_tolerance` は時刻（日 or 秒, 境界で `JulianDate2` 差分へ）。`root_tolerance_seconds`（EngineConfig）を目標の 1/10 以下で渡す（accuracy.md §2.1）。
- 失敗は `EclipseError::RootNotBracketed` / `SolverDidNotConverge`（api-draft §3.5, From<SolverError>）。

## 数式・アルゴリズムの出典
- **合の定義**: Meeus, *Astronomical Algorithms* (2nd ed.), **Ch.54「Eclipses」冒頭** および Ch.49（朔）。日食判定の伝統定義は黄経合（`λ_moon = λ_sun`）。NASA/Espenak のベッセル要素は最大食付近を扱うため赤道座標ベースだが、**朔の検出は黄経合で行い、最接近（最大食）は別途**（ISSUE-023）という分担が一般的。採用座標を実装コメントで明記（conventions §10）。
- **連続関数化**: `f(t) = (λ_moon(t) − λ_sun(t))` を `[-π, π)` 正規化したものは合付近で連続だが、窓全体で折返しが入りうる。窓内では `Δλ` は単調に増加（月が太陽を追い越す）するため、**窓中心の値を基準に unwrap して連続化**（conventions §2「±π 折返しを除いた連続関数」, architecture §12）。出典: 一般的な root-finding 前処理（NR の連続化）。
- **角速度の符号**: 月の黄経角速度 ≈ 13.2°/日 ≫ 太陽 ≈ 0.99°/日 → `Δλ` は単調増加。これを符号変化検出・ブラケット方向の根拠とする（Meeus Ch.47 月運動 / Ch.25 太陽運動）。

## 単位 / 時刻系 / 座標系
- 時刻系: 入力窓・出力ともに **TtInstant**（conventions §6）。Brent の独立変数は窓内オフセット（日 or 秒, `JulianDate2` 差分で橋渡し, ISSUE-008 注記）。
- 角度: ラジアン。`Δλ`/`Δα` の正規化は**連続化用に専用**（two_pi でも signed でもなく unwrap, conventions §2「用途ごとに正規化関数を分け混在させない」）。
- 座標系: 見かけ地心（ISSUE-015）。黄経合は黄道 of date、赤経合は赤道 CIRS/of date。返す `separation` は角距離（acos 由来は ISSUE-018 と整合, [-1,1] クランプ）。

## アルゴリズム概要
1. 窓 `[t0, t1]`（TT）を等間隔で粗走査（刻みは平均月運動から「窓内に符号変化が1回だけ入る」よう設定し、根拠コメント化）。
2. 各サンプルで `f(t) = Δλ_continuous(t)` を評価（ISSUE-015 で見かけ位置）。窓中心基準で unwrap し連続化。
3. 隣接サンプルで符号変化（`f(t_i)·f(t_{i+1}) < 0`）を検出 → ブラケット `[t_i, t_{i+1}]` 確定。
4. Brent（ISSUE-008）でゼロ点を `x_tolerance` まで精解。**Newton 単独は使わない**（conventions §11）。
5. 合時刻で月-太陽角距離 `separation` を `acos(clamp(û_moon·û_sun, -1, 1))` で算出（accuracy.md §2.2 クランプ）。
- 数値安定性: 符号変化が窓内で 0 回（窓不足＝ISSUE-016 のマージン不良 → `RootNotBracketed` で顕在化）、複数回（粗走査刻みが粗すぎ）を異常として検出。acos 引数クランプ必須。Newton 禁止（conventions §11）。

## 受け入れテスト
accuracy.md テストレベル **L2/L3（時刻・位置）境界の合検出**、補助 L1（連続化）。
- **DE 差分（第一義, accuracy.md §3.1）**: DE440 で同一の合関数を解いた朔時刻と、解析暦（Standard）の朔時刻の差を 1900–2100 で測定。差は暦残差由来（§4 層分解）として記録。
- **NASA カタログ整合（第二義, data-sources §4.1）**: 既知日食の朔/合時刻と整合（ΔT・座標慣習を揃えて比較。絶対基準にしない, accuracy.md §3.1）。
- **MockEphemeris 人工ケース（accuracy.md §3.1）**: 既知の合時刻を持つ人工配置（線形に動く月・太陽）で、解析解とソルバ結果を厳密比較（オラクル＝人工配置の数式, 実装からコピーしない）。
- 連続化テスト（L1）: 窓が ±π 境界をまたぐケースで unwrap 後の `f` が連続・単調であること、Brent がブラケットを失わないこと。
- 異常系: 窓に符号変化なし → `RootNotBracketed`。`max_iterations` 不足 → `SolverDidNotConverge`。
- 収束: root_tolerance を目標の 1/10（≤0.15s 相当, accuracy.md §2.1）に設定し達成を確認。

## 許容誤差
- accuracy.md §2.1「solver 収束 0.05″」「root_tolerance を目標の 1/10 以下」。合時刻の独立変数許容は **≤ root_tolerance_seconds（目標 ±1.5s の 1/10 = ≤0.15s 相当）**。
- 感度: 角速度 ≈0.5″/s（§2.1）。合時刻 0.15s 誤差 ≈ 0.075″ ≪ 0.05″ バジェットに収まるよう、呼出側は `root_tolerance_seconds` を厳しめに設定可。
- 注意: 合時刻は**最大食時刻ではない**（最大食はベッセル面の最接近, ISSUE-023）。本 Issue の許容は合のゼロ点求解そのもの。最大食精度は ISSUE-023 の責務。
- 許容を「通すためだけに拡大しない」（conventions §11）。

## 実装メモ
- 連続化は「窓中心の `Δλ` を基準に `±2π` を引き去る」方式が単純。窓が小さい（ISSUE-016 で ±1日程度）ので折返しは高々1回。
- 黄経合 vs 赤経合の選択は後段（ISSUE-018 棄却・ISSUE-023 種別）の要求に合わせる。**黄経合を既定**とし、赤経合は照合用（NASA 表記突合）。レビューで確定。
- 粗走査刻みは月運動（13.2°/日）から「窓内で `Δλ` が単調 1 通過」を保証する値に。刻み過大で根の見落とし＝偽陰性なので慎重に（ISSUE-016 と同じく偽陰性は致命）。
- 速度が要る場合（収束加速）でも Newton 単独は禁止。Brent の超線形収束で十分（ISSUE-008）。
- レビュー重点: 連続化の正しさ、ブラケット保証、Newton 非混入、acos クランプ、合≠最大食の混同防止。
