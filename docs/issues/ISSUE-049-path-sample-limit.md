---
milestone: M9
---

# ISSUE-049: `EclipseEngine::path` のサンプル数上限（ハング防止・入力検証）

## 背景（実測された欠陥）

ISSUE-048（`umbra path` CLI）の実装レビューが発見した。`PathOptions::sample_interval_seconds` が
**正の有限値でも極端に小さい**（例 `1e-300`）と、`trace_central` / `trace_penumbral_limits` /
`trace_terminator_lunes` / `build_partial_limit` の走査回数が `span / interval` ＝事実上無限になり、
**`path()` が返らない**（エラーも進捗も出ないまま停止しない）。

ISSUE-048 では **CLI 境界**で `check_path_sample_count` を入れて塞いだが、
**`EclipseEngine::path` は公開 API であり、CLI を経由しない任意の呼び出し元が同じ入力でハングする**。
本 issue でエンジン側に正本の検証を置く。

## 目的

`path()` が**必ず有限時間で返る**ようにする。過大なサンプル数を要求されたら、走査を始める前に
**エラーを返す**（黙ってクランプしない＝conventions §11「誤差・制限を隠さない」）。

## 非目的

- サンプリング方式そのものの変更（等間隔サンプルのまま）。
- `interval` が**非正**のときの挙動変更。これは既に定義済み（「非正なら始点のみ」＝
  `trace_penumbral_limits` 等のループ規約）で、ハングしない。本 issue は**正の極小値**が対象。
- CLI 側の検証の削除（早期・親切なエラーとして残す。定数はエンジンの正本を参照して二重管理しない）。

## 公開インターフェース

```rust
// umbra-eclipse
pub const MAX_PATH_SAMPLES: f64 = 100_000.0;

pub enum EclipseError {
    // 追加
    PathIntervalTooSmall { interval_seconds: f64, estimated_samples: f64 },
}
```

## 確定仕様

1. **検査位置**: `path()` の冒頭、いかなる走査も始める前。
2. **検査対象の span**: `path()` が実際に走査する区間の**最大**を用いる。
   - 中心食（`central_begin`/`central_end` 両 `Some`）なら U1〜U4。
   - 部分食域（`partial_begin`/`partial_end` 両 `Some` かつ `include_limits`）なら P1〜P4。
   - 両方該当するなら**長い方**（P1〜P4 ⊇ U1〜U4）。どちらも該当しないなら走査しないので**検査しない**。
3. **判定**: `estimated_samples = span_seconds / interval_seconds` が `MAX_PATH_SAMPLES` を
   **超える**（厳密 `>`）なら `Err(EclipseError::PathIntervalTooSmall { .. })`。
   `interval_seconds` が**非正**のときは検査しない（既定義の「始点のみ」挙動に委ねる・ハングしない）。
   `interval_seconds` が **NaN** のときは `estimated_samples` も NaN となり比較は false なので素通りするが、
   NaN interval では `t_sec` が NaN になりループの `>=` 比較が false のまま**ハングしうる**ため、
   **NaN は明示的に弾く**（`estimated_samples` に NaN を載せて返す）。
4. **黙ったクランプはしない**: `interval` を上限に丸めて続行すると、利用者が要求した分解能と
   異なる結果を**無言で**返すことになる（誤差を隠す）。必ずエラーにする。
5. **CLI（ISSUE-048）**: 既存の `check_path_sample_count` は残すが、上限定数は
   `umbra_eclipse::MAX_PATH_SAMPLES` を参照し、CLI 側に数値を重複定義しない。

## 受け入れテスト戦略

FAST（合成 bessel・実エンジン不要）:
- 極小 `interval`（上限を超える）で `path()` が**即座に** `PathIntervalTooSmall` を返す
  （返ること自体がハングしない証拠）。payload の `interval_seconds`・`estimated_samples` を検証。
- 上限ちょうど（`span / interval == MAX_PATH_SAMPLES`）は**成功**（厳密 `>` の固定）。
- `interval` が NaN で `PathIntervalTooSmall`（ハングしない）。
- `interval` が非正（0・負）は**従来どおり成功**（始点のみ・エラーにしない＝既定義挙動の回帰）。
- 既定 `PathOptions::default()` は成功（上限が通常利用を制限しない）。
- 中心食でない・`include_limits=false` で走査区間が無い場合に検査が発火しない。

## 依存

ISSUE-045（path 本体）、ISSUE-048（CLI・上限定数の参照元を切り替える）。

## 完了状況（2026-09-14・strict）

**完了**。`EclipseEngine::path` が**必ず有限時間で返る**ようになった。

### 実装

- `umbra-eclipse` に `MAX_PATH_SAMPLES: f64 = 100_000.0` と
  `EclipseError::PathIntervalTooSmall { interval_seconds, estimated_samples }` を追加。
- `check_path_sample_limit(eclipse, options)` を `path()` 冒頭（走査前）で呼ぶ。
  検査対象の span は `path()` が実際に走査する最長区間＝
  中心食なら U1〜U4、部分食域を組むなら P1〜P4、両方なら長い方、どちらも無ければ**検査しない**。
- **NaN は明示的に弾く**（比較が常に false で素通りし、走査側のループ条件も false のままハングするため）。
- **非正は従来どおり成功**（「始点のみ」の既定義挙動。3 つの走査ループすべてが
  `t_sec >= span || interval <= 0.0` で 1 回で抜けることを実装レビューが確認）。
- **CLI（ISSUE-048）は数値を重複定義せず** `pub use umbra_eclipse::MAX_PATH_SAMPLES;` で再輸出。
  CLI 側の検査は「より親切に早期に弾く」ためのもので、**安全性の正本はエンジン側**。

### 実装レビュー（別サブエージェント）の指摘

- **critical 1 件**: 新関数を `trace_central` の直前に挿入したため、`trace_central` に付いていた
  `#[allow(clippy::type_complexity)]` が**新関数に付け替わり**、`trace_central` が無防備になった
  （`clippy -D warnings` ゲートを壊す）。属性を `trace_central` へ戻して修正。
- guard の網羅性（到達しうる 3 つの走査ループすべてが検査対象 span に収まること）、
  ±inf・subnormal・span=0 の挙動、CLI の定数参照は**いずれも問題なし**と確認された。

### テスト設計での発見（記録）

境界テスト（`span / interval == MAX_PATH_SAMPLES` ちょうど）は、当初 fixture の span = 7200 s では
**構成不能**だった。除算は `interval` に対し単調なので、商がちょうど `100000.0` に丸まる `interval` の帯は
**連続した 1 本**であり、span = 7200 ではその帯が**空**（探索窓を広げても無意味）。
span = 5400 s の fixture（U1〜U4 = ±0.75 h）に変えて `5400 / 0.054 == 100000.0` を厳密に満たす構成へ是正した。
「境界を近似で緩める」のではなく**fixture を変える**のが正しい対処である、という轍として記録する。

### 検証

`crates/umbra-eclipse/tests/path_limits.rs` に **9 テスト**（同ファイル計 37・全通過）:
極小正 `interval` で即エラー（**返ること自体がハングしない証拠**）／上限ちょうどは成功・1 ULP 下は失敗
（厳密 `>` の固定）／NaN はエラー／非正（0・負）は従来どおり成功で点数 1／既定は成功／
走査区間が無ければ検査しない／**P1〜P4 が使われること**（U1〜U4 は通るが P1〜P4 で落ちる構成）／
`include_limits=false` でも中心線の走査は検査される。

fmt / clippy `-D warnings` / 全テスト通過（Docker 内）。
