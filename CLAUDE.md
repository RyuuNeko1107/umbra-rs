# CLAUDE.md

`umbra-rs` — experimental pure-Rust solar eclipse prediction engine (AI-driven design,
implementation, and verification). Design docs live under `docs/`; verification runs in
Docker (`docker compose -p umbra-rs run --rm rust ...`).

## Implementation workflow (tdd-workflow.md — 遵守必須)

実装は `tdd-workflow.md` を**遵守する**。standard/strict の issue では役割をサブエージェントで分離する（追認禁止）:

1. **工程0**: モード判定（trivial / standard / strict）を最終報告に明記。
2. **テスト設計・作成**: テスト担当サブエージェントに委譲（**確定仕様・公開IF・関連既存テストのみ**渡す。実装案・予定差分は渡さない）。
3. **red 確認**: 実装前にテストが想定どおりの理由で失敗することを確認。
4. **実装**: メインエージェント。
5. **実装レビュー**: 実装者とは別のサブエージェント（**作成者の自己評価・結論を渡さない**）。
6. **全テスト実行**: Docker 内で fmt / clippy -D warnings / test。
7. **テスト結果レビュー＋ミューテーション**（strict は `cargo mutants`。生存変異は列挙・許容判断）。

trivial（typo/コメント/設定値/機械的リネーム）のみサブエージェント分離を省略可。**実験・検証は全て Docker 内**で実行する。

## Agent skills

### Issue tracker

Issues/PRDs live as GitHub issues (`RyuuNeko1107/umbra-rs`) via the `gh` CLI. See `docs/agents/issue-tracker.md`.

### Triage labels

Five canonical triage labels with default strings (`needs-triage` / `needs-info` / `ready-for-agent` / `ready-for-human` / `wontfix`). See `docs/agents/triage-labels.md`.

### Domain docs

Single-context: `CONTEXT.md` + `docs/adr/` at the repo root. See `docs/agents/domain.md`.
