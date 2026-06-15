//! `xtask` バイナリの統合テスト（ISSUE-046）。
//!
//! `main` の終了コード対応づけ（`Ok → SUCCESS` / `Err → FAILURE`）を、実際にコンパイル済み
//! バイナリを起動して検証する。`run` のディスパッチ論理は lib 単体テストで網羅済みだが、
//! `main` 自体（引数収集と `ExitCode` 写像）はここでのみ実行経路に乗る。

use std::process::Command;

/// コンパイル済み `xtask` バイナリ（cargo が `CARGO_BIN_EXE_<name>` を提供）。
fn xtask() -> Command {
    Command::new(env!("CARGO_BIN_EXE_xtask"))
}

/// 未知サブコマンドは非 0 終了し、診断を stderr へ出す。
#[test]
fn unknown_subcommand_exits_nonzero_with_stderr() {
    let output = xtask()
        .arg("bogus-cmd")
        .output()
        .expect("xtask binary runs");
    assert!(
        !output.status.success(),
        "unknown subcommand must exit non-zero"
    );
    assert!(
        !output.stderr.is_empty(),
        "error path should write a diagnostic to stderr"
    );
}

/// `--help` は 0 終了し、使用法を stdout へ出す。
#[test]
fn help_exits_zero_with_usage() {
    let output = xtask().arg("--help").output().expect("xtask binary runs");
    assert!(output.status.success(), "--help must exit zero");
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(
        stdout.contains("USAGE") || stdout.contains("xtask"),
        "usage text should be printed: {stdout}"
    );
}

/// 既知サブコマンドだが生成本体未実装 → 非 0 終了（NotImplemented）。
#[test]
fn not_implemented_subcommand_exits_nonzero() {
    let output = xtask()
        .arg("verify-data")
        .output()
        .expect("xtask binary runs");
    assert!(
        !output.status.success(),
        "not-yet-implemented subcommand must exit non-zero"
    );
}
