//! `xtask` バイナリの統合テスト（ISSUE-046）。
//!
//! `main` の終了コード対応づけ（`Ok → SUCCESS` / `Err → FAILURE`）を、実際にコンパイル済み
//! バイナリを起動して検証する。`run` のディスパッチ論理は lib 単体テストで網羅済みだが、
//! `main` 自体（引数収集と `ExitCode` 写像）はここでのみ実行経路に乗る。

use std::path::PathBuf;
use std::process::Command;

/// コンパイル済み `xtask` バイナリ（cargo が `CARGO_BIN_EXE_<name>` を提供）。
fn xtask() -> Command {
    Command::new(env!("CARGO_BIN_EXE_xtask"))
}

/// リポジトリルート（`crates/xtask` の 2 階層上）。`generate`/`verify` は
/// `data/`・`generated/` をルート相対で参照するため、ここを cwd にして実行する。
fn repo_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..")
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

/// 未実装データセット（vsop87）の生成 → 非 0 終了（NotImplemented）。
#[test]
fn not_implemented_dataset_exits_nonzero() {
    let output = xtask()
        .args(["generate-coefficients", "--dataset", "vsop87"])
        .output()
        .expect("xtask binary runs");
    assert!(
        !output.status.success(),
        "not-yet-implemented dataset must exit non-zero"
    );
}

/// コミット済み章動係数が一次原データから決定的に再生成できる（生成→packed→checksum 経路の
/// end-to-end 回帰）。リポジトリルートを cwd にして実 `data/`・`generated/` を参照する。
#[test]
fn verify_generated_nutation_succeeds() {
    let output = xtask()
        .current_dir(repo_root())
        .args(["verify-generated", "--dataset", "nutation-iau2000a"])
        .output()
        .expect("xtask binary runs");
    assert!(
        output.status.success(),
        "committed nutation artifact must verify against source; stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
}

/// 原データが見つからない作業ディレクトリでは verify-generated は **失敗**する
/// （`verify_against_disk` が無条件 Ok を返さないことの保証＝検証経路が実在する）。
#[test]
fn verify_generated_fails_without_source() {
    let empty = std::env::temp_dir().join(format!("umbra_xtask_no_source_{}", std::process::id()));
    std::fs::create_dir_all(&empty).expect("create temp cwd");
    let output = xtask()
        .current_dir(&empty)
        .args(["verify-generated", "--dataset", "nutation-iau2000a"])
        .output()
        .expect("xtask binary runs");
    let _ = std::fs::remove_dir_all(&empty);
    assert!(
        !output.status.success(),
        "verify must fail when source data is absent (not unconditionally Ok)"
    );
}
