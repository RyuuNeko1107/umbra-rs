//! `xtask` バイナリ入口（ISSUE-046）。ロジックは `xtask` ライブラリ（`run`）に集約し、
//! ここは引数収集と終了コード対応づけのみ。

use std::process::ExitCode;

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    match xtask::run(&args) {
        Ok(()) => ExitCode::SUCCESS,
        Err(err) => {
            eprintln!("xtask: {err}");
            ExitCode::FAILURE
        }
    }
}
