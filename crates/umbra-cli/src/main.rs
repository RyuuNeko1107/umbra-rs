//! `umbra` CLI バイナリ（ISSUE-031/032）。引数解釈は [`umbra_cli`]、本体は薄いディスパッチ。
//!
//! サブコマンド `search`（ISSUE-031）/ `local`（ISSUE-032）。path / bessel / inspect / validate は後続 issue。

use std::process::ExitCode;

use clap::Parser;
use umbra_cli::{Cli, Command};

fn main() -> ExitCode {
    let cli = Cli::parse();
    let result = match cli.command {
        Command::Search(args) => umbra_cli::run_search(&args),
        Command::Local(args) => umbra_cli::run_local(&args),
    };
    match result {
        Ok(output) => {
            print!("{output}");
            ExitCode::SUCCESS
        }
        Err(error) => {
            eprintln!("error: {error}");
            ExitCode::FAILURE
        }
    }
}
