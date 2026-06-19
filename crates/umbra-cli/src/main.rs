//! `umbra` CLI バイナリ（ISSUE-031）。引数解釈は [`umbra_cli`]、本体は薄いディスパッチ。
//!
//! サブコマンド `search`（S31a 実装済）。local / path / bessel / inspect / validate は後続 issue。

use std::process::ExitCode;

use clap::Parser;
use umbra_cli::{Cli, Command};

fn main() -> ExitCode {
    let cli = Cli::parse();
    match cli.command {
        Command::Search(args) => match umbra_cli::run_search(&args) {
            Ok(output) => {
                print!("{output}");
                ExitCode::SUCCESS
            }
            Err(error) => {
                eprintln!("error: {error}");
                ExitCode::FAILURE
            }
        },
    }
}
