//! umbra-rs 内部ビルドタスク `xtask`（ISSUE-046）。
//!
//! 係数生成（ISSUE-033/034/040）・生成物検証・データ期限/ライセンス検査のサブコマンド分岐と、
//! それらが共有する決定的プリミティブ（checksum・packed 直列化・データセット識別）を提供する。
//! 生成ロジック本体は各 Issue が `Dataset` 配下に実装する（本 crate は骨子）。
//!
//! `docs/architecture.md` §11 / `docs/data-sources.md` §5/§6 / `docs/accuracy.md` §5。

pub mod checksum;
pub mod dataset;
pub mod differential;
pub mod elp;
pub mod eop;
pub mod error;
pub mod nutation;
pub mod packed;
pub mod spk;
pub mod validate;
pub mod vsop87;

use dataset::Dataset;
use error::XtaskError;

const USAGE: &str = "\
umbra-rs xtask — internal build tasks (ISSUE-046)

USAGE:
    cargo xtask <SUBCOMMAND> [--dataset <vsop87|elp2000-82b|elp-mpp02|nutation-iau2000a|eop-c04|all>]

SUBCOMMANDS:
    generate-coefficients   一次データ→packed 係数を決定的に生成（033/034/040）
    verify-generated        コミット済み generated と再生成の checksum 差分を検査
    verify-data             EOP/閏秒/ΔT の valid_to 期限・checksum を検査
    check-licenses          cargo-deny + データ provenance/NOTICE 整合を検査
    validate                ゴールデン照合を実エンジンで実走し誤差レポートを出力（ISSUE-030）
                            [--accuracy <standard|reference>] [--format <text|json>]
    differential            解析暦×JPL DE の 2 エンジンで誤差を層分解（暦層/幾何層）出力（ISSUE-030・SLOW）
                            [--accuracy <standard|reference>] [--format <text|json>]
    fetch-de440s            JPL DE440s SPK を NAIF から data/spk/ へ取得し SHA-256 照合（ISSUE-036・非同梱）
    verify-de440s           取得済み data/spk/de440s.bsp の SHA-256 整合を検査（DL 不要）";

/// 記録済み checksum と再生成物の checksum を比較し、不一致なら [`XtaskError::ChecksumMismatch`]。
/// `verify-generated` の中核（1 バイトの差異でも fail する決定的検査）。
pub fn compare_checksum(
    dataset: &str,
    stored: &str,
    regenerated_bytes: &[u8],
) -> Result<(), XtaskError> {
    let regenerated = checksum::sha256_hex(regenerated_bytes);
    if regenerated == stored {
        Ok(())
    } else {
        Err(XtaskError::ChecksumMismatch {
            dataset: dataset.to_string(),
            stored: stored.to_string(),
            regenerated,
        })
    }
}

/// `--dataset <value>` の値を取り出すヘルパ。指定が無い／`all` は `None`（=全データセット意）、
/// 値欠落は [`XtaskError::MissingArgument`]、未知値は [`XtaskError::UnknownDataset`]。
pub fn dataset_arg(args: &[String]) -> Result<Option<Dataset>, XtaskError> {
    match args.iter().position(|arg| arg == "--dataset") {
        None => Ok(None),
        Some(flag_index) => {
            let value = args
                .get(flag_index + 1)
                .ok_or_else(|| XtaskError::MissingArgument("--dataset".to_string()))?;
            if value == "all" {
                Ok(None)
            } else {
                Ok(Some(Dataset::from_arg(value)?))
            }
        }
    }
}

/// サブコマンドを解釈して実行する（プログラム名を除いた引数列）。
/// 未知のサブコマンド/データセットは非 Ok を返し、`main` が非 0 終了に対応づける。
pub fn run(args: &[String]) -> Result<(), XtaskError> {
    let subcommand = match args.first() {
        None => {
            println!("{USAGE}");
            return Ok(());
        }
        Some(first) if first == "--help" || first == "-h" => {
            println!("{USAGE}");
            return Ok(());
        }
        Some(first) => first.as_str(),
    };

    match subcommand {
        "generate-coefficients" => match dataset_arg(args)? {
            // 章動（ISSUE-040）・VSOP87D（ISSUE-033）は実装済み。ELP/MPP02 は ISSUE-034。
            Some(Dataset::NutationIau2000a) => {
                let m = nutation::generate_to_disk()?;
                println!("generated nutation-iau2000a (checksum {})", m.checksum);
                Ok(())
            }
            Some(Dataset::Vsop87) => {
                let m = vsop87::generate_to_disk()?;
                println!("generated vsop87 (checksum {})", m.checksum);
                Ok(())
            }
            Some(Dataset::Elp200082b) => {
                let m = elp::generate_to_disk()?;
                println!("generated elp2000-82b (checksum {})", m.checksum);
                Ok(())
            }
            Some(Dataset::ElpMpp02) => Err(XtaskError::NotImplemented(
                "generate-coefficients (dataset: elp-mpp02) — 将来の MPP02 アップグレード候補"
                    .to_string(),
            )),
            Some(Dataset::EopC04) => {
                let m = eop::generate_to_disk()?;
                println!("generated eop-c04 (checksum {})", m.checksum);
                Ok(())
            }
            None => {
                // all（実装済みデータセット: nutation + vsop87 + elp2000-82b + eop-c04）。
                let n = nutation::generate_to_disk()?;
                let v = vsop87::generate_to_disk()?;
                let e = elp::generate_to_disk()?;
                let o = eop::generate_to_disk()?;
                println!(
                    "generated nutation-iau2000a ({}) + vsop87 ({}) + elp2000-82b ({}) + eop-c04 ({})",
                    n.checksum, v.checksum, e.checksum, o.checksum
                );
                Ok(())
            }
        },
        "verify-generated" => match dataset_arg(args)? {
            Some(Dataset::NutationIau2000a) => nutation::verify_against_disk(),
            Some(Dataset::Vsop87) => vsop87::verify_against_disk(),
            Some(Dataset::Elp200082b) => elp::verify_against_disk(),
            Some(Dataset::ElpMpp02) => Err(XtaskError::NotImplemented(
                "verify-generated (dataset: elp-mpp02) — 将来の MPP02 アップグレード候補"
                    .to_string(),
            )),
            Some(Dataset::EopC04) => eop::verify_against_disk(),
            None => {
                nutation::verify_against_disk()?;
                vsop87::verify_against_disk()?;
                elp::verify_against_disk()?;
                eop::verify_against_disk()
            }
        },
        "verify-data" => Err(XtaskError::NotImplemented(
            "verify-data — EOP/閏秒/ΔT 期限検査は ISSUE-007 データ管理で実装".to_string(),
        )),
        "check-licenses" => Err(XtaskError::NotImplemented(
            "check-licenses — cargo-deny + provenance 整合は CI ゲートで実装".to_string(),
        )),
        "validate" => validate::run_validate(args),
        "differential" => differential::run_differential(args),
        "fetch-de440s" => spk::fetch_de440s(),
        "verify-de440s" => spk::verify_de440s(),
        other => Err(XtaskError::UnknownSubcommand(other.to_string())),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::checksum::sha256_hex;
    use crate::packed::pack_f64_le;

    /// &str スライス → Vec<String>（run/dataset_arg の引数列構築ヘルパ）。
    fn args(parts: &[&str]) -> Vec<String> {
        parts.iter().map(|s| s.to_string()).collect()
    }

    // ------------------------------------------------------------------
    // compare_checksum: 1 バイト改変で必ず fail する決定的差分検出。
    // ------------------------------------------------------------------

    /// 記録済み checksum が再生成バイト列の SHA-256 と一致すれば Ok。
    #[test]
    fn compare_checksum_accepts_matching() {
        let bytes = pack_f64_le_bytes();
        let stored = sha256_hex(&bytes);
        assert!(compare_checksum("nutation-iau2000a", &stored, &bytes).is_ok());
    }

    /// バイト列を 1 バイト改変すると ChecksumMismatch（acceptance: 1 バイトで非 0）。
    /// dataset 名・stored・regenerated の各フィールドが正しく埋まることも確認。
    #[test]
    fn compare_checksum_detects_single_byte_change() {
        let bytes = pack_f64_le_bytes();
        let stored = sha256_hex(&bytes); // 正しい（改変前の）checksum
        let mut tampered = bytes.clone();
        tampered[0] ^= 0x01; // 先頭 1 バイトのみ反転

        let err = compare_checksum("nutation-iau2000a", &stored, &tampered)
            .expect_err("1-byte change must mismatch");
        match err {
            XtaskError::ChecksumMismatch {
                dataset,
                stored: s,
                regenerated,
            } => {
                assert_eq!(dataset, "nutation-iau2000a", "dataset name propagated");
                assert_eq!(s, stored, "stored checksum propagated verbatim");
                assert_eq!(
                    regenerated,
                    sha256_hex(&tampered),
                    "regenerated = sha256 of the (tampered) bytes"
                );
                assert_ne!(s, regenerated, "stored and regenerated must differ");
            }
            other => panic!("expected ChecksumMismatch, got {other:?}"),
        }
    }

    /// stored 側を改変しても（バイト列は正規でも）ChecksumMismatch になる。
    #[test]
    fn compare_checksum_detects_wrong_stored() {
        let bytes = pack_f64_le_bytes();
        let wrong_stored = "0000000000000000000000000000000000000000000000000000000000000000";
        let err = compare_checksum("vsop87", wrong_stored, &bytes)
            .expect_err("wrong stored must mismatch");
        assert!(
            matches!(err, XtaskError::ChecksumMismatch { .. }),
            "expected ChecksumMismatch, got {err:?}"
        );
    }

    /// 代表バイト列（packed f64 経由・非自明な内容）。
    fn pack_f64_le_bytes() -> Vec<u8> {
        pack_f64_le(&[1.0, -2.5, 0.0, f64::MIN_POSITIVE, 1e300])
    }

    // ------------------------------------------------------------------
    // run: サブコマンドディスパッチ。
    // ------------------------------------------------------------------

    /// 引数なしと --help は使用法表示で Ok(())。
    #[test]
    fn run_no_args_and_help_are_ok() {
        assert!(run(&[]).is_ok(), "no args prints usage and succeeds");
        assert!(run(&args(&["--help"])).is_ok(), "--help succeeds");
        assert!(run(&args(&["-h"])).is_ok(), "-h succeeds");
    }

    /// generate-coefficients は既知データセットでも生成本体は別 Issue → NotImplemented。
    #[test]
    fn run_generate_coefficients_known_dataset_is_not_implemented() {
        // elp-mpp02 は未実装（ISSUE-034）。nutation/vsop87 は実装済みのため別データセットで検証。
        let err = run(&args(&["generate-coefficients", "--dataset", "elp-mpp02"]))
            .expect_err("elp-mpp02 generation lives in ISSUE-034");
        let message = err.to_string();
        assert!(
            matches!(err, XtaskError::NotImplemented(_)),
            "expected NotImplemented, got {err:?}"
        );
        // メッセージは対象データセット名を含む。
        assert!(
            message.contains("elp-mpp02"),
            "NotImplemented message should name the dataset: {message}"
        );
    }

    /// verify-generated（未実装データセット elp-mpp02）は NotImplemented。
    #[test]
    fn run_verify_generated_is_not_implemented() {
        let err = run(&args(&["verify-generated", "--dataset", "elp-mpp02"]))
            .expect_err("verify-generated body not yet implemented");
        assert!(
            matches!(err, XtaskError::NotImplemented(_)),
            "expected NotImplemented, got {err:?}"
        );
    }

    /// verify-data（--dataset 不要）は NotImplemented。
    #[test]
    fn run_verify_data_is_not_implemented() {
        let err = run(&args(&["verify-data"])).expect_err("verify-data not yet implemented");
        assert!(
            matches!(err, XtaskError::NotImplemented(_)),
            "expected NotImplemented, got {err:?}"
        );
    }

    /// check-licenses（--dataset 不要）は NotImplemented。
    #[test]
    fn run_check_licenses_is_not_implemented() {
        let err = run(&args(&["check-licenses"])).expect_err("check-licenses not yet implemented");
        assert!(
            matches!(err, XtaskError::NotImplemented(_)),
            "expected NotImplemented, got {err:?}"
        );
    }

    /// 未知サブコマンドは UnknownSubcommand。
    #[test]
    fn run_unknown_subcommand() {
        let err = run(&args(&["bogus-cmd"])).expect_err("unknown subcommand must error");
        assert!(
            matches!(err, XtaskError::UnknownSubcommand(_)),
            "expected UnknownSubcommand, got {err:?}"
        );
    }

    /// 既知サブコマンド + 未知 --dataset 値は UnknownDataset。
    #[test]
    fn run_known_subcommand_unknown_dataset() {
        let err = run(&args(&["generate-coefficients", "--dataset", "bogus"]))
            .expect_err("unknown dataset must error");
        assert!(
            matches!(err, XtaskError::UnknownDataset(_)),
            "expected UnknownDataset, got {err:?}"
        );
    }

    // ------------------------------------------------------------------
    // dataset_arg: --dataset <value> の解析。
    // ------------------------------------------------------------------

    /// --dataset vsop87 → Some(Vsop87)。
    #[test]
    fn dataset_arg_parses_specific() {
        assert_eq!(
            dataset_arg(&args(&["--dataset", "vsop87"])).unwrap(),
            Some(Dataset::Vsop87)
        );
    }

    /// --dataset all → None（全データセット意）。
    #[test]
    fn dataset_arg_all_means_none() {
        assert_eq!(dataset_arg(&args(&["--dataset", "all"])).unwrap(), None);
    }

    /// --dataset 未指定 → None。
    #[test]
    fn dataset_arg_absent_means_none() {
        assert_eq!(dataset_arg(&args(&["verify-data"])).unwrap(), None);
        assert_eq!(dataset_arg(&[]).unwrap(), None);
    }

    /// --dataset が引数列末尾で値が無ければ MissingArgument。
    #[test]
    fn dataset_arg_missing_value() {
        let err = dataset_arg(&args(&["generate-coefficients", "--dataset"]))
            .expect_err("missing value must error");
        assert!(
            matches!(err, XtaskError::MissingArgument(_)),
            "expected MissingArgument, got {err:?}"
        );
    }

    /// --dataset bogus → UnknownDataset。
    #[test]
    fn dataset_arg_unknown_value() {
        let err =
            dataset_arg(&args(&["--dataset", "bogus"])).expect_err("unknown value must error");
        assert!(
            matches!(err, XtaskError::UnknownDataset(_)),
            "expected UnknownDataset, got {err:?}"
        );
    }
}
