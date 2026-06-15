//! 係数データセットの識別（ISSUE-046 骨子）。
//!
//! 各データセットの **生成ロジック本体**は別 Issue（VSOP87=033 / ELP-MPP02=034 /
//! 章動 IAU2000A=040）。本モジュールはサブコマンド分岐で使う識別子の解析のみを担い、
//! 生成は当該 Issue が実装するまで [`XtaskError::NotImplemented`] を返す。

use crate::error::XtaskError;

/// 係数データセット種別。
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Dataset {
    /// VSOP87D 太陽（地球）系列（ISSUE-033）。
    Vsop87,
    /// ELP/MPP02 月（ISSUE-034）。
    ElpMpp02,
    /// IAU2000A 章動 nut00a（ISSUE-040）。
    NutationIau2000a,
}

impl Dataset {
    /// `--dataset` 引数の文字列から解析する。未知の値は [`XtaskError::UnknownDataset`]。
    pub fn from_arg(arg: &str) -> Result<Dataset, XtaskError> {
        match arg {
            "vsop87" => Ok(Dataset::Vsop87),
            "elp-mpp02" => Ok(Dataset::ElpMpp02),
            "nutation-iau2000a" => Ok(Dataset::NutationIau2000a),
            other => Err(XtaskError::UnknownDataset(other.to_string())),
        }
    }

    /// `--dataset` の正準文字列表現（[`Dataset::from_arg`] と往復する）。
    pub fn as_arg(self) -> &'static str {
        match self {
            Dataset::Vsop87 => "vsop87",
            Dataset::ElpMpp02 => "elp-mpp02",
            Dataset::NutationIau2000a => "nutation-iau2000a",
        }
    }

    /// 全データセット（`--dataset all` 展開用）。
    pub fn all() -> &'static [Dataset] {
        &[
            Dataset::Vsop87,
            Dataset::ElpMpp02,
            Dataset::NutationIau2000a,
        ]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// 3 種の正準引数文字列がそれぞれ対応する Dataset に解析される。
    #[test]
    fn from_arg_parses_canonical_strings() {
        assert_eq!(Dataset::from_arg("vsop87").unwrap(), Dataset::Vsop87);
        assert_eq!(Dataset::from_arg("elp-mpp02").unwrap(), Dataset::ElpMpp02);
        assert_eq!(
            Dataset::from_arg("nutation-iau2000a").unwrap(),
            Dataset::NutationIau2000a
        );
    }

    /// as_arg → from_arg の往復で同一 Dataset に戻る（正準文字列の一貫性）。
    #[test]
    fn as_arg_round_trips_through_from_arg() {
        for &d in Dataset::all() {
            assert_eq!(
                Dataset::from_arg(d.as_arg()).unwrap(),
                d,
                "round-trip failed for {d:?} via {:?}",
                d.as_arg()
            );
        }
    }

    /// 未知の文字列は UnknownDataset。
    #[test]
    fn from_arg_rejects_unknown() {
        let err = Dataset::from_arg("bogus").expect_err("unknown must error");
        assert!(
            matches!(err, XtaskError::UnknownDataset(_)),
            "expected UnknownDataset, got {err:?}"
        );
    }

    /// all() は 3 種をちょうど重複なく列挙する。
    #[test]
    fn all_contains_three_distinct_datasets() {
        let all = Dataset::all();
        assert_eq!(all.len(), 3, "exactly three datasets");
        assert!(all.contains(&Dataset::Vsop87));
        assert!(all.contains(&Dataset::ElpMpp02));
        assert!(all.contains(&Dataset::NutationIau2000a));
        // 重複なし: 各ペアが相異なる。
        for i in 0..all.len() {
            for j in (i + 1)..all.len() {
                assert_ne!(all[i], all[j], "duplicate dataset at {i},{j}");
            }
        }
    }
}
