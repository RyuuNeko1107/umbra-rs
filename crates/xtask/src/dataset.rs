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
    /// ELP2000-82B 月（ISSUE-034、現行採用）。
    Elp200082b,
    /// ELP/MPP02 月（将来のアップグレード候補。原データ未入手・未実装）。
    ElpMpp02,
    /// IAU2000A 章動 nut00a（ISSUE-040）。
    NutationIau2000a,
    /// IERS EOP 14 C04 地球姿勢（UT1−UTC・極運動, ISSUE-007 EOP）。
    EopC04,
}

impl Dataset {
    /// `--dataset` 引数の文字列から解析する。未知の値は [`XtaskError::UnknownDataset`]。
    pub fn from_arg(arg: &str) -> Result<Dataset, XtaskError> {
        match arg {
            "vsop87" => Ok(Dataset::Vsop87),
            "elp2000-82b" => Ok(Dataset::Elp200082b),
            "elp-mpp02" => Ok(Dataset::ElpMpp02),
            "nutation-iau2000a" => Ok(Dataset::NutationIau2000a),
            "eop-c04" => Ok(Dataset::EopC04),
            other => Err(XtaskError::UnknownDataset(other.to_string())),
        }
    }

    /// `--dataset` の正準文字列表現（[`Dataset::from_arg`] と往復する）。
    pub fn as_arg(self) -> &'static str {
        match self {
            Dataset::Vsop87 => "vsop87",
            Dataset::Elp200082b => "elp2000-82b",
            Dataset::ElpMpp02 => "elp-mpp02",
            Dataset::NutationIau2000a => "nutation-iau2000a",
            Dataset::EopC04 => "eop-c04",
        }
    }

    /// 全データセット（`--dataset all` 展開用）。
    pub fn all() -> &'static [Dataset] {
        &[
            Dataset::Vsop87,
            Dataset::Elp200082b,
            Dataset::ElpMpp02,
            Dataset::NutationIau2000a,
            Dataset::EopC04,
        ]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// 正準引数文字列がそれぞれ対応する Dataset に解析される。
    #[test]
    fn from_arg_parses_canonical_strings() {
        assert_eq!(Dataset::from_arg("vsop87").unwrap(), Dataset::Vsop87);
        assert_eq!(
            Dataset::from_arg("elp2000-82b").unwrap(),
            Dataset::Elp200082b
        );
        assert_eq!(Dataset::from_arg("elp-mpp02").unwrap(), Dataset::ElpMpp02);
        assert_eq!(
            Dataset::from_arg("nutation-iau2000a").unwrap(),
            Dataset::NutationIau2000a
        );
        assert_eq!(Dataset::from_arg("eop-c04").unwrap(), Dataset::EopC04);
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

    /// all() は 5 種をちょうど重複なく列挙する。
    #[test]
    fn all_contains_distinct_datasets() {
        let all = Dataset::all();
        assert_eq!(all.len(), 5, "exactly five datasets");
        assert!(all.contains(&Dataset::Vsop87));
        assert!(all.contains(&Dataset::Elp200082b));
        assert!(all.contains(&Dataset::ElpMpp02));
        assert!(all.contains(&Dataset::NutationIau2000a));
        assert!(all.contains(&Dataset::EopC04));
        // 重複なし: 各ペアが相異なる。
        for i in 0..all.len() {
            for j in (i + 1)..all.len() {
                assert_ne!(all[i], all[j], "duplicate dataset at {i},{j}");
            }
        }
    }
}
