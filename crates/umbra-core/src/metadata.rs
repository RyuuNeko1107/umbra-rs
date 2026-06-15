//! データセットの出所・完全性メタデータ（`docs/architecture.md` §11 / `docs/data-sources.md` §5）。
//!
//! 係数・EOP・閏秒・ΔT など、外部一次データから生成・取り込みした全データセットに
//! 付随する provenance（出典・取得年代）と完全性（checksum）の共有型。xtask（生成側,
//! ISSUE-033/034/040/046）とライブラリ本体（`TimeData` 等, ISSUE-042）で共有する純粋型
//! （確定 B3）。serde 付与は公開シリアライズ API を担う ISSUE-042 で行う。

/// 全データセット共通の出所・完全性メタデータ。
///
/// 全フィールドは非空であることが provenance の最低要件（`docs/data-sources.md` §0/§5）。
/// `license` 欄は GPL 派生物を含まないこと（§0）。`checksum` は生成物バイト列の決定的ハッシュ。
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DataSetMetadata {
    /// データセット名（例: `"nutation-iau2000a"`）。
    pub name: String,
    /// 版識別子（例: `"IAU2000A"` / `"VSOP87D"`）。
    pub version: String,
    /// 一次配布元・論文引用（`docs/data-sources.md` §2/§3）。
    pub source: String,
    /// ライセンス区分（GPL 派生物は不可・§0）。
    pub license: String,
    /// 対応年代の下端。
    pub valid_from: String,
    /// 対応年代の上端。
    pub valid_to: String,
    /// 生成物バイト列の決定的ハッシュ（16 進小文字）。
    pub checksum: String,
}

impl DataSetMetadata {
    /// provenance の最低要件: 全フィールドが非空であること（`docs/data-sources.md` §0/§5）。
    pub fn has_complete_provenance(&self) -> bool {
        [
            &self.name,
            &self.version,
            &self.source,
            &self.license,
            &self.valid_from,
            &self.valid_to,
            &self.checksum,
        ]
        .iter()
        .all(|field| !field.is_empty())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// 全フィールド非空の代表メタデータ（章動 IAU2000A を模した完全 provenance）。
    fn complete() -> DataSetMetadata {
        DataSetMetadata {
            name: "nutation-iau2000a".to_string(),
            version: "IAU2000A".to_string(),
            source: "IERS Conventions 2010 / nut00a".to_string(),
            license: "public-domain".to_string(),
            valid_from: "1900-01-01".to_string(),
            valid_to: "2100-01-01".to_string(),
            checksum: "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
                .to_string(),
        }
    }

    /// 全フィールド非空なら provenance 完全（最低要件を満たす）。
    #[test]
    fn complete_metadata_has_provenance() {
        assert!(complete().has_complete_provenance());
    }

    /// どのフィールドであっても 1 つでも空文字列なら provenance 不完全。
    /// 各フィールドを順に空にしてループ検証し、欠落検出が特定フィールドに
    /// 偏っていない（全フィールドを検査している）ことを保証する。
    #[test]
    fn any_empty_field_breaks_provenance() {
        // (フィールド名, 当該フィールドのみ空にする変異関数) の列。
        type FieldMutator = fn(&mut DataSetMetadata);
        let mutators: &[(&str, FieldMutator)] = &[
            ("name", |m| m.name = String::new()),
            ("version", |m| m.version = String::new()),
            ("source", |m| m.source = String::new()),
            ("license", |m| m.license = String::new()),
            ("valid_from", |m| m.valid_from = String::new()),
            ("valid_to", |m| m.valid_to = String::new()),
            ("checksum", |m| m.checksum = String::new()),
        ];
        for (field, mutate) in mutators {
            let mut m = complete();
            mutate(&mut m);
            assert!(
                !m.has_complete_provenance(),
                "empty `{field}` should make provenance incomplete"
            );
        }
    }
}
