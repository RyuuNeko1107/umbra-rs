//! IAU 2000_R06 章動係数の取り込み（ISSUE-040）。
//!
//! IERS Conventions (2010) Chapter 5 の章動表 `tab5.3a.txt`（黄経 Δψ）/ `tab5.3b.txt`
//! （黄道傾斜 Δε）をパースし、packed 形式（消費側 ISSUE-035 と byte-for-byte 契約）へ
//! 直列化する。フォーマットの正本は `data/coefficient-source/nutation/PROVENANCE.md`。
//!
//! 各データ行 = `i  <sin 係数>  <cos 係数>  l l' F D Ω L_Me L_Ve L_E L_Ma L_J L_Sa L_U L_Ne p_A`。
//! 位置規約（両表共通）: 第2列 = sin の係数、第3列 = cos の係数（単位 µas）。
//! ブロック `j = 0`（定数項）/ `j = 1`（×t 項）。評価式は ISSUE-035。

use crate::checksum::sha256_hex;
use crate::error::XtaskError;
use crate::packed::{pack_f64_le, unpack_f64_le};
use std::fs;
use umbra_core::DataSetMetadata;

/// 原データ配置（`cargo xtask` をリポジトリルートから実行する前提。`docs/accuracy.md` §5）。
const SOURCE_DIR: &str = "data/coefficient-source/nutation";
/// 生成物の出力先。
const GENERATED_DIR: &str = "generated/nutation";
/// packed 係数ファイル名。
const PACKED_NAME: &str = "nut_iau2000_r06.bin";

/// 基本引数乗数の個数（luni-solar 5 + planetary 9）。
pub const N_ARGUMENTS: usize = 14;

/// 各ブロックの項数（IAU 2000_R06, IERS Conventions 2010 tab5.3a/b の宣言値）。
pub const N_PSI_CONSTANT: usize = 1320;
/// Δψ の ×t 項数。
pub const N_PSI_RATE: usize = 38;
/// Δε の定数項数。
pub const N_EPS_CONSTANT: usize = 1037;
/// Δε の ×t 項数。
pub const N_EPS_RATE: usize = 19;

/// 章動級数 1 項: 14 基本引数の整数乗数 + sin/cos 振幅（µas）。
#[derive(Clone, Debug, PartialEq)]
pub struct NutationTerm {
    /// 基本引数の整数乗数 `[l, l', F, D, Ω, L_Me, L_Ve, L_E, L_Ma, L_J, L_Sa, L_U, L_Ne, p_A]`。
    pub multipliers: [i32; N_ARGUMENTS],
    /// sin(ARG) の係数（µas）。
    pub sin_amp_uas: f64,
    /// cos(ARG) の係数（µas）。
    pub cos_amp_uas: f64,
}

/// IAU 2000_R06 章動モデル全体（4 ブロック）。
#[derive(Clone, Debug, PartialEq)]
pub struct NutationModel {
    /// Δψ 定数項（j=0）。
    pub psi_constant: Vec<NutationTerm>,
    /// Δψ ×t 項（j=1）。
    pub psi_rate: Vec<NutationTerm>,
    /// Δε 定数項（j=0）。
    pub eps_constant: Vec<NutationTerm>,
    /// Δε ×t 項（j=1）。
    pub eps_rate: Vec<NutationTerm>,
}

/// packed 1 項あたりの f64 数（14 乗数 + sin + cos）。
const TERM_STRIDE: usize = N_ARGUMENTS + 2;
/// 健全性のための項数上限（実モデルは数千。これを超える packed ヘッダは破損とみなす）。
const MAX_BLOCK_COUNT: f64 = 1_000_000.0;

/// ブロック見出し `j = K  Number of terms = N` から宣言項数 N を取り出す。
/// それ以外の行は `None`（`tab5.3b` の "Number  of terms" 二重空白も許容）。
fn parse_declared_count(line: &str) -> Option<usize> {
    let trimmed = line.trim_start();
    if !trimmed.starts_with('j') || !trimmed.contains("terms") {
        return None;
    }
    // 最後の '=' の後ろが宣言項数（`j = K ... terms = N` の N）。
    trimmed.rsplit('=').next()?.trim().parse::<usize>().ok()
}

/// データ行 `i  <sin>  <cos>  m0 .. m13` を 1 項にパースする。
/// 罫線・列見出し・空行は（先頭が整数 index でない・列数不足で）`None`。
fn parse_data_row(line: &str) -> Option<NutationTerm> {
    let tokens: Vec<&str> = line.split_whitespace().collect();
    if tokens.len() < 3 + N_ARGUMENTS {
        return None;
    }
    // 先頭は項番号（整数）。列見出し "i A_i ..." はここで弾かれる。
    tokens[0].parse::<i64>().ok()?;
    let sin_amp_uas = tokens[1].parse::<f64>().ok()?;
    let cos_amp_uas = tokens[2].parse::<f64>().ok()?;
    let mut multipliers = [0i32; N_ARGUMENTS];
    for (k, slot) in multipliers.iter_mut().enumerate() {
        *slot = tokens[3 + k].parse::<i32>().ok()?;
    }
    Some(NutationTerm {
        multipliers,
        sin_amp_uas,
        cos_amp_uas,
    })
}

/// tab5.3a/b 形式のテキストを `(j=0 ブロック, j=1 ブロック)` にパースする。
/// `j = N Number of terms = K` 宣言と実際の行数が食い違えば [`XtaskError::MalformedSource`]。
pub fn parse_table(text: &str) -> Result<(Vec<NutationTerm>, Vec<NutationTerm>), XtaskError> {
    let mut blocks: Vec<(usize, Vec<NutationTerm>)> = Vec::new();
    let mut current: Option<(usize, Vec<NutationTerm>)> = None;
    for line in text.lines() {
        if let Some(declared) = parse_declared_count(line) {
            if let Some(finished) = current.take() {
                blocks.push(finished);
            }
            current = Some((declared, Vec::new()));
        } else if let Some(term) = parse_data_row(line) {
            match current.as_mut() {
                Some((_, terms)) => terms.push(term),
                None => {
                    return Err(XtaskError::MalformedSource(
                        "data row encountered before any 'j = ...' block header".to_string(),
                    ))
                }
            }
        }
    }
    if let Some(finished) = current.take() {
        blocks.push(finished);
    }
    if blocks.len() != 2 {
        return Err(XtaskError::MalformedSource(format!(
            "expected 2 blocks (j=0, j=1), found {}",
            blocks.len()
        )));
    }
    for (declared, terms) in &blocks {
        if terms.len() != *declared {
            return Err(XtaskError::MalformedSource(format!(
                "block declares {} terms but parsed {}",
                declared,
                terms.len()
            )));
        }
    }
    let mut blocks = blocks.into_iter();
    let constant = blocks.next().expect("checked len == 2").1;
    let rate = blocks.next().expect("checked len == 2").1;
    Ok((constant, rate))
}

/// 黄経・黄道傾斜の 2 テキストから完全な章動モデルを構成する。
pub fn parse_model(
    longitude_text: &str,
    obliquity_text: &str,
) -> Result<NutationModel, XtaskError> {
    let (psi_constant, psi_rate) = parse_table(longitude_text)?;
    let (eps_constant, eps_rate) = parse_table(obliquity_text)?;
    Ok(NutationModel {
        psi_constant,
        psi_rate,
        eps_constant,
        eps_rate,
    })
}

/// モデルを packed 形式（flat little-endian f64）へ直列化する（決定的）。
/// レイアウト: `[n_psi0, n_psi1, n_eps0, n_eps1, <各ブロックの項>...]`、
/// 各項 = `[m0..m13 (14), sin_amp_uas, cos_amp_uas]`（16 f64）。
pub fn pack_model(model: &NutationModel) -> Vec<u8> {
    let blocks: [&Vec<NutationTerm>; 4] = [
        &model.psi_constant,
        &model.psi_rate,
        &model.eps_constant,
        &model.eps_rate,
    ];
    let total: usize = blocks.iter().map(|b| b.len()).sum();
    let mut values: Vec<f64> = Vec::with_capacity(4 + total * TERM_STRIDE);
    for block in blocks {
        // 項数は数千以下で i32 に収まり、i32→f64 は無損失。
        values.push(f64::from(
            i32::try_from(block.len()).expect("block length fits in i32"),
        ));
    }
    for block in blocks {
        for term in block {
            for &multiplier in &term.multipliers {
                values.push(f64::from(multiplier));
            }
            values.push(term.sin_amp_uas);
            values.push(term.cos_amp_uas);
        }
    }
    pack_f64_le(&values)
}

/// 検証済み f64 をブロック項数へ変換する（負・非整数・非有限・過大は破損）。
#[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
fn f64_to_count(value: f64) -> Result<usize, XtaskError> {
    if !value.is_finite() || value < 0.0 || value.fract() != 0.0 || value > MAX_BLOCK_COUNT {
        return Err(XtaskError::MalformedPacked(format!(
            "invalid block count {value}"
        )));
    }
    Ok(value as usize) // 検証済み: 0..=MAX_BLOCK_COUNT の整数
}

/// 検証済み f64 を整数乗数へ変換する（非整数・範囲外は破損）。
#[allow(clippy::cast_possible_truncation)]
fn f64_to_multiplier(value: f64) -> Result<i32, XtaskError> {
    if !value.is_finite()
        || value.fract() != 0.0
        || value < f64::from(i32::MIN)
        || value > f64::from(i32::MAX)
    {
        return Err(XtaskError::MalformedPacked(format!(
            "invalid integer multiplier {value}"
        )));
    }
    Ok(value as i32) // 検証済み: i32 範囲の整数
}

/// `values[*cursor]` を境界チェック付きで取り出し、カーソルを進める。
/// 範囲外は [`XtaskError::MalformedPacked`]（呼び出し側の不変条件が崩れても panic しない）。
fn take(values: &[f64], cursor: &mut usize) -> Result<f64, XtaskError> {
    let value = values.get(*cursor).copied().ok_or_else(|| {
        XtaskError::MalformedPacked("packed nutation data ended mid-term".to_string())
    })?;
    *cursor += 1;
    Ok(value)
}

/// `values[*cursor..]` から `n` 項を読み出し、カーソルを進める。
fn read_terms(
    values: &[f64],
    cursor: &mut usize,
    n: usize,
) -> Result<Vec<NutationTerm>, XtaskError> {
    let mut block = Vec::with_capacity(n);
    for _ in 0..n {
        let mut multipliers = [0i32; N_ARGUMENTS];
        for slot in multipliers.iter_mut() {
            *slot = f64_to_multiplier(take(values, cursor)?)?;
        }
        let sin_amp_uas = take(values, cursor)?;
        let cos_amp_uas = take(values, cursor)?;
        block.push(NutationTerm {
            multipliers,
            sin_amp_uas,
            cos_amp_uas,
        });
    }
    Ok(block)
}

/// packed バイト列を [`NutationModel`] へ復元する。
pub fn unpack_model(bytes: &[u8]) -> Result<NutationModel, XtaskError> {
    let values = unpack_f64_le(bytes)?;
    if values.len() < 4 {
        return Err(XtaskError::MalformedPacked(format!(
            "packed nutation needs at least 4 header f64, got {}",
            values.len()
        )));
    }
    let counts = [
        f64_to_count(values[0])?,
        f64_to_count(values[1])?,
        f64_to_count(values[2])?,
        f64_to_count(values[3])?,
    ];
    let total: usize = counts.iter().sum();
    let expected = 4 + total
        .checked_mul(TERM_STRIDE)
        .ok_or_else(|| XtaskError::MalformedPacked("term count overflow".to_string()))?;
    if values.len() != expected {
        return Err(XtaskError::MalformedPacked(format!(
            "packed length {} does not match header counts {:?} (expected {} f64)",
            values.len(),
            counts,
            expected
        )));
    }
    let mut cursor = 4usize;
    let psi_constant = read_terms(&values, &mut cursor, counts[0])?;
    let psi_rate = read_terms(&values, &mut cursor, counts[1])?;
    let eps_constant = read_terms(&values, &mut cursor, counts[2])?;
    let eps_rate = read_terms(&values, &mut cursor, counts[3])?;
    Ok(NutationModel {
        psi_constant,
        psi_rate,
        eps_constant,
        eps_rate,
    })
}

/// 章動原データ（2 テキスト）から packed バイト列と [`DataSetMetadata`] を構成する（純関数・決定的）。
/// `checksum` は packed バイト列の SHA-256。
pub fn build_artifact(
    longitude_text: &str,
    obliquity_text: &str,
) -> Result<(Vec<u8>, DataSetMetadata), XtaskError> {
    let model = parse_model(longitude_text, obliquity_text)?;
    let bytes = pack_model(&model);
    let checksum = sha256_hex(&bytes);
    let metadata = DataSetMetadata {
        name: "nutation-iau2000a".to_string(),
        version: "IAU 2000_R06".to_string(),
        source: "IERS Conventions (2010) Ch.5 tab5.3a (Δψ) / tab5.3b (Δε); \
                 Mathews-Herring-Buffett (2002), IAU2006 adj. Capitaine-Wallace-Chapront (2003)"
            .to_string(),
        license: "IERS published scientific data (attribution); SOFA reference-only, not ported"
            .to_string(),
        valid_from: "1900-01-01".to_string(),
        valid_to: "2100-01-01".to_string(),
        checksum,
    };
    Ok((bytes, metadata))
}

/// `SOURCE_DIR` 配下のファイルを読む（I/O エラーは [`XtaskError::Io`]）。
fn read_source(file_name: &str) -> Result<String, XtaskError> {
    let path = format!("{SOURCE_DIR}/{file_name}");
    fs::read_to_string(&path).map_err(|source| XtaskError::Io { path, source })
}

/// `metadata.txt` の人間可読レンダリング。
fn render_metadata(metadata: &DataSetMetadata) -> String {
    format!(
        "name = {}\nversion = {}\nsource = {}\nlicense = {}\nvalid_from = {}\nvalid_to = {}\nchecksum = {}\n",
        metadata.name,
        metadata.version,
        metadata.source,
        metadata.license,
        metadata.valid_from,
        metadata.valid_to,
        metadata.checksum,
    )
}

/// 一次原データを読み、packed 係数・metadata・NOTICE を `GENERATED_DIR` へ決定的に書き出す。
/// `cargo xtask generate-coefficients --dataset nutation-iau2000a`（リポジトリルートで実行）。
pub fn generate_to_disk() -> Result<DataSetMetadata, XtaskError> {
    let longitude = read_source("tab5.3a.txt")?;
    let obliquity = read_source("tab5.3b.txt")?;
    let (bytes, metadata) = build_artifact(&longitude, &obliquity)?;

    fs::create_dir_all(GENERATED_DIR).map_err(|source| XtaskError::Io {
        path: GENERATED_DIR.to_string(),
        source,
    })?;
    let write = |name: &str, content: &[u8]| -> Result<(), XtaskError> {
        let path = format!("{GENERATED_DIR}/{name}");
        fs::write(&path, content).map_err(|source| XtaskError::Io { path, source })
    };
    write(PACKED_NAME, &bytes)?;
    write("metadata.txt", render_metadata(&metadata).as_bytes())?;
    write(
        "NOTICE.md",
        format!(
            "# Generated nutation coefficients (IAU 2000_R06)\n\n\
             Generated by `cargo xtask generate-coefficients --dataset nutation-iau2000a` from\n\
             `{SOURCE_DIR}/` (see PROVENANCE.md). Do not edit by hand; regenerate and verify with\n\
             `cargo xtask verify-generated --dataset nutation-iau2000a`.\n\n\
             {}\n",
            render_metadata(&metadata)
        )
        .as_bytes(),
    )?;
    Ok(metadata)
}

/// コミット済み packed 係数が一次原データから決定的に再生成できることを検証する。
/// 差分があれば [`XtaskError::ChecksumMismatch`]（CI fail）。
pub fn verify_against_disk() -> Result<(), XtaskError> {
    let longitude = read_source("tab5.3a.txt")?;
    let obliquity = read_source("tab5.3b.txt")?;
    let (regenerated, _) = build_artifact(&longitude, &obliquity)?;
    let committed_path = format!("{GENERATED_DIR}/{PACKED_NAME}");
    let committed = fs::read(&committed_path).map_err(|source| XtaskError::Io {
        path: committed_path,
        source,
    })?;
    crate::compare_checksum("nutation-iau2000a", &sha256_hex(&committed), &regenerated)
}

#[cfg(test)]
mod tests {
    use super::*;

    // ------------------------------------------------------------------
    // 原データ（コミット済み実ファイル）。manifest 相対 include_str! で
    // cwd 非依存に読み込む（src→xtask→crates→repo root）。
    // 出典: data/coefficient-source/nutation/PROVENANCE.md（IERS Conventions
    // 2010, Chapter 5, tab5.3a Δψ / tab5.3b Δε, IAU 2000_R06）。
    // ------------------------------------------------------------------
    const TAB_5_3A: &str = include_str!("../../../data/coefficient-source/nutation/tab5.3a.txt");
    const TAB_5_3B: &str = include_str!("../../../data/coefficient-source/nutation/tab5.3b.txt");

    /// スポット値検証用ヘルパ: 1 項の multipliers / sin / cos を厳密一致で確認。
    /// 振幅は表の 10 進が f64 へ一意変換されるため厳密 `==`（リテラル比較）でよい。
    #[track_caller]
    fn assert_term(term: &NutationTerm, multipliers: [i32; N_ARGUMENTS], sin: f64, cos: f64) {
        assert_eq!(term.multipliers, multipliers, "multipliers mismatch");
        assert_eq!(term.sin_amp_uas, sin, "sin_amp_uas mismatch (col2=sin)");
        assert_eq!(term.cos_amp_uas, cos, "cos_amp_uas mismatch (col3=cos)");
    }

    // ------------------------------------------------------------------
    // (1) 項数: parse_model の 4 ブロック長が宣言値と一致。
    //     宣言: Δψ j=0=1320 / j=1=38, Δε j=0=1037 / j=1=19
    //     (tab5.3a:19,1345 / tab5.3b:19,1062 のブロック見出し)。
    // ------------------------------------------------------------------
    #[test]
    fn parse_model_block_lengths_match_declared_counts() {
        let model = parse_model(TAB_5_3A, TAB_5_3B).expect("real IERS tables parse");
        assert_eq!(model.psi_constant.len(), N_PSI_CONSTANT, "Δψ j=0 = 1320");
        assert_eq!(model.psi_rate.len(), N_PSI_RATE, "Δψ j=1 = 38");
        assert_eq!(model.eps_constant.len(), N_EPS_CONSTANT, "Δε j=0 = 1037");
        assert_eq!(model.eps_rate.len(), N_EPS_RATE, "Δε j=1 = 19");
    }

    // ------------------------------------------------------------------
    // (2) スポット値: 原データから検証済みの 4 項を厳密一致で確認。
    //     列ズレ・sin/cos 取り違え・lunisolar/planetary 取り違えを潰す。
    // ------------------------------------------------------------------

    /// Δψ 定数 第1項（tab5.3a:23）:
    ///   `1   -17206424.18   3338.60   0 0 0 0 1 0 0 0 0 0 0 0 0 0`
    #[test]
    fn psi_constant_first_term_spot_value() {
        let model = parse_model(TAB_5_3A, TAB_5_3B).expect("real IERS tables parse");
        assert_term(
            &model.psi_constant[0],
            [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            -17206424.18,
            3338.60,
        );
    }

    /// Δψ 定数 第37項（tab5.3a:59）: 惑星乗数を含む混合項。
    ///   `37   -308.40   512.30   0 0 1 -1 1 0 0 -1 0 -2 5 0 0 0`
    /// 惑星乗数 [.., L_E=-1, L_Ma=0, L_J=-2, L_Sa=5, ..] の桁・符号・列対応を検証。
    #[test]
    fn psi_constant_planetary_term_37_spot_value() {
        let model = parse_model(TAB_5_3A, TAB_5_3B).expect("real IERS tables parse");
        assert_term(
            &model.psi_constant[36],
            [0, 0, 1, -1, 1, 0, 0, -1, 0, -2, 5, 0, 0, 0],
            -308.40,
            512.30,
        );
    }

    /// Δψ ×t 第1項（tab5.3a:1349, i=1321）:
    ///   `1321   -17418.82   2.89   0 0 0 0 1 0 0 0 0 0 0 0 0 0`
    #[test]
    fn psi_rate_first_term_spot_value() {
        let model = parse_model(TAB_5_3A, TAB_5_3B).expect("real IERS tables parse");
        assert_term(
            &model.psi_rate[0],
            [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            -17418.82,
            2.89,
        );
    }

    /// Δε 定数 第1項（tab5.3b:23）: 位置規約 col2=sin / col3=cos を検証。
    ///   `1   1537.70   9205233.10   0 0 0 0 1 0 0 0 0 0 0 0 0 0`
    /// （tab5.3b 見出しは `B"_i B_i` = sin が先・cos が後 → col2=sin=1537.70,
    ///   col3=cos=9205233.10）。sin/cos 取り違えがあれば即 fail。
    #[test]
    fn eps_constant_first_term_spot_value_and_column_convention() {
        let model = parse_model(TAB_5_3A, TAB_5_3B).expect("real IERS tables parse");
        assert_term(
            &model.eps_constant[0],
            [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            1537.70,
            9205233.10,
        );
    }

    // ------------------------------------------------------------------
    // (3) parse_table 単体: 小さな合成テキスト（罫線/ヘッダ含む）で
    //     (j=0=2項, j=1=1項) を返し、各値が一致。
    // ------------------------------------------------------------------

    /// 罫線・ヘッダ・空行を含む最小の合成テーブル（宣言と実行数が一致）。
    /// データ行は `i  col2(sin)  col3(cos)  m0..m13`。
    fn synthetic_table_2_then_1() -> String {
        "\
----------------------------------------
j = 0  Number of terms = 2
----------------------------------------
    i        A_i             A\"_i     l    l'   F    D    Om  L_Me L_Ve  L_E L_Ma  L_J L_Sa  L_U L_Ne  p_A
----------------------------------------
    1       -17206424.18     3338.60    0    0    0    0    1    0    0    0    0    0    0    0    0    0
    2        -308.40         512.30     0    0    1   -1    1    0    0   -1    0   -2    5    0    0    0

----------------------------------------
j = 1  Number of terms = 1
----------------------------------------
    i        A'_i            A\"'_i     l    l'   F    D    Om  L_Me L_Ve  L_E L_Ma  L_J L_Sa  L_U L_Ne  p_A
----------------------------------------
 1321       -17418.82           2.89    0    0    0    0    1    0    0    0    0    0    0    0    0    0
"
        .to_string()
    }

    #[test]
    fn parse_table_synthetic_returns_two_then_one_block() {
        let (j0, j1) = parse_table(&synthetic_table_2_then_1()).expect("synthetic table parses");
        assert_eq!(j0.len(), 2, "j=0 block has 2 terms");
        assert_eq!(j1.len(), 1, "j=1 block has 1 term");

        assert_term(
            &j0[0],
            [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            -17206424.18,
            3338.60,
        );
        assert_term(
            &j0[1],
            [0, 0, 1, -1, 1, 0, 0, -1, 0, -2, 5, 0, 0, 0],
            -308.40,
            512.30,
        );
        assert_term(
            &j1[0],
            [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            -17418.82,
            2.89,
        );
    }

    // ------------------------------------------------------------------
    // (4) 項数不整合: 宣言 `= 3` だが 2 データ行 → MalformedSource。
    // ------------------------------------------------------------------
    #[test]
    fn parse_table_rejects_term_count_mismatch() {
        let text = "\
----------------------------------------
j = 0  Number of terms = 3
----------------------------------------
    i        A_i             A\"_i     l    l'   F    D    Om  L_Me L_Ve  L_E L_Ma  L_J L_Sa  L_U L_Ne  p_A
----------------------------------------
    1       -17206424.18     3338.60    0    0    0    0    1    0    0    0    0    0    0    0    0    0
    2        -308.40         512.30     0    0    1   -1    1    0    0   -1    0   -2    5    0    0    0

----------------------------------------
j = 1  Number of terms = 1
----------------------------------------
    i        A'_i            A\"'_i     l    l'   F    D    Om  L_Me L_Ve  L_E L_Ma  L_J L_Sa  L_U L_Ne  p_A
----------------------------------------
 1321       -17418.82           2.89    0    0    0    0    1    0    0    0    0    0    0    0    0    0
";
        let err = parse_table(text).expect_err("declared 3 but only 2 rows must error");
        assert!(
            matches!(err, XtaskError::MalformedSource(_)),
            "expected MalformedSource, got {err:?}"
        );
    }

    // ------------------------------------------------------------------
    // (5) packed 往復: parse_model → pack_model → unpack_model が元と完全一致。
    // ------------------------------------------------------------------
    #[test]
    fn pack_unpack_round_trips_full_model() {
        let model = parse_model(TAB_5_3A, TAB_5_3B).expect("real IERS tables parse");
        let restored = unpack_model(&pack_model(&model)).expect("packed model round-trips");
        assert_eq!(restored, model, "round-tripped model must equal original");
    }

    // ------------------------------------------------------------------
    // (6) pack_model 決定性・長さ。
    //     長さ = 8 * (4 + 16 * (1320+38+1037+19)) = 8 * 38628 バイト。
    // ------------------------------------------------------------------
    #[test]
    fn pack_model_is_deterministic_and_has_exact_length() {
        let model = parse_model(TAB_5_3A, TAB_5_3B).expect("real IERS tables parse");
        let a = pack_model(&model);
        let b = pack_model(&model);
        assert_eq!(a, b, "pack_model must be deterministic");

        let n_terms = N_PSI_CONSTANT + N_PSI_RATE + N_EPS_CONSTANT + N_EPS_RATE;
        let expected_f64 = 4 + 16 * n_terms; // 4 ヘッダ件数 + 各項 16 f64
        assert_eq!(
            a.len(),
            8 * expected_f64,
            "packed length = 8 * (4 + 16 * total_terms)"
        );
    }

    // ------------------------------------------------------------------
    // (7) unpack 異常: 壊れた packed バイト列で Err（パニックしない）。
    //     長さ検査 / 件数整合のどちらに委ねるかは未確定なので緩く検証。
    // ------------------------------------------------------------------
    #[test]
    fn unpack_rejects_non_f64_aligned_bytes() {
        // 8 の倍数でない長さ → 必ず Err（MalformedPacked か MalformedSource）。
        let bytes = vec![0u8; 13];
        let err = unpack_model(&bytes).expect_err("non-8-multiple must error");
        assert!(
            matches!(
                err,
                XtaskError::MalformedPacked(_) | XtaskError::MalformedSource(_)
            ),
            "expected Malformed*, got {err:?}"
        );
    }

    #[test]
    fn unpack_rejects_header_data_length_mismatch() {
        // ヘッダで巨大な件数を宣言するがデータが伴わない（8 の倍数長は満たす）。
        // [n_psi0=1e9, n_psi1=0, n_eps0=0, n_eps1=0] の 4 f64 のみ → 件数とデータ長が矛盾。
        let mut values = vec![1.0e9_f64, 0.0, 0.0, 0.0];
        let _ = &mut values; // 明示
        let bytes = pack_f64_le(&values);
        let err =
            unpack_model(&bytes).expect_err("header count vs data length mismatch must error");
        assert!(
            matches!(
                err,
                XtaskError::MalformedPacked(_) | XtaskError::MalformedSource(_)
            ),
            "expected Malformed*, got {err:?}"
        );
    }

    // ------------------------------------------------------------------
    // (8) build_artifact（純関数）: packed バイト・DataSetMetadata の整合。
    // ------------------------------------------------------------------

    /// checksum は packed バイトの SHA-256 と一致し、metadata は全フィールド非空・決定的。
    #[test]
    fn build_artifact_checksum_and_metadata_are_consistent() {
        let (bytes, metadata) = build_artifact(TAB_5_3A, TAB_5_3B).expect("artifact builds");
        // packed バイトは pack_model と同一（純関数の決定性）。
        let model = parse_model(TAB_5_3A, TAB_5_3B).unwrap();
        assert_eq!(bytes, pack_model(&model), "artifact bytes == pack_model");
        // checksum は packed バイトの SHA-256。
        assert_eq!(
            metadata.checksum,
            sha256_hex(&bytes),
            "checksum = sha256(bytes)"
        );
        // provenance 完全（全フィールド非空）。
        assert!(
            metadata.has_complete_provenance(),
            "metadata must have complete provenance: {metadata:?}"
        );
        // 版・名称は R06 章動を示す。
        assert_eq!(metadata.name, "nutation-iau2000a");
        assert_eq!(metadata.version, "IAU 2000_R06");
    }

    /// build_artifact は決定的（同一入力 → 同一バイト・同一 checksum）。
    #[test]
    fn build_artifact_is_deterministic() {
        let (bytes_a, meta_a) = build_artifact(TAB_5_3A, TAB_5_3B).unwrap();
        let (bytes_b, meta_b) = build_artifact(TAB_5_3A, TAB_5_3B).unwrap();
        assert_eq!(bytes_a, bytes_b, "bytes deterministic");
        assert_eq!(meta_a.checksum, meta_b.checksum, "checksum deterministic");
    }

    // ------------------------------------------------------------------
    // (9) 検証ヘルパの境界（f64_to_count / f64_to_multiplier / parse_declared_count）。
    //     各条件の境界（< vs <=、> vs >=、|| の各節）を独立に踏む。
    // ------------------------------------------------------------------

    #[test]
    fn f64_to_count_accepts_valid_and_boundary() {
        assert_eq!(f64_to_count(0.0).unwrap(), 0, "zero is a valid count");
        assert_eq!(f64_to_count(1320.0).unwrap(), 1320, "typical count");
        // 上限ちょうどは受理（`>` 境界、`>=` 改変で fail させる）。
        assert_eq!(
            f64_to_count(MAX_BLOCK_COUNT).unwrap(),
            1_000_000,
            "exactly MAX_BLOCK_COUNT is accepted"
        );
    }

    #[test]
    fn f64_to_count_rejects_invalid() {
        // 上限超過（`>` を `==`/`>=` に改変すると 1 ずれて生存するため境界+1 で踏む）。
        assert!(f64_to_count(MAX_BLOCK_COUNT + 1.0).is_err(), "above max");
        assert!(f64_to_count(-1.0).is_err(), "negative");
        assert!(f64_to_count(1.5).is_err(), "non-integer");
        assert!(f64_to_count(f64::NAN).is_err(), "NaN");
        assert!(f64_to_count(f64::INFINITY).is_err(), "infinity");
    }

    #[test]
    fn f64_to_multiplier_accepts_valid_and_boundary() {
        assert_eq!(f64_to_multiplier(0.0).unwrap(), 0);
        assert_eq!(f64_to_multiplier(5.0).unwrap(), 5);
        assert_eq!(f64_to_multiplier(-2.0).unwrap(), -2);
        // i32 範囲端ちょうどは受理（`<`/`>` 境界）。
        assert_eq!(f64_to_multiplier(f64::from(i32::MAX)).unwrap(), i32::MAX);
        assert_eq!(f64_to_multiplier(f64::from(i32::MIN)).unwrap(), i32::MIN);
    }

    #[test]
    fn f64_to_multiplier_rejects_invalid() {
        assert!(
            f64_to_multiplier(f64::from(i32::MAX) + 1.0).is_err(),
            "above i32::MAX"
        );
        assert!(
            f64_to_multiplier(f64::from(i32::MIN) - 1.0).is_err(),
            "below i32::MIN"
        );
        assert!(f64_to_multiplier(1.5).is_err(), "non-integer");
        assert!(f64_to_multiplier(f64::NAN).is_err(), "NaN");
        assert!(f64_to_multiplier(f64::INFINITY).is_err(), "infinity");
    }

    #[test]
    fn parse_declared_count_recognizes_only_real_headers() {
        // 実ヘッダ（単一空白・二重空白の両方）。
        assert_eq!(
            parse_declared_count("j = 0  Number of terms = 1320"),
            Some(1320)
        );
        assert_eq!(
            parse_declared_count("j = 0  Number  of terms = 1037"),
            Some(1037),
            "double-space 'Number  of terms' (tab5.3b) is accepted"
        );
        // 'j' で始まるが "terms" を含まない行は見出しでない（`||`→`&&` 改変を殺す）。
        assert_eq!(parse_declared_count("j stuff = 7"), None);
        // "terms" を含むが 'j' で始まらない行も見出しでない。
        assert_eq!(parse_declared_count("Total terms = 5"), None);
        // データ行・罫線は None。
        assert_eq!(
            parse_declared_count("    1   -17206424.18  3338.60  0"),
            None
        );
    }

    // ------------------------------------------------------------------
    // (10) unpack の空モデル・ヘッダ長境界。
    // ------------------------------------------------------------------

    /// 4 ヘッダ件数がすべて 0 の packed は「空モデル」として Ok（`len() < 4` 境界）。
    #[test]
    fn unpack_all_zero_counts_is_empty_model() {
        let bytes = pack_f64_le(&[0.0, 0.0, 0.0, 0.0]);
        let model = unpack_model(&bytes).expect("4 zero counts → empty model");
        assert!(model.psi_constant.is_empty());
        assert!(model.psi_rate.is_empty());
        assert!(model.eps_constant.is_empty());
        assert!(model.eps_rate.is_empty());
    }

    // ------------------------------------------------------------------
    // (11) Δε ×t 第1項スポット値（4 ブロック目を独立に検証）。
    //   tab5.3b j=1 第1項。round-trip では潰せない列ズレを独立確認。
    // ------------------------------------------------------------------
    #[test]
    fn eps_rate_first_term_spot_value() {
        let (_, eps_rate) = parse_table(TAB_5_3B).expect("tab5.3b parses");
        assert_eq!(eps_rate.len(), N_EPS_RATE);
        // tab5.3b j=1 第1項（tab5.3b:1066, i=1038）: `0.20  883.03  0 0 0 0 1 ...`
        // col2=sin=0.20, col3=cos=883.03, Ω 主項。
        assert_term(
            &eps_rate[0],
            [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            0.20,
            883.03,
        );
    }

    // ------------------------------------------------------------------
    // (12) render_metadata: 全フィールドを行として出力する（mutation: 本体 → "" を殺す）。
    // ------------------------------------------------------------------
    #[test]
    fn render_metadata_includes_all_fields() {
        let (_, metadata) = build_artifact(TAB_5_3A, TAB_5_3B).unwrap();
        let rendered = render_metadata(&metadata);
        assert!(rendered.contains("name = nutation-iau2000a"), "{rendered}");
        assert!(rendered.contains("version = IAU 2000_R06"), "{rendered}");
        assert!(
            rendered.contains(&format!("checksum = {}", metadata.checksum)),
            "{rendered}"
        );
    }
}
