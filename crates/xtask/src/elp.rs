//! ELP2000-82B 月係数の取り込み（ISSUE-034）。
//!
//! IMCCE 配布の 36 ファイル `ELP1`..`ELP36` をパースし、packed 形式（消費側 ISSUE-014 と
//! byte-for-byte 契約）へ直列化する。フォーマットの正本は
//! `data/coefficient-source/elp2000-82b/PROVENANCE.md`。
//!
//! 2 形式: 主問題（ELP1-3）= `4 整数(D,l',l,F 乗数) + 7 実数(A=coef1, 偏微分 B1..B6=coef2..7)`
//! （評価器 ISSUE-014 は coef1..6 を DE200/LE200 フィット補正に使う。coef7 は未使用だが忠実保存）;
//! 摂動（ELP4-36）= `K 整数乗数 + φ(度) + A + 周期`（period 不使用）。係数の解釈・引数定義・
//! sin/cos・×T の適用（評価式）は ISSUE-014。本 Issue は乗数群 + 係数を忠実に取り込む。

use crate::checksum::sha256_hex;
use crate::error::XtaskError;
use crate::packed::{pack_f64_le, unpack_f64_le};
use std::fs;
use umbra_core::DataSetMetadata;

/// ELP2000-82B のファイル数（ELP1..ELP36）。
pub const N_FILES: usize = 36;

/// 原データ配置（`cargo xtask` をリポジトリルートから実行する前提。`docs/accuracy.md` §5）。
const SOURCE_DIR: &str = "data/coefficient-source/elp2000-82b";
/// 生成物の出力先。
const GENERATED_DIR: &str = "generated/elp2000-82b";
/// packed 係数ファイル名。
const PACKED_NAME: &str = "elp2000-82b_moon.bin";
/// 健全性のための上限（カウント。実データ最大 14328 項/系列）。
const MAX_COUNT: f64 = 1_000_000.0;

/// 1 項: 整数乗数 + 係数。主問題は coefficients=\[A, B1..B6\]（7 個）、摂動は \[φ(度), A\]。
#[derive(Clone, Debug, PartialEq)]
pub struct Elp82bTerm {
    /// 引数の整数乗数（主問題=4, 摂動=系列ごとに異なる K）。
    pub multipliers: Vec<i32>,
    /// 係数（主問題=\[A, B1..B6\]（7 個）、摂動=\[phase_deg, A\]）。
    pub coefficients: Vec<f64>,
}

/// 1 ファイル（系列）。`file` は 1..=36（変数 L/B/R・摂動群・×T・引数集合は ISSUE-014 が file から導く）。
#[derive(Clone, Debug, PartialEq)]
pub struct Elp82bSeries {
    /// ファイル番号（1..=36）。
    pub file: u8,
    /// 項列。
    pub terms: Vec<Elp82bTerm>,
}

/// ELP2000-82B モデル全体（36 系列、ファイル順）。
#[derive(Clone, Debug, PartialEq)]
pub struct Elp82bModel {
    /// 系列列（ELP1..ELP36 の順）。
    pub series: Vec<Elp82bSeries>,
}

/// 系列ごとの整数乗数の個数 K（Fortran 固定幅 i3 で並ぶ）。
/// 主問題(1-3)=4（D,l',l,F）、地球形状(4-9)=5、惑星表1/2(10-21)=11、その他摂動(22-36)=5。
/// 出典: ELP2000-82B 文書 / PROVENANCE.md「行フォーマット」。
fn multiplier_count(file_number: u8) -> usize {
    match file_number {
        1..=3 => 4,
        4..=9 => 5,
        10..=21 => 11,
        _ => 5,
    }
}

/// 1 ファイルのテキストを系列にパースする。`file_number` が主問題(1-3)/摂動(4-36)を決める。
/// 先頭行は見出し（スキップ）。各項行は **Fortran 固定幅 i3**（3 文字幅）で K 個の整数乗数が並び、
/// 続く空白区切りの小数列が係数（主問題=\[A, B1..B6\]、摂動=\[φ, A\]）。i3 は隣接フィールドが空白なしで
/// 接する（例 ELP10 の "4-11" = " 4" の続きに "-11"）ため、空白分割ではなく固定幅で切る。
pub fn parse_series(file_number: u8, text: &str) -> Result<Elp82bSeries, XtaskError> {
    let is_main = file_number <= 3;
    let k = multiplier_count(file_number);
    let mut terms = Vec::new();
    // 先頭行は見出し（"MAIN PROBLEM..." 等）。以降が項行。
    for line in text.lines().skip(1) {
        if line.trim().is_empty() {
            continue;
        }
        // 乗数フィールドは ASCII の固定幅 i3。先頭 K*3 文字が乗数群。
        if line.len() < k * 3 {
            return Err(XtaskError::MalformedSource(format!(
                "ELP{file_number}: line shorter than {} multiplier columns: {line:?}",
                k * 3
            )));
        }
        let mut multipliers = Vec::with_capacity(k);
        for i in 0..k {
            let field = line[i * 3..i * 3 + 3].trim();
            let m = field.parse::<i32>().map_err(|_| {
                XtaskError::MalformedSource(format!(
                    "ELP{file_number}: non-integer multiplier field {field:?}: {line:?}"
                ))
            })?;
            multipliers.push(m);
        }
        // 乗数群の後ろは空白区切りの小数列。主問題=[A, B1..B6]（7 実数）、摂動=[φ, A]。
        let floats: Vec<&str> = line[k * 3..].split_whitespace().collect();
        let coeff = |i: usize| -> Result<f64, XtaskError> {
            floats
                .get(i)
                .and_then(|t| t.parse::<f64>().ok())
                .ok_or_else(|| {
                    XtaskError::MalformedSource(format!(
                        "ELP{file_number}: bad coefficient: {line:?}"
                    ))
                })
        };
        // 主問題は振幅 A＋偏微分 6 列（`4i3,2x,f13.5,6(2x,f10.2)`）の全 7 実数を忠実保存。
        // 評価器(ISSUE-014)は coef(1..6) を DE200/LE200 フィット補正に使う（coef(7) は未使用）。
        // 摂動は [φ(度), A]（period は不使用）。
        let coefficients = if is_main {
            (0..7)
                .map(coeff)
                .collect::<Result<Vec<f64>, XtaskError>>()?
        } else {
            vec![coeff(0)?, coeff(1)?]
        };
        terms.push(Elp82bTerm {
            multipliers,
            coefficients,
        });
    }
    Ok(Elp82bSeries {
        file: file_number,
        terms,
    })
}

/// モデルを packed（flat little-endian f64）へ直列化する（決定的）。
/// レイアウト: `[n_series, <各系列 = [file, n_mult, n_coeff, n_terms, <各項 = n_mult 乗数 + n_coeff 係数>]>...]`。
pub fn pack_model(model: &Elp82bModel) -> Vec<u8> {
    let mut values: Vec<f64> = Vec::new();
    values.push(f64::from(
        i32::try_from(model.series.len()).expect("series count fits in i32"),
    ));
    for series in &model.series {
        let n_mult = series.terms.first().map_or(0, |t| t.multipliers.len());
        let n_coeff = series.terms.first().map_or(0, |t| t.coefficients.len());
        values.push(f64::from(series.file));
        values.push(f64::from(
            i32::try_from(n_mult).expect("n_mult fits in i32"),
        ));
        values.push(f64::from(
            i32::try_from(n_coeff).expect("n_coeff fits in i32"),
        ));
        values.push(f64::from(
            i32::try_from(series.terms.len()).expect("term count fits in i32"),
        ));
        for term in &series.terms {
            for &m in &term.multipliers {
                values.push(f64::from(m));
            }
            for &c in &term.coefficients {
                values.push(c);
            }
        }
    }
    pack_f64_le(&values)
}

/// 検証済み f64 を非負カウントへ。
#[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
fn f64_to_count(value: f64) -> Result<usize, XtaskError> {
    if !value.is_finite() || value < 0.0 || value.fract() != 0.0 || value > MAX_COUNT {
        return Err(XtaskError::MalformedPacked(format!(
            "invalid count {value}"
        )));
    }
    Ok(value as usize)
}

/// 検証済み f64 を u8 識別子（file 番号）へ。
#[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
fn f64_to_u8(value: f64) -> Result<u8, XtaskError> {
    if !value.is_finite() || value < 0.0 || value.fract() != 0.0 || value > 255.0 {
        return Err(XtaskError::MalformedPacked(format!("invalid u8 {value}")));
    }
    Ok(value as u8)
}

/// 検証済み f64 を整数乗数 i32 へ。
#[allow(clippy::cast_possible_truncation)]
fn f64_to_i32(value: f64) -> Result<i32, XtaskError> {
    if !value.is_finite()
        || value.fract() != 0.0
        || value < f64::from(i32::MIN)
        || value > f64::from(i32::MAX)
    {
        return Err(XtaskError::MalformedPacked(format!(
            "invalid multiplier {value}"
        )));
    }
    Ok(value as i32)
}

/// `values[*cursor]` を境界チェック付きで取り出す。
fn take(values: &[f64], cursor: &mut usize) -> Result<f64, XtaskError> {
    let value = values
        .get(*cursor)
        .copied()
        .ok_or_else(|| XtaskError::MalformedPacked("ELP packed ended mid-record".to_string()))?;
    *cursor += 1;
    Ok(value)
}

/// packed バイト列を [`Elp82bModel`] へ復元する。
pub fn unpack_model(bytes: &[u8]) -> Result<Elp82bModel, XtaskError> {
    let values = unpack_f64_le(bytes)?;
    let mut cursor = 0usize;
    let n_series = f64_to_count(take(&values, &mut cursor)?)?;
    let mut series = Vec::with_capacity(n_series);
    for _ in 0..n_series {
        let file = f64_to_u8(take(&values, &mut cursor)?)?;
        let n_mult = f64_to_count(take(&values, &mut cursor)?)?;
        let n_coeff = f64_to_count(take(&values, &mut cursor)?)?;
        let n_terms = f64_to_count(take(&values, &mut cursor)?)?;
        let mut terms = Vec::with_capacity(n_terms);
        for _ in 0..n_terms {
            let mut multipliers = Vec::with_capacity(n_mult);
            for _ in 0..n_mult {
                multipliers.push(f64_to_i32(take(&values, &mut cursor)?)?);
            }
            let mut coefficients = Vec::with_capacity(n_coeff);
            for _ in 0..n_coeff {
                coefficients.push(take(&values, &mut cursor)?);
            }
            terms.push(Elp82bTerm {
                multipliers,
                coefficients,
            });
        }
        series.push(Elp82bSeries { file, terms });
    }
    if cursor != values.len() {
        return Err(XtaskError::MalformedPacked(format!(
            "ELP packed has {} trailing f64 after {n_series} series",
            values.len() - cursor
        )));
    }
    Ok(Elp82bModel { series })
}

/// 36 ファイル `(file_number, text)` から packed バイト列と [`DataSetMetadata`] を構成する（純関数・決定的）。
pub fn build_artifact(files: &[(u8, String)]) -> Result<(Vec<u8>, DataSetMetadata), XtaskError> {
    let series = files
        .iter()
        .map(|(file_number, text)| parse_series(*file_number, text))
        .collect::<Result<Vec<_>, _>>()?;
    let model = Elp82bModel { series };
    let bytes = pack_model(&model);
    let checksum = sha256_hex(&bytes);
    let metadata = DataSetMetadata {
        name: "elp2000-82b-moon".to_string(),
        version: "ELP2000-82B".to_string(),
        source: "IMCCE ELP2000-82B (ELP1..ELP36); Chapront-Touzé & Chapront (1983/1988)"
            .to_string(),
        license: "IMCCE scientific data (attribution); regenerated from primary, GPL not used"
            .to_string(),
        valid_from: "1900-01-01".to_string(),
        valid_to: "2100-01-01".to_string(),
        checksum,
    };
    Ok((bytes, metadata))
}

/// `SOURCE_DIR` 配下の 36 ファイルを `(file_number, text)` として読む。
fn read_sources() -> Result<Vec<(u8, String)>, XtaskError> {
    (1..=N_FILES)
        .map(|n| {
            let file_number = u8::try_from(n).expect("1..=36 fits u8");
            let path = format!("{SOURCE_DIR}/ELP{n}");
            let text =
                fs::read_to_string(&path).map_err(|source| XtaskError::Io { path, source })?;
            Ok((file_number, text))
        })
        .collect()
}

/// `metadata.txt` の人間可読レンダリング。
fn render_metadata(m: &DataSetMetadata) -> String {
    format!(
        "name = {}\nversion = {}\nsource = {}\nlicense = {}\nvalid_from = {}\nvalid_to = {}\nchecksum = {}\n",
        m.name, m.version, m.source, m.license, m.valid_from, m.valid_to, m.checksum,
    )
}

/// 一次原データ（ELP1..ELP36）を読み、packed・metadata・NOTICE を `generated/elp2000-82b/` へ書き出す。
pub fn generate_to_disk() -> Result<DataSetMetadata, XtaskError> {
    let files = read_sources()?;
    let (bytes, metadata) = build_artifact(&files)?;
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
            "# Generated ELP2000-82B Moon coefficients\n\n\
             Generated by `cargo xtask generate-coefficients --dataset elp2000-82b` from\n\
             `{SOURCE_DIR}/ELP1..ELP36` (see PROVENANCE.md). Do not edit by hand; regenerate and\n\
             verify with `cargo xtask verify-generated --dataset elp2000-82b`.\n\n{}\n",
            render_metadata(&metadata)
        )
        .as_bytes(),
    )?;
    Ok(metadata)
}

/// コミット済み packed 係数が一次原データから決定的に再生成できることを検証する。
pub fn verify_against_disk() -> Result<(), XtaskError> {
    let files = read_sources()?;
    let (regenerated, _) = build_artifact(&files)?;
    let committed_path = format!("{GENERATED_DIR}/{PACKED_NAME}");
    let committed = fs::read(&committed_path).map_err(|source| XtaskError::Io {
        path: committed_path,
        source,
    })?;
    crate::compare_checksum("elp2000-82b", &sha256_hex(&committed), &regenerated)
}

#[cfg(test)]
// ELP2000-82B の乗数（整数）・係数（A / φ,A）は原データの 10 進が f64 へ一意変換されるため、
// スポット値はリテラルの厳密 `==` で照合する（`clippy::float_cmp` を mod 全体で許容）。
#[allow(clippy::float_cmp)]
mod tests {
    use super::*;
    use crate::checksum::sha256_hex;

    // ------------------------------------------------------------------
    // 原データ（コミット済み実ファイル ELP1..ELP36）。manifest 相対 include_str! で
    // cwd 非依存に読み込む（src→xtask→crates→repo root）。
    // 出典: data/coefficient-source/elp2000-82b/PROVENANCE.md
    //   （IMCCE, ELP2000-82B, Chapront-Touzé & Chapront 1983/1988）。
    // 各ファイル: 先頭 1 行が見出し（スキップ）、以降が項行。
    //   主問題 ELP1-3: `i1 i2 i3 i4   A   B1..B6` → 乗数=[i1..i4], 係数=[A, B1..B6]（7 個）。
    //   摂動  ELP4-36: `m1..mK   φ   A   period`  → 乗数=[m1..mK], 係数=[φ, A]。
    // ------------------------------------------------------------------
    const ELP_FILES: [&str; N_FILES] = [
        include_str!("../../../data/coefficient-source/elp2000-82b/ELP1"),
        include_str!("../../../data/coefficient-source/elp2000-82b/ELP2"),
        include_str!("../../../data/coefficient-source/elp2000-82b/ELP3"),
        include_str!("../../../data/coefficient-source/elp2000-82b/ELP4"),
        include_str!("../../../data/coefficient-source/elp2000-82b/ELP5"),
        include_str!("../../../data/coefficient-source/elp2000-82b/ELP6"),
        include_str!("../../../data/coefficient-source/elp2000-82b/ELP7"),
        include_str!("../../../data/coefficient-source/elp2000-82b/ELP8"),
        include_str!("../../../data/coefficient-source/elp2000-82b/ELP9"),
        include_str!("../../../data/coefficient-source/elp2000-82b/ELP10"),
        include_str!("../../../data/coefficient-source/elp2000-82b/ELP11"),
        include_str!("../../../data/coefficient-source/elp2000-82b/ELP12"),
        include_str!("../../../data/coefficient-source/elp2000-82b/ELP13"),
        include_str!("../../../data/coefficient-source/elp2000-82b/ELP14"),
        include_str!("../../../data/coefficient-source/elp2000-82b/ELP15"),
        include_str!("../../../data/coefficient-source/elp2000-82b/ELP16"),
        include_str!("../../../data/coefficient-source/elp2000-82b/ELP17"),
        include_str!("../../../data/coefficient-source/elp2000-82b/ELP18"),
        include_str!("../../../data/coefficient-source/elp2000-82b/ELP19"),
        include_str!("../../../data/coefficient-source/elp2000-82b/ELP20"),
        include_str!("../../../data/coefficient-source/elp2000-82b/ELP21"),
        include_str!("../../../data/coefficient-source/elp2000-82b/ELP22"),
        include_str!("../../../data/coefficient-source/elp2000-82b/ELP23"),
        include_str!("../../../data/coefficient-source/elp2000-82b/ELP24"),
        include_str!("../../../data/coefficient-source/elp2000-82b/ELP25"),
        include_str!("../../../data/coefficient-source/elp2000-82b/ELP26"),
        include_str!("../../../data/coefficient-source/elp2000-82b/ELP27"),
        include_str!("../../../data/coefficient-source/elp2000-82b/ELP28"),
        include_str!("../../../data/coefficient-source/elp2000-82b/ELP29"),
        include_str!("../../../data/coefficient-source/elp2000-82b/ELP30"),
        include_str!("../../../data/coefficient-source/elp2000-82b/ELP31"),
        include_str!("../../../data/coefficient-source/elp2000-82b/ELP32"),
        include_str!("../../../data/coefficient-source/elp2000-82b/ELP33"),
        include_str!("../../../data/coefficient-source/elp2000-82b/ELP34"),
        include_str!("../../../data/coefficient-source/elp2000-82b/ELP35"),
        include_str!("../../../data/coefficient-source/elp2000-82b/ELP36"),
    ];

    /// `build_artifact` 用の `(file_number, text)` ×36 を実ファイルから構築する。
    /// file_number は 1..=36（ファイル順）。
    fn all_files() -> Vec<(u8, String)> {
        ELP_FILES
            .iter()
            .enumerate()
            .map(|(i, text)| (u8::try_from(i + 1).unwrap(), (*text).to_string()))
            .collect()
    }

    /// 実ファイルから全 36 系列をパースしてモデルを組む（build_artifact を介さない経路）。
    fn parse_all() -> Elp82bModel {
        let series = ELP_FILES
            .iter()
            .enumerate()
            .map(|(i, text)| {
                let file_number = u8::try_from(i + 1).unwrap();
                parse_series(file_number, text)
                    .unwrap_or_else(|e| panic!("ELP{file_number} parses: {e:?}"))
            })
            .collect();
        Elp82bModel { series }
    }

    /// 各ファイルの期待項数 = 「実ファイル行数 − 1（見出し）」。
    /// PROVENANCE.md の系列マップと一致（Σ=37872）。出典: 各 ELP ファイルの行数。
    const EXPECTED_TERMS: [usize; N_FILES] = [
        1023, 918, 704, // ELP1-3 主問題（計 2645）
        347, 316, 237, 14, 11, 8, // ELP4-9 地球形状摂動
        14328, 5233, 6631, 4384, 833, 1715, // ELP10-15 惑星摂動 表1
        170, 150, 114, 226, 188, 169, // ELP16-21 惑星摂動 表2
        3, 2, 2, 6, 4, 5, // ELP22-27 潮汐
        20, 12, 14, // ELP28-30 月形状摂動
        11, 4, 10, // ELP31-33 相対論摂動
        28, 13, 19, // ELP34-36 惑星摂動(太陽離心率)
    ];

    /// 総項数（主問題 2645 + 摂動 35227）。
    const TOTAL_TERMS: usize = 37872;

    // ==================================================================
    // (1) スポット値（主問題 + 摂動）: 各ファイル第1項の乗数・係数を厳密一致で確認。
    //     「最初の小数点トークンまでが整数乗数」規約・列ズレ・B列/period 誤読を潰す。
    //     主問題=乗数4個+A1個、摂動=乗数K個+[φ,A]2個。
    // ==================================================================

    /// ELP1（主問題・経度 sine）第1項（ELP1:2）:
    ///   `  0  0  0  2     -411.60287      168.48   -18433.81     -121.62        0.40       -0.18        0.00`
    /// 乗数=[0,0,0,2]（4 個・D,l',l,F）、係数=coef(1..7)=[A, B1..B6]（全 7 floats を file 順で保持）。
    /// elp82b_1（DE200/LE200 適合補正）が coef(1..6) を消費するため、A だけでなく続く B 偏微分列も全て取り込む。
    #[test]
    fn elp1_main_first_term_spot_value() {
        let s = parse_series(1, ELP_FILES[0]).expect("ELP1 parses");
        let t = &s.terms[0];
        assert_eq!(t.multipliers, vec![0, 0, 0, 2], "ELP1 first multipliers");
        assert_eq!(
            t.coefficients,
            vec![-411.60287, 168.48, -18433.81, -121.62, 0.40, -0.18, 0.00],
            "ELP1 first coeff = coef(1..7) = [A, B1..B6]"
        );
        assert_eq!(t.multipliers.len(), 4, "main problem: 4 multipliers");
        assert_eq!(
            t.coefficients.len(),
            7,
            "main problem: 7 coefficients (A, B1..B6)"
        );
    }

    /// ELP3（主問題・距離 cosine）第1項（ELP3:2、距離定数 ~385000 km）:
    ///   `  0  0  0  0   385000.52719    -7992.63      -11.06    21578.08       -4.53       11.39       -0.06`
    /// 乗数=[0,0,0,0]、係数=coef(1..7)=[A, B1..B6]。主問題で全乗数 0・大きな振幅・全 7 floats 取込みの境界を独立確認。
    #[test]
    fn elp3_main_first_term_spot_value() {
        let s = parse_series(3, ELP_FILES[2]).expect("ELP3 parses");
        let t = &s.terms[0];
        assert_eq!(t.multipliers, vec![0, 0, 0, 0], "ELP3 first multipliers");
        assert_eq!(
            t.coefficients,
            vec![
                385000.52719,
                -7992.63,
                -11.06,
                21578.08,
                -4.53,
                11.39,
                -0.06
            ],
            "ELP3 first coeff = coef(1..7) = [A, B1..B6]"
        );
        assert_eq!(t.multipliers.len(), 4, "main problem: 4 multipliers");
        assert_eq!(
            t.coefficients.len(),
            7,
            "main problem: 7 coefficients (A, B1..B6)"
        );
    }

    /// ELP4（地球形状摂動・経度）第1項（ELP4:2）:
    ///   `  0  0  0  0  1 270.00000   0.00003     0.075`
    /// 乗数=[0,0,0,0,1]（5 個）、係数=[φ=270.00000, A=0.00003]（period=0.075 不使用）。
    #[test]
    fn elp4_perturbation_first_term_spot_value() {
        let s = parse_series(4, ELP_FILES[3]).expect("ELP4 parses");
        let t = &s.terms[0];
        assert_eq!(t.multipliers, vec![0, 0, 0, 0, 1], "ELP4 first multipliers");
        assert_eq!(
            t.coefficients,
            vec![270.00000, 0.00003],
            "ELP4 first coeff = [phi, A]"
        );
        assert_eq!(t.multipliers.len(), 5, "earth-figure pert: K=5 multipliers");
        assert_eq!(
            t.coefficients.len(),
            2,
            "perturbation: 2 coefficients (phi, A)"
        );
    }

    /// ELP10（惑星摂動 表1・経度）第1項（ELP10:2）:
    ///   `  0  0  0  0  0  0  0  0  0  0  2 359.99831   0.00020     0.037`
    /// 乗数=[0,…,0,2]（11 個）、係数=[φ=359.99831, A=0.00020]。K=11 の系列で乗数個数を検証。
    #[test]
    fn elp10_perturbation_first_term_spot_value() {
        let s = parse_series(10, ELP_FILES[9]).expect("ELP10 parses");
        let t = &s.terms[0];
        assert_eq!(
            t.multipliers,
            vec![0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2],
            "ELP10 first multipliers (11)"
        );
        assert_eq!(
            t.coefficients,
            vec![359.99831, 0.00020],
            "ELP10 first coeff = [phi, A]"
        );
        assert_eq!(
            t.multipliers.len(),
            11,
            "planetary table1: K=11 multipliers"
        );
        assert_eq!(t.coefficients.len(), 2, "perturbation: 2 coefficients");
    }

    /// ELP22（潮汐・経度）第1項（ELP22:2）:
    ///   `  0  1  1 -1 -1 192.93665   0.00004     0.076`
    /// 乗数=[0,1,1,-1,-1]（5 個、負の乗数を含む）、係数=[φ=192.93665, A=0.00004]。
    /// 負の整数乗数の符号取り違えを潰す。
    #[test]
    fn elp22_perturbation_first_term_spot_value() {
        let s = parse_series(22, ELP_FILES[21]).expect("ELP22 parses");
        let t = &s.terms[0];
        assert_eq!(
            t.multipliers,
            vec![0, 1, 1, -1, -1],
            "ELP22 first multipliers (with negatives)"
        );
        assert_eq!(
            t.coefficients,
            vec![192.93665, 0.00004],
            "ELP22 first coeff = [phi, A]"
        );
        assert_eq!(t.multipliers.len(), 5, "tidal: K=5 multipliers");
        assert_eq!(t.coefficients.len(), 2, "perturbation: 2 coefficients");
    }

    /// ELP34（惑星摂動・太陽離心率・経度）第1項（ELP34:2）:
    ///   `  0  0  1 -2  0   0.00000   0.00007     0.039`
    /// 乗数=[0,0,1,-2,0]（5 個）、係数=[φ=0.00000, A=0.00007]。
    /// φ=0 の項で「φ も小数トークンとして A の前に拾う」ことを確認（φ を読み飛ばさない）。
    #[test]
    fn elp34_perturbation_first_term_spot_value() {
        let s = parse_series(34, ELP_FILES[33]).expect("ELP34 parses");
        let t = &s.terms[0];
        assert_eq!(
            t.multipliers,
            vec![0, 0, 1, -2, 0],
            "ELP34 first multipliers"
        );
        assert_eq!(
            t.coefficients,
            vec![0.00000, 0.00007],
            "ELP34 first coeff = [phi=0, A]"
        );
        assert_eq!(t.multipliers.len(), 5, "solar-ecc pert: K=5 multipliers");
        assert_eq!(t.coefficients.len(), 2, "perturbation: 2 coefficients");
    }

    // ==================================================================
    // (2) 項数: parse_series の terms.len() が各ファイルの「行数 − 1」と一致。
    //     主要ファイルのリテラル（ELP1=1023, ELP4=347, ELP10=14328, ELP22=3）も独立確認。
    // ==================================================================

    /// 主要 4 ファイルの terms.len() をリテラルで確認（EXPECTED_TERMS とは独立の固定値）。
    #[test]
    fn parse_series_term_counts_for_key_files() {
        assert_eq!(parse_series(1, ELP_FILES[0]).unwrap().terms.len(), 1023);
        assert_eq!(parse_series(4, ELP_FILES[3]).unwrap().terms.len(), 347);
        assert_eq!(parse_series(10, ELP_FILES[9]).unwrap().terms.len(), 14328);
        assert_eq!(parse_series(22, ELP_FILES[21]).unwrap().terms.len(), 3);
    }

    /// 全 36 ファイルの terms.len() が EXPECTED_TERMS（行数 − 1）と一致し、
    /// file 番号が i+1（ファイル順）であること。
    #[test]
    fn parse_series_all_files_term_counts_and_file_numbers() {
        for (i, text) in ELP_FILES.iter().enumerate() {
            let file_number = u8::try_from(i + 1).unwrap();
            let s = parse_series(file_number, text).expect("file parses");
            assert_eq!(s.file, file_number, "ELP{file_number}: series.file");
            assert_eq!(
                s.terms.len(),
                EXPECTED_TERMS[i],
                "ELP{file_number}: term count (lines - 1 header)"
            );
        }
    }

    /// build_artifact を介してモデル化したとき: series.len()==36、series[i].file==i+1、
    /// 総項数 Σ == 37872（主問題 2645 + 摂動 35227）。
    #[test]
    fn full_model_has_36_series_in_order_and_total_37872_terms() {
        let model = parse_all();
        assert_eq!(model.series.len(), N_FILES, "36 series");
        assert_eq!(N_FILES, 36, "N_FILES == 36");
        for (i, s) in model.series.iter().enumerate() {
            assert_eq!(
                s.file,
                u8::try_from(i + 1).unwrap(),
                "series[{i}].file == i+1 (file order)"
            );
        }
        let total: usize = model.series.iter().map(|s| s.terms.len()).sum();
        assert_eq!(total, TOTAL_TERMS, "total terms = 37872");
        // EXPECTED_TERMS の和も 37872（表自体の整合）。
        assert_eq!(EXPECTED_TERMS.iter().sum::<usize>(), TOTAL_TERMS);

        let main_total: usize = model.series[0..3].iter().map(|s| s.terms.len()).sum();
        assert_eq!(main_total, 2645, "main problem (ELP1-3) = 2645 terms");
        let pert_total: usize = model.series[3..].iter().map(|s| s.terms.len()).sum();
        assert_eq!(pert_total, 35227, "perturbations (ELP4-36) = 35227 terms");
    }

    // ==================================================================
    // (3) packed 往復: 全 36 ファイル → Elp82bModel → pack_model → unpack_model が
    //     元モデルと完全一致（assert_eq!）。
    // ==================================================================
    #[test]
    fn pack_unpack_round_trips_full_model() {
        let model = parse_all();
        let restored = unpack_model(&pack_model(&model)).expect("packed model round-trips");
        assert_eq!(restored, model, "round-tripped model must equal original");
    }

    /// build_artifact のバイト列も unpack_model でモデルへ戻せる（同じ packed 契約）。
    #[test]
    fn build_artifact_bytes_unpack_to_full_model() {
        let model = parse_all();
        let (bytes, _) = build_artifact(&all_files()).expect("artifact builds");
        let restored = unpack_model(&bytes).expect("artifact bytes unpack");
        assert_eq!(restored, model, "artifact bytes == packed full model");
    }

    // ==================================================================
    // (4a) pack 決定性・長さ。
    //   レイアウト: [n_series, <系列 = [file, n_mult, n_coeff, n_terms,
    //               <項 = n_mult 乗数 + n_coeff 係数>]>...]。
    //   f64 数 = 1 + Σ_series(4 + n_terms*(n_mult + n_coeff))。
    //   主問題 n_mult=4,n_coeff=7; 摂動 n_mult=K,n_coeff=2。
    // ==================================================================
    #[test]
    fn pack_model_is_deterministic_and_has_exact_length() {
        let model = parse_all();
        let a = pack_model(&model);
        let b = pack_model(&model);
        assert_eq!(a, b, "pack_model must be deterministic");

        // モデルからレイアウトどおりに期待 f64 数を計算（実装の n_mult/n_coeff に依存させる）。
        let mut expected_f64 = 1usize; // n_series
        for s in &model.series {
            let n_mult = s.terms.first().map_or(0, |t| t.multipliers.len());
            let n_coeff = s.terms.first().map_or(0, |t| t.coefficients.len());
            // 系列内の全項で n_mult/n_coeff が一定であることも検証（パース規約の不変条件）。
            for t in &s.terms {
                assert_eq!(t.multipliers.len(), n_mult, "uniform n_mult in series");
                assert_eq!(t.coefficients.len(), n_coeff, "uniform n_coeff in series");
            }
            // 主問題(1-3)=4+7、摂動(4-36)=K+2 の対応も確認。
            if s.file <= 3 {
                assert_eq!(n_mult, 4, "main problem n_mult=4");
                assert_eq!(n_coeff, 7, "main problem n_coeff=7");
            } else {
                assert_eq!(n_coeff, 2, "perturbation n_coeff=2");
            }
            expected_f64 += 4 + s.terms.len() * (n_mult + n_coeff);
        }
        assert_eq!(
            a.len(),
            8 * expected_f64,
            "packed length = 8 * (1 + Σ(4 + n_terms*(n_mult+n_coeff)))"
        );
    }

    // ==================================================================
    // (4b) 合成テキスト: 主問題形式 1 項 / 摂動形式 1 項を parse_series が正しく解釈。
    //      乗数が空（先頭トークンが小数）の異常行は MalformedSource（パース規約の境界）。
    // ==================================================================

    /// 合成主問題（1 項）: 見出し行 + `i1 i2 i3 i4   A   B1..B6`。
    /// 4 整数 i3 までが乗数、続く 7 floats coef(1..7)=[A, B1..B6] が全て係数。
    #[test]
    fn parse_series_main_synthetic_single_term() {
        let text = "MAIN PROBLEM. LONGITUDE(SINE)\n  1 -2  3 -4     12.50000      1.0   2.0   3.0   4.0   5.0   6.0\n";
        let s = parse_series(1, text).expect("synthetic main parses");
        assert_eq!(s.file, 1);
        assert_eq!(s.terms.len(), 1, "one term row");
        assert_eq!(s.terms[0].multipliers, vec![1, -2, 3, -4]);
        assert_eq!(
            s.terms[0].coefficients,
            vec![12.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
            "coef(1..7) = [A, B1..B6]"
        );
    }

    /// 合成摂動（1 項）: 見出し行 + `m1..m5   φ   A   period`。
    /// φ と A の 2 つが係数、period は無視される。
    #[test]
    fn parse_series_perturbation_synthetic_single_term() {
        let text = "EARTH FIGURE PERTURBATIONS. LONGITUDE\n  0  1 -1  2 -3 123.45000   0.06789     0.500\n";
        let s = parse_series(4, text).expect("synthetic perturbation parses");
        assert_eq!(s.file, 4);
        assert_eq!(s.terms.len(), 1, "one term row");
        assert_eq!(s.terms[0].multipliers, vec![0, 1, -1, 2, -3]);
        assert_eq!(
            s.terms[0].coefficients,
            vec![123.45, 0.06789],
            "[phi, A]; period ignored"
        );
    }

    /// 乗数が空（最初のトークンが小数）の異常項行は MalformedSource。
    /// 「最初の小数トークンの手前までが整数乗数」規約で乗数 0 個になる退化を弾く。
    #[test]
    fn parse_series_rejects_term_with_no_integer_multipliers() {
        let text = "MAIN PROBLEM. LONGITUDE(SINE)\n   12.50000   1.0   2.0   3.0\n";
        let err = parse_series(1, text).expect_err("leading decimal (no multipliers) must error");
        assert!(
            matches!(err, XtaskError::MalformedSource(_)),
            "expected MalformedSource for empty multipliers, got {err:?}"
        );
    }

    // ==================================================================
    // (5) build_artifact（純関数）: packed バイト・DataSetMetadata の整合。
    // ==================================================================
    #[test]
    fn build_artifact_checksum_and_metadata_are_consistent() {
        let (bytes, metadata) = build_artifact(&all_files()).expect("artifact builds");
        // packed バイトは pack_model（全モデル）と同一（純関数の決定性）。
        let model = parse_all();
        assert_eq!(bytes, pack_model(&model), "artifact bytes == pack_model");
        // checksum は packed バイトの SHA-256。
        assert_eq!(
            metadata.checksum,
            sha256_hex(&bytes),
            "checksum = sha256(packed bytes)"
        );
        // provenance 完全（全フィールド非空）。
        assert!(
            metadata.has_complete_provenance(),
            "metadata must have complete provenance: {metadata:?}"
        );
        // 版識別子は ELP2000-82B 系。
        assert!(
            metadata.version.contains("ELP2000-82B"),
            "version should denote ELP2000-82B: {:?}",
            metadata.version
        );
        // 名称は月/elp を示す。
        let name_lower = metadata.name.to_ascii_lowercase();
        assert!(
            name_lower.contains("elp") || name_lower.contains("moon"),
            "name should reference elp/moon: {:?}",
            metadata.name
        );
    }

    /// build_artifact は決定的（同一入力 → 同一バイト・同一 checksum）。
    #[test]
    fn build_artifact_is_deterministic() {
        let (bytes_a, meta_a) = build_artifact(&all_files()).expect("artifact builds");
        let (bytes_b, meta_b) = build_artifact(&all_files()).expect("artifact builds");
        assert_eq!(bytes_a, bytes_b, "bytes deterministic");
        assert_eq!(meta_a.checksum, meta_b.checksum, "checksum deterministic");
    }

    // ==================================================================
    // (6) unpack 異常: 壊れた packed バイト列で Err（パニックしない）。
    //     長さ検査 / 件数整合のどちらに委ねるかは未確定なので緩く照合。
    // ==================================================================

    /// 8 の倍数でない長さ → 必ず Err（Malformed*）。
    #[test]
    fn unpack_rejects_non_f64_aligned_bytes() {
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

    /// ヘッダで巨大な系列数を宣言するがデータが伴わない（8 の倍数長は満たす）。
    /// [n_series = 1e9] の 1 f64 のみ → 件数とデータ長が矛盾。
    #[test]
    fn unpack_rejects_header_data_length_mismatch() {
        let bytes = crate::packed::pack_f64_le(&[1.0e9_f64]);
        let err = unpack_model(&bytes)
            .expect_err("header series count vs data length mismatch must error");
        assert!(
            matches!(
                err,
                XtaskError::MalformedPacked(_) | XtaskError::MalformedSource(_)
            ),
            "expected Malformed*, got {err:?}"
        );
    }

    // ==================================================================
    // (7) 検証ヘルパ・パース規約の境界（ミューテーション堅牢化）。
    // ==================================================================

    #[test]
    fn f64_to_count_accepts_valid_and_boundary() {
        assert_eq!(f64_to_count(0.0).unwrap(), 0);
        assert_eq!(f64_to_count(14328.0).unwrap(), 14328); // 実データ最大項数/系列
        assert_eq!(f64_to_count(MAX_COUNT).unwrap(), 1_000_000); // 上限ちょうど受理
    }

    #[test]
    fn f64_to_count_rejects_invalid() {
        assert!(f64_to_count(MAX_COUNT + 1.0).is_err(), "above max");
        assert!(f64_to_count(-1.0).is_err(), "negative");
        assert!(f64_to_count(1.5).is_err(), "non-integer");
        assert!(f64_to_count(f64::NAN).is_err(), "NaN");
        assert!(f64_to_count(f64::INFINITY).is_err(), "infinity");
    }

    #[test]
    fn f64_to_u8_accepts_valid_and_boundary() {
        assert_eq!(f64_to_u8(0.0).unwrap(), 0);
        assert_eq!(f64_to_u8(36.0).unwrap(), 36); // 最大 file 番号
        assert_eq!(f64_to_u8(255.0).unwrap(), 255); // 上限ちょうど受理
    }

    #[test]
    fn f64_to_u8_rejects_invalid() {
        assert!(f64_to_u8(256.0).is_err(), "above 255");
        assert!(f64_to_u8(-1.0).is_err(), "negative");
        assert!(f64_to_u8(1.5).is_err(), "non-integer");
        assert!(f64_to_u8(f64::NAN).is_err(), "NaN");
        assert!(f64_to_u8(f64::INFINITY).is_err(), "infinity");
    }

    /// f64_to_i32: 0・典型（負を含む）・i32 両端ちょうどは受理。範囲外/非整数/非有限は Err。
    /// 乗数は負値を取りうる（ELP22 等）ため負の典型値も確認する。
    #[test]
    fn f64_to_i32_accepts_valid_and_boundary() {
        assert_eq!(f64_to_i32(0.0).unwrap(), 0);
        assert_eq!(f64_to_i32(-69.0).unwrap(), -69); // 実データ範囲の負乗数
        assert_eq!(f64_to_i32(51.0).unwrap(), 51);
        assert_eq!(f64_to_i32(f64::from(i32::MAX)).unwrap(), i32::MAX); // 上端ちょうど
        assert_eq!(f64_to_i32(f64::from(i32::MIN)).unwrap(), i32::MIN); // 下端ちょうど
    }

    #[test]
    fn f64_to_i32_rejects_invalid() {
        assert!(f64_to_i32(f64::from(i32::MAX) + 1.0).is_err(), "above MAX");
        assert!(f64_to_i32(f64::from(i32::MIN) - 1.0).is_err(), "below MIN");
        assert!(f64_to_i32(1.5).is_err(), "non-integer");
        assert!(f64_to_i32(f64::NAN).is_err(), "NaN");
        assert!(f64_to_i32(f64::INFINITY).is_err(), "infinity");
    }

    /// 乗数列が固定幅 K*3 文字に満たない項行は MalformedSource（パニックしない）。
    /// 摂動 K=5 → 15 文字必要。12 文字行（i3 フィールド 4 個のみ）で境界を検証する。
    /// `line.len() < k*3` の `<`/`*` 改変は、ここでスライス越境 panic または早期 Err で捕捉される。
    #[test]
    fn parse_series_rejects_line_shorter_than_multiplier_columns() {
        // ヘッダ + 12 文字の項行（"  0  0  0  1"）。K=5 では乗数 15 文字に 3 文字不足。
        let text = "EARTH FIGURE PERTURBATIONS. LONGITUDE\n  0  0  0  1\n";
        let err = parse_series(4, text).expect_err("short multiplier columns must error");
        assert!(
            matches!(err, XtaskError::MalformedSource(_)),
            "expected MalformedSource for short line, got {err:?}"
        );
    }

    /// render_metadata は全フィールドを行として出力する（本体→"" / 固定文字列 変異を殺す）。
    #[test]
    fn render_metadata_includes_all_fields() {
        let (_, metadata) = build_artifact(&all_files()).unwrap();
        let rendered = render_metadata(&metadata);
        assert!(rendered.contains("name = elp2000-82b-moon"), "{rendered}");
        assert!(rendered.contains("version = ELP2000-82B"), "{rendered}");
        assert!(
            rendered.contains(&format!("checksum = {}", metadata.checksum)),
            "{rendered}"
        );
    }
}
