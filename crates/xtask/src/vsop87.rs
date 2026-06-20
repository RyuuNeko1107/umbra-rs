//! VSOP87D 地球係数の取り込み（ISSUE-033）。
//!
//! IMCCE 配布の `VSOP87D.ear`（版D・地球・黄道 of date 球面 heliocentric）をパースし、
//! packed 形式（消費側 ISSUE-013 と byte-for-byte 契約）へ直列化する。フォーマットの正本は
//! `data/coefficient-source/vsop87/PROVENANCE.md`。
//!
//! 各項は `amplitude·cos(phase + frequency·T)`（T = ユリウス千年 from J2000 TDB）。
//! セクション = 変数 v（1=L 黄経, 2=B 黄緯, 3=R 動径）× べき α（T^0..T^5）。評価は ISSUE-013。

use crate::checksum::sha256_hex;
use crate::error::XtaskError;
use crate::packed::{pack_f64_le, unpack_f64_le};
use std::fs;
use umbra_core::DataSetMetadata;

/// VSOP87 のべき次数の上限（T^0..T^5）。
pub const N_POWERS: usize = 6;

/// 原データ配置（`cargo xtask` をリポジトリルートから実行する前提。`docs/accuracy.md` §5）。
const SOURCE_DIR: &str = "data/coefficient-source/vsop87";
/// 原データファイル名（版D・地球＝黄道 of date 球面 L,B,R）。
const SOURCE_FILE_D: &str = "VSOP87D.ear";
/// packed 係数ファイル名（版D）。
const PACKED_NAME_D: &str = "vsop87d_earth.bin";
/// 原データファイル名（版A・地球＝黄道 J2000 直交 X,Y,Z。M10 太陽フレーム修正・ISSUE-013/035）。
const SOURCE_FILE_A: &str = "VSOP87A.ear";
/// packed 係数ファイル名（版A）。
const PACKED_NAME_A: &str = "vsop87a_earth.bin";
/// 生成物の出力先。
const GENERATED_DIR: &str = "generated/vsop87";

/// VSOP87 級数 1 項: `amplitude·cos(phase + frequency·T)`。
#[derive(Clone, Debug, PartialEq)]
pub struct Vsop87Term {
    /// 振幅 A（L,B は無次元 rad / R は AU）。
    pub amplitude: f64,
    /// 位相 B \[rad\]。
    pub phase: f64,
    /// 振動数 C \[rad/千年\]。
    pub frequency: f64,
}

/// 1 セクション（変数 v × べき α）の項列。
#[derive(Clone, Debug, PartialEq)]
pub struct Vsop87Series {
    /// 変数（1=L, 2=B, 3=R）。
    pub variable: u8,
    /// T のべき（0..=5）。
    pub power: u8,
    /// 項列。
    pub terms: Vec<Vsop87Term>,
}

/// 地球の VSOP87D モデル全体（ファイル出現順のセクション列）。
#[derive(Clone, Debug, PartialEq)]
pub struct Vsop87Earth {
    /// セクション列（L T0..T5, B T0..T4, R T0..T5）。
    pub series: Vec<Vsop87Series>,
}

/// 健全性のための上限（セクション項数。実データは最大 559）。
const MAX_TERMS: f64 = 100_000.0;

/// セクション見出し行か（`VARIABLE` と `TERMS` を含む）。
fn is_header(line: &str) -> bool {
    line.contains("VARIABLE") && line.contains("TERMS")
}

/// セクション見出しを `(variable, power, declared_terms)` にパースし、版（`version_letter`）・
/// 地球を検査する（B4(c)）。版D=黄道 of date 球面 L,B,R、版A=黄道 J2000 直交 X,Y,Z。
fn parse_header(line: &str, version_letter: char) -> Result<(u8, u8, usize), XtaskError> {
    let tokens: Vec<&str> = line.split_whitespace().collect();
    let malformed = |m: &str| XtaskError::MalformedSource(format!("VSOP87 header: {m}: {line:?}"));

    // 版・地球の検査（EMB/他版/他天体の取り違えは 6.4″/月の系統誤差・PROVENANCE.md §B4(c)）。
    let version_idx = tokens
        .iter()
        .position(|t| *t == "VERSION")
        .ok_or_else(|| malformed("missing VERSION"))?;
    let version = tokens
        .get(version_idx + 1)
        .ok_or_else(|| malformed("missing version id"))?;
    if !version.starts_with(version_letter) {
        return Err(malformed(&format!("not VSOP87 version {version_letter}")));
    }
    let body = tokens
        .get(version_idx + 2)
        .ok_or_else(|| malformed("missing body"))?;
    if *body != "EARTH" {
        return Err(malformed("body is not EARTH"));
    }

    let var_idx = tokens
        .iter()
        .position(|t| *t == "VARIABLE")
        .ok_or_else(|| malformed("missing VARIABLE"))?;
    let variable: u8 = tokens
        .get(var_idx + 1)
        .and_then(|t| t.parse().ok())
        .ok_or_else(|| malformed("bad variable index"))?;

    let power: u8 = tokens
        .iter()
        .find_map(|t| t.strip_prefix("*T**").and_then(|p| p.parse().ok()))
        .ok_or_else(|| malformed("missing *T** power"))?;

    let terms_idx = tokens
        .iter()
        .position(|t| *t == "TERMS")
        .ok_or_else(|| malformed("missing TERMS"))?;
    let declared: usize = terms_idx
        .checked_sub(1)
        .and_then(|i| tokens.get(i))
        .and_then(|t| t.parse().ok())
        .ok_or_else(|| malformed("bad declared term count"))?;

    Ok((variable, power, declared))
}

/// 項行 `<id> <rank> <12 整数乗数>  S  K  A  B  C` を 1 項にパースする。末尾 3 トークン = A,B,C。
/// 見出し・空行は（先頭が整数 id でない・トークン不足で）`None`。
fn parse_term_line(line: &str) -> Option<Vsop87Term> {
    let tokens: Vec<&str> = line.split_whitespace().collect();
    if tokens.len() < 3 {
        return None;
    }
    // 先頭は項識別子（整数）。見出し行（"VSOP87"...）はここで弾かれる。
    tokens[0].parse::<i64>().ok()?;
    let n = tokens.len();
    let amplitude = tokens[n - 3].parse::<f64>().ok()?;
    let phase = tokens[n - 2].parse::<f64>().ok()?;
    let frequency = tokens[n - 1].parse::<f64>().ok()?;
    Some(Vsop87Term {
        amplitude,
        phase,
        frequency,
    })
}

/// `VSOP87<version>.ear` テキストをパースする。版・地球であることを検査し（B4(c)）、各セクション
/// 見出しの宣言項数と実際の項数が食い違えば [`XtaskError::MalformedSource`]。
pub fn parse_earth(text: &str, version_letter: char) -> Result<Vsop87Earth, XtaskError> {
    let mut series: Vec<Vsop87Series> = Vec::new();
    let mut current: Option<(u8, u8, usize, Vec<Vsop87Term>)> = None;
    let finish = |sec: (u8, u8, usize, Vec<Vsop87Term>)| -> Result<Vsop87Series, XtaskError> {
        let (variable, power, declared, terms) = sec;
        if terms.len() != declared {
            return Err(XtaskError::MalformedSource(format!(
                "VSOP87 section (var={variable}, T**{power}) declares {declared} terms but parsed {}",
                terms.len()
            )));
        }
        Ok(Vsop87Series {
            variable,
            power,
            terms,
        })
    };
    for line in text.lines() {
        if is_header(line) {
            if let Some(sec) = current.take() {
                series.push(finish(sec)?);
            }
            let (variable, power, declared) = parse_header(line, version_letter)?;
            current = Some((variable, power, declared, Vec::new()));
        } else if let Some(term) = parse_term_line(line) {
            match current.as_mut() {
                Some((_, _, _, terms)) => terms.push(term),
                None => {
                    return Err(XtaskError::MalformedSource(
                        "VSOP87 term row before any section header".to_string(),
                    ))
                }
            }
        }
    }
    if let Some(sec) = current.take() {
        series.push(finish(sec)?);
    }
    if series.is_empty() {
        return Err(XtaskError::MalformedSource(
            "VSOP87: no sections parsed".to_string(),
        ));
    }
    Ok(Vsop87Earth { series })
}

/// モデルを packed（flat little-endian f64）へ直列化する（決定的）。
/// レイアウト: `[n_sections, <各セクション = [variable, power, n_terms, <各項 amp,phase,freq>]>...]`。
pub fn pack_earth(model: &Vsop87Earth) -> Vec<u8> {
    let total_terms: usize = model.series.iter().map(|s| s.terms.len()).sum();
    let mut values: Vec<f64> = Vec::with_capacity(1 + model.series.len() * 3 + total_terms * 3);
    values.push(f64::from(
        i32::try_from(model.series.len()).expect("section count fits in i32"),
    ));
    for section in &model.series {
        values.push(f64::from(section.variable));
        values.push(f64::from(section.power));
        values.push(f64::from(
            i32::try_from(section.terms.len()).expect("term count fits in i32"),
        ));
        for term in &section.terms {
            values.push(term.amplitude);
            values.push(term.phase);
            values.push(term.frequency);
        }
    }
    pack_f64_le(&values)
}

/// 検証済み f64 を非負カウントへ変換する（負・非整数・非有限・過大は破損）。
#[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
fn f64_to_count(value: f64) -> Result<usize, XtaskError> {
    if !value.is_finite() || value < 0.0 || value.fract() != 0.0 || value > MAX_TERMS {
        return Err(XtaskError::MalformedPacked(format!(
            "invalid count {value}"
        )));
    }
    Ok(value as usize)
}

/// 検証済み f64 を変数/べき識別子（u8）へ変換する。
#[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
fn f64_to_u8(value: f64) -> Result<u8, XtaskError> {
    if !value.is_finite() || value < 0.0 || value.fract() != 0.0 || value > 255.0 {
        return Err(XtaskError::MalformedPacked(format!(
            "invalid u8 field {value}"
        )));
    }
    Ok(value as u8)
}

/// `values[*cursor]` を境界チェック付きで取り出す。
fn take(values: &[f64], cursor: &mut usize) -> Result<f64, XtaskError> {
    let value = values
        .get(*cursor)
        .copied()
        .ok_or_else(|| XtaskError::MalformedPacked("VSOP87 packed ended mid-record".to_string()))?;
    *cursor += 1;
    Ok(value)
}

/// packed バイト列を [`Vsop87Earth`] へ復元する。
pub fn unpack_earth(bytes: &[u8]) -> Result<Vsop87Earth, XtaskError> {
    let values = unpack_f64_le(bytes)?;
    let mut cursor = 0usize;
    let n_sections = f64_to_count(take(&values, &mut cursor)?)?;
    let mut series = Vec::with_capacity(n_sections);
    for _ in 0..n_sections {
        let variable = f64_to_u8(take(&values, &mut cursor)?)?;
        let power = f64_to_u8(take(&values, &mut cursor)?)?;
        let n_terms = f64_to_count(take(&values, &mut cursor)?)?;
        let mut terms = Vec::with_capacity(n_terms);
        for _ in 0..n_terms {
            let amplitude = take(&values, &mut cursor)?;
            let phase = take(&values, &mut cursor)?;
            let frequency = take(&values, &mut cursor)?;
            terms.push(Vsop87Term {
                amplitude,
                phase,
                frequency,
            });
        }
        series.push(Vsop87Series {
            variable,
            power,
            terms,
        });
    }
    if cursor != values.len() {
        return Err(XtaskError::MalformedPacked(format!(
            "VSOP87 packed has {} trailing f64 after {n_sections} sections",
            values.len() - cursor
        )));
    }
    Ok(Vsop87Earth { series })
}

/// 原データテキストから packed バイト列と [`DataSetMetadata`] を構成する（純関数・決定的）。
/// `version_letter`='D'（黄道 of date 球面 L,B,R）or 'A'（黄道 J2000 直交 X,Y,Z）。
pub fn build_artifact(
    text: &str,
    version_letter: char,
) -> Result<(Vec<u8>, DataSetMetadata), XtaskError> {
    let model = parse_earth(text, version_letter)?;
    let bytes = pack_earth(&model);
    let checksum = sha256_hex(&bytes);
    let frame = if version_letter == 'A' {
        "version A, ecliptic J2000 rectangular X,Y,Z"
    } else {
        "version D, ecliptic of date spherical L,B,R"
    };
    let metadata = DataSetMetadata {
        name: "vsop87-earth".to_string(),
        version: format!("VSOP87{version_letter}"),
        source: format!(
            "IMCCE VSOP87{version_letter}.ear (Earth, {frame}); Bretagnon & Francou (1988), A&A 202, 309"
        ),
        license: "IMCCE scientific data (attribution); regenerated from primary, GPL not used"
            .to_string(),
        valid_from: "1900-01-01".to_string(),
        valid_to: "2100-01-01".to_string(),
        checksum,
    };
    Ok((bytes, metadata))
}

/// `SOURCE_DIR` 配下のファイルを読む。
fn read_source(file_name: &str) -> Result<String, XtaskError> {
    let path = format!("{SOURCE_DIR}/{file_name}");
    fs::read_to_string(&path).map_err(|source| XtaskError::Io { path, source })
}

/// `metadata.txt` の人間可読レンダリング。
fn render_metadata(m: &DataSetMetadata) -> String {
    format!(
        "name = {}\nversion = {}\nsource = {}\nlicense = {}\nvalid_from = {}\nvalid_to = {}\nchecksum = {}\n",
        m.name, m.version, m.source, m.license, m.valid_from, m.valid_to, m.checksum,
    )
}

/// 1 版（D or A）の packed 係数・metadata・NOTICE を `generated/vsop87/` へ書き出す。
fn generate_one(
    source_file: &str,
    packed_name: &str,
    version_letter: char,
    meta_name: &str,
    notice_name: &str,
) -> Result<DataSetMetadata, XtaskError> {
    let text = read_source(source_file)?;
    let (bytes, metadata) = build_artifact(&text, version_letter)?;
    fs::create_dir_all(GENERATED_DIR).map_err(|source| XtaskError::Io {
        path: GENERATED_DIR.to_string(),
        source,
    })?;
    let write = |name: &str, content: &[u8]| -> Result<(), XtaskError> {
        let path = format!("{GENERATED_DIR}/{name}");
        fs::write(&path, content).map_err(|source| XtaskError::Io { path, source })
    };
    write(packed_name, &bytes)?;
    write(meta_name, render_metadata(&metadata).as_bytes())?;
    write(
        notice_name,
        format!(
            "# Generated VSOP87{version_letter} Earth coefficients\n\n\
             Generated by `cargo xtask generate-coefficients --dataset vsop87` from\n\
             `{SOURCE_DIR}/{source_file}` (see PROVENANCE.md). Do not edit by hand; regenerate and\n\
             verify with `cargo xtask verify-generated --dataset vsop87`.\n\n{}\n",
            render_metadata(&metadata)
        )
        .as_bytes(),
    )?;
    Ok(metadata)
}

/// 一次原データを読み、版D（黄道 of date L,B,R）と版A（黄道 J2000 X,Y,Z）の packed 係数・
/// metadata・NOTICE を `generated/vsop87/` へ書き出す。版A は M10 太陽フレーム修正（ISSUE-013/035）。
pub fn generate_to_disk() -> Result<DataSetMetadata, XtaskError> {
    let d = generate_one(
        SOURCE_FILE_D,
        PACKED_NAME_D,
        'D',
        "metadata.txt",
        "NOTICE.md",
    )?;
    let a = generate_one(
        SOURCE_FILE_A,
        PACKED_NAME_A,
        'A',
        "metadata_a.txt",
        "NOTICE_a.md",
    )?;
    println!("  vsop87a (checksum {})", a.checksum);
    Ok(d)
}

/// 1 版のコミット済み packed が一次原データから決定的に再生成できることを検証する。
fn verify_one(
    source_file: &str,
    packed_name: &str,
    version_letter: char,
) -> Result<(), XtaskError> {
    let text = read_source(source_file)?;
    let (regenerated, _) = build_artifact(&text, version_letter)?;
    let committed_path = format!("{GENERATED_DIR}/{packed_name}");
    let committed = fs::read(&committed_path).map_err(|source| XtaskError::Io {
        path: committed_path,
        source,
    })?;
    crate::compare_checksum("vsop87", &sha256_hex(&committed), &regenerated)
}

/// コミット済み packed 係数（版D・版A）が一次原データから決定的に再生成できることを検証する。
pub fn verify_against_disk() -> Result<(), XtaskError> {
    verify_one(SOURCE_FILE_D, PACKED_NAME_D, 'D')?;
    verify_one(SOURCE_FILE_A, PACKED_NAME_A, 'A')
}

#[cfg(test)]
// VSOP87D の振幅・位相・振動数は原データの 10 進が f64 へ一意変換されるため、スポット値は
// リテラルの厳密 `==` で照合する（`clippy::float_cmp` を mod 全体で許容）。
#[allow(clippy::float_cmp)]
mod tests {
    use super::*;
    use crate::checksum::sha256_hex;

    // ------------------------------------------------------------------
    // 原データ（コミット済み実ファイル）。manifest 相対 include_str! で
    // cwd 非依存に読み込む（src→xtask→crates→repo root）。
    // 出典: data/coefficient-source/vsop87/PROVENANCE.md（IMCCE, VSOP87
    // 版D・地球, Bretagnon & Francou 1988）。
    // ------------------------------------------------------------------
    const VSOP87D_EAR: &str = include_str!("../../../data/coefficient-source/vsop87/VSOP87D.ear");
    /// 版A（黄道 J2000 直交 X,Y,Z・M10 太陽フレーム修正）。出典は同 PROVENANCE.md。
    const VSOP87A_EAR: &str = include_str!("../../../data/coefficient-source/vsop87/VSOP87A.ear");

    /// 地球の宣言項数表（PROVENANCE.md §系列の構造）。
    /// 各要素 = (variable, power, declared_terms)。出現順 = L T0..T5, B T0..T4, R T0..T5。
    const DECLARED_SECTIONS: &[(u8, u8, usize)] = &[
        // L（variable=1）: T0..T5
        (1, 0, 559),
        (1, 1, 341),
        (1, 2, 142),
        (1, 3, 22),
        (1, 4, 11),
        (1, 5, 5),
        // B（variable=2）: T0..T4（T5 なし）
        (2, 0, 184),
        (2, 1, 99),
        (2, 2, 49),
        (2, 3, 11),
        (2, 4, 5),
        // R（variable=3）: T0..T5
        (3, 0, 526),
        (3, 1, 292),
        (3, 2, 139),
        (3, 3, 27),
        (3, 4, 10),
        (3, 5, 3),
    ];

    /// 全項数の合計（pack 長計算に使用）。Σ declared = 2425。
    const TOTAL_TERMS: usize = 2425;

    /// スポット値検証用ヘルパ: 1 項の amplitude / phase / frequency を厳密一致で確認。
    /// 末尾 3 トークン = A(amplitude), B(phase), C(frequency)。S,K の誤読・列ズレを潰す。
    #[track_caller]
    fn assert_term(term: &Vsop87Term, amplitude: f64, phase: f64, frequency: f64) {
        assert_eq!(term.amplitude, amplitude, "amplitude (A) mismatch");
        assert_eq!(term.phase, phase, "phase (B) mismatch");
        assert_eq!(term.frequency, frequency, "frequency (C) mismatch");
    }

    // ==================================================================
    // (1) セクション構造: 17 セクション、各 (variable, power, terms.len())
    //     が宣言値と一致。順序（L0..5, B0..4, R0..5）も検証する。
    //     参照: VSOP87D.ear の "VARIABLE v (LBR) *T**α N TERMS" 見出し
    //     17 行（grep 行 1,561,903,...,2439）。
    // ==================================================================
    #[test]
    fn parse_earth_has_17_sections_in_declared_order_with_declared_counts() {
        let earth = parse_earth(VSOP87D_EAR, 'D').expect("real VSOP87D.ear parses");
        assert_eq!(earth.series.len(), 17, "L 6 + B 5 + R 6 = 17 sections");
        assert_eq!(
            earth.series.len(),
            DECLARED_SECTIONS.len(),
            "section count matches declaration table"
        );
        for (i, ((var, pow, n), series)) in DECLARED_SECTIONS
            .iter()
            .copied()
            .zip(earth.series.iter())
            .enumerate()
        {
            assert_eq!(series.variable, var, "section {i}: variable");
            assert_eq!(series.power, pow, "section {i}: power");
            assert_eq!(
                series.terms.len(),
                n,
                "section {i} (var={var}, pow={pow}): declared term count"
            );
        }
    }

    /// 先頭セクションの不変条件を独立に確認（順序付けの起点）。
    #[test]
    fn parse_earth_first_section_is_l_t0() {
        let earth = parse_earth(VSOP87D_EAR, 'D').expect("real VSOP87D.ear parses");
        let first = &earth.series[0];
        assert_eq!(first.variable, 1, "first section variable = L(1)");
        assert_eq!(first.power, 0, "first section power = T^0");
        assert_eq!(first.terms.len(), 559, "L T0 = 559 terms");
    }

    // ==================================================================
    // (2) スポット値: series[0]（L T0）の第1項・第2項、および R T0 先頭項を
    //     原ファイルから厳密一致で確認。末尾 3 トークン = A,B,C の対応、
    //     および S,K（先行 2 トークン）の誤読・列ズレ・A/B/C 取り違えを潰す。
    // ==================================================================

    /// L T0 第1項（VSOP87D.ear:2、定数項）:
    ///   ` 4310    1  0 0 0 0 0 0 0 0 0 0 0 0  0.00000000000  1.75347045673  1.75347045673 0.00000000000  0.00000000000`
    /// 末尾 3 トークン = A=1.75347045673, B=0.0, C=0.0。
    #[test]
    fn l_t0_first_term_spot_value() {
        let earth = parse_earth(VSOP87D_EAR, 'D').expect("real VSOP87D.ear parses");
        assert_term(&earth.series[0].terms[0], 1.75347045673, 0.0, 0.0);
    }

    /// L T0 第2項（VSOP87D.ear:3）:
    ///   末尾 3 トークン = A=0.03341656456, B=4.66925680417, C=6283.07584999140。
    /// 直前の S=-0.00748171065, K=-0.03256824823 を誤って拾わないことを保証する。
    #[test]
    fn l_t0_second_term_spot_value() {
        let earth = parse_earth(VSOP87D_EAR, 'D').expect("real VSOP87D.ear parses");
        assert_term(
            &earth.series[0].terms[1],
            0.03341656456,
            4.66925680417,
            6283.07584999140,
        );
    }

    /// R T0 先頭セクション先頭項（VSOP87D.ear:1441、動径定数 ~1.00014 AU）:
    ///   ` 4330    1  0 ... 0  0.00000000000  1.00013988799  1.00013988799 0.00000000000  0.00000000000`
    /// 末尾 3 トークン = A=1.00013988799, B=0.0, C=0.0。
    /// R セクション（variable=3, power=0）の起点を独立に検証する。
    #[test]
    fn r_t0_first_term_spot_value() {
        let earth = parse_earth(VSOP87D_EAR, 'D').expect("real VSOP87D.ear parses");
        // R T0 = series[11]（L 6 + B 5 = 11 セクションの次）。
        let r_t0 = &earth.series[11];
        assert_eq!(r_t0.variable, 3, "series[11] is R(3)");
        assert_eq!(r_t0.power, 0, "series[11] is T^0");
        assert_term(&r_t0.terms[0], 1.00013988799, 0.0, 0.0);
    }

    // ==================================================================
    // (3) 版D・地球チェック（B4(c)）: 版≠D または body≠EARTH を弾く。
    //     合成見出し（VERSION A / VENUS）＋数項で MalformedSource。
    // ==================================================================

    /// 1 セクション・宣言項数と実項数が一致する最小の合成テキスト（版・天体を差し替え可能）。
    /// 見出し: `VSOP87 VERSION {ver}  {body}  VARIABLE 1 (LBR)  *T**0  {n} TERMS  ...`。
    fn synthetic_section(version: &str, body: &str, declared: usize, n_rows: usize) -> String {
        let mut s = format!(
            " VSOP87 VERSION {version}    {body}     VARIABLE 1 (LBR)       *T**0    {declared} TERMS    HELIOCENTRIC DYNAMICAL ECLIPTIC AND EQUINOX OF THE DATE\n"
        );
        for i in 1..=n_rows {
            // <id> <rank> <12 整数乗数>  S  K  A  B  C。末尾 3 = A,B,C。
            s.push_str(&format!(
                " 4310    {i}  0  0  0  0  0  0  0  0  0  0  0  0  0.00000000000     0.00000000000     0.50000000000 1.00000000000    2.00000000000 \n"
            ));
        }
        s
    }

    /// 版D・地球の合成テキストは（構造が整っていれば）パースできる
    /// ＝ チェックが「版≠D / body≠EARTH」だけを弾き、正常系を巻き込まないことの対照。
    #[test]
    fn parse_earth_accepts_synthetic_version_d_earth() {
        let text = synthetic_section("D4", "EARTH", 1, 1);
        let earth = parse_earth(&text, 'D').expect("version D + EARTH synthetic parses");
        assert_eq!(earth.series.len(), 1);
        assert_eq!(earth.series[0].variable, 1);
        assert_eq!(earth.series[0].power, 0);
        // 末尾 3 トークン A,B,C が拾えていること。
        assert_term(&earth.series[0].terms[0], 0.5, 1.0, 2.0);
    }

    /// 版が D でない（VERSION A）見出しは MalformedSource。
    #[test]
    fn parse_earth_rejects_non_version_d() {
        let text = synthetic_section("A", "EARTH", 1, 1);
        let err = parse_earth(&text, 'D').expect_err("version A must be rejected");
        assert!(
            matches!(err, XtaskError::MalformedSource(_)),
            "expected MalformedSource for non-D version, got {err:?}"
        );
    }

    /// 天体が EARTH でない（VENUS）見出しは MalformedSource。
    /// EMB/他天体取り違えは 6.4″/月オーダーの系統誤差（PROVENANCE.md §B4(c)）。
    #[test]
    fn parse_earth_rejects_non_earth_body() {
        let text = synthetic_section("D4", "VENUS", 1, 1);
        let err = parse_earth(&text, 'D').expect_err("VENUS body must be rejected");
        assert!(
            matches!(err, XtaskError::MalformedSource(_)),
            "expected MalformedSource for non-EARTH body, got {err:?}"
        );
    }

    // ==================================================================
    // (4) 項数不整合: 宣言 `*T**0  3 TERMS` だが 2 項のみ → MalformedSource。
    // ==================================================================
    #[test]
    fn parse_earth_rejects_term_count_mismatch() {
        // 宣言 3、実 2 行。
        let text = synthetic_section("D4", "EARTH", 3, 2);
        let err = parse_earth(&text, 'D').expect_err("declared 3 but 2 rows must error");
        assert!(
            matches!(err, XtaskError::MalformedSource(_)),
            "expected MalformedSource for term-count mismatch, got {err:?}"
        );
    }

    // ==================================================================
    // (5) packed 往復: parse_earth → pack_earth → unpack_earth が元と完全一致。
    // ==================================================================
    #[test]
    fn pack_unpack_round_trips_full_model() {
        let earth = parse_earth(VSOP87D_EAR, 'D').expect("real VSOP87D.ear parses");
        let restored = unpack_earth(&pack_earth(&earth)).expect("packed model round-trips");
        assert_eq!(restored, earth, "round-tripped model must equal original");
    }

    // ==================================================================
    // (6) pack 決定性・長さ。
    //     レイアウト: [n_sections, <各 = variable, power, n_terms, <各項 amp,phase,freq>>...]。
    //     f64 数 = 1 + Σ_sections(3 + 3*n_terms) = 1 + 17*3 + 3*2425 = 7327。
    //     バイト長 = 8 * 7327 = 58616。
    // ==================================================================
    #[test]
    fn pack_earth_is_deterministic_and_has_exact_length() {
        let earth = parse_earth(VSOP87D_EAR, 'D').expect("real VSOP87D.ear parses");
        let a = pack_earth(&earth);
        let b = pack_earth(&earth);
        assert_eq!(a, b, "pack_earth must be deterministic");

        // 1 (n_sections) + Σ(3 ヘッダ + 3*n_terms) = 1 + 17*3 + 3*TOTAL_TERMS。
        let expected_f64 = 1 + DECLARED_SECTIONS.len() * 3 + 3 * TOTAL_TERMS;
        assert_eq!(expected_f64, 7327, "f64 count = 1 + 17*3 + 3*2425");
        assert_eq!(
            a.len(),
            8 * expected_f64,
            "packed length = 8 * (1 + Σ(3 + 3*n_terms))"
        );
    }

    // ==================================================================
    // (7) build_artifact（純関数）: packed バイト・DataSetMetadata の整合。
    // ==================================================================
    #[test]
    fn build_artifact_checksum_and_metadata_are_consistent() {
        let (bytes, metadata) = build_artifact(VSOP87D_EAR, 'D').expect("artifact builds");
        // packed バイトは pack_earth と同一（純関数の決定性）。
        let earth = parse_earth(VSOP87D_EAR, 'D').expect("real VSOP87D.ear parses");
        assert_eq!(bytes, pack_earth(&earth), "artifact bytes == pack_earth");
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
        // 版識別子は D 系（"D" を含む）、名称は earth/vsop87 を示す。
        assert!(
            metadata.version.contains('D'),
            "version should denote VSOP87 version D: {:?}",
            metadata.version
        );
        let name_lower = metadata.name.to_ascii_lowercase();
        assert!(
            name_lower.contains("vsop87") || name_lower.contains("earth"),
            "name should reference vsop87/earth: {:?}",
            metadata.name
        );
    }

    /// build_artifact は決定的（同一入力 → 同一バイト・同一 checksum）。
    #[test]
    fn build_artifact_is_deterministic() {
        let (bytes_a, meta_a) = build_artifact(VSOP87D_EAR, 'D').expect("artifact builds");
        let (bytes_b, meta_b) = build_artifact(VSOP87D_EAR, 'D').expect("artifact builds");
        assert_eq!(bytes_a, bytes_b, "bytes deterministic");
        assert_eq!(meta_a.checksum, meta_b.checksum, "checksum deterministic");
    }

    // ==================================================================
    // (8) unpack 異常: 壊れた packed バイト列で Err（パニックしない）。
    //     長さ検査 / 件数整合のどちらに委ねるかは未確定なので緩く照合。
    // ==================================================================

    /// 8 の倍数でない長さ → 必ず Err（Malformed*）。
    #[test]
    fn unpack_rejects_non_f64_aligned_bytes() {
        let bytes = vec![0u8; 13];
        let err = unpack_earth(&bytes).expect_err("non-8-multiple must error");
        assert!(
            matches!(
                err,
                XtaskError::MalformedPacked(_) | XtaskError::MalformedSource(_)
            ),
            "expected Malformed*, got {err:?}"
        );
    }

    /// ヘッダで巨大なセクション数を宣言するがデータが伴わない（8 の倍数長は満たす）。
    /// [n_sections = 1e9] の 1 f64 のみ → 件数とデータ長が矛盾。
    #[test]
    fn unpack_rejects_header_data_length_mismatch() {
        let bytes = crate::packed::pack_f64_le(&[1.0e9_f64]);
        let err = unpack_earth(&bytes)
            .expect_err("header section count vs data length mismatch must error");
        assert!(
            matches!(
                err,
                XtaskError::MalformedPacked(_) | XtaskError::MalformedSource(_)
            ),
            "expected Malformed*, got {err:?}"
        );
    }

    // ==================================================================
    // 検証ヘルパ・判定の境界（ミューテーション堅牢化）。
    // ==================================================================

    #[test]
    fn f64_to_count_accepts_valid_and_boundary() {
        assert_eq!(f64_to_count(0.0).unwrap(), 0);
        assert_eq!(f64_to_count(2425.0).unwrap(), 2425);
        assert_eq!(f64_to_count(MAX_TERMS).unwrap(), 100_000); // 上限ちょうど受理
    }

    #[test]
    fn f64_to_count_rejects_invalid() {
        assert!(f64_to_count(MAX_TERMS + 1.0).is_err(), "above max");
        assert!(f64_to_count(-1.0).is_err(), "negative");
        assert!(f64_to_count(1.5).is_err(), "non-integer");
        assert!(f64_to_count(f64::NAN).is_err(), "NaN");
        assert!(f64_to_count(f64::INFINITY).is_err(), "infinity");
    }

    #[test]
    fn f64_to_u8_accepts_valid_and_boundary() {
        assert_eq!(f64_to_u8(0.0).unwrap(), 0);
        assert_eq!(f64_to_u8(3.0).unwrap(), 3);
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

    /// is_header は VARIABLE と TERMS の**両方**を含む行のみ見出しと判定する（`&&`→`||` を殺す）。
    #[test]
    fn is_header_requires_both_variable_and_terms() {
        assert!(is_header(
            " VSOP87 VERSION D4  EARTH  VARIABLE 1 (LBR)  *T**0  559 TERMS  ..."
        ));
        assert!(!is_header(
            " something with VARIABLE but no t-e-r-m-s keyword"
        ));
        assert!(!is_header(" 5 TERMS but no var-keyword here"));
        assert!(!is_header(" 4310    1  0  0  0.0  1.0  2.0")); // 項行
    }

    /// parse_term_line: 末尾3トークン=A,B,C。3 トークン未満は None（`<3` 境界を殺す）。
    #[test]
    fn parse_term_line_edge_cases() {
        // 現実的な項行（先頭 id + 多数トークン、末尾3=A,B,C）。
        let t =
            parse_term_line(" 4310 2 0 0 1 0 0 0 0 0 0 0 0 0 -0.1 -0.2 0.5 1.5 6283.0").unwrap();
        assert_eq!((t.amplitude, t.phase, t.frequency), (0.5, 1.5, 6283.0));
        // ちょうど3トークン（id=10, B/C…ここでは A=10,B=1.5,C=2.5）は Some（`<=3` を殺す）。
        assert!(parse_term_line("10 1.5 2.5").is_some());
        // 2 トークンは None（`==3` 改変はここで underflow→panic で捕捉、`<` は None）。
        assert!(parse_term_line("10 1.5").is_none());
        // 先頭が非整数の行（見出し等）は None。
        assert!(parse_term_line("VSOP87 VERSION D4 1.0 2.0 3.0").is_none());
    }

    /// render_metadata は全フィールドを行として出力する（本体→"" 変異を殺す）。
    #[test]
    fn render_metadata_includes_all_fields() {
        let (_, metadata) = build_artifact(VSOP87D_EAR, 'D').unwrap();
        let rendered = render_metadata(&metadata);
        assert!(rendered.contains("name = vsop87-earth"), "{rendered}");
        assert!(rendered.contains("version = VSOP87D"), "{rendered}");
        assert!(
            rendered.contains(&format!("checksum = {}", metadata.checksum)),
            "{rendered}"
        );
    }

    /// 最終セクション R T5（series[16]）の構造を独立に確認（順序付けの終端）。
    #[test]
    fn last_section_is_r_t5() {
        let earth = parse_earth(VSOP87D_EAR, 'D').expect("real VSOP87D.ear parses");
        let last = earth.series.last().unwrap();
        assert_eq!(last.variable, 3, "last section variable = R(3)");
        assert_eq!(last.power, 5, "last section power = T^5");
        assert_eq!(last.terms.len(), 3, "R T5 = 3 terms");
    }

    // ==================================================================
    // 版A（黄道 J2000 直交 X,Y,Z）— M10 太陽フレーム修正（ISSUE-013/035）。
    //   構造は版D と同形式（A,B,C 項）。版判別・XYZ 18 セクション・round-trip を検証。
    //   VSOP87A.ear ヘッダ宣言項数（PROVENANCE.md）:
    //     X(1): 843,491,204,18,15,6  Y(2): 854,496,202,17,15,6  Z(3): 178,120,53,12,6,2
    // ==================================================================

    /// 版A・地球の宣言項数表（出現順 X T0..T5, Y T0..T5, Z T0..T5 = 18 セクション）。
    const DECLARED_SECTIONS_A: &[(u8, u8, usize)] = &[
        (1, 0, 843),
        (1, 1, 491),
        (1, 2, 204),
        (1, 3, 18),
        (1, 4, 15),
        (1, 5, 6),
        (2, 0, 854),
        (2, 1, 496),
        (2, 2, 202),
        (2, 3, 17),
        (2, 4, 15),
        (2, 5, 6),
        (3, 0, 178),
        (3, 1, 120),
        (3, 2, 53),
        (3, 3, 12),
        (3, 4, 6),
        (3, 5, 2),
    ];

    /// 版A は 18 セクション、各 (variable, power, terms.len()) が宣言値と出現順で一致。
    #[test]
    fn parse_earth_version_a_has_18_sections_in_declared_order() {
        let earth = parse_earth(VSOP87A_EAR, 'A').expect("real VSOP87A.ear parses");
        assert_eq!(earth.series.len(), 18, "X 6 + Y 6 + Z 6 = 18 sections");
        for (got, &(var, pow, n)) in earth.series.iter().zip(DECLARED_SECTIONS_A) {
            assert_eq!(got.variable, var, "variable order");
            assert_eq!(got.power, pow, "power order");
            assert_eq!(got.terms.len(), n, "term count (var={var}, T**{pow})");
        }
    }

    /// 版A の packed round-trip（parse → pack → unpack が元と完全一致）。
    #[test]
    fn version_a_packed_roundtrips() {
        let earth = parse_earth(VSOP87A_EAR, 'A').expect("real VSOP87A.ear parses");
        let bytes = pack_earth(&earth);
        let restored = unpack_earth(&bytes).expect("unpack");
        assert_eq!(restored, earth, "round-trip identity");
    }

    /// 版A の build_artifact メタデータ（version=VSOP87A・source に version A 明記）。
    #[test]
    fn build_artifact_version_a_metadata() {
        let (bytes, metadata) = build_artifact(VSOP87A_EAR, 'A').expect("artifact builds");
        assert_eq!(metadata.version, "VSOP87A");
        assert!(
            metadata.source.contains("VSOP87A.ear") && metadata.source.contains("version A"),
            "source names version A: {:?}",
            metadata.source
        );
        assert_eq!(
            metadata.checksum,
            sha256_hex(&bytes),
            "checksum matches bytes"
        );
    }

    /// 版判別: 版D を 'A' で、版A を 'D' で読むと MalformedSource（取り違え検出）。
    #[test]
    fn version_letter_mismatch_is_rejected() {
        assert!(matches!(
            parse_earth(VSOP87D_EAR, 'A'),
            Err(XtaskError::MalformedSource(_))
        ));
        assert!(matches!(
            parse_earth(VSOP87A_EAR, 'D'),
            Err(XtaskError::MalformedSource(_))
        ));
    }
}
