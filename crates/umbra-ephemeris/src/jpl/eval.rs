//! SPK type 2（Chebyshev 位置）評価（ISSUE-036 S2・feature `jpl`）。
//!
//! type 2 セグメントは末尾 4 double のディレクトリ [INIT, INTLEN, RSIZE, N] と、その手前に並ぶ
//! N 個の論理レコード（各 RSIZE double）から成る。各レコードは [MID, RADIUS, Cx[ncoeff],
//! Cy[ncoeff], Cz[ncoeff]]（RSIZE = 2 + 3·ncoeff）。ET から被覆レコードを特定し、正規化時刻
//! τ=(et−MID)/RADIUS ∈ [−1,1] で Chebyshev 級数を評価して位置を、その微分/RADIUS で速度を得る。
//!
//! 仕様: NAIF SPK Required Reading（type 2）。アドレスは 1 始まり DAF ワード（8 バイト）。
//! 本スライスは type 2（位置）のみ。原点/フレームはセグメント native（center/frame）のまま返す
//! （body 差・フレーム変換は S3 の `JplEphemeris` 層）。

// S2（評価）は単体で完結し、非テストビルドからの呼び出しは S3（JplEphemeris）で配線される。
// それまで crate 内利用者が無く dead_code 警告となるためモジュール単位で許可する（配線後に解除）。
#![allow(dead_code)]

use crate::ephemeris::EphemerisError;
use crate::jpl::daf::SpkSegment;

/// DAF ワード長（バイト）= f64 1 個。
const WORD_BYTES: usize = 8;
/// type 2 ディレクトリの double 数（末尾 [INIT, INTLEN, RSIZE, N]）。
const DIRECTORY_DOUBLES: usize = 4;
/// 1 レコードの非係数 double 数（MID, RADIUS）。
const RECORD_HEADER_DOUBLES: usize = 2;

/// セグメント native（center 原点・frame）での状態。位置[km]・速度[km/s]。
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct SegmentState {
    /// 位置 \[km\]（X, Y, Z）。
    pub position: [f64; 3],
    /// 速度 \[km/s\]（dX, dY, dZ）。type 2 は Chebyshev 微分から導出。
    pub velocity: [f64; 3],
}

fn malformed(msg: impl Into<String>) -> EphemerisError {
    EphemerisError::MalformedSpk(msg.into())
}

/// 1 始まり DAF ワード番地 `word` の LE f64 を読む。範囲外は MalformedSpk。
fn read_word(bytes: &[u8], word: usize) -> Result<f64, EphemerisError> {
    if word < 1 {
        return Err(malformed(format!("DAF word address {word} < 1")));
    }
    let off = (word - 1) * WORD_BYTES;
    let slice = bytes
        .get(off..off + WORD_BYTES)
        .ok_or_else(|| malformed(format!("DAF word {word} out of range (byte {off})")))?;
    Ok(f64::from_le_bytes(
        slice.try_into().expect("slice is 8 bytes"),
    ))
}

/// 整数値として格納された double（RSIZE/N）を正の usize へ。非有限・非整数・<1 は MalformedSpk。
#[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
fn as_positive_count(value: f64, what: &str) -> Result<usize, EphemerisError> {
    if !value.is_finite() || value < 1.0 || value.fract() != 0.0 {
        return Err(malformed(format!(
            "{what} is not a positive integer: {value}"
        )));
    }
    Ok(value as usize)
}

/// Chebyshev 級数 Σ C_k T_k(τ) と微分 Σ C_k T_k'(τ) を返す（第 1 種, T_0=1,T_1=τ）。
/// 微分は T_k' = k·U_{k-1}（U_0=1, U_1=2τ, U_{k+1}=2τU_k−U_{k-1}）。
#[allow(clippy::cast_precision_loss)]
fn chebyshev_value_and_deriv(coeffs: &[f64], tau: f64) -> (f64, f64) {
    if coeffs.is_empty() {
        return (0.0, 0.0);
    }
    // 値: T 漸化。
    let mut t_km1 = 1.0; // T_0
    let mut t_k = tau; // T_1
    let mut value = coeffs[0]; // C_0·T_0
    if coeffs.len() > 1 {
        value += coeffs[1] * t_k;
    }
    // ncoeff<2 のとき skip(2) は空（添字スライス [2..] は len<2 でパニックするため使わない）。
    for &c in coeffs.iter().skip(2) {
        let t_kp1 = 2.0 * tau * t_k - t_km1;
        value += c * t_kp1;
        t_km1 = t_k;
        t_k = t_kp1;
    }
    // 微分: U 漸化（T_k' = k·U_{k-1}, T_0'=0）。
    let mut deriv = 0.0;
    let mut u_km1 = 0.0; // U_{-1}
    let mut u_k = 1.0; // U_0（k=1 のとき U_{k-1}=U_0）
    for (k, &c) in coeffs.iter().enumerate() {
        if k == 0 {
            continue; // T_0' = 0
        }
        deriv += c * (k as f64) * u_k; // u_k が U_{k-1} を保持
        let u_kp1 = 2.0 * tau * u_k - u_km1;
        u_km1 = u_k;
        u_k = u_kp1;
    }
    (value, deriv)
}

/// SPK type 2 セグメントを `et`（TDB 秒, J2000 基準）で評価する（S2: Chebyshev 評価）。
///
/// `segment.data_type` が 2 でない場合・`et` がセグメント被覆外・データ構造不整合は
/// [`EphemerisError`]（範囲外は `OutOfSupportedRange`, 構造不整合は `MalformedSpk`）。
#[allow(
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::cast_precision_loss
)]
pub(crate) fn eval_type2(
    bytes: &[u8],
    segment: &SpkSegment,
    et: f64,
) -> Result<SegmentState, EphemerisError> {
    // 本スライスは type 2（Chebyshev 位置）のみ対応。
    if segment.data_type != 2 {
        return Err(malformed(format!(
            "unsupported SPK data type {} (eval_type2 handles type 2 only)",
            segment.data_type
        )));
    }

    let start_addr = usize::try_from(segment.start_addr)
        .map_err(|_| malformed(format!("negative start_addr {}", segment.start_addr)))?;
    let end_addr = usize::try_from(segment.end_addr)
        .map_err(|_| malformed(format!("negative end_addr {}", segment.end_addr)))?;
    if start_addr < 1 || end_addr < start_addr + DIRECTORY_DOUBLES {
        return Err(malformed(format!(
            "segment address range too small (start={start_addr}, end={end_addr})"
        )));
    }

    // 末尾 4 double のディレクトリ [INIT, INTLEN, RSIZE, N]（ワード end_addr-3..end_addr）。
    let init = read_word(bytes, end_addr - 3)?;
    let intlen = read_word(bytes, end_addr - 2)?;
    let rsize = as_positive_count(read_word(bytes, end_addr - 1)?, "RSIZE")?;
    let n = as_positive_count(read_word(bytes, end_addr)?, "N")?;
    if !intlen.is_finite() || intlen <= 0.0 {
        return Err(malformed(format!("non-positive INTLEN {intlen}")));
    }
    // RSIZE = 2 + 3·ncoeff（X/Y/Z 同数係数）。
    if rsize < RECORD_HEADER_DOUBLES + 3 || (rsize - RECORD_HEADER_DOUBLES) % 3 != 0 {
        return Err(malformed(format!("RSIZE {rsize} is not 2 + 3·ncoeff")));
    }
    let ncoeff = (rsize - RECORD_HEADER_DOUBLES) / 3;

    // N・RSIZE と宣言アドレス範囲の整合（破損 N による隣接領域の誤読を防ぐ）。
    // データ領域 = N レコード（各 RSIZE）+ ディレクトリ 4 = N·RSIZE + 4 ワード。
    // よって end_addr = start_addr + N·RSIZE + 4 − 1。checked 演算でオーバーフローも弾く
    // （以降の record_word = start_addr + m·RSIZE は m<N ゆえ本検証通過後は安全）。
    let expected_end = n
        .checked_mul(rsize)
        .and_then(|x| x.checked_add(start_addr))
        .and_then(|x| x.checked_add(DIRECTORY_DOUBLES))
        .and_then(|x| x.checked_sub(1))
        .ok_or_else(|| malformed("segment size arithmetic overflow"))?;
    if end_addr != expected_end {
        return Err(malformed(format!(
            "segment size mismatch: end_addr={end_addr} but N={n}·RSIZE={rsize}+{DIRECTORY_DOUBLES} implies {expected_end}"
        )));
    }

    // 被覆区間 [INIT, INIT + N·INTLEN]。外は範囲外。
    let coverage_end = init + (n as f64) * intlen;
    if et < init || et > coverage_end {
        return Err(EphemerisError::OutOfSupportedRange);
    }

    // レコード索引 m = floor((et−INIT)/INTLEN) を [0, N−1] にクランプ（被覆末尾は m=N−1）。
    let m_f = ((et - init) / intlen).floor();
    let m = m_f.clamp(0.0, (n - 1) as f64) as usize;

    // レコード m の先頭ワード（MID）。各レコードは RSIZE ワード。
    let record_word = start_addr + m * rsize;
    let mid = read_word(bytes, record_word)?;
    let radius = read_word(bytes, record_word + 1)?;
    if !radius.is_finite() || radius <= 0.0 {
        return Err(malformed(format!("non-positive RADIUS {radius}")));
    }

    // 係数: X[ncoeff], Y[ncoeff], Z[ncoeff]（MID, RADIUS の直後）。
    let coeff_word0 = record_word + RECORD_HEADER_DOUBLES;
    let tau = (et - mid) / radius;
    let mut position = [0.0_f64; 3];
    let mut velocity = [0.0_f64; 3];
    let mut coeffs = vec![0.0_f64; ncoeff];
    for axis in 0..3 {
        let base = coeff_word0 + axis * ncoeff;
        for (j, c) in coeffs.iter_mut().enumerate() {
            *c = read_word(bytes, base + j)?;
        }
        let (value, deriv) = chebyshev_value_and_deriv(&coeffs, tau);
        position[axis] = value;
        // dX/dt = (dX/dτ)·(dτ/dt) = deriv / radius。
        velocity[axis] = deriv / radius;
    }

    Ok(SegmentState { position, velocity })
}

#[cfg(test)]
mod tests {
    // 位置/速度と SPICE 基準を並列添字で比較する（st.position[i] と want[i] を同 i で参照）。
    #![allow(clippy::needless_range_loop)]

    use super::*;
    use crate::jpl::daf::{parse_spk_segments, SpkSegment};

    // ------------------------------------------------------------------
    // 仕様出典: NAIF "SPK Required Reading"（type 2: Chebyshev 位置）。
    //   type 2 セグメントのデータ領域（DAF ワード start_addr..=end_addr）は
    //     [ レコード0 ][ レコード1 ]…[ レコードN-1 ][ INIT, INTLEN, RSIZE, N ]
    //   末尾 4 double がディレクトリ。各レコード（RSIZE double）:
    //     [ MID, RADIUS, Cx[0..ncoeff], Cy[0..ncoeff], Cz[0..ncoeff] ]
    //   RSIZE = 2 + 3*ncoeff。被覆区間 = [INIT, INIT + N*INTLEN]。
    //   レコード索引 m = clamp(floor((et-INIT)/INTLEN), 0, N-1)。
    //   正規化時刻 τ = (et - MID)/RADIUS ∈ [-1,1]。
    //   X(τ)=Σ C_k T_k(τ)（T_0=1,T_1=τ,T_{k+1}=2τT_k-T_{k-1}）。位置[km]。
    //   速度 dX/dt = (Σ C_k T_k'(τ))/RADIUS [km/s]（T_k'=k·U_{k-1}）。
    //
    //   評価器の契約: bytes を DAF ワードアドレスで参照し、データ領域
    //   （ワード start_addr..=end_addr）のみを読む。ワード w（1 始まり）は
    //   バイト offset (w-1)*8 の LE f64。よって合成テストはデータ領域だけを
    //   正しいアドレスに置いた Vec<u8> を作れば良い（有効なファイルレコード不要）。
    // ------------------------------------------------------------------

    /// f64 1 個分のバイト数（= DAF ワード長）。
    const WORD: usize = 8;

    /// 許容誤差つき比較（数値丸めのみを許す厳密寄りの相対/絶対許容）。
    /// `|a-b| <= tol*(1+|b|)` で 0 近傍は絶対、大きい値は相対に振る舞う。
    fn approx(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() <= tol * (1.0 + b.abs())
    }

    /// 位置・速度の各成分を許容誤差つきで照合する。
    fn assert_state_close(actual: &SegmentState, pos: [f64; 3], vel: [f64; 3], tol: f64) {
        for i in 0..3 {
            assert!(
                approx(actual.position[i], pos[i], tol),
                "position[{i}]: got {} want {} (tol {tol})",
                actual.position[i],
                pos[i]
            );
            assert!(
                approx(actual.velocity[i], vel[i], tol),
                "velocity[{i}]: got {} want {} (tol {tol})",
                actual.velocity[i],
                vel[i]
            );
        }
    }

    /// 1 レコードの記述（合成セグメント組立用）。X/Y/Z の Chebyshev 係数は
    /// 各成分が同数（ncoeff）の係数を持つ。
    #[derive(Clone)]
    struct Record {
        mid: f64,
        radius: f64,
        cx: Vec<f64>,
        cy: Vec<f64>,
        cz: Vec<f64>,
    }

    /// type 2 セグメントを合成して `(bytes, segment)` を返す。
    ///
    /// データ領域を `start_addr`（1 始まり DAF ワード番地）に置き、その手前は
    /// 任意のパディングで埋める（評価器がデータ領域外を読まないことの検証にも使う）。
    /// 全レコードの ncoeff（= cx/cy/cz の長さ）は一致している前提。
    fn build_segment(
        init: f64,
        intlen: f64,
        records: &[Record],
        start_addr: i32,
        data_type: i32,
    ) -> (Vec<u8>, SpkSegment) {
        let ncoeff = records[0].cx.len();
        // 各成分の係数数は一致している前提（RSIZE = 2 + 3*ncoeff）。
        for r in records {
            assert_eq!(r.cx.len(), ncoeff, "cx の ncoeff 不一致");
            assert_eq!(r.cy.len(), ncoeff, "cy の ncoeff 不一致");
            assert_eq!(r.cz.len(), ncoeff, "cz の ncoeff 不一致");
        }
        let rsize = 2 + 3 * ncoeff; // double 数
        let n = records.len();
        // データ領域の総 double 数 = N*RSIZE + 4（末尾ディレクトリ）。
        let data_doubles = n * rsize + 4;

        // データ領域を double 列として組み立てる。
        let mut data: Vec<f64> = Vec::with_capacity(data_doubles);
        for r in records {
            data.push(r.mid);
            data.push(r.radius);
            data.extend_from_slice(&r.cx);
            data.extend_from_slice(&r.cy);
            data.extend_from_slice(&r.cz);
        }
        // ディレクトリ [INIT, INTLEN, RSIZE, N]。
        data.push(init);
        data.push(intlen);
        #[allow(clippy::cast_precision_loss)]
        {
            data.push(rsize as f64);
            data.push(n as f64);
        }
        assert_eq!(data.len(), data_doubles);

        // start_addr（1 始まりワード）→ データ先頭のバイト offset = (start_addr-1)*8。
        let start_word = usize::try_from(start_addr).expect("start_addr 非負");
        let data_byte_off = (start_word - 1) * WORD;
        let total_bytes = data_byte_off + data_doubles * WORD;
        let mut bytes = vec![0u8; total_bytes];
        // データ領域手前のパディングは「読まれてはいけない」毒値で埋める。
        for b in bytes.iter_mut().take(data_byte_off) {
            *b = 0xAB;
        }
        for (i, v) in data.iter().enumerate() {
            let off = data_byte_off + i * WORD;
            bytes[off..off + WORD].copy_from_slice(&v.to_le_bytes());
        }

        // end_addr = データ末尾ワード（1 始まり, 含む）。
        let end_word = start_word + data_doubles - 1;
        let segment = SpkSegment {
            start_et: init,
            #[allow(clippy::cast_precision_loss)]
            end_et: init + (n as f64) * intlen,
            target: 10,
            center: 0,
            frame: 1,
            data_type,
            start_addr,
            end_addr: i32::try_from(end_word).expect("end_addr 範囲内"),
        };
        (bytes, segment)
    }

    /// 単一レコードの type 2 セグメント（data_type=2, start_addr=1）を合成する近道。
    fn single_record_segment(
        init: f64,
        intlen: f64,
        mid: f64,
        radius: f64,
        cx: Vec<f64>,
        cy: Vec<f64>,
        cz: Vec<f64>,
    ) -> (Vec<u8>, SpkSegment) {
        let rec = Record {
            mid,
            radius,
            cx,
            cy,
            cz,
        };
        build_segment(init, intlen, &[rec], 1, 2)
    }

    /// 第 1 種 Chebyshev 多項式 T_k(τ) を漸化で評価（テスト独立実装）。
    fn cheb_t(coeffs: &[f64], tau: f64) -> f64 {
        if coeffs.is_empty() {
            return 0.0;
        }
        let mut tkm1 = 1.0; // T_0
        let mut tk = tau; // T_1
        let mut sum = coeffs[0] * tkm1;
        if coeffs.len() > 1 {
            sum += coeffs[1] * tk;
        }
        for &c in &coeffs[2..] {
            let tkp1 = 2.0 * tau * tk - tkm1;
            sum += c * tkp1;
            tkm1 = tk;
            tk = tkp1;
        }
        sum
    }

    /// dX/dτ = Σ C_k T_k'(τ)（T_k' = k·U_{k-1}, U_0=1, U_1=2τ, …）。テスト独立実装。
    fn cheb_dt(coeffs: &[f64], tau: f64) -> f64 {
        // U_{-1}=0, U_0=1 とすると T_k' = k·U_{k-1}。
        let mut sum = 0.0;
        if coeffs.len() <= 1 {
            return 0.0; // T_0'=0
        }
        // U 多項式列を回しながら T_k'（k>=1）を加算。
        let mut ukm1 = 0.0; // U_{-1}
        let mut uk = 1.0; // U_0
        for (k, &c) in coeffs.iter().enumerate() {
            if k == 0 {
                continue; // T_0' = 0
            }
            #[allow(clippy::cast_precision_loss)]
            let tkp = (k as f64) * uk; // T_k' = k·U_{k-1}（uk が U_{k-1} を保持）
            sum += c * tkp;
            let ukp1 = 2.0 * tau * uk - ukm1;
            ukm1 = uk;
            uk = ukp1;
        }
        sum
    }

    // ==================================================================
    // A. 合成セグメント（厳密値・外部依存なし）
    // ==================================================================

    /// 観点: 定数係数 Cx=[c,0,…] は全 τ で X=c・速度0。X/Y/Z 別々の定数も成立。
    #[test]
    fn constant_record_returns_constant_position_and_zero_velocity() {
        let init = 100.0;
        let intlen = 20.0;
        let mid = 110.0;
        let radius = 10.0; // 区間 [100, 120]
                           // ncoeff=2（2 個目は 0）で定数を表現。X/Y/Z で異なる定数。
        let (bytes, seg) = single_record_segment(
            init,
            intlen,
            mid,
            radius,
            vec![3.5, 0.0],
            vec![-7.0, 0.0],
            vec![42.0, 0.0],
        );
        // 区間内の複数の et で常に同じ位置・速度0。
        for et in [100.0, 105.0, 110.0, 119.999] {
            let st = eval_type2(&bytes, &seg, et).expect("被覆内 et は成功すべき");
            assert_state_close(&st, [3.5, -7.0, 42.0], [0.0, 0.0, 0.0], 1e-12);
        }
    }

    /// 観点: ncoeff=1（係数 1 個のみ）でも定数として評価できる（RSIZE=2+3=5 の最小レコード）。
    #[test]
    fn single_coefficient_is_constant() {
        let (bytes, seg) =
            single_record_segment(0.0, 10.0, 5.0, 5.0, vec![11.0], vec![22.0], vec![33.0]);
        let st = eval_type2(&bytes, &seg, 3.0).expect("被覆内 et は成功すべき");
        assert_state_close(&st, [11.0, 22.0, 33.0], [0.0, 0.0, 0.0], 1e-12);
    }

    /// 観点: 1 次 Cx=[c0,c1] → X=c0+c1·τ・速度=c1/RADIUS。τ を手計算して検算。
    #[test]
    fn linear_record_matches_hand_computed_value() {
        let init = 0.0;
        let intlen = 100.0;
        let mid = 50.0;
        let radius = 50.0; // 区間 [0,100]
        let (c0, c1) = (10.0, 4.0);
        let (bytes, seg) = single_record_segment(
            init,
            intlen,
            mid,
            radius,
            vec![c0, c1],
            vec![0.0, 0.0],
            vec![0.0, 0.0],
        );
        let et = 75.0; // τ = (75-50)/50 = 0.5
        let tau = (et - mid) / radius;
        let want_x = c0 + c1 * tau; // 10 + 4*0.5 = 12
        let want_vx = c1 / radius; // 4/50 = 0.08
        let st = eval_type2(&bytes, &seg, et).expect("被覆内 et は成功すべき");
        assert_state_close(&st, [want_x, 0.0, 0.0], [want_vx, 0.0, 0.0], 1e-12);
    }

    /// 観点: 2 次 Cx=[c0,c1,c2] → X=c0+c1·τ+c2·(2τ²−1)・速度=(c1+4c2·τ)/RADIUS。
    /// X,Y,Z で異なる係数を与え、独立評価（テスト側の cheb_t/cheb_dt）と一致すること。
    #[test]
    fn quadratic_record_matches_independent_evaluation() {
        let init = -50.0;
        let intlen = 100.0;
        let mid = 0.0;
        let radius = 50.0; // 区間 [-50,50]
        let cx = vec![1.0, 2.0, 3.0];
        let cy = vec![-4.0, 0.5, -1.5];
        let cz = vec![7.0, -2.0, 0.25];
        let (bytes, seg) = single_record_segment(
            init,
            intlen,
            mid,
            radius,
            cx.clone(),
            cy.clone(),
            cz.clone(),
        );
        let et = 12.5; // τ = 0.25
        let tau = (et - mid) / radius;
        // 期待位置（明示の閉形式）と独立評価の両方で検算。
        let want_x = cx[0] + cx[1] * tau + cx[2] * (2.0 * tau * tau - 1.0);
        let want_vx = (cx[1] + 4.0 * cx[2] * tau) / radius;
        let want_pos = [cheb_t(&cx, tau), cheb_t(&cy, tau), cheb_t(&cz, tau)];
        let want_vel = [
            cheb_dt(&cx, tau) / radius,
            cheb_dt(&cy, tau) / radius,
            cheb_dt(&cz, tau) / radius,
        ];
        // 閉形式と独立評価が整合していることをテスト内で自己検証。
        assert!(approx(want_pos[0], want_x, 1e-12));
        assert!(approx(want_vel[0], want_vx, 1e-12));

        let st = eval_type2(&bytes, &seg, et).expect("被覆内 et は成功すべき");
        assert_state_close(&st, want_pos, want_vel, 1e-9);
    }

    /// 観点: 高次（5 次）係数でも独立評価と一致（Chebyshev 漸化の安定性）。
    #[test]
    fn high_degree_record_matches_independent_evaluation() {
        let init = 0.0;
        let intlen = 2.0;
        let mid = 1.0;
        let radius = 1.0; // 区間 [0,2]
        let cx = vec![0.3, -1.1, 0.7, 0.2, -0.05, 0.9];
        let cy = vec![1.0, 0.0, -2.0, 0.5, 0.1, -0.3];
        let cz = vec![-0.5, 0.25, 0.125, -0.0625, 0.03125, 0.015625];
        let (bytes, seg) = single_record_segment(
            init,
            intlen,
            mid,
            radius,
            cx.clone(),
            cy.clone(),
            cz.clone(),
        );
        let et = 1.7; // τ = 0.7
        let tau = (et - mid) / radius;
        let want_pos = [cheb_t(&cx, tau), cheb_t(&cy, tau), cheb_t(&cz, tau)];
        let want_vel = [
            cheb_dt(&cx, tau) / radius,
            cheb_dt(&cy, tau) / radius,
            cheb_dt(&cz, tau) / radius,
        ];
        let st = eval_type2(&bytes, &seg, et).expect("被覆内 et は成功すべき");
        assert_state_close(&st, want_pos, want_vel, 1e-9);
    }

    /// 観点: 複数レコード（N=3, INTLEN 区間）で各 et が正しいレコードを選択。
    /// レコード毎に異なる定数を持たせ、被覆内 et が属レコードの値を返すこと。
    #[test]
    fn selects_correct_record_among_multiple() {
        let init = 0.0;
        let intlen = 10.0; // レコード0:[0,10], 1:[10,20], 2:[20,30]
                           // 各レコードは自区間の中点・半長を持つ。定数係数で識別。
        let records = vec![
            Record {
                mid: 5.0,
                radius: 5.0,
                cx: vec![100.0, 0.0],
                cy: vec![0.0, 0.0],
                cz: vec![0.0, 0.0],
            },
            Record {
                mid: 15.0,
                radius: 5.0,
                cx: vec![200.0, 0.0],
                cy: vec![0.0, 0.0],
                cz: vec![0.0, 0.0],
            },
            Record {
                mid: 25.0,
                radius: 5.0,
                cx: vec![300.0, 0.0],
                cy: vec![0.0, 0.0],
                cz: vec![0.0, 0.0],
            },
        ];
        let (bytes, seg) = build_segment(init, intlen, &records, 1, 2);

        // 各区間内部の et → 対応する定数。
        let cases = [(3.0, 100.0), (12.0, 200.0), (27.0, 300.0)];
        for (et, want) in cases {
            let st = eval_type2(&bytes, &seg, et).expect("被覆内 et は成功すべき");
            assert_state_close(&st, [want, 0.0, 0.0], [0.0, 0.0, 0.0], 1e-12);
        }
    }

    /// 観点: レコード境界の連続性。隣接レコードの係数を「境界時刻で同じ位置」を返すよう
    /// 整合させ、その境界 et を評価したとき（どちらのレコードが選ばれても）同一位置になること。
    /// レコード0:[0,10] は X=t を、レコード1:[10,20] も X=t を表す 1 次多項式で構成する。
    #[test]
    fn record_boundary_is_continuous() {
        let init = 0.0;
        let intlen = 10.0;
        // X = et を表現する。レコード i では X(τ)=mid + radius·τ（τ=(et-mid)/radius）。
        // よって cx=[mid, radius]（1 次）。両レコードとも X=et を返すので境界で連続。
        let records = vec![
            Record {
                mid: 5.0,
                radius: 5.0,
                cx: vec![5.0, 5.0], // X = 5 + 5τ = et
                cy: vec![0.0, 0.0],
                cz: vec![0.0, 0.0],
            },
            Record {
                mid: 15.0,
                radius: 5.0,
                cx: vec![15.0, 5.0], // X = 15 + 5τ = et
                cy: vec![0.0, 0.0],
                cz: vec![0.0, 0.0],
            },
        ];
        let (bytes, seg) = build_segment(init, intlen, &records, 1, 2);

        // 境界 et=10.0 直前/直後で X≈et かつ連続（速度は 5/5=1）。
        let st_lo = eval_type2(&bytes, &seg, 9.999).expect("被覆内");
        let st_hi = eval_type2(&bytes, &seg, 10.001).expect("被覆内");
        assert!(approx(st_lo.position[0], 9.999, 1e-9));
        assert!(approx(st_hi.position[0], 10.001, 1e-9));
        // 境界点そのもの（どちらのレコードでも X=et）。
        let st_b = eval_type2(&bytes, &seg, 10.0).expect("被覆内");
        assert!(approx(st_b.position[0], 10.0, 1e-9));
        // 速度はどちらの側でも一致（1.0 km/s 相当）。
        assert!(approx(st_lo.velocity[0], 1.0, 1e-9));
        assert!(approx(st_hi.velocity[0], 1.0, 1e-9));
    }

    /// 観点: start_addr のオフセット。データ領域を非ゼロ番地（先頭に毒値パディング）に
    /// 置いても正しく読む（DAF ワードアドレス計算の検証）。
    #[test]
    fn nonzero_start_addr_is_read_correctly() {
        let init = 0.0;
        let intlen = 100.0;
        let mid = 50.0;
        let radius = 50.0;
        let rec = Record {
            mid,
            radius,
            cx: vec![10.0, 4.0],
            cy: vec![1.0, 0.0],
            cz: vec![2.0, 0.0],
        };
        // start_addr=129（先頭 128 ワード=1024 バイトのパディングを毒値で置く）。
        let (bytes, seg) = build_segment(init, intlen, &[rec], 129, 2);
        let et = 75.0; // τ=0.5
        let st = eval_type2(&bytes, &seg, et).expect("被覆内 et は成功すべき");
        // 先頭の毒値を読んでいれば値は壊れる。正しいデータ領域を読めば X=12。
        assert_state_close(&st, [12.0, 1.0, 2.0], [4.0 / 50.0, 0.0, 0.0], 1e-9);
    }

    /// 観点: 同一点を両側区間の端で評価して一致すること（被覆末尾 et=INIT+N·INTLEN は
    /// m=N-1 にクランプ）。被覆末尾の連続評価が破綻しないこと。
    #[test]
    fn coverage_end_clamps_to_last_record() {
        let init = 0.0;
        let intlen = 10.0;
        // 単一レコード [0,10]、X=et（cx=[5,5]）。被覆末尾 et=10.0 は m=0 にクランプ。
        let (bytes, seg) = single_record_segment(
            init,
            intlen,
            5.0,
            5.0,
            vec![5.0, 5.0],
            vec![0.0, 0.0],
            vec![0.0, 0.0],
        );
        let st = eval_type2(&bytes, &seg, 10.0).expect("被覆末尾 et は m=N-1 で成功すべき");
        assert!(approx(st.position[0], 10.0, 1e-9), "被覆末尾でも X=et");
    }

    // ==================================================================
    // A. 異常系
    // ==================================================================

    /// 観点: data_type≠2（例 3）は MalformedSpk（メッセージ非依存・変種のみ）。
    #[test]
    fn wrong_data_type_is_malformed() {
        let (bytes, seg) = build_segment(
            0.0,
            100.0,
            &[Record {
                mid: 50.0,
                radius: 50.0,
                cx: vec![1.0, 0.0],
                cy: vec![0.0, 0.0],
                cz: vec![0.0, 0.0],
            }],
            1,
            3,
        ); // data_type=3
        let err = eval_type2(&bytes, &seg, 50.0).unwrap_err();
        assert!(matches!(err, EphemerisError::MalformedSpk(_)));
    }

    /// 観点: et が被覆開始未満（INIT 未満）は OutOfSupportedRange。
    #[test]
    fn et_before_coverage_is_out_of_range() {
        let (bytes, seg) = single_record_segment(
            0.0,
            100.0,
            50.0,
            50.0,
            vec![1.0, 0.0],
            vec![0.0, 0.0],
            vec![0.0, 0.0],
        );
        // start_et=0.0。少し手前。
        let err = eval_type2(&bytes, &seg, -1.0).unwrap_err();
        assert!(matches!(err, EphemerisError::OutOfSupportedRange));
    }

    /// 観点: et が被覆末尾超過（INIT+N·INTLEN 超）は OutOfSupportedRange。
    #[test]
    fn et_after_coverage_is_out_of_range() {
        let (bytes, seg) = single_record_segment(
            0.0,
            100.0,
            50.0,
            50.0,
            vec![1.0, 0.0],
            vec![0.0, 0.0],
            vec![0.0, 0.0],
        );
        // end_et=100.0。少し超過。
        let err = eval_type2(&bytes, &seg, 100.0 + 1.0).unwrap_err();
        assert!(matches!(err, EphemerisError::OutOfSupportedRange));
    }

    /// 1 始まり DAF ワード番地 `word` に LE f64 を書き込む（破損注入用）。
    fn put_word(bytes: &mut [u8], word: usize, v: f64) {
        let off = (word - 1) * WORD;
        bytes[off..off + WORD].copy_from_slice(&v.to_le_bytes());
    }

    /// 観点: ディレクトリの INTLEN が非正（0）なら MalformedSpk。
    /// 健全なセグメントを作り、ディレクトリの INTLEN（ワード end_addr-2）を 0.0 に上書きして
    /// 被覆内のはずだった et を評価しても、構造不整合として拒否されること（パニックしない）。
    #[test]
    fn rejects_nonpositive_intlen() {
        let (mut bytes, seg) = single_record_segment(
            0.0,
            100.0,
            50.0,
            50.0,
            vec![1.0, 0.0],
            vec![0.0, 0.0],
            vec![0.0, 0.0],
        );
        let end = usize::try_from(seg.end_addr).unwrap();
        // ディレクトリ [INIT, INTLEN, RSIZE, N] の INTLEN はワード end_addr-2。
        put_word(&mut bytes, end - 2, 0.0);
        // 健全なら被覆内だった et（初期値付近）でも MalformedSpk。
        let err = eval_type2(&bytes, &seg, 50.0).unwrap_err();
        assert!(matches!(err, EphemerisError::MalformedSpk(_)));
    }

    /// 観点: レコードの RADIUS が非正（0）なら MalformedSpk。
    /// 健全なセグメントを作り、レコード0 の RADIUS（ワード start_addr+1）を 0.0 に上書きして
    /// そのレコードが選ばれる et を評価しても、τ 正規化の分母不正として拒否されること。
    #[test]
    fn rejects_nonpositive_radius() {
        let (mut bytes, seg) = single_record_segment(
            0.0,
            100.0,
            50.0,
            50.0,
            vec![1.0, 0.0],
            vec![0.0, 0.0],
            vec![0.0, 0.0],
        );
        let start = usize::try_from(seg.start_addr).unwrap();
        // レコード0 = [MID, RADIUS, …]。RADIUS はワード start_addr+1。
        put_word(&mut bytes, start + 1, 0.0);
        // レコード0 が選ばれる被覆内 et でも MalformedSpk。
        let err = eval_type2(&bytes, &seg, 50.0).unwrap_err();
        assert!(matches!(err, EphemerisError::MalformedSpk(_)));
    }

    // ==================================================================
    // B. 実 DE440s 受入（ゲート・SPICE オラクル）
    //   data/spk/de440s.bsp 存在時のみ実行（CI 非同梱・ISSUE-036）。
    //   パス: CARGO_MANIFEST_DIR (crates/umbra-ephemeris)/../../data/spk/de440s.bsp。
    //
    // オラクル基準値の生成（出典・逐語転記）:
    //   spiceypy 8.1.2 / SPICE toolkit CSPICE_N0067（Docker python:3.12-slim, pip install spiceypy）
    //   import spiceypy as sp; sp.furnsh('de440s.bsp')
    //   sp.spkgeo(target, et, 'J2000', obs)  ->  (state[0:3]=位置km, state[3:6]=速度km/s, lt)
    //   spkgeo は body 差・J2000 フレームでセグメント native と一致する。
    //
    //   Sun(10) wrt SSB(0):  sp.spkgeo(10, et, 'J2000', 0)
    //     et=0.0:
    //       [-1067706.8053809535, -396036.18479594623, -138065.18428688092,
    //         0.009312571926520472, -0.01170150612817771, -0.005251266205200356]
    //     et=750000000.0:
    //       [-1249290.2370180925, -324838.0277918544, -106018.13729732925,
    //         0.007071092471312386, -0.012387616579323522, -0.005423569393921095]
    //     et=-1000000000.0:
    //       [571997.7425724803, -213216.87848410784, -98740.71783691116,
    //         0.005839325981523217, 0.007990608420087036, 0.0033003327192824682]
    //   Moon(301) wrt EMB(3):  sp.spkgeo(301, et, 'J2000', 3)
    //     et=0.0:
    //       [-288065.17234541546, -263476.06800028845, -75177.79740766216,
    //         0.6357121052811876, -0.6579943294710526, -0.29766442157325324]
    //     et=750000000.0:
    //       [-201891.21293137703, 298470.0858977215, 169054.0330397744,
    //        -0.8392452377606863, -0.42264621266089697, -0.19040420324384058]
    // ==================================================================

    /// DE440s から (target,center) 一致の先頭セグメントを取得する。不在なら None。
    fn find_segment(segs: &[SpkSegment], target: i32, center: i32) -> Option<SpkSegment> {
        segs.iter()
            .copied()
            .find(|s| s.target == target && s.center == center)
    }

    /// 実 DE440s の Sun(10) wrt SSB(0) セグメントを既知 ET で評価し、SPICE と突合する。
    /// 位置は数 m（相対 ~1e-10）、速度は mm/s 級まで一致すべき。
    #[test]
    fn de440s_sun_matches_spice_oracle() {
        let path = concat!(env!("CARGO_MANIFEST_DIR"), "/../../data/spk/de440s.bsp");
        let bytes = match std::fs::read(path) {
            Ok(b) => b,
            Err(e) => {
                eprintln!(
                    "skip de440s_sun_matches_spice_oracle: {path} を読めない（{e}）。\
                     実ファイルは CI 非同梱（ISSUE-036）。"
                );
                return;
            }
        };
        let segs = parse_spk_segments(&bytes).expect("実 DE440s は解析成功すべき");
        let seg = find_segment(&segs, 10, 0).expect("Sun(10) wrt SSB(0) セグメントがあるはず");

        // (et, SPICE state[6]) の逐語転記（出典は上記コメント）。
        let cases: [(f64, [f64; 6]); 3] = [
            (
                0.0,
                [
                    -1067706.8053809535,
                    -396036.18479594623,
                    -138065.18428688092,
                    0.009312571926520472,
                    -0.01170150612817771,
                    -0.005251266205200356,
                ],
            ),
            (
                750000000.0,
                [
                    -1249290.2370180925,
                    -324838.0277918544,
                    -106018.13729732925,
                    0.007071092471312386,
                    -0.012387616579323522,
                    -0.005423569393921095,
                ],
            ),
            (
                -1000000000.0,
                [
                    571997.7425724803,
                    -213216.87848410784,
                    -98740.71783691116,
                    0.005839325981523217,
                    0.007990608420087036,
                    0.0033003327192824682,
                ],
            ),
        ];

        for (et, want) in cases {
            let st = eval_type2(&bytes, &seg, et)
                .unwrap_or_else(|e| panic!("et={et} の評価に失敗: {e:?}"));
            // 位置: 絶対 1e-2 km = 10 m 以内。
            for i in 0..3 {
                assert!(
                    (st.position[i] - want[i]).abs() < 1.0e-2,
                    "et={et} position[{i}]: got {} want {} (Δ={:.3e} km)",
                    st.position[i],
                    want[i],
                    (st.position[i] - want[i]).abs()
                );
            }
            // 速度: 絶対 1e-6 km/s = 1 mm/s 以内。
            for i in 0..3 {
                assert!(
                    (st.velocity[i] - want[3 + i]).abs() < 1.0e-6,
                    "et={et} velocity[{i}]: got {} want {} (Δ={:.3e} km/s)",
                    st.velocity[i],
                    want[3 + i],
                    (st.velocity[i] - want[3 + i]).abs()
                );
            }
        }
    }

    /// 実 DE440s の Moon(301) wrt EMB(3) セグメントを既知 ET で評価し、SPICE と突合する。
    #[test]
    fn de440s_moon_matches_spice_oracle() {
        let path = concat!(env!("CARGO_MANIFEST_DIR"), "/../../data/spk/de440s.bsp");
        let bytes = match std::fs::read(path) {
            Ok(b) => b,
            Err(e) => {
                eprintln!(
                    "skip de440s_moon_matches_spice_oracle: {path} を読めない（{e}）。\
                     実ファイルは CI 非同梱（ISSUE-036）。"
                );
                return;
            }
        };
        let segs = parse_spk_segments(&bytes).expect("実 DE440s は解析成功すべき");
        let seg = find_segment(&segs, 301, 3).expect("Moon(301) wrt EMB(3) セグメントがあるはず");

        let cases: [(f64, [f64; 6]); 2] = [
            (
                0.0,
                [
                    -288065.17234541546,
                    -263476.06800028845,
                    -75177.79740766216,
                    0.6357121052811876,
                    -0.6579943294710526,
                    -0.29766442157325324,
                ],
            ),
            (
                750000000.0,
                [
                    -201891.21293137703,
                    298470.0858977215,
                    169054.0330397744,
                    -0.8392452377606863,
                    -0.42264621266089697,
                    -0.19040420324384058,
                ],
            ),
        ];

        for (et, want) in cases {
            let st = eval_type2(&bytes, &seg, et)
                .unwrap_or_else(|e| panic!("et={et} の評価に失敗: {e:?}"));
            for i in 0..3 {
                assert!(
                    (st.position[i] - want[i]).abs() < 1.0e-2,
                    "et={et} position[{i}]: got {} want {} (Δ={:.3e} km)",
                    st.position[i],
                    want[i],
                    (st.position[i] - want[i]).abs()
                );
            }
            for i in 0..3 {
                assert!(
                    (st.velocity[i] - want[3 + i]).abs() < 1.0e-6,
                    "et={et} velocity[{i}]: got {} want {} (Δ={:.3e} km/s)",
                    st.velocity[i],
                    want[3 + i],
                    (st.velocity[i] - want[3 + i]).abs()
                );
            }
        }
    }
}
