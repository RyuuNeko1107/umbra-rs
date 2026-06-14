//! 3×3 行列と座標回転（`docs/issues/ISSUE-003`、`docs/conventions.md` §5）。
//!
//! 回転は **SOFA 準拠の座標（受動）回転** R1/R2/R3。座標系を角度 θ だけ回す向きで、
//! フレーム変換（GCRS→CIRS→TIRS→ITRS 等、`docs/algorithms/02-frames.md`）に用いる。
//!
//! ```text
//! R1(θ) = [[1, 0, 0], [0, cosθ, sinθ], [0, -sinθ, cosθ]]
//! R2(θ) = [[cosθ, 0, -sinθ], [0, 1, 0], [sinθ, 0, cosθ]]
//! R3(θ) = [[cosθ, sinθ, 0], [-sinθ, cosθ, 0], [0, 0, 1]]
//! ```

use crate::vector::Vector3;

/// 行優先の 3×3 行列。
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Matrix3 {
    /// 行 `[row][col]`。
    pub rows: [[f64; 3]; 3],
}

impl Matrix3 {
    /// 単位行列。
    pub const IDENTITY: Matrix3 = Matrix3 {
        rows: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
    };

    /// 行から構築。
    pub const fn new(rows: [[f64; 3]; 3]) -> Self {
        Matrix3 { rows }
    }

    /// x 軸まわりの座標回転 R1(θ)。
    pub fn rotation_x(theta: f64) -> Matrix3 {
        let (s, c) = theta.sin_cos();
        Matrix3::new([[1.0, 0.0, 0.0], [0.0, c, s], [0.0, -s, c]])
    }

    /// y 軸まわりの座標回転 R2(θ)。
    pub fn rotation_y(theta: f64) -> Matrix3 {
        let (s, c) = theta.sin_cos();
        Matrix3::new([[c, 0.0, -s], [0.0, 1.0, 0.0], [s, 0.0, c]])
    }

    /// z 軸まわりの座標回転 R3(θ)。
    pub fn rotation_z(theta: f64) -> Matrix3 {
        let (s, c) = theta.sin_cos();
        Matrix3::new([[c, s, 0.0], [-s, c, 0.0], [0.0, 0.0, 1.0]])
    }

    /// 行列 × ベクトル。
    pub fn mul_vec(&self, v: Vector3) -> Vector3 {
        let r = &self.rows;
        Vector3::new(
            r[0][0] * v.x + r[0][1] * v.y + r[0][2] * v.z,
            r[1][0] * v.x + r[1][1] * v.y + r[1][2] * v.z,
            r[2][0] * v.x + r[2][1] * v.y + r[2][2] * v.z,
        )
    }

    /// 行列 × 行列（self を後段に適用: `self * other`）。
    pub fn mul_mat(&self, other: &Matrix3) -> Matrix3 {
        let mut out = [[0.0; 3]; 3];
        for (i, row) in out.iter_mut().enumerate() {
            for (j, cell) in row.iter_mut().enumerate() {
                *cell = self.rows[i][0] * other.rows[0][j]
                    + self.rows[i][1] * other.rows[1][j]
                    + self.rows[i][2] * other.rows[2][j];
            }
        }
        Matrix3::new(out)
    }

    /// 転置（回転行列では逆回転に等しい）。
    pub fn transpose(&self) -> Matrix3 {
        let r = &self.rows;
        Matrix3::new([
            [r[0][0], r[1][0], r[2][0]],
            [r[0][1], r[1][1], r[2][1]],
            [r[0][2], r[1][2], r[2][2]],
        ])
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use core::f64::consts::PI;

    fn vclose(a: Vector3, b: Vector3) -> bool {
        (a.x - b.x).abs() < 1e-12 && (a.y - b.y).abs() < 1e-12 && (a.z - b.z).abs() < 1e-12
    }

    #[test]
    fn identity_is_neutral() {
        let v = Vector3::new(1.0, 2.0, 3.0);
        assert_eq!(Matrix3::IDENTITY.mul_vec(v), v);
    }

    #[test]
    fn rotation_z_is_coordinate_rotation() {
        // R3(+90°) で x 軸ベクトルの座標は (0, -1, 0) になる（受動回転）。
        let r = Matrix3::rotation_z(PI / 2.0);
        assert!(vclose(
            r.mul_vec(Vector3::new(1.0, 0.0, 0.0)),
            Vector3::new(0.0, -1.0, 0.0)
        ));
    }

    #[test]
    fn rotation_x_is_coordinate_rotation() {
        let r = Matrix3::rotation_x(PI / 2.0);
        assert!(vclose(
            r.mul_vec(Vector3::new(0.0, 1.0, 0.0)),
            Vector3::new(0.0, 0.0, -1.0)
        ));
    }

    #[test]
    fn transpose_is_inverse_rotation() {
        let r = Matrix3::rotation_y(0.7);
        let prod = r.mul_mat(&r.transpose());
        for i in 0..3 {
            for j in 0..3 {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!((prod.rows[i][j] - expected).abs() < 1e-12);
            }
        }
    }

    #[test]
    fn composed_z_rotations_add_angles() {
        let a = Matrix3::rotation_z(0.3);
        let b = Matrix3::rotation_z(0.4);
        let combined = a.mul_mat(&b);
        let direct = Matrix3::rotation_z(0.7);
        for i in 0..3 {
            for j in 0..3 {
                assert!((combined.rows[i][j] - direct.rows[i][j]).abs() < 1e-12);
            }
        }
    }
}
