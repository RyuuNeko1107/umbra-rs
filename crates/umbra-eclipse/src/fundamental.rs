//! ベッセル基本面の基底（`docs/issues/ISSUE-020`、`docs/algorithms/06-besselian.md`）。
//!
//! 影軸（太陽方向 ẑ）に直交する基本面に、右手系の直交基底を張る:
//! - ẑ = 影軸方向（太陽向き。赤道座標 ICRS/CIRS で与える）
//! - x̂ = (天の北極 × ẑ) を正規化（基本面と赤道面の交線方向＝東）
//! - ŷ = ẑ × x̂（基本面内で天の北側へ射影した方向）
//!
//! 極端な配置（ẑ が天の北極と平行）でも壊れないよう、退化時は `DegenerateGeometry`。

use umbra_core::{UnitVector3, Vector3};

use crate::error::EclipseError;

/// ベッセル基本面の右手系直交基底。
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FundamentalPlaneBasis {
    /// x 軸（東方向、基本面と赤道面の交線）。
    pub x_axis: UnitVector3,
    /// y 軸（基本面内の天の北方向）。
    pub y_axis: UnitVector3,
    /// z 軸（影軸＝太陽方向）。
    pub z_axis: UnitVector3,
}

/// 影軸方向 `z_axis`（太陽向き、赤道座標）から基本面基底を構成する。
pub fn fundamental_plane_basis(z_axis: UnitVector3) -> Result<FundamentalPlaneBasis, EclipseError> {
    let z = z_axis.get();
    let celestial_north = Vector3::new(0.0, 0.0, 1.0);
    // x̂ = 北極 × ẑ（赤道面内＝東）。ẑ が極と平行なら零ベクトル → 退化。
    let x_axis = celestial_north
        .cross(z)
        .normalized()
        .ok_or(EclipseError::DegenerateGeometry)?;
    // ŷ = ẑ × x̂（基本面内の北側）。直交単位ベクトルの外積ゆえ単位。
    let y_axis = z
        .cross(x_axis.get())
        .normalized()
        .ok_or(EclipseError::DegenerateGeometry)?;
    Ok(FundamentalPlaneBasis {
        x_axis,
        y_axis,
        z_axis,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn unit(x: f64, y: f64, z: f64) -> UnitVector3 {
        Vector3::new(x, y, z).normalized().unwrap()
    }

    fn dot(a: Vector3, b: Vector3) -> f64 {
        a.dot(b)
    }

    #[test]
    fn axis_aligned_basis_is_canonical() {
        // ẑ=+x（太陽方向）→ x̂=北極×ẑ=(0,1,0), ŷ=ẑ×x̂=(0,0,1)。
        let b = fundamental_plane_basis(unit(1.0, 0.0, 0.0)).unwrap();
        assert_eq!(b.x_axis.get(), Vector3::new(0.0, 1.0, 0.0));
        assert_eq!(b.y_axis.get(), Vector3::new(0.0, 0.0, 1.0));
        assert_eq!(b.z_axis.get(), Vector3::new(1.0, 0.0, 0.0));
    }

    #[test]
    fn basis_is_orthonormal_and_right_handed() {
        let z = unit(0.3, -0.5, 0.8); // 任意の傾いた軸
        let b = fundamental_plane_basis(z).unwrap();
        let (x, y, zz) = (b.x_axis.get(), b.y_axis.get(), b.z_axis.get());
        // 単位長。
        for v in [x, y, zz] {
            assert!((v.norm() - 1.0).abs() < 1e-12);
        }
        // 直交。
        assert!(dot(x, y).abs() < 1e-12);
        assert!(dot(y, zz).abs() < 1e-12);
        assert!(dot(zz, x).abs() < 1e-12);
        // 右手系: x̂ × ŷ = ẑ。
        let cross = x.cross(y);
        assert!((cross - zz).norm() < 1e-12);
    }

    #[test]
    fn x_axis_lies_in_equatorial_plane() {
        // x̂ = 北極 × ẑ なので必ず赤道面内（z 成分 0）。
        let b = fundamental_plane_basis(unit(0.3, -0.5, 0.8)).unwrap();
        assert!(b.x_axis.get().z.abs() < 1e-12);
    }

    #[test]
    fn y_axis_points_toward_celestial_north() {
        // 軸が赤道近く（低赤緯）なら ŷ は北半球側（z 成分 > 0）。
        let b = fundamental_plane_basis(unit(1.0, 0.0, 0.1)).unwrap();
        assert!(b.y_axis.get().z > 0.0);
    }

    #[test]
    fn degenerate_when_axis_parallel_to_pole() {
        assert_eq!(
            fundamental_plane_basis(unit(0.0, 0.0, 1.0)).unwrap_err(),
            EclipseError::DegenerateGeometry
        );
    }
}
