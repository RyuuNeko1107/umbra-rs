//! 3次元ベクトル（右手系。`docs/conventions.md` §5）。
//!
//! 公開は最小限に保つ（`docs/api-draft.md` §1.4）。フレームの区別は上位層の型で行う。

/// 3次元ベクトル。
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Vector3 {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl Vector3 {
    /// 零ベクトル。
    pub const ZERO: Vector3 = Vector3 {
        x: 0.0,
        y: 0.0,
        z: 0.0,
    };

    /// 構築。
    pub const fn new(x: f64, y: f64, z: f64) -> Self {
        Vector3 { x, y, z }
    }

    /// 内積。
    pub fn dot(self, other: Vector3) -> f64 {
        self.x * other.x + self.y * other.y + self.z * other.z
    }

    /// 外積（右手系）。
    pub fn cross(self, other: Vector3) -> Vector3 {
        Vector3 {
            x: self.y * other.z - self.z * other.y,
            y: self.z * other.x - self.x * other.z,
            z: self.x * other.y - self.y * other.x,
        }
    }

    /// ノルム（ユークリッド長）。
    pub fn norm(self) -> f64 {
        self.dot(self).sqrt()
    }

    /// スカラ倍。
    pub fn scale(self, s: f64) -> Vector3 {
        Vector3 {
            x: self.x * s,
            y: self.y * s,
            z: self.z * s,
        }
    }
}

impl core::ops::Add for Vector3 {
    type Output = Vector3;
    fn add(self, o: Vector3) -> Vector3 {
        Vector3::new(self.x + o.x, self.y + o.y, self.z + o.z)
    }
}

impl core::ops::Sub for Vector3 {
    type Output = Vector3;
    fn sub(self, o: Vector3) -> Vector3 {
        Vector3::new(self.x - o.x, self.y - o.y, self.z - o.z)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dot_and_norm() {
        let a = Vector3::new(3.0, 4.0, 0.0);
        assert_eq!(a.norm(), 5.0);
        assert_eq!(a.dot(Vector3::new(1.0, 0.0, 0.0)), 3.0);
    }

    #[test]
    fn cross_is_right_handed() {
        let x = Vector3::new(1.0, 0.0, 0.0);
        let y = Vector3::new(0.0, 1.0, 0.0);
        // x × y = z
        assert_eq!(x.cross(y), Vector3::new(0.0, 0.0, 1.0));
        // y × x = -z（反交換性）
        assert_eq!(y.cross(x), Vector3::new(0.0, 0.0, -1.0));
    }

    #[test]
    fn add_sub_scale() {
        let a = Vector3::new(1.0, 2.0, 3.0);
        let b = Vector3::new(4.0, 5.0, 6.0);
        assert_eq!(a + b, Vector3::new(5.0, 7.0, 9.0));
        assert_eq!(b - a, Vector3::new(3.0, 3.0, 3.0));
        assert_eq!(a.scale(2.0), Vector3::new(2.0, 4.0, 6.0));
    }
}
