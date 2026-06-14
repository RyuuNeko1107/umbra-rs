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

impl Vector3 {
    /// 正規化して [`UnitVector3`] にする。零ベクトルは `None`。
    pub fn normalized(self) -> Option<UnitVector3> {
        UnitVector3::from_vector(self)
    }
}

/// 単位ベクトル（生成時に正規化済みであることを型で保証）。
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct UnitVector3(Vector3);

impl UnitVector3 {
    /// ベクトルを正規化して構築。零ベクトルは `None`。
    pub fn from_vector(v: Vector3) -> Option<Self> {
        let n = v.norm();
        if n == 0.0 || !n.is_finite() {
            None
        } else {
            Some(UnitVector3(v.scale(1.0 / n)))
        }
    }

    /// 内側の [`Vector3`]（長さ 1）。
    pub fn get(self) -> Vector3 {
        self.0
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
    fn dot_uses_every_component() {
        // 全成分が非ゼロかつ非対称: 各項の符号・積の取り違えを区別できる。
        let a = Vector3::new(1.0, 2.0, 3.0);
        let b = Vector3::new(4.0, 5.0, 6.0);
        assert_eq!(a.dot(b), 4.0 + 10.0 + 18.0); // 32
    }

    #[test]
    fn norm_uses_every_component() {
        // (1,2,2) の各成分が二乗和に寄与（z=0 にしない）。
        assert_eq!(Vector3::new(1.0, 2.0, 2.0).norm(), 3.0);
    }

    #[test]
    fn cross_all_nonzero_and_orthogonal() {
        let a = Vector3::new(1.0, 2.0, 3.0);
        let b = Vector3::new(4.0, 5.0, 6.0);
        // (2*6-3*5, 3*4-1*6, 1*5-2*4) = (-3, 6, -3)
        assert_eq!(a.cross(b), Vector3::new(-3.0, 6.0, -3.0));
        assert_eq!(b.cross(a), Vector3::new(3.0, -6.0, 3.0)); // 反交換
        assert_eq!(a.cross(b).dot(a), 0.0); // 直交
        assert_eq!(a.cross(b).dot(b), 0.0);
    }

    #[test]
    fn cross_is_right_handed() {
        let x = Vector3::new(1.0, 0.0, 0.0);
        let y = Vector3::new(0.0, 1.0, 0.0);
        assert_eq!(x.cross(y), Vector3::new(0.0, 0.0, 1.0)); // x × y = z
    }

    #[test]
    fn add_sub() {
        let a = Vector3::new(1.0, 2.0, 3.0);
        let b = Vector3::new(4.0, 5.0, 6.0);
        assert_eq!(a + b, Vector3::new(5.0, 7.0, 9.0));
        assert_eq!(b - a, Vector3::new(3.0, 3.0, 3.0));
    }

    #[test]
    fn scale_factor_differs_from_every_component() {
        // 係数 10 はどの成分とも一致しないので * と + を区別できる。
        let a = Vector3::new(1.0, 2.0, 3.0);
        assert_eq!(a.scale(10.0), Vector3::new(10.0, 20.0, 30.0));
    }

    #[test]
    fn normalized_has_unit_length_and_direction() {
        let u = Vector3::new(3.0, 4.0, 0.0).normalized().unwrap();
        assert!((u.get().norm() - 1.0).abs() < 1e-12);
        // 方向保存: (3,4,0)/5
        assert!((u.get().x - 0.6).abs() < 1e-12 && (u.get().y - 0.8).abs() < 1e-12);
    }

    #[test]
    fn normalized_zero_is_none() {
        assert!(Vector3::ZERO.normalized().is_none());
    }
}
