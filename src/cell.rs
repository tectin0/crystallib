use cgmath::{vec3, Matrix3, Matrix4};

#[cfg_attr(feature = "serde", derive(serde::Deserialize, serde::Serialize))]
#[cfg_attr(feature = "serde", serde(default))]
#[derive(Default, Debug, Clone, PartialEq)]
pub struct Cell {
    pub a: f64,
    pub b: f64,
    pub c: f64,
    pub alpha: f64,
    pub beta: f64,
    pub gamma: f64,
    pub volume: f64,
    pub space_group: String,
    pub space_group_number: u8,
}

impl Cell {
    /// Returns the transformation matrix for the unit cell with Y-up coordinates.
    pub fn transformation_matrix_y_up(&self) -> Matrix4<f32> {
        let x = self.a as f32;
        let y = self.c as f32;
        let z = self.b as f32;
        let angle_yz = (self.alpha as f32).to_radians();
        let angle_xz = (self.gamma as f32).to_radians();
        let angle_yx = (self.beta as f32).to_radians();

        let x0 = x;
        let x1 = 0.0;
        let x2 = 0.0;

        let y0 = y * (angle_yx.cos());
        let y1 = y * (angle_yx.sin());
        let y2 = 0.0;

        let z0 = z * (angle_xz.cos());
        let z1 = z * (angle_yz.cos() - angle_xz.cos() * angle_yx.cos()) / angle_yx.sin();
        let z2 = z
            * (angle_xz.sin().powi(2)
                - (angle_yz.cos() - angle_xz.cos() * angle_yx.cos()).powi(2)
                    / angle_yx.sin().powi(2))
            .sqrt();

        let c0 = vec3(x0, x1, x2);
        let c1 = vec3(y0, y1, y2);
        let c2 = vec3(z0, z1, z2);

        let m = Matrix3::from_cols(c0, c1, c2);

        Matrix4::from(m)
    }
}
