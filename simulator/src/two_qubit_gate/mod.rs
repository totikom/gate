use num_complex::Complex32;

pub mod gates;
mod ops;

#[derive(Clone, Copy, PartialEq)]
pub struct TwoQubitGate(pub [[Complex32; 4]; 4]);
impl TwoQubitGate {
    pub fn distanse(&self, rhs: &Self) -> f32 {
        let mut sum = 0.0;

        for i in 0..4 {
            for j in 0..4 {
                sum += (self.0[i][j] - rhs.0[i][j]).norm_sqr();
            }
        }
        sum.sqrt()
    }

    pub fn determinant(&self) -> Complex32 {
        self.0[0][0]
            * (self.0[1][1] * self.0[2][2] * self.0[3][3]
                + self.0[1][2] * self.0[2][3] * self.0[3][1]
                + self.0[1][3] * self.0[2][1] * self.0[3][2]
                - self.0[1][3] * self.0[2][2] * self.0[3][1]
                - self.0[1][2] * self.0[2][1] * self.0[3][3]
                - self.0[1][1] * self.0[2][3] * self.0[3][2])
            - self.0[1][0]
                * (self.0[0][1] * self.0[2][2] * self.0[3][3]
                    + self.0[0][2] * self.0[2][3] * self.0[3][1]
                    + self.0[0][3] * self.0[2][1] * self.0[3][2]
                    - self.0[0][3] * self.0[2][2] * self.0[3][1]
                    - self.0[0][2] * self.0[2][1] * self.0[3][3]
                    - self.0[0][1] * self.0[2][3] * self.0[3][2])
            + self.0[2][0]
                * (self.0[0][1] * self.0[1][2] * self.0[3][3]
                    + self.0[0][2] * self.0[1][3] * self.0[3][1]
                    + self.0[0][3] * self.0[1][1] * self.0[3][2]
                    - self.0[0][3] * self.0[1][2] * self.0[3][1]
                    - self.0[0][2] * self.0[1][1] * self.0[3][3]
                    - self.0[0][1] * self.0[1][3] * self.0[3][2])
            - self.0[3][0]
                * (self.0[0][1] * self.0[1][2] * self.0[2][3]
                    + self.0[0][2] * self.0[1][3] * self.0[2][1]
                    + self.0[0][3] * self.0[1][1] * self.0[2][2]
                    - self.0[0][3] * self.0[1][2] * self.0[2][1]
                    - self.0[0][2] * self.0[1][1] * self.0[2][3]
                    - self.0[0][1] * self.0[1][3] * self.0[2][2])
    }

    pub fn h_conj(&self) -> Self {
        let mut result = [[Complex32::new(0.0,0.0);4];4];
        for (i, row) in result.iter_mut().enumerate() {
            for (j, val) in row.iter_mut().enumerate() {
                *val = self.0[j][i].conj();
            }
        }
        Self(result)
    }
}
