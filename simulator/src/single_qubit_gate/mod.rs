use num_complex::Complex32;

pub mod consts;
mod ops;

#[derive(Clone, Copy, PartialEq)]
pub struct SingleQubitGate(pub [[Complex32; 2]; 2]);

impl SingleQubitGate {
    pub fn distanse(&self, rhs: &Self) -> f32 {
        let mut sum = 0.0;

        for i in 0..2 {
            for j in 0..2 {
                sum += (self.0[i][j] - rhs.0[i][j]).norm_sqr();
            }
        }
        sum.sqrt()
    }

    pub fn determinant(&self) -> Complex32 {
        self.0[0][0] * self.0[1][1] - self.0[0][1] * self.0[1][0]
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use consts::*;
    use pretty_assertions::assert_eq;

    #[test]
    fn involution() {
        assert_eq!(X * X, I);
        assert_eq!(Y * Y, I);
        assert_eq!(Z * Z, I);
        dbg!(I - H * H);
        assert!(I.distanse(&(H * H)) < 1e-4);
    }
}
