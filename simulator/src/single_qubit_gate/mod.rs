use num_complex::Complex32;

pub mod gates;
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

    pub fn h_conj(&self) -> Self {
        let mut result = [[Complex32::new(0.0, 0.0); 2]; 2];
        for (i, row) in result.iter_mut().enumerate() {
            for (j, val) in row.iter_mut().enumerate() {
                *val = self.0[j][i].conj();
            }
        }
        Self(result)
    }

    pub fn pow(self, pow: usize) -> SingleQubitGate {
        let mut value = SingleQubitGate([
            [Complex32::new(1.0, 0.0), Complex32::new(0.0, 0.0)],
            [Complex32::new(0.0, 0.0), Complex32::new(1.0, 0.0)],
        ]);

        for _ in 0..pow {
            value = value * self;
        }

        value
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use gates::*;
    use pretty_assertions::assert_eq;

    #[test]
    fn involution() {
        assert_eq!(X * X, I);
        assert_eq!(Y * Y, I);
        assert_eq!(Z * Z, I);
        dbg!(I - H * H);
        assert!(I.distanse(&(H * H)) < 1e-4);
        assert_eq!(X * X.h_conj(), I);
        assert_eq!(Y * Y.h_conj(), I);
        assert_eq!(Z * Z.h_conj(), I);
    }
}
