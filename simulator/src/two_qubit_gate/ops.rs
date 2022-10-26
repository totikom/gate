use num_complex::Complex32;

use super::TwoQubitGate;

impl std::ops::Mul for TwoQubitGate {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self {
        let mut result = [[Complex32::new(0.0, 0.0); 4]; 4];

        for (i, row) in result.iter_mut().enumerate() {
            for (j, elem) in row.iter_mut().enumerate() {
                *elem = self.0[i]
                    .iter()
                    .enumerate()
                    .map(|(k, val)| *val * rhs.0[k][j])
                    .reduce(|x, y| x + y)
                    .unwrap();
            }
        }
        Self(result)
    }
}

impl std::ops::Add for TwoQubitGate {
    type Output = Self;

    fn add(self, rhs: Self) -> Self {
        let mut result = [[Complex32::new(0.0, 0.0); 4]; 4];

        for i in 0..4 {
            for j in 0..4 {
                result[i][j] = self.0[i][j] + rhs.0[i][j]
            }
        }
        Self(result)
    }
}

impl std::ops::Sub for TwoQubitGate {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self {
        let mut result = [[Complex32::new(0.0, 0.0); 4]; 4];

        for i in 0..4 {
            for j in 0..4 {
                result[i][j] = self.0[i][j] - rhs.0[i][j]
            }
        }

        Self(result)
    }
}

impl std::ops::Mul<Complex32> for TwoQubitGate {
    type Output = Self;

    fn mul(self, rhs: Complex32) -> Self {
        let mut result = [[Complex32::new(0.0, 0.0); 4]; 4];

        for i in 0..4 {
            for j in 0..4 {
                result[i][j] = self.0[i][j] * rhs
            }
        }

        Self(result)
    }
}
