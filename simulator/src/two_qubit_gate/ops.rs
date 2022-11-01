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

        for (i, row) in result.iter_mut().enumerate() {
            for (j, val) in row.iter_mut().enumerate() {
                *val = self.0[i][j] + rhs.0[i][j]
            }
        }

        Self(result)
    }
}

impl std::ops::Sub for TwoQubitGate {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self {
        let mut result = [[Complex32::new(0.0, 0.0); 4]; 4];

        for (i, row) in result.iter_mut().enumerate() {
            for (j, val) in row.iter_mut().enumerate() {
                *val = self.0[i][j] - rhs.0[i][j]
            }
        }

        Self(result)
    }
}

impl std::ops::Mul<Complex32> for TwoQubitGate {
    type Output = Self;

    fn mul(self, rhs: Complex32) -> Self {
        let mut result = [[Complex32::new(0.0, 0.0); 4]; 4];

        for (i, row) in result.iter_mut().enumerate() {
            for (j, val) in row.iter_mut().enumerate() {
                *val = self.0[i][j] * rhs
            }
        }

        Self(result)
    }
}

impl std::fmt::Display for TwoQubitGate {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for row in self.0.iter() {
            let mut value_iter = row.iter();
            let value = value_iter.next().unwrap();
            write!(f, "{:.4}{:+.4}i", value.re, value.im)?;
            for value in value_iter {
                write!(f, "\t{:.4}{:+.4}i", value.re, value.im)?;
            }
            writeln!(f)?;
        }
        Ok(())
    }
}

impl std::fmt::Debug for TwoQubitGate {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "\n{}", self)
    }
}
