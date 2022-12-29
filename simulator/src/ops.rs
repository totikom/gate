use super::State;
use std::{fmt, ops};
use num_complex::Complex32;

impl ops::Add for State {
    type Output = Self;

    fn add(self, rhs: Self) -> Self {
        let result: Vec<_> = self
            .0
            .into_iter()
            .zip(rhs.0.iter())
            .map(|(x, y)| x + y)
            .collect();

        Self(result)
    }
}

impl ops::Mul<f32> for State {
    type Output = Self;

    fn mul(self, rhs: f32) -> Self {
        let result: Vec<_> = self.0.into_iter().map(|x| x * rhs).collect();

        Self(result)
    }
}

impl ops::Mul<Complex32> for State {
    type Output = Self;

    fn mul(self, rhs: Complex32) -> Self {
        let result: Vec<_> = self.0.into_iter().map(|x| x * rhs).collect();

        Self(result)
    }
}

impl fmt::Display for State {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for (idx, amp) in self.0.iter().enumerate() {
            let len = self.0.len().ilog2() as usize;
            writeln!(f, "|{:0len$b}> ({}{:+}i)", idx, amp.re, amp.im)?;
        }
        Ok(())
    }
}

impl fmt::Debug for State {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self)
    }
}
