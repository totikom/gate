use num_complex::Complex32;
use std::f32::consts::TAU;

mod block;
mod rand;

use super::single_qubit_gate::consts::arbitrary_unitary_matrix;
use super::{SingleQubitGate, TwoQubitGate};
use block::Block;
use rand::{Rand, K};

pub struct RandomCircuitIter {
    n_qubits: u64,
    rand: Rand,
}

impl RandomCircuitIter {
    pub fn new(seed: u64, n_qubits: u64) -> Self {
        let rand = Rand::new(seed);
        Self { rand, n_qubits }
    }
}

impl Iterator for RandomCircuitIter {
    type Item = Block;

    fn next(&mut self) -> Option<Self::Item> {
        let parity = self.rand.next().unwrap();
        if parity % 2 == 1 {
            let gate = random_2x2_unitary(&mut self.rand);
            let qubit_idx = self.rand.next().unwrap() % self.n_qubits;
            Some(Block::SingleQubitGate { gate, qubit_idx })
        } else {
            let gate = random_4x4_unitary(&mut self.rand);

            let control_qubit_idx = self.rand.next().unwrap() % self.n_qubits;
            let mut target_qubit_idx = self.rand.next().unwrap() % self.n_qubits;

            while target_qubit_idx == control_qubit_idx {
                target_qubit_idx = self.rand.next().unwrap() % self.n_qubits;
            }

            Some(Block::TwoQubitGate {
                gate,
                control_qubit_idx,
                target_qubit_idx,
            })
        }
    }
}

fn random_2x2_unitary(rand: &mut Rand) -> SingleQubitGate {
    let a = rand.next().unwrap() as f32 / (K as f32) * TAU;
    let b = rand.next().unwrap() as f32 / (K as f32) * TAU;
    let c = rand.next().unwrap() as f32 / (K as f32) * TAU;
    let d = rand.next().unwrap() as f32 / (K as f32) * TAU;

    arbitrary_unitary_matrix(a, b, c, d)
}

fn random_4x4_unitary(rand: &mut Rand) -> TwoQubitGate {
    let u1 = random_2x2_unitary(rand);
    let u2 = random_2x2_unitary(rand);
    let u3 = random_2x2_unitary(rand);
    let u4 = random_2x2_unitary(rand);
    let u5 = random_2x2_unitary(rand);
    let u6 = random_2x2_unitary(rand);

    construct_type_1_matrix(&u1, &u2)
        * construct_type_2_matrix(&u3)
        * construct_type_1_matrix(&u4, &u5)
        * construct_type_2_matrix(&u6)
}

fn construct_type_1_matrix(u1: &SingleQubitGate, u2: &SingleQubitGate) -> TwoQubitGate {
    let mut result = [[Complex32::new(0.0, 0.0); 4]; 4];

    for (row, u1_row) in result.iter_mut().zip(u1.0.iter()) {
        row[0..2].copy_from_slice(u1_row);
    }
    for (row, u2_row) in result.iter_mut().skip(2).zip(u2.0.iter()) {
        row[2..4].copy_from_slice(u2_row);
    }

    TwoQubitGate(result)
}

fn construct_type_2_matrix(u: &SingleQubitGate) -> TwoQubitGate {
    let mut result = [[Complex32::new(0.0, 0.0); 4]; 4];

    for (row, u_row) in result.iter_mut().skip(1).zip(u.0.iter()) {
        row[1..3].copy_from_slice(u_row);
    }
    result[0][0] = Complex32::new(1.0, 0.0);
    result[3][3] = Complex32::new(1.0, 0.0);

    TwoQubitGate(result)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::single_qubit_gate::consts::*;
    use pretty_assertions::assert_eq;

    #[test]
    fn type_1_test() {
        let expected_gate = TwoQubitGate([
            [
                Complex32::new(0.0, 0.0),
                Complex32::new(1.0, 0.0),
                Complex32::new(0.0, 0.0),
                Complex32::new(0.0, 0.0),
            ],
            [
                Complex32::new(1.0, 0.0),
                Complex32::new(0.0, 0.0),
                Complex32::new(0.0, 0.0),
                Complex32::new(0.0, 0.0),
            ],
            [
                Complex32::new(0.0, 0.0),
                Complex32::new(0.0, 0.0),
                Complex32::new(0.0, 0.0),
                Complex32::new(1.0, 0.0),
            ],
            [
                Complex32::new(0.0, 0.0),
                Complex32::new(0.0, 0.0),
                Complex32::new(1.0, 0.0),
                Complex32::new(0.0, 0.0),
            ],
        ]);

        assert_eq!(construct_type_1_matrix(&X, &X), expected_gate);
    }

    #[test]
    fn type_2_test() {
        let expected_gate = TwoQubitGate([
            [
                Complex32::new(1.0, 0.0),
                Complex32::new(0.0, 0.0),
                Complex32::new(0.0, 0.0),
                Complex32::new(0.0, 0.0),
            ],
            [
                Complex32::new(0.0, 0.0),
                Complex32::new(0.0, 0.0),
                Complex32::new(1.0, 0.0),
                Complex32::new(0.0, 0.0),
            ],
            [
                Complex32::new(0.0, 0.0),
                Complex32::new(1.0, 0.0),
                Complex32::new(0.0, 0.0),
                Complex32::new(0.0, 0.0),
            ],
            [
                Complex32::new(0.0, 0.0),
                Complex32::new(0.0, 0.0),
                Complex32::new(0.0, 0.0),
                Complex32::new(1.0, 0.0),
            ],
        ]);

        assert_eq!(construct_type_2_matrix(&X), expected_gate);
    }

    #[test]
    fn unitary_2x2() {
        let mut rand = Rand::new(0);
        for _ in 0..10 {
            assert!((random_2x2_unitary(&mut rand).determinant().norm() -1.0).abs() < 1e-4);
        }
    }

    #[test]
    fn unitary_4x4() {
        let mut rand = Rand::new(0);
        for _ in 0..10 {
            let matrix = random_4x4_unitary(&mut rand);
            assert!((matrix.determinant().norm() -1.0).abs() < 1e-4);
        }
    }
}
