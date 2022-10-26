use std::f32::consts::TAU;

mod rand;

use super::single_qubit_gate::consts::arbitrary_unitary_matrix;
use super::SingleQubitGate;
use rand::{Rand, K};

pub struct RandomCircuitIter {
    n_qubits: usize,
    rand: Rand,
}

impl RandomCircuitIter {
    pub fn new(seed: u64, n_qubits: usize) -> Self {
        let rand = Rand::new(seed);
        Self { rand, n_qubits }
    }
}

fn random_2x2_unitary(rand: &mut Rand) -> SingleQubitGate {
    let a = rand.next().unwrap() as f32 / (K as f32) * TAU;
    let b = rand.next().unwrap() as f32 / (K as f32) * TAU;
    let c = rand.next().unwrap() as f32 / (K as f32) * TAU;
    let d = rand.next().unwrap() as f32 / (K as f32) * TAU;

    arbitrary_unitary_matrix(a, b, c, d)
}
