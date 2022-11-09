#![feature(int_log)]
#[cfg(feature = "indicatif")]
use indicatif::{ProgressIterator, ProgressStyle};
use num_complex::Complex32;

use std::fmt;

pub mod block;
pub mod random_circuit;
pub mod single_qubit_gate;
pub mod two_qubit_gate;

pub use block::Block;
pub use single_qubit_gate::SingleQubitGate;
pub use two_qubit_gate::TwoQubitGate;

use single_qubit_gate::gates::{H, T};
use two_qubit_gate::gates::{controlled_u, CNOT};

#[derive(Clone, PartialEq)]
pub struct State(Vec<Complex32>);

impl State {
    pub fn new(state: Vec<Complex32>) -> Self {
        Self(state)
    }

    pub fn from_bit_str(state_str: &str) -> Result<Self, std::num::ParseIntError> {
        let mut state = vec![Complex32::new(0.0, 0.0); 2_usize.pow(state_str.len() as u32)];
        let index = usize::from_str_radix(state_str, 2)?;

        state[index] = Complex32::new(1.0, 0.0);
        Ok(Self(state))
    }

    pub fn apply_single_qubit_gate(
        &mut self,
        index: usize,
        gate: &SingleQubitGate,
        temp_state: &mut State,
    ) {
        let result = &mut temp_state.0;
        let chunk_len = if result.len() <= num_cpus::get() {
            result.len()
        } else {
            result.len() / num_cpus::get()
        };

        crossbeam::scope(|scope| {
            for (thread_idx, chunk) in result.chunks_mut(chunk_len).enumerate() {
                let old_state = &self.0;
                scope.spawn(move |_| {
                    let offset = thread_idx * chunk_len;
                    for (i, val) in chunk.iter_mut().enumerate() {
                        let i = i + offset;
                        let i_negated = i ^ (1 << index);

                        *val = if (i & 1 << index) == 0 {
                            gate.0[0][0] * old_state[i] + gate.0[0][1] * old_state[i_negated]
                        } else {
                            gate.0[1][1] * old_state[i] + gate.0[1][0] * old_state[i_negated]
                        }
                    }
                });
            }
        })
        .unwrap();

        std::mem::swap(&mut self.0, &mut temp_state.0);
    }

    pub fn apply_two_qubit_gate(
        &mut self,
        control_index: usize,
        target_index: usize,
        gate: &TwoQubitGate,
        temp_state: &mut State,
    ) {
        let result = &mut temp_state.0;

        let chunk_len = if result.len() <= num_cpus::get() {
            result.len()
        } else {
            result.len() / num_cpus::get()
        };

        crossbeam::scope(|scope| {
            for (thread_idx, chunk) in result.chunks_mut(chunk_len).enumerate() {
                let old_state = &self.0;
                scope.spawn(move |_| {
                    let offset = thread_idx * chunk_len;

                    for (i, val) in chunk.iter_mut().enumerate() {
                        //dbg!(i);
                        let i = i + offset;
                        let gate_index = if control_index == 0 {
                            (i & 1 << target_index) >> target_index | (i & 1 << control_index) << 1
                        } else {
                            (i & 1 << target_index) >> target_index
                                | (i & 1 << control_index) >> (control_index - 1)
                        };
                        //dbg!(gate_index);

                        *val = gate.0[gate_index][0]
                            * old_state[(i & !(1 << target_index)) & !(1 << control_index)]
                            + gate.0[gate_index][1]
                                * old_state[(i | (1 << target_index)) & !(1 << control_index)]
                            + gate.0[gate_index][2]
                                * old_state[(i & !(1 << target_index)) | (1 << control_index)]
                            + gate.0[gate_index][3]
                                * old_state[(i | (1 << target_index)) | (1 << control_index)];
                    }
                });
            }
        })
        .unwrap();

        std::mem::swap(&mut self.0, &mut temp_state.0);
    }

    pub fn diffuse(&mut self, index: usize) {
        self.0[2_usize.pow(index as u32)] = -self.0[2_usize.pow(index as u32)];
    }

    pub fn apply_toffoli_gate(
        &mut self,
        control_1: usize,
        control_2: usize,
        target: usize,
        temp_state: &mut State,
    ) {
        self.apply_single_qubit_gate(target, &H, temp_state);
        self.apply_two_qubit_gate(control_2, target, &CNOT, temp_state);
        self.apply_single_qubit_gate(target, &T.h_conj(), temp_state);
        self.apply_two_qubit_gate(control_1, target, &CNOT, temp_state);
        self.apply_single_qubit_gate(target, &T, temp_state);
        self.apply_two_qubit_gate(control_2, target, &CNOT, temp_state);
        self.apply_single_qubit_gate(target, &T.h_conj(), temp_state);
        self.apply_two_qubit_gate(control_1, target, &CNOT, temp_state);
        self.apply_single_qubit_gate(control_2, &T, temp_state);
        self.apply_single_qubit_gate(target, &T, temp_state);
        self.apply_two_qubit_gate(control_1, control_2, &CNOT, temp_state);
        self.apply_single_qubit_gate(target, &H, temp_state);
        self.apply_single_qubit_gate(control_1, &T, temp_state);
        self.apply_single_qubit_gate(control_2, &T.h_conj(), temp_state);
        self.apply_two_qubit_gate(control_1, control_2, &CNOT, temp_state);
    }

    pub fn apply_controlled_controlled_gate(
        &mut self,
        control_1: usize,
        control_2: usize,
        target: usize,
        sqr_root_u: &SingleQubitGate,
        temp_state: &mut State,
    ) {
        let controlled_root_u = controlled_u(sqr_root_u);
        self.apply_two_qubit_gate(control_2, target, &controlled_root_u, temp_state);
        self.apply_two_qubit_gate(control_1, control_2, &CNOT, temp_state);
        self.apply_two_qubit_gate(control_2, target, &controlled_root_u.h_conj(), temp_state);
        self.apply_two_qubit_gate(control_1, control_2, &CNOT, temp_state);
        self.apply_two_qubit_gate(control_1, target, &controlled_root_u, temp_state);
    }

    pub fn apply_n_controlled_gate(
        mut self,
        controllers: Vec<usize>,
        ancillas: Vec<usize>,
        target: usize,
        u: &SingleQubitGate,
        temp_state: &mut State,
    ) {
        assert_eq!(controllers.len(), ancillas.len() + 1);
        assert!(controllers.len() > 2);

        let controlled_u = controlled_u(u);

        self.apply_toffoli_gate(controllers[0], controllers[1], ancillas[0], temp_state);

        for (ancilla_pair, control) in ancillas.windows(2).zip(controllers.iter().skip(2)) {
            self.apply_toffoli_gate(*control, ancilla_pair[0], ancilla_pair[1], temp_state);
        }

        self.apply_two_qubit_gate(
            ancillas[ancillas.len() - 1],
            target,
            &controlled_u,
            temp_state,
        );

        for (ancilla_pair, control) in ancillas.windows(2).rev().zip(controllers.iter().rev()) {
            self.apply_toffoli_gate(*control, ancilla_pair[0], ancilla_pair[1], temp_state);
        }
    }

    pub fn norm(&self) -> f32 {
        self.0
            .iter()
            .map(|x| x.norm())
            .reduce(|x, y| x + y)
            .unwrap()
    }

    pub fn distance(&self, rhs: &Self) -> f32 {
        self.0
            .iter()
            .zip(rhs.0.iter())
            .map(|(x, y)| (x - y).norm())
            .reduce(|x, y| x + y)
            .unwrap()
    }

    pub fn evaluate_circuit<I>(mut self, circuit: I, temp_state: &mut State)
    where
        I: Iterator<Item = Block>,
    {
        for gate in circuit {
            self.evaluate_gate(&gate, temp_state);
        }
    }

    #[cfg(feature = "indicatif")]
    pub fn evaluate_circuit_progress_bar<I>(&mut self, circuit: I, temp_state: &mut State)
    where
        I: Iterator<Item = Block> + std::iter::ExactSizeIterator,
    {
        for gate in circuit.progress().with_style(
            ProgressStyle::with_template(
                "[{elapsed_precise}] {wide_bar:.cyan/blue} {pos:>7}/{len:7} {msg}",
            )
            .unwrap()
            .progress_chars("##-"),
        ) {
            self.evaluate_gate(&gate, temp_state);
        }
    }

    fn evaluate_gate<'a>(&mut self, gate: &'a Block, temp_state: &mut State) {
        match gate {
            Block::SingleQubitGate { gate, qubit_idx } => {
                self.apply_single_qubit_gate(*qubit_idx as usize, gate, temp_state)
            }
            Block::TwoQubitGate {
                gate,
                control_qubit_idx,
                target_qubit_idx,
            } => self.apply_two_qubit_gate(
                *control_qubit_idx as usize,
                *target_qubit_idx as usize,
                gate,
                temp_state,
            ),
        }
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::single_qubit_gate::gates::*;
    use crate::two_qubit_gate::gates::*;
    use std::f32::consts::{FRAC_PI_3, PI, SQRT_2};

    mod single_qubit_gate {
        use super::*;
        use pretty_assertions::assert_eq;

        #[test]
        fn single_qubit_x() {
            let qubit = State(vec![Complex32::new(1.0, 0.0), Complex32::new(0.0, 0.0)]);

            let result = qubit.apply_single_qubit_gate(0, &X);

            let expected_qubit = State(vec![Complex32::new(0.0, 0.0), Complex32::new(1.0, 0.0)]);

            assert_eq!(result, expected_qubit);
        }

        #[test]
        fn two_qubits_x() {
            let state = State(vec![
                Complex32::new(1.0, 0.0),
                Complex32::new(0.0, 0.0),
                Complex32::new(0.0, 0.0),
                Complex32::new(0.0, 0.0),
            ]);

            let result = state.clone().apply_single_qubit_gate(1, &X);

            let expected_state = State(vec![
                Complex32::new(0.0, 0.0),
                Complex32::new(0.0, 0.0),
                Complex32::new(1.0, 0.0),
                Complex32::new(0.0, 0.0),
            ]);

            assert_eq!(result, expected_state);

            let result = state.apply_single_qubit_gate(0, &X);

            let expected_state = State(vec![
                Complex32::new(0.0, 0.0),
                Complex32::new(1.0, 0.0),
                Complex32::new(0.0, 0.0),
                Complex32::new(0.0, 0.0),
            ]);

            assert_eq!(result, expected_state);
        }

        #[test]
        fn three_qubits_x() {
            let state = State(vec![
                Complex32::new(1.0, 0.0),
                Complex32::new(0.0, 0.0),
                Complex32::new(0.0, 0.0),
                Complex32::new(0.0, 0.0),
                Complex32::new(0.0, 0.0),
                Complex32::new(0.0, 0.0),
                Complex32::new(0.0, 0.0),
                Complex32::new(0.0, 0.0),
            ]);

            let result = state.clone().apply_single_qubit_gate(0, &X);

            let expected_state = State(vec![
                Complex32::new(0.0, 0.0),
                Complex32::new(1.0, 0.0),
                Complex32::new(0.0, 0.0),
                Complex32::new(0.0, 0.0),
                Complex32::new(0.0, 0.0),
                Complex32::new(0.0, 0.0),
                Complex32::new(0.0, 0.0),
                Complex32::new(0.0, 0.0),
            ]);

            assert_eq!(result, expected_state);

            let result = state.clone().apply_single_qubit_gate(1, &X);

            let expected_state = State(vec![
                Complex32::new(0.0, 0.0),
                Complex32::new(0.0, 0.0),
                Complex32::new(1.0, 0.0),
                Complex32::new(0.0, 0.0),
                Complex32::new(0.0, 0.0),
                Complex32::new(0.0, 0.0),
                Complex32::new(0.0, 0.0),
                Complex32::new(0.0, 0.0),
            ]);

            assert_eq!(result, expected_state);

            let result = state.apply_single_qubit_gate(2, &X);

            let expected_state = State(vec![
                Complex32::new(0.0, 0.0),
                Complex32::new(0.0, 0.0),
                Complex32::new(0.0, 0.0),
                Complex32::new(0.0, 0.0),
                Complex32::new(1.0, 0.0),
                Complex32::new(0.0, 0.0),
                Complex32::new(0.0, 0.0),
                Complex32::new(0.0, 0.0),
            ]);

            assert_eq!(result, expected_state);
        }
        #[test]
        fn single_qubit_z() {
            let qubit = State(vec![Complex32::new(1.0, 0.0), Complex32::new(0.0, 0.0)]);

            let result = qubit.apply_single_qubit_gate(0, &Z);

            let expected_qubit = State(vec![Complex32::new(1.0, 0.0), Complex32::new(0.0, 0.0)]);

            assert_eq!(result, expected_qubit);

            let qubit = State(vec![Complex32::new(0.0, 0.0), Complex32::new(1.0, 0.0)]);

            let result = qubit.apply_single_qubit_gate(0, &Z);

            let expected_qubit = State(vec![Complex32::new(0.0, 0.0), Complex32::new(-1.0, 0.0)]);

            assert_eq!(result, expected_qubit);
        }

        #[test]
        fn two_qubits_z() {
            let state = State(vec![
                Complex32::new(0.0, 0.0),
                Complex32::new(1.0, 0.0),
                Complex32::new(0.0, 0.0),
                Complex32::new(0.0, 0.0),
            ]);

            let result = state.clone().apply_single_qubit_gate(1, &Z);

            let expected_state = State(vec![
                Complex32::new(0.0, 0.0),
                Complex32::new(1.0, 0.0),
                Complex32::new(0.0, 0.0),
                Complex32::new(0.0, 0.0),
            ]);

            assert_eq!(result, expected_state);

            let result = state.apply_single_qubit_gate(0, &Z);

            let expected_state = State(vec![
                Complex32::new(0.0, 0.0),
                Complex32::new(-1.0, 0.0),
                Complex32::new(0.0, 0.0),
                Complex32::new(0.0, 0.0),
            ]);

            assert_eq!(result, expected_state);
        }
    }

    mod two_qubit_gate {
        use super::*;
        use pretty_assertions::assert_eq;

        #[test]
        fn cnot() {
            let state = State::from_bit_str("00").unwrap();

            let result = state.apply_two_qubit_gate(0, 1, &CNOT);

            let expected_state = State::from_bit_str("00").unwrap();

            assert_eq!(result, expected_state);

            let state = State::from_bit_str("01").unwrap();

            let result = state.apply_two_qubit_gate(0, 1, &CNOT);

            let expected_state = State::from_bit_str("11").unwrap();

            assert_eq!(result, expected_state);

            let state = State::from_bit_str("10").unwrap();

            let result = state.apply_two_qubit_gate(0, 1, &CNOT);

            let expected_state = State::from_bit_str("10").unwrap();

            assert_eq!(result, expected_state);

            let state = State::from_bit_str("11").unwrap();

            let result = state.apply_two_qubit_gate(0, 1, &CNOT);

            let expected_state = State::from_bit_str("01").unwrap();

            assert_eq!(result, expected_state);
        }

        #[test]
        fn suren() {
            let state = State(vec![
                Complex32::new(1.0, 0.0),
                Complex32::new(0.0, 0.0),
                Complex32::new(0.0, 0.0),
                Complex32::new(0.0, 0.0),
            ]);

            let suren_gate = suren_gate(0.0, 1.0, 2.0, 3.0);

            let result = state.apply_two_qubit_gate(1, 0, &suren_gate);

            let expected_state = State(vec![
                Complex32::new(1.0, 0.0),
                Complex32::new(0.0, 0.0),
                Complex32::new(0.0, 0.0),
                Complex32::new(0.0, 0.0),
            ]);

            assert_eq!(result, expected_state);

            let state = State(vec![
                Complex32::new(0.0, 0.0),
                Complex32::new(1.0, 0.0),
                Complex32::new(0.0, 0.0),
                Complex32::new(0.0, 0.0),
            ]);

            let result = state.apply_two_qubit_gate(1, 0, &suren_gate);

            let expected_state = State(vec![
                Complex32::new(0.0, 0.0),
                Complex32::new(1.0_f32.cos(), 1.0_f32.sin()),
                Complex32::new(0.0, 0.0),
                Complex32::new(0.0, 0.0),
            ]);

            assert_eq!(result, expected_state);

            let state = State(vec![
                Complex32::new(0.0, 0.0),
                Complex32::new(0.0, 0.0),
                Complex32::new(1.0, 0.0),
                Complex32::new(0.0, 0.0),
            ]);

            let result = state.apply_two_qubit_gate(1, 0, &suren_gate);

            let expected_state = State(vec![
                Complex32::new(0.0, 0.0),
                Complex32::new(0.0, 0.0),
                Complex32::new(2.0_f32.cos(), 2.0_f32.sin()),
                Complex32::new(0.0, 0.0),
            ]);

            assert_eq!(result, expected_state);

            let state = State(vec![
                Complex32::new(0.0, 0.0),
                Complex32::new(0.0, 0.0),
                Complex32::new(0.0, 0.0),
                Complex32::new(1.0, 0.0),
            ]);

            let result = state.apply_two_qubit_gate(1, 0, &suren_gate);

            let expected_state = State(vec![
                Complex32::new(0.0, 0.0),
                Complex32::new(0.0, 0.0),
                Complex32::new(0.0, 0.0),
                Complex32::new(3.0_f32.cos(), 3.0_f32.sin()),
            ]);

            assert_eq!(result, expected_state);
        }

        #[test]
        #[ignore]
        fn test_suren_circuit() {
            let state = State(vec![
                Complex32::new(1.0, 0.0),
                Complex32::new(0.0, 0.0),
                Complex32::new(0.0, 0.0),
                Complex32::new(0.0, 0.0),
                Complex32::new(0.0, 0.0),
                Complex32::new(0.0, 0.0),
                Complex32::new(0.0, 0.0),
                Complex32::new(0.0, 0.0),
                Complex32::new(0.0, 0.0),
                Complex32::new(0.0, 0.0),
                Complex32::new(0.0, 0.0),
                Complex32::new(0.0, 0.0),
                Complex32::new(0.0, 0.0),
                Complex32::new(0.0, 0.0),
                Complex32::new(0.0, 0.0),
                Complex32::new(0.0, 0.0),
            ]);

            let r_y = r_y(FRAC_PI_3);
            let r_x = r_x(PI * SQRT_2);
            let suren_gate = suren_gate(0.0, 1.0, 2.0, 3.0);

            let result = state
                .apply_two_qubit_gate(0, 1, &CNOT)
                .apply_single_qubit_gate(2, &Y)
                .apply_two_qubit_gate(1, 2, &CNOT)
                .apply_single_qubit_gate(2, &r_x)
                .apply_single_qubit_gate(3, &r_y)
                .apply_two_qubit_gate(0, 3, &CNOT)
                .apply_single_qubit_gate(1, &X)
                .apply_single_qubit_gate(3, &T)
                .apply_single_qubit_gate(0, &H)
                .apply_single_qubit_gate(3, &H)
                .apply_two_qubit_gate(3, 0, &suren_gate)
                .apply_single_qubit_gate(3, &H)
                .apply_single_qubit_gate(0, &H);

            println!("{}", &result);
        }
    }

    mod toffoli {
        use super::*;

        #[test]
        fn state_000() {
            let state = State::from_bit_str("000").unwrap();

            let result = state.apply_toffoli_gate(0, 1, 2);

            let expected_state = State::from_bit_str("000").unwrap();

            assert!(result.distance(&expected_state) < 1e-4);
        }

        #[test]
        fn state_001() {
            let state = State::from_bit_str("001").unwrap();

            let result = state.apply_toffoli_gate(0, 1, 2);

            let expected_state = State::from_bit_str("001").unwrap();

            assert!(dbg!(result).distance(&dbg!(expected_state)) < 1e-4);
        }

        #[test]
        fn state_010() {
            let state = State::from_bit_str("010").unwrap();

            let result = state.apply_toffoli_gate(0, 1, 2);

            let expected_state = State::from_bit_str("010").unwrap();

            assert!(dbg!(result).distance(&dbg!(expected_state)) < 1e-4);
        }

        #[test]
        fn state_011() {
            let state = State::from_bit_str("011").unwrap();

            let result = state.apply_toffoli_gate(0, 1, 2);

            let expected_state = State::from_bit_str("111").unwrap();

            assert!(dbg!(result).distance(&dbg!(expected_state)) < 1e-4);
        }

        #[test]
        fn state_100() {
            let state = State::from_bit_str("100").unwrap();

            let result = state.apply_toffoli_gate(0, 1, 2);

            let expected_state = State::from_bit_str("100").unwrap();

            assert!(dbg!(result).distance(&dbg!(expected_state)) < 1e-4);
        }

        #[test]
        fn state_101() {
            let state = State::from_bit_str("101").unwrap();

            let result = state.apply_toffoli_gate(0, 1, 2);

            let expected_state = State::from_bit_str("101").unwrap();

            assert!(result.distance(&expected_state) < 1e-4);
        }

        #[test]
        fn state_110() {
            let state = State::from_bit_str("110").unwrap();

            let result = state.apply_toffoli_gate(0, 1, 2);

            let expected_state = State::from_bit_str("110").unwrap();

            assert!(result.distance(&expected_state) < 1e-4);
        }

        #[test]
        fn state_111() {
            let state = State::from_bit_str("111").unwrap();

            let result = state.apply_toffoli_gate(0, 1, 2);

            let expected_state = State::from_bit_str("011").unwrap();

            assert!(result.distance(&expected_state) < 1e-4);
        }
    }

    mod ccnot {
        use super::*;

        #[test]
        fn state_000() {
            let state = State::from_bit_str("000").unwrap();

            let result = state.apply_controlled_controlled_gate(0, 1, 2, &SX);

            let expected_state = State::from_bit_str("000").unwrap();

            assert!(result.distance(&expected_state) < 1e-4);
        }

        #[test]
        fn state_001() {
            let state = State::from_bit_str("001").unwrap();

            let result = state.apply_controlled_controlled_gate(0, 1, 2, &SX);

            let expected_state = State::from_bit_str("001").unwrap();

            assert!(dbg!(result).distance(&dbg!(expected_state)) < 1e-4);
        }

        #[test]
        fn state_010() {
            let state = State::from_bit_str("010").unwrap();

            let result = state.apply_controlled_controlled_gate(0, 1, 2, &SX);

            let expected_state = State::from_bit_str("010").unwrap();

            assert!(dbg!(result).distance(&dbg!(expected_state)) < 1e-4);
        }

        #[test]
        fn state_011() {
            let state = State::from_bit_str("011").unwrap();

            let result = state.apply_controlled_controlled_gate(0, 1, 2, &SX);

            let expected_state = State::from_bit_str("111").unwrap();

            assert!(dbg!(result).distance(&dbg!(expected_state)) < 1e-4);
        }

        #[test]
        fn state_100() {
            let state = State::from_bit_str("100").unwrap();

            let result = state.apply_controlled_controlled_gate(0, 1, 2, &SX);

            let expected_state = State::from_bit_str("100").unwrap();

            assert!(dbg!(result).distance(&dbg!(expected_state)) < 1e-4);
        }

        #[test]
        fn state_101() {
            let state = State::from_bit_str("101").unwrap();

            let result = state.apply_controlled_controlled_gate(0, 1, 2, &SX);

            let expected_state = State::from_bit_str("101").unwrap();

            assert!(result.distance(&expected_state) < 1e-4);
        }

        #[test]
        fn state_110() {
            let state = State::from_bit_str("110").unwrap();

            let result = state.apply_controlled_controlled_gate(0, 1, 2, &SX);

            let expected_state = State::from_bit_str("110").unwrap();

            assert!(result.distance(&expected_state) < 1e-4);
        }

        #[test]
        fn state_111() {
            let state = State::from_bit_str("111").unwrap();

            let result = state.apply_controlled_controlled_gate(0, 1, 2, &SX);

            let expected_state = State::from_bit_str("011").unwrap();

            assert!(result.distance(&expected_state) < 1e-4);
        }
    }

    mod str {
        use super::*;
        use pretty_assertions::assert_eq;

        #[test]
        fn state_000() {
            let result = State::from_bit_str("000").unwrap();

            let expected_state = State(vec![
                Complex32::new(1.0, 0.0),
                Complex32::new(0.0, 0.0),
                Complex32::new(0.0, 0.0),
                Complex32::new(0.0, 0.0),
                Complex32::new(0.0, 0.0),
                Complex32::new(0.0, 0.0),
                Complex32::new(0.0, 0.0),
                Complex32::new(0.0, 0.0),
            ]);

            assert_eq!(result, expected_state);
        }

        #[test]
        fn state_001() {
            let result = State::from_bit_str("001").unwrap();

            let expected_state = State(vec![
                Complex32::new(0.0, 0.0),
                Complex32::new(1.0, 0.0),
                Complex32::new(0.0, 0.0),
                Complex32::new(0.0, 0.0),
                Complex32::new(0.0, 0.0),
                Complex32::new(0.0, 0.0),
                Complex32::new(0.0, 0.0),
                Complex32::new(0.0, 0.0),
            ]);

            assert_eq!(result, expected_state);
        }

        #[test]
        fn state_10() {
            let result = State::from_bit_str("10").unwrap();

            let expected_state = State(vec![
                Complex32::new(0.0, 0.0),
                Complex32::new(0.0, 0.0),
                Complex32::new(1.0, 0.0),
                Complex32::new(0.0, 0.0),
            ]);

            assert_eq!(result, expected_state);
        }
    }

    #[test]
    fn grover() {
        let mut state = State::from_bit_str("0000000").unwrap();
        //let mut state = State::from_br("00_0_0000").unwrap();

        for i in 0..5 {
            state = state.apply_single_qubit_gate(i, &H);
        }

        for _ in 0..2_usize.pow(2) {
            state = state.apply_n_controlled_gate(vec![0, 1, 2], vec![5, 6], 3, &X);

            for i in 0..5 {
                state = state.apply_single_qubit_gate(i, &H);
            }
            state = state.diffuse(2);
            for i in 0..5 {
                state = state.apply_single_qubit_gate(i, &H);
            }
        }

        println!("{}", state);
        todo!();
    }
}
