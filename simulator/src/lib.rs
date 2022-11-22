#![feature(int_log)]
#[cfg(feature = "indicatif")]
use indicatif::{ProgressIterator, ProgressStyle};
use num_complex::Complex32;
use std::f32::consts::FRAC_PI_4;

pub mod block;
mod ops;
pub mod random_circuit;
pub mod single_qubit_gate;
pub mod two_qubit_gate;

pub use block::Block;
pub use single_qubit_gate::SingleQubitGate;
pub use two_qubit_gate::TwoQubitGate;

use single_qubit_gate::gates::{H, T, X, Z};
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

    pub fn grover_diffusion<T>(
        &mut self,
        qubits_to_be_diffused: T,
        ancillas: T,
        temp_state: &mut State,
    ) where
        T: std::ops::Deref<Target = [usize]>,
    {
        assert_eq!(qubits_to_be_diffused.len(), ancillas.len() + 2);
        assert!(qubits_to_be_diffused.len() > 2);

        for &idx in qubits_to_be_diffused.iter() {
            self.apply_single_qubit_gate(idx, &H, temp_state);
            self.apply_single_qubit_gate(idx, &X, temp_state);
        }

        self.apply_n_controlled_gate(
            &qubits_to_be_diffused[1..],
            &ancillas,
            qubits_to_be_diffused[0],
            &Z,
            temp_state,
        );

        for &idx in qubits_to_be_diffused.iter() {
            self.apply_single_qubit_gate(idx, &X, temp_state);
            self.apply_single_qubit_gate(idx, &H, temp_state);
        }
    }

    pub fn apply_grover_algorithm<T, O>(
        &mut self,
        search_qubits: T,
        ancillas: T,
        M: usize,
        oracle: O,
        temp_state: &mut State,
    ) where
        T: std::ops::Deref<Target = [usize]>,
        O: Iterator<Item = Block<T>> + Clone,
    {
        let N = 2_u32.pow(search_qubits.len() as u32);
        let n = (FRAC_PI_4 * (N as f32 / M as f32).sqrt()).floor() as usize;

        for i in search_qubits.deref() {
            self.apply_single_qubit_gate(*i, &H, temp_state);
        }

        for _ in 0..n {
            self.evaluate_circuit(oracle.clone(), temp_state);

            self.grover_diffusion(search_qubits.deref(), ancillas.deref(), temp_state);
        }
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

    pub fn apply_n_controlled_gate<T>(
        &mut self,
        controllers: T,
        ancillas: T,
        target: usize,
        u: &SingleQubitGate,
        temp_state: &mut State,
    ) where
        T: std::ops::Deref<Target = [usize]>,
    {
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
        self.apply_toffoli_gate(controllers[0], controllers[1], ancillas[0], temp_state);
    }

    pub fn norm(&self) -> f32 {
        self.scalar_product(self).norm()
    }

    pub fn distance(&self, rhs: &Self) -> f32 {
        self.0
            .iter()
            .zip(rhs.0.iter())
            .map(|(x, y)| (x - y).norm())
            .reduce(|x, y| x + y)
            .unwrap()
    }

    pub fn scalar_product(&self, rhs: &Self) -> Complex32 {
        self.0
            .iter()
            .zip(rhs.0.iter())
            .map(|(&x, &y)| x.conj() * y)
            .reduce(|x, y| x + y)
            .unwrap()
    }

    pub fn evaluate_circuit<I, T>(&mut self, circuit: I, temp_state: &mut State)
    where
        I: Iterator<Item = Block<T>>,
        T: std::ops::Deref<Target = [usize]>,
    {
        for gate in circuit {
            self.evaluate_gate(&gate, temp_state);
        }
    }

    #[cfg(feature = "indicatif")]
    pub fn evaluate_circuit_progress_bar<I, T>(&mut self, circuit: I, temp_state: &mut State)
    where
        I: Iterator<Item = Block<T>> + std::iter::ExactSizeIterator,
        T: std::ops::Deref<Target = [usize]>,
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

    fn evaluate_gate<'a, T>(&mut self, gate: &'a Block<T>, temp_state: &mut State)
    where
        T: std::ops::Deref<Target = [usize]>,
    {
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
            Block::ToffoliGate {
                control_0_qubit_idx,
                control_1_qubit_idx,
                target_qubit_idx,
            } => self.apply_toffoli_gate(
                *control_0_qubit_idx as usize,
                *control_1_qubit_idx as usize,
                *target_qubit_idx as usize,
                temp_state,
            ),
            Block::CCGate {
                control_0_qubit_idx,
                control_1_qubit_idx,
                target_qubit_idx,
                root_gate,
            } => self.apply_controlled_controlled_gate(
                *control_0_qubit_idx as usize,
                *control_1_qubit_idx as usize,
                *target_qubit_idx as usize,
                root_gate,
                temp_state,
            ),
            Block::NCGate {
                controllers,
                ancillas,
                target,
                gate,
            } => self.apply_n_controlled_gate(
                controllers.deref(),
                ancillas.deref(),
                *target,
                gate,
                temp_state,
            ),
            Block::GroverDiffusion {
                diffusion_qubits,
                ancillas,
            } => self.grover_diffusion(diffusion_qubits.deref(), ancillas.deref(), temp_state),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::single_qubit_gate::gates::*;
    use crate::two_qubit_gate::gates::*;
    use num_complex::ComplexFloat;
    use std::f32::consts::{FRAC_1_SQRT_2, FRAC_PI_3, PI, SQRT_2};

    mod single_qubit_gate {
        use super::*;
        use pretty_assertions::assert_eq;

        #[test]
        fn single_qubit_x() {
            let mut state = State::from_bit_str("0").unwrap();
            let mut temp_state = State::from_bit_str("0").unwrap();

            state.apply_single_qubit_gate(0, &X, &mut temp_state);

            let expected_state = State::from_bit_str("1").unwrap();

            assert_eq!(state, expected_state);
        }

        #[test]
        fn two_qubits_x() {
            let mut state = State::from_bit_str("00").unwrap();
            let mut temp_state = State::from_bit_str("00").unwrap();

            state.apply_single_qubit_gate(0, &X, &mut temp_state);

            let expected_state = State::from_bit_str("01").unwrap();

            assert_eq!(state, expected_state);

            let mut state = State::from_bit_str("10").unwrap();

            state.apply_single_qubit_gate(0, &X, &mut temp_state);

            let expected_state = State::from_bit_str("11").unwrap();

            assert_eq!(state, expected_state);
        }

        #[test]
        fn three_qubits_x() {
            let mut state = State::from_bit_str("000").unwrap();
            let mut temp_state = State::from_bit_str("000").unwrap();

            state.apply_single_qubit_gate(0, &X, &mut temp_state);

            let expected_state = State::from_bit_str("001").unwrap();

            assert_eq!(state, expected_state);

            let mut state = State::from_bit_str("010").unwrap();

            state.apply_single_qubit_gate(0, &X, &mut temp_state);

            let expected_state = State::from_bit_str("011").unwrap();

            assert_eq!(state, expected_state);

            let mut state = State::from_bit_str("000").unwrap();

            state.apply_single_qubit_gate(1, &X, &mut temp_state);

            let expected_state = State::from_bit_str("010").unwrap();

            assert_eq!(state, expected_state);

            let mut state = State::from_bit_str("010").unwrap();

            state.apply_single_qubit_gate(1, &X, &mut temp_state);

            let expected_state = State::from_bit_str("000").unwrap();

            assert_eq!(state, expected_state);
        }

        #[test]
        fn single_qubit_z() {
            let mut qubit = State(vec![Complex32::new(1.0, 0.0), Complex32::new(0.0, 0.0)]);
            let mut tmp_qubit = State(vec![Complex32::new(1.0, 0.0), Complex32::new(0.0, 0.0)]);

            qubit.apply_single_qubit_gate(0, &Z, &mut tmp_qubit);

            let expected_qubit = State(vec![Complex32::new(1.0, 0.0), Complex32::new(0.0, 0.0)]);

            assert_eq!(qubit, expected_qubit);

            let mut qubit = State(vec![Complex32::new(0.0, 0.0), Complex32::new(1.0, 0.0)]);

            qubit.apply_single_qubit_gate(0, &Z, &mut tmp_qubit);

            let expected_qubit = State(vec![Complex32::new(0.0, 0.0), Complex32::new(-1.0, 0.0)]);

            assert_eq!(qubit, expected_qubit);
        }

        #[test]
        fn two_qubits_z() {
            let mut state = State(vec![
                Complex32::new(0.0, 0.0),
                Complex32::new(1.0, 0.0),
                Complex32::new(0.0, 0.0),
                Complex32::new(0.0, 0.0),
            ]);
            let mut tmp_state = State(vec![
                Complex32::new(0.0, 0.0),
                Complex32::new(1.0, 0.0),
                Complex32::new(0.0, 0.0),
                Complex32::new(0.0, 0.0),
            ]);

            state.apply_single_qubit_gate(1, &Z, &mut tmp_state);

            let expected_state = State(vec![
                Complex32::new(0.0, 0.0),
                Complex32::new(1.0, 0.0),
                Complex32::new(0.0, 0.0),
                Complex32::new(0.0, 0.0),
            ]);

            assert_eq!(state, expected_state);

            let mut state = State(vec![
                Complex32::new(0.0, 0.0),
                Complex32::new(1.0, 0.0),
                Complex32::new(0.0, 0.0),
                Complex32::new(0.0, 0.0),
            ]);

            state.apply_single_qubit_gate(0, &Z, &mut tmp_state);

            let expected_state = State(vec![
                Complex32::new(0.0, 0.0),
                Complex32::new(-1.0, 0.0),
                Complex32::new(0.0, 0.0),
                Complex32::new(0.0, 0.0),
            ]);

            assert_eq!(state, expected_state);
        }
    }

    mod two_qubit_gate {
        use super::*;
        use pretty_assertions::assert_eq;

        #[test]
        fn cnot() {
            let mut state = State::from_bit_str("00").unwrap();
            let mut temp_state = State::from_bit_str("00").unwrap();

            state.apply_two_qubit_gate(0, 1, &CNOT, &mut temp_state);

            let expected_state = State::from_bit_str("00").unwrap();

            assert_eq!(state, expected_state);

            let mut state = State::from_bit_str("01").unwrap();

            state.apply_two_qubit_gate(0, 1, &CNOT, &mut temp_state);

            let expected_state = State::from_bit_str("11").unwrap();

            assert_eq!(state, expected_state);

            let mut state = State::from_bit_str("10").unwrap();

            state.apply_two_qubit_gate(0, 1, &CNOT, &mut temp_state);

            let expected_state = State::from_bit_str("10").unwrap();

            assert_eq!(state, expected_state);

            let mut state = State::from_bit_str("11").unwrap();

            state.apply_two_qubit_gate(0, 1, &CNOT, &mut temp_state);

            let expected_state = State::from_bit_str("01").unwrap();

            assert_eq!(state, expected_state);
        }

        #[test]
        fn suren() {
            let mut state = State::from_bit_str("00").unwrap();
            let mut temp_state = State::from_bit_str("00").unwrap();

            let suren_gate = suren_gate(0.0, 1.0, 2.0, 3.0);

            state.apply_two_qubit_gate(1, 0, &suren_gate, &mut temp_state);

            let expected_state = State(vec![
                Complex32::new(1.0, 0.0),
                Complex32::new(0.0, 0.0),
                Complex32::new(0.0, 0.0),
                Complex32::new(0.0, 0.0),
            ]);

            assert_eq!(state, expected_state);

            let mut state = State(vec![
                Complex32::new(0.0, 0.0),
                Complex32::new(1.0, 0.0),
                Complex32::new(0.0, 0.0),
                Complex32::new(0.0, 0.0),
            ]);

            state.apply_two_qubit_gate(1, 0, &suren_gate, &mut temp_state);

            let expected_state = State(vec![
                Complex32::new(0.0, 0.0),
                Complex32::new(1.0_f32.cos(), 1.0_f32.sin()),
                Complex32::new(0.0, 0.0),
                Complex32::new(0.0, 0.0),
            ]);

            assert_eq!(state, expected_state);

            let mut state = State(vec![
                Complex32::new(0.0, 0.0),
                Complex32::new(0.0, 0.0),
                Complex32::new(1.0, 0.0),
                Complex32::new(0.0, 0.0),
            ]);

            state.apply_two_qubit_gate(1, 0, &suren_gate, &mut temp_state);

            let expected_state = State(vec![
                Complex32::new(0.0, 0.0),
                Complex32::new(0.0, 0.0),
                Complex32::new(2.0_f32.cos(), 2.0_f32.sin()),
                Complex32::new(0.0, 0.0),
            ]);

            assert_eq!(state, expected_state);

            let mut state = State(vec![
                Complex32::new(0.0, 0.0),
                Complex32::new(0.0, 0.0),
                Complex32::new(0.0, 0.0),
                Complex32::new(1.0, 0.0),
            ]);

            state.apply_two_qubit_gate(1, 0, &suren_gate, &mut temp_state);

            let expected_state = State(vec![
                Complex32::new(0.0, 0.0),
                Complex32::new(0.0, 0.0),
                Complex32::new(0.0, 0.0),
                Complex32::new(3.0_f32.cos(), 3.0_f32.sin()),
            ]);

            assert_eq!(state, expected_state);
        }

        #[test]
        #[ignore]
        fn test_suren_circuit() {
            let mut state = State::from_bit_str("0000").unwrap();
            let mut temp_state = State::from_bit_str("0000").unwrap();

            let r_y = r_y(FRAC_PI_3);
            let r_x = r_x(PI * SQRT_2);
            let suren_gate = suren_gate(0.0, 1.0, 2.0, 3.0);

            state.apply_two_qubit_gate(0, 1, &CNOT, &mut temp_state);
            state.apply_single_qubit_gate(2, &Y, &mut temp_state);
            state.apply_two_qubit_gate(1, 2, &CNOT, &mut temp_state);
            state.apply_single_qubit_gate(2, &r_x, &mut temp_state);
            state.apply_single_qubit_gate(3, &r_y, &mut temp_state);
            state.apply_two_qubit_gate(0, 3, &CNOT, &mut temp_state);
            state.apply_single_qubit_gate(1, &X, &mut temp_state);
            state.apply_single_qubit_gate(3, &T, &mut temp_state);
            state.apply_single_qubit_gate(0, &H, &mut temp_state);
            state.apply_single_qubit_gate(3, &H, &mut temp_state);
            state.apply_two_qubit_gate(3, 0, &suren_gate, &mut temp_state);
            state.apply_single_qubit_gate(3, &H, &mut temp_state);
            state.apply_single_qubit_gate(0, &H, &mut temp_state);

            println!("{}", &state);
        }
    }

    #[test]
    fn toffoli() {
        let mut temp_state = State::from_bit_str("000").unwrap();
        for i in 0..7 {
            let control_0 = i & 1;
            let control_1 = (i & 2) >> 1;
            let target = (i & 4) >> 2;
            let mut state =
                State::from_bit_str(&format!("{}{}{}", target, control_1, control_0)).unwrap();
            let control_0_bool = control_0 == 1;
            let control_1_bool = control_1 == 1;
            let target_bool = target == 1;

            let expected_result = u32::from(target_bool ^ (control_0_bool & control_1_bool));

            let expected_state =
                State::from_bit_str(&format!("{}{}{}", expected_result, control_1, control_0))
                    .unwrap();

            state.apply_toffoli_gate(0, 1, 2, &mut temp_state);
            assert!(state.distance(&expected_state) < 1e-4);
        }
    }

    #[test]
    fn ccnot() {
        let mut temp_state = State::from_bit_str("000").unwrap();
        for i in 0..7 {
            let control_0 = i & 1;
            let control_1 = (i & 2) >> 1;
            let target = (i & 4) >> 2;
            let mut state =
                State::from_bit_str(&format!("{}{}{}", target, control_1, control_0)).unwrap();
            let control_0_bool = control_0 == 1;
            let control_1_bool = control_1 == 1;
            let target_bool = target == 1;

            let expected_result = u32::from(target_bool ^ (control_0_bool & control_1_bool));

            let expected_state =
                State::from_bit_str(&format!("{}{}{}", expected_result, control_1, control_0))
                    .unwrap();

            state.apply_controlled_controlled_gate(0, 1, 2, &SX, &mut temp_state);
            assert!(state.distance(&expected_state) < 1e-4);
        }
    }

    #[test]
    fn cccnot() {
        let mut temp_state = State::from_bit_str("000000").unwrap();
        for i in 0..15 {
            let control_0 = i & 1;
            let control_1 = (i & 2) >> 1;
            let control_2 = (i & 4) >> 2;
            let target = (i & 8) >> 3;
            let mut state = State::from_bit_str(&dbg!(format!(
                "00{}{}{}{}",
                target, control_2, control_1, control_0
            )))
            .unwrap();
            let control_0_bool = control_0 == 1;
            let control_1_bool = control_1 == 1;
            let control_2_bool = control_2 == 1;
            let target_bool = target == 1;

            let expected_result =
                u32::from(target_bool ^ (control_0_bool & control_1_bool & control_2_bool));

            let expected_state = State::from_bit_str(&dbg!(format!(
                "00{}{}{}{}",
                expected_result, control_2, control_1, control_0
            )))
            .unwrap();

            state.apply_n_controlled_gate(vec![0, 1, 2], vec![4, 5], 3, &X, &mut temp_state);
            assert!(dbg!(state).distance(&dbg!(expected_state)) < 1e-4);
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
    fn grover_naive() {
        let expected_state = State::from_bit_str("001111").unwrap();
        let mut state = State::from_bit_str("000000").unwrap();
        let mut temp_state = State::from_bit_str("000000").unwrap();
        //                                  00_0000
        //                                  2 ancillas
        //                                  4 qubits

        let n = (FRAC_PI_4 * (2_u32.pow(4) as f32).sqrt()).floor() as usize;

        for i in 0..4 {
            state.apply_single_qubit_gate(i, &H, &mut temp_state);
        }

        for i in 0..n {
            state.apply_n_controlled_gate(vec![0, 1, 2], vec![4, 5], 3, &Z, &mut temp_state);

            state.grover_diffusion(vec![0, 1, 2, 3], vec![4, 5], &mut temp_state);
            println!("{} {}", i, state.scalar_product(&expected_state).abs());
        }

        let theta = (2.0 * 15.0.sqrt() / 16.0).asin();
        let probabity = ((n as f32 + 0.5) * theta).sin().powi(2);

        assert!(
            (dbg!(state.scalar_product(&expected_state).abs().powi(2)) - dbg!(probabity)).abs()
                <= 1e-3
        );
    }

    #[test]
    fn grover() {
        let expected_state = State::from_bit_str("001111").unwrap();
        let mut state = State::from_bit_str("000000").unwrap();
        let mut temp_state = State::from_bit_str("000000").unwrap();
        //                                  00_0000
        //                                  2 ancillas
        //                                  4 qubits

        let oracle = vec![Block::NCGate {
            controllers: vec![0, 1, 2],
            ancillas: vec![4, 5],
            target: 3,
            gate: Z,
        }];

        state.apply_grover_algorithm(
            vec![0, 1, 2, 3],
            vec![4, 5],
            1,
            oracle.into_iter(),
            &mut temp_state,
        );

        let N = 2.0_f32.powi(4);
        let M = 1.0;
        let k = (FRAC_PI_4 * (N / M).sqrt()).floor();
        let probabity = expected_probability(N, M, k);

        assert!(
            (dbg!(state.scalar_product(&expected_state).abs().powi(2)) - dbg!(probabity)).abs()
                <= 1e-3
        );
    }

    #[test]
    #[ignore]
    fn grover_2_answers() {
        let expected_state_0 = State::from_bit_str("00011110").unwrap();
        let expected_state_1 = State::from_bit_str("00011111").unwrap();
        let expected_state = (expected_state_1 + expected_state_0) * FRAC_1_SQRT_2;

        let mut state = State::from_bit_str("00000000").unwrap();
        let mut temp_state = State::from_bit_str("00000000").unwrap();
        //                                  000_00000
        //                                  3 ancillas
        //                                  5 qubits

        let oracle = vec![Block::NCGate {
            controllers: vec![1, 2, 3],
            ancillas: vec![5, 6],
            target: 4,
            gate: Z,
        }];

        state.apply_grover_algorithm(
            vec![0, 1, 2, 3, 4],
            vec![5, 6, 7],
            2,
            oracle.into_iter(),
            &mut temp_state,
        );

        //dbg!(&state);
        //dbg!(&expected_state);

        let N = 2.0_f32.powi(5);
        let M = 2.0;
        let k = (FRAC_PI_4 * (N / M).sqrt()).floor();
        let probabity = expected_probability(N, M, k);

        for i in 0..k as usize + 5 {
            print!("{} ", expected_probability(N, M, i as f32));
        }
        println!();

        assert!(
            (dbg!(state.scalar_product(&expected_state).abs().powi(2)) - dbg!(probabity)).abs()
                <= 1e-3
        );
    }

    fn expected_probability(n_variants: f32, m_answers: f32, k_iterations: f32) -> f32 {
        let theta = (2.0 * (m_answers * (n_variants - m_answers).sqrt()) / n_variants).asin();
        ((2.0 * k_iterations + 1.0) / 2.0 * theta).sin().powi(2)
    }
}
