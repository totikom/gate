#[cfg(feature = "indicatif")]
use indicatif::{ProgressIterator, ProgressStyle};
use num_complex::Complex32;
use std::f32::consts::{FRAC_PI_4, TAU};

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
        let state_str = state_str.replace(&[' ', '_'], "");
        let mut state = vec![Complex32::new(0.0, 0.0); 2_usize.pow(state_str.len() as u32)];
        let index = usize::from_str_radix(&state_str, 2)?;

        state[index] = Complex32::new(1.0, 0.0);
        Ok(Self(state))
    }

    pub fn apply_single_qubit_gate(
        &mut self,
        index: usize,
        gate: &SingleQubitGate,
        temp_state: &mut State,
    ) {
        debug_assert_eq!(self.0.len(), temp_state.0.len());
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
        debug_assert_eq!(self.0.len(), temp_state.0.len());
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

    pub fn apply_phase_estimation_algorithm<U>(
        &mut self,
        estimation_start_idx: usize,
        estimation_end_idx: usize,
        additional_ancilla_idx: usize,
        additional_ancilla_idx2: usize,
        operator: U,
        temp_state: &mut State,
    ) where
        U: Iterator<Item = Block<Vec<usize>>> + Clone + std::iter::ExactSizeIterator,
    {
        for (logical_idx, i) in (estimation_start_idx..=estimation_end_idx)
            .into_iter()
            .enumerate()
        {
            self.apply_single_qubit_gate(i, &H, temp_state);
            let controlled_u = Self::convert_u(
                operator.clone(),
                i,
                additional_ancilla_idx,
                additional_ancilla_idx2,
            );
            for _ in 0..2_usize.pow(logical_idx as u32) {
                self.evaluate_circuit(controlled_u.clone(), temp_state);
            }
        }

        self.apply_inv_qft(estimation_start_idx, estimation_end_idx, temp_state);
    }

    pub fn apply_qft(&mut self, start_idx: usize, end_idx: usize, temp_state: &mut State) {
        assert!(start_idx <= end_idx);
        for abs_idx in start_idx..end_idx {
            self.apply_single_qubit_gate(abs_idx, &H, temp_state);
            for j in 1..=(end_idx - abs_idx) {
                let angle = TAU / 2.0_f32.powi(j as i32 + 1);
                let r = SingleQubitGate([
                    [Complex32::new(1.0, 0.0), Complex32::new(0.0, 0.0)],
                    [
                        Complex32::new(0.0, 0.0),
                        Complex32::new(angle.cos(), angle.sin()),
                    ],
                ]);
                let r = controlled_u(&r);
                self.apply_two_qubit_gate(abs_idx, abs_idx + j, &r, temp_state);
            }
        }
        self.apply_single_qubit_gate(end_idx, &H, temp_state);

        for i in start_idx..((end_idx - start_idx + 1) / 2) {
            self.swap_qubits(start_idx + i, end_idx - i, temp_state);
        }
    }

    pub fn apply_inv_qft(&mut self, start_idx: usize, end_idx: usize, temp_state: &mut State) {
        assert!(start_idx <= end_idx);
        for i in start_idx..((end_idx - start_idx + 1) / 2) {
            self.swap_qubits(start_idx + i, end_idx - i, temp_state);
        }

        self.apply_single_qubit_gate(end_idx, &H, temp_state);

        for abs_idx in (start_idx..end_idx).into_iter().rev() {
            for j in (1..=(end_idx - abs_idx)).into_iter().rev() {
                let angle = -TAU / 2.0_f32.powi(j as i32 + 1);
                let r = SingleQubitGate([
                    [Complex32::new(1.0, 0.0), Complex32::new(0.0, 0.0)],
                    [
                        Complex32::new(0.0, 0.0),
                        Complex32::new(angle.cos(), angle.sin()),
                    ],
                ]);
                let r = controlled_u(&r);
                self.apply_two_qubit_gate(abs_idx, abs_idx + j, &r, temp_state);
            }
            self.apply_single_qubit_gate(abs_idx, &H, temp_state);
        }
    }

    pub fn swap_qubits(&mut self, first: usize, second: usize, temp_state: &mut State) {
        self.apply_two_qubit_gate(first, second, &CNOT, temp_state);
        self.apply_two_qubit_gate(second, first, &CNOT, temp_state);
        self.apply_two_qubit_gate(first, second, &CNOT, temp_state);
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

    fn convert_u<I>(
        u: I,
        control_idx: usize,
        additional_ancilla_idx: usize,
        additional_ancilla_idx2: usize,
    ) -> impl Iterator<Item = Block<Vec<usize>>> + std::iter::ExactSizeIterator + Clone
    where
        I: Iterator<Item = Block<Vec<usize>>> + std::iter::ExactSizeIterator + Clone,
    {
        u.map(move |b| match b {
            Block::SingleQubitGate { gate, qubit_idx } => {
                let controlled_u = controlled_u(&gate);
                Block::TwoQubitGate {
                    gate: controlled_u,
                    control_qubit_idx: control_idx,
                    target_qubit_idx: qubit_idx,
                }
            }
            Block::TwoQubitGate { .. } => {
                unimplemented!();
            }
            Block::ToffoliGate {
                control_0_qubit_idx,
                control_1_qubit_idx,
                target_qubit_idx,
            } => Block::NCGate {
                controllers: vec![control_0_qubit_idx, control_1_qubit_idx, control_idx],
                ancillas: vec![additional_ancilla_idx, additional_ancilla_idx2],
                target: target_qubit_idx,
                gate: X,
            },
            Block::CCGate {
                control_0_qubit_idx,
                control_1_qubit_idx,
                target_qubit_idx,
                root_gate,
            } => Block::NCGate {
                controllers: vec![control_0_qubit_idx, control_1_qubit_idx, control_idx],
                ancillas: vec![additional_ancilla_idx],
                target: target_qubit_idx,
                gate: root_gate * root_gate,
            },
            Block::NCGate {
                mut controllers,
                mut ancillas,
                target,
                gate,
            } => {
                controllers.push(control_idx);
                ancillas.push(additional_ancilla_idx);
                Block::NCGate {
                    controllers,
                    ancillas,
                    target: target,
                    gate: gate,
                }
            }
            Block::GroverDiffusion { .. } => {
                unimplemented!();
            }
        })
    }

    pub fn kron(&self, state: &State) -> Self {
        let mut new_state = vec![Complex32::new(0.0, 0.0); self.0.len() * state.0.len()];
        for (self_idx, a) in self.0.iter().enumerate() {
            for (state_idx, b) in state.0.iter().enumerate() {
                new_state[self_idx * state.0.len() + state_idx] = a * b;
            }
        }
        Self(new_state)
    }

    pub fn prob_reduce(&self, n_qubits: usize) -> Vec<f32> {
        let mut probs = vec![0.0; 2_usize.pow(n_qubits as u32)];
        for (idx, val) in self.0.iter().enumerate() {
            probs[idx % 2_usize.pow(n_qubits as u32)] += val.norm() * val.norm();
        }
        probs
    }
}

#[cfg(test)]
mod tests;
