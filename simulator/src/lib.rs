#![feature(int_log)]
use num_complex::Complex32;
use std::fmt;

pub mod single_qubit_gate;
pub mod two_qubit_gate;
pub mod random_circuit;

pub use single_qubit_gate::*;
pub use two_qubit_gate::*;

#[derive(Debug, Clone, PartialEq)]
pub struct State(Vec<Complex32>);

impl State {
    pub fn apply_single_qubit_gate(&self, index: usize, gate: &SingleQubitGate) -> Self {
        let mut result: Vec<Complex32> = vec![Complex32::new(0.0, 0.0); self.0.len()];

        for i in 0..self.0.len() {
            let i_negated = i ^ (1 << index);

            result[i] = if (i & 1 << index) == 0 {
                gate.0[0][0] * self.0[i] + gate.0[0][1] * self.0[i_negated]
            } else {
                gate.0[1][1] * self.0[i] + gate.0[1][0] * self.0[i_negated]
            }
        }
        Self(result)
    }

    pub fn apply_two_qubit_gate(
        &self,
        first_index: usize,
        second_index: usize,
        gate: &TwoQubitGate,
    ) -> Self {
        let mut result: Vec<Complex32> = vec![Complex32::new(0.0, 0.0); self.0.len()];

        for i in 0..self.0.len() {
            //dbg!(i);
            let gate_index = if second_index == 0 {
                (i & 1 << first_index) >> first_index | (i & 1 << second_index) << 1
            } else {
                (i & 1 << first_index) >> first_index
                    | (i & 1 << second_index) >> (second_index - 1)
            };
            //dbg!(gate_index);

            result[i] = gate.0[gate_index][0]
                * self.0[(i & !(1 << first_index)) & !(1 << second_index)]
                + gate.0[gate_index][1] * self.0[(i | (1 << first_index)) & !(1 << second_index)]
                + gate.0[gate_index][2] * self.0[(i & !(1 << first_index)) | (1 << second_index)]
                + gate.0[gate_index][3] * self.0[(i | (1 << first_index)) | (1 << second_index)];
        }
        Self(result)
    }
}

impl fmt::Display for State {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for (idx, amp) in self.0.iter().enumerate() {
            let len = self.0.len().ilog2() as usize;
            writeln!(
                f,
                "{} ({}{:+}i)",
                format!("|{:0len$b}>", idx),
                amp.re,
                amp.im
            )?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f32::consts::{FRAC_PI_2, FRAC_PI_3, FRAC_PI_4, FRAC_PI_8, PI, SQRT_2};

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

            let result = state.apply_single_qubit_gate(1, &X);

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

            let result = state.apply_single_qubit_gate(0, &X);

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

            let result = state.apply_single_qubit_gate(1, &X);

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

            let result = state.apply_single_qubit_gate(1, &Z);

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
            let state = State(vec![
                Complex32::new(1.0, 0.0),
                Complex32::new(0.0, 0.0),
                Complex32::new(0.0, 0.0),
                Complex32::new(0.0, 0.0),
            ]);

            let result = state.apply_two_qubit_gate(0, 1, &CNOT);

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

            let result = state.apply_two_qubit_gate(0, 1, &CNOT);

            let expected_state = State(vec![
                Complex32::new(0.0, 0.0),
                Complex32::new(1.0, 0.0),
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

            let result = state.apply_two_qubit_gate(0, 1, &CNOT);

            let expected_state = State(vec![
                Complex32::new(0.0, 0.0),
                Complex32::new(0.0, 0.0),
                Complex32::new(0.0, 0.0),
                Complex32::new(1.0, 0.0),
            ]);

            assert_eq!(result, expected_state);

            let state = State(vec![
                Complex32::new(0.0, 0.0),
                Complex32::new(0.0, 0.0),
                Complex32::new(0.0, 0.0),
                Complex32::new(1.0, 0.0),
            ]);

            let result = state.apply_two_qubit_gate(0, 1, &CNOT);

            let expected_state = State(vec![
                Complex32::new(0.0, 0.0),
                Complex32::new(0.0, 0.0),
                Complex32::new(1.0, 0.0),
                Complex32::new(0.0, 0.0),
            ]);

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

            let result = state.apply_two_qubit_gate(0, 1, &suren_gate);

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

            let result = state.apply_two_qubit_gate(0, 1, &suren_gate);

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

            let result = state.apply_two_qubit_gate(0, 1, &suren_gate);

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

            let result = state.apply_two_qubit_gate(0, 1, &suren_gate);

            let expected_state = State(vec![
                Complex32::new(0.0, 0.0),
                Complex32::new(0.0, 0.0),
                Complex32::new(0.0, 0.0),
                Complex32::new(3.0_f32.cos(), 3.0_f32.sin()),
            ]);

            assert_eq!(result, expected_state);
        }

        #[test]
        fn test_suren_gate() {
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
                .apply_single_qubit_gate(3, &r_y)
                .apply_single_qubit_gate(2, &Y)
                .apply_two_qubit_gate(2, 1, &CNOT)
                .apply_single_qubit_gate(2, &r_x)
                .apply_two_qubit_gate(3, 0, &CNOT)
                .apply_single_qubit_gate(3, &T)
                .apply_single_qubit_gate(1, &X)
                .apply_single_qubit_gate(0, &H)
                .apply_single_qubit_gate(3, &H)
                .apply_two_qubit_gate(0, 3, &suren_gate)
                .apply_single_qubit_gate(0, &H)
                .apply_single_qubit_gate(3, &H);

            println!("{}", &result);

            todo!();
        }
    }
}
