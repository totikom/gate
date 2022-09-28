use num_complex::{Complex, Complex64};

type Amp = Complex<f32>;

#[derive(Debug, Clone, PartialEq)]
struct State(Vec<Amp>);

#[derive(Debug, Clone, PartialEq)]
struct SingleQubitGate([[Amp; 2]; 2]);

impl State {
    pub fn apply_single_qubit_gate(&self, index: usize, gate: &SingleQubitGate) -> Self {
        let mut result: Vec<Amp> = vec![Complex::new(0.0, 0.0); self.0.len()];

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
}

#[cfg(test)]
mod tests {
    use super::*;
    use pretty_assertions::assert_eq;

    #[test]
    fn single_qubit() {
        let qubit = State(vec![Complex::new(1.0, 0.0), Complex::new(0.0, 0.0)]);

        let matrix = SingleQubitGate([
            [Complex::new(0.0, 0.0), Complex::new(1.0, 0.0)],
            [Complex::new(1.0, 0.0), Complex::new(0.0, 0.0)],
        ]);

        let result = qubit.apply_single_qubit_gate(0, &matrix);

        let expected_qubit = State(vec![Complex::new(0.0, 0.0), Complex::new(1.0, 0.0)]);

        assert_eq!(result, expected_qubit);
    }

    #[test]
    fn two_qubits() {
        let state = State(vec![
            Complex::new(1.0, 0.0),
            Complex::new(0.0, 0.0),
            Complex::new(0.0, 0.0),
            Complex::new(0.0, 0.0),
        ]);

        let matrix = SingleQubitGate([
            [Complex::new(0.0, 0.0), Complex::new(1.0, 0.0)],
            [Complex::new(1.0, 0.0), Complex::new(0.0, 0.0)],
        ]);

        let result = state.apply_single_qubit_gate(1, &matrix);

        let expected_state = State(vec![
            Complex::new(0.0, 0.0),
            Complex::new(0.0, 0.0),
            Complex::new(1.0, 0.0),
            Complex::new(0.0, 0.0),
        ]);

        assert_eq!(result, expected_state);

        let result = state.apply_single_qubit_gate(0, &matrix);

        let expected_state = State(vec![
            Complex::new(0.0, 0.0),
            Complex::new(1.0, 0.0),
            Complex::new(0.0, 0.0),
            Complex::new(0.0, 0.0),
        ]);

        assert_eq!(result, expected_state);
    }

    #[test]
    fn three_qubits() {
        let state = State(vec![
            Complex::new(1.0, 0.0),
            Complex::new(0.0, 0.0),
            Complex::new(0.0, 0.0),
            Complex::new(0.0, 0.0),
            Complex::new(0.0, 0.0),
            Complex::new(0.0, 0.0),
            Complex::new(0.0, 0.0),
            Complex::new(0.0, 0.0),
        ]);

        let matrix = SingleQubitGate([
            [Complex::new(0.0, 0.0), Complex::new(1.0, 0.0)],
            [Complex::new(1.0, 0.0), Complex::new(0.0, 0.0)],
        ]);

        let result = state.apply_single_qubit_gate(0, &matrix);

        let expected_state = State(vec![
            Complex::new(0.0, 0.0),
            Complex::new(1.0, 0.0),
            Complex::new(0.0, 0.0),
            Complex::new(0.0, 0.0),
            Complex::new(0.0, 0.0),
            Complex::new(0.0, 0.0),
            Complex::new(0.0, 0.0),
            Complex::new(0.0, 0.0),
        ]);

        assert_eq!(result, expected_state);

        let result = state.apply_single_qubit_gate(1, &matrix);

        let expected_state = State(vec![
            Complex::new(0.0, 0.0),
            Complex::new(0.0, 0.0),
            Complex::new(1.0, 0.0),
            Complex::new(0.0, 0.0),
            Complex::new(0.0, 0.0),
            Complex::new(0.0, 0.0),
            Complex::new(0.0, 0.0),
            Complex::new(0.0, 0.0),
        ]);

        assert_eq!(result, expected_state);

        let result = state.apply_single_qubit_gate(2, &matrix);

        let expected_state = State(vec![
            Complex::new(0.0, 0.0),
            Complex::new(0.0, 0.0),
            Complex::new(0.0, 0.0),
            Complex::new(0.0, 0.0),
            Complex::new(1.0, 0.0),
            Complex::new(0.0, 0.0),
            Complex::new(0.0, 0.0),
            Complex::new(0.0, 0.0),
        ]);

        assert_eq!(result, expected_state);
    }
}
