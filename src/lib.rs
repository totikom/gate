use ndarray::Array1;
use num_complex::{Complex, Complex64};

type Amp = Complex<f32>;

#[derive(Debug, Clone, PartialEq)]
struct State(Array1<Amp>);

#[derive(Debug, Clone, PartialEq)]
struct SingleQubitGate([[Amp; 2]; 2]);

impl State {
    pub fn apply_single_qubit_gate_to_qubit(&self, index: usize, gate: &SingleQubitGate) -> Self {
        let mut result = Array1::<Amp>::zeros(self.0.dim());

        for i in 0..self.0.dim() {
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
    use ndarray::array;
    use pretty_assertions::assert_eq;

    #[test]
    fn single_qubit() {
        let qubit = State(array![Complex::new(1.0, 0.0), Complex::new(0.0, 0.0)]);

        let matrix = SingleQubitGate([
            [Complex::new(0.0, 0.0), Complex::new(1.0, 0.0)],
            [Complex::new(1.0, 0.0), Complex::new(0.0, 0.0)],
        ]);

        let result = qubit.apply_single_qubit_gate_to_qubit(0, &matrix);

        let expected_qubit = State(array![Complex::new(0.0, 0.0), Complex::new(1.0, 0.0)]);

        assert_eq!(result, expected_qubit);
    }

    #[test]
    fn two_qubits() {
        let state = State(array![
            Complex::new(1.0, 0.0),
            Complex::new(0.0, 0.0),
            Complex::new(0.0, 0.0),
            Complex::new(0.0, 0.0)
        ]);

        let matrix = SingleQubitGate([
            [Complex::new(0.0, 0.0), Complex::new(1.0, 0.0)],
            [Complex::new(1.0, 0.0), Complex::new(0.0, 0.0)],
        ]);

        let result = state.apply_single_qubit_gate_to_qubit(1, &matrix);

        let expected_state = State(array![
            Complex::new(0.0, 0.0),
            Complex::new(0.0, 0.0),
            Complex::new(1.0, 0.0),
            Complex::new(0.0, 0.0)
        ]);

        assert_eq!(result, expected_state);

        let result = state.apply_single_qubit_gate_to_qubit(0, &matrix);

        let expected_state = State(array![
            Complex::new(0.0, 0.0),
            Complex::new(1.0, 0.0),
            Complex::new(0.0, 0.0),
            Complex::new(0.0, 0.0)
        ]);

        assert_eq!(result, expected_state);
    }

    #[test]
    fn three_qubits() {
        let state = State(array![
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

        let result = state.apply_single_qubit_gate_to_qubit(0, &matrix);

        let expected_state = State(array![
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

        let result = state.apply_single_qubit_gate_to_qubit(1, &matrix);

        let expected_state = State(array![
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

        let result = state.apply_single_qubit_gate_to_qubit(2, &matrix);

        let expected_state = State(array![
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
