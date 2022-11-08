use num_complex::Complex32;

use super::TwoQubitGate;
use crate::SingleQubitGate;

pub const CNOT: TwoQubitGate = TwoQubitGate([
    [
        Complex32::new(1.0, 0.0),
        Complex32::new(0.0, 0.0),
        Complex32::new(0.0, 0.0),
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
    [
        Complex32::new(0.0, 0.0),
        Complex32::new(0.0, 0.0),
        Complex32::new(1.0, 0.0),
        Complex32::new(0.0, 0.0),
    ],
]);

pub fn suren_gate(phi_1: f32, phi_2: f32, phi_3: f32, phi_4: f32) -> TwoQubitGate {
    TwoQubitGate([
        [
            Complex32::new(phi_1.cos(), phi_1.sin()),
            Complex32::new(0.0, 0.0),
            Complex32::new(0.0, 0.0),
            Complex32::new(0.0, 0.0),
        ],
        [
            Complex32::new(0.0, 0.0),
            Complex32::new(phi_2.cos(), phi_2.sin()),
            Complex32::new(0.0, 0.0),
            Complex32::new(0.0, 0.0),
        ],
        [
            Complex32::new(0.0, 0.0),
            Complex32::new(0.0, 0.0),
            Complex32::new(phi_3.cos(), phi_3.sin()),
            Complex32::new(0.0, 0.0),
        ],
        [
            Complex32::new(0.0, 0.0),
            Complex32::new(0.0, 0.0),
            Complex32::new(0.0, 0.0),
            Complex32::new(phi_4.cos(), phi_4.sin()),
        ],
    ])
}

pub fn controlled_u(u: &SingleQubitGate) -> TwoQubitGate {
    let mut gate = [[Complex32::new(0.0, 0.0); 4]; 4];
    gate[0][0] = Complex32::new(1.0, 0.0);
    gate[1][1] = Complex32::new(1.0, 0.0);

    for (controlled_u_row, u_row) in gate.iter_mut().skip(2).zip(u.0.iter()) {
        controlled_u_row[2..4].copy_from_slice(u_row);
    }
    TwoQubitGate(gate)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::single_qubit_gate::gates::{X, Z};

    use pretty_assertions::assert_eq;

    #[test]
    fn cnot() {
        let result = controlled_u(&X);

        assert_eq!(result, CNOT);
    }

    #[test]
    fn cz() {
        let result = controlled_u(&Z);

        let mut expected_result = [[Complex32::new(0.0, 0.0); 4]; 4];

        for i in 0..3 {
            expected_result[i][i] = Complex32::new(1.0, 0.0);
        }
        expected_result[3][3] = Complex32::new(-1.0, 0.0);

        let expected_result = TwoQubitGate(expected_result);

        assert_eq!(result, expected_result);
    }
}
