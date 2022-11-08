use super::{SingleQubitGate, TwoQubitGate};
use std::fmt;

#[derive(Debug, PartialEq, Copy, Clone)]
pub enum Block {
    SingleQubitGate {
        qubit_idx: u64,
        gate: SingleQubitGate,
    },
    TwoQubitGate {
        control_qubit_idx: u64,
        target_qubit_idx: u64,
        gate: TwoQubitGate,
    },
}

impl fmt::Display for Block {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Block::SingleQubitGate { qubit_idx, gate } => {
                writeln!(f, "Single qubit gate on {} qubit", qubit_idx)?;
                write!(f, "{}", gate)
            }
            Block::TwoQubitGate {
                control_qubit_idx,
                target_qubit_idx,
                gate,
            } => {
                writeln!(
                    f,
                    "Two qubit gate on {} and {} qubits",
                    control_qubit_idx, target_qubit_idx
                )?;
                write!(f, "{}", gate)
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::single_qubit_gate::gates::*;
    use crate::two_qubit_gate::gates::*;
    use pretty_assertions::assert_eq;

    const X_STRING: &str = "Single qubit gate on 1 qubit\n\
                            0.0000+0.0000i\t1.0000+0.0000i\n\
                            1.0000+0.0000i\t0.0000+0.0000i\n";
    const CNOT_STRING: &str = "Two qubit gate on 0 and 1 qubits\n\
        1.0000+0.0000i\t0.0000+0.0000i\t0.0000+0.0000i\t0.0000+0.0000i\n\
        0.0000+0.0000i\t1.0000+0.0000i\t0.0000+0.0000i\t0.0000+0.0000i\n\
        0.0000+0.0000i\t0.0000+0.0000i\t0.0000+0.0000i\t1.0000+0.0000i\n\
        0.0000+0.0000i\t0.0000+0.0000i\t1.0000+0.0000i\t0.0000+0.0000i\n";

    #[test]
    fn format_single() {
        let result = format!(
            "{}",
            Block::SingleQubitGate {
                qubit_idx: 1,
                gate: X
            }
        );
        assert_eq!(&result, X_STRING);
    }

    #[test]
    fn format_two() {
        let result = format!(
            "{}",
            Block::TwoQubitGate {
                control_qubit_idx: 0,
                target_qubit_idx: 1,
                gate: CNOT,
            }
        );
        assert_eq!(&result, CNOT_STRING);
    }
}
