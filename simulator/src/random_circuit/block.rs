use super::{SingleQubitGate, TwoQubitGate};
use std::fmt;

#[derive(Debug, PartialEq, Copy, Clone)]
pub enum Block {
    SingleQubitGate {
        qubit_idx: usize,
        gate: SingleQubitGate,
    },
    TwoQubitGate {
        control_qubit_idx: usize,
        target_qubit_idx: usize,
        gate: TwoQubitGate,
    },
}

impl fmt::Display for Block {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Block::SingleQubitGate { qubit_idx, gate } => {
                writeln!(f, "Single qubit gate on {} qubit", qubit_idx)?;
                for row in gate.0.iter() {
                    let mut value_iter = row.iter();
                    let value = value_iter.next().unwrap();
                    write!(f, "{}{:+}i", value.re, value.im)?;
                    for value in value_iter {
                        write!(f, "\t{}{:+}i", value.re, value.im)?;
                    }
                    writeln!(f, "")?;
                }
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
                for row in gate.0.iter() {
                    let mut value_iter = row.iter();
                    let value = value_iter.next().unwrap();
                    write!(f, "{}{:+}i", value.re, value.im)?;
                    for value in value_iter {
                        write!(f, "\t{}{:+}i", value.re, value.im)?;
                    }
                    writeln!(f, "")?;
                }
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::single_qubit_gate::consts::*;
    use crate::two_qubit_gate::consts::*;
    use pretty_assertions::assert_eq;

    const X_STRING: &str = "Single qubit gate on 1 qubit\n0+0i\t1+0i\n1+0i\t0+0i\n";
    const CNOT_STRING: &str = "Two qubit gate on 0 and 1 qubits\n1+0i\t0+0i\t0+0i\t0+0i\n0+0i\t1+0i\t0+0i\t0+0i\n0+0i\t0+0i\t0+0i\t1+0i\n0+0i\t0+0i\t1+0i\t0+0i\n";

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
