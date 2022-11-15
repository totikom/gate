use super::{SingleQubitGate, TwoQubitGate};
use std::fmt;

#[derive(Debug, PartialEq, Copy, Clone)]
pub enum Block<T>
where
    T: std::ops::Deref<Target = [usize]>,
{
    SingleQubitGate {
        qubit_idx: u64,
        gate: SingleQubitGate,
    },
    TwoQubitGate {
        control_qubit_idx: u64,
        target_qubit_idx: u64,
        gate: TwoQubitGate,
    },
    ToffoliGate {
        control_0_qubit_idx: u64,
        control_1_qubit_idx: u64,
        target_qubit_idx: u64,
    },
    CCGate {
        control_0_qubit_idx: u64,
        control_1_qubit_idx: u64,
        target_qubit_idx: u64,
        root_gate: SingleQubitGate,
    },
    NCGate {
        controllers: T,
        ancillas: T,
        target: usize,
        gate: SingleQubitGate,
    },
    GroverDiffusion {
        diffusion_qubits: T,
        ancillas: T,
    },
}

impl<T: fmt::Debug + std::ops::Deref<Target = [usize]>> fmt::Display for Block<T> {
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
            Block::ToffoliGate {
                control_0_qubit_idx,
                control_1_qubit_idx,
                target_qubit_idx,
            } => {
                writeln!(
                    f,
                    "Toffoli gate controlled by {} and {} qubits, target on {}",
                    control_0_qubit_idx, control_1_qubit_idx, target_qubit_idx,
                )
            }
            Block::CCGate {
                control_0_qubit_idx,
                control_1_qubit_idx,
                target_qubit_idx,
                root_gate,
            } => {
                writeln!(
                    f,
                    "Twice-controlled qubit gate on {} qubit, controlled by {} and {}",
                    target_qubit_idx, control_0_qubit_idx, control_1_qubit_idx
                )?;
                write!(f, "{}", root_gate)
            }
            Block::NCGate {
                controllers,
                ancillas,
                target,
                gate,
            } => {
                writeln!(
                    f,
                    "N-controlled gate on {} qubit, controlled by {:?}, ancillas are: {:?}",
                    target, controllers, ancillas,
                )?;
                write!(f, "{}", gate)
            }
            Block::GroverDiffusion {
                diffusion_qubits,
                ancillas,
            } => {
                writeln!(
                    f,
                    "Grover diffusion {:?} qubits, ancillas are: {:?}",
                    diffusion_qubits, ancillas,
                )
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
            Block::<Vec<usize>>::SingleQubitGate {
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
            Block::<Vec<usize>>::TwoQubitGate {
                control_qubit_idx: 0,
                target_qubit_idx: 1,
                gate: CNOT,
            }
        );
        assert_eq!(&result, CNOT_STRING);
    }
}
