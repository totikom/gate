use num_complex::Complex32;
use simulator::random_circuit::{Block, RandomCircuitIter};
use simulator::State;

fn main() {
    let seed = 0;
    let n_qubits = 4;
    let n_gates = 10;
    let mut circuit_builder = RandomCircuitIter::new(seed, n_qubits);

    let mut state = vec![Complex32::new(0.0, 0.0); 1 << n_qubits];
    state[0] = Complex32::new(1.0, 0.0);
    let mut state = State::new(state);
    println!(
        "Seed =  {}, number of qubits = {}, number of gates {}",
        seed, n_qubits, n_gates
    );
    for _ in 0..n_gates {
        let gate = circuit_builder.next().unwrap();
        //println!("{}", &gate);
        match gate {
            Block::SingleQubitGate { gate, qubit_idx } => {
                state = state.apply_single_qubit_gate(qubit_idx as usize, &gate);
            }
            Block::TwoQubitGate {
                gate,
                control_qubit_idx,
                target_qubit_idx,
            } => {
                state = state.apply_two_qubit_gate(
                    target_qubit_idx as usize,
                    control_qubit_idx as usize,
                    &gate,
                );
            }
        }
    }

    println!("{}", state);
}
