use clap::Parser;
use indicatif::{ProgressIterator, ProgressStyle};
use num_complex::Complex32;
use simulator::random_circuit::{Block, RandomCircuitIter};
use simulator::State;

#[derive(Debug, Parser)]
#[command(version, about, long_about = None)]
struct Args {
    #[arg(short, long, value_name = "NUM", default_value_t = 0)]
    seed: u64,
    #[arg(short, long, value_name = "NUM", default_value_t = 3)]
    qubits: u64,
    #[arg(short, long, value_name = "NUM", default_value_t = 10)]
    gates: usize,
    #[arg(short, long)]
    drop_final_state: bool,
}
fn main() {
    let cli = Args::parse();

    let mut circuit_builder = RandomCircuitIter::new(cli.seed, cli.qubits);

    let mut state = vec![Complex32::new(0.0, 0.0); 1 << cli.qubits];

    state[0] = Complex32::new(1.0, 0.0);
    let mut state = State::new(state);

    println!(
        "Seed =  {}, number of qubits = {}, number of gates {}",
        cli.seed, cli.qubits, cli.gates
    );

    for _ in (0..cli.gates).progress().with_style(
        ProgressStyle::with_template(
            "[{elapsed_precise}] {wide_bar:.cyan/blue} {pos:>7}/{len:7} {msg}",
        )
        .unwrap()
        .progress_chars("##-"),
    ) {
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

    if !cli.drop_final_state {
        println!("{}", state);
    }
}
