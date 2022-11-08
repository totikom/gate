use clap::Parser;
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

    let circuit_builder = RandomCircuitIter::new(cli.seed, cli.qubits, cli.gates);

    let mut state = vec![Complex32::new(0.0, 0.0); 1 << cli.qubits];

    state[0] = Complex32::new(1.0, 0.0);
    let mut state = State::new(state);

    println!(
        "Seed =  {}, number of qubits = {}, number of gates {}",
        cli.seed, cli.qubits, cli.gates
    );

    state = state.evaluate_circuit_progress_bar(circuit_builder);

    if !cli.drop_final_state {
        println!("{}", state);
    }
}
