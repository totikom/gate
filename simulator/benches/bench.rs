use criterion::{black_box, criterion_group, criterion_main, BatchSize, BenchmarkId, Criterion};
use num_complex::Complex32;
use simulator::{random_circuit::RandomCircuitIter, State};
use std::fmt;

#[derive(Copy, Clone, Debug)]
struct Args {
    seed: u64,
    qubits: u64,
    gates: usize,
}
impl fmt::Display for Args {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(
            f,
            "Seed =  {}, number of qubits = {}, number of gates {}",
            self.seed, self.qubits, self.gates
        )
    }
}

fn random_circuit(c: &mut Criterion) {
    let args = Args {
        seed: 0,
        qubits: 4,
        gates: 10,
    };
    let mut group = c.benchmark_group("rand");

    group.measurement_time(std::time::Duration::new(10,0)).sample_size(100);

    group.bench_with_input(
        BenchmarkId::new("random_circuit", "with 4 qubits and 10 gates"),
        &args,
        |b, &args| {
            let mut state = vec![Complex32::new(0.0, 0.0); 1 << args.qubits];
            state[0] = Complex32::new(1.0, 0.0);
            let tmp_state = vec![Complex32::new(0.0, 0.0); 1 << args.qubits];

            let state = State::new(state);
            let mut temp_state = State::new(tmp_state);

            b.iter_batched(
                || {
                    (
                        state.clone(),
                        RandomCircuitIter::new(args.seed, args.qubits, args.gates),
                    )
                },
                |(state, circuit_builder)| {
                    black_box(state).evaluate_circuit(black_box(circuit_builder), &mut temp_state)
                },
                BatchSize::SmallInput,
            )
        },
    );
}

criterion_group!(benches, random_circuit);
criterion_main!(benches);
