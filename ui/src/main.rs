use simulator::random_circuit::RandomCircuitIter;

fn main() {
    let seed = 0;
    let n_quibits = 4;
    let mut circuit_builder = RandomCircuitIter::new(seed, n_quibits);

    println!("Seed =  {}, number of qubits = {}", seed, n_quibits);
    for _ in 0..10 {
        println!("{}", circuit_builder.next().unwrap());
    }
}
