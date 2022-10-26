use simulator::random_circuit::Rand;

fn main() {
    let mut rand = Rand::new(0);

    for _ in 0..10 {
        println!("{}", rand.next().unwrap());
    }
}
