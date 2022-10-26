const A: u64 = 1_103_515_245;
const B: u64 = 12_345;
const K: u64 = 32_767;

#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub struct Rand {
    next: u64,
}

impl Rand {
    pub fn new(seed: u64) -> Self {
        Self { next: seed }
    }
}

impl Iterator for Rand {
    type Item = u64;

    fn next(&mut self) -> Option<Self::Item> {
        self.next = (A * self.next + B) % K;
        Some(self.next)
    }
}

