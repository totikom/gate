use num_complex::Complex32;

pub mod consts;
mod ops;

#[derive(Clone, Copy, PartialEq)]
pub struct TwoQubitGate(pub [[Complex32; 4]; 4]);
