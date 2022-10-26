use num_complex::Complex32;

pub mod consts;
mod ops;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TwoQubitGate(pub [[Complex32; 4]; 4]);

pub const CNOT: TwoQubitGate = TwoQubitGate([
    [
        Complex32::new(1.0, 0.0),
        Complex32::new(0.0, 0.0),
        Complex32::new(0.0, 0.0),
        Complex32::new(0.0, 0.0),
    ],
    [
        Complex32::new(0.0, 0.0),
        Complex32::new(1.0, 0.0),
        Complex32::new(0.0, 0.0),
        Complex32::new(0.0, 0.0),
    ],
    [
        Complex32::new(0.0, 0.0),
        Complex32::new(0.0, 0.0),
        Complex32::new(0.0, 0.0),
        Complex32::new(1.0, 0.0),
    ],
    [
        Complex32::new(0.0, 0.0),
        Complex32::new(0.0, 0.0),
        Complex32::new(1.0, 0.0),
        Complex32::new(0.0, 0.0),
    ],
]);

pub fn suren_gate(phi_1: f32, phi_2: f32, phi_3: f32, phi_4: f32) -> TwoQubitGate {
    TwoQubitGate([
        [
            Complex32::new(phi_1.cos(), phi_1.sin()),
            Complex32::new(0.0, 0.0),
            Complex32::new(0.0, 0.0),
            Complex32::new(0.0, 0.0),
        ],
        [
            Complex32::new(0.0, 0.0),
            Complex32::new(phi_2.cos(), phi_2.sin()),
            Complex32::new(0.0, 0.0),
            Complex32::new(0.0, 0.0),
        ],
        [
            Complex32::new(0.0, 0.0),
            Complex32::new(0.0, 0.0),
            Complex32::new(phi_3.cos(), phi_3.sin()),
            Complex32::new(0.0, 0.0),
        ],
        [
            Complex32::new(0.0, 0.0),
            Complex32::new(0.0, 0.0),
            Complex32::new(0.0, 0.0),
            Complex32::new(phi_4.cos(), phi_4.sin()),
        ],
    ])
}
