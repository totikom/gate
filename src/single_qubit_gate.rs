use num_complex::{Complex32};
use std::f32::consts::FRAC_1_SQRT_2;

#[derive(Debug, Clone, PartialEq)]
pub struct SingleQubitGate(pub [[Complex32; 2]; 2]);

pub const X: SingleQubitGate = SingleQubitGate([
    [Complex32::new(0.0, 0.0), Complex32::new(1.0, 0.0)],
    [Complex32::new(1.0, 0.0), Complex32::new(0.0, 0.0)],
]);

pub const Y: SingleQubitGate = SingleQubitGate([
    [Complex32::new(0.0, 0.0), Complex32::new(0.0, -1.0)],
    [Complex32::new(0.0, 1.0), Complex32::new(0.0, 0.0)],
]);

pub const Z: SingleQubitGate = SingleQubitGate([
    [Complex32::new(1.0, 0.0), Complex32::new(0.0, 0.0)],
    [Complex32::new(0.0, 0.0), Complex32::new(-1.0, 0.0)],
]);

pub const T: SingleQubitGate = SingleQubitGate([
    [Complex32::new(1.0, 0.0), Complex32::new(0.0, 0.0)],
    [
        Complex32::new(0.0, 0.0),
        Complex32::new(FRAC_1_SQRT_2, FRAC_1_SQRT_2),
    ],
]);

pub const H: SingleQubitGate = SingleQubitGate([
    [
        Complex32::new(FRAC_1_SQRT_2, 0.0),
        Complex32::new(FRAC_1_SQRT_2, 0.0),
    ],
    [
        Complex32::new(FRAC_1_SQRT_2, 0.0),
        Complex32::new(-FRAC_1_SQRT_2, 0.0),
    ],
]);

pub fn r_x(theta: f32) -> SingleQubitGate {
    let arg = theta / 2.0;
    SingleQubitGate([
        [Complex32::new(arg.cos(), 0.0), Complex32::new(0.0, -arg.sin())],
        [Complex32::new(0.0, arg.sin()), Complex32::new(arg.cos(), 0.0)],
    ])
}

pub fn r_y(theta: f32) -> SingleQubitGate {
    let arg = theta / 2.0;
    SingleQubitGate([
        [Complex32::new(arg.cos(), 0.0), Complex32::new(-arg.sin(), 0.0)],
        [Complex32::new(arg.sin(), 0.0), Complex32::new(arg.cos(), 0.0)],
    ])
}

pub fn r_z(theta: f32) -> SingleQubitGate {
    let arg = theta / 2.0;
    SingleQubitGate([
        [Complex32::new(arg.cos(), -arg.sin()), Complex32::new(0.0, 0.0)],
        [Complex32::new(0.0, 0.0), Complex32::new(arg.cos(), arg.sin())],
    ])
}
