use num_complex::Complex32;
use std::f32::consts::FRAC_1_SQRT_2;

use super::SingleQubitGate;

pub const I: SingleQubitGate = SingleQubitGate([
    [Complex32::new(1.0, 0.0), Complex32::new(0.0, 0.0)],
    [Complex32::new(0.0, 0.0), Complex32::new(1.0, 0.0)],
]);

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

pub const SX: SingleQubitGate = SingleQubitGate([
    [Complex32::new(0.5, 0.5), Complex32::new(0.5, -0.5)],
    [Complex32::new(0.5, -0.5), Complex32::new(0.5, 0.5)],
]);

pub const SZ: SingleQubitGate = SingleQubitGate([
    [Complex32::new(1.0, 0.0), Complex32::new(0.0, 0.0)],
    [Complex32::new(0.0, 0.0), Complex32::new(FRAC_1_SQRT_2, FRAC_1_SQRT_2)],
]);

pub fn r_x(theta: f32) -> SingleQubitGate {
    let arg = theta / 2.0;
    SingleQubitGate([
        [
            Complex32::new(arg.cos(), 0.0),
            Complex32::new(0.0, -arg.sin()),
        ],
        [
            Complex32::new(0.0, arg.sin()),
            Complex32::new(arg.cos(), 0.0),
        ],
    ])
}

pub fn r_y(theta: f32) -> SingleQubitGate {
    let arg = theta / 2.0;
    SingleQubitGate([
        [
            Complex32::new(arg.cos(), 0.0),
            Complex32::new(-arg.sin(), 0.0),
        ],
        [
            Complex32::new(arg.sin(), 0.0),
            Complex32::new(arg.cos(), 0.0),
        ],
    ])
}

pub fn r_z(theta: f32) -> SingleQubitGate {
    let arg = theta / 2.0;
    SingleQubitGate([
        [
            Complex32::new(arg.cos(), -arg.sin()),
            Complex32::new(0.0, 0.0),
        ],
        [
            Complex32::new(0.0, 0.0),
            Complex32::new(arg.cos(), arg.sin()),
        ],
    ])
}

pub fn arbitrary_unitary_matrix(a: f32, b: f32, c: f32, d: f32) -> SingleQubitGate {
    let r1 = r_z(b);
    let r2 = r_y(c);
    let r3 = r_z(d);

    r1 * r2 * r3 * Complex32::new(a.cos(), a.sin())
}
