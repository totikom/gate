use super::*;
use crate::single_qubit_gate::gates::*;
use crate::two_qubit_gate::gates::*;
use num_complex::ComplexFloat;
use std::f32::consts::{FRAC_1_SQRT_2, FRAC_PI_3, PI, SQRT_2};

mod single_qubit_gate {
    use super::*;
    use pretty_assertions::assert_eq;

    #[test]
    fn single_qubit_x() {
        let mut state = State::from_bit_str("0").unwrap();
        let mut temp_state = State::from_bit_str("0").unwrap();

        state.apply_single_qubit_gate(0, &X, &mut temp_state);

        let expected_state = State::from_bit_str("1").unwrap();

        assert_eq!(state, expected_state);
    }

    #[test]
    fn two_qubits_x() {
        let mut state = State::from_bit_str("00").unwrap();
        let mut temp_state = State::from_bit_str("00").unwrap();

        state.apply_single_qubit_gate(0, &X, &mut temp_state);

        let expected_state = State::from_bit_str("01").unwrap();

        assert_eq!(state, expected_state);

        let mut state = State::from_bit_str("10").unwrap();

        state.apply_single_qubit_gate(0, &X, &mut temp_state);

        let expected_state = State::from_bit_str("11").unwrap();

        assert_eq!(state, expected_state);
    }

    #[test]
    fn three_qubits_x() {
        let mut state = State::from_bit_str("000").unwrap();
        let mut temp_state = State::from_bit_str("000").unwrap();

        state.apply_single_qubit_gate(0, &X, &mut temp_state);

        let expected_state = State::from_bit_str("001").unwrap();

        assert_eq!(state, expected_state);

        let mut state = State::from_bit_str("010").unwrap();

        state.apply_single_qubit_gate(0, &X, &mut temp_state);

        let expected_state = State::from_bit_str("011").unwrap();

        assert_eq!(state, expected_state);

        let mut state = State::from_bit_str("000").unwrap();

        state.apply_single_qubit_gate(1, &X, &mut temp_state);

        let expected_state = State::from_bit_str("010").unwrap();

        assert_eq!(state, expected_state);

        let mut state = State::from_bit_str("010").unwrap();

        state.apply_single_qubit_gate(1, &X, &mut temp_state);

        let expected_state = State::from_bit_str("000").unwrap();

        assert_eq!(state, expected_state);
    }

    #[test]
    fn single_qubit_z() {
        let mut qubit = State(vec![Complex32::new(1.0, 0.0), Complex32::new(0.0, 0.0)]);
        let mut tmp_qubit = State(vec![Complex32::new(1.0, 0.0), Complex32::new(0.0, 0.0)]);

        qubit.apply_single_qubit_gate(0, &Z, &mut tmp_qubit);

        let expected_qubit = State(vec![Complex32::new(1.0, 0.0), Complex32::new(0.0, 0.0)]);

        assert_eq!(qubit, expected_qubit);

        let mut qubit = State(vec![Complex32::new(0.0, 0.0), Complex32::new(1.0, 0.0)]);

        qubit.apply_single_qubit_gate(0, &Z, &mut tmp_qubit);

        let expected_qubit = State(vec![Complex32::new(0.0, 0.0), Complex32::new(-1.0, 0.0)]);

        assert_eq!(qubit, expected_qubit);
    }

    #[test]
    fn two_qubits_z() {
        let mut state = State(vec![
            Complex32::new(0.0, 0.0),
            Complex32::new(1.0, 0.0),
            Complex32::new(0.0, 0.0),
            Complex32::new(0.0, 0.0),
        ]);
        let mut tmp_state = State(vec![
            Complex32::new(0.0, 0.0),
            Complex32::new(1.0, 0.0),
            Complex32::new(0.0, 0.0),
            Complex32::new(0.0, 0.0),
        ]);

        state.apply_single_qubit_gate(1, &Z, &mut tmp_state);

        let expected_state = State(vec![
            Complex32::new(0.0, 0.0),
            Complex32::new(1.0, 0.0),
            Complex32::new(0.0, 0.0),
            Complex32::new(0.0, 0.0),
        ]);

        assert_eq!(state, expected_state);

        let mut state = State(vec![
            Complex32::new(0.0, 0.0),
            Complex32::new(1.0, 0.0),
            Complex32::new(0.0, 0.0),
            Complex32::new(0.0, 0.0),
        ]);

        state.apply_single_qubit_gate(0, &Z, &mut tmp_state);

        let expected_state = State(vec![
            Complex32::new(0.0, 0.0),
            Complex32::new(-1.0, 0.0),
            Complex32::new(0.0, 0.0),
            Complex32::new(0.0, 0.0),
        ]);

        assert_eq!(state, expected_state);
    }
}

mod two_qubit_gate {
    use super::*;
    use pretty_assertions::assert_eq;

    #[test]
    fn cnot() {
        let mut state = State::from_bit_str("00").unwrap();
        let mut temp_state = State::from_bit_str("00").unwrap();

        state.apply_two_qubit_gate(0, 1, &CNOT, &mut temp_state);

        let expected_state = State::from_bit_str("00").unwrap();

        assert_eq!(state, expected_state);

        let mut state = State::from_bit_str("01").unwrap();

        state.apply_two_qubit_gate(0, 1, &CNOT, &mut temp_state);

        let expected_state = State::from_bit_str("11").unwrap();

        assert_eq!(state, expected_state);

        let mut state = State::from_bit_str("10").unwrap();

        state.apply_two_qubit_gate(0, 1, &CNOT, &mut temp_state);

        let expected_state = State::from_bit_str("10").unwrap();

        assert_eq!(state, expected_state);

        let mut state = State::from_bit_str("11").unwrap();

        state.apply_two_qubit_gate(0, 1, &CNOT, &mut temp_state);

        let expected_state = State::from_bit_str("01").unwrap();

        assert_eq!(state, expected_state);
    }

    #[test]
    fn suren() {
        let mut state = State::from_bit_str("00").unwrap();
        let mut temp_state = State::from_bit_str("00").unwrap();

        let suren_gate = suren_gate(0.0, 1.0, 2.0, 3.0);

        state.apply_two_qubit_gate(1, 0, &suren_gate, &mut temp_state);

        let expected_state = State(vec![
            Complex32::new(1.0, 0.0),
            Complex32::new(0.0, 0.0),
            Complex32::new(0.0, 0.0),
            Complex32::new(0.0, 0.0),
        ]);

        assert_eq!(state, expected_state);

        let mut state = State(vec![
            Complex32::new(0.0, 0.0),
            Complex32::new(1.0, 0.0),
            Complex32::new(0.0, 0.0),
            Complex32::new(0.0, 0.0),
        ]);

        state.apply_two_qubit_gate(1, 0, &suren_gate, &mut temp_state);

        let expected_state = State(vec![
            Complex32::new(0.0, 0.0),
            Complex32::new(1.0_f32.cos(), 1.0_f32.sin()),
            Complex32::new(0.0, 0.0),
            Complex32::new(0.0, 0.0),
        ]);

        assert_eq!(state, expected_state);

        let mut state = State(vec![
            Complex32::new(0.0, 0.0),
            Complex32::new(0.0, 0.0),
            Complex32::new(1.0, 0.0),
            Complex32::new(0.0, 0.0),
        ]);

        state.apply_two_qubit_gate(1, 0, &suren_gate, &mut temp_state);

        let expected_state = State(vec![
            Complex32::new(0.0, 0.0),
            Complex32::new(0.0, 0.0),
            Complex32::new(2.0_f32.cos(), 2.0_f32.sin()),
            Complex32::new(0.0, 0.0),
        ]);

        assert_eq!(state, expected_state);

        let mut state = State(vec![
            Complex32::new(0.0, 0.0),
            Complex32::new(0.0, 0.0),
            Complex32::new(0.0, 0.0),
            Complex32::new(1.0, 0.0),
        ]);

        state.apply_two_qubit_gate(1, 0, &suren_gate, &mut temp_state);

        let expected_state = State(vec![
            Complex32::new(0.0, 0.0),
            Complex32::new(0.0, 0.0),
            Complex32::new(0.0, 0.0),
            Complex32::new(3.0_f32.cos(), 3.0_f32.sin()),
        ]);

        assert_eq!(state, expected_state);
    }
}

#[test]
fn toffoli() {
    let mut temp_state = State::from_bit_str("000").unwrap();
    for i in 0..7 {
        let control_0 = i & 1;
        let control_1 = (i & 2) >> 1;
        let target = (i & 4) >> 2;
        let mut state =
            State::from_bit_str(&format!("{}{}{}", target, control_1, control_0)).unwrap();
        let control_0_bool = control_0 == 1;
        let control_1_bool = control_1 == 1;
        let target_bool = target == 1;

        let expected_result = u32::from(target_bool ^ (control_0_bool & control_1_bool));

        let expected_state =
            State::from_bit_str(&format!("{}{}{}", expected_result, control_1, control_0)).unwrap();

        state.apply_toffoli_gate(0, 1, 2, &mut temp_state);
        assert!(state.distance(&expected_state) < 1e-4);
    }
}

#[test]
fn ccnot() {
    let mut temp_state = State::from_bit_str("000").unwrap();
    for i in 0..7 {
        let control_0 = i & 1;
        let control_1 = (i & 2) >> 1;
        let target = (i & 4) >> 2;
        let mut state =
            State::from_bit_str(&format!("{}{}{}", target, control_1, control_0)).unwrap();
        let control_0_bool = control_0 == 1;
        let control_1_bool = control_1 == 1;
        let target_bool = target == 1;

        let expected_result = u32::from(target_bool ^ (control_0_bool & control_1_bool));

        let expected_state =
            State::from_bit_str(&format!("{}{}{}", expected_result, control_1, control_0)).unwrap();

        state.apply_controlled_controlled_gate(0, 1, 2, &SX, &mut temp_state);
        assert!(state.distance(&expected_state) < 1e-4);
    }
}

#[test]
fn cccnot() {
    let mut temp_state = State::from_bit_str("000000").unwrap();
    for i in 0..15 {
        let control_0 = i & 1;
        let control_1 = (i & 2) >> 1;
        let control_2 = (i & 4) >> 2;
        let target = (i & 8) >> 3;
        let mut state = State::from_bit_str(&dbg!(format!(
            "00{}{}{}{}",
            target, control_2, control_1, control_0
        )))
        .unwrap();
        let control_0_bool = control_0 == 1;
        let control_1_bool = control_1 == 1;
        let control_2_bool = control_2 == 1;
        let target_bool = target == 1;

        let expected_result =
            u32::from(target_bool ^ (control_0_bool & control_1_bool & control_2_bool));

        let expected_state = State::from_bit_str(&dbg!(format!(
            "00{}{}{}{}",
            expected_result, control_2, control_1, control_0
        )))
        .unwrap();

        state.apply_n_controlled_gate(vec![0, 1, 2], vec![4, 5], 3, &X, &mut temp_state);
        assert!(dbg!(state).distance(&dbg!(expected_state)) < 1e-4);
    }
}

mod str {
    use super::*;
    use pretty_assertions::assert_eq;

    #[test]
    fn state_000() {
        let result = State::from_bit_str("000").unwrap();

        let expected_state = State(vec![
            Complex32::new(1.0, 0.0),
            Complex32::new(0.0, 0.0),
            Complex32::new(0.0, 0.0),
            Complex32::new(0.0, 0.0),
            Complex32::new(0.0, 0.0),
            Complex32::new(0.0, 0.0),
            Complex32::new(0.0, 0.0),
            Complex32::new(0.0, 0.0),
        ]);

        assert_eq!(result, expected_state);
    }

    #[test]
    fn state_001() {
        let result = State::from_bit_str("001").unwrap();

        let expected_state = State(vec![
            Complex32::new(0.0, 0.0),
            Complex32::new(1.0, 0.0),
            Complex32::new(0.0, 0.0),
            Complex32::new(0.0, 0.0),
            Complex32::new(0.0, 0.0),
            Complex32::new(0.0, 0.0),
            Complex32::new(0.0, 0.0),
            Complex32::new(0.0, 0.0),
        ]);

        assert_eq!(result, expected_state);
    }

    #[test]
    fn state_00_1() {
        let result = State::from_bit_str("00_1").unwrap();

        let expected_state = State(vec![
            Complex32::new(0.0, 0.0),
            Complex32::new(1.0, 0.0),
            Complex32::new(0.0, 0.0),
            Complex32::new(0.0, 0.0),
            Complex32::new(0.0, 0.0),
            Complex32::new(0.0, 0.0),
            Complex32::new(0.0, 0.0),
            Complex32::new(0.0, 0.0),
        ]);

        assert_eq!(result, expected_state);
    }

    #[test]
    fn state_10() {
        let result = State::from_bit_str("10").unwrap();

        let expected_state = State(vec![
            Complex32::new(0.0, 0.0),
            Complex32::new(0.0, 0.0),
            Complex32::new(1.0, 0.0),
            Complex32::new(0.0, 0.0),
        ]);

        assert_eq!(result, expected_state);
    }
}

mod grover {
    use super::*;

    #[test]
    fn grover_naive() {
        let expected_state = State::from_bit_str("001111").unwrap();
        let mut state = State::from_bit_str("000000").unwrap();
        let mut temp_state = State::from_bit_str("000000").unwrap();
        //                                  00_0000
        //                                  2 ancillas
        //                                  4 qubits

        let n = (FRAC_PI_4 * (2_u32.pow(4) as f32).sqrt()).floor() as usize;

        for i in 0..4 {
            state.apply_single_qubit_gate(i, &H, &mut temp_state);
        }

        for i in 0..n {
            state.apply_n_controlled_gate(vec![0, 1, 2], vec![4, 5], 3, &Z, &mut temp_state);

            state.grover_diffusion(vec![0, 1, 2, 3], vec![4, 5], &mut temp_state);
            println!("{} {}", i, state.scalar_product(&expected_state).abs());
        }

        let theta = (2.0 * 15.0.sqrt() / 16.0).asin();
        let probabity = ((n as f32 + 0.5) * theta).sin().powi(2);

        assert!(
            (dbg!(state.scalar_product(&expected_state).abs().powi(2)) - dbg!(probabity)).abs()
                <= 1e-3
        );
    }

    #[test]
    fn grover() {
        let expected_state = State::from_bit_str("00_1_111").unwrap();
        let mut state = State::from_bit_str("00_0_000").unwrap();
        let mut temp_state = State::from_bit_str("00_0_000").unwrap();
        //                                  00_0000
        //                                  2 ancillas
        //                                  4 qubits

        let oracle = vec![Block::NCGate {
            controllers: vec![0, 1, 2],
            ancillas: vec![4, 5],
            target: 3,
            gate: Z,
        }];

        state.apply_grover_algorithm(
            vec![0, 1, 2, 3],
            vec![4, 5],
            1,
            oracle.into_iter(),
            &mut temp_state,
        );

        let N = 2.0_f32.powi(4);
        let M = 1.0;
        let k = (FRAC_PI_4 * (N / M).sqrt()).floor();
        let probabity = expected_probability(N, M, k);

        assert!(
            (dbg!(state.scalar_product(&expected_state).abs().powi(2)) - dbg!(probabity)).abs()
                <= 1e-3
        );
    }

    #[test]
    fn grover_2_answers() {
        let expected_state_0 = State::from_bit_str("000_1_111_0").unwrap();
        let expected_state_1 = State::from_bit_str("000_1_111_1").unwrap();
        let expected_state = (expected_state_1 + expected_state_0) * FRAC_1_SQRT_2;

        let mut state = State::from_bit_str("000_0_000_0").unwrap();
        let mut temp_state = State::from_bit_str("000_0_000_0").unwrap();
        //                                  00_0000
        //                                  2 ancillas
        //                                  4 qubits

        let oracle = vec![Block::NCGate {
            controllers: vec![1, 2, 3],
            ancillas: vec![5, 6],
            target: 4,
            gate: Z,
        }];

        state.apply_grover_algorithm(
            vec![0, 1, 2, 3, 4],
            vec![5, 6, 7],
            2,
            oracle.into_iter(),
            &mut temp_state,
        );

        let N = 2.0_f32.powi(5);
        let M = 2.0;
        let k = (FRAC_PI_4 * (N / M).sqrt()).floor();
        let probabity = expected_probability(N, M, k);

        for i in 0..10 {
            print!("{} ", expected_probability(N, M, i as f32));
        }

        println!();

        assert!(
            (dbg!(state.scalar_product(&expected_state).abs().powi(2)) - dbg!(probabity)).abs()
                <= 1e-3
        );
    }

    #[test]
    fn grover_2_answers_6q() {
        let expected_state_0 = State::from_bit_str("000000_1_1111_0").unwrap();
        let expected_state_1 = State::from_bit_str("000000_1_1111_1").unwrap();
        let expected_state = (expected_state_1 + expected_state_0) * FRAC_1_SQRT_2;

        let mut state = State::from_bit_str("000000_0_0000_0").unwrap();
        let mut temp_state = State::from_bit_str("000000_0_0000_0").unwrap();

        let oracle = vec![Block::NCGate {
            controllers: vec![1, 2, 3, 4],
            ancillas: vec![6, 7, 8],
            target: 5,
            gate: Z,
        }];

        state.apply_grover_algorithm(
            vec![0, 1, 2, 3, 4, 5],
            vec![6, 7, 8, 9],
            2,
            oracle.into_iter(),
            &mut temp_state,
        );

        let N = 2.0_f32.powi(6);
        let M = 2.0;
        let k = (FRAC_PI_4 * (N / M).sqrt()).floor();
        let probabity = expected_probability(N, M, k);

        for i in 0..10 {
            print!("{} ", expected_probability(N, M, i as f32));
        }

        println!();

        assert!(
            (dbg!(state.scalar_product(&expected_state).abs().powi(2)) - dbg!(probabity)).abs()
                <= 1e-3
        );
    }

    #[test]
    fn grover_4_answers() {
        let expected_state_0 = State::from_bit_str("0000_1_111_00").unwrap();
        let expected_state_1 = State::from_bit_str("0000_1_111_01").unwrap();
        let expected_state_2 = State::from_bit_str("0000_1_111_10").unwrap();
        let expected_state_3 = State::from_bit_str("0000_1_111_11").unwrap();
        let expected_state =
            (expected_state_0 + expected_state_1 + expected_state_2 + expected_state_3) * 0.5;

        let mut state = State::from_bit_str("0000_0_000_00").unwrap();
        let mut temp_state = State::from_bit_str("0000_0_000_00").unwrap();
        //                                  00_0000
        //                                  2 ancillas
        //                                  4 qubits

        let oracle = vec![Block::NCGate {
            controllers: vec![2, 3, 4],
            ancillas: vec![6, 7],
            target: 5,
            gate: Z,
        }];

        state.apply_grover_algorithm(
            vec![0, 1, 2, 3, 4, 5],
            vec![6, 7, 8, 9],
            4,
            oracle.into_iter(),
            &mut temp_state,
        );

        let N = 2.0_f32.powi(6);
        let M = 4.0;
        let k = (FRAC_PI_4 * (N / M).sqrt()).floor();
        let probabity = expected_probability(N, M, k);

        for i in 0..k as usize + 5 {
            print!("{} ", expected_probability(N, M, i as f32));
        }

        println!();

        assert!(
            (dbg!(state.scalar_product(&expected_state).abs().powi(2)) - dbg!(probabity)).abs()
                <= 1e-3
        );
    }

    fn expected_probability(n_variants: f32, m_answers: f32, k_iterations: f32) -> f32 {
        let theta = (2.0 * (m_answers * (n_variants - m_answers)).sqrt() / n_variants).asin();
        ((2.0 * k_iterations + 1.0) / 2.0 * theta).sin().powi(2)
    }
}

mod qft {
    use super::*;

    #[test]
    fn two_qubit_in_out() {
        let mut tmp = State::from_bit_str("00").unwrap();
        for i in 0..4 {
            let control_0 = i & 1;
            let control_1 = (i & 2) >> 1;
            let mut state = State::from_bit_str(&format!("{}{}", control_0, control_1)).unwrap();
            let expected_state = state.clone();

            state.apply_qft(0, 1, &mut tmp);
            state.apply_inv_qft(0, 1, &mut tmp);
            assert!(dbg!(state).distance(&dbg!(expected_state)) < 1e-4);
        }
    }

    #[test]
    fn three_qubit_in_out() {
        let mut tmp = State::from_bit_str("000").unwrap();
        for i in 0..8 {
            let digits = to_digits(i, 3);

            let digits_str = to_digit_str(&digits);
            println!("{}", &digits_str);

            let mut state = State::from_bit_str(&digits_str).unwrap();
            let expected_state = state.clone();

            state.apply_qft(0, 2, &mut tmp);
            state.apply_inv_qft(0, 2, &mut tmp);

            assert!(dbg!(state).distance(&dbg!(expected_state)) < 1e-4);
        }
    }

    #[test]
    fn four_qubit_in_out() {
        let mut tmp = State::from_bit_str("0000").unwrap();
        for i in 0..16 {
            let digits = to_digits(i, 4);

            let digits_str = to_digit_str(&digits);
            println!("{}", &digits_str);

            let mut state = State::from_bit_str(&digits_str).unwrap();
            let expected_state = state.clone();

            state.apply_qft(0, 3, &mut tmp);
            state.apply_inv_qft(0, 3, &mut tmp);

            assert!(dbg!(state).distance(&dbg!(expected_state)) < 1e-4);
        }
    }

    #[test]
    fn one_qubit() {
        let mut tmp = State::from_bit_str("0").unwrap();
        for i in 0..2 {
            let digits = to_digits(i, 1);

            let digits_str = to_digit_str(&digits);
            println!("{}", &digits_str);

            let mut state = State::from_bit_str(&digits_str).unwrap();

            state.apply_qft(0, 0, &mut tmp);

            let phases = to_phases(&digits);
            dbg!(&phases);

            let expected_states = to_states(&phases);

            let expected_state = to_state(expected_states);

            assert!(dbg!(state).distance(&dbg!(expected_state)) < 1e-4);
        }
    }

    #[test]
    fn one_qubit_inv() {
        let mut tmp = State::from_bit_str("0").unwrap();
        for i in 0..2 {
            let digits = to_digits(i, 1);

            let digits_str = to_digit_str(&digits);
            println!("{}", &digits_str);

            let expected_state = State::from_bit_str(&digits_str).unwrap();

            let phases = to_phases(&digits);
            dbg!(&phases);

            let states = to_states(&phases);

            let mut state = to_state(states);

            state.apply_qft(0, 0, &mut tmp);

            assert!(dbg!(state).distance(&dbg!(expected_state)) < 1e-4);
        }
    }

    #[test]
    fn two_qubit() {
        let mut tmp = State::from_bit_str("00").unwrap();
        for i in 1..4 {
            let digits = to_digits(i, 2);

            let digits_str = to_digit_str(&digits);
            println!("{}", &digits_str);

            let mut state = State::from_bit_str(&digits_str).unwrap();

            state.apply_qft(0, 1, &mut tmp);

            let phases = to_phases(&digits);
            dbg!(&phases);

            let expected_states = to_states(&phases);
            dbg!(&expected_states);

            let expected_state = to_state(expected_states);

            assert!(dbg!(state).distance(&dbg!(expected_state)) < 1e-4);
        }
    }

    #[test]
    fn two_qubit_inv() {
        let mut tmp = State::from_bit_str("00").unwrap();
        for i in 0..4 {
            let digits = to_digits(i, 2);

            let digits_str = to_digit_str(&digits);
            println!("{}", &digits_str);

            let expected_state = State::from_bit_str(&digits_str).unwrap();

            let phases = to_phases(&digits);
            dbg!(&phases);

            let states = to_states(&phases);

            let mut state = to_state(states);

            state.apply_inv_qft(0, 1, &mut tmp);

            assert!(dbg!(state).distance(&dbg!(expected_state)) < 1e-4);
        }
    }

    #[test]
    fn three_qubit() {
        let mut tmp = State::from_bit_str("000").unwrap();
        for i in 0..8 {
            let digits = to_digits(i, 3);

            let digits_str = to_digit_str(&digits);
            println!("{}", &digits_str);

            let mut state = State::from_bit_str(&digits_str).unwrap();

            state.apply_qft(0, 2, &mut tmp);

            let phases = to_phases(&digits);
            dbg!(&phases);

            let expected_states = to_states(&phases);
            dbg!(&expected_states);

            let expected_state = to_state(expected_states);

            assert!(dbg!(state).distance(&dbg!(expected_state)) < 1e-4);
        }
    }

    #[test]
    fn four_qubit() {
        let mut tmp = State::from_bit_str("0000").unwrap();
        for i in 0..16 {
            let digits = to_digits(i, 4);

            let digits_str = to_digit_str(&digits);
            println!("{}", &digits_str);

            let mut state = State::from_bit_str(&digits_str).unwrap();

            state.apply_qft(0, 3, &mut tmp);

            let phases = to_phases(&digits);
            dbg!(&phases);

            let expected_states = to_states(&phases);

            let expected_state = to_state(expected_states);

            assert!(dbg!(state).distance(&dbg!(expected_state)) < 1e-4);
        }
    }

    #[test]
    fn four_qubit_inv() {
        let mut tmp = State::from_bit_str("0000").unwrap();
        for i in 0..16 {
            let digits = to_digits(i, 4);

            let digits_str = to_digit_str(&digits);
            println!("{}", &digits_str);

            let expected_state = State::from_bit_str(&digits_str).unwrap();

            let phases = to_phases(&digits);
            dbg!(&phases);

            let states = to_states(&phases);

            let mut state = to_state(states);

            state.apply_inv_qft(0, 3, &mut tmp);

            assert!(dbg!(state).distance(&dbg!(expected_state)) < 1e-4);
        }
    }

    fn to_digits(num: usize, size: usize) -> Vec<u8> {
        (0..size)
            .into_iter()
            .map(|i| ((num & (1 << i)) >> i) as u8)
            .collect()
    }

    fn to_digit_str(digits: &Vec<u8>) -> String {
        digits.iter().rev().map(|i| format!("{}", i)).collect()
    }

    fn to_phases(digits: &Vec<u8>) -> Vec<f32> {
        let mut phase = 0.0;
        let mut phases = Vec::new();
        for &digit in digits.iter().rev() {
            phase = (phase + digit as f32) * 0.5;
            phases.push(phase);
        }
        phases
    }

    fn to_states(phases: &Vec<f32>) -> Vec<State> {
        phases
            .iter()
            .map(|phase| {
                let phase = TAU * phase;
                (State::from_bit_str("0").unwrap()
                    + State::from_bit_str("1").unwrap() * Complex32::new(phase.cos(), phase.sin()))
                    * FRAC_1_SQRT_2
            })
            .collect()
    }

    fn to_state(states: Vec<State>) -> State {
        states
            .into_iter()
            .reduce(|accum, state| state.kron(&accum))
            .unwrap()
    }
}

#[test]
fn kron() {
    let state1 = State::from_bit_str("0").unwrap();
    let state2 = State::from_bit_str("1").unwrap();
    let state_res = State::from_bit_str("01").unwrap();
    assert_eq!(state1.kron(&state2), state_res);

    let state1 = State::from_bit_str("1").unwrap();
    let state2 = State::from_bit_str("1").unwrap();
    let state_res = State::from_bit_str("11").unwrap();
    assert_eq!(state1.kron(&state2), state_res);

    let state1 = State::from_bit_str("1").unwrap();
    let state2 = State::from_bit_str("0").unwrap();
    let state_res = State::from_bit_str("10").unwrap();
    assert_eq!(state1.kron(&state2), state_res);

    let state1 = State::from_bit_str("000").unwrap();
    let state2 = State::from_bit_str("00").unwrap();
    let state_res = State::from_bit_str("000_00").unwrap();
    dbg!(state1.0.len());
    dbg!(state2.0.len());
    let result = state1.kron(&state2);
    assert_eq!(result.0.len(), state_res.0.len());
    assert_eq!(result, state_res);

    let state1 = State::from_bit_str("100").unwrap();
    let state2 = State::from_bit_str("01").unwrap();
    let state_res = State::from_bit_str("100_01").unwrap();
    dbg!(state1.0.len());
    dbg!(state2.0.len());
    let result = state1.kron(&state2);
    assert_eq!(result.0.len(), state_res.0.len());
    assert_eq!(result, state_res);
}

mod phase_estimation {
    use super::*;
    use pretty_assertions::assert_eq;

    #[test]
    #[ignore]// too long
    fn three_qubit_op_1() {
        let expected_phase = 0.5042053881557941;

        let eigen_vec = vec![
            Complex32::new(0.06327676, -0.11673236),
            Complex32::new(-0.00090105, 0.35390646),
            Complex32::new(0.23446311, -0.13134822),
            Complex32::new(-0.3327006, 0.29969659),
            Complex32::new(-0.22010026, -0.12466702),
            Complex32::new(0.69412754, 0.),
            Complex32::new(0.0316633, -0.10823434),
            Complex32::new(-0.15695969, 0.03511926),
        ];
        let eig = State::new(eigen_vec);

        for n_estimation in 1..6 {
            let estimation_zeroes = std::iter::repeat("0")
                .take(n_estimation)
                .collect::<String>();
            let estimation_state = State::from_bit_str(&estimation_zeroes).unwrap();
            let mut temp_state =
                State::from_bit_str(&format!("00_000_{}", estimation_zeroes)).unwrap();
            //                                               estimation
            //                                           eig
            //                                         ancilla

            let circuit: Vec<Block<Vec<usize>>> = vec![
                // 1 layer
                Block::SingleQubitGate {
                    gate: H,
                    qubit_idx: n_estimation + 0,
                },
                Block::SingleQubitGate {
                    gate: H,
                    qubit_idx: n_estimation + 1,
                },
                Block::SingleQubitGate {
                    gate: H,
                    qubit_idx: n_estimation + 2,
                },
                // 2 layer
                Block::SingleQubitGate {
                    gate: X,
                    qubit_idx: n_estimation + 0,
                },
                Block::SingleQubitGate {
                    gate: SX,
                    qubit_idx: n_estimation + 1,
                },
                Block::SingleQubitGate {
                    gate: T,
                    qubit_idx: n_estimation + 2,
                },
                // 3 layer
                Block::ToffoliGate {
                    control_0_qubit_idx: n_estimation + 0,
                    control_1_qubit_idx: n_estimation + 1,
                    target_qubit_idx: n_estimation + 2,
                },
            ];

            //dbg!(&eig);
            let ancilla = State::from_bit_str("00").unwrap();

            let initial_state = eig.kron(&estimation_state);
            let mut initial_state = ancilla.kron(&initial_state);
            //dbg!(&initial_state);

            initial_state.apply_phase_estimation_algorithm(
                0,
                n_estimation - 1,
                n_estimation + 3,
                n_estimation + 3 + 1,
                circuit.into_iter(),
                &mut temp_state,
            );
            //dbg!(&initial_state);

            let probs = initial_state.prob_reduce(n_estimation);
            //dbg!(&probs);

            let value = probs
                .iter()
                .enumerate()
                .max_by(|(_, a), (_, b)| a.total_cmp(b))
                .map(|(index, _)| index)
                .unwrap();

            let phase = value as f32 / 2.0.powi(n_estimation as i32);
            let delta = (phase - expected_phase).abs();

            dbg!(phase);
            dbg!(expected_phase);
            assert!(delta < 2.0.powi(-(n_estimation as i32 + 1)));
        }
    }

    #[test]
    #[ignore]
    fn three_qubit_op_2() {
        let expected_phase = 0.3068456619416857;

        let eigen_vec = vec![
            Complex32::new(0.20775065, -0.10035436),
            Complex32::new(-0.32302179, 0.29197462),
            Complex32::new(-0.17216106, -0.13002442),
            Complex32::new(0.69480461, 0.),
            Complex32::new(0.07218861, 0.01968729),
            Complex32::new(0.13538043, 0.105178),
            Complex32::new(0.30927723, 0.0552938),
            Complex32::new(-0.30467237, -0.03678417),
        ];
        let eig = State::new(eigen_vec);

        for n_estimation in 1..6 {
            //let n_estimation = 5;
            let estimation_zeroes = std::iter::repeat("0")
                .take(n_estimation)
                .collect::<String>();
            let estimation_state = State::from_bit_str(&estimation_zeroes).unwrap();
            let mut temp_state =
                State::from_bit_str(&format!("00_000_{}", estimation_zeroes)).unwrap();
            //                                               estimation
            //                                           eig
            //                                         ancilla

            let circuit: Vec<Block<Vec<usize>>> = vec![
                // 1 layer
                Block::SingleQubitGate {
                    gate: H,
                    qubit_idx: n_estimation + 0,
                },
                Block::SingleQubitGate {
                    gate: H,
                    qubit_idx: n_estimation + 1,
                },
                Block::SingleQubitGate {
                    gate: H,
                    qubit_idx: n_estimation + 2,
                },
                // 2 layer
                Block::SingleQubitGate {
                    gate: X,
                    qubit_idx: n_estimation + 0,
                },
                Block::SingleQubitGate {
                    gate: SX,
                    qubit_idx: n_estimation + 1,
                },
                Block::SingleQubitGate {
                    gate: T,
                    qubit_idx: n_estimation + 2,
                },
                // 3 layer
                Block::ToffoliGate {
                    control_0_qubit_idx: n_estimation + 0,
                    control_1_qubit_idx: n_estimation + 1,
                    target_qubit_idx: n_estimation + 2,
                },
            ];

            //dbg!(&eig);
            let ancilla = State::from_bit_str("00").unwrap();

            let initial_state = eig.kron(&estimation_state);
            let mut initial_state = ancilla.kron(&initial_state);
            //dbg!(&initial_state);

            initial_state.apply_phase_estimation_algorithm(
                0,
                n_estimation - 1,     // last estimation qubit
                n_estimation + 3,     //ancilla 1
                n_estimation + 3 + 1, //ancilla 2
                circuit.into_iter(),
                &mut temp_state,
            );
            //dbg!(&initial_state);

            let probs = initial_state.prob_reduce(n_estimation);
            //dbg!(&probs);

            let value = probs
                .iter()
                .enumerate()
                .max_by(|(_, a), (_, b)| a.total_cmp(b))
                .map(|(index, _)| index)
                .unwrap();
            //dbg!(&value);

            let mut prob_value_pairs: Vec<(usize, f32)> =
                probs.iter().enumerate().map(|(a, b)| (a, *b)).collect();
            prob_value_pairs.sort_by(|(_, val), (_, val2)| val.partial_cmp(val2).unwrap());
            //dbg!(prob_value_pairs);

            let phase = value as f32 / 2.0.powi(n_estimation as i32);
            let delta = (phase - expected_phase).abs();

            dbg!(phase);
            dbg!(expected_phase);
            assert!(delta < 2.0.powi(-(n_estimation as i32)));
        }
    }

    #[test]
    //#[ignore]
    fn one_qubit_X_1() {
        let expected_phase = 0.0;

        let eigen_vec = vec![
            Complex32::new(FRAC_1_SQRT_2, 0.),
            Complex32::new(FRAC_1_SQRT_2, 0.),
        ];
        let eig = State::new(eigen_vec);

        for n_estimation in 1..10 {
            let estimation_zeroes = std::iter::repeat("0")
                .take(n_estimation)
                .collect::<String>();
            let estimation_state = State::from_bit_str(&estimation_zeroes).unwrap();
            let mut temp_state =
                State::from_bit_str(&format!("00_0_{}", estimation_zeroes)).unwrap();
            //                                               estimation
            //                                           eig
            //                                         ancilla

            let circuit: Vec<Block<Vec<usize>>> = vec![
                // 1 layer
                Block::SingleQubitGate {
                    gate: H,
                    qubit_idx: n_estimation + 0,
                },
            ];

            //dbg!(&eig);
            let ancilla = State::from_bit_str("00").unwrap();

            let initial_state = eig.kron(&estimation_state);
            let mut initial_state = ancilla.kron(&initial_state);
            //dbg!(&initial_state);

            initial_state.apply_phase_estimation_algorithm(
                0,                    // estimation_start
                n_estimation - 1,     // estimation_end
                n_estimation + 1,     // ancilla 1
                n_estimation + 1 + 1, //ancilla 2
                circuit.into_iter(),
                &mut temp_state,
            );
            //dbg!(&initial_state);

            let probs = initial_state.prob_reduce(n_estimation);
            //dbg!(&probs);

            let value = probs
                .iter()
                .enumerate()
                .max_by(|(_, a), (_, b)| a.total_cmp(b))
                .map(|(index, _)| index)
                .unwrap();

            let phase = value as f32 / 2.0.powi(n_estimation as i32);
            let delta = (phase - expected_phase).abs();

            dbg!(phase);
            dbg!(expected_phase);
            assert!(delta < 2.0.powi(-(n_estimation as i32 + 1)));
        }
    }

    #[test]
    //#[ignore]
    fn one_qubit_X_2() {
        let expected_phase = 0.5;

        let eigen_vec = vec![
            Complex32::new(FRAC_1_SQRT_2, 0.),
            Complex32::new(-FRAC_1_SQRT_2, 0.),
        ];
        let eig = State::new(eigen_vec);

        //for n_estimation in 1..10 {
        let n_estimation = 2;
            let estimation_zeroes = std::iter::repeat("0")
                .take(n_estimation)
                .collect::<String>();
            let estimation_state = State::from_bit_str(&estimation_zeroes).unwrap();
            let mut temp_state =
                State::from_bit_str(&format!("00_0_{}", estimation_zeroes)).unwrap();
            //                                               estimation
            //                                           eig
            //                                         ancilla

            let circuit: Vec<Block<Vec<usize>>> = vec![
                // 1 layer
                Block::SingleQubitGate {
                    gate: H,
                    qubit_idx: n_estimation + 0,
                },
            ];

            //dbg!(&eig);
            let ancilla = State::from_bit_str("00").unwrap();

            let initial_state = eig.kron(&estimation_state);
            let mut initial_state = ancilla.kron(&initial_state);
            //dbg!(&initial_state);

            initial_state.apply_phase_estimation_algorithm(
                0,                    // estimation_start
                n_estimation - 1,     // estimation_end
                n_estimation + 1,     // ancilla 1
                n_estimation + 1 + 1, //ancilla 2
                circuit.into_iter(),
                &mut temp_state,
            );
            //dbg!(&initial_state);

            let probs = initial_state.prob_reduce(n_estimation);
            //dbg!(&probs);

            let value = probs
                .iter()
                .enumerate()
                .max_by(|(_, a), (_, b)| a.total_cmp(b))
                .map(|(index, _)| index)
                .unwrap();

            let phase = value as f32 / 2.0.powi(n_estimation as i32);
            let delta = (phase - expected_phase).abs();

            dbg!(phase);
            dbg!(expected_phase);
            //assert!(delta < 2.0.powi(-(n_estimation as i32 + 1)));
        //}
        todo!();
    }

    #[test]
    //#[ignore]
    fn one_qubit_op_2() {
        let expected_phase = 1.0;

        let eigen_vec = vec![
            Complex32::new(-0.38268343, 0.),
            Complex32::new(0.92387953, 0.),
        ];
        let eig = State::new(eigen_vec);

        for n_estimation in 1..6 {
            let estimation_zeroes = std::iter::repeat("0")
                .take(n_estimation)
                .collect::<String>();
            let estimation_state = State::from_bit_str(&estimation_zeroes).unwrap();
            let mut temp_state =
                State::from_bit_str(&format!("00_0_{}", estimation_zeroes)).unwrap();
            //                                               estimation
            //                                           eig
            //                                         ancilla

            let circuit: Vec<Block<Vec<usize>>> = vec![
                // 1 layer
                Block::SingleQubitGate {
                    gate: H,
                    qubit_idx: n_estimation + 0,
                },
            ];

            //dbg!(&eig);
            let ancilla = State::from_bit_str("00").unwrap();

            let initial_state = eig.kron(&estimation_state);
            let mut initial_state = ancilla.kron(&initial_state);
            //dbg!(&initial_state);

            initial_state.apply_phase_estimation_algorithm(
                0,                    // estimation_start
                n_estimation - 1,     // estimation_end
                n_estimation + 1,     // ancilla 1
                n_estimation + 1 + 1, //ancilla 2
                circuit.into_iter(),
                &mut temp_state,
            );
            //dbg!(&initial_state);

            let probs = initial_state.prob_reduce(n_estimation);
            //dbg!(&probs);

            let value = probs
                .iter()
                .enumerate()
                .max_by(|(_, a), (_, b)| a.total_cmp(b))
                .map(|(index, _)| index)
                .unwrap();

            let phase = value as f32 / 2.0.powi(n_estimation as i32);
            let delta = (phase - expected_phase).abs();

            dbg!(phase);
            dbg!(expected_phase);
            assert!(delta <= 2.0.powi(-(n_estimation as i32)));
        }
    }

    #[test]
    #[ignore]
    fn four_qubit_op_1() {
        let expected_phase = 0.413154598834622;

        let eigen_vec = vec![
            Complex32::new(-0.03542531, -0.14944166),
            Complex32::new(-0.07025411, 0.05198089),
            Complex32::new(0.12899156, 0.15777689),
            Complex32::new(-0.14034703, -0.22655787),
            Complex32::new(-0.32270646, -0.01983551),
            Complex32::new(0.04293359, 0.06855498),
            Complex32::new(0.44154969, 0.),
            Complex32::new(0.22142507, 0.07148536),
            Complex32::new(-0.11146848, 0.04651739),
            Complex32::new(0.21147643, 0.18366894),
            Complex32::new(-0.09914681, -0.27989018),
            Complex32::new(-0.28971536, 0.06116787),
            Complex32::new(0.03876883, 0.08057337),
            Complex32::new(0.37341838, -0.03809144),
            Complex32::new(0.11467434, 0.12426874),
            Complex32::new(0.00216074, -0.22284398),
        ];
        let eig = State::new(eigen_vec);

        for n_estimation in 1..6 {
            //let n_estimation = 5;
            let estimation_zeroes = std::iter::repeat("0")
                .take(n_estimation)
                .collect::<String>();
            let estimation_state = State::from_bit_str(&estimation_zeroes).unwrap();
            let mut temp_state =
                State::from_bit_str(&format!("00_0000_{}", estimation_zeroes)).unwrap();
            //                                               estimation
            //                                           eig
            //                                         ancilla

            let circuit: Vec<Block<Vec<usize>>> = vec![
                // 1 layer
                Block::SingleQubitGate {
                    gate: H,
                    qubit_idx: n_estimation + 0,
                },
                Block::SingleQubitGate {
                    gate: H,
                    qubit_idx: n_estimation + 1,
                },
                Block::SingleQubitGate {
                    gate: H,
                    qubit_idx: n_estimation + 2,
                },
                Block::SingleQubitGate {
                    gate: H,
                    qubit_idx: n_estimation + 3,
                },
                // 2 layer
                Block::SingleQubitGate {
                    gate: X,
                    qubit_idx: n_estimation + 0,
                },
                Block::SingleQubitGate {
                    gate: SX,
                    qubit_idx: n_estimation + 1,
                },
                Block::SingleQubitGate {
                    gate: T,
                    qubit_idx: n_estimation + 2,
                },
                Block::SingleQubitGate {
                    gate: X,
                    qubit_idx: n_estimation + 3,
                },
                // 3 layer
                Block::ToffoliGate {
                    control_0_qubit_idx: n_estimation + 1,
                    control_1_qubit_idx: n_estimation + 2,
                    target_qubit_idx: n_estimation + 3,
                },
                // 4 layer
                Block::ToffoliGate {
                    control_0_qubit_idx: n_estimation + 0,
                    control_1_qubit_idx: n_estimation + 1,
                    target_qubit_idx: n_estimation + 2,
                },
            ];

            //dbg!(&eig);
            let ancilla = State::from_bit_str("00").unwrap();

            let initial_state = eig.kron(&estimation_state);
            let mut initial_state = ancilla.kron(&initial_state);
            //dbg!(&initial_state);

            initial_state.apply_phase_estimation_algorithm(
                0,
                n_estimation - 1,
                n_estimation + 4,
                n_estimation + 4 + 1,
                circuit.into_iter(),
                &mut temp_state,
            );
            //dbg!(&initial_state);

            let probs = initial_state.prob_reduce(n_estimation);
            //dbg!(&probs);

            let value = probs
                .iter()
                .enumerate()
                .max_by(|(_, a), (_, b)| a.total_cmp(b))
                .map(|(index, _)| index)
                .unwrap();
            //dbg!(&value);

            let mut prob_value_pairs: Vec<(usize, f32)> =
                probs.iter().enumerate().map(|(a, b)| (a, *b)).collect();
            prob_value_pairs.sort_by(|(_, val), (_, val2)| val.partial_cmp(val2).unwrap());
            //dbg!(prob_value_pairs);

            let phase = value as f32 / 2.0.powi(n_estimation as i32);
            let delta = (phase - expected_phase).abs();

            dbg!(phase);
            dbg!(expected_phase);
            //assert!(delta < 2.0.powi(-(n_estimation as i32)));
        }
        todo!();
    }
}
