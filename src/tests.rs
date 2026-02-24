use super::*;

mod qubitlayer_tests {
    use super::*;

    #[test]
    fn test_add_qubit_layers_owned() {
        let mut lhs = QubitLayer::new(1);
        lhs.main[1] = Complex::new(2.0, 0.0);

        let mut rhs = QubitLayer::new(1);
        rhs.main[1] = Complex::new(3.0, 0.0);

        let sum = lhs + rhs;
        assert_eq!(Complex::new(2.0, 0.0), sum.main[0]);
        assert_eq!(Complex::new(5.0, 0.0), sum.main[1]);
    }

    #[test]
    fn test_add_qubit_layers_borrowed() {
        let mut lhs = QubitLayer::new(1);
        lhs.main[0] = Complex::new(0.5, 0.0);
        lhs.main[1] = Complex::new(1.5, 0.0);

        let mut rhs = QubitLayer::new(1);
        rhs.main[0] = Complex::new(1.5, 0.0);
        rhs.main[1] = Complex::new(0.5, 0.0);

        let sum = &lhs + &rhs;
        assert_eq!(Complex::new(2.0, 0.0), sum.main[0]);
        assert_eq!(Complex::new(2.0, 0.0), sum.main[1]);
    }

    #[test]
    fn test_add_assign_qubit_layers_borrowed() {
        let mut lhs = QubitLayer::new(1);
        lhs.main[0] = Complex::new(0.5, 0.0);
        lhs.main[1] = Complex::new(1.0, 0.0);

        let mut rhs = QubitLayer::new(1);
        rhs.main[0] = Complex::new(1.5, 0.0);
        rhs.main[1] = Complex::new(2.0, 0.0);

        lhs += &rhs;
        assert_eq!(Complex::new(2.0, 0.0), lhs.main[0]);
        assert_eq!(Complex::new(3.0, 0.0), lhs.main[1]);
    }

    #[test]
    fn test_add_assign_qubit_layers_owned() {
        let mut lhs = QubitLayer::new(1);
        lhs.main[1] = Complex::new(2.0, 0.0);

        let mut rhs = QubitLayer::new(1);
        rhs.main[1] = Complex::new(3.0, 0.0);

        lhs += rhs;
        assert_eq!(Complex::new(2.0, 0.0), lhs.main[0]);
        assert_eq!(Complex::new(5.0, 0.0), lhs.main[1]);
    }

    #[test]
    #[should_panic(expected = "Cannot add QubitLayers with different numbers of qubits")]
    fn test_add_qubit_layers_size_mismatch_panics() {
        let lhs = QubitLayer::new(1);
        let rhs = QubitLayer::new(2);
        let _ = lhs + rhs;
    }

    #[test]
    fn test_execute_shots() {
        let instructions = vec![(QuantumOp::Hadamard, 0), (QuantumOp::PauliX, 1)];

        let mut accumulated_layer = QubitLayer::new(3);
        for _ in 0..4 {
            let mut shot_layer = QubitLayer::new(3);
            let result = shot_layer.execute_noiseless(instructions.clone());
            assert!(result.is_ok());
            accumulated_layer += &shot_layer;
        }
        accumulated_layer.scale_amplitudes(4);

        let measured = accumulated_layer.measure_qubits();
        assert_eq!(0.5, (measured[0] * 10.0).round() / 10.0);
        assert_eq!(1.0, (measured[1] * 10.0).round() / 10.0);
        assert_eq!(0.0, (measured[2] * 10.0).round() / 10.0);
    }

    #[test]
    fn test_execute_shots_zero() {
        let mut q_layer = QubitLayer::new(3);
        let instructions = vec![(QuantumOp::Hadamard, 0)];
        let noise = NoiseModel {
            gate_error_prob: 0.0,
            readout_flip_prob: 0.0,
        };

        let result = q_layer.execute_noisy_shots(&instructions, 0, noise);
        assert!(result.is_ok());

        let measured = q_layer.measure_qubits();
        assert_eq!(0.0, measured[0]);
        assert_eq!(0.0, measured[1]);
        assert_eq!(0.0, measured[2]);
    }

    #[test]
    fn test_execute_shots_failed_execute() {
        let mut q_layer = QubitLayer::new(3);
        let instructions = vec![(QuantumOp::PauliX, 10)];

        let result = q_layer.execute_noiseless(instructions);
        assert!(result.is_err());
    }

    #[test]
    fn test_execute_noisy_shots_zero_noise() {
        let instructions = vec![(QuantumOp::Hadamard, 0), (QuantumOp::PauliX, 1)];
        let noise = NoiseModel {
            gate_error_prob: 0.0,
            readout_flip_prob: 0.0,
        };

        let mut accumulated_layer = QubitLayer::new(3);
        let result = accumulated_layer.execute_noisy_shots(&instructions, 3, noise);
        assert!(result.is_ok());

        let measured = accumulated_layer.measure_qubits();
        assert_eq!(0.5, (measured[0] * 10.0).round() / 10.0);
        assert_eq!(1.0, (measured[1] * 10.0).round() / 10.0);
        assert_eq!(0.0, (measured[2] * 10.0).round() / 10.0);
    }

    #[test]
    fn test_execute_noisy_shots_readout_flip_full() {
        let instructions: Vec<(QuantumOp, TargetQubit)> = vec![];
        let noise = NoiseModel {
            gate_error_prob: 0.0,
            readout_flip_prob: 1.0,
        };

        let mut accumulated_layer = QubitLayer::new(3);
        let result = accumulated_layer.execute_noisy_shots(&instructions, 2, noise);
        assert!(result.is_ok());

        let measured = accumulated_layer.measure_qubits();
        assert_eq!(1.0, measured[0]);
        assert_eq!(1.0, measured[1]);
        assert_eq!(1.0, measured[2]);
    }

    #[test]
    fn test_execute_noisy_shots_invalid_noise() {
        let mut q_layer = QubitLayer::new(1);
        let instructions = vec![(QuantumOp::PauliX, 0)];
        let noise = NoiseModel {
            gate_error_prob: 1.1,
            readout_flip_prob: 0.0,
        };

        let result = q_layer.execute_noisy_shots(&instructions, 1, noise);
        assert!(result.is_err());
    }

    #[test]
    fn test_random_executions() {
        let mut q_layer: QubitLayer = QubitLayer::new(3);
        let instructions = vec![
            (QuantumOp::Hadamard, 0),
            (QuantumOp::Hadamard, 1),
            (QuantumOp::Hadamard, 2),
            (QuantumOp::Hadamard, 0),
            (QuantumOp::Hadamard, 1),
            (QuantumOp::Hadamard, 2),
        ];
        let _ = q_layer.execute_noiseless(instructions);
        for it in 0..q_layer.get_num_qubits() {
            assert_eq!(
                0.0,
                (q_layer.measure_qubits()[it as usize] * 10.0).round() / 10.0
            );
        }
    }

    #[test]
    fn test_measure_qubits() {
        let mut q_layer: QubitLayer = QubitLayer::new(3);
        for it in 0..q_layer.get_num_qubits() {
            assert_eq!(0.0, q_layer.measure_qubits()[it as usize]);
        }

        let instructions: Vec<(QuantumOp, TargetQubit)> = vec![];
        let _ = q_layer.execute_noiseless(instructions);

        for it in 0..q_layer.get_num_qubits() {
            assert_eq!(0.0, q_layer.measure_qubits()[it as usize]);
        }
    }

    #[test]
    fn test_spins_on_superposition() {
        let mut q_layer: QubitLayer = QubitLayer::new(3);
        let instructions = vec![
            (QuantumOp::Hadamard, 0),
            (QuantumOp::Hadamard, 1),
            (QuantumOp::Hadamard, 2),
            (QuantumOp::PauliX, 0),
            (QuantumOp::PauliY, 1),
            (QuantumOp::PauliZ, 2),
        ];
        let _ = q_layer.execute_noiseless(instructions);
        for it in 0..q_layer.get_num_qubits() {
            assert_eq!(
                0.5,
                (q_layer.measure_qubits()[it as usize] * 10.0).round() / 10.0
            );
        }
    }

    #[test]
    fn test_failed_execute() {
        let mut q_layer: QubitLayer = QubitLayer::new(10);
        let instructions = vec![(QuantumOp::PauliX, 10)]; // index goes up to 9

        let result: Result<(), String> = q_layer.execute_noiseless(instructions);
        assert!(result.is_err());

        let result: Result<(), String> =
            q_layer.execute_noiseless(vec![(QuantumOp::Hadamard, 2112)]);
        assert!(result.is_err());
    }

    #[test]
    fn test_execute() {
        let mut q_layer: QubitLayer = QubitLayer::new(10);
        let instructions = vec![
            (QuantumOp::PauliX, 0),
            (QuantumOp::PauliY, 1),
            (QuantumOp::PauliZ, 2),
        ];

        if let Err(e) = q_layer.execute_noiseless(instructions) {
            panic!("Should not panic!. Error: {e}");
        }

        assert_eq!(1.0, q_layer.measure_qubits()[0].round());
    }

    #[test]
    fn test_get_num_qubits() {
        let num_qubits = 10;
        let q_layer: QubitLayer = QubitLayer::new(num_qubits);
        assert_eq!(num_qubits, q_layer.get_num_qubits());
    }

    #[test]
    fn test_get_mem_usage() {
        let num_qubits = 20;
        let q_layer: QubitLayer = QubitLayer::new(num_qubits);

        let expected = (8_u64 * 2_u64) * (2_u64.pow(num_qubits));
        assert_eq!(expected, q_layer.get_mem_usage());
    }

    #[test]
    fn test_display_trait_print() {
        let q_layer: QubitLayer = QubitLayer::new(2);
        let expected = "1+0i 0+0i 0+0i 0+0i";
        println!("{}", q_layer);
        assert_eq!(expected, format!("{}", q_layer));
    }

    #[test]
    fn test_debug_trait_print() {
        let q_layer: QubitLayer = QubitLayer::new(2);
        let expected = "|0>\t -> 1+0i\n|1>\t -> 0+0i\n|10>\t -> 0+0i\n|11>\t -> 0+0i\n";
        println!("{:?}", q_layer);
        assert_eq!(expected, format!("{:?}", q_layer));
    }

    #[test]
    fn test_hadamard_simple() {
        let hadamard_const = 1.0 / std::f64::consts::SQRT_2;
        let num_qubits = 3;
        let mut q_layer: QubitLayer = QubitLayer::new(num_qubits);
        q_layer.hadamard(0);
        let mut test_vec = vec![Complex::new(0.0, 0.0); 2_usize.pow(num_qubits)];
        test_vec[0] = Complex::new(hadamard_const, 0.0);
        test_vec[1] = Complex::new(hadamard_const, 0.0);
        assert_eq!(test_vec, q_layer.main);

        let mut q_layer: QubitLayer = QubitLayer::new(num_qubits);
        q_layer.hadamard(2);
        let mut test_vec = vec![Complex::new(0.0, 0.0); 2_usize.pow(num_qubits)];
        test_vec[0] = Complex::new(hadamard_const, 0.0);
        test_vec[4] = Complex::new(hadamard_const, 0.0);
        assert_eq!(test_vec, q_layer.main);
        q_layer = QubitLayer::new(num_qubits);
        q_layer.hadamard(0);
        q_layer.hadamard(1);
        q_layer.hadamard(2);
        let test_vec = vec![Complex::new(pow(hadamard_const, 3), 0.0); 2_usize.pow(num_qubits)];
        assert_eq!(test_vec, q_layer.main);
    }

    #[test]
    fn test_pauli_z_simple() {
        let num_qubits = 3;
        let mut q_layer: QubitLayer = QubitLayer::new(num_qubits);
        q_layer.pauli_z(0);
        let mut test_vec = vec![Complex::new(0.0, 0.0); 2_usize.pow(num_qubits)];
        test_vec[0] = Complex::new(1.0, 0.0);
        assert_eq!(test_vec, q_layer.main);

        let mut q_layer: QubitLayer = QubitLayer::new(num_qubits);
        q_layer.pauli_z(2);
        let mut test_vec = vec![Complex::new(0.0, 0.0); 2_usize.pow(num_qubits)];
        test_vec[0] = Complex::new(1.0, 0.0);
        assert_eq!(test_vec, q_layer.main);
        q_layer = QubitLayer::new(num_qubits);
        q_layer.pauli_z(0);
        q_layer.pauli_z(1);
        q_layer.pauli_z(2);
        let mut test_vec = vec![Complex::new(0.0, 0.0); 2_usize.pow(num_qubits)];
        test_vec[0] = Complex::new(1.0, 0.0);
        assert_eq!(test_vec, q_layer.main);
    }

    #[test]
    fn test_pauli_y_simple() {
        let num_qubits = 3;
        let mut q_layer: QubitLayer = QubitLayer::new(num_qubits);
        q_layer.pauli_y(0);
        let mut test_vec = vec![Complex::new(0.0, 0.0); 2_usize.pow(num_qubits)];
        test_vec[1] = Complex::new(0.0, 1.0);
        assert_eq!(test_vec, q_layer.main);

        let mut q_layer: QubitLayer = QubitLayer::new(num_qubits);
        q_layer.pauli_y(2);
        let mut test_vec = vec![Complex::new(0.0, 0.0); 2_usize.pow(num_qubits)];
        test_vec[4] = Complex::new(0.0, 1.0);
        assert_eq!(test_vec, q_layer.main);
        q_layer = QubitLayer::new(num_qubits);
        q_layer.pauli_y(0);
        q_layer.pauli_y(1);
        q_layer.pauli_y(2);
        let mut test_vec = vec![Complex::new(0.0, 0.0); 2_usize.pow(num_qubits)];
        test_vec[7] = Complex::new(0.0, -1.0);
        assert_eq!(test_vec, q_layer.main);
    }

    #[test]
    fn test_pauli_x_simple() {
        let num_qubits = 3;
        let mut q_layer: QubitLayer = QubitLayer::new(num_qubits);
        q_layer.pauli_x(0);
        let mut test_vec = vec![Complex::new(0.0, 0.0); 2_usize.pow(num_qubits)];
        test_vec[1] = Complex::new(1.0, 0.0);
        assert_eq!(test_vec, q_layer.main);

        let mut q_layer: QubitLayer = QubitLayer::new(num_qubits);
        q_layer.pauli_x(2);
        let mut test_vec = vec![Complex::new(0.0, 0.0); 2_usize.pow(num_qubits)];
        test_vec[4] = Complex::new(1.0, 0.0);
        assert_eq!(test_vec, q_layer.main);
        q_layer = QubitLayer::new(num_qubits);
        q_layer.pauli_x(0);
        q_layer.pauli_x(1);
        q_layer.pauli_x(2);
        let mut test_vec = vec![Complex::new(0.0, 0.0); 2_usize.pow(num_qubits)];
        test_vec[7] = Complex::new(1.0, 0.0);
        assert_eq!(test_vec, q_layer.main);
    }

    #[test]
    fn test_s_gate_simple() {
        let hadamard_const = 1.0 / std::f64::consts::SQRT_2;
        let mut q_layer: QubitLayer = QubitLayer::new(1);

        let result = q_layer.execute_noiseless(vec![(QuantumOp::Hadamard, 0), (QuantumOp::S, 0)]);
        assert!(result.is_ok());

        let expected = vec![
            Complex::new(hadamard_const, 0.0),
            Complex::new(0.0, hadamard_const),
        ];
        assert_eq!(expected, q_layer.main);
    }

    #[test]
    fn test_t_gate_simple() {
        let hadamard_const = 1.0 / std::f64::consts::SQRT_2;
        let mut q_layer: QubitLayer = QubitLayer::new(1);

        let result = q_layer.execute_noiseless(vec![(QuantumOp::Hadamard, 0), (QuantumOp::T, 0)]);
        assert!(result.is_ok());

        let expected = vec![
            Complex::new(hadamard_const, 0.0),
            Complex::from_polar(hadamard_const, std::f64::consts::FRAC_PI_4),
        ];
        assert_eq!(expected, q_layer.main);
    }

    #[test]
    fn test_controlled_x_simple() {
        let num_qubits = 2;
        let mut q_layer: QubitLayer = QubitLayer::new(num_qubits);

        q_layer.pauli_x(1);
        q_layer.controlled_x(1, 0);

        let measured = q_layer.measure_qubits();
        assert_eq!(1.0, measured[0]);
        assert_eq!(1.0, measured[1]);
    }

    #[test]
    fn test_controlled_z_simple() {
        let num_qubits = 2;
        let mut q_layer: QubitLayer = QubitLayer::new(num_qubits);

        q_layer.pauli_x(0);
        q_layer.pauli_x(1);
        q_layer.controlled_z(1, 0);

        let measured = q_layer.measure_qubits();
        assert_eq!(1.0, measured[0]);
        assert_eq!(1.0, measured[1]);
    }

    #[test]
    fn test_toffoli_simple() {
        let num_qubits = 3;
        let mut q_layer: QubitLayer = QubitLayer::new(num_qubits);

        q_layer.pauli_x(0);
        q_layer.pauli_x(1);
        q_layer.toffoli(0, 1, 2);

        let measured = q_layer.measure_qubits();
        assert_eq!(1.0, measured[0]);
        assert_eq!(1.0, measured[1]);
        assert_eq!(1.0, measured[2]);
    }

    #[test]
    fn test_toffoli_no_flip_when_not_all_controls_active() {
        let num_qubits = 3;
        let mut q_layer: QubitLayer = QubitLayer::new(num_qubits);

        q_layer.pauli_x(0);
        q_layer.toffoli(0, 1, 2);

        let measured = q_layer.measure_qubits();
        assert_eq!(1.0, measured[0]);
        assert_eq!(0.0, measured[1]);
        assert_eq!(0.0, measured[2]);
    }

    #[test]
    fn test_execute_noiseless_toffoli_instruction() {
        let mut q_layer: QubitLayer = QubitLayer::new(3);

        let prep = vec![(QuantumOp::PauliX, 0), (QuantumOp::PauliX, 1)];
        let prep_result = q_layer.execute_noiseless(prep);
        assert!(prep_result.is_ok());

        let toffoli_result = q_layer.execute_noiseless(vec![(TwoCtrlQubitOp::Toffoli, 0, 1, 2)]);
        assert!(toffoli_result.is_ok());

        let measured = q_layer.measure_qubits();
        assert_eq!(1.0, measured[0]);
        assert_eq!(1.0, measured[1]);
        assert_eq!(1.0, measured[2]);
    }

    #[test]
    fn test_execute_noiseless_toffoli_out_of_range() {
        let mut q_layer: QubitLayer = QubitLayer::new(3);

        let result = q_layer.execute_noiseless(vec![(TwoCtrlQubitOp::Toffoli, 0, 7, 2)]);
        assert!(result.is_err());
    }

    #[test]
    fn test_execute_controlled_x_file() {
        let qasm = std::fs::read_to_string("qasm_files/Controllex_X_5_qubits.openqasm")
            .expect("should read Controllex_X_5_qubits.openqasm");

        let parsed = openq3_parser::parse(&qasm).expect("parser should parse cx qasm");
        let mut q_layer = QubitLayer::new(parsed.num_qubits);

        let result = q_layer.execute_noiseless(parsed.ops.clone());
        assert!(result.is_ok());

        let measured = q_layer.measure_qubits();
        assert_eq!(parsed.num_qubits as usize, measured.len());
        for value in measured {
            assert_eq!(0.5, (value * 10.0).round() / 10.0);
        }
    }

    #[test]
    fn test_execute_controlled_z_file() {
        let qasm = std::fs::read_to_string("qasm_files/Controllex_Z_5_qubits.openqasm")
            .expect("should read Controllex_Z_5_qubits.openqasm");

        let parsed = openq3_parser::parse(&qasm).expect("parser should parse cz qasm");
        let mut q_layer = QubitLayer::new(parsed.num_qubits);

        let result = q_layer.execute_noiseless(parsed.ops.clone());
        assert!(result.is_ok());

        let measured = q_layer.measure_qubits();
        assert_eq!(parsed.num_qubits as usize, measured.len());
        for value in measured {
            assert_eq!(0.5, (value * 10.0).round() / 10.0);
        }
    }

    #[test]
    fn test_execute_grovers_file() {
        let qasm = std::fs::read_to_string("qasm_files/Grovers_5_qubits.openqasm")
            .expect("should read Grovers_5_qubits.openqasm");

        let parsed = openq3_parser::parse(&qasm).expect("parser should parse Grover qasm");
        let mut q_layer = QubitLayer::new(parsed.num_qubits);

        let result = q_layer.execute_noiseless(parsed.ops.clone());
        assert!(result.is_ok());

        let measured = q_layer.measure_qubits();
        assert_eq!(parsed.num_qubits as usize, measured.len());
        assert!(measured.iter().all(|value| (0.0..=1.0).contains(value)));
    }
}

mod parser_tests {
    use super::SingleCtrlQubitOp;
    use super::TwoCtrlQubitOp;
    use super::*;
    use std::fs;

    #[test]
    fn parse_valid() {
        let qasm = "qreg q[3];\nx q[0];\ny q[1];\nz q[2];\nh q[0];\ns q[1];\nt q[2];";

        let parsed = openq3_parser::parse(qasm).expect("parser should parse valid qasm");

        assert_eq!(3, parsed.num_qubits);
        assert_eq!(6, parsed.ops.len());
        assert_eq!(QInstruct::Single((QuantumOp::PauliX, 0)), parsed.ops[0]);
        assert_eq!(QInstruct::Single((QuantumOp::PauliY, 1)), parsed.ops[1]);
        assert_eq!(QInstruct::Single((QuantumOp::PauliZ, 2)), parsed.ops[2]);
        assert_eq!(QInstruct::Single((QuantumOp::Hadamard, 0)), parsed.ops[3]);
        assert_eq!(QInstruct::Single((QuantumOp::S, 1)), parsed.ops[4]);
        assert_eq!(QInstruct::Single((QuantumOp::T, 2)), parsed.ops[5]);
    }

    #[test]
    fn parse_controlled_x_ops() {
        let qasm = fs::read_to_string("qasm_files/Controllex_X_5_qubits.openqasm")
            .expect("should read Controllex_X_5_qubits.openqasm");

        let parsed = openq3_parser::parse(&qasm).expect("parser should parse cx qasm");

        assert_eq!(5, parsed.num_qubits);
        assert_eq!(7, parsed.ops.len());
        assert_eq!(QInstruct::Single((QuantumOp::Hadamard, 0)), parsed.ops[0]);
        assert_eq!(QInstruct::Single((QuantumOp::Hadamard, 1)), parsed.ops[1]);
        assert_eq!(QInstruct::Single((QuantumOp::Hadamard, 2)), parsed.ops[2]);
        assert_eq!(QInstruct::Single((QuantumOp::Hadamard, 3)), parsed.ops[3]);
        assert_eq!(QInstruct::Single((QuantumOp::Hadamard, 4)), parsed.ops[4]);
        assert_eq!(
            QInstruct::SingleCtrl((SingleCtrlQubitOp::ControlledX, 0, 1)),
            parsed.ops[5]
        );
        assert_eq!(
            QInstruct::SingleCtrl((SingleCtrlQubitOp::ControlledX, 2, 3)),
            parsed.ops[6]
        );
    }

    #[test]
    fn parse_controlled_z_ops() {
        let qasm = fs::read_to_string("qasm_files/Controllex_Z_5_qubits.openqasm")
            .expect("should read Controllex_Z_5_qubits.openqasm");

        let parsed = openq3_parser::parse(&qasm).expect("parser should parse cz qasm");

        assert_eq!(5, parsed.num_qubits);
        assert_eq!(7, parsed.ops.len());
        assert_eq!(QInstruct::Single((QuantumOp::Hadamard, 0)), parsed.ops[0]);
        assert_eq!(QInstruct::Single((QuantumOp::Hadamard, 1)), parsed.ops[1]);
        assert_eq!(QInstruct::Single((QuantumOp::Hadamard, 2)), parsed.ops[2]);
        assert_eq!(QInstruct::Single((QuantumOp::Hadamard, 3)), parsed.ops[3]);
        assert_eq!(QInstruct::Single((QuantumOp::Hadamard, 4)), parsed.ops[4]);
        assert_eq!(
            QInstruct::SingleCtrl((SingleCtrlQubitOp::ControlledZ, 0, 1)),
            parsed.ops[5]
        );
        assert_eq!(
            QInstruct::SingleCtrl((SingleCtrlQubitOp::ControlledZ, 2, 3)),
            parsed.ops[6]
        );
    }

    #[test]
    fn parse_grovers_file() {
        let qasm = fs::read_to_string("qasm_files/Grovers_5_qubits.openqasm")
            .expect("should read Grovers_5_qubits.openqasm");

        let parsed = openq3_parser::parse(&qasm).expect("parser should parse Grover qasm");

        assert!(parsed.num_qubits == 5);
        assert!(!parsed.ops.is_empty());
        assert_eq!(QInstruct::Single((QuantumOp::Hadamard, 0)), parsed.ops[0]);

        let has_toffoli = parsed
            .ops
            .iter()
            .any(|op| matches!(op, QInstruct::TwoCtrl((TwoCtrlQubitOp::Toffoli, _, _, _))));
        assert!(has_toffoli);
    }

    #[test]
    fn parse_multiops_3_file() {
        let qasm = fs::read_to_string("qasm_files/MultiOps_3_qubits.openqasm")
            .expect("should read MultiOps_3_qubits.openqasm");

        let parsed = openq3_parser::parse(&qasm).expect("parser should parse MultiOps_3 qasm");

        assert_eq!(3, parsed.num_qubits);
        assert_eq!(9, parsed.ops.len());
        assert!(parsed.ops.iter().any(|op| matches!(
            op,
            QInstruct::SingleCtrl((SingleCtrlQubitOp::ControlledX, _, _))
        )));
        assert!(parsed.ops.iter().any(|op| matches!(
            op,
            QInstruct::SingleCtrl((SingleCtrlQubitOp::ControlledZ, _, _))
        )));
    }

    #[test]
    fn parse_multiops_4_file() {
        let qasm = fs::read_to_string("qasm_files/MultiOps_4_qubits.openqasm")
            .expect("should read MultiOps_4_qubits.openqasm");

        let parsed = openq3_parser::parse(&qasm).expect("parser should parse MultiOps_4 qasm");

        assert_eq!(4, parsed.num_qubits);
        assert_eq!(10, parsed.ops.len());
        assert!(parsed
            .ops
            .iter()
            .any(|op| matches!(op, QInstruct::Single((QuantumOp::PauliY, 0)))));
        assert!(parsed.ops.iter().any(|op| matches!(
            op,
            QInstruct::SingleCtrl((SingleCtrlQubitOp::ControlledX, 0, 2))
        )));
        assert!(parsed.ops.iter().any(|op| matches!(
            op,
            QInstruct::SingleCtrl((SingleCtrlQubitOp::ControlledZ, 1, 3))
        )));
    }

    #[test]
    fn parse_multiops_5_file() {
        let qasm = fs::read_to_string("qasm_files/MultiOps_5_qubits.openqasm")
            .expect("should read MultiOps_5_qubits.openqasm");

        let parsed = openq3_parser::parse(&qasm).expect("parser should parse MultiOps_5 qasm");

        assert_eq!(5, parsed.num_qubits);
        assert_eq!(13, parsed.ops.len());

        let cx_count = parsed
            .ops
            .iter()
            .filter(|op| {
                matches!(
                    op,
                    QInstruct::SingleCtrl((SingleCtrlQubitOp::ControlledX, _, _))
                )
            })
            .count();
        let cz_count = parsed
            .ops
            .iter()
            .filter(|op| {
                matches!(
                    op,
                    QInstruct::SingleCtrl((SingleCtrlQubitOp::ControlledZ, _, _))
                )
            })
            .count();

        assert_eq!(2, cx_count);
        assert_eq!(2, cz_count);
    }

    #[test]
    fn skips_unknown_op() {
        let qasm = "qreg q[2];\nfoo q[1];\nx q[0];";

        let parsed = openq3_parser::parse(qasm).expect("parser should skip unknown op");

        assert_eq!(2, parsed.num_qubits);
        assert_eq!(1, parsed.ops.len());
        assert_eq!(QInstruct::Single((QuantumOp::PauliX, 0)), parsed.ops[0]);
    }

    #[test]
    fn invalid_num_qubits() {
        let qasm = "qreg q[];\nx q[0];";

        let parsed = openq3_parser::parse(qasm);
        assert!(parsed.is_err());
    }

    #[test]
    fn invalid_target_qubit() {
        let qasm = "qreg q[2];\nx q[a];";

        let parsed = openq3_parser::parse(qasm);
        assert!(parsed.is_err());
    }
}
