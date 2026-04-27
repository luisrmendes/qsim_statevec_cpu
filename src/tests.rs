use super::*;

// More extensive functionality tests on quantum gate operations
mod openqasm_tests {
    use super::*;

    fn parse_qasm_file_to_ops(qasm_path: &str) -> (u32, QInstructs) {
        let circuit = oq3_semantics::circuit::parse_circuit_file(qasm_path)
            .expect("parser should parse OpenQASM file");

        let ops: QInstructs = circuit
            .gates
            .iter()
            .map(|g| map_gate_to_instruction(&g.name, &g.qubits))
            .collect::<Result<_, _>>()
            .expect("gate mapping should succeed");

        (circuit.num_qubits, ops)
    }

    fn map_gate_to_instruction(gate_name: &str, qubits: &[u32]) -> Result<QInstruct, String> {
        match (gate_name, qubits) {
            ("x", [q]) => Ok(QInstruct::Single((QuantumOp::PauliX, *q))),
            ("y", [q]) => Ok(QInstruct::Single((QuantumOp::PauliY, *q))),
            ("z", [q]) => Ok(QInstruct::Single((QuantumOp::PauliZ, *q))),
            ("h", [q]) => Ok(QInstruct::Single((QuantumOp::Hadamard, *q))),
            ("s", [q]) => Ok(QInstruct::Single((QuantumOp::S, *q))),
            ("t", [q]) => Ok(QInstruct::Single((QuantumOp::T, *q))),
            ("sx", [q]) => Ok(QInstruct::Single((QuantumOp::SX, *q))),
            ("sy", [q]) => Ok(QInstruct::Single((QuantumOp::SY, *q))),
            ("cx", [c, t]) => Ok(QInstruct::SingleCtrl((
                SingleCtrlQubitOp::ControlledX,
                *c,
                *t,
            ))),
            ("cy", [c, t]) => Ok(QInstruct::SingleCtrl((
                SingleCtrlQubitOp::ControlledY,
                *c,
                *t,
            ))),
            ("cz", [c, t]) => Ok(QInstruct::SingleCtrl((
                SingleCtrlQubitOp::ControlledZ,
                *c,
                *t,
            ))),
            ("ccx", [c1, c2, t]) => Ok(QInstruct::TwoCtrl((TwoCtrlQubitOp::Toffoli, *c1, *c2, *t))),
            _ => Err(format!(
                "Unsupported gate or arity: {gate_name} with {} operands",
                qubits.len()
            )),
        }
    }

    #[test]
    fn misc1_3_qubits() {
        let (num_qubits, ops) = parse_qasm_file_to_ops("qasm_files/misc1_3_qubits.qasm");
        let mut q_layer = QubitLayer::new(num_qubits);

        let result = q_layer.execute_noiseless(&ops);
        assert!(result.is_ok());

        let measured = q_layer.measure_qubits();
        assert_eq!(round_to(measured[0], 2), 0.5);
        assert_eq!(round_to(measured[1], 2), 1.0);
        assert_eq!(round_to(measured[2], 2), 0.5);
    }

    #[test]
    fn misc1_4_qubits() {
        let (num_qubits, ops) = parse_qasm_file_to_ops("qasm_files/misc1_4_qubits.qasm");
        let mut q_layer = QubitLayer::new(num_qubits);

        let result = q_layer.execute_noiseless(&ops);
        assert!(result.is_ok());

        let measured = q_layer.measure_qubits();
        assert_eq!(round_to(measured[0], 2), 0.5);
        assert_eq!(round_to(measured[1], 2), 0.5);
        assert_eq!(round_to(measured[2], 2), 1.0);
        assert_eq!(round_to(measured[3], 2), 0.5);
    }

    #[test]
    fn misc1_5_qubits() {
        let (num_qubits, ops) = parse_qasm_file_to_ops("qasm_files/misc1_5_qubits.qasm");
        let mut q_layer = QubitLayer::new(num_qubits);

        let result = q_layer.execute_noiseless(&ops);
        assert!(result.is_ok());

        let measured = q_layer.measure_qubits();
        assert_eq!(round_to(measured[0], 2), 0.5);
        assert_eq!(round_to(measured[1], 2), 0.5);
        assert_eq!(round_to(measured[2], 2), 0.5);
        assert_eq!(round_to(measured[3], 2), 1.0);
        assert_eq!(round_to(measured[4], 2), 0.5);
    }

    #[test]
    fn misc2_3_qubits() {
        let (num_qubits, ops) = parse_qasm_file_to_ops("qasm_files/misc2_3_qubits.qasm");
        let mut q_layer = QubitLayer::new(num_qubits);

        let result = q_layer.execute_noiseless(&ops);
        assert!(result.is_ok());

        let measured = q_layer.measure_qubits();
        assert_eq!(measured[0], 1.0);
        assert_eq!(measured[1], 0.0);
        assert_eq!(measured[2], 0.0);
    }

    #[test]
    fn ctrl_x_1_5_qubits() {
        let (num_qubits, ops) = parse_qasm_file_to_ops("qasm_files/ctrl_x_1_5_qubits.qasm");
        let mut q_layer = QubitLayer::new(num_qubits);

        let result = q_layer.execute_noiseless(&ops);
        assert!(result.is_ok());

        let measured = q_layer.measure_qubits();
        assert_eq!(num_qubits as usize, measured.len());
        for value in measured {
            assert_eq!(0.5, round_to(value, 2));
        }
    }

    #[test]
    fn ctrl_z_1_5_qubits() {
        let (num_qubits, ops) = parse_qasm_file_to_ops("qasm_files/ctrl_z_1_5_qubits.qasm");
        let mut q_layer = QubitLayer::new(num_qubits);

        let result = q_layer.execute_noiseless(&ops);
        assert!(result.is_ok());

        let measured = q_layer.measure_qubits();
        assert_eq!(num_qubits as usize, measured.len());
        for value in measured {
            assert_eq!(0.5, (value * 10.0).round() / 10.0);
        }
    }

    fn round_to(x: f64, places: u32) -> f64 {
        let factor = 10_f64.powi(places as i32);
        (x * factor).round() / factor
    }
}

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
            let result = shot_layer.execute_noiseless(&instructions);
            assert!(result.is_ok());
            accumulated_layer += &shot_layer;
        }
        accumulated_layer /= 4;

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
        assert!(result.is_err());

        let measured = q_layer.measure_qubits();
        assert_eq!(0.0, measured[0]);
        assert_eq!(0.0, measured[1]);
        assert_eq!(0.0, measured[2]);
    }

    #[test]
    fn test_execute_shots_failed_execute() {
        let mut q_layer = QubitLayer::new(3);
        let instructions = vec![(QuantumOp::PauliX, 10)];

        let result = q_layer.execute_noiseless(&instructions);
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
        let _ = q_layer.execute_noiseless(&instructions);
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
        let _ = q_layer.execute_noiseless(&instructions);

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
        let _ = q_layer.execute_noiseless(&instructions);
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

        let result: Result<(), String> = q_layer.execute_noiseless(&instructions);
        assert!(result.is_err());

        let result: Result<(), String> = q_layer.execute_noiseless(&[(QuantumOp::Hadamard, 2112)]);
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

        if let Err(e) = q_layer.execute_noiseless(&instructions) {
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

        let result = q_layer.execute_noiseless(&[(QuantumOp::Hadamard, 0), (QuantumOp::S, 0)]);
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

        let result = q_layer.execute_noiseless(&[(QuantumOp::Hadamard, 0), (QuantumOp::T, 0)]);
        assert!(result.is_ok());

        let expected = vec![
            Complex::new(hadamard_const, 0.0),
            Complex::from_polar(hadamard_const, std::f64::consts::FRAC_PI_4),
        ];
        assert_eq!(expected, q_layer.main);
    }

    #[test]
    fn test_sqrt_pauli_x_simple() {
        let mut q_layer: QubitLayer = QubitLayer::new(1);

        let result = q_layer.execute_noiseless(&[(QuantumOp::SX, 0)]);
        assert!(result.is_ok());

        let expected = vec![Complex::new(0.5, 0.5), Complex::new(0.5, -0.5)];
        assert_eq!(expected, q_layer.main);
    }

    #[test]
    fn test_sqrt_pauli_x_squared_equals_pauli_x() {
        let mut q_layer: QubitLayer = QubitLayer::new(1);

        let result = q_layer.execute_noiseless(&[(QuantumOp::SX, 0), (QuantumOp::SX, 0)]);
        assert!(result.is_ok());

        let results = q_layer.measure_qubits();
        assert_eq!(results[0], 1.0);
    }

    #[test]
    fn test_sqrt_pauli_y_simple() {
        let mut q_layer: QubitLayer = QubitLayer::new(1);

        let result = q_layer.execute_noiseless(&[(QuantumOp::SY, 0)]);
        assert!(result.is_ok());

        let expected = vec![Complex::new(0.5, 0.5), Complex::new(0.5, 0.5)];
        assert_eq!(expected, q_layer.main);
    }

    #[test]
    fn test_sqrt_pauli_y_squared_equals_pauli_y() {
        let mut q_layer: QubitLayer = QubitLayer::new(1);

        let result = q_layer.execute_noiseless(&[(QuantumOp::SY, 0), (QuantumOp::SY, 0)]);
        assert!(result.is_ok());

        let expected = vec![Complex::new(0.0, 0.0), Complex::new(0.0, 1.0)];
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
    fn test_controlled_y_simple() {
        let num_qubits = 2;
        let mut q_layer: QubitLayer = QubitLayer::new(num_qubits);

        q_layer.pauli_x(1);
        q_layer.controlled_y(1, 0);

        let expected = vec![
            Complex::new(0.0, 0.0),
            Complex::new(0.0, 0.0),
            Complex::new(0.0, 0.0),
            Complex::new(0.0, 1.0),
        ];
        assert_eq!(expected, q_layer.main);
    }

    #[test]
    fn test_controlled_y_no_change_when_control_inactive() {
        let num_qubits = 2;
        let mut q_layer: QubitLayer = QubitLayer::new(num_qubits);

        q_layer.pauli_x(0);
        q_layer.controlled_y(1, 0);

        let expected = vec![
            Complex::new(0.0, 0.0),
            Complex::new(1.0, 0.0),
            Complex::new(0.0, 0.0),
            Complex::new(0.0, 0.0),
        ];
        assert_eq!(expected, q_layer.main);
    }

    #[test]
    fn test_execute_noiseless_controlled_y_instruction() {
        let mut q_layer: QubitLayer = QubitLayer::new(2);

        let prep = vec![(QuantumOp::PauliX, 1)];
        let prep_result = q_layer.execute_noiseless(&prep);
        assert!(prep_result.is_ok());

        let cy_result = q_layer.execute_noiseless(&[(SingleCtrlQubitOp::ControlledY, 1, 0)]);
        assert!(cy_result.is_ok());

        let expected = vec![
            Complex::new(0.0, 0.0),
            Complex::new(0.0, 0.0),
            Complex::new(0.0, 0.0),
            Complex::new(0.0, 1.0),
        ];
        assert_eq!(expected, q_layer.main);
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

        let prep = vec![(QuantumOp::PauliX, 1), (QuantumOp::PauliX, 2)];
        let prep_result = q_layer.execute_noiseless(&prep);
        assert!(prep_result.is_ok());

        let toffoli_result = q_layer.execute_noiseless(&[(TwoCtrlQubitOp::Toffoli, 1, 2, 0)]);
        assert!(toffoli_result.is_ok());

        let measured = q_layer.measure_qubits();
        assert_eq!(1.0, measured[0]);
        assert_eq!(1.0, measured[1]);
        assert_eq!(1.0, measured[2]);
    }

    #[test]
    fn test_execute_noiseless_toffoli_out_of_range() {
        let mut q_layer: QubitLayer = QubitLayer::new(3);

        let result = q_layer.execute_noiseless(&[(TwoCtrlQubitOp::Toffoli, 0, 7, 2)]);
        assert!(result.is_err());
    }
}
