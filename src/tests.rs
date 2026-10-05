mod qsim_statevec_cpu_tests {
    use crate::*;

    #[test]
    fn test_execute_shots() {
        let instructions = vec![
            QSimGate::Single {
                op: SingleQubitOp::Hadamard,
                target: 0,
            },
            QSimGate::Single {
                op: SingleQubitOp::PauliX,
                target: 1,
            },
        ];

        let mut accumulated_layer = QubitLayer::new(3);
        for _ in 0..4 {
            let mut shot_layer = QubitLayer::new(3);
            let result = shot_layer.execute_instructions(instructions.clone());
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
    fn test_execute_zero_shots_returns_error() {
        let circuit = QSimCircuit {
            num_qubits: 3,
            gates: vec![QSimGate::Single {
                op: SingleQubitOp::Hadamard,
                target: 0,
            }],
        };
        let noise = NoiseModel {
            gate_error_prob: 0.0,
            readout_flip_prob: 0.0,
        };

        let result = execute_noisy_shots(circuit, 0, noise);
        assert!(result.is_err());
    }

    #[test]
    fn test_execute_shots_failed_execute() {
        let mut q_layer = QubitLayer::new(3);
        let instructions = vec![QSimGate::Single {
            op: SingleQubitOp::PauliX,
            target: 10,
        }];

        let result = q_layer.execute_instructions(instructions);
        assert!(result.is_err());
    }

    #[test]
    fn test_execute_noisy_shots_zero_noise() {
        let circuit = QSimCircuit {
            num_qubits: 3,
            gates: vec![
                QSimGate::Single {
                    op: SingleQubitOp::Hadamard,
                    target: 0,
                },
                QSimGate::Single {
                    op: SingleQubitOp::PauliX,
                    target: 1,
                },
            ],
        };
        let noise = NoiseModel {
            gate_error_prob: 0.0,
            readout_flip_prob: 0.0,
        };

        let measured = execute_noisy_shots(circuit, 3, noise).expect("shots should execute");
        assert_eq!(0.5, (measured[0] * 10.0).round() / 10.0);
        assert_eq!(1.0, (measured[1] * 10.0).round() / 10.0);
        assert_eq!(0.0, (measured[2] * 10.0).round() / 10.0);
    }

    #[test]
    fn test_execute_noisy_shots_readout_flip_full() {
        let circuit = QSimCircuit {
            num_qubits: 3,
            gates: vec![],
        };
        let noise = NoiseModel {
            gate_error_prob: 0.0,
            readout_flip_prob: 1.0,
        };

        let measured = execute_noisy_shots(circuit, 2, noise).expect("shots should execute");

        // Full flip probability turns 0% probabilities into 100% probabilities
        assert_eq!(1.0, measured[0]);
        assert_eq!(1.0, measured[1]);
        assert_eq!(1.0, measured[2]);
    }

    #[test]
    fn test_execute_noisy_shots_invalid_noise() {
        let circuit = QSimCircuit {
            num_qubits: 1,
            gates: vec![QSimGate::Single {
                op: SingleQubitOp::PauliX,
                target: 0,
            }],
        };
        let noise = NoiseModel {
            gate_error_prob: 1.1,
            readout_flip_prob: 0.0,
        };

        let result = execute_noisy_shots(circuit, 1, noise);
        assert!(result.is_err());
    }

    #[test]
    fn test_random_executions() {
        let mut q_layer: QubitLayer = QubitLayer::new(3);
        let instructions = vec![
            QSimGate::Single {
                op: SingleQubitOp::Hadamard,
                target: 0,
            },
            QSimGate::Single {
                op: SingleQubitOp::Hadamard,
                target: 1,
            },
            QSimGate::Single {
                op: SingleQubitOp::Hadamard,
                target: 2,
            },
            QSimGate::Single {
                op: SingleQubitOp::Hadamard,
                target: 0,
            },
            QSimGate::Single {
                op: SingleQubitOp::Hadamard,
                target: 1,
            },
            QSimGate::Single {
                op: SingleQubitOp::Hadamard,
                target: 2,
            },
        ];
        let _ = q_layer.execute_instructions(instructions);
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

        let instructions: Vec<QSimGate> = vec![];
        let _ = q_layer.execute_instructions(instructions);

        for it in 0..q_layer.get_num_qubits() {
            assert_eq!(0.0, q_layer.measure_qubits()[it as usize]);
        }
    }

    #[test]
    fn test_spins_on_superposition() {
        let mut q_layer: QubitLayer = QubitLayer::new(3);
        let instructions = vec![
            QSimGate::Single {
                op: SingleQubitOp::Hadamard,
                target: 0,
            },
            QSimGate::Single {
                op: SingleQubitOp::Hadamard,
                target: 1,
            },
            QSimGate::Single {
                op: SingleQubitOp::Hadamard,
                target: 2,
            },
            QSimGate::Single {
                op: SingleQubitOp::PauliX,
                target: 0,
            },
            QSimGate::Single {
                op: SingleQubitOp::PauliY,
                target: 1,
            },
            QSimGate::Single {
                op: SingleQubitOp::PauliZ,
                target: 2,
            },
        ];
        let _ = q_layer.execute_instructions(instructions);
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
        let instructions = vec![QSimGate::Single {
            op: SingleQubitOp::PauliX,
            target: 10,
        }]; // index goes up to 9

        let result: Result<(), String> = q_layer.execute_instructions(instructions);
        assert!(result.is_err());

        let result: Result<(), String> = q_layer.execute_instructions(vec![QSimGate::Single {
            op: SingleQubitOp::Hadamard,
            target: 2112,
        }]);
        assert!(result.is_err());
    }

    #[test]
    fn test_execute() {
        let mut q_layer: QubitLayer = QubitLayer::new(10);
        let instructions = vec![
            QSimGate::Single {
                op: SingleQubitOp::PauliX,
                target: 0,
            },
            QSimGate::Single {
                op: SingleQubitOp::PauliY,
                target: 1,
            },
            QSimGate::Single {
                op: SingleQubitOp::PauliZ,
                target: 2,
            },
        ];

        if let Err(e) = q_layer.execute_instructions(instructions) {
            panic!("Should not panic!. Error: {e}");
        }

        assert_eq!(1.0, q_layer.measure_qubits()[0].round());
    }
}
