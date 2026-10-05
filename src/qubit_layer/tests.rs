use super::*;

use num::pow;
use num::Complex;

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

        let result = q_layer.execute_instructions(vec![
            QSimGate::Single {
                op: SingleQubitOp::Hadamard,
                target: 0,
            },
            QSimGate::Single {
                op: SingleQubitOp::S,
                target: 0,
            },
        ]);

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

        let result = q_layer.execute_instructions(vec![
            QSimGate::Single {
                op: SingleQubitOp::Hadamard,
                target: 0,
            },
            QSimGate::Single {
                op: SingleQubitOp::T,
                target: 0,
            },
        ]);
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

        let result = q_layer.execute_instructions(vec![QSimGate::Single {
            op: SingleQubitOp::SX,
            target: 0,
        }]);
        assert!(result.is_ok());

        let expected = vec![Complex::new(0.5, 0.5), Complex::new(0.5, -0.5)];
        assert_eq!(expected, q_layer.main);
    }

    #[test]
    fn test_sqrt_pauli_x_squared_equals_pauli_x() {
        let mut q_layer: QubitLayer = QubitLayer::new(1);

        let result = q_layer.execute_instructions(vec![
            QSimGate::Single {
                op: SingleQubitOp::SX,
                target: 0,
            },
            QSimGate::Single {
                op: SingleQubitOp::SX,
                target: 0,
            },
        ]);
        assert!(result.is_ok());

        let results = q_layer.measure_qubits();
        assert_eq!(results[0], 1.0);
    }

    #[test]
    fn test_sqrt_pauli_y_simple() {
        let mut q_layer: QubitLayer = QubitLayer::new(1);

        let result = q_layer.execute_instructions(vec![QSimGate::Single {
            op: SingleQubitOp::SY,
            target: 0,
        }]);
        assert!(result.is_ok());

        let expected = vec![Complex::new(0.5, 0.5), Complex::new(0.5, 0.5)];
        assert_eq!(expected, q_layer.main);
    }

    #[test]
    fn test_sqrt_pauli_y_squared_equals_pauli_y() {
        let mut q_layer: QubitLayer = QubitLayer::new(1);

        let result = q_layer.execute_instructions(vec![
            QSimGate::Single {
                op: SingleQubitOp::SY,
                target: 0,
            },
            QSimGate::Single {
                op: SingleQubitOp::SY,
                target: 0,
            },
        ]);
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

        let prep = vec![QSimGate::Single {
            op: SingleQubitOp::PauliX,
            target: 1,
        }];
        let prep_result = q_layer.execute_instructions(prep);
        assert!(prep_result.is_ok());

        let cy_result = q_layer.execute_instructions(vec![QSimGate::SingleCtrl {
            op: SingleCtrlQubitOp::ControlledY,
            control: 1,
            target: 0,
        }]);
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

        let prep = vec![
            QSimGate::Single {
                op: SingleQubitOp::PauliX,
                target: 1,
            },
            QSimGate::Single {
                op: SingleQubitOp::PauliX,
                target: 2,
            },
        ];
        let prep_result = q_layer.execute_instructions(prep);
        assert!(prep_result.is_ok());

        let toffoli_result = q_layer.execute_instructions(vec![QSimGate::TwoCtrl {
            op: TwoCtrlQubitOp::Toffoli,
            controls: [1, 2],
            target: 0,
        }]);
        assert!(toffoli_result.is_ok());

        let measured = q_layer.measure_qubits();
        assert_eq!(1.0, measured[0]);
        assert_eq!(1.0, measured[1]);
        assert_eq!(1.0, measured[2]);
    }

    #[test]
    fn test_execute_noiseless_toffoli_out_of_range() {
        let mut q_layer: QubitLayer = QubitLayer::new(3);

        let result = q_layer.execute_instructions(vec![QSimGate::TwoCtrl {
            op: TwoCtrlQubitOp::Toffoli,
            controls: [0, 7],
            target: 2,
        }]);
        assert!(result.is_err());
        assert_eq!(
            result.unwrap_err(),
            "Control qubit 7 is out of range. Size of layer is 3"
        );
    }

    #[test]
    fn test_execute_toffoli_with_same_control_qubits() {
        let mut q_layer: QubitLayer = QubitLayer::new(2);

        let result = q_layer.execute_instructions(vec![QSimGate::TwoCtrl {
            op: TwoCtrlQubitOp::Toffoli,
            controls: [0, 1],
            target: 1,
        }]);
        assert!(result.is_err());
        assert_eq!(
            result.unwrap_err(),
            "Target qubit and control qubit 2 are the same: 1"
        );
    }
}
