use super::*;

#[test]
fn test_measure_qubits() {
    let mut q_layer: QubitLayer = QubitLayer::new(3);
    for it in 0..q_layer.get_num_qubits() {
        assert_eq!(0.0, q_layer.measure_qubits()[it as usize]);
    }

    let _ = q_layer.execute(vec![]);

    for it in 0..q_layer.get_num_qubits() {
        assert_eq!(0.0, q_layer.measure_qubits()[it as usize]);
    }
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
    let _ = q_layer.execute(instructions);
    for it in 0..q_layer.get_num_qubits() {
        assert_eq!(
            0.0,
            (q_layer.measure_qubits()[it as usize] * 10.0).round() / 10.0
        );
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
    let _ = q_layer.execute(instructions);
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

    let result: Result<(), String> = q_layer.execute(instructions);
    assert!(result.is_err());

    let result: Result<(), String> = q_layer.execute(vec![(QuantumOp::Hadamard, 2112)]);
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

    if let Err(e) = q_layer.execute(instructions) {
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
fn test_display_trait_print() {
    let q_layer: QubitLayer = QubitLayer::new(2);
    let expected = "1+0i 0+0i 0+0i 0+0i";
    println!("{}", q_layer);
    assert_eq!(expected, format!("{}", q_layer));
}

#[test]
fn test_debug_trait_print() {
    let q_layer: QubitLayer = QubitLayer::new(2);
    let expected = "state 0 -> 1+0i\nstate 1 -> 0+0i\nstate 10 -> 0+0i\nstate 11 -> 0+0i\n";
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
