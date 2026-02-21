use qsim_statevec_cpu::QuantumOp;

fn main() {
    let instructions = vec![
        (QuantumOp::Hadamard, 0),
        (QuantumOp::PauliX, 1),
        (QuantumOp::PauliZ, 2),
    ];

    let res = qsim_statevec_cpu::execute_shots(instructions, 3, 1000).unwrap();

    let measured_qubits = res.measure_qubits();

    println!("Measured qubits: {measured_qubits:#?}");
    println!("{res:?}");
}
