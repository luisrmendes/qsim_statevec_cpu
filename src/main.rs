use qsim_statevec_cpu::{QuantumOp, QubitLayer};

fn main() {
    let mut layer = QubitLayer::new(5);
    let instructions = vec![
        (QuantumOp::PauliX, 0),
        (QuantumOp::Hadamard, 1),
        (QuantumOp::S, 2),
        (QuantumOp::T, 3),
        (QuantumOp::PauliZ, 4),
        (QuantumOp::SX, 0),
        (QuantumOp::SY, 1),
    ];
    if let Err(error) = layer.execute_noiseless(&instructions) {
        eprintln!("Failed to execute instructions: {error}");
        return;
    }
    // println!("{layer:?}");

    let res = layer.measure_qubits();
    println!("{res:?}");

}
