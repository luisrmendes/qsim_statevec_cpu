use qsim_statevec_cpu::{QuantumOp, QubitLayer};

fn main() {
    // let file_path = "qasm_files/Grovers_5_qubits.openqasm";
    // let shots = 1;

    // let file_contents = match fs::read_to_string(file_path) {
    //     Ok(contents) => contents,
    //     Err(error) => {
    //         eprintln!("Failed to read QASM file '{file_path}': {error}");
    //         return;
    //     }
    // };

    // let parsed = match openq3_parser::parse(&file_contents) {
    //     Ok(parsed) => parsed,
    //     Err(error) => {
    //         eprintln!("Failed to parse QASM file '{file_path}': {error}");
    //         return;
    //     }
    // };

    // let noise_model = NoiseModel {
    //     gate_error_prob: 0.02,
    //     readout_flip_prob: 0.01,
    // };

    // let mut noisy_layer = QubitLayer::new(parsed.num_qubits);
    // let result = noisy_layer.execute_noisy_shots(&parsed.ops, shots, noise_model);
    // if let Err(error) = result {
    //     eprintln!("Failed to execute parsed instructions with noise: {error}");
    //     return;
    // }
    // println!(
    //     "Average noisy measured qubits: {:#?}",
    //     noisy_layer.measure_qubits()
    // );

    let mut layer = QubitLayer::new(2);
    let instructions = vec![
        (QuantumOp::Hadamard, 0),
        (QuantumOp::PauliY, 0),
        (QuantumOp::Hadamard, 1),
        // (QuantumOp::S, 1),
        (QuantumOp::T, 0),
        (QuantumOp::T, 1),
    ];
    if let Err(error) = layer.execute_noiseless(&instructions) {
        eprintln!("Failed to execute instructions: {error}");
        return;
    }
    println!("{layer:#?}");
}
