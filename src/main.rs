use qsim_statevec_cpu::{openq3_parser, execute_shots, execute_shots_noisy, NoiseModel};
use std::fs;

fn main() {
    let file_path = "qasm_files/Grovers_5_qubits.openqasm";
    let shots = 1;

    let file_contents = match fs::read_to_string(file_path) {
        Ok(contents) => contents,
        Err(error) => {
            eprintln!("Failed to read QASM file '{}': {error}", file_path);
            return;
        }
    };

    let parsed = match openq3_parser::parse(&file_contents) {
        Ok(parsed) => parsed,
        Err(error) => {
            eprintln!("Failed to parse QASM file '{}': {error}", file_path);
            return;
        }
    };

    let layer = match execute_shots(parsed.ops.clone(), parsed.num_qubits, shots) {
        Ok(layer) => layer,
        Err(error) => {
            eprintln!("Failed to execute parsed instructions: {error}");
            return;
        }
    };
    println!("Average measured qubits: {:#?}", layer.measure_qubits());

    let noise_model = NoiseModel {
        gate_error_prob: 0.02,
        readout_flip_prob: 0.01,
    };

    let noisy_layer = match execute_shots_noisy(parsed.ops, parsed.num_qubits, shots, noise_model)
    {
        Ok(layer) => layer,
        Err(error) => {
            eprintln!("Failed to execute parsed instructions with noise: {error}");
            return;
        }
    };
    println!("Average noisy measured qubits: {:#?}", noisy_layer.measure_qubits());
}
