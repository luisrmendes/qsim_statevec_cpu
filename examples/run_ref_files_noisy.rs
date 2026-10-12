mod common;

use qsim_statevec_cpu::{execute_noisy_shots, NoiseModel};

// TODO: Compare output to Qiskit's
fn main() {
    for (name, circuit) in common::reference_circuits() {
        let results = execute_noisy_shots(
            circuit,
            100,
            NoiseModel {
                gate_error_prob: 0.01,
                readout_flip_prob: 0.01,
            },
        )
        .unwrap_or_else(|error| panic!("failed to execute {name}: {error}"));

        let measured_qubits: Vec<String> = results.iter().map(f64::to_string).collect();
        println!("{name}={}", measured_qubits.join(","));
    }
}
