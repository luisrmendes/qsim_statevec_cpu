mod common;

use qsim_statevec_cpu::execute_noiseless;

// Prints `NAME=P0,P1,...` (probability of each qubit being 1) for every reference circuit.
// Run through functional_tests/run_functional_tests.py, which compares the output against Qiskit.
fn main() {
    for (name, circuit) in common::reference_circuits() {
        let results = execute_noiseless(circuit)
            .unwrap_or_else(|error| panic!("failed to execute {name}: {error}"));

        let measured_qubits: Vec<String> = results.iter().map(f64::to_string).collect();
        println!("{name}={}", measured_qubits.join(","));
    }
}
