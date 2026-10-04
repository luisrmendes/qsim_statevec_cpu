use oq3_circuit::{parse_circuit_file, GateApplication};
use qsim_statevec_cpu::{QInstruct, QubitLayer, SingleCtrlQubitOp, SingleQubitOp, TwoCtrlQubitOp};
use std::fs;
use std::path::Path;

// Map instructions gathered from file to qsim_statevec_cpu's `QInstruct` representation.
fn to_instruction(gate: &GateApplication) -> QInstruct {
    match (gate.name.as_str(), gate.qubits.as_slice()) {
        ("x", &[t]) => (SingleQubitOp::PauliX, t).into(),
        ("y", &[t]) => (SingleQubitOp::PauliY, t).into(),
        ("z", &[t]) => (SingleQubitOp::PauliZ, t).into(),
        ("h", &[t]) => (SingleQubitOp::Hadamard, t).into(),
        ("s", &[t]) => (SingleQubitOp::S, t).into(),
        ("t", &[t]) => (SingleQubitOp::T, t).into(),
        ("sx", &[t]) => (SingleQubitOp::SX, t).into(),
        ("sy", &[t]) => (SingleQubitOp::SY, t).into(),
        ("cx", &[c, t]) => (SingleCtrlQubitOp::ControlledX, c, t).into(),
        ("cy", &[c, t]) => (SingleCtrlQubitOp::ControlledY, c, t).into(),
        ("cz", &[c, t]) => (SingleCtrlQubitOp::ControlledZ, c, t).into(),
        ("ccx", &[c0, c1, t]) => (TwoCtrlQubitOp::Toffoli, c0, c1, t).into(),
        _ => panic!("unsupported gate: {gate:?}"),
    }
}

// Prints `NAME=P0,P1,...` (probability of each qubit being 1) for every reference circuit.
// Run through functional_tests/run_functional_tests.py, which compares the output against Qiskit.
fn main() {
    let circuits_dir = Path::new("functional_tests/reference_qasm");
    let entries = fs::read_dir(circuits_dir)
        .unwrap_or_else(|error| panic!("failed to read {}: {error}", circuits_dir.display()));

    for path in entries.filter_map(Result::ok).map(|entry| entry.path()) {
        if path.extension().is_none_or(|ext| ext != "openqasm") {
            continue;
        }
        let circuit = parse_circuit_file(&path)
            .unwrap_or_else(|error| panic!("failed to parse {}: {error}", path.display()));
        let instructions: Vec<QInstruct> = circuit.gates.iter().map(to_instruction).collect();

        let mut layer = QubitLayer::new(circuit.num_qubits);
        layer
            .execute_noiseless(&instructions)
            .unwrap_or_else(|error| panic!("failed to execute {}: {error}", path.display()));

        let probabilities: Vec<String> =
            layer.measure_qubits().iter().map(f64::to_string).collect();
        let name = path.file_stem().unwrap_or_default().to_string_lossy();
        println!("{name}={}", probabilities.join(","));
    }
}
