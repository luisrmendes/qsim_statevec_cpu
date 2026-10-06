use oq3_circuit::{parse_circuit_file, Gate, ParsedCircuit};
use qsim_statevec_cpu::{
    execute_noisy_shots, QSimCircuit, QSimGate, SingleCtrlQubitOp, SingleQubitOp, TwoCtrlQubitOp,
};
use std::fs;
use std::path::Path;

// TODO: Compare output to Qiskit's
fn main() {
    let circuits_dir = Path::new("functional_tests/reference_qasm");
    let entries = fs::read_dir(circuits_dir)
        .unwrap_or_else(|error| panic!("failed to read {}: {error}", circuits_dir.display()));

    for path in entries.filter_map(Result::ok).map(|entry| entry.path()) {
        if path.extension().is_none_or(|ext| ext != "openqasm") {
            continue;
        }
        let circuit: ParsedCircuit = parse_circuit_file(&path)
            .unwrap_or_else(|error| panic!("failed to parse {}: {error}", path.display()));

        let circuit: QSimCircuit = convert_parsed_circuit_into_qsim_circuit(circuit);

        let results = execute_noisy_shots(
            circuit,
            100,
            qsim_statevec_cpu::NoiseModel {
                gate_error_prob: 0.01,
                readout_flip_prob: 0.01,
            },
        )
        .unwrap_or_else(|error| panic!("failed to execute {}: {error}", path.display()));

        let measured_qubits: Vec<String> = results.iter().map(f64::to_string).collect();

        let name = path.file_stem().unwrap_or_default().to_string_lossy();
        println!("{name}={}", measured_qubits.join(","));
    }
}

// Convert ParsedCircuit into QSimCircuit
fn convert_parsed_circuit_into_qsim_circuit(circuit: ParsedCircuit) -> QSimCircuit {
    fn to_qsim_gate(gate: Gate) -> QSimGate {
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
    let gates = circuit.gates.into_iter().map(to_qsim_gate).collect();

    QSimCircuit {
        num_qubits: circuit.num_qubits,
        gates,
    }
}
