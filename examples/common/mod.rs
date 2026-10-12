use oq3_circuit::{parse_circuit_file, Gate, ParsedCircuit};
use qsim_statevec_cpu::{QSimCircuit, QSimGate, SingleCtrlQubitOp, SingleQubitOp, TwoCtrlQubitOp};
use std::fs;
use std::path::{Path, PathBuf};

// Anchored on the crate root so the examples run from any working directory
const REFERENCE_DIR: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/qasm_files/reference");

// Parses every reference circuit, sorted by file name, as (name, circuit) pairs
pub fn reference_circuits() -> Vec<(String, QSimCircuit)> {
    let circuits_dir = Path::new(REFERENCE_DIR);
    let entries = fs::read_dir(circuits_dir)
        .unwrap_or_else(|error| panic!("failed to read {}: {error}", circuits_dir.display()));

    let mut paths: Vec<PathBuf> = entries
        .filter_map(Result::ok)
        .map(|entry| entry.path())
        .filter(|path| path.extension().is_some_and(|ext| ext == "qasm"))
        .collect();
    paths.sort();

    paths
        .into_iter()
        .map(|path| {
            let circuit: ParsedCircuit = parse_circuit_file(&path)
                .unwrap_or_else(|error| panic!("failed to parse {}: {error}", path.display()));
            let name = path
                .file_stem()
                .unwrap_or_default()
                .to_string_lossy()
                .into_owned();
            (name, convert_parsed_circuit_into_qsim_circuit(circuit))
        })
        .collect()
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
