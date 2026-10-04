use oq3_circuit::{parse_circuit_file, GateApplication};
use qsim_statevec_cpu::{QInstruct, QubitLayer, SingleCtrlQubitOp, SingleQubitOp, TwoCtrlQubitOp};
use std::fs;
use std::path::Path;

fn format_probability(value: f64) -> String {
    let rounded_int = value.round();
    if (value - rounded_int).abs() < 1e-12 {
        return (rounded_int as i64).to_string();
    }

    let rounded = (value * 1_000_000_000_000.0).round() / 1_000_000_000_000.0;
    let mut output = format!("{rounded:.12}");
    while output.contains('.') && output.ends_with('0') {
        output.pop();
    }
    if output.ends_with('.') {
        output.pop();
    }
    output
}

fn format_measurements(measured: &[f64]) -> String {
    measured
        .iter()
        .map(|&value| format_probability(value))
        .collect::<Vec<_>>()
        .join(",")
}

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

fn main() {
    let circuits_dir = Path::new("functional_tests/reference_qasm");
    let mut qasm_files = fs::read_dir(circuits_dir)
        .unwrap_or_else(|error| panic!("failed to read {}: {error}", circuits_dir.display()))
        .filter_map(Result::ok)
        .map(|entry| entry.path())
        .filter(|path| path.extension().is_some_and(|ext| ext == "openqasm"))
        .collect::<Vec<_>>();

    qasm_files.sort();

    assert!(
        !qasm_files.is_empty(),
        "no .openqasm files found in {}",
        circuits_dir.display()
    );

    for qasm_file in qasm_files {
        let circuit = parse_circuit_file(&qasm_file)
            .unwrap_or_else(|error| panic!("failed to parse {}: {error}", qasm_file.display()));

        let mut layer = QubitLayer::new(circuit.num_qubits);

        let instructions: Vec<QInstruct> = circuit.gates.iter().map(to_instruction).collect();

        layer
            .execute_noiseless(&instructions)
            .unwrap_or_else(|error| panic!("{qasm_file:?} should execute without errors: {error}"));

        let measured = layer.measure_qubits();
        let name = qasm_file.file_stem().unwrap_or_default().to_string_lossy();
        println!("{name}={}", format_measurements(&measured));
    }
}
