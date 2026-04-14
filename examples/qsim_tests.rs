use qsim_statevec_cpu::{QInstruct, QuantumOp, QubitLayer, SingleCtrlQubitOp, TwoCtrlQubitOp};
use serde::Deserialize;

#[derive(Deserialize)]
struct CircuitSuite {
    cases: Vec<CircuitCase>,
}

#[derive(Deserialize)]
struct CircuitCase {
    name: String,
    num_qubits: u32,
    instructions: Vec<Instruction>,
}

#[derive(Deserialize)]
struct Instruction {
    op: String,
    qubits: Vec<u32>,
}

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

fn run_case(name: &str, num_qubits: u32, instructions: Vec<QInstruct>) {
    let mut layer = QubitLayer::new(num_qubits);
    layer
        .execute_noiseless(&instructions)
        .unwrap_or_else(|error| panic!("{name} should execute without errors: {error}"));

    let measured = layer.measure_qubits();
    println!("{name}={}", format_measurements(&measured));
}

fn build_instruction(instruction: Instruction) -> Result<QInstruct, String> {
    match (instruction.op.as_str(), instruction.qubits.as_slice()) {
        ("x", [target]) => Ok((QuantumOp::PauliX, *target).into()),
        ("z", [target]) => Ok((QuantumOp::PauliZ, *target).into()),
        ("h", [target]) => Ok((QuantumOp::Hadamard, *target).into()),
        ("s", [target]) => Ok((QuantumOp::S, *target).into()),
        ("t", [target]) => Ok((QuantumOp::T, *target).into()),
        ("sx", [target]) => Ok((QuantumOp::SX, *target).into()),
        ("sy", [target]) => Ok((QuantumOp::SY, *target).into()),
        ("cx", [control, target]) => Ok((SingleCtrlQubitOp::ControlledX, *control, *target).into()),
        ("cz", [control, target]) => Ok((SingleCtrlQubitOp::ControlledZ, *control, *target).into()),
        ("ccx", [control1, control2, target]) => {
            Ok((TwoCtrlQubitOp::Toffoli, *control1, *control2, *target).into())
        }
        _ => Err(format!(
            "unsupported instruction op='{}' qubits={:?}",
            instruction.op, instruction.qubits
        )),
    }
}

fn main() {
    let suite: CircuitSuite =
        serde_json::from_str(include_str!("../functional_tests/reference_circuits.json"))
            .expect("reference_circuits.json must be valid JSON");

    for case in suite.cases {
        let instructions = case
            .instructions
            .into_iter()
            .map(build_instruction)
            .collect::<Result<Vec<_>, _>>()
            .unwrap_or_else(|error| panic!("{} has invalid instructions: {error}", case.name));

        run_case(&case.name, case.num_qubits, instructions);
    }
}
