use oq3_semantics::circuit::parse_circuit_file;
use qsim_statevec_cpu::{
    QInstruct, QInstructs, QuantumOp, QubitLayer, SingleCtrlQubitOp, TwoCtrlQubitOp,
};

fn main() {
    let qasm_path = "qasm_files/openqasm3_simple.qasm";

    let circuit = match parse_circuit_file(qasm_path) {
        Ok(c) => c,
        Err(error) => {
            eprintln!("Failed to parse OpenQASM file '{qasm_path}': {error}");
            return;
        }
    };

    let ops: QInstructs = match circuit
        .gates
        .iter()
        .map(|g| map_gate_to_instruction(&g.name, &g.qubits))
        .collect::<Result<_, _>>()
    {
        Ok(ops) => ops,
        Err(error) => {
            eprintln!("Failed to map gates to instructions: {error}");
            return;
        }
    };

    let mut layer = QubitLayer::new(circuit.num_qubits);
    if let Err(error) = layer.execute_noiseless(&ops) {
        eprintln!("Failed to execute parsed OpenQASM instructions: {error}");
        return;
    }

    println!(
        "Executed '{}' with {} qubits",
        qasm_path, circuit.num_qubits
    );
    println!("Measured probabilities: {:?}", layer.measure_qubits());
}

fn map_gate_to_instruction(gate_name: &str, qubits: &[u32]) -> Result<QInstruct, String> {
    match (gate_name, qubits) {
        ("x", [q]) => Ok(QInstruct::Single((QuantumOp::PauliX, *q))),
        ("y", [q]) => Ok(QInstruct::Single((QuantumOp::PauliY, *q))),
        ("z", [q]) => Ok(QInstruct::Single((QuantumOp::PauliZ, *q))),
        ("h", [q]) => Ok(QInstruct::Single((QuantumOp::Hadamard, *q))),
        ("s", [q]) => Ok(QInstruct::Single((QuantumOp::S, *q))),
        ("t", [q]) => Ok(QInstruct::Single((QuantumOp::T, *q))),
        ("sx", [q]) => Ok(QInstruct::Single((QuantumOp::SX, *q))),
        ("sy", [q]) => Ok(QInstruct::Single((QuantumOp::SY, *q))),
        ("cx", [c, t]) => Ok(QInstruct::SingleCtrl((
            SingleCtrlQubitOp::ControlledX,
            *c,
            *t,
        ))),
        ("cy", [c, t]) => Ok(QInstruct::SingleCtrl((
            SingleCtrlQubitOp::ControlledY,
            *c,
            *t,
        ))),
        ("cz", [c, t]) => Ok(QInstruct::SingleCtrl((
            SingleCtrlQubitOp::ControlledZ,
            *c,
            *t,
        ))),
        ("ccx", [c1, c2, t]) => Ok(QInstruct::TwoCtrl((TwoCtrlQubitOp::Toffoli, *c1, *c2, *t))),
        _ => Err(format!(
            "Unsupported gate or arity: {gate_name} with {} operands",
            qubits.len()
        )),
    }
}
