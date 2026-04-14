from pathlib import Path

from qiskit import QuantumCircuit
from qiskit.quantum_info import Statevector


def format_probability(value: float) -> str:
    rounded = round(value)
    if abs(value - rounded) < 1e-12:
        return str(int(rounded))
    return f"{value:.12f}".rstrip("0").rstrip(".")


def format_measurements(qc: QuantumCircuit, num_qubits: int) -> str:
    probabilities = Statevector.from_instruction(qc).probabilities_dict()
    qubit_probs = [
        sum(float(prob) for basis, prob in probabilities.items() if str(basis)[-1 - qubit] == "1")
        for qubit in range(num_qubits)
    ]
    return ",".join(format_probability(value) for value in qubit_probs)


circuits_dir = Path(__file__).with_name("reference_qasm")
qasm_files = sorted(circuits_dir.glob("*.openqasm"))

if not qasm_files:
    raise FileNotFoundError(f"No QASM files found in {circuits_dir}")

for qasm_file in qasm_files:
    qc = QuantumCircuit.from_qasm_file(str(qasm_file))
    print(f"{qasm_file.stem}={format_measurements(qc, qc.num_qubits)}")
