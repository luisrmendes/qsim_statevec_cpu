import json
from pathlib import Path

from qiskit import QuantumCircuit
from qiskit.circuit.library import YGate
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


def apply_instruction(qc: QuantumCircuit, instruction: dict) -> None:
    op = instruction["op"]
    qubits = instruction["qubits"]

    if op == "x" and len(qubits) == 1:
        qc.x(qubits[0])
    elif op == "z" and len(qubits) == 1:
        qc.z(qubits[0])
    elif op == "h" and len(qubits) == 1:
        qc.h(qubits[0])
    elif op == "s" and len(qubits) == 1:
        qc.s(qubits[0])
    elif op == "t" and len(qubits) == 1:
        qc.t(qubits[0])
    elif op == "sx" and len(qubits) == 1:
        qc.sx(qubits[0])
    elif op == "sy" and len(qubits) == 1:
        qc.append(YGate().power(0.5), [qubits[0]])
    elif op == "cx" and len(qubits) == 2:
        qc.cx(qubits[0], qubits[1])
    elif op == "cz" and len(qubits) == 2:
        qc.cz(qubits[0], qubits[1])
    elif op == "ccx" and len(qubits) == 3:
        qc.ccx(qubits[0], qubits[1], qubits[2])
    else:
        raise ValueError(f"Unsupported instruction: {instruction}")


circuits_file = Path(__file__).with_name("reference_circuits.json")
with circuits_file.open("r", encoding="utf-8") as handle:
    circuits_data = json.load(handle)

for case in circuits_data["cases"]:
    qc = QuantumCircuit(case["num_qubits"])
    for instruction in case["instructions"]:
        apply_instruction(qc, instruction)
    print(f"{case['name']}={format_measurements(qc, case['num_qubits'])}")
