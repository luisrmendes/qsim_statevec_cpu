use qsim_statevec_cpu::{NoiseModel, QuantumOp};

fn main() {
    let instructions = vec![
        (QuantumOp::Hadamard, 0),
        (QuantumOp::PauliX, 1),
        (QuantumOp::PauliZ, 2),
    ];

    let noise_model = NoiseModel {
        gate_error_prob: 0.02,
        readout_flip_prob: 0.01,
    };

    let avg_noisy_layer =
        qsim_statevec_cpu::execute_shots_noisy(instructions, 3, 1000000, noise_model).unwrap();
    let avg = avg_noisy_layer.measure_qubits();

    println!("Average noisy measured qubits: {avg:#?}");
}
