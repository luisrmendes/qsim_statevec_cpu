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

    let noisy_shots = qsim_statevec_cpu::execute_shots_noisy(instructions, 3, 1000, noise_model)
        .unwrap();

    let mut avg = vec![0.0; 3];
    for shot in &noisy_shots {
        for (index, value) in shot.iter().enumerate() {
            avg[index] += *value;
        }
    }
    for value in &mut avg {
        *value /= noisy_shots.len() as f64;
    }

    println!("Average noisy measured qubits: {avg:#?}");
}
