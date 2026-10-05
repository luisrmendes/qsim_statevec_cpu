//! Quantum circuit simulator
//!
//! Provides an abstraction for quantum circuit simulations.
//! Uses the state vector simulation method running on CPU using main system memory.
//! Memory consumption is 8 * 2 * 2<sup>`num_qubits`</sup> bytes. For example, simulating 25 qubits costs ~537 MB.
//!
//! # Example
//!
//! ```
//! use qsim_statevec_cpu::{QSimGate, QubitLayer, SingleQubitOp};
//!
//! let mut q_layer: QubitLayer = QubitLayer::new(20);
//!
//! let instructions = vec![
//!     QSimGate::Single { op: SingleQubitOp::PauliX, target: 0 },
//!     QSimGate::Single { op: SingleQubitOp::PauliY, target: 1 },
//!     QSimGate::Single { op: SingleQubitOp::PauliZ, target: 2 },
//!     QSimGate::Single { op: SingleQubitOp::Hadamard, target: 3 },
//! ];
//!
//! if let Err(e) = q_layer.execute_instructions(instructions) {
//!     panic!("Failed to execute instructions! Error: {e}");
//! }
//!
//! let measured_qubits = q_layer.measure_qubits();
//! println!("{:?}", measured_qubits);
//!
//! // Check if the first qubit has flipped to 1 due to the Pauli X operation
//! assert_eq!(1.0, measured_qubits[0].round()); // measures might come with floating-point precision loss
//!
//! ```

pub mod qubit_layer;
pub mod types;

#[cfg(test)]
mod tests;

pub use qubit_layer::QubitLayer;
use rand::RngExt;
pub use types::*;

/// Executes multiple shots with stochastic noise.
/// Injects random errors in each axis of the qubit by inserting Pauli gates in the circuit gate list.
/// Randomly flips the measured qubit probabilities.
///
/// - `gate_error_prob`: after each gate, applies a random Pauli error (`X`, `Y`, or `Z`) on the same target qubit.
/// - `readout_flip_prob`: before measurement, applies a stochastic bit-flip (`X`) per qubit.
///
/// Returns the accumulated noisy layer averaged by the number of shots.
///
/// # Errors
/// Returns error if operation target qubit is out of range or if noise probabilities are outside `[0.0, 1.0]`.
pub fn execute_noisy_shots(
    circuit: QSimCircuit,
    shots: u32,
    noise_model: NoiseModel,
) -> Result<MeasuredQubits, String> {
    if !noise_model.is_valid() {
        return Err("Noise probabilities must be in the range [0.0, 1.0]".to_owned());
    }
    if shots == 0 {
        return Err("Number of shots must be greater than 0".to_owned());
    }

    let mut accumulated_results: MeasuredQubits = vec![0.0; circuit.num_qubits as usize].into();

    // Inject errors in the circuit as gates
    let mut rng = rand::rng();
    let mut gates_with_noise: Vec<QSimGate> = Vec::with_capacity(circuit.gates.len());
    for it in circuit.gates {
        gates_with_noise.push(it);

        // Affect one of qubit's 3 axis via Pauli gate application
        if rng.random::<f64>() < noise_model.gate_error_prob {
            let pauli = match rng.random_range(0..3) {
                0 => SingleQubitOp::PauliX,
                1 => SingleQubitOp::PauliY,
                _ => SingleQubitOp::PauliZ,
            };
            gates_with_noise.push(QSimGate::Single {
                op: pauli,
                target: it.target(),
            });
        }
    }

    for _ in 0..shots {
        let mut results = execute_noiseless(QSimCircuit {
            num_qubits: circuit.num_qubits,
            gates: gates_with_noise.clone(),
        })?;

        // Randomly flips the probability of the result
        for p in results.iter_mut() {
            if rng.random::<f64>() < noise_model.readout_flip_prob {
                *p = 1.0 - *p;
            }
        }

        accumulated_results += &results;
    }

    if shots > 0 {
        accumulated_results /= shots.into();
    }

    Ok(accumulated_results)
}

/// Executes noiseless quantum circuit.
/// Receives a ParsedCircuit, converts to QSim gate types and executes on a QubitLayer.
///
/// # Examples
/// ```
/// use qsim_statevec_cpu::{execute_noiseless, QSimCircuit, QSimGate, SingleQubitOp};
///
/// let circuit = QSimCircuit {
///     num_qubits: 2,
///     gates: vec![
///         QSimGate::Single { op: SingleQubitOp::PauliX, target: 0 },
///         QSimGate::Single { op: SingleQubitOp::PauliX, target: 1 },
///     ],
/// };
///
/// let results = execute_noiseless(circuit).expect("circuit should execute");
///
/// // qubits 0 and 1 must be 1.0
/// assert_eq!(results[0], 1.0);
/// assert_eq!(results[1], 1.0);
/// ```
///
/// # Errors
/// Returns error if operation target qubit is out of range.
pub fn execute_noiseless(circuit: QSimCircuit) -> Result<MeasuredQubits, String> {
    let mut layer = QubitLayer::new(circuit.num_qubits);

    layer.execute_instructions(circuit.gates)?;

    Ok(layer.measure_qubits())
}
