//! Quantum circuit simulator
//!
//! Provides an abstraction for quantum circuit simulations.
//! Uses the state vector simulation method running on CPU using main system memory.
//! Memory consumption is 8 * 2 * 2<sup>`num_qubits`</sup> bytes. For example, simulating 25 qubits costs ~537 MB.
//!
//! # Example
//!
//! ```
//! use qsim_statevec_cpu::{QubitLayer, QuantumOp};
//!
//! let mut q_layer: QubitLayer = QubitLayer::new(20);
//!
//! let instructions = vec![
//!     (QuantumOp::PauliX, 0),
//!     (QuantumOp::PauliY, 1),
//!     (QuantumOp::PauliZ, 2),
//!     (QuantumOp::Hadamard, 3),
//! ];
//!
//! if let Err(e) = q_layer.execute_noiseless(&instructions) {
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

pub mod openq3_parser;
pub mod qubit_layer;
pub mod types;

pub use qubit_layer::QubitLayer;
pub use types::*;
