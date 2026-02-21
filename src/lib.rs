//! Quantum circuit simulator
//!
//! Provides an abstraction for quantum circuit simulations.
//! Uses the state vector simulation method.
//! Memory consumption is 2 * 8 * 2<sup>`num_qubits`</sup> bytes. For example, simulating 25 qubits cost ~537 MB.
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
//! if let Err(e) = q_layer.execute(instructions) {
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

use num::pow;
use num::Complex;
use serde::de;
use serde::de::{VariantAccess, Visitor};
use serde::{Deserialize, Deserializer, Serialize, Serializer};
use std::fmt;
use std::fmt::Write;

/// Supported quantum operations, equivalent to quantum gates in a circuit.
/// Operations with 'Par' suffix are experimental multi-threaded implementations, not guaranteed to improve performance.
#[derive(Clone, PartialEq, Debug)]
pub enum QuantumOp {
    PauliX,
    PauliY,
    PauliZ,
    Hadamard,
}

impl Serialize for QuantumOp {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        match *self {
            QuantumOp::PauliX => serializer.serialize_unit_variant("QuantumOp", 0, "PauliX"),
            QuantumOp::PauliY => serializer.serialize_unit_variant("QuantumOp", 1, "PauliY"),
            QuantumOp::PauliZ => serializer.serialize_unit_variant("QuantumOp", 1, "PauliZ"),
            QuantumOp::Hadamard => serializer.serialize_unit_variant("QuantumOp", 1, "Hadamard"),
        }
    }
}

impl<'de> Deserialize<'de> for QuantumOp {
    fn deserialize<D>(deserializer: D) -> Result<QuantumOp, D::Error>
    where
        D: Deserializer<'de>,
    {
        pub enum Field {
            PauliX,
            PauliY,
            PauliZ,
            Hadamard,
        }
        impl<'de> Deserialize<'de> for Field {
            fn deserialize<D>(deserializer: D) -> Result<Field, D::Error>
            where
                D: Deserializer<'de>,
            {
                struct FieldVisitor;

                impl Visitor<'_> for FieldVisitor {
                    type Value = Field;

                    fn expecting(&self, formatter: &mut fmt::Formatter) -> fmt::Result {
                        formatter.write_str("`PauliX`, `PauliY`, `PauliZ`, `Hadamard`")
                    }

                    fn visit_str<E>(self, value: &str) -> Result<Field, E>
                    where
                        E: de::Error,
                    {
                        match value {
                            "PauliX" => Ok(Field::PauliX),
                            "PauliY" => Ok(Field::PauliY),
                            "PauliZ" => Ok(Field::PauliZ),
                            "Hadamard" => Ok(Field::Hadamard),

                            _ => Err(de::Error::unknown_variant(
                                value,
                                &["PauliX", "PauliY", "PauliZ", "Hadamard"],
                            )),
                        }
                    }
                }

                deserializer.deserialize_identifier(FieldVisitor)
            }
        }

        struct MyEnumVisitor;

        impl<'de> Visitor<'de> for MyEnumVisitor {
            type Value = QuantumOp;

            fn expecting(&self, formatter: &mut fmt::Formatter) -> fmt::Result {
                formatter.write_str("struct QuantumOp")
            }

            fn visit_enum<A>(self, data: A) -> Result<QuantumOp, A::Error>
            where
                A: de::EnumAccess<'de>,
            {
                let (field, variant) = data.variant::<Field>()?;
                match field {
                    Field::PauliX => variant.unit_variant().map(|()| QuantumOp::PauliX),
                    Field::PauliY => variant.unit_variant().map(|()| QuantumOp::PauliY),
                    Field::PauliZ => variant.unit_variant().map(|()| QuantumOp::PauliX),
                    Field::Hadamard => variant.unit_variant().map(|()| QuantumOp::Hadamard),
                }
            }
        }

        deserializer.deserialize_enum(
            "QuantumOp",
            &["PauliX", "PauliY", "PauliZ", "Hadamard"],
            MyEnumVisitor,
        )
    }
}

pub type TargetQubit = u32;
pub type MeasuredQubits = Vec<f64>;

/// Quantum Assembly parser module.
/// Supports a simple subset of `OpenQASM` 3.0 (<https://openqasm.com/versions/3.0/index.html>)
pub mod qasm_parser {
    use crate::QuantumOp;
    use crate::TargetQubit;

    pub struct ParsedInstruct {
        pub num_qubits: u32,
        pub ops: Vec<(QuantumOp, TargetQubit)>,
    }

    /// Parses the contents of a qasm file
    ///
    /// # Errors
    /// Returns error if encounters semantic errors in the qasm file contents
    pub fn parse(file_contents: &str) -> Result<ParsedInstruct, String> {
        // create a vector of strings split by newline
        let mut lines: Vec<String> = file_contents
            .split('\n')
            .map(std::borrow::ToOwned::to_owned)
            .collect();

        // remove all lines until qreg
        let mut remove_delim = 0;
        for line in &lines {
            if line.contains("qreg") {
                break;
            }
            remove_delim += 1;
        }

        // Parse the number of qubits
        let num_qubits: String = lines[0].chars().filter(|&c| c.is_numeric()).collect();
        let Ok(num_qubits) = num_qubits.parse::<u32>() else {
            return Err("Failed to parse the number of qubits!".to_owned());
        };

        // remove all lines before and including qreg
        lines.drain(0..=remove_delim);

        // Filter each newline
        let mut parsed_instructions: Vec<(QuantumOp, TargetQubit)> = vec![];
        for line in &lines {
            let operation: &str = match line.split_whitespace().next() {
                Some(operation) => operation,
                None => continue,
            };

            // parse qubit target list
            let target_qubits: TargetQubit = match operation {
                // fetch the only target qubit after the op string
                "x" | "y" | "z" | "h" => {
                    let filtered_line: String = line
                        .chars()
                        .filter(|&c| c.is_numeric() || c == ' ')
                        .collect();
                    let filtered_line: Vec<String> = filtered_line
                        .split(' ')
                        .map(std::borrow::ToOwned::to_owned)
                        .collect();
                    match filtered_line[1].parse::<TargetQubit>() {
                        Ok(x) => x,
                        Err(_) => return Err("Failed to parse the target qubit!".to_owned()),
                    }
                }
                // TODO: fetch two target qubits after the op string cx and cz
                // TODO: fetch three target qubits after the op string ccx
                _ => {
                    // trace!("Skipping unknown operation {}", other);
                    continue;
                }
            };

            // parse operation codes
            let operation: QuantumOp = match operation {
                "x" => QuantumOp::PauliX,
                "y" => QuantumOp::PauliY,
                "z" => QuantumOp::PauliZ,
                "h" => QuantumOp::Hadamard,
                other => return Err(format!("Operation Code {other} not recognized!")),
            };

            parsed_instructions.push((operation, target_qubits));
        }

        Ok(ParsedInstruct {
            num_qubits,
            ops: parsed_instructions,
        })
    }
}

/// The main abstraction of quantum circuit simulation.
/// Contains the complex values of each possible state.
#[derive(Clone, PartialEq)]
pub struct QubitLayer {
    main: Vec<Complex<f64>>,
    parity: Vec<Complex<f64>>,
}

impl QubitLayer {
    /// Executes multiple quantum assembly instructions.
    /// Receives a vector containing pairs of (`QuantumOp`, `TargetQubit`).
    ///
    /// # Examples
    /// ```
    /// use qsim_statevec_cpu::{QubitLayer, QuantumOp};
    ///
    /// let mut q_layer = QubitLayer::new(2);
    /// let instructions = vec![(QuantumOp::PauliX, 0), (QuantumOp::PauliX, 1)];
    /// q_layer.execute(instructions);
    ///
    /// // qubits 0 and 1 must be 1.0
    /// assert_eq!(q_layer.measure_qubits()[0], 1.0);
    /// assert_eq!(q_layer.measure_qubits()[1], 1.0);
    /// ```
    ///
    /// # Errors
    /// If operation target qubit is out of range.
    pub fn execute(
        &mut self,
        quantum_instructions: Vec<(QuantumOp, TargetQubit)>,
    ) -> Result<(), String> {
        for (op, target_qubit) in quantum_instructions {
            if target_qubit >= self.get_num_qubits() {
                return Err(format!(
                    "Target qubit {target_qubit:?} is out of range. Size of layer is {}",
                    self.get_num_qubits()
                )
                .to_owned());
            }
            match op {
                QuantumOp::PauliX => {
                    self.pauli_x(target_qubit);
                }
                QuantumOp::PauliY => {
                    self.pauli_y(target_qubit);
                }
                QuantumOp::PauliZ => {
                    self.pauli_z(target_qubit);
                }
                QuantumOp::Hadamard => {
                    self.hadamard(target_qubit);
                }
            }
        }
        Ok(())
    }

    /// Returns the number of qubits represented in the `QubitLayer`.  
    /// ```
    /// use qsim_statevec_cpu::QubitLayer;
    ///
    /// let num_qubits = 20;
    /// let q_layer = QubitLayer::new(num_qubits);
    /// assert_eq!(num_qubits, 20);
    /// ```
    #[must_use]
    pub fn get_num_qubits(&self) -> u32 {
        self.main.len().ilog2()
    }

    /// Returns the results of the operations performed in the `QubitLayer`.
    /// Equivalent to collapsing qubits to obtain its state.
    /// # Examples
    /// ```
    /// use qsim_statevec_cpu::{QubitLayer, QuantumOp};
    ///
    /// let mut q_layer = QubitLayer::new(20);
    /// q_layer.execute(vec![(QuantumOp::Hadamard, 0)]);
    /// println!("{:?}", q_layer.measure_qubits());
    ///
    /// ```
    #[must_use]
    pub fn measure_qubits(&self) -> MeasuredQubits {
        let num_qubits = self.get_num_qubits();
        let mut measured_qubits: Vec<f64> = vec![0.0; num_qubits as usize];

        for index_main in 0..self.main.len() {
            if self.main[index_main] == Complex::new(0.0, 0.0) {
                continue;
            }
            for (index_measured_qubits, value) in measured_qubits.iter_mut().enumerate() {
                // check if the state has a bit in common with the measured_qubit index
                // does not matter which it is, thats why >= 1
                if (index_main & Self::mask(index_measured_qubits)) > 0 {
                    *value += pow(self.main[index_main].norm(), 2);
                }
            }
        }
        measured_qubits
    }

    /// Creates a new `QubitLayer` representing `num_qubits` qubits.  
    /// # Examples
    /// ```
    /// use qsim_statevec_cpu::QubitLayer;
    ///
    /// let q_layer = QubitLayer::new(20);
    /// ```
    #[must_use]
    pub fn new(num_qubits: u32) -> Self {
        let mut main = vec![Complex::new(0.0, 0.0); 2_usize.pow(num_qubits)];
        main[0] = Complex::new(1.0, 0.0);

        Self {
            main,
            parity: vec![Complex::new(0.0, 0.0); 2_usize.pow(num_qubits)],
        }
    }

    fn hadamard(&mut self, target_qubit: u32) {
        let hadamard_const = 1.0 / std::f64::consts::SQRT_2;
        for state in 0..self.main.len() {
            if self.main[state] != Complex::new(0.0, 0.0) {
                if state & Self::mask(target_qubit as usize) != 0 {
                    self.parity[state] -= hadamard_const * self.main[state];
                } else {
                    self.parity[state] += hadamard_const * self.main[state];
                }
            }
        }
        for state in 0..self.main.len() {
            if self.main[state] != Complex::new(0.0, 0.0) {
                let target_state: usize = state ^ Self::mask(target_qubit as usize);
                self.parity[target_state] += hadamard_const * self.main[state];
            }
        }
        self.reset_parity_layer();
    }

    fn pauli_z(&mut self, target_qubit: u32) {
        for state in 0..self.main.len() {
            if self.main[state] != Complex::new(0.0, 0.0) {
                if state & Self::mask(target_qubit as usize) != 0 {
                    self.parity[state] = -self.main[state];
                } else {
                    self.parity[state] = self.main[state];
                }
            }
        }

        self.reset_parity_layer();
    }

    fn pauli_y(&mut self, target_qubit: u32) {
        for state in 0..self.main.len() {
            if self.main[state] != Complex::new(0.0, 0.0) {
                let target_state: usize = state ^ Self::mask(target_qubit as usize);
                // if |0>, scalar 1i applies to |1>
                // if |1>, scalar -1i
                // TODO: probabily room for optimization here
                if target_state & Self::mask(target_qubit as usize) != 0 {
                    self.parity[target_state] = self.main[state] * Complex::new(0.0, 1.0);
                } else {
                    self.parity[target_state] = self.main[state] * Complex::new(0.0, -1.0);
                }
            }
        }
        self.reset_parity_layer();
    }

    fn pauli_x(&mut self, target_qubit: u32) {
        for state in 0..self.main.len() {
            if self.main[state] != Complex::new(0.0, 0.0) {
                let mut target_state: usize = state;
                target_state ^= Self::mask(target_qubit as usize); // flip bit 0
                self.parity[target_state] = self.main[state];
            }
        }

        self.reset_parity_layer();
    }

    fn reset_parity_layer(&mut self) {
        // clone parity qubit layer to qubit layer
        // self.main = self.parity.clone();
        self.main.clone_from(&self.parity);

        // reset parity qubit layer with 0
        self.parity
            .iter_mut()
            .map(|x| *x = Complex::new(0.0, 0.0))
            .count();
    }

    fn mask(position: usize) -> usize {
        0x1usize << position
    }
}

impl fmt::Debug for QubitLayer {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut output = String::new();
        for index_main in 0..self.main.len() {
            writeln!(
                &mut output,
                "state {:b} -> {}",
                index_main, self.main[index_main]
            )?;
        }
        write!(f, "{output}")
    }
}

impl fmt::Display for QubitLayer {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut output = String::new();
        for state in &self.main {
            let str = format!("{output} {state}");
            output = str;
        }
        output.remove(0);
        write!(f, "{output}")
    }
}

#[cfg(test)]
mod tests;
