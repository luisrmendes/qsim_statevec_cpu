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
//! if let Err(e) = q_layer.execute_noiseless(instructions) {
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
use rand::Rng;
use std::fmt;
use std::fmt::Write;
use std::ops::Add;
use std::ops::AddAssign;

/// Supported quantum operations, equivalent to quantum gates in a circuit.
/// Operations with 'Par' suffix are experimental multi-threaded implementations, not guaranteed to improve performance.
#[derive(Clone, PartialEq, Debug)]
pub enum SingleQubitOp {
    PauliX,
    PauliY,
    PauliZ,
    Hadamard,
    S,
    T,
}

pub type QuantumOp = SingleQubitOp;

#[derive(Clone, PartialEq, Debug)]
pub enum SingleCtrlQubitOp {
    ControlledX,
    ControlledZ,
}

#[derive(Clone, PartialEq, Debug)]
pub enum TwoCtrlQubitOp {
    Toffoli,
}

pub type SingleQubitInstruct = (SingleQubitOp, TargetQubit);
pub type SingleCtrlQubitInstruct = (SingleCtrlQubitOp, TargetQubit, CtrlQubit);
pub type TwoCtrlQubitInstruct = (TwoCtrlQubitOp, TargetQubit, CtrlQubit, CtrlQubit);

#[derive(Clone, PartialEq, Debug)]
pub enum QInstruct {
    Single(SingleQubitInstruct),
    SingleCtrl(SingleCtrlQubitInstruct),
    TwoCtrl(TwoCtrlQubitInstruct),
}

impl From<(QuantumOp, TargetQubit)> for QInstruct {
    fn from(value: (QuantumOp, TargetQubit)) -> Self {
        QInstruct::Single(value)
    }
}

impl From<(SingleCtrlQubitOp, TargetQubit, CtrlQubit)> for QInstruct {
    fn from(value: (SingleCtrlQubitOp, TargetQubit, CtrlQubit)) -> Self {
        QInstruct::SingleCtrl(value)
    }
}

impl From<(TwoCtrlQubitOp, TargetQubit, CtrlQubit, CtrlQubit)> for QInstruct {
    fn from(value: (TwoCtrlQubitOp, TargetQubit, CtrlQubit, CtrlQubit)) -> Self {
        QInstruct::TwoCtrl(value)
    }
}

pub type QInstructs = Vec<QInstruct>;
pub type TargetQubit = u32;
pub type CtrlQubit = u32;
pub type MeasuredQubits = Vec<f64>;

#[derive(Clone, Copy, Debug)]
pub struct NoiseModel {
    pub gate_error_prob: f64,
    pub readout_flip_prob: f64,
}

impl NoiseModel {
    fn is_valid(self) -> bool {
        (0.0..=1.0).contains(&self.gate_error_prob) && (0.0..=1.0).contains(&self.readout_flip_prob)
    }
}

/// The main abstraction of quantum circuit simulation.
/// Contains the complex values of each possible state.
#[derive(Clone, PartialEq)]
pub struct QubitLayer {
    main: Vec<Complex<f64>>,
    parity: Vec<Complex<f64>>,
    num_qubits: u32,
}

impl QubitLayer {
    /// Executes multiple shots with stochastic noise.
    ///
    /// - `gate_error_prob`: after each gate, applies a random Pauli error (`X`, `Y`, or `Z`) on the same target qubit.
    /// - `readout_flip_prob`: before measurement, applies a stochastic bit-flip (`X`) per qubit.
    ///
    /// Returns the accumulated noisy layer averaged by the number of shots.
    ///
    /// # Errors
    /// Returns error if operation target qubit is out of range or if noise probabilities are outside `[0.0, 1.0]`.
    pub fn execute_noisy_shots<T>(
        &mut self,
        quantum_instructions: &[T],
        shots: u32,
        noise_model: NoiseModel,
    ) -> Result<(), String>
    where
        T: Clone + Into<QInstruct>,
    {
        if !noise_model.is_valid() {
            return Err("Noise probabilities must be in the range [0.0, 1.0]".to_owned());
        }

        let mut rng = rand::thread_rng();
        let mut accumulated_qubit_layer = QubitLayer::new(self.get_num_qubits());

        for _ in 0..shots {
            let mut qubit_layer = QubitLayer::new(self.get_num_qubits());

            for instruction in quantum_instructions.iter().cloned() {
                let target_qubit = qubit_layer.execute_instruction(instruction.into())?;

                if rng.gen::<f64>() < noise_model.gate_error_prob {
                    match rng.gen_range(0..3) {
                        0 => qubit_layer.pauli_x(target_qubit),
                        1 => qubit_layer.pauli_y(target_qubit),
                        _ => qubit_layer.pauli_z(target_qubit),
                    }
                }
            }

            for qubit in 0..self.get_num_qubits() {
                if rng.gen::<f64>() < noise_model.readout_flip_prob {
                    qubit_layer.pauli_x(qubit);
                }
            }

            accumulated_qubit_layer += &qubit_layer;
        }

        if shots > 0 {
            accumulated_qubit_layer.scale_amplitudes(shots);
        }

        self.main = accumulated_qubit_layer.main;
        Ok(())
    }

    /// Executes multiple quantum assembly instructions.
    /// Receives a vector containing pairs of (`QuantumOp`, `TargetQubit`).
    ///
    /// # Examples
    /// ```
    /// use qsim_statevec_cpu::{QubitLayer, QuantumOp};
    ///
    /// let mut q_layer = QubitLayer::new(2);
    /// let instructions = vec![(QuantumOp::PauliX, 0), (QuantumOp::PauliX, 1)];
    /// q_layer.execute_noiseless(instructions);
    ///
    /// // qubits 0 and 1 must be 1.0
    /// assert_eq!(q_layer.measure_qubits()[0], 1.0);
    /// assert_eq!(q_layer.measure_qubits()[1], 1.0);
    /// ```
    ///
    /// # Errors
    /// If operation target qubit is out of range.
    pub fn execute_noiseless<I, T>(&mut self, quantum_instructions: I) -> Result<(), String>
    where
        I: IntoIterator<Item = T>,
        T: Into<QInstruct>,
    {
        for instruction in quantum_instructions {
            let _ = self.execute_instruction(instruction.into())?;
        }
        Ok(())
    }

    fn execute_instruction(&mut self, instruction: QInstruct) -> Result<TargetQubit, String> {
        match instruction {
            QInstruct::Single((op, target_qubit)) => {
                if target_qubit >= self.get_num_qubits() {
                    return Err(format!(
                        "Target qubit {target_qubit:?} is out of range. Size of layer is {}",
                        self.get_num_qubits()
                    ));
                }

                match op {
                    QuantumOp::PauliX => self.pauli_x(target_qubit),
                    QuantumOp::PauliY => self.pauli_y(target_qubit),
                    QuantumOp::PauliZ => self.pauli_z(target_qubit),
                    QuantumOp::Hadamard => self.hadamard(target_qubit),
                    QuantumOp::S => self.s_gate(target_qubit),
                    QuantumOp::T => self.t_gate(target_qubit),
                }

                Ok(target_qubit)
            }
            QInstruct::SingleCtrl((op, target_qubit, control_qubit)) => {
                if target_qubit >= self.get_num_qubits() {
                    return Err(format!(
                        "Target qubit {target_qubit:?} is out of range. Size of layer is {}",
                        self.get_num_qubits()
                    ));
                }
                if control_qubit >= self.get_num_qubits() {
                    return Err(format!(
                        "Control qubit {control_qubit:?} is out of range. Size of layer is {}",
                        self.get_num_qubits()
                    ));
                }

                match op {
                    SingleCtrlQubitOp::ControlledX => {
                        self.controlled_x(control_qubit, target_qubit);
                    }
                    SingleCtrlQubitOp::ControlledZ => {
                        self.controlled_z(control_qubit, target_qubit);
                    }
                }

                Ok(target_qubit)
            }
            QInstruct::TwoCtrl((op, target_qubit, control_qubit1, control_qubit2)) => {
                if target_qubit >= self.get_num_qubits() {
                    return Err(format!(
                        "Target qubit {target_qubit:?} is out of range. Size of layer is {}",
                        self.get_num_qubits()
                    ));
                }
                if control_qubit1 >= self.get_num_qubits() {
                    return Err(format!(
                        "Control qubit {control_qubit1:?} is out of range. Size of layer is {}",
                        self.get_num_qubits()
                    ));
                }
                if control_qubit2 >= self.get_num_qubits() {
                    return Err(format!(
                        "Control qubit {control_qubit2:?} is out of range. Size of layer is {}",
                        self.get_num_qubits()
                    ));
                }

                match op {
                    TwoCtrlQubitOp::Toffoli => {
                        self.toffoli(control_qubit1, control_qubit2, target_qubit);
                    }
                }

                Ok(target_qubit)
            }
        }
    }

    /// Returns the estimated memory usage in bytes (`8 * 2 * 2^num_qubits`).
    #[must_use]
    pub fn get_mem_usage(&self) -> u64 {
        (8_u64 * 2_u64) * (2_u64.pow(self.num_qubits))
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
        self.num_qubits
    }

    /// Returns the results of the operations performed in the `QubitLayer`.
    /// Equivalent to collapsing qubits to obtain its state.
    /// # Examples
    /// ```
    /// use qsim_statevec_cpu::{QubitLayer, QuantumOp};
    ///
    /// let mut q_layer = QubitLayer::new(20);
    /// q_layer.execute_noiseless(vec![(QuantumOp::Hadamard, 0)]);
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
            num_qubits,
        }
    }

    fn toffoli(&mut self, control_qubit1: u32, control_qubit2: u32, target_qubit: u32) {
        for state in 0..self.main.len() {
            if self.main[state] != Complex::new(0.0, 0.0) {
                if state & Self::mask(control_qubit1 as usize) != 0
                    && state & Self::mask(control_qubit2 as usize) != 0
                {
                    let target_state: usize = state ^ Self::mask(target_qubit as usize);
                    self.parity[target_state] = self.main[state];
                } else {
                    self.parity[state] = self.main[state];
                }
            }
        }

        self.reset_parity_layer();
    }

    fn controlled_z(&mut self, control_qubit: u32, target_qubit: u32) {
        for state in 0..self.main.len() {
            if self.main[state] != Complex::new(0.0, 0.0) {
                if state & Self::mask(control_qubit as usize) != 0
                    && state & Self::mask(target_qubit as usize) != 0
                {
                    self.parity[state] = -self.main[state];
                } else {
                    self.parity[state] = self.main[state];
                }
            }
        }

        self.reset_parity_layer();
    }

    fn controlled_x(&mut self, control_qubit: u32, target_qubit: u32) {
        for state in 0..self.main.len() {
            if self.main[state] != Complex::new(0.0, 0.0) {
                if state & Self::mask(control_qubit as usize) != 0 {
                    let target_state: usize = state ^ Self::mask(target_qubit as usize);
                    self.parity[target_state] = self.main[state];
                } else {
                    self.parity[state] = self.main[state];
                }
            }
        }

        self.reset_parity_layer();
    }

    fn s_gate(&mut self, target_qubit: u32) {
        for state in 0..self.main.len() {
            if self.main[state] != Complex::new(0.0, 0.0) {
                if state & Self::mask(target_qubit as usize) != 0 {
                    self.parity[state] = Complex::new(0.0, 1.0) * self.main[state];
                } else {
                    self.parity[state] = self.main[state];
                }
            }
        }

        self.reset_parity_layer();
    }

    fn t_gate(&mut self, target_qubit: u32) {
        let t_const = Complex::from_polar(1.0, std::f64::consts::FRAC_PI_4);
        for state in 0..self.main.len() {
            if self.main[state] != Complex::new(0.0, 0.0) {
                if state & Self::mask(target_qubit as usize) != 0 {
                    self.parity[state] = t_const * self.main[state];
                } else {
                    self.parity[state] = self.main[state];
                }
            }
        }

        self.reset_parity_layer();
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

    /// Scales the amplitudes of the `QubitLayer` by a factor of `scale`.
    fn scale_amplitudes(&mut self, scale: u32) {
        let scale = f64::from(scale);
        for amplitude in &mut self.main {
            *amplitude /= scale;
        }

        for amplitude in &mut self.parity {
            *amplitude /= scale;
        }
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
                "|{:b}>\t -> {}",
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

impl Add for QubitLayer {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        assert_eq!(
            self.get_num_qubits(),
            rhs.get_num_qubits(),
            "Cannot add QubitLayers with different numbers of qubits"
        );

        let main = self
            .main
            .into_iter()
            .zip(rhs.main)
            .map(|(lhs, rhs)| lhs + rhs)
            .collect();

        let parity = self
            .parity
            .into_iter()
            .zip(rhs.parity)
            .map(|(lhs, rhs)| lhs + rhs)
            .collect();

        Self {
            main,
            parity,
            num_qubits: self.num_qubits,
        }
    }
}

impl Add<&QubitLayer> for &QubitLayer {
    type Output = QubitLayer;

    fn add(self, rhs: &QubitLayer) -> Self::Output {
        assert_eq!(
            self.main.len(),
            rhs.main.len(),
            "Cannot add QubitLayers with different numbers of qubits"
        );

        let main = self
            .main
            .iter()
            .zip(rhs.main.iter())
            .map(|(lhs, rhs)| *lhs + *rhs)
            .collect();

        let parity = self
            .parity
            .iter()
            .zip(rhs.parity.iter())
            .map(|(lhs, rhs)| *lhs + *rhs)
            .collect();

        QubitLayer {
            main,
            parity,
            num_qubits: self.num_qubits,
        }
    }
}

impl AddAssign<&QubitLayer> for QubitLayer {
    fn add_assign(&mut self, rhs: &QubitLayer) {
        assert_eq!(
            self.main.len(),
            rhs.main.len(),
            "Cannot add QubitLayers with different numbers of qubits"
        );

        for (lhs, rhs_value) in self.main.iter_mut().zip(rhs.main.iter()) {
            *lhs += *rhs_value;
        }

        for (lhs, rhs_value) in self.parity.iter_mut().zip(rhs.parity.iter()) {
            *lhs += *rhs_value;
        }
    }
}

impl AddAssign for QubitLayer {
    fn add_assign(&mut self, rhs: Self) {
        *self += &rhs;
    }
}

/// Quantum Assembly parser.
/// Supports a simple subset of `OpenQASM` 3.0 (<https://openqasm.com/versions/3.0/index.html>)
pub mod openq3_parser {
    use crate::QInstruct;
    use crate::QInstructs;
    use crate::QuantumOp;
    use crate::SingleCtrlQubitOp;
    use crate::TargetQubit;
    use crate::TwoCtrlQubitOp;

    pub struct ParsedInstruct {
        pub num_qubits: u32,
        pub ops: QInstructs,
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

        // find qreg declaration line
        let Some(remove_delim) = lines.iter().position(|line| line.contains("qreg")) else {
            return Err("Failed to parse the number of qubits!".to_owned());
        };

        // Parse the number of qubits
        let num_qubits: String = lines[remove_delim]
            .chars()
            .filter(|&c| c.is_numeric())
            .collect();
        let Ok(num_qubits) = num_qubits.parse::<u32>() else {
            return Err("Failed to parse the number of qubits!".to_owned());
        };

        // remove all lines before and including qreg
        lines.drain(0..=remove_delim);

        // Filter each newline
        let mut parsed_instructions: QInstructs = vec![];
        for line in &lines {
            let operation: &str = match line.split_whitespace().next() {
                Some(operation) => operation,
                None => continue,
            };

            let filtered_line: String = line
                .chars()
                .map(|c| if c.is_ascii_digit() { c } else { ' ' })
                .collect();
            let qubits: Vec<u32> = filtered_line
                .split_whitespace()
                .filter_map(|x| x.parse::<u32>().ok())
                .collect();

            if operation == "cx" || operation == "cz" {
                if qubits.len() != 2 {
                    return Err("Failed to parse control and target qubits!".to_owned());
                }

                let op = if operation == "cx" {
                    SingleCtrlQubitOp::ControlledX
                } else {
                    SingleCtrlQubitOp::ControlledZ
                };

                parsed_instructions.push(QInstruct::SingleCtrl((op, qubits[0], qubits[1])));
                continue;
            }

            if operation == "ccx" {
                if qubits.len() != 3 {
                    return Err("Failed to parse two controls and target qubits!".to_owned());
                }

                parsed_instructions.push(QInstruct::TwoCtrl((
                    TwoCtrlQubitOp::Toffoli,
                    qubits[0],
                    qubits[1],
                    qubits[2],
                )));
                continue;
            }

            // parse qubit target list
            let target_qubits: TargetQubit = match operation {
                // fetch the only target qubit after the op string
                "x" | "y" | "z" | "h" | "s" | "t" => {
                    let Some(&target) = qubits.first() else {
                        return Err("Failed to parse the target qubit!".to_owned());
                    };
                    target
                }

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
                "s" => QuantumOp::S,
                "t" => QuantumOp::T,
                other => return Err(format!("Operation Code {other} not recognized!")),
            };

            parsed_instructions.push(QInstruct::Single((operation, target_qubits)));
        }

        Ok(ParsedInstruct {
            num_qubits,
            ops: parsed_instructions,
        })
    }
}

#[cfg(test)]
mod tests;
