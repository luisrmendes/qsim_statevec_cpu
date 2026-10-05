use num::pow;
use num::Complex;
use std::fmt;
use std::fmt::Write;
use std::ops::Add;
use std::ops::AddAssign;
use std::ops::DivAssign;

use crate::types::*;

/// The main abstraction of quantum circuit simulation.
/// Contains the complex values of each possible state.
#[derive(Clone, PartialEq)]
pub struct QubitLayer {
    main: Vec<Complex<f64>>,
    parity: Vec<Complex<f64>>,
    num_qubits: u32,
}

impl QubitLayer {
    pub fn execute_instructions(&mut self, instructions: Vec<QSimGate>) -> Result<(), String> {
        for gate in instructions {
            match gate {
                QSimGate::Single { op, target } => {
                    if target >= self.get_num_qubits() {
                        return Err(format!(
                            "Target qubit {target:?} is out of range. Size of layer is {}",
                            self.get_num_qubits()
                        ));
                    }

                    match op {
                        SingleQubitOp::PauliX => self.pauli_x(target),
                        SingleQubitOp::PauliY => self.pauli_y(target),
                        SingleQubitOp::PauliZ => self.pauli_z(target),
                        SingleQubitOp::Hadamard => self.hadamard(target),
                        SingleQubitOp::S => self.s_gate(target),
                        SingleQubitOp::T => self.t_gate(target),
                        SingleQubitOp::SX => self.sqrt_pauli_x(target),
                        SingleQubitOp::SY => self.sqrt_pauli_y(target),
                    }
                }
                QSimGate::SingleCtrl {
                    op,
                    control,
                    target,
                } => {
                    if control >= self.get_num_qubits() {
                        return Err(format!(
                            "Control qubit {control:?} is out of range. Size of layer is {}",
                            self.get_num_qubits()
                        ));
                    }
                    if control >= self.get_num_qubits() {
                        return Err(format!(
                            "Target qubit {target:?} is out of range. Size of layer is {}",
                            self.get_num_qubits()
                        ));
                    }
                    if target == control {
                        return Err(format!(
                            "Target qubit and control qubit are the same: {target:?}"
                        ));
                    }

                    match op {
                        SingleCtrlQubitOp::ControlledX => {
                            self.controlled_x(control, target);
                        }
                        SingleCtrlQubitOp::ControlledY => {
                            self.controlled_y(control, target);
                        }
                        SingleCtrlQubitOp::ControlledZ => {
                            self.controlled_z(control, target);
                        }
                    }
                }
                QSimGate::TwoCtrl {
                    op,
                    controls,
                    target,
                } => {
                    if controls[0] >= self.get_num_qubits() {
                        return Err(format!(
                            "Control qubit {:?} is out of range. Size of layer is {}",
                            controls[0],
                            self.get_num_qubits()
                        ));
                    }
                    if controls[1] >= self.get_num_qubits() {
                        return Err(format!(
                            "Control qubit {:?} is out of range. Size of layer is {}",
                            controls[1],
                            self.get_num_qubits()
                        ));
                    }
                    if target >= self.get_num_qubits() {
                        return Err(format!(
                            "Target qubit {target:?} is out of range. Size of layer is {}",
                            self.get_num_qubits()
                        ));
                    }
                    if target == controls[0] {
                        return Err(format!(
                            "Target qubit and control qubit 1 are the same: {target:?}"
                        ));
                    }
                    if target == controls[1] {
                        return Err(format!(
                            "Target qubit and control qubit 2 are the same: {target:?}"
                        ));
                    }
                    if controls[0] == controls[1] {
                        return Err(format!(
                            "Control qubit 1 and 2 are the same: {:?}",
                            controls[0]
                        ));
                    }

                    match op {
                        TwoCtrlQubitOp::Toffoli => {
                            self.toffoli(controls[0], controls[1], target);
                        }
                    }
                }
            }
        }
        Ok(())
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
    /// use qsim_statevec_cpu::{QSimGate, QubitLayer, SingleQubitOp};
    ///
    /// let mut q_layer = QubitLayer::new(20);
    /// q_layer.execute_instructions(vec![QSimGate::Single { op: SingleQubitOp::Hadamard, target: 0 }]);
    /// println!("{:?}", q_layer.measure_qubits());
    ///
    /// ```
    #[must_use]
    pub fn measure_qubits(&self) -> MeasuredQubits {
        let num_qubits = self.get_num_qubits();
        let mut measured_qubits: MeasuredQubits = vec![0.0; num_qubits as usize].into();

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

    fn sqrt_pauli_x(&mut self, target_qubit: u32) {
        let const_same_state = Complex::new(0.5, 0.5);
        let const_flipped_state = Complex::new(0.5, -0.5);

        for state in 0..self.main.len() {
            if self.main[state] != Complex::new(0.0, 0.0) {
                let target_state: usize = state ^ Self::mask(target_qubit as usize);

                // |0> and |1> components both contribute with:
                // (1+i)/2 to the same basis index and (1-i)/2 to the flipped index.
                self.parity[state] += const_same_state * self.main[state];
                self.parity[target_state] += const_flipped_state * self.main[state];
            }
        }

        self.reset_parity_layer();
    }

    fn sqrt_pauli_y(&mut self, target_qubit: u32) {
        let const_same_state = Complex::new(0.5, 0.5);
        let const_zero_to_one = Complex::new(0.5, 0.5);
        let const_one_to_zero = Complex::new(-0.5, -0.5);

        for state in 0..self.main.len() {
            if self.main[state] != Complex::new(0.0, 0.0) {
                let target_state: usize = state ^ Self::mask(target_qubit as usize);

                self.parity[state] += const_same_state * self.main[state];

                if state & Self::mask(target_qubit as usize) == 0 {
                    self.parity[target_state] += const_zero_to_one * self.main[state];
                } else {
                    self.parity[target_state] += const_one_to_zero * self.main[state];
                }
            }
        }

        self.reset_parity_layer();
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

    fn controlled_y(&mut self, control_qubit: u32, target_qubit: u32) {
        for state in 0..self.main.len() {
            if self.main[state] != Complex::new(0.0, 0.0) {
                if state & Self::mask(control_qubit as usize) != 0 {
                    let target_state: usize = state ^ Self::mask(target_qubit as usize);
                    if state & Self::mask(target_qubit as usize) == 0 {
                        self.parity[target_state] = self.main[state] * Complex::new(0.0, 1.0);
                    } else {
                        self.parity[target_state] = self.main[state] * Complex::new(0.0, -1.0);
                    }
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

impl DivAssign<u32> for QubitLayer {
    fn div_assign(&mut self, rhs: u32) {
        assert_ne!(rhs, 0, "Cannot divide QubitLayer by zero");

        let scale = f64::from(rhs);
        for amplitude in &mut self.main {
            *amplitude /= scale;
        }

        for amplitude in &mut self.parity {
            *amplitude /= scale;
        }
    }
}

#[cfg(test)]
mod tests;
