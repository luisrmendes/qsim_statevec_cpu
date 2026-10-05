use std::ops::{AddAssign, Deref, DerefMut, DivAssign};

/// Supported quantum operations, equivalent to quantum gates in a circuit.
#[derive(Clone, Copy, PartialEq, Debug)]
pub enum SingleQubitOp {
    PauliX,
    PauliY,
    PauliZ,
    Hadamard,
    S,
    T,
    SX,
    SY,
}

#[derive(Clone, Copy, PartialEq, Debug)]
pub enum SingleCtrlQubitOp {
    ControlledX,
    ControlledY,
    ControlledZ,
}

#[derive(Clone, Copy, PartialEq, Debug)]
pub enum TwoCtrlQubitOp {
    Toffoli,
}

#[derive(Clone, Copy, PartialEq, Debug)]
pub enum QSimGate {
    Single {
        op: SingleQubitOp,
        target: TargetQubit,
    },
    SingleCtrl {
        op: SingleCtrlQubitOp,
        control: CtrlQubit,
        target: TargetQubit,
    },
    TwoCtrl {
        op: TwoCtrlQubitOp,
        controls: [CtrlQubit; 2],
        target: TargetQubit,
    },
}

impl QSimGate {
    /// The qubit the gate acts on (for controlled gates, the target, not the controls).
    pub fn target(&self) -> TargetQubit {
        match self {
            QSimGate::Single { target, .. }
            | QSimGate::SingleCtrl { target, .. }
            | QSimGate::TwoCtrl { target, .. } => *target,
        }
    }
}

impl From<(SingleQubitOp, TargetQubit)> for QSimGate {
    fn from((op, target): (SingleQubitOp, TargetQubit)) -> Self {
        QSimGate::Single { op, target }
    }
}

impl From<(SingleCtrlQubitOp, CtrlQubit, TargetQubit)> for QSimGate {
    fn from((op, control, target): (SingleCtrlQubitOp, CtrlQubit, TargetQubit)) -> Self {
        QSimGate::SingleCtrl {
            op,
            control,
            target,
        }
    }
}

impl From<(TwoCtrlQubitOp, CtrlQubit, CtrlQubit, TargetQubit)> for QSimGate {
    fn from(
        (op, control1, control2, target): (TwoCtrlQubitOp, CtrlQubit, CtrlQubit, TargetQubit),
    ) -> Self {
        QSimGate::TwoCtrl {
            op,
            controls: [control1, control2],
            target,
        }
    }
}

#[derive(Debug)]
pub struct QSimCircuit {
    pub num_qubits: u32,
    pub gates: Vec<QSimGate>,
}

pub type QSimGates = Vec<QSimGate>;
pub type TargetQubit = u32;
pub type CtrlQubit = u32;

/// Probability of each qubit being measured as 1.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct MeasuredQubits(pub Vec<f64>);

impl From<Vec<f64>> for MeasuredQubits {
    fn from(v: Vec<f64>) -> Self {
        Self(v)
    }
}

impl AddAssign<&MeasuredQubits> for MeasuredQubits {
    fn add_assign(&mut self, other: &MeasuredQubits) {
        assert_eq!(
            self.len(),
            other.len(),
            "MeasuredQubits must have the same length to be added"
        );
        for (a, b) in self.iter_mut().zip(other.iter()) {
            *a += b;
        }
    }
}

// For averaging over shots: `accumulated /= shots as f64`
impl DivAssign<f64> for MeasuredQubits {
    fn div_assign(&mut self, divisor: f64) {
        self.iter_mut().for_each(|p| *p /= divisor);
    }
}

// Let it behave like a Vec<f64> everywhere else: indexing, .len(), .iter(), slices
impl Deref for MeasuredQubits {
    type Target = Vec<f64>;
    fn deref(&self) -> &Vec<f64> {
        &self.0
    }
}
impl DerefMut for MeasuredQubits {
    fn deref_mut(&mut self) -> &mut Vec<f64> {
        &mut self.0
    }
}

#[derive(Clone, Copy, Debug)]
pub struct NoiseModel {
    pub gate_error_prob: f64,
    pub readout_flip_prob: f64,
}

impl NoiseModel {
    pub fn is_valid(self) -> bool {
        (0.0..=1.0).contains(&self.gate_error_prob) && (0.0..=1.0).contains(&self.readout_flip_prob)
    }
}
