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
    SX,
    SY,
}

pub type QuantumOp = SingleQubitOp;

#[derive(Clone, PartialEq, Debug)]
pub enum SingleCtrlQubitOp {
    ControlledX,
    ControlledY,
    ControlledZ,
}

#[derive(Clone, PartialEq, Debug)]
pub enum TwoCtrlQubitOp {
    Toffoli,
}

#[derive(Clone, PartialEq, Debug)]
pub enum QInstruct {
    Single((SingleQubitOp, TargetQubit)),
    SingleCtrl((SingleCtrlQubitOp, CtrlQubit, TargetQubit)),
    TwoCtrl((TwoCtrlQubitOp, CtrlQubit, CtrlQubit, TargetQubit)),
}

impl From<(QuantumOp, TargetQubit)> for QInstruct {
    fn from(value: (QuantumOp, TargetQubit)) -> Self {
        QInstruct::Single(value)
    }
}

impl From<(SingleCtrlQubitOp, CtrlQubit, TargetQubit)> for QInstruct {
    fn from(value: (SingleCtrlQubitOp, CtrlQubit, TargetQubit)) -> Self {
        QInstruct::SingleCtrl(value)
    }
}

impl From<(TwoCtrlQubitOp, CtrlQubit, CtrlQubit, TargetQubit)> for QInstruct {
    fn from(value: (TwoCtrlQubitOp, CtrlQubit, CtrlQubit, TargetQubit)) -> Self {
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
    pub fn is_valid(self) -> bool {
        (0.0..=1.0).contains(&self.gate_error_prob) && (0.0..=1.0).contains(&self.readout_flip_prob)
    }
}
