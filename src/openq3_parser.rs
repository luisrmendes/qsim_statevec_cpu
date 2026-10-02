//! Quantum Assembly parser.
//! Supports a simple subset of `OpenQASM` 3.0 (<https://openqasm.com/versions/3.0/index.html>)

use crate::types::QInstruct;
use crate::types::QInstructs;
use crate::types::QuantumOp;
use crate::types::SingleCtrlQubitOp;
use crate::types::TargetQubit;
use crate::types::TwoCtrlQubitOp;

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

        if operation == "cx" || operation == "cy" || operation == "cz" {
            if qubits.len() != 2 {
                return Err("Failed to parse control and target qubits!".to_owned());
            }

            let op = if operation == "cx" {
                SingleCtrlQubitOp::ControlledX
            } else if operation == "cy" {
                SingleCtrlQubitOp::ControlledY
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
            "x" | "y" | "z" | "h" | "s" | "t" | "sx" | "sy" => {
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
            "sx" => QuantumOp::SX,
            "sy" => QuantumOp::SY,
            other => return Err(format!("Operation Code {other} not recognized!")),
        };

        parsed_instructions.push(QInstruct::Single((operation, target_qubits)));
    }

    Ok(ParsedInstruct {
        num_qubits,
        ops: parsed_instructions,
    })
}
