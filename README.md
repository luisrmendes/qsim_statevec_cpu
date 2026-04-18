# qsim_statevec_cpu

<!-- Workflow status badges -->
[![CI](https://github.com/luisrmendes/qsim_statevec_cpu/actions/workflows/ci.yml/badge.svg?event=push&job=code_checks)](https://github.com/luisrmendes/qsim_statevec_cpu/actions/workflows/ci.yml)
[![License](https://img.shields.io/badge/license-MIT_OR_Apache--2.0-blue.svg)](https://github.com/luisrmendes/qsim_statevec_cpu#license)
[![Cargo](https://img.shields.io/crates/v/quantum_state_sim.svg)](https://docs.rs/quantum_state_sim/0.1.0/quantum_state_sim/)
[![Rust 1.70](https://img.shields.io/badge/rust-1.70-lightgray.svg)](
https://www.rust-lang.org)
[![Documentation](https://docs.rs/quantum_state_sim/badge.svg)](https://docs.rs/quantum_state_sim/0.1.0/quantum_state_sim/)

This provides a quantum simulation abstraction tool to simulate quantum circuits.  
Uses the state vector simulation method.  
Memory consumption is 8 * 2 * 2<sup>`num_qubits`</sup> bytes. For example, simulating 25 qubits costs ~537 MB.  

The following gate operations are implemented:

- Pauli X
- Pauli Y
- Pauli Z
- Hadamard
- T
- S
- Controlled X
- Controlled Y
- Controlled Z
- Toffoli
- Sqrt Pauli X
- Sqrt Pauli Y

## Dev Notes

### Git hooks

This repository includes a pre-commit hook in `.githooks/pre-commit` that runs:

- `cargo fmt -- --check`
- `cargo check`
- `cargo test -q`

To enable it:

```bash
git config core.hooksPath .githooks
chmod +x .githooks/pre-commit
```
