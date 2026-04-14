use qsim_statevec_cpu::openq3_parser;
use qsim_statevec_cpu::QubitLayer;
use std::fs;
use std::path::Path;

fn format_probability(value: f64) -> String {
    let rounded_int = value.round();
    if (value - rounded_int).abs() < 1e-12 {
        return (rounded_int as i64).to_string();
    }

    let rounded = (value * 1_000_000_000_000.0).round() / 1_000_000_000_000.0;
    let mut output = format!("{rounded:.12}");
    while output.contains('.') && output.ends_with('0') {
        output.pop();
    }
    if output.ends_with('.') {
        output.pop();
    }
    output
}

fn format_measurements(measured: &[f64]) -> String {
    measured
        .iter()
        .map(|&value| format_probability(value))
        .collect::<Vec<_>>()
        .join(",")
}

fn run_case(name: &str, qasm_contents: &str) {
    let parsed = openq3_parser::parse(qasm_contents)
        .unwrap_or_else(|error| panic!("{name} should parse without errors: {error}"));

    let num_qubits = parsed.num_qubits;
    let instructions = parsed.ops;

    let mut layer = QubitLayer::new(num_qubits);
    layer
        .execute_noiseless(&instructions)
        .unwrap_or_else(|error| panic!("{name} should execute without errors: {error}"));

    let measured = layer.measure_qubits();
    println!("{name}={}", format_measurements(&measured));
}

fn main() {
    let circuits_dir = Path::new("functional_tests/reference_qasm");
    let mut qasm_files = fs::read_dir(circuits_dir)
        .unwrap_or_else(|error| panic!("failed to read {}: {error}", circuits_dir.display()))
        .filter_map(Result::ok)
        .map(|entry| entry.path())
        .filter(|path| path.extension().is_some_and(|ext| ext == "openqasm"))
        .collect::<Vec<_>>();

    qasm_files.sort();

    assert!(
        !qasm_files.is_empty(),
        "no .openqasm files found in {}",
        circuits_dir.display()
    );

    for qasm_file in qasm_files {
        let name = qasm_file
            .file_stem()
            .and_then(|stem| stem.to_str())
            .unwrap_or_else(|| panic!("invalid file name: {}", qasm_file.display()));
        let qasm_contents = fs::read_to_string(&qasm_file)
            .unwrap_or_else(|error| panic!("failed to read {}: {error}", qasm_file.display()));

        run_case(name, &qasm_contents);
    }
}
