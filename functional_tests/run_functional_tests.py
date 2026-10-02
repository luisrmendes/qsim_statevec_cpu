#!/usr/bin/env python3
'''
This file defines functional tests comparing the qsim_statevec_cpu simulator lib against Qiskit's backend simulator
'''
from __future__ import annotations

import math
import subprocess
import sys
from pathlib import Path

# Both sides print probabilities rounded to 12 decimal places
ABS_TOLERANCE = 1e-9

USE_COLOR = sys.stdout.isatty()
GREEN = "\033[92m" if USE_COLOR else ""
RED = "\033[91m" if USE_COLOR else ""
BOLD = "\033[1m" if USE_COLOR else ""
RESET = "\033[0m" if USE_COLOR else ""


def run_command(command: list[str], cwd: Path) -> list[str]:
    try:
        completed = subprocess.run(command, cwd=cwd, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as error:
        print(f"Command failed with exit code {error.returncode}: {' '.join(command)}", file=sys.stderr)
        print(error.stderr, file=sys.stderr)
        raise SystemExit(1) from error
    return [line for line in completed.stdout.splitlines() if line]


def parse_results(lines: list[str], source: str) -> dict[str, str]:
    '''Parses lines in the format NAME=RESULTS into a {NAME: RESULTS} dict'''
    results: dict[str, str] = {}
    for line in lines:
        name, sep, result = line.partition("=")
        if not sep:
            raise SystemExit(f"Malformed {source} output line (expected NAME=RESULTS): {line!r}")
        if name in results:
            raise SystemExit(f"Duplicate test name in {source} output: {name}")
        results[name] = result
    return results


def results_match(qiskit_result: str, qsim_result: str) -> bool:
    qiskit_values = [float(value) for value in qiskit_result.split(",")]
    qsim_values = [float(value) for value in qsim_result.split(",")]
    return len(qiskit_values) == len(qsim_values) and all(
        math.isclose(a, b, rel_tol=0.0, abs_tol=ABS_TOLERANCE) for a, b in zip(qiskit_values, qsim_values)
    )


def has_qiskit(python_bin: Path) -> bool:
    if not python_bin.is_file():
        return False
    return subprocess.run([str(python_bin), "-c", "import qiskit"], capture_output=True, check=False).returncode == 0


def main() -> int:
    script_dir = Path(__file__).resolve().parent
    repo_root = script_dir.parent
    venv_dir = script_dir / ".venv"
    python_bin = venv_dir / "bin" / "python"
    python_script = script_dir / "qiskit_tests.py"
    requirements = script_dir / "requirements.txt"

    if not has_qiskit(python_bin):
        print("Qiskit virtual environment not found or incomplete. Setting it up...")
        subprocess.run([sys.executable, "-m", "venv", str(venv_dir)], check=True)
        subprocess.run([str(python_bin), "-m", "pip", "install", "-r", str(requirements)], check=True)

    qsim = parse_results(
        run_command(["cargo", "run", "--quiet", "--example", "qsim_tests"], repo_root), "qsim_statevec_cpu"
    )
    qiskit = parse_results(run_command([str(python_bin), str(python_script)], repo_root), "Qiskit")

    failing_tests: list[tuple[str, str, str]] = []
    print("Test results (Qiskit vs qsim_statevec_cpu):\n")
    for name in sorted(qiskit.keys() | qsim.keys()):
        qiskit_result = qiskit.get(name, "<missing>")
        qsim_result = qsim.get(name, "<missing>")
        print(f"{BOLD}{name}{RESET}")
        print(f"  Qiskit:            {qiskit_result}")
        print(f"  qsim_statevec_cpu: {qsim_result}")
        if name in qiskit and name in qsim and results_match(qiskit_result, qsim_result):
            print(f"  Result:            {GREEN}PASS{RESET}")
        else:
            print(f"  Result:            {RED}FAIL{RESET}")
            failing_tests.append((name, qiskit_result, qsim_result))

    if not failing_tests:
        print("\nAll tests passed: qsim_statevec_cpu output matches Qiskit for all reference circuits.")
        return 0

    print(f"\n{len(failing_tests)} test(s) failed. Failing tests summary:")
    for name, qiskit_result, qsim_result in failing_tests:
        print(f"{BOLD}{name}{RESET}")
        print(f"  Qiskit:            {qiskit_result}")
        print(f"  qsim_statevec_cpu: {qsim_result}")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
