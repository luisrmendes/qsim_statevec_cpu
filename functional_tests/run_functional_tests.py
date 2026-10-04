#!/usr/bin/env python3
# Qiskit ships without type stubs
# pyright: reportMissingTypeStubs=false
"""Compare qsim_statevec_cpu against Qiskit on every circuit in reference_qasm/.

Runs

Runs with any Python 3: on first use, Qiskit is installed into functional_tests/.venv
and the script re-executes itself with that interpreter.
"""

import math
import os
import subprocess
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent
CIRCUITS_DIR = SCRIPT_DIR / "reference_qasm"
VENV_PYTHON = SCRIPT_DIR / ".venv" / "bin" / "python"
ABS_TOLERANCE = 1e-9

try:
    from qiskit import qasm3
    from qiskit.quantum_info import Statevector
except ImportError:
    if os.environ.get("QSIM_FT_BOOTSTRAPPED"):
        raise
    print("Setting up the Qiskit virtual environment...", file=sys.stderr)
    subprocess.run([sys.executable, "-m", "venv", SCRIPT_DIR / ".venv"], check=True)
    subprocess.run(
        [
            VENV_PYTHON,
            "-m",
            "pip",
            "install",
            "-q",
            "-r",
            SCRIPT_DIR / "requirements.txt",
        ],
        check=True,
    )
    os.environ["QSIM_FT_BOOTSTRAPPED"] = "1"
    os.execv(VENV_PYTHON, [VENV_PYTHON, __file__, *sys.argv[1:]])

USE_COLOR = sys.stdout.isatty()
PASS = "\033[92mPASS\033[0m" if USE_COLOR else "PASS"
FAIL = "\033[91mFAIL\033[0m" if USE_COLOR else "FAIL"


def qiskit_results() -> dict[str, list[float]]:
    """P(qubit = 1) for each qubit of each reference circuit, computed by Qiskit."""
    results: dict[str, list[float]] = {}
    for path in sorted(CIRCUITS_DIR.glob("*.openqasm")):
        state = Statevector(qasm3.load(path))  # pyright: ignore[reportUnknownMemberType, reportUnknownArgumentType]
        num_qubits = state.num_qubits or 0  # only None for non-qubit dimensions
        results[path.stem] = [float(state.probabilities([q])[1]) for q in range(num_qubits)]
    return results


def qsim_results() -> dict[str, list[float]]:
    """Same as qiskit_results, from the qsim_tests example (lines of NAME=P0,P1,...)."""
    command = ["cargo", "run", "--quiet", "--example", "run_ref_files"]
    completed = subprocess.run(
        command, cwd=REPO_ROOT, capture_output=True, text=True, check=False
    )
    if completed.returncode != 0:
        sys.exit(
            f"{' '.join(command)} failed with exit code {completed.returncode}:\n{completed.stderr}"
        )
    results: dict[str, list[float]] = {}
    for line in completed.stdout.splitlines():
        name, _, values = line.partition("=")
        results[name] = [float(v) for v in values.split(",")]
    return results


def matches(expected: list[float] | None, actual: list[float] | None) -> bool:
    return (
        expected is not None
        and actual is not None
        and len(expected) == len(actual)
        and all(
            math.isclose(e, a, rel_tol=0.0, abs_tol=ABS_TOLERANCE)
            for e, a in zip(expected, actual)
        )
    )


def fmt(values: list[float] | None) -> str:
    return "<missing>" if values is None else ",".join(f"{v:.6g}" for v in values)


def main() -> int:
    expected = qiskit_results()
    if not expected:
        sys.exit(f"No .openqasm files found in {CIRCUITS_DIR}")
    actual = qsim_results()

    failures = 0
    for name in sorted(expected.keys() | actual.keys()):
        want, got = expected.get(name), actual.get(name)
        if matches(want, got):
            print(f"{PASS}  {name}")
        else:
            failures += 1
            print(
                f"{FAIL}  {name}\n        Qiskit: {fmt(want)}\n        qsim:   {fmt(got)}"
            )

    total = len(expected.keys() | actual.keys())
    print(f"\n{total - failures}/{total} circuits match Qiskit.")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
