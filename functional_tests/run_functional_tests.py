#!/usr/bin/env python3
from __future__ import annotations

import difflib
import subprocess
import sys
from pathlib import Path


def remove_empty_lines(text: str) -> str:
    return "\n".join(line for line in text.splitlines() if line != "")


def run_command(command: list[str], cwd: Path) -> str:
    completed = subprocess.run(
        command,
        cwd=cwd,
        check=True,
        capture_output=True,
        text=True,
    )
    return remove_empty_lines(completed.stdout)


def main() -> int:
    script_dir = Path(__file__).resolve().parent
    repo_root = script_dir.parent
    python_bin = script_dir / ".venv" / "bin" / "python"
    python_script = script_dir / "qiskit_tests.py"

    if not python_bin.exists() or not python_bin.is_file() or not python_bin.stat().st_mode & 0o111:
        print("ERROR: qiskit_tests virtual environment is missing or not executable.", file=sys.stderr)
        print("Create it first and install qiskit inside qiskit_tests/.venv.", file=sys.stderr)
        return 1

    qsim_lines = run_command(["cargo", "run", "--quiet", "--example", "qsim_tests"], repo_root).splitlines()
    qiskit_lines = run_command([str(python_bin), str(python_script)], repo_root).splitlines()

    all_passed = True
    print("Test results (Qiskit vs qsim_statevec_cpu):\n")
    GREEN = "\033[92m"
    RED = "\033[91m"
    BOLD = "\033[1m"
    RESET = "\033[0m"
    for qiskit_line, qsim_line in zip(qiskit_lines, qsim_lines):
        # Both lines are in the format: NAME=RESULTS
        qiskit_name, qiskit_result = qiskit_line.split("=", 1)
        _, qsim_result = qsim_line.split("=", 1)
        print(f"{BOLD}{qiskit_name}{RESET}")
        print(f"  Qiskit:            {qiskit_result}")
        print(f"  qsim_statevec_cpu: {qsim_result}")
        if qiskit_result == qsim_result:
            print(f"  Result:           {GREEN}PASS{RESET}")
        else:
            print(f"  Result:           {RED}FAIL{RESET}")
            all_passed = False

    if all_passed:
        print("\nAll tests passed: qsim_statevec_cpu output matches Qiskit for all reference circuits.")
        return 0
    else:
        print("\nSome tests failed: see above for mismatches.", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())