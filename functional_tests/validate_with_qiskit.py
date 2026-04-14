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

    qsim_output = run_command(["cargo", "run", "--quiet", "--example", "qsim_tests"], repo_root)
    qiskit_output = run_command([str(python_bin), str(python_script)], repo_root)

    print("qsim_statevec_cpu output:")
    print(qsim_output)
    print()
    print("Qiskit output:")
    print(qiskit_output)

    if qsim_output != qiskit_output:
        print()
        print("FAIL: Rust simulator output does not match Qiskit.", file=sys.stderr)
        qsim_lines = qsim_output.splitlines()
        qiskit_lines = qiskit_output.splitlines()
        diff = difflib.unified_diff(qsim_lines, qiskit_lines, fromfile="qsim_output", tofile="qiskit_output", lineterm="")
        for line in diff:
            print(line, file=sys.stderr)
        return 1

    print()
    print("PASS: qsim_statevec_cpu simulator output matches Qiskit for all reference circuits.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())