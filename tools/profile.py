"""Print measured error and timing tables for the four demo programs."""

import argparse
from fractions import Fraction
from pathlib import Path
import re
import statistics
import subprocess

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--build-dir", type=Path, default=Path("build"))
parser.add_argument("--runs", type=int, default=21)
args = parser.parse_args()
if args.runs < 1:
    parser.error("--runs must be positive")

programs = ["double_rref", "double_lu", "double_lu_pivoted", "rational_rref"]
expected = Fraction(987, 8573324)
results = {}
for name in programs:
    executable = (args.build_dir / name).resolve()
    if not executable.exists():
        executable = executable.with_suffix(".exe")
    times = []
    reference = None
    for run in range(args.runs + 3):
        output = subprocess.check_output([str(executable)], text=True)
        entries = re.findall(r"^c([0-4])\s+(\S+)", output, re.MULTILINE)
        if len(entries) != 5:
            raise RuntimeError(f"Missing coefficients from {name}")
        values = [Fraction(value) if name == "rational_rref" else Fraction(float(value))
                  for _, value in entries]
        if reference is not None and values != reference:
            raise RuntimeError(f"Coefficient results changed between runs for {name}")
        reference = values
        if run >= 3:
            times.append(float(re.search(r"Solve time \(us\): (\S+)", output)[1]))
    errors = [abs(value - expected) for value in reference]
    residual = 0.0 if name == "rational_rref" else float(
        re.search(r"Relative equation residual:\s+(\S+)", output)[1])
    results[name] = (errors, statistics.median(times), residual)

print("| Coefficient | Double: row operations | Double: LU | Double: pivoted LU | Rational: row operations |")
print("| --- | ---: | ---: | ---: | ---: |")
for row in range(5):
    cells = ["0" if results[name][0][row] == 0 else f"{float(results[name][0][row]):.9e}"
             for name in programs]
    print(f"| `c{row}` | " + " | ".join(f"`{cell}`" for cell in cells) + " |")

print("\n| Method | Median solve time (microseconds) | Largest relative coefficient error | Relative equation residual |")
print("| --- | ---: | ---: | ---: |")
for name in programs:
    errors, median, residual = results[name]
    relative = float(max(errors) / abs(expected))
    print(f"| `{name}` | {median:.3f} | `{relative:.6e}` | `{residual:.6e}` |")
