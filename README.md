# LUDemo

Four small programs solve the same five polynomial-fitting equations:

| Program | Number type | Method |
| --- | --- | --- |
| `double_rref` | `double` | Row operations with `a.augment(b).rref()` |
| `double_lu` | `double` | `a.lu()`, forward substitution, and back substitution |
| `double_lu_pivoted` | `double` | `a.luPartialPivoting()`, then substitution using `P*b` |
| `rational_rref` | `RationalNumber<ArbitraryInteger>` | Row operations with `a.augment(b).rref()` |

## The system

```text
p(x) = c0 + c1*x + c2*x^2 + c3*x^3 + c4*x^4
x    = 1000, 1001, 1002, 1003, 1004
q    = 987/8573324
p(x) = q * (1 + x + x^2 + x^3 + x^4)
```

Each matrix row is `[1, x, x^2, x^3, x^4]`; the right-hand side is `p(x)`.
The exact answer is `c0 = c1 = c2 = c3 = c4 = q`. The five distinct sample
locations make the solution unique. This is degree-four interpolation, with
five equations for five unknowns and no added measurement noise.

The large, nearby sample locations make this Vandermonde system poorly
conditioned: tiny changes in the sample values can cause large changes in the
coefficients. The double programs round the input fractions and introduce
further rounding during elimination. The rational program constructs and solves
the system exactly. ArbitraryInteger prevents fixed-size integer overflow.
See the [NumPy polynomial-fitting notes](https://numpy.org/doc/stable/reference/generated/numpy.polyfit.html#notes)
for why poorly centered sample locations make coefficient recovery sensitive.

## Rounding-error comparison

Absolute coefficient error is `abs(computed coefficient - 987/8573324)`.
These measured results use AppleClang 21, an Apple M4, and macOS 26.6.2 with a Release build;
floating-point results can vary with compiler and platform. The table compares
the computed doubles against the exact fraction, with errors rounded for display.

| Coefficient | Double: row operations | Double: LU | Double: pivoted LU | Rational: row operations |
| --- | ---: | ---: | ---: | ---: |
| `c0` | `1.219680132e+03` | `1.219680113e+03` | `1.219680113e+03` | `0` |
| `c1` | `4.867382290e+00` | `4.867382215e+00` | `4.867382215e+00` | `0` |
| `c2` | `7.284099935e-03` | `7.284099824e-03` | `7.284099824e-03` | `0` |
| `c3` | `4.844773495e-06` | `4.844773421e-06` | `4.844773421e-06` | `0` |
| `c4` | `1.208375718e-09` | `1.208375699e-09` | `1.208375699e-09` | `0` |

The true value of every coefficient is approximately `0.0001151245421262511`.
All three double methods produced approximately `1219.68` for `c0`, even though their
relative equation residuals, measured against the stored double equations,
were below `1e-15`. The rational method recovered
`987/8573324` for every coefficient, with exactly zero error and residual.

The double programs print relative coefficient errors and the relative residual
`max(abs(a*c - b)) / max(abs(b))`. A small residual does not guarantee accurate
coefficients in this system. The result reflects a poorly conditioned problem
and simple elimination; double arithmetic is not always inaccurate.
Partial pivoting improves the choice of pivots, but does not materially improve
the coefficients in this example. The system is highly sensitive to the input
rounding that has already occurred. At the displayed precision, plain and
pivoted LU have the same errors.

## Partial pivoting

```cpp
const auto [P, L, U] = a.luPartialPivoting(); // P*a = L*U
const auto permutedB = P * b;
```

At each elimination step, the library chooses the largest absolute value in
the remaining column and swaps that row into place. The demo then solves
`L*y = P*b`, followed by `U*c = y`, using the same simple substitution loops as
the original LU program. `double_lu` retains the original method for comparison.

## Solve-time profile

Measured on the same Release build as the error table. Each median covers 21
process runs after three discarded warm-up runs. An internal steady clock times
only the solve and result extraction, excluding input construction, correctness
checks, console output, and process startup. These tiny-system timings are
approximate and platform-dependent.

| Method | Median solve time (microseconds) | Largest relative coefficient error | Relative equation residual |
| --- | ---: | ---: | ---: |
| `double_rref` | 1.208 | `1.059444e+07` | `3.817735e-16` |
| `double_lu` | 0.917 | `1.059444e+07` | `1.272578e-16` |
| `double_lu_pivoted` | 1.375 | `1.059444e+07` | `1.272578e-16` |
| `rational_rref` | 3245.000 | `0.000000e+00` | `0.000000e+00` |

After building, regenerate both Markdown tables with Python 3:

```bash
python3 tools/profile.py --build-dir build --runs 21
```

## Build and run

Install Git, CMake 3.16 or newer, and a C++17 compiler. CMake fetches
MatrixClassDemo `v3.1.0`, ArbitraryInteger `v4.0.0`, and RationalNumber `v5.0.0`.
RationalNumber is private: your GitHub account needs access, and Git must be able
to clone it using SSH. The first configuration needs internet access.

### Mac

Install Apple's command-line tools with `xcode-select --install` if needed.
From this repository:

```bash
cmake -S . -B build -G "Unix Makefiles" -DCMAKE_BUILD_TYPE=Release
cmake --build build
./build/double_rref
./build/double_lu
./build/double_lu_pivoted
./build/rational_rref
ctest --test-dir build --output-on-failure
```

### Windows with Git Bash

Install MinGW-w64 and Ninja, and add them and CMake to `PATH`. Reopen Git Bash:

```bash
cmake -S . -B build -G Ninja -DCMAKE_CXX_COMPILER=g++ -DCMAKE_BUILD_TYPE=Release
cmake --build build
./build/double_rref.exe
./build/double_lu.exe
./build/double_lu_pivoted.exe
./build/rational_rref.exe
ctest --test-dir build --output-on-failure
```

Each program is one source file in `programs/`. CTest runs all four: the double
programs check for finite answers and small equation residuals, while the
rational program requires exact coefficients and an exactly zero residual.

## Versions

| LUDemo | MatrixClassDemo | Programs |
| --- | --- | --- |
| `v0.1.0` | `v3.0.0` | Double row reduction, double LU, and rational row reduction |
| `v0.2.0` | `v3.1.0` | Adds double LU with partial pivoting and reproducible error/time tables |
