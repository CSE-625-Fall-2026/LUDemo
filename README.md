# LUDemo

Small C++ programs compare row operations and LU decomposition. LUDemo is
maintained directly on `main`, without numbered releases.

The first four programs solve the same five polynomial-fitting equations.
The fifth demonstrates partial pivoting with a sensor-calibration system.

| Program | Number type | Method |
| --- | --- | --- |
| `double_rref` | `double` | Row operations with `a.augment(b).rref()` |
| `double_lu` | `double` | `a.lu()`, forward substitution, and back substitution |
| `double_lu_pivoted` | `double` | `a.luPartialPivoting()`, then substitution using `P*b` |
| `rational_rref` | `RationalNumber<ArbitraryInteger>` | Row operations with `a.augment(b).rref()` |
| `sensor_calibration` | `double` | Plain LU and partial-pivoting LU on the same 5×5 system |

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

## A 5×5 sensor-calibration example

Suppose five measurement channels each pick up a mixture of five unknown
signal voltages. The calibration matrix `A` contains dimensionless gains;
`b` contains the measured voltages. We solve `A*x = b` to recover the signals.
These gains and voltages are plausible for a small analog measurement system:
the readings fit within ordinary [±10 V measurement ranges](https://www.ni.com/en/shop/hardware-portfolio/daq-sensor-io/voltage.html).
This is a deliberately constructed teaching example, not recorded sensor data.
No measurement noise is added, so we can isolate the arithmetic error.

```text
      [ 0.1  0.3  0.2  0.1  0.4 ]        [ 3.7 ]
      [ 0.3  0.9  0.1  0.2  0.1 ]        [ 3.7 ]
A  =  [ 0.5  0.2  0.8  0.1  0.3 ]   b =  [ 5.2 ] volts
      [ 0.2  0.7  0.3  0.9  0.1 ]        [ 6.6 ]
      [ 0.4  0.1  0.2  0.3  0.8 ]        [ 6.4 ]

Exact solution for these decimal values: x = [1, 2, 3, 4, 5] volts
```

The `sensor_calibration` program passes the same `double` inputs to both
methods and uses the same forward- and back-substitution function. This run
used AppleClang 21 on an Apple M4, macOS 26.6.2, with a Release build.
All solution values and absolute errors below are in volts.

| Signal | Expected | Plain LU | Absolute error | LU with partial pivoting | Absolute error |
| --- | ---: | ---: | ---: | ---: | ---: |
| `x1` | 1 | -6.8473677618669901 | `7.847368e+00` | 0.99999999999999922 | `7.771561e-16` |
| `x2` | 2 | 4.6314031645305391 | `2.631403e+00` | 2.0000000000000000 | `0` |
| `x3` | 3 | 3.5471965933286032 | `5.471966e-01` | 3.0000000000000000 | `0` |
| `x4` | 4 | 3.7700496806245560 | `2.299503e-01` | 3.9999999999999996 | `4.440892e-16` |
| `x5` | 5 | 4.7721788502484026 | `2.278211e-01` | 5.0000000000000000 | `0` |

| Diagnostic | Plain LU | LU with partial pivoting |
| --- | ---: | ---: |
| Second pivot | `1.887379e-16` | `0.78` |
| Largest absolute solution error (V) | `7.847368e+00` | `7.771561e-16` |
| Largest absolute equation residual, `max(abs(A*x - b))` (V) | `3.050987e+00` | `8.881784e-16` |

The first two rows begin with proportional pairs: `[0.1, 0.3]` and
`[0.3, 0.9]`. Without a row swap, the second pivot would be zero in exact
arithmetic. In this run, rounding leaves a tiny nonzero pivot, and dividing by
it produces huge multipliers that amplify arithmetic errors. The full 5×5
matrix still has a unique solution. Partial pivoting chooses better rows and
recovers the known voltages to nearly machine precision.

Results vary with compiler and platform; plain LU may instead encounter an
exactly zero pivot and report that a row swap is required. Run
`./build/sensor_calibration` to see both solutions on your machine.

## Polynomial solve-time profile

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

Each polynomial program prints its own solve time and error measurements.

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
./build/sensor_calibration
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
./build/sensor_calibration.exe
ctest --test-dir build --output-on-failure
```

Each program is one source file in `programs/`. CTest runs all five. The double
polynomial programs check for finite answers and small equation residuals;
the rational program requires exact coefficients and an exactly zero residual.
The sensor example checks that the pivoted solution and equation residual are
accurate within `1e-10`. Its plain LU result is printed for comparison.
