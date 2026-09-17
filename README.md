# LUDemo

Five small C++ programs compare ways to solve matrix equations.
Everything is kept on `main`.

| Program | Method |
| --- | --- |
| `double_rref` | Row operations with `double` |
| `double_lu` | LU decomposition with `double` |
| `double_lu_pivoted` | LU with partial pivoting and `double` |
| `rational_rref` | Exact row operations with `RationalNumber<ArbitraryInteger>` |
| `sensor_calibration` | Plain LU and pivoted LU on a second 5×5 system |

## Polynomial example

The first four programs fit this polynomial using five equations:

```text
p(x) = c0 + c1*x + c2*x^2 + c3*x^3 + c4*x^4
x = 1000, 1001, 1002, 1003, 1004
q = 987/8573324
p(x) = q * (1 + x + x^2 + x^3 + x^4)
```

The correct answer is `q` for all five coefficients, about `0.0001151245421262511`.
The large, closely spaced x values make this problem sensitive to rounding.
The rational program uses exact fractions throughout the calculation.

This table shows absolute error, or `abs(computed value - q)`.

| Coefficient | Double: row operations | Double: LU | Double: pivoted LU | Rational: row operations |
| --- | ---: | ---: | ---: | ---: |
| `c0` | `1.219680132e+03` | `1.219680113e+03` | `1.219680113e+03` | `0` |
| `c1` | `4.867382290e+00` | `4.867382215e+00` | `4.867382215e+00` | `0` |
| `c2` | `7.284099935e-03` | `7.284099824e-03` | `7.284099824e-03` | `0` |
| `c3` | `4.844773495e-06` | `4.844773421e-06` | `4.844773421e-06` | `0` |
| `c4` | `1.208375718e-09` | `1.208375699e-09` | `1.208375699e-09` | `0` |

All three double methods returned about `1219.68` for `c0`.
Partial pivoting does not fix the input rounding in this example.
The rational method returned the exact answer for every coefficient.

## Sensor example

Five measurement channels mix five unknown signal voltages.
We solve `A*x = b` to recover those voltages.
This is a made-up teaching example with realistic values and no measurement noise.

```text
      [ 0.1  0.3  0.2  0.1  0.4 ]        [ 3.7 ]
      [ 0.3  0.9  0.1  0.2  0.1 ]        [ 3.7 ]
A  =  [ 0.5  0.2  0.8  0.1  0.3 ]   b =  [ 5.2 ] volts
      [ 0.2  0.7  0.3  0.9  0.1 ]        [ 6.6 ]
      [ 0.4  0.1  0.2  0.3  0.8 ]        [ 6.4 ]

Correct answer: x = [1, 2, 3, 4, 5] volts
```

Both methods use the same inputs. Values and errors below are in volts.

| Signal | Expected | Plain LU | Absolute error | Pivoted LU | Absolute error |
| --- | ---: | ---: | ---: | ---: | ---: |
| `x1` | 1 | -6.8473677618669901 | `7.847368e+00` | 0.99999999999999922 | `7.771561e-16` |
| `x2` | 2 | 4.6314031645305391 | `2.631403e+00` | 2.0000000000000000 | `0` |
| `x3` | 3 | 3.5471965933286032 | `5.471966e-01` | 3.0000000000000000 | `0` |
| `x4` | 4 | 3.7700496806245560 | `2.299503e-01` | 3.9999999999999996 | `4.440892e-16` |
| `x5` | 5 | 4.7721788502484026 | `2.278211e-01` | 5.0000000000000000 | `0` |

Plain LU divides by a tiny second pivot, about `1.89e-16`, and loses accuracy.
Partial pivoting swaps rows to choose larger pivots and gives an accurate answer.
On other platforms, plain LU may hit a zero pivot and report that a row swap is needed.

## Timing

These results were measured with AppleClang 21 on an Apple M4 running macOS
26.6.2, using a Release build. Results and timings can vary by platform.

Times below are medians of 21 runs after three warm-up runs.
Only the solve is timed, including result extraction. Input setup and printing are excluded.
The residual measures how closely the answer satisfies the stored equations.
A small residual does not guarantee accurate coefficients in the polynomial example.

| Program | Solve time (microseconds) | Largest relative coefficient error | Relative equation residual |
| --- | ---: | ---: | ---: |
| `double_rref` | 1.208 | `1.059444e+07` | `3.817735e-16` |
| `double_lu` | 0.917 | `1.059444e+07` | `1.272578e-16` |
| `double_lu_pivoted` | 1.375 | `1.059444e+07` | `1.272578e-16` |
| `rational_rref` | 3245.000 | `0` | `0` |

## Build and run

You need Git, CMake 3.16 or newer, and a C++17 compiler.
CMake downloads MatrixClassDemo, ArbitraryInteger, and RationalNumber.
The first build needs internet access and GitHub SSH access to the private RationalNumber repository.

### Mac

Install Apple's command-line tools with `xcode-select --install` if needed.
Run these commands from the repository folder:

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

Install MinGW-w64 and Ninja. Add them and CMake to `PATH`, then reopen Git Bash.

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

Each program is one file in `programs/`. CTest runs all five.
The rational program checks for exact answers. The double polynomial programs
check equation residuals. The sensor program checks the pivoted answer.
