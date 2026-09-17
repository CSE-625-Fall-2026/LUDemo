# LUDemo

Three small programs solve the same five polynomial-fitting equations:

| Program | Number type | Method |
| --- | --- | --- |
| `double_rref` | `double` | Row operations with `a.augment(b).rref()` |
| `double_lu` | `double` | `a.lu()`, forward substitution, and back substitution |
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
These measured results use AppleClang 21 on macOS with the build commands below;
floating-point results can vary with compiler and platform. The table compares
the computed doubles against the exact fraction, with errors rounded for display.

| Coefficient | Double: row operations | Double: LU | Rational: row operations |
| --- | ---: | ---: | ---: |
| `c0` | `1.219680132e+03` | `1.219680113e+03` | `0` |
| `c1` | `4.867382290e+00` | `4.867382215e+00` | `0` |
| `c2` | `7.284099935e-03` | `7.284099824e-03` | `0` |
| `c3` | `4.844773495e-06` | `4.844773421e-06` | `0` |
| `c4` | `1.208375718e-09` | `1.208375699e-09` | `0` |

The true value of every coefficient is approximately `0.0001151245421262511`.
Both double methods produced approximately `1219.68` for `c0`, even though their
relative equation residuals, measured against the stored double equations,
were below `1e-15`. The rational method recovered
`987/8573324` for every coefficient, with exactly zero error and residual.

The double programs print relative coefficient errors and the relative residual
`max(abs(a*c - b)) / max(abs(b))`. A small residual does not guarantee accurate
coefficients in this system. The result reflects a poorly conditioned problem
and simple elimination; double arithmetic is not always inaccurate.
MatrixClassDemo's LU implementation does not pivot.

## Build and run

Install Git, CMake 3.16 or newer, and a C++17 compiler. CMake fetches
MatrixClassDemo `v3.0.0`, ArbitraryInteger `v4.0.0`, and RationalNumber `v5.0.0`.
RationalNumber is private: your GitHub account needs access, and Git must be able
to clone it using SSH. The first configuration needs internet access.

### Mac

Install Apple's command-line tools with `xcode-select --install` if needed.
From this repository:

```bash
cmake -S . -B build -G "Unix Makefiles"
cmake --build build
./build/double_rref
./build/double_lu
./build/rational_rref
ctest --test-dir build --output-on-failure
```

### Windows with Git Bash

Install MinGW-w64 and Ninja, and add them and CMake to `PATH`. Reopen Git Bash:

```bash
cmake -S . -B build -G Ninja -DCMAKE_CXX_COMPILER=g++
cmake --build build
./build/double_rref.exe
./build/double_lu.exe
./build/rational_rref.exe
ctest --test-dir build --output-on-failure
```

Each program is one source file in `programs/`. CTest runs all three: the double
programs check for finite answers and small equation residuals, while the
rational program requires exact coefficients and an exactly zero residual.
