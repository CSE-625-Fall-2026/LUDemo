#include <matrix/Matrix.hpp>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>

using Matrix = matrix::Matrix<double>;

int main() {
    const int n = 5;
    const double expected = 987.0 / 8573324.0;
    Matrix a(n, n), b(n, 1);

    // Fit c0 + c1*x + c2*x^2 + c3*x^3 + c4*x^4 at x = 1000,...,1004.
    // In exact arithmetic, every coefficient is 987/8573324.
    for (int row = 0; row < n; ++row) {
        const double x = 1000 + row;
        double power = 1;
        double sum = 0;
        for (int col = 0; col < n; ++col) {
            a.set(row, col, power);
            sum += power;
            power *= x;
        }
        b.set(row, 0, (987.0 * sum) / 8573324.0);
    }

    // Solve by row operations on [a | b].
    const auto reduced = a.augment(b).rref();
    Matrix solution(n, 1);
    for (int row = 0; row < n; ++row) {
        solution.set(row, 0, reduced.get(row, n));
        for (int col = 0; col < n; ++col) {
            if (reduced.get(row, col) != (row == col ? 1.0 : 0.0)) {
                std::cerr << "Row reduction did not produce an identity block.\n";
                return 1;
            }
        }
    }

    std::cout << std::scientific << std::setprecision(16)
              << "double with row reduction\n"
              << "Expected every coefficient: 987/8573324 = " << expected << "\n\n"
              << "Coefficient     Computed value          Relative error\n";
    double maxError = 0;
    for (int row = 0; row < n; ++row) {
        const double value = solution.get(row, 0);
        if (!std::isfinite(value)) {
            std::cerr << "Non-finite coefficient.\n";
            return 1;
        }
        const double error = std::abs((value - expected) / expected);
        maxError = std::max(maxError, error);
        std::cout << "c" << row << "              " << value << "      " << error << '\n';
    }

    // A small equation residual does not guarantee accurate coefficients.
    const auto fitted = a * solution;
    double maxResidual = 0, maxB = 0;
    for (int row = 0; row < n; ++row) {
        maxResidual = std::max(maxResidual, std::abs(fitted.get(row, 0) - b.get(row, 0)));
        maxB = std::max(maxB, std::abs(b.get(row, 0)));
    }
    const double relativeResidual = maxResidual / maxB;
    std::cout << "\nLargest relative coefficient error: " << maxError
              << "\nRelative equation residual:          " << relativeResidual << '\n';
    return std::isfinite(relativeResidual) && relativeResidual < 1e-10 ? 0 : 1;
}
