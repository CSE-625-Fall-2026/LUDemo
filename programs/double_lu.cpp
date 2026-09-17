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

    // The same five polynomial equations as in the row-reduction programs.
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

    // a = L*U. First solve L*y = b, then U*c = y.
    const auto [L, U] = a.lu();
    Matrix y(n, 1), solution(n, 1);

    // Forward substitution.
    for (int row = 0; row < n; ++row) {
        double value = b.get(row, 0);
        for (int col = 0; col < row; ++col) {
            value -= L.get(row, col) * y.get(col, 0);
        }
        y.set(row, 0, value / L.get(row, row));
    }

    // Back substitution.
    for (int row = n - 1; row >= 0; --row) {
        double value = y.get(row, 0);
        for (int col = row + 1; col < n; ++col) {
            value -= U.get(row, col) * solution.get(col, 0);
        }
        solution.set(row, 0, value / U.get(row, row));
    }

    std::cout << std::scientific << std::setprecision(16)
              << "double with LU decomposition\n"
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
