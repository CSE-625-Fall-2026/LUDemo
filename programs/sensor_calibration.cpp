#include <matrix/Matrix.hpp>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>

using Matrix = matrix::Matrix<double>;

Matrix substitute(const Matrix& L, const Matrix& U, const Matrix& b) {
    const int n = b.getNumRows();
    Matrix y(n, 1), x(n, 1);

    // Forward substitution: L*y = b.
    for (int row = 0; row < n; ++row) {
        double value = b.get(row, 0);
        for (int col = 0; col < row; ++col) {
            value -= L.get(row, col) * y.get(col, 0);
        }
        y.set(row, 0, value / L.get(row, row));
    }

    // Back substitution: U*x = y.
    for (int row = n - 1; row >= 0; --row) {
        if (U.get(row, row) == 0) {
            throw std::domain_error("Back substitution requires a nonzero pivot");
        }
        double value = y.get(row, 0);
        for (int col = row + 1; col < n; ++col) {
            value -= U.get(row, col) * x.get(col, 0);
        }
        x.set(row, 0, value / U.get(row, row));
    }
    return x;
}

bool report(const Matrix& a, const Matrix& b, const Matrix& x) {
    const auto fitted = a * x;
    double maxError = 0, maxResidual = 0;
    bool finite = true;
    std::cout << "Signal   Computed voltage       Absolute error (V)\n";
    for (int row = 0; row < x.getNumRows(); ++row) {
        const double value = x.get(row, 0);
        const double error = std::abs(value - (row + 1));
        const double residual = std::abs(fitted.get(row, 0) - b.get(row, 0));
        finite = finite && std::isfinite(value) && std::isfinite(residual);
        maxError = std::max(maxError, error);
        maxResidual = std::max(maxResidual, residual);
        std::cout << "x" << row + 1 << "       " << value << "   " << error << '\n';
    }
    std::cout << "Largest absolute solution error (V): " << maxError
              << "\nLargest absolute equation residual (V): " << maxResidual << "\n\n";
    return finite && maxError < 1e-10 && maxResidual < 1e-10;
}

int main() {
    // Synthetic sensor calibration: dimensionless gains and voltages in volts.
    const double gains[5][5] = {
        {0.1, 0.3, 0.2, 0.1, 0.4},
        {0.3, 0.9, 0.1, 0.2, 0.1},
        {0.5, 0.2, 0.8, 0.1, 0.3},
        {0.2, 0.7, 0.3, 0.9, 0.1},
        {0.4, 0.1, 0.2, 0.3, 0.8}
    };
    const double readings[5] = {3.7, 3.7, 5.2, 6.6, 6.4};
    Matrix a(5, 5), b(5, 1);
    for (int row = 0; row < 5; ++row) {
        for (int col = 0; col < 5; ++col) {
            a.set(row, col, gains[row][col]);
        }
        b.set(row, 0, readings[row]);
    }

    std::cout << std::scientific << std::setprecision(16)
              << "Expected voltages: 1, 2, 3, 4, 5 V\n\n";
    // Depending on rounding, plain LU either loses accuracy or hits a zero pivot.
    try {
        const auto [L, U] = a.lu();
        std::cout << "Plain LU\nSecond pivot: " << U.get(1, 1) << '\n';
        report(a, b, substitute(L, U, b));
    } catch (const std::domain_error& error) {
        std::cout << "Plain LU could not solve the system: " << error.what() << "\n\n";
    }

    const auto [P, L, U] = a.luPartialPivoting();
    std::cout << "LU with partial pivoting\nSecond pivot: " << U.get(1, 1) << '\n';
    const auto x = substitute(L, U, P * b);
    // Check the pivoted answer; the inaccurate plain result is the demonstration.
    return report(a, b, x) ? 0 : 1;
}
