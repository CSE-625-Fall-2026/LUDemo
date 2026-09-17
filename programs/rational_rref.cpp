#include <matrix/Matrix.hpp>
#include "ArbitraryInteger.hpp"
#include "RationalNumber.hpp"

#include <iostream>
#include <chrono>

using Rational = RationalNumber<ArbitraryInteger>;
using Matrix = matrix::Matrix<Rational>;

int main() {
    const int n = 5;
    const Rational expected(987, 8573324);
    Matrix a(n, n), b(n, 1);

    // The same five equations, constructed entirely with exact fractions.
    for (int row = 0; row < n; ++row) {
        const Rational x(1000 + row);
        Rational power(1);
        Rational sum(0);
        for (int col = 0; col < n; ++col) {
            a.set(row, col, power);
            sum += power;
            power *= x;
        }
        b.set(row, 0, Rational(987) * sum / Rational(8573324));
    }

    // Solve by row operations on [a | b].
    const auto start = std::chrono::steady_clock::now();
    const auto reduced = a.augment(b).rref();
    Matrix solution(n, 1);
    for (int row = 0; row < n; ++row) {
        solution.set(row, 0, reduced.get(row, n));
    }
    const double solveTime = std::chrono::duration<double, std::micro>(
        std::chrono::steady_clock::now() - start).count();

    std::cout << "RationalNumber<ArbitraryInteger> with row reduction\n"
              << "Expected every coefficient: " << expected << "\n\n"
              << "Coefficient   Computed value\n";
    for (int row = 0; row < n; ++row) {
        const Rational value = solution.get(row, 0);
        std::cout << "c" << row << "            " << value << '\n';
        if (value != expected) {
            std::cerr << "Exact coefficient check failed.\n";
            return 1;
        }
    }

    if (a * solution != b) {
        std::cerr << "Original equation check failed.\n";
        return 1;
    }
    std::cout << "\nAll five coefficients are exact; the equation residual is exactly zero.\n";
    std::cout << "Solve time (us): " << solveTime << '\n';
    return 0;
}
