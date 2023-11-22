#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

// Function to perform Gauss-Seidel iteration
void gaussSeidelIteration(const std::vector<std::vector<double>>& coefficients,
                          const std::vector<double>& constants,
                          std::vector<double>& initialGuess,
                          double tol = 1e-10,
                          int maxIterations = 1000,
                          int decimalPoints = 6) {
    int n = coefficients.size();
    std::vector<double> x = initialGuess;
    std::vector<double> xNew(n, 0.0);

    std::cout << std::fixed << std::setprecision(decimalPoints);

    for (int iteration = 0; iteration < maxIterations; ++iteration) {
        double error = 0.0;

        for (int i = 0; i < n; ++i) {
            double sum_val = 0.0;

            for (int j = 0; j < n; ++j) {
                if (j != i) {
                    sum_val += coefficients[i][j] * x[j];
                }
            }

            xNew[i] = (constants[i] - sum_val) / coefficients[i][i];
            error += std::abs(xNew[i] - x[i]);
            x[i] = xNew[i];
        }

        std::cout << "Iteration " << iteration + 1 << ": ";
        for (double val : x) {
            std::cout << val << " ";
        }
        std::cout << std::endl;

        if (error < tol) {
            std::cout << "Converged after " << iteration + 1 << " iterations\n";
            return;
        }
    }

    std::cout << "Did not converge within the specified tolerance and maximum iterations.\n";
}

int main() {
    // Define the coefficient matrix and constants vector
    std::vector<std::vector<double>> coefficients = {{{27,6,-1},
                                                     {6, 15, 2},
                                                     {1, 1, 54}}};
    std::vector<double> constants = {85,72,110};

    // Initial guess
    std::vector<double> initialGuess = {0, 0, 0};

    // Solve the system of linear equations using Gauss-Seidel iteration
    int decimalPoints = 3; // Change this value to the desired decimal points
    gaussSeidelIteration(coefficients, constants, initialGuess, 1e-10, 1000, decimalPoints);

    return 0;
}
