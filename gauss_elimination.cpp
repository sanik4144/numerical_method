#include <iostream>
#include <vector>
#include <iomanip>

// Function to perform Gauss Elimination
std::vector<double> gaussElimination(const std::vector<std::vector<double>>& coefficients,
                                     const std::vector<double>& constants,
                                     int decimalPoints = 6) {
    int n = coefficients.size();
    std::vector<std::vector<double>> augmentedMatrix(n, std::vector<double>(n + 1, 0.0));

    // Construct the augmented matrix
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            augmentedMatrix[i][j] = coefficients[i][j];
        }
        augmentedMatrix[i][n] = constants[i];
    }

    std::cout << std::fixed << std::setprecision(decimalPoints);

    // Perform forward elimination
    for (int i = 0; i < n - 1; ++i) {
        for (int k = i + 1; k < n; ++k) {
            double factor = augmentedMatrix[k][i] / augmentedMatrix[i][i];
            for (int j = i; j < n + 1; ++j) {
                augmentedMatrix[k][j] -= factor * augmentedMatrix[i][j];
            }
        }
    }

    // Perform back-substitution
    std::vector<double> solution(n, 0.0);
    for (int i = n - 1; i >= 0; --i) {
        solution[i] = augmentedMatrix[i][n];
        for (int j = i + 1; j < n; ++j) {
            solution[i] -= augmentedMatrix[i][j] * solution[j];
        }
        solution[i] /= augmentedMatrix[i][i];
    }

    return solution;
}

int main() {
    // Define the coefficient matrix and constants vector
    std::vector<std::vector<double>> coefficients = {{{1,-1,2},
                                                     {1,2,3},
                                                     {3,-4,-5}}};
    std::vector<double> constants = {3,5,-13};

    // Solve the system of linear equations using Gauss Elimination
    int decimalPoints = 4; // Change this value to the desired decimal points
    std::vector<double> solution = gaussElimination(coefficients, constants, decimalPoints);

    // Display the solution
    std::cout << "Solution: ";
    for (double val : solution) {
        std::cout << val << " ";
    }
    std::cout << std::endl;

    return 0;
}
