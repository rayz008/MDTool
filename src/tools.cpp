#include <vector>
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <Eigen/Dense>
#include "system.h"
#include "tools.h"

std::vector<double> savitzkyGolay(const std::vector<double>& y,
                                  int window_size, int order,
                                  int deriv, double rate) {
    // Input validation
    if (window_size <= 0 || window_size % 2 == 0) {
        throw std::invalid_argument("window_size must be a positive odd number");
    }
    if (order < 0) {
        throw std::invalid_argument("order must be non-negative");
    }
    if (window_size < order + 2) {
        throw std::invalid_argument("window_size is too small for the polynomial order");
    }
    if (deriv < 0 || deriv > order) {
        throw std::invalid_argument("deriv must be between 0 and order");
    }
    if (y.size() < static_cast<size_t>(window_size)) {
        throw std::invalid_argument("Input data is too short for the specified window size");
    }

    int half_window = (window_size - 1) / 2;

    // Build Vandermonde matrix B
    // B[i][j] = k^j where k = i - half_window
    Eigen::MatrixXd B(window_size, order + 1);
    for (int i = 0; i < window_size; ++i) {
        int k = i - half_window;
        for (int j = 0; j <= order; ++j) {
            B(i, j) = std::pow(k, j);
        }
    }

    // Compute pseudoinverse and extract the derivative row
    Eigen::MatrixXd B_pinv = B.completeOrthogonalDecomposition().pseudoInverse();
    Eigen::VectorXd coeffs = B_pinv.row(deriv);

    // Apply derivative and rate scaling
    double deriv_factor = std::pow(rate, deriv) * factorial(deriv);
    coeffs *= deriv_factor;

    // Pad the signal at the extremes (matching Python implementation exactly)
    std::vector<double> y_padded;
    y_padded.reserve(y.size() + 2 * half_window);

    // First padding: firstvals = y[0] - abs(y[1:half_window+1][::-1] - y[0])
    for (int i = half_window; i >= 1; --i) {
        if (i < static_cast<int>(y.size())) {
            double val = y[0] - std::abs(y[i] - y[0]);
            y_padded.push_back(val);
        } else {
            y_padded.push_back(y[0]);
        }
    }

    // Original data
    y_padded.insert(y_padded.end(), y.begin(), y.end());

    // Last padding: lastvals = y[-1] + abs(y[-half_window-1:-1][::-1] - y[-1])
    size_t n = y.size();
    for (int i = 1; i <= half_window; ++i) {
        if (i < static_cast<int>(n)) {
            double val = y[n-1] + std::abs(y[n-1-i] - y[n-1]);
            y_padded.push_back(val);
        } else {
            y_padded.push_back(y[n-1]);
        }
    }

    // Apply convolution (equivalent to np.convolve(m[::-1], y, mode='valid'))
    std::vector<double> result;
    result.reserve(y.size());

    for (size_t i = half_window; i < y_padded.size() - half_window; ++i) {
        double sum = 0.0;
        for (int j = 0; j < window_size; ++j) {
            // Reverse the coefficients (m[::-1] in Python)
            sum += coeffs[window_size - 1 - j] * y_padded[i - half_window + j];
        }
        result.push_back(sum);
    }

    return result;
}

// Simple factorial function (already in your code, but included for completeness)
int factorial(int n) {
    if (n <= 1) return 1;
    return n * factorial(n - 1);
}

// Convenience function for just smoothing (deriv=0, rate=1)
std::vector<double> smoothData(const std::vector<double>& data,
                               int window_size, int order) {
    return savitzkyGolay(data, window_size, order, 0, 1.0);
}
