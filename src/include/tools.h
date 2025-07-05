#ifndef TOOLS_H
#define TOOLS_H
#include <vector>
#include <string>
#include "system.h"

/**
 * @brief Applies Savitzky Golay smoothing to input array
 *
 * @param[in,out] y The input vector for smoothing
 * @param[in] window_size
 * @param[in] order
 * @param[in] deriv
 * @param[in] rate
 */
std::vector<double> savitzkyGolay(const std::vector<double>& y,
                                  int window_size, int order,
                                  int deriv, double rate);

/**
 * @brief Computes n factorial
 */
int factorial(int n);

/**
 * @brief Applies smoothing to input array
 *
 * @param[in,out] data The input vector for smoothing
 * @param[in] window_size
 * @param[in] order
 */
std::vector<double> smoothData(const std::vector<double>& data,
                               int window_size, int order);

#endif
