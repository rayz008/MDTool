#ifndef TOOLS_H
#define TOOLS_H

#include <vector>
#include <string>
#include "system.h"

std::vector<double> savitzkyGolay(const std::vector<double>& y,
                                  int window_size, int order,
                                  int deriv, double rate);

int factorial(int n);

std::vector<double> smoothData(const std::vector<double>& data,
                               int window_size, int order);

#endif
