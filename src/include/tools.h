#ifndef TOOLS_H
#define TOOLS_H

#include <vector>
#include <string>
#include "system.h"

std::vector<std::pair<int, int>> countPairs(
    const System& sys, const std::string &atomA, const std::string &atomB,
    int &numA, int &numB);

#endif
