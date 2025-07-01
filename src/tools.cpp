#include <vector>
#include "system.h"

std::vector<std::pair<int, int>> calculate_couples(
    const System& sys,
    const std::string &atomA,
    const std::string &atomB,
    int &numA,
    int &numB)
{
    std::vector<std::pair<int, int>> vec;

    for(int a{0}; a < sys.natoms; a++) {
        if(sys.atoms[a] == atomA) {
            numA++;
        }
        else if(sys.atoms[a] == atomB) {
            numB++;
        }
    }

    if(atomA == atomB) {
        numB = numA;
    }

    vec.reserve(numA * numB);
    // reduce
    for(int a{0}; a < sys.natoms; a++) {
        if(sys.atoms[a] == atomA) {
            for(int b{0}; b < sys.natoms; b++) {
                if(sys.atoms[b] == atomB) {
                    vec.push_back(std::make_pair(a, b));
                }
            }
        }
    }
    return vec;
}
