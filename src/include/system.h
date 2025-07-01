#ifndef SYSTEM_H
#define SYSTEM_H

#include <string>
#include "settings.h"

// system stores the information of the atoms and coordinates as an array
struct System {
    int nframes{0};
    int natoms{0};
    double boxvol{0};
    // stores pointers to position of the atom array and the coord array
    std::string* atoms{nullptr};
    double* coords{nullptr};
    double* box{nullptr};
    double* boxmat{nullptr};
    double* boxinv{nullptr};
    
    System() = default;
    ~System();
    System(const System& other)              = delete;
    System& operator = (const System& other) = delete;
    System(System&& other)                   = delete;
    System& operator = (System&& other)      = delete;

    //initialize the system with the atom array and the coord array
    void AllocateMemory();
    void ReadXYZ(const std::string &filename);
    void UpdateNPTBox(int frame);
    void UpdateNVTBox(Settings& setting);
    void UpdateBoxInv();
};

#endif
