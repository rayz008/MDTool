#ifndef SYSTEM_H
#define SYSTEM_H

#include <string>
#include "settings.h"

// system stores the information of the atoms and coordinates as an array
struct System {
    bool traj_allocated = false;
    bool box_allocated = false;
    bool fixed_volume;

    int nframes{0};
    int natoms{0};
    double box_volume{0};
    // stores pointers to position of the atom array and the coord array
    std::string* atoms{nullptr};
    double* coords{nullptr};
    double* boxes{nullptr};
    double* box_matrix{nullptr};
    double* box_inverse{nullptr};
    
    System() = default;
    ~System();
    System(const System& other)              = delete;
    System& operator = (const System& other) = delete;
    System(System&& other)                   = delete;
    System& operator = (System&& other)      = delete;

    //initialize the system with the atom array and the coord array
    void allocateTrajectoryMemory();
    void allocateBoxMemory();
    void readXYZ(const std::string &filename);
    void readBoxFromXYZ(const std::string &trajectory_file_name);
    void readBoxFromFile(const std::string &box_file_name);
    void updateBoxInverse();
    void updateBoxInformation(int frame);
};

#endif
