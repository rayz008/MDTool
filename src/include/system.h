#ifndef SYSTEM_H
#define SYSTEM_H
#include <string>
#include "settings.h"

/**
 * @struct System
 * @brief System constructed from trajectory and box information
 * 
 * The System struct stores atomic and periodic boundary condition info from trajectory
 * and box input files stored in JSON Setting. 
 *
 * @note Currently only handles xyz input. Can read box information from xyz or separate file.
 */
struct System {
    bool traj_allocated = false;
    bool box_allocated = false;
    bool fixed_volume;                  // if pbc is fixed

    int nframes{0};
    int natoms{0};
    double box_volume{0};               // Volume of box at requested frame
    std::string* atoms{nullptr};
    double* coords{nullptr};
    double* boxes{nullptr};
    double* box_matrix{nullptr};        // Matrix of box at requested frame
    double* box_inverse{nullptr};       // Inverse matrix of box at requested frame
    
    System() = default;
    ~System();
    System(const System& other)              = delete;
    System& operator = (const System& other) = delete;
    System(System&& other)                   = delete;
    System& operator = (System&& other)      = delete;

    /**
     * @brief Allocates memory for trajectory information
     */
    void allocateTrajectoryMemory();

    /**
     * @brief Allocates memory for box information
     */
    void allocateBoxMemory();

    /**
     * @brief Reads atoms and coordinate information from filename
     *
     * @param[in] filename The xyz trajectory file name. 
     */
    void readXYZ(const std::string &filename);

    /**
     * @brief Reads box information from Settings parameter
     *
     * @param[in] settings The JSON setting information. 
     */
    void readBox(const Settings& settings);

    /**
     * @brief Reads box information from comment lines in xyz trajectory file
     */
    void readBoxFromXYZ(const std::string &trajectory_file_name);

    /**
     * @brief Reads box information from a separate box file
     */
    void readBoxFromFile(const std::string &box_file_name);

    /**
     * @brief Calculates inverse matrix of the periodic boundary box
     */
    void updateBoxInverse();

    /**
     * @brief Updates box information in current frame.
     *
     * Updates the box matrix, box volume and box inverse matrix for current frame
     * @param[in] frame Current frame index for box information update
     */
    void updateBoxInformation(int frame);
};

#endif
