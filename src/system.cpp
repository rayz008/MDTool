/**
 * @file system.cpp
 * @brief System constructed from trajectory and box information
 * 
 * Atomic system built from MD trajectory information files.
 * Currently only reads xyz trajectory.
 */

#include <fstream>
#include <sstream>
#include <stdexcept>
#include <vector>
#include <iostream>
#include <math.h>
#include <algorithm>
#include "system.h"

System::~System() {
    delete[] atoms;
    delete[] coords;
    delete[] boxes;
}

void System::allocateTrajectoryMemory() {
    if (traj_allocated) {
        delete[] atoms;
        delete[] coords;
    }

    // TODO: assumes each frame contains same atoms in same sequence
    atoms  = new std::string[natoms]; 
    coords = new double[nframes * natoms * 3];

    traj_allocated = true;
}

void System::allocateBoxMemory() {
    if (box_allocated) {
        delete[] boxes;
    }

    if (!traj_allocated) {
        throw std::runtime_error("Box memory allocated before trajectory memory allocation.");
    }

    if (fixed_volume) {
        boxes = new double[6];
    } else {
        boxes = new double[nframes * 6];
    }

    box_matrix  = new double[9];
    box_inverse = new double[9];

    box_allocated = true;
}

void System::readXYZ(const std::string &trajectory_file_name) {
    if (atoms || coords) {
        throw std::logic_error("Coordinates already in System instance.");
        return;
    }

    if (trajectory_file_name.empty()) {
        throw std::runtime_error("Trajectory file is not specified, please include an xyz file.");
    }
    
    std::ifstream file(trajectory_file_name);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open trajectory file, please check if it exists.");
    }

    std::cout << "Parsing trajectory file: " + trajectory_file_name << std::endl;

    // set natoms and nframes
    std::string line;
    std::istringstream iss;
    iss.clear();

    std::getline(file, line);
    iss.str(line);
    iss >> natoms;

    nframes++;
    while(std::getline(file, line)) {
        nframes++;
    }
    nframes /= (natoms + 2);

    file.clear();
    file.seekg(0);
    allocateTrajectoryMemory();

    // set atoms and coords
    for (int i = 0; i < nframes; i++) {
        std::getline(file, line);
        std::getline(file, line);

        for (int j = 0; j < natoms; j++) {
            std::istringstream tempiss;
            std::string dummy_atom;
            tempiss.clear();

            std::getline(file, line);
            tempiss.str(line);

            if (i == 0) {
                tempiss >> atoms[j];
            } else {
                tempiss >> dummy_atom;
            }
            
            tempiss.clear();
            tempiss.str(line);
            tempiss >> dummy_atom
                    >> coords[i * natoms * 3 + j * 3]
                    >> coords[i * natoms * 3 + j * 3 + 1]
                    >> coords[i * natoms * 3 + j * 3 + 2];

            if (dummy_atom != atoms[j] && i != 0) {
                std::cerr << "atomname" << dummy_atom << "  atoms[j]" << atoms[j] << std::endl;
                std::cerr << "Frame " << i << " has different atom name at index " << j << std::endl;
            }
        }
    }

    file.close();
    
    std::cout << "Trajectory file parsed successfully!" << std::endl;

}

void System::readBox(const Settings& settings) {
    bool box_read = false;
    
    // read from separate box file if specified
    if (!settings.box_infile.empty()) {
        try {
            std::cout << "Attempting to read box from file: " << settings.box_infile << std::endl;
            readBoxFromFile(settings.box_infile);
            box_read = true;
            std::cout << "Successfully read box information from separate file." << std::endl;
        } catch (const std::exception& e) {
            std::cerr << "Failed to read box from file '" << settings.box_infile 
                      << "': " << e.what() << std::endl;
        }
    }
    
    // read from XYZ trajectory file if no separate box file specified
    if (!box_read) {
        try {
            std::cout << "Attempting to read box from XYZ file: " << settings.traj_infile << std::endl;
            readBoxFromXYZ(settings.traj_infile);
            box_read = true;
            std::cout << "Successfully read box information from XYZ file." << std::endl;
        } catch (const std::exception& e) {
            std::cerr << "Failed to read box from XYZ file: " << e.what() << std::endl;
        }
    }
    
    if (!box_read) {
        std::cerr << "No box information found. Please verify box information in trajectory or in box file" << std::endl;
    }
    
    std::cout << "Box setup complete!" << std::endl;
}

void System::readBoxFromXYZ(const std::string &trajectory_file_name) {
    if (trajectory_file_name.empty()) {
        throw std::runtime_error("Trajectory file is not specified, please include an xyz file.");
    }
    
    std::ifstream file(trajectory_file_name);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open trajectory file, please check if it exists.");
    }

    std::vector<std::vector<double>> box_data;
    std::string line;

    for (int i = 0; i < nframes; i++) {
        getline(file, line);
        getline(file, line);
        std::istringstream iss(line);

        // box_format: -1: unknown, 3: orthogonal, 6: triclinic
        int box_format = -1;
        std::vector<double> temp_params;
        double value;

        while (iss >> value) {
            temp_params.push_back(value);
        }

        if (box_format == -1) {
            if (temp_params.size() == 3) {
                box_format = 3;
            } else if (temp_params.size() == 6) {
                box_format = 6;
            } else {
                throw std::runtime_error(std::string("Box can either take 3 parameters (a, b, c) for orthorhombic box")
                                        + " or 6 parameters (a, b, c, A, B, C) for triclinic box in a line.");
            }
        } else {
            if (static_cast<int>(temp_params.size()) != box_format) {
                throw std::runtime_error("Box format in consistent to the first frame!");
            }
        }

        std::vector<double> box_params(6);
        if (box_format == 3) {
            box_params[0] = temp_params[0];
            box_params[1] = temp_params[1];
            box_params[2] = temp_params[2];
            box_params[3] = 90;
            box_params[4] = 90;
            box_params[5] = 90;
        } else {
            box_params = temp_params;
        }

        //validateBoxParams(box_params); //TODO
        box_data.push_back(box_params);

        for (int j = 0; j < natoms; j++) {
            getline(file, line);
        }
    }

    file.close();

    if (box_data.empty()) {
        throw std::runtime_error("No valid box data found in file.");
    }

    fixed_volume = (box_data.size() == 1);
    allocateBoxMemory();

    if (fixed_volume) {
        for (int i = 0; i < 6; i++) {
            boxes[i] = box_data[0][i];
        }
    } else {
        if (static_cast<int>(box_data.size()) != nframes) {
            throw std::logic_error("Box entries not matching trajectory frame numbers");
        }

        for (int frame = 0; frame < nframes; frame++) {
            for (int i = 0; i < 6; i++) {
                boxes[frame * 6 + i] = box_data[frame][i];
            }
        }
    }

}

void System::readBoxFromFile(const std::string &box_file_name) {
    if (box_file_name.empty()) {
        throw std::runtime_error("Box file is not specified, please include a box file.");
    }
    
    std::ifstream file(box_file_name);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open box file, please check if it exists.");
    }

    std::vector<std::vector<double>> box_data;
    std::string line;
    // box_format: -1: unknown, 3: orthogonal, 6: triclinic
    int box_format = -1;

    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') {
            continue;
        }

        std::istringstream iss(line);
        std::vector<double> temp_params;
        double value;

        while (iss >> value) {
            temp_params.push_back(value);
        }

        if (box_format == -1) {
            if (temp_params.size() == 3) {
                box_format = 3;
            } else if (temp_params.size() == 6) {
                box_format = 6;
            } else {
                throw std::runtime_error(std::string("Box can either take 3 parameters (a, b, c) for orthorhombic box")
                                        + " or 6 parameters (a, b, c, A, B, C) for triclinic box in a line.");
            }
        } else {
            if (static_cast<int>(temp_params.size()) != box_format) {
                throw std::runtime_error("Box file format in consistent to the first line!");
            }
        }

        std::vector<double> box_params(6);
        if (box_format == 3) {
            box_params[0] = temp_params[0];
            box_params[1] = temp_params[1];
            box_params[2] = temp_params[2];
            box_params[3] = 90;
            box_params[4] = 90;
            box_params[5] = 90;
        } else {
            box_params = temp_params;
        }

        //validateBoxParams(box_params); //TODO
        box_data.push_back(box_params);

    }

    file.close();

    if (box_data.empty()) {
        throw std::runtime_error("No valid box data found in file.");
    }

    fixed_volume = (box_data.size() == 1);
    allocateBoxMemory();

    if (fixed_volume) {
        for (int i = 0; i < 6; i++) {
            boxes[i] = box_data[0][i];
        }
    } else {
        if (static_cast<int>(box_data.size()) != nframes) {
            throw std::logic_error("Box entries not matching trajectory frame numbers");
        }

        for (int frame = 0; frame < nframes; frame++) {
            for (int i = 0; i < 6; i++) {
                boxes[frame * 6 + i] = box_data[frame][i];
            }
        }
    }

}

void System::updateBoxInverse() {
    for (int i = 0; i < 3; i++) {
         for (int j = 0; j < 3; j++) {
             box_inverse[i * 3 + j] = ((box_matrix[3 * ((j + 1) % 3) + (i + 1) % 3]  * 
                                   box_matrix[3 * ((j + 2) % 3) + (i + 2) % 3]) -
                                  (box_matrix[3 * ((j + 1) % 3) + (i + 2) % 3]  *
                                   box_matrix[3 * ((j + 2) % 3) + (i + 1) % 3])) 
                                  / box_volume;
         }
    }

    return;
}

void System::updateBoxInformation(int frame) {
    double PI = 3.141592653589793238L;
    double radian_to_degree = PI / 180;

    // update box information only if volume not fixed or first frame if fixed volume
    if (!fixed_volume || frame == 0) {
        // update box_matrix
        box_matrix[0] = boxes[6 * frame];
        box_matrix[1] = boxes[6 * frame + 1] * cos(boxes[6 * frame + 5] * radian_to_degree);
        box_matrix[2] = boxes[6 * frame + 2] * cos(boxes[6 * frame + 4] * radian_to_degree);
        box_matrix[3] = 0;
        box_matrix[4] = boxes[6 * frame + 1] * sin(boxes[6 * frame + 5] * radian_to_degree);
        box_matrix[5] = (boxes[6 * frame + 1] * boxes[6 * frame + 2] * cos(boxes[6 * frame + 3]
                      * radian_to_degree) - box_matrix[1] * box_matrix[2]) / box_matrix[4];
        box_matrix[6] = 0;
        box_matrix[7] = 0;
        box_matrix[8] = sqrt(boxes[6 * frame + 2] * boxes[6 * frame + 2] 
                      - box_matrix[2] * box_matrix[2] - box_matrix[5] * box_matrix[5]);

        // update box_volume
        box_volume = 0;
        for (int i = 0; i < 3; i++) {
            box_volume += box_matrix[i] *
                         (box_matrix[3 + (i + 1) % 3] * box_matrix[6 + (i + 2) % 3] -
                          box_matrix[3 + (i + 2) % 3] * box_matrix[6 + (i + 1) % 3]);
        }

        if (box_volume <= 0) {
            throw std::logic_error("PBC box volume should be positive!");
        }

        // update box_inverse
        updateBoxInverse();
        
    }
}

