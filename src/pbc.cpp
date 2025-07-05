/**
 * @file pbc.cpp
 * @brief Minimum image conversion to distance vectors based on periodic 
 * boundary conditions given orthorhombic or triclinic boxes.
 */

#include <math.h>
#include "system.h"

void pbcOrthorhombic(double& dx, double& dy, double& dz, System& sys) {
    dx -= rint(dx / sys.box_matrix[0]) * sys.box_matrix[0];
    dy -= rint(dy / sys.box_matrix[4]) * sys.box_matrix[4];
    dz -= rint(dz / sys.box_matrix[8]) * sys.box_matrix[8];
}

void pbcTriclinic(double& dx, double& dy, double& dz, System& sys) {

    double ds[3] = {0, 0, 0};
    for (int i = 0; i < 3; i++) {
        ds[i] = rint(sys.box_inverse[i * 3]     * dx
                   + sys.box_inverse[i * 3 + 1] * dy
                   + sys.box_inverse[i * 3 + 2] * dz);
    }
    dx -= sys.box_matrix[0] * ds[0] + sys.box_matrix[1] * ds[1] + sys.box_matrix[2] * ds[2];
    dy -= sys.box_matrix[3] * ds[0] + sys.box_matrix[4] * ds[1] + sys.box_matrix[5] * ds[2];
    dz -= sys.box_matrix[6] * ds[0] + sys.box_matrix[7] * ds[1] + sys.box_matrix[8] * ds[2];
}
