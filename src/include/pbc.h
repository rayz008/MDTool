#ifndef PBC_H
#define PBC_H
#include "settings.h"

/**
 * @brief Applies minimum image conversion for orthorhombic pbc box
 *
 * @param[in,out] dx The x-distance between target atom pair
 * @param[in,out] dy The y-distance between target atom pair
 * @param[in,out] dz The z-distance between target atom pair
 * @param[in] sys System containing coordinates and box information
 */
void pbcOrthorhombic(double& dx, double& dy, double& dz, System& sys);

/**
 * @brief Applies minimum image conversion for any triclinic pbc box
 *
 * @param[in,out] dx The x-distance between target atom pair
 * @param[in,out] dy The y-distance between target atom pair
 * @param[in,out] dz The z-distance between target atom pair
 * @param[in] sys System containing coordinates and box information
 */
void pbcTriclinic(double& dx, double& dy, double& dz, System& sys);

#endif
