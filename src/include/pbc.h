#ifndef PBC_H
#define PBC_H

#include "settings.h"

void PBCOrthorhombic(
    double& dx,
    double& dy,
    double& dz,
    System& sys);
void PBCTriclinic(
    double& dx,
    double& dy,
    double& dz,
    System& sys);

#endif
