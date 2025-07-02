#ifndef PBC_H
#define PBC_H

#include "settings.h"

void pbcOrthorhombic(double& dx, double& dy, double& dz, System& sys);
void pbcTriclinic(double& dx, double& dy, double& dz, System& sys);

#endif
