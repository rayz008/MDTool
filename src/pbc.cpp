#include <math.h>
#include "system.h"

// stores box and boxinverse for a frame
// update box

void PBCOrthorhombic(
    double& dx,
    double& dy,
    double& dz,
    System& sys)
{
    dx -= rint(dx / sys.boxmat[0]) * sys.boxmat[0];
    dy -= rint(dy / sys.boxmat[4]) * sys.boxmat[4];
    dz -= rint(dz / sys.boxmat[8]) * sys.boxmat[8];
}

void PBCTriclinic(
    double& dx,
    double& dy,
    double& dz,
    System& sys)
{
    double ds[3] = {0, 0, 0};
    for(int i = 0; i < 3; i++) 
    {
        ds[i] = rint(sys.boxinv[i * 3]     * dx
                   + sys.boxinv[i * 3 + 1] * dy
                   + sys.boxinv[i * 3 + 2] * dz);
    }
    dx -= sys.boxmat[0] * ds[0] + sys.boxmat[1] * ds[1] + sys.boxmat[2] * ds[2];
    dy -= sys.boxmat[3] * ds[0] + sys.boxmat[4] * ds[1] + sys.boxmat[5] * ds[2];
    dz -= sys.boxmat[6] * ds[0] + sys.boxmat[7] * ds[1] + sys.boxmat[8] * ds[2];
}
