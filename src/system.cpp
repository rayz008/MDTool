#include <fstream>
#include <sstream>
#include <stdexcept>
#include <vector>
#include <iostream>
#include <math.h>
#include <thread>
#include <algorithm>

#include "system.h"
#include "settings.h"


System::~System()
{
    delete[] atoms;
    delete[] coords;
}

void System::AllocateMemory()
{
    atoms = new std::string[natoms];
    coords = new double[nframes*natoms*3];
    box = new double[nframes*6];
    boxmat = new double[9];
    boxinv = new double[9];

    if(!atoms || !coords)
        throw std::bad_alloc();
}


void System::ReadXYZ(const std::string &filename) 
{

    if(atoms || coords)
    {
        throw std::logic_error(
            "Coordinates already in System instance.");
        return;
    }
    std::ifstream file(filename);
    if(file.is_open())
    {
        std::string line;
        //get number of atoms
        std::getline(file, line);
        std::istringstream iss;
        iss.clear();
        iss.str(line);
        iss >> natoms;

        //get number of frames
        nframes++;
        while(getline(file, line)) 
        {
            nframes++;
        }
        nframes /= (natoms+2);
        file.clear();
        file.seekg(0);

        // settings atomtype array and coords array=[x0,y1,z1,x2,y2,z2,....]
        // box array = [a1, b1, c1, A1, B1, C1, a2, b2, c2, ...]
        AllocateMemory();


        for(int i=0; i < nframes; i++) 
        {
            getline(file, line);
            
            // get box information
            getline(file, line);
            if(line.size() != 0)
            {
                line = line.substr(line.find("CELL(abcABC):")+13, line.find("Step:")-line.find("CELL(abcABC):")-13);
                iss.clear();
                iss.str(line);
                double box_tmp;
                int k=0;
                while(iss >> box_tmp){
                    box[6*i+k] = box_tmp;
                    k++;
                }
            }
           
            for(int j=0; j < natoms; j++) 
            {
                // write in atomtypes and coords
                getline(file, line);
                iss.clear();
                iss.str(line);
                iss >> atoms[j]
                    >> coords[i*natoms*3 + j*3]
                    >> coords[i*natoms*3 + j*3 + 1]
                    >> coords[i*natoms*3 + j*3 + 2];
            }
        }

        file.close();

    } 
    else 
    {
        throw std::logic_error("Unexpected file format.");
    }
 
    return;
}

void System::UpdateBoxInv()
{
    for(int i = 0; i < 3; i++) {
         for(int j = 0; j < 3; j++){
             boxinv[i*3+j] = ((boxmat[3*((j + 1)%3) + (i + 1)%3] *
                               boxmat[3*((j + 2)%3) + (i + 2)%3]) -
                              (boxmat[3*((j + 1)%3) + (i + 2)%3] *
                               boxmat[3*((j + 2)%3) + (i + 1)%3])) 
                              /boxvol;
         }
    }

    return;
}

void System::UpdateNPTBox(int nframe)
{
    double PI = 3.141592653589793238L;
    //update box matrix
    boxmat[0] = box[6*nframe];
    boxmat[1] = box[6*nframe+1]*cos(box[6*nframe+5]*PI/180);
    boxmat[2] = box[6*nframe+2]*cos(box[6*nframe+4]*PI/180);
    boxmat[3] = 0;
    boxmat[4] = box[6*nframe+1]*sin(box[6*nframe+5]*PI/180);
    boxmat[5] = ( box[6*nframe+1] * box[6*nframe+2] * cos(box[6*nframe+3]*PI/180)
                - boxmat[1] * boxmat[2])/boxmat[4];
    boxmat[6] = 0;
    boxmat[7] = 0;
    boxmat[8] = sqrt(box[6*nframe+2]*box[6*nframe+2] 
                  - boxmat[2]*boxmat[2] 
                  - boxmat[5]*boxmat[5]);

    //update box volume
    boxvol = 0;
    for(int i{0}; i < 3; i++) 
    {
        boxvol += boxmat[i] *
                    (boxmat[3+(i+1)%3] * boxmat[6+(i+2)%3] -
                     boxmat[3+(i+2)%3] * boxmat[6+(i+1)%3]);
    }

    if(boxvol == 0) 
    {
        throw std::logic_error(
            "Determinant for NPT inputting box is 0, more dimension required!");
    }

    //update box inverse matrix
    UpdateBoxInv();

    return;
}

// add options to read with Cell abcABC
// detects NPT and NVT by the comment line specified in the xyz trajectory
void System::UpdateNVTBox(Settings& settings)
{
    std::string line;
    line = settings.box;

    line.erase(std::remove(line.begin(), line.end(), '['), line.end());
    line.erase(std::remove(line.begin(), line.end(), ']'), line.end());
    
    std::istringstream iss;
    iss.clear();
    iss.str(line);
    iss.ignore(256, ' ');
    iss.ignore(256, ' ');
    
    int i{0};
    double box_tmp;
    // box formated like box = [[xa xb xc] [ya yb yc] [za zb zc]] for a b c as the three edges
    while(iss >> box_tmp) 
    {
        boxmat[i] = box_tmp;
        i++;
    }
    
    //update box volume
    boxvol = 0;
    for(int i{0}; i < 3; i++) 
    {
        boxvol += boxmat[i] *
                    (boxmat[3+(i+1)%3] * boxmat[6+(i+2)%3] -
                     boxmat[3+(i+2)%3] * boxmat[6+(i+1)%3]);
    }
    
    if(boxvol == 0) 
    {
        throw std::logic_error(
            "Determinant for NVT inputting box is 0, more dimension required!");
    }

    //update box inverse matrix
    UpdateBoxInv();

    return;

}


