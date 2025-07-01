// Radial Distribution Functions
//#include <cstdio>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdexcept>
//#include <string>
#include <math.h>
//#include <utility>
//#include <vector>
#include <algorithm>

#include "settings.h"
#include "system.h"
#include "tools.h"
#include "pbc.h"


constexpr double PI = 3.141592653589793238L;


void calculate_rdf(System& sys, Settings& settings) {
    // get error message for wrong r_max and r_min;
    if(settings.r_max <= settings.r_min) {
        throw std::logic_error(
            "Maximum radius is smaller than minimum radius!");
    }
    if(settings.bins <= 0) {
        throw std::logic_error(
            "Bin number must be larger than 0!");
    }

    /*
    if(std::none_of(sys.atoms, sys.atoms + sys.natoms,
        [&](std::string str) {return (str == atomA || str == atomB);})) {
        throw std::logic_error(
            "Selected atoms not found.");
    }
    */

    double dx, dy, dz, dAB;
    double dr = (settings.r_max - settings.r_min)/settings.bins;
    double* g = new double[settings.bins];
    double* incre_g = new double[settings.bins * settings.increment];
    double* minAB = new double[settings.increment];

    int numA{0}, numB{0};
    std::vector<std::pair<int, int>> couples{
        calculate_couples(sys, settings.atom1, settings.atom2, numA, numB)};
    double factor = numB*numA*4*PI*dr;

    for(int i{0}; i < settings.bins; i++)
        g[i] = 0;


    for(int i{0}; i < settings.bins * settings.increment; i++)
        incre_g[i] = 0;


    if(settings.NVT){
            sys.UpdateNVTBox(settings);
    }

    for(int frame{0}; frame < sys.nframes; frame++) {
        
        if(!settings.NVT){
            sys.UpdateNPTBox(frame);
        }

        int count = 0;
        int check = couples[0].first;

        for(std::pair<int, int> &couple : couples) {
            dx = sys.coords[frame*sys.natoms*3 + couple.first*3]
                 - sys.coords[frame*sys.natoms*3 + couple.second*3];
            dy = sys.coords[frame*sys.natoms*3 + couple.first*3 + 1]
                 - sys.coords[frame*sys.natoms*3 + couple.second*3 + 1];
            dz = sys.coords[frame*sys.natoms*3 + couple.first*3 + 2]
                 - sys.coords[frame*sys.natoms*3 + couple.second*3 + 2];

            PBCTriclinic(dx, dy, dz, sys);

            dAB = sqrt(dx*dx + dy*dy + dz*dz);

            // check if should include this distance in RDF
            if(dAB < settings.r_max && dAB > settings.r_min) {
                int layer = static_cast<int>((dAB-settings.r_min)/dr);
                g[layer] += sys.boxvol; 
            }

            // check and include distance in minAB array for i-RDF
            auto RefreshMinAB = [&]()
            {
                if(     dAB > settings.r_min
                    &&  dAB < settings.r_max 
                    &&  count < settings.increment) {
                    minAB[count] = dAB;
                    count++;
                }
                else if(     dAB > settings.r_min 
                         &&  dAB < settings.r_max
                         &&  count >= settings.increment){
                    int max_index = 0;
                    for(int k = 1; k < settings.increment; k++){
                        if(minAB[k] > minAB[max_index]){
                            max_index = k;
                        }
                    }
                    if(minAB[max_index] > dAB){
                        minAB[max_index] = dAB;
                    }
                }
            };

            if(settings.increment) {
                if(check == couple.first){
                    RefreshMinAB();
                }
                else{
                    // if changed atom1, load sorted minAB array to incre_g
                    std::sort(minAB, minAB+settings.increment);
                    for(int i{0}; i<settings.increment; i++){
                        int layer = static_cast<int>((minAB[i]
                                                     -settings.r_min)/dr);
                        incre_g[layer+i*settings.bins] += sys.boxvol; 
                    }

                    // zero the temp minAB array
                    for(int i{0}; i < settings.increment; i++){
                        minAB[i] = 0;
                    }
                }
                check = couple.first;
            }
            

        }

        //load the last sorted minAB to incre_g
        std::sort(minAB, minAB+settings.increment);
        for(int i{0}; i<settings.increment; i++){
            int layer = static_cast<int>((minAB[i]-settings.r_min)/dr);
            incre_g[layer+i*settings.bins] += sys.boxvol; 

        }
       
        // zero the temp minAB array
        for(int i{0}; i < settings.increment; i++){
            minAB[i] = 0;
        }

    }

    for(int i{0}; i < settings.bins; i++){
        std::cout << incre_g[i] << std::endl;
    }


    // write rdf output file
    std::ofstream rdffile;
    rdffile.open(settings.outfile);
    rdffile << settings.bins << "  " << settings.atom1 << "  " << settings.atom2 << "\n";
    rdffile << "distance:\tRDF value:\n";
    for(int i{0}; i < settings.bins; i++) {
        //final calculation
        double r = settings.r_min + i*dr;
        if(i != 0) {
            g[i] = g[i]/(factor*sys.nframes*r*r);
        }
        else{
            g[i] = 0;   // set the rdf for 0 radius as 0;
        }
        rdffile << std::fixed << std::setprecision(5) << settings.r_min+i*dr << "\t"
                << std::fixed << std::setprecision(8) << g[i] << "\n";
    }
    rdffile.close();



    // write i-rdf output file
    if(settings.increfile.size() != 0 && settings.increment) {
        std::ofstream irdffile;
        irdffile.open(settings.increfile);
        irdffile << settings.bins << "  " << settings.atom1 << "  " << settings.atom2 << "\n";
        for(int i{0}; i < settings.increment; i++) {
            irdffile << "iRDF: " << i << "\n";
            irdffile << "distance:\tRDF value:\n";
            for(int j{0}; j < settings.bins; j++) {
                double r = settings.r_min + j*dr;
                if(j != 0) {
                    incre_g[j+i*settings.bins] = incre_g[j+i*settings.bins]
                                                      /(factor*sys.nframes*r*r);
                }
                else{
                    incre_g[j+i*settings.bins] = 0;   // set the rdf for 0 radius as 0;
                }
                irdffile << std::fixed << std::setprecision(5) << settings.r_min+j*dr << "\t"
                         << std::fixed << std::setprecision(8) << incre_g[j+i*settings.bins] << "\n";
            }
            irdffile << "\n";
        }
    }

    delete [] g;
    delete [] incre_g;
    delete [] minAB;
}


int main(int argc, char** argv){
    if(argc != 2)
    {
        printf("Example: %s settings.ini\n", argv[0]);
        return 1;
    }

    const char* settings_file(argv[1]);

    //read all basic settings from settings.ini file
    struct Settings settings(settings_file);
    struct System sys;

    // read and check stored system info
    try {
        sys.ReadXYZ(settings.infile);
    } catch(const std::bad_alloc &e) {
        std::cerr << "Memory allocation failed." << std::endl;
        return 1;
    } catch(const std::logic_error &e) {
        std::cerr << e.what() << std::endl;
        return 1;
    }


    // calculate g info
    try {
        calculate_rdf(sys, settings);
    } catch(const std::logic_error &e) {
        std::cerr << e.what() << std::endl;
    }

    return 0;

}
