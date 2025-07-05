// Radial Distribution Functions
#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <math.h>
#include <algorithm>

#include "settings.h"
#include "system.h"
#include "tools.h"
#include "pbc.h"
#include "rdf.h"

constexpr double PI = 3.141592653589793238L;

void RDFCalculator::initializeVectors(const Settings& settings) {
    // Calculate bin width
    dr_ = (settings.r_max - settings.r_min) / settings.bins;

    // Initialize RDF vector
    g_.clear();
    g_.resize(settings.bins, 0.0);

    // Initialize incremental RDF vectors if needed
    if (settings.increments > 0) {
        incre_g_.clear();
        incre_g_.resize(settings.bins * settings.increments, 0.0);

        minAB_.clear();
        minAB_.resize(settings.increments, 0.0);
    } else {
        // Clear these vectors if increments is 0
        incre_g_.clear();
        minAB_.clear();
    }
}

void RDFCalculator::compute(System& sys, const Settings& settings) {
    initializeVectors(settings);
    num_A_ = 0;
    num_B_ = 0;
    pairs_ = generatePairs(sys, settings.atomA, settings.atomB, num_A_, num_B_);
    factor_ = num_A_ * num_B_ * 4 * PI * dr_;

    for (int frame = 0; frame < sys.nframes; frame++) {
        sys.updateBoxInformation(frame);
        calculateRDF(sys, settings, frame);

        if (settings.increments > 0) {
            calculateIncrementalRDF(sys, settings, frame);
        }
    }

    normalizeRDF(settings, sys.nframes);
    if (settings.increments > 0) {
        normalizeIncrementalRDF(settings, sys.nframes);
    }
    //smoothRDF(settings);

    writeRDFOutput(settings);
    if (!settings.irdf_outfile.empty() && settings.increments > 0) {
        writeIncrementalRDFOutput(settings);
    }
}

void RDFCalculator::calculateRDF(System& sys, const Settings& settings, int frame) {
    for(std::pair<int, int> &pair : pairs_) {
        double dx = sys.coords[frame*sys.natoms * 3 + pair.first * 3]
                  - sys.coords[frame*sys.natoms * 3 + pair.second * 3];
        double dy = sys.coords[frame*sys.natoms * 3 + pair.first * 3 + 1]
                  - sys.coords[frame*sys.natoms * 3 + pair.second * 3 + 1];
        double dz = sys.coords[frame*sys.natoms * 3 + pair.first * 3 + 2]
                  - sys.coords[frame*sys.natoms * 3 + pair.second * 3 + 2];

        pbcTriclinic(dx, dy, dz, sys);

        double dAB = sqrt(dx * dx + dy * dy + dz * dz);

        if(dAB < settings.r_max && dAB >= settings.r_min) {
            int layer = static_cast<int>((dAB - settings.r_min) / dr_);
            if (layer >= 0 && layer < settings.bins) {
                std::cout << sys.box_volume << std::endl;
                g_[layer] += sys.box_volume; 
            }
        }
    }

}

void RDFCalculator::calculateIncrementalRDF(System& sys, const Settings& settings, int frame) {
    int count = 0;
    int check = pairs_.empty() ? -1 : pairs_[0].first;

    for(std::pair<int, int> &pair : pairs_) {
        double dx = sys.coords[frame*sys.natoms * 3 + pair.first * 3]
                  - sys.coords[frame*sys.natoms * 3 + pair.second * 3];
        double dy = sys.coords[frame*sys.natoms * 3 + pair.first * 3 + 1]
                  - sys.coords[frame*sys.natoms * 3 + pair.second * 3 + 1];
        double dz = sys.coords[frame*sys.natoms * 3 + pair.first * 3 + 2]
                  - sys.coords[frame*sys.natoms * 3 + pair.second * 3 + 2];

        pbcTriclinic(dx, dy, dz, sys); // TODO: add usage for pbcOrth

        double dAB = sqrt(dx * dx + dy * dy + dz * dz);

        if (check == pair.first) {
            refreshMinAB(settings, dAB, count);
        } else {
            // if changed atom1, load sorted minAB array to incre_g
            std::sort(minAB_.begin(), minAB_.end());
            for(int i = 0; i < settings.increments; i++) {
                int layer = static_cast<int>((minAB_[i] -settings.r_min) / dr_);
                if (layer >= 0 && layer < settings.bins) {
                    incre_g_[layer + i*settings.bins] += sys.box_volume; 
                }
            }

            for(int i = 0; i < settings.increments; i++) {
                minAB_[i] = 0;
            }
        }
        check = pair.first;
    }

    std::sort(minAB_.begin(), minAB_.end());
    for (int i = 0; i < settings.increments; i++) {
        int layer = static_cast<int>((minAB_[i] - settings.r_min) / dr_);
        if (layer >= 0 && layer < settings.bins) {
            incre_g_[layer + i*settings.bins] += sys.box_volume; 
        }
    }
    
    for(int i = 0; i < settings.increments; i++) {
        minAB_[i] = 0.0;
    }
}


void RDFCalculator::normalizeRDF(const Settings& settings, int nframes) {
    for (int i = 0; i < settings.bins; i++) {
        double r = settings.r_min + i * dr_;
        if (i != 0) {
            g_[i] = g_[i] / (factor_ * nframes * r * r);
        } else {
            g_[i] = 0;
        }
    }
}

void RDFCalculator::normalizeIncrementalRDF(const Settings& settings, int nframes) {
    for (int i = 0; i < settings.increments; i++) {
        for (int j = 0; j < settings.bins; j++) {
            double r = settings.r_min + j * dr_;
            if (j != 0) {
                incre_g_[j + i * settings.bins] = incre_g_[j + i * settings.bins]
                                                  / (factor_ * nframes * r * r);
            } else {
                incre_g_[j + i * settings.bins] = 0;
            }
        }
    }
}


void RDFCalculator::writeRDFOutput(const Settings& settings) const {
    std::ofstream rdffile(settings.rdf_outfile);
    if (!rdffile.is_open()) {
        throw std::runtime_error("Failed to write RDF output file");
    }

    rdffile << settings.bins << "  " << settings.atomA << "  " << settings.atomB << "\n";
    rdffile << "distance:\tRDF value:\n";

    for (int i = 0; i < settings.bins; i++) {
        double r = settings.r_min + i * dr_;
        rdffile << std::fixed << std::setprecision(5) << r << "\t"
                << std::fixed << std::setprecision(8) << g_[i] << "\n";
    }

    rdffile.close();
}

void RDFCalculator::writeIncrementalRDFOutput(const Settings& settings) const {
    std::ofstream irdffile(settings.irdf_outfile);
    if (!irdffile.is_open()) {
        throw std::runtime_error("Failed to write incremental RDF output file");
    }

    irdffile << settings.bins << "  " << settings.atomA << "  " << settings.atomB << "\n";

    for (int i = 0; i < settings.increments; i++) {
        irdffile << "iRDF: " << i << "\n";
        irdffile << "distance:\tRDF value:\n";

        for (int j = 0; j < settings.bins; j++) {
            double r = settings.r_min + j * dr_;
            irdffile << std::fixed << std::setprecision(5) << r << "\t"
                     << std::fixed << std::setprecision(8) << incre_g_[j + i * settings.bins] << "\n";
        }
    }

    irdffile.close();
}


void RDFCalculator::refreshMinAB(const Settings& settings, double dAB, int& count) {
    if (dAB > settings.r_min && dAB < settings.r_max && count < settings.increments) {
        minAB_[count] = dAB;
        count++;
    } else if (dAB > settings.r_min && dAB < settings.r_max && count >= settings.increments) {
        int max_index = 0;
        for (int k = 1; k < settings.increments; k++) {
            if (minAB_[k] > minAB_[max_index]) {
                max_index = k;
            }
        }
        if (minAB_[max_index] > dAB) {
            minAB_[max_index] = dAB;
        }
    }
}

std::vector<std::pair<int, int>> RDFCalculator::generatePairs(
    const System& sys, const std::string &atomA, const std::string &atomB,
    int &numA, int &numB) {
    std::vector<std::pair<int, int>> pair_vector;

    for (int a = 0; a < sys.natoms; a++) {
        if (sys.atoms[a] == atomA) {
            numA++;
        } else if(sys.atoms[a] == atomB) {
            numB++;
        }
    }

    if (atomA == atomB) {
        numB = numA;
    }

    pair_vector.reserve(numA * numB);
    for (int a = 0; a < sys.natoms; a++) {
        if (sys.atoms[a] == atomA) {
            for (int b = 0; b < sys.natoms; b++) {
                if (sys.atoms[b] == atomB) {
                    pair_vector.push_back(std::make_pair(a, b));
                    //std::cout << "pairs: " << a << " , " << b << std::endl;
                }
            }
        }
    }
    return pair_vector;
}

