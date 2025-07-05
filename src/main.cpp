#include <iostream>
#include <stdexcept>
#include "settings.h"
#include "system.h"
#include "rdf.h"

int main(int argc, char** argv){
    if(argc != 2) {
        std::cerr << "Usage: " << argv[0] << "<settings.json>" << std::endl;
        return 1;
    }

    const char* settings_file(argv[1]);

    try {
        // Load settings and system info
        Settings settings(settings_file);
        System sys;
        sys.readXYZ(settings.traj_infile);
        sys.readBox(settings);

        // Compute (i)RDFs
        RDFCalculator calculator;
        std::cout << "Starting RDF calculation ..." << std::endl;
        calculator.compute(sys, settings);
        std::cout << "RDF calculation completed successfully!" << std::endl;

    } catch(const std::bad_alloc &e) {
        std::cerr << "Memory allocation failed." << std::endl;
        return 1;
    } catch(const std::logic_error &e) {
        std::cerr << e.what() << std::endl;
        return 1;
    }

    return 0;
}
