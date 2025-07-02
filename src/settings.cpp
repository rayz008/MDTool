#include <fstream>
#include <stdexcept>

#include "json.hpp"
#include "settings.h"

using json = nlohmann::json;

/**
 * Settings handles parsing of JSON configuration files 
 *
 */

Settings::Settings(const char* settingfile) {
    readSettings(settingfile);
}

void Settings::readSettings(const char* settingfile) {
    std::ifstream file(settingfile);

    if (!file.is_open()) {
        throw std::runtime_error("Failed to open setting file: " + 
                                 std::string(settingfile) + ", check filename!");
    }

    json settingconfig;
    try {
        file >> settingconfig;
    } catch (const json::parse_error& e) {
        throw std::runtime_error("Failed to parse JSON file: " + std::string(e.what()));
    }
    file.close();

    try {
        // Required configuration fields
        traj_infile = settingconfig.at("trajectory_input").get<std::string>();
        atom1 = settingconfig.at("atom_type_1").get<std::string>();
        atom2 = settingconfig.at("atom_type_2").get<std::string>();

        // Optional configuration fields
        box_infile = settingconfig.value("box_input", std::string(""));
        rdf_outfile = settingconfig.value("rdf_output", std::string("radial.dat"));
        r_min = settingconfig.value("r_min", 0.0);
        r_max = settingconfig.value("r_max", 10.0);
        bins = settingconfig.value("bins", 200);
        increments = settingconfig.value("increment", 0);
        irdf_outfile = settingconfig.value("irdf_output", std::string(""));

        validateSettings();

    } catch (const json::out_of_range& e) {
        throw std::runtime_error("Missing required configuration: " + std::string(e.what()));
    } catch (const json::type_error& e) {
        throw std::runtime_error("Invalid type for configuration field: " + std::string(e.what()));
    }
}

void Settings::validateSettings() const {
    if (r_max <= r_min) {
        throw std::logic_error("r_max must be greater than r_min");
    }
    if (bins <= 0) {
        throw std::logic_error("Number of bins must be positive");
    }
    if (increments < 0) {
        throw std::runtime_error("Number of increments must be non-negative");
    }
    if (traj_infile.empty()) {
        throw std::runtime_error("Trajectory input file need to be specified");
    }
    if (atom1.empty() || atom2.empty()) {
        throw std::runtime_error("Both atom types should be specified");
    }
}

