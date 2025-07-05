#ifndef SETTINGS_H
#define SETTINGS_H

#include <vector>
#include <string>

struct Settings {
    // files
    std::string traj_infile;
    std::string box_infile;
    std::string rdf_outfile;
    std::string irdf_outfile;
    // atoms
    std::string atomA;
    std::string atomB;
    // parameters
    double r_min, r_max;
    int bins, increments;

    void readSettings(const char* filename);
    void validateSettings() const;

    Settings(const char* filename);
    ~Settings() = default;
    Settings(const Settings& other)              = delete;
    Settings& operator = (const Settings& other) = delete;
    Settings(Settings&& other)                   = delete;
    Settings& operator = (Settings&& other)      = delete;
};

#endif
