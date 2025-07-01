#ifndef SETTINGS_H
#define SETTINGS_H

#include <vector>
#include <string>

struct Settings {
    //files
    std::string infile;
    std::string outfile;
    std::string increfile;   // new
    //system
    std::string atom1;
    std::string atom2;
    std::string box;
    bool NVT{false};
    //settings
    double r_min{0}, r_max{0};
    int bins{0};
    int increment{0};        // new

    void ReadSettings(const char* filename);

    Settings(const char* filename);
    ~Settings() = default;
    Settings(const Settings& other)              = delete;
    Settings& operator = (const Settings& other) = delete;
    Settings(Settings&& other)                   = delete;
    Settings& operator = (Settings&& other)      = delete;
};

#endif
