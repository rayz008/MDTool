#ifndef SETTINGS_H
#define SETTINGS_H
#include <string>

/**
 * @struct Settings
 * @brief Processes parameters used for MD trajectory analysis from JSON input
 */
struct Settings {
    // files
    std::string traj_infile;
    std::string box_infile;
    std::string rdf_outfile;
    std::string irdf_outfile;
    // atoms
    std::string atomA;
    std::string atomB;
    // parameters for rdf and i-rdf
    double r_min, r_max;
    int bins, increments;

    /**
     * @brief Reads setting information from JSON file
     */
    void readSettings(const char* filename);

    /**
     * @brief Validates setting parameters
     */
    void validateSettings() const;

    Settings(const char* filename);
    ~Settings() = default;
    Settings(const Settings& other)              = delete;
    Settings& operator = (const Settings& other) = delete;
    Settings(Settings&& other)                   = delete;
    Settings& operator = (Settings&& other)      = delete;
};

#endif
