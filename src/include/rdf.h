#ifndef RDF_H
#define RDF_H
#include <vector>
#include <utility>
#include <string>
#include <memory>

/**
 * @class RDFCalculator
 * @brief Calculates radial distribution functions from MD trajectories
 *
 * The RDFCalculator class processes atomic coordinates from trajectory files
 * to compute standard RDFs and incremental RDFs (iRDFs). It handles:
 * - Pair generation between specified atom types
 * - Distance calculations with periodic boundary conditions
 * - Histogram binning and normalization
 * - Output file generation
 *
 * @note System struct assumes that trajectory frames contain atoms in consistent order
 */
class RDFCalculator {
public:
    RDFCalculator() = default;
    ~RDFCalculator() = default;

    /**
     * @brief Compute RDF and iRDF based on setting information
     */
    void compute(System& sys, const Settings& settings);

private:
    double dr_;                                   // Bin width for RDF calculation
    int num_A_;                                   // Number of atoms of type A
    int num_B_;                                   // Number of atoms of type B
    double factor_;                               // Normalization factor
    std::vector<std::pair<int, int>> pairs_;      // Vector of atom pairs to analyze
    
    std::vector<double> g_;                 // RDF histogram data
    std::vector<double> incre_g_;           // Incremental RDF histogram data
    std::vector<double> minAB_;             // Array for storing minimum distances

    /**
     * @brief Initialize necessary vectors g_, incre_g_ and minAB_ based on settings
     */
    void initializeVectors(const Settings& settings);

    /**
     * @brief Radial distribution function calculation
     *
     * @param[in] sys System information
     * @param[in] settings Setting information
     * @param[in] frame Current frame index
     */
    void calculateRDF(System& sys, const Settings& settings, int frame);

    /**
     * @brief Incremental radial distribution function calculation
     *
     * @param[in] sys System information
     * @param[in] settings Setting information
     * @param[in] frame Current frame index
     */
    void calculateIncrementalRDF(System& sys, const Settings& settings, int frame);

    /**
     * @brief Normalize RDF g_
     */
    void normalizeRDF(const Settings& settings, int nframes);

    /**
     * @brief Normalize iRDF incre_g_
     */
    void normalizeIncrementalRDF(const Settings& settings, int nframes);

    /**
     * @brief Write g_ to output file
     */
    void writeRDFOutput(const Settings& settings) const;

    /**
     * @brief Write incre_g_ to output file
     */
    void writeIncrementalRDFOutput(const Settings& settings) const;

    /**
     * @brief Renews minAB array for iRDF calculation
     */
    void refreshMinAB(const Settings& settings, double dAB, int& count);

    /**
     * @brief Generates all index pairs for atomA and atomB from given system atom list
     *
     * @param[in] sys System information
     * @param[in] atomA Name of atomA
     * @param[in] atomB Name of atomB
     * @param[in,out] num_A The number of atom A contained in system
     * @param[in,out] num_B The number of atom B contained in system
     */
    std::vector<std::pair<int, int>> generatePairs(
        const System& sys, 
        const std::string& atomA, 
        const std::string& atomB,
        int& numA, 
        int& numB);
};


#endif // RDF_H
