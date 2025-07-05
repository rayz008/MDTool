#ifndef RDF_H
#define RDF_H

#include <vector>
#include <utility>
#include <string>
#include <memory>

class RDFCalculator {
public:
    RDFCalculator() = default;

    ~RDFCalculator() = default;

    void compute(System& sys, const Settings& settings);

private:
    double dr_;                                   ///< Bin width for RDF calculation
    int num_A_;                                   ///< Number of atoms of type A
    int num_B_;                                   ///< Number of atoms of type B
    double factor_;                               ///< Normalization factor
    std::vector<std::pair<int, int>> pairs_;      ///< Vector of atom pairs to analyze
    
    std::vector<double> g_;                 ///< RDF histogram data
    std::vector<double> incre_g_;           ///< Incremental RDF histogram data
    std::vector<double> minAB_;             ///< Array for storing minimum distances

    void initializeVectors(const Settings& settings);

    void calculateRDF(System& sys, const Settings& settings, int frame);

    void calculateIncrementalRDF(System& sys, const Settings& settings, int frame);

    void normalizeRDF(const Settings& settings, int nframes);

    void normalizeIncrementalRDF(const Settings& settings, int nframes);

    void writeRDFOutput(const Settings& settings) const;

    void writeIncrementalRDFOutput(const Settings& settings) const;

    void refreshMinAB(const Settings& settings, double dAB, int& count);

    std::vector<std::pair<int, int>> generatePairs(
        const System& sys, 
        const std::string& atomA, 
        const std::string& atomB,
        int& numA, 
        int& numB);
};


#endif // RDF_H
