#ifndef ECOEVOLITY_MATRIX_HPP
#define ECOEVOLITY_MATRIX_HPP

#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "debug.hpp"
#include "assert.hpp"
#include "error.hpp"

class BiallelicPatternProbabilityMatrix {
    public:
        // Constructor
        BiallelicPatternProbabilityMatrix() { }
        BiallelicPatternProbabilityMatrix(unsigned int allele_count);
        BiallelicPatternProbabilityMatrix(
                const BiallelicPatternProbabilityMatrix& matrix);
        // Construct a leaf likelihood
        BiallelicPatternProbabilityMatrix(
                unsigned int allele_count,
                unsigned int red_allele_count);
        // Construct a top-of-branch likelihood
        BiallelicPatternProbabilityMatrix(
                unsigned int allele_count,
                const std::vector<double>& pattern_probs);
        // Destructor
        // ~AlleleProbabilityMatrix();
        
        //Methods
        const double& get_pattern_probability(
                unsigned int allele_count,
                unsigned int red_allele_count) const;
        void set_pattern_probability(
                unsigned int allele_count,
                unsigned int red_allele_count,
                double probability);
        const unsigned int& get_allele_count() const;
        const std::vector<double>& get_pattern_prob_matrix() const;

        void resize(unsigned int allele_count);
        void reset(unsigned int allele_count);
        void copy(const BiallelicPatternProbabilityMatrix& m);

        std::string to_string() const;

    private:
        unsigned int allele_count_ = 0;
        std::vector<double> pattern_prob_matrix_;
};

#endif
