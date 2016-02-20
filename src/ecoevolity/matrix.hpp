/******************************************************************************
 * Copyright (C) 2016 Jamie R. Oaks.
 *
 * This file is part of Ecoevolity.
 *
 * Ecoevolity is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * Ecoevolity is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with Ecoevolity.  If not, see <http://www.gnu.org/licenses/>.
******************************************************************************/

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
