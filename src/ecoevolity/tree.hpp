/******************************************************************************
 * Copyright (C) 2015-2016 Jamie R. Oaks.
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

#ifndef ECOEVOLITY_TREE_HPP
#define ECOEVOLITY_TREE_HPP

#include "data.hpp"
#include "node.hpp"
#include "error.hpp"
#include "assert.hpp"

class PopulationTree {
    private:
        BiallelicData data_;
        PopulationNode * root_;
        double u = 1.0;
        double v = 1.0;
        std::vector<double> pattern_probs_;

        void compute_pattern_probability(unsigned int pattern_index);
        void compute_pattern_probabilities();

    public:
        PopulationTree(
                const std::string path, 
                const char population_name_delimiter = '_',
                const bool population_name_is_prefix = true,
                const bool genotypes_are_diploid = true,
                const bool markers_are_dominant = false,
                const bool validate = true);
        ~PopulationTree () { }
};

#endif
