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

#ifndef ECOEVOLITY_GENERAL_TREE_HPP
#define ECOEVOLITY_GENERAL_TREE_HPP

#include "tree.hpp"
#include "error.hpp"
#include "assert.hpp"


class GeneralPopulationTree : public BasePopulationTree {
    protected:
        std::vector< std::shared_ptr<PositiveRealParameter> > node_heights_;

        void init_node_heights();

    public:
        GeneralPopulationTree() { }
        GeneralPopulationTree(
                std::string path, 
                char population_name_delimiter = ' ',
                bool population_name_is_prefix = true,
                bool genotypes_are_diploid = true,
                bool markers_are_dominant = false,
                bool constant_sites_removed = true,
                bool validate = true,
                bool strict_on_constant_sites = false,
                bool strict_on_missing_sites = false,
                bool strict_on_triallelic_sites = true,
                double ploidy = 2.0,
                bool store_seq_loci_info = false
                ) : BasePopulationTree(
                    path,
                    population_name_delimiter,
                    population_name_is_prefix,
                    genotypes_are_diploid,
                    markers_are_dominant,
                    constant_sites_removed,
                    validate,
                    strict_on_constant_sites,
                    strict_on_missing_sites,
                    strict_on_triallelic_sites,
                    ploidy,
                    store_seq_loci_info) {
            this->init_node_heights();
        }

        GeneralPopulationTree(
                std::shared_ptr<PopulationNode> root,
                unsigned int number_of_loci = 10000,
                unsigned int length_of_loci = 1,
                bool validate_data = false
                ) : BasePopulationTree(
                    root,
                    number_of_loci,
                    length_of_loci,
                    validate_data) {
            this->init_node_heights();
        }

        const std::vector< std::shared_ptr<PositiveRealParameter> >& get_node_height_parameters() const {
            return this->node_heights_;
        }

        unsigned int get_number_of_node_heights() const {
            return this->node_heights_.size();
        }

        std::vector<double> get_node_heights() const {
            std::vector<double> heights (this->node_heights_.size());
            for (unsigned int i = 0; i < this->node_heights_.size(); ++i) {
                heights.at(i) = this->node_heights_.at(i)->get_value();
            }
            return heights;
        }
};

#endif
