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

#ifndef ECOEVOLITY_OPTIONS_HPP
#define ECOEVOLITY_OPTIONS_HPP


class EcoevolityOptions {

    public:
        EcoevolityOptions() { }

		enum ModelPrior {
            fixed = 1,
            uniform = 2,
            dpp = 3,
            pyp = 4,
        };

        enum ModelOperator {
            none = 1,
            rj = 2,
            gibbs_dpp = 3,
            gibbs_pyp = 4,
        };

		enum TreePrior {
            uniform_root_and_betas = 1,
        };

		enum TreeSpace {
            generalized     = 1,
            bifurcating     = 2,
        };


        static bool is_tree_prior(const std::string & s) {
            try {
                EcoevolityOptions::get_tree_prior(s);
            }
            catch (...) {
                return false;
            }
            return true;
        }
        static std::string get_tree_prior_string(const EcoevolityOptions::TreePrior tp) {
            if (tp == EcoevolityOptions::TreePrior::uniform_root_and_betas) {
                return "uniform_root_and_betas";
            }
            throw "EcoevolityOptions: Invalid tree prior";
        }
        static EcoevolityOptions::TreePrior get_tree_prior(const std::string & s) {
            if (s == EcoevolityOptions::get_tree_prior_string(
                        EcoevolityOptions::TreePrior::uniform_root_and_betas)) {
                return EcoevolityOptions::TreePrior::uniform_root_and_betas;
            }
            throw "EcoevolityOptions: Invalid tree prior string";
        }

        static bool is_tree_space(const std::string & s) {
            try {
                EcoevolityOptions::get_tree_space(s);
            }
            catch (...) {
                return false;
            }
            return true;
        }
        static std::string get_tree_space_string(const EcoevolityOptions::TreeSpace ts) {
            if (ts == EcoevolityOptions::TreeSpace::generalized) {
                return "generalized";
            }
            if (ts == EcoevolityOptions::TreeSpace::bifurcating) {
                return "bifurcating";
            }
            throw "EcoevolityOptions: Invalid tree space";
        }
        static EcoevolityOptions::TreeSpace get_tree_space(const std::string & s) {
            if (s == EcoevolityOptions::get_tree_space_string(
                        EcoevolityOptions::TreeSpace::generalized)) {
                return EcoevolityOptions::TreeSpace::generalized;
            }
            if (s == EcoevolityOptions::get_tree_space_string(
                        EcoevolityOptions::TreeSpace::bifurcating)) {
                return EcoevolityOptions::TreeSpace::bifurcating;
            }
            throw "EcoevolityOptions: Invalid tree space string";
        }
};

#endif
