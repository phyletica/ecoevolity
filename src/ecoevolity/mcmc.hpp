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

#ifndef ECOEVOLITY_MCMC_HPP
#define ECOEVOLITY_MCMC_HPP

#include "assert.hpp"
#include "error.hpp"


template<class TreeType>
inline void mcmc(
        RandomNumberGenerator & rng,
        TreeType & tree,
        GeneralTreeOperatorSchedule<TreeType> & operator_schedule,
        unsigned int chain_length,
        unsigned int sample_frequency,
        unsigned int number_of_moves_per_generation,
        const std::string & tree_log_path,
        const std::string & state_log_path,
        const std::string & operator_log_path,
        unsigned int nthreads) {
}

#endif
