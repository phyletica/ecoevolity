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
        std::ostream & tree_log_stream,
        std::ostream & state_log_stream,
        std::ostream & operator_log_stream,
        std::ostream & std_output_stream,
        std::ostream & std_error_stream,
        unsigned int logging_precision = 18,
        unsigned int nthreads = 1) {
    if (! tree_log_stream.is_open()) {
        throw EcoevolityError("mcmc: tree log stream is not open");
    }
    if (! state_log_stream.is_open()) {
        throw EcoevolityError("mcmc: state log stream is not open");
    }
    if (! operator_log_stream.is_open()) {
        throw EcoevolityError("mcmc: operator log stream is not open");
    }
    if (! std_output_stream.is_open()) {
        throw EcoevolityError("mcmc: std output is not open");
    }
    if (! std_error_stream.is_open()) {
        throw EcoevolityError("mcmc: std error is not open");
    }
    
    tree_log_stream.precision(logging_precision);
    state_log_stream.precision(logging_precision);
    operator_log_stream.precision(logging_precision);

    write_state_log_header(state_log_stream);
    write_state_log_header(std_output_stream, true);

    this->make_trees_dirty();
    this->compute_log_likelihood_and_prior(true);
    if (this->get_log_likelihood() == -std::numeric_limits<double>::infinity()) {
        std::ostringstream message;
        message << "\n"
                << "\n#######################################################################\n"
                <<   "###############################  ERROR  ###############################\n"
                <<   "The initial model state has a probability density of zero.\n"
                <<   "#######################################################################\n";
        throw EcoevolityError(message.str());
    }
    if (std::isnan(this->get_log_likelihood())) {
        std::ostringstream message;
        message << "\n"
                << "\n#######################################################################\n"
                <<   "###############################  ERROR  ###############################\n"
                <<   "The initial model state has a NAN log likelihood.\n"
                <<   "#######################################################################\n";
        throw EcoevolityError(message.str());
    }
    this->log_state(state_log_stream, 0);
    this->log_state(std::cout, 0, true);

    unsigned int gen;
    unsigned int gen_of_last_state_log = 0;
    unsigned int gen_of_last_operator_log = 0;
    for (gen = 0; gen < chain_length; ++gen) {
        OperatorInterface& op = this->operator_schedule_.draw_operator(rng);
        op.operate(rng, this, this->get_number_of_threads());

        if ((gen + 1) % sample_frequency == 0) {
            this->log_state(state_log_stream, gen + 1);
            gen_of_last_state_log = gen;
            // Log every 10th sample to std out
            if ((gen + 1) % (sample_frequency * 10) == 0) {
                this->log_state(std::cout, gen + 1, true);
                // Log operator performance every 100 samples
                if ((gen + 1) % (sample_frequency * 100) == 0) {
                    operator_log_stream << "generation " << gen + 1 << ":\n";
                    this->operator_schedule_.write_operator_rates(operator_log_stream);
                    gen_of_last_operator_log = gen;
                }
            }
        }

        // Check if the chain has gone astray
        ECOEVOLITY_DEBUG(
            if ((gen + 1) % (sample_frequency * 100) == 0) {
                double chain_ln_likelihood = this->log_likelihood_.get_value();
                this->make_trees_dirty();
                this->compute_log_likelihood_and_prior(true);
                double expected_ln_likelihood = this->log_likelihood_.get_value();
                if (fabs(chain_ln_likelihood - expected_ln_likelihood) > 1e-6) {
                    std::ostringstream message;
                    message << "ERROR: The likelihood at generation "
                            << gen
                            << " is " << chain_ln_likelihood
                            << "; expected " << expected_ln_likelihood
                            << "; last operator " << op.get_name();
                    state_log_stream.close();
                    operator_log_stream.close();
                    throw EcoevolityError(message.str());
                }
            }
        )
    }
    // Make sure last generation is reported
    if (gen > (gen_of_last_state_log + 1)) {
        this->log_state(state_log_stream, gen + 1);
        this->log_state(std::cout, gen + 1, true);
    }
    if (gen > (gen_of_last_operator_log + 1)) {
        operator_log_stream << "generation " << gen + 1 << ":\n";
        this->operator_schedule_.write_operator_rates(operator_log_stream);
    }
    std::cout << "\nOperator stats:\n";
    this->operator_schedule_.write_operator_rates(std::cout);
    std::cout << "\n";
    state_log_stream.close();
    operator_log_stream.close();
}

#endif
