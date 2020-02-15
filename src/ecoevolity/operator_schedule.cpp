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

#include "operator_schedule.hpp"
#include "operator.hpp"


OperatorSchedule::OperatorSchedule(const CollectionSettings& collection_settings) {
    const OperatorScheduleSettings & settings = collection_settings.get_operator_schedule_settings();

    this->turn_off_auto_optimize();
    if (settings.auto_optimizing()) {
        this->turn_on_auto_optimize();
    }
    this->set_auto_optimize_delay(settings.get_auto_optimize_delay());

    if (settings.get_model_operator_settings().get_weight() > 0.0) {
        if (collection_settings.get_model_operator() == EcoevolityOptions::ModelOperator::gibbs_dpp) {
            this->add_operator(std::make_shared<DirichletProcessGibbsSampler>(
                    settings.get_model_operator_settings().get_weight(),
                    settings.get_model_operator_settings().get_number_of_auxiliary_categories()));
        }
        else if (collection_settings.get_model_operator() == EcoevolityOptions::ModelOperator::gibbs_pyp) {
            this->add_operator(std::make_shared<PitmanYorProcessGibbsSampler>(
                    settings.get_model_operator_settings().get_weight(),
                    settings.get_model_operator_settings().get_number_of_auxiliary_categories()));
        }
        else if (collection_settings.get_model_operator() == EcoevolityOptions::ModelOperator::rj) {
            this->add_operator(std::make_shared<ReversibleJumpSampler>(
                    settings.get_model_operator_settings().get_weight()));
        }
        else if (collection_settings.get_model_operator() == EcoevolityOptions::ModelOperator::none) {
            // Do not add a model operator
        }
        else {
            std::ostringstream message;
            message << "ERROR: Unexpected EcoevolityOptions::ModelOperator \'"
                    << (int)collection_settings.get_model_operator()
                    << "\'\n";
            throw EcoevolityError(message.str());
        }
    }

    if (settings.get_concentration_scaler_settings().get_weight() > 0.0) {
        this->add_operator(std::make_shared<ConcentrationScaler>(
                settings.get_concentration_scaler_settings().get_weight(),
                settings.get_concentration_scaler_settings().get_scale()
                ));
    }

    if (settings.get_discount_scaler_settings().get_weight() > 0.0) {
        this->add_operator(std::make_shared<DiscountMixer>(
                settings.get_discount_scaler_settings().get_weight(),
                settings.get_discount_scaler_settings().get_scale()
                ));
    }

    if (settings.get_discount_mover_settings().get_weight() > 0.0) {
        this->add_operator(std::make_shared<DiscountMover>(
                settings.get_discount_mover_settings().get_weight(),
                settings.get_discount_mover_settings().get_window()
                ));
    }

    if (settings.get_time_size_rate_mixer_settings().get_weight() > 0.0) {
        this->add_operator(std::make_shared<TimeSizeRateMixer>(
                settings.get_time_size_rate_mixer_settings().get_weight(),
                settings.get_time_size_rate_mixer_settings().get_scale()
                ));
    }

    if (settings.get_time_root_size_mixer_settings().get_weight() > 0.0) {
        this->add_operator(std::make_shared<TimeRootSizeMixer>(
                settings.get_time_root_size_mixer_settings().get_weight(),
                settings.get_time_root_size_mixer_settings().get_scale()
                ));
    }

    if (settings.get_time_size_rate_scaler_settings().get_weight() > 0.0) {
        this->add_operator(std::make_shared<TimeSizeRateScaler>(
                settings.get_time_size_rate_scaler_settings().get_weight(),
                settings.get_time_size_rate_scaler_settings().get_scale()
                ));
    }

    if (settings.get_event_time_scaler_settings().get_weight() > 0.0) {
        this->add_operator(std::make_shared<EventTimeScaler>(
                settings.get_event_time_scaler_settings().get_weight(),
                settings.get_event_time_scaler_settings().get_scale()
                ));
    }

    for (unsigned int i = 0; i < collection_settings.get_number_of_comparisons(); ++i) {
        auto comp_settings = collection_settings.get_comparison_setting(i);

        if (comp_settings.get_operator_settings().get_time_size_rate_mixer_settings().get_weight() > 0.0) {
            this->add_operator(std::make_shared<TimeSizeRateMixer>(
                    i,
                    comp_settings.get_operator_settings().get_time_size_rate_mixer_settings().get_weight(),
                    comp_settings.get_operator_settings().get_time_size_rate_mixer_settings().get_scale()
                    ));
        }

        if (comp_settings.get_operator_settings().get_time_size_rate_scaler_settings().get_weight() > 0.0) {
            this->add_operator(std::make_shared<TimeSizeRateScaler>(
                    i,
                    comp_settings.get_operator_settings().get_time_size_rate_scaler_settings().get_weight(),
                    comp_settings.get_operator_settings().get_time_size_rate_scaler_settings().get_scale()
                    ));
        }

        if (comp_settings.get_operator_settings().get_event_time_scaler_settings().get_weight() > 0.0) {
            this->add_operator(std::make_shared<EventTimeScaler>(
                    i,
                    comp_settings.get_operator_settings().get_event_time_scaler_settings().get_weight(),
                    comp_settings.get_operator_settings().get_event_time_scaler_settings().get_scale()
                    ));
        }

        if (comp_settings.get_operator_settings().get_mutation_rate_scaler_settings().get_weight() > 0.0) {
            this->add_operator(std::make_shared<MutationRateScaler>(
                    i,
                    comp_settings.get_operator_settings().get_mutation_rate_scaler_settings().get_weight(),
                    comp_settings.get_operator_settings().get_mutation_rate_scaler_settings().get_scale()
                    ));
        }

        if (comp_settings.get_operator_settings().get_root_population_size_scaler_settings().get_weight() > 0.0) {
            this->add_operator(std::make_shared<RootPopulationSizeScaler>(
                    i,
                    comp_settings.get_operator_settings().get_root_population_size_scaler_settings().get_weight(),
                    comp_settings.get_operator_settings().get_root_population_size_scaler_settings().get_scale()
                    ));
        }

        if (comp_settings.get_operator_settings().get_leaf_population_size_scaler_settings().get_weight() > 0.0) {
            this->add_operator(std::make_shared<LeafPopulationSizeScaler>(
                    i,
                    comp_settings.get_operator_settings().get_leaf_population_size_scaler_settings().get_weight(),
                    comp_settings.get_operator_settings().get_leaf_population_size_scaler_settings().get_scale()
                    ));
        }

        if (comp_settings.get_operator_settings().get_freq_mover_settings().get_weight() > 0.0) {
            this->add_operator(std::make_shared<FreqMover>(
                    i,
                    comp_settings.get_operator_settings().get_freq_mover_settings().get_weight(),
                    comp_settings.get_operator_settings().get_freq_mover_settings().get_window()
                    ));
        }

        if (comp_settings.get_operator_settings().get_time_root_size_mixer_settings().get_weight() > 0.0) {
            this->add_operator(std::make_shared<TimeRootSizeMixer>(
                    i,
                    comp_settings.get_operator_settings().get_time_root_size_mixer_settings().get_weight(),
                    comp_settings.get_operator_settings().get_time_root_size_mixer_settings().get_scale()
                    ));
        }
    }
}

OperatorSchedule::OperatorSchedule(
        const RelativeRootCollectionSettings& collection_settings) {
    const OperatorScheduleSettings & settings = collection_settings.get_operator_schedule_settings();

    this->turn_off_auto_optimize();
    if (settings.auto_optimizing()) {
        this->turn_on_auto_optimize();
    }
    this->set_auto_optimize_delay(settings.get_auto_optimize_delay());

    if (settings.get_model_operator_settings().get_weight() > 0.0) {
        if (collection_settings.get_model_operator() == EcoevolityOptions::ModelOperator::gibbs_dpp) {
            this->add_operator(std::make_shared<DirichletProcessGibbsSampler>(
                    settings.get_model_operator_settings().get_weight(),
                    settings.get_model_operator_settings().get_number_of_auxiliary_categories()));
        }
        else if (collection_settings.get_model_operator() == EcoevolityOptions::ModelOperator::gibbs_pyp) {
            this->add_operator(std::make_shared<PitmanYorProcessGibbsSampler>(
                    settings.get_model_operator_settings().get_weight(),
                    settings.get_model_operator_settings().get_number_of_auxiliary_categories()));
        }
        else if (collection_settings.get_model_operator() == EcoevolityOptions::ModelOperator::rj) {
            this->add_operator(std::make_shared<ReversibleJumpSampler>(
                    settings.get_model_operator_settings().get_weight()));
        }
        else if (collection_settings.get_model_operator() == EcoevolityOptions::ModelOperator::none) {
            // Do not add a model operator
        }
        else {
            std::ostringstream message;
            message << "ERROR: Unexpected EcoevolityOptions::ModelOperator \'"
                    << (int)collection_settings.get_model_operator()
                    << "\'\n";
            throw EcoevolityError(message.str());
        }
    }

    if (settings.get_concentration_scaler_settings().get_weight() > 0.0) {
        this->add_operator(std::make_shared<ConcentrationScaler>(
                settings.get_concentration_scaler_settings().get_weight(),
                settings.get_concentration_scaler_settings().get_scale()
                ));
    }

    if (settings.get_discount_scaler_settings().get_weight() > 0.0) {
        this->add_operator(std::make_shared<DiscountMixer>(
                settings.get_discount_scaler_settings().get_weight(),
                settings.get_discount_scaler_settings().get_scale()
                ));
    }

    if (settings.get_discount_mover_settings().get_weight() > 0.0) {
        this->add_operator(std::make_shared<DiscountMover>(
                settings.get_discount_mover_settings().get_weight(),
                settings.get_discount_mover_settings().get_window()
                ));
    }

    if (settings.get_time_size_rate_mixer_settings().get_weight() > 0.0) {
        this->add_operator(std::make_shared<TimeSizeRateMixer>(
                settings.get_time_size_rate_mixer_settings().get_weight(),
                settings.get_time_size_rate_mixer_settings().get_scale()
                ));
    }

    if (settings.get_time_root_size_mixer_settings().get_weight() > 0.0) {
        this->add_operator(std::make_shared<TimeRootSizeMixer>(
                settings.get_time_root_size_mixer_settings().get_weight(),
                settings.get_time_root_size_mixer_settings().get_scale()
                ));
    }

    if (settings.get_time_size_rate_scaler_settings().get_weight() > 0.0) {
        this->add_operator(std::make_shared<TimeSizeRateScaler>(
                settings.get_time_size_rate_scaler_settings().get_weight(),
                settings.get_time_size_rate_scaler_settings().get_scale()
                ));
    }

    if (settings.get_event_time_scaler_settings().get_weight() > 0.0) {
        this->add_operator(std::make_shared<EventTimeScaler>(
                settings.get_event_time_scaler_settings().get_weight(),
                settings.get_event_time_scaler_settings().get_scale()
                ));
    }

    for (unsigned int i = 0; i < collection_settings.get_number_of_comparisons(); ++i) {
        auto comp_settings = collection_settings.get_comparison_setting(i);
        if (comp_settings.get_operator_settings().get_time_size_rate_mixer_settings().get_weight() > 0.0) {
            this->add_operator(std::make_shared<TimeSizeRateMixer>(
                    i,
                    comp_settings.get_operator_settings().get_time_size_rate_mixer_settings().get_weight(),
                    comp_settings.get_operator_settings().get_time_size_rate_mixer_settings().get_scale()
                    ));
        }

        if (comp_settings.get_operator_settings().get_time_size_rate_scaler_settings().get_weight() > 0.0) {
            this->add_operator(std::make_shared<TimeSizeRateScaler>(
                    i,
                    comp_settings.get_operator_settings().get_time_size_rate_scaler_settings().get_weight(),
                    comp_settings.get_operator_settings().get_time_size_rate_scaler_settings().get_scale()
                    ));
        }

        if (comp_settings.get_operator_settings().get_event_time_scaler_settings().get_weight() > 0.0) {
            this->add_operator(std::make_shared<EventTimeScaler>(
                    i,
                    comp_settings.get_operator_settings().get_event_time_scaler_settings().get_weight(),
                    comp_settings.get_operator_settings().get_event_time_scaler_settings().get_scale()
                    ));
        }

        if (comp_settings.get_operator_settings().get_mutation_rate_scaler_settings().get_weight() > 0.0) {
            this->add_operator(std::make_shared<MutationRateScaler>(
                    i,
                    comp_settings.get_operator_settings().get_mutation_rate_scaler_settings().get_weight(),
                    comp_settings.get_operator_settings().get_mutation_rate_scaler_settings().get_scale()
                    ));
        }

        if (comp_settings.get_operator_settings().get_root_population_size_scaler_settings().get_weight() > 0.0) {
            this->add_operator(std::make_shared<RootPopulationSizeScaler>(
                    i,
                    comp_settings.get_operator_settings().get_root_population_size_scaler_settings().get_weight(),
                    comp_settings.get_operator_settings().get_root_population_size_scaler_settings().get_scale()
                    ));
        }

        if (comp_settings.get_operator_settings().get_leaf_population_size_scaler_settings().get_weight() > 0.0) {
            this->add_operator(std::make_shared<LeafPopulationSizeScaler>(
                    i,
                    comp_settings.get_operator_settings().get_leaf_population_size_scaler_settings().get_weight(),
                    comp_settings.get_operator_settings().get_leaf_population_size_scaler_settings().get_scale()
                    ));
        }

        if (comp_settings.get_operator_settings().get_freq_mover_settings().get_weight() > 0.0) {
            this->add_operator(std::make_shared<FreqMover>(
                    i,
                    comp_settings.get_operator_settings().get_freq_mover_settings().get_weight(),
                    comp_settings.get_operator_settings().get_freq_mover_settings().get_window()
                    ));
        }

        if (comp_settings.get_operator_settings().get_time_root_size_mixer_settings().get_weight() > 0.0) {
            this->add_operator(std::make_shared<TimeRootSizeMixer>(
                    i,
                    comp_settings.get_operator_settings().get_time_root_size_mixer_settings().get_weight(),
                    comp_settings.get_operator_settings().get_time_root_size_mixer_settings().get_scale()
                    ));
        }
    }
}

OperatorSchedule::OperatorSchedule(
        const DirichletCollectionSettings& collection_settings) {
    const OperatorScheduleSettings & settings = collection_settings.get_operator_schedule_settings();

    this->turn_off_auto_optimize();
    if (settings.auto_optimizing()) {
        this->turn_on_auto_optimize();
    }
    this->set_auto_optimize_delay(settings.get_auto_optimize_delay());

    if (settings.get_model_operator_settings().get_weight() > 0.0) {
        if (collection_settings.get_model_operator() == EcoevolityOptions::ModelOperator::gibbs_dpp) {
            this->add_operator(std::make_shared<DirichletProcessGibbsSampler>(
                    settings.get_model_operator_settings().get_weight(),
                    settings.get_model_operator_settings().get_number_of_auxiliary_categories()));
        }
        else if (collection_settings.get_model_operator() == EcoevolityOptions::ModelOperator::gibbs_pyp) {
            this->add_operator(std::make_shared<PitmanYorProcessGibbsSampler>(
                    settings.get_model_operator_settings().get_weight(),
                    settings.get_model_operator_settings().get_number_of_auxiliary_categories()));
        }
        else if (collection_settings.get_model_operator() == EcoevolityOptions::ModelOperator::rj) {
            this->add_operator(std::make_shared<ReversibleJumpSampler>(
                    settings.get_model_operator_settings().get_weight()));
        }
        else if (collection_settings.get_model_operator() == EcoevolityOptions::ModelOperator::none) {
            // Do not add a model operator
        }
        else {
            std::ostringstream message;
            message << "ERROR: Unexpected EcoevolityOptions::ModelOperator \'"
                    << (int)collection_settings.get_model_operator()
                    << "\'\n";
            throw EcoevolityError(message.str());
        }
    }

    if (settings.get_concentration_scaler_settings().get_weight() > 0.0) {
        this->add_operator(std::make_shared<ConcentrationScaler>(
                settings.get_concentration_scaler_settings().get_weight(),
                settings.get_concentration_scaler_settings().get_scale()
                ));
    }

    if (settings.get_discount_scaler_settings().get_weight() > 0.0) {
        this->add_operator(std::make_shared<DiscountMixer>(
                settings.get_discount_scaler_settings().get_weight(),
                settings.get_discount_scaler_settings().get_scale()
                ));
    }

    if (settings.get_discount_mover_settings().get_weight() > 0.0) {
        this->add_operator(std::make_shared<DiscountMover>(
                settings.get_discount_mover_settings().get_weight(),
                settings.get_discount_mover_settings().get_window()
                ));
    }

    if (settings.get_time_size_rate_mixer_settings().get_weight() > 0.0) {
        this->add_operator(std::make_shared<TimeMeanSizeRateMixer>(
                settings.get_time_size_rate_mixer_settings().get_weight(),
                settings.get_time_size_rate_mixer_settings().get_scale()
                ));
    }

    if (settings.get_time_size_rate_scaler_settings().get_weight() > 0.0) {
        this->add_operator(std::make_shared<TimeMeanSizeRateScaler>(
                settings.get_time_size_rate_scaler_settings().get_weight(),
                settings.get_time_size_rate_scaler_settings().get_scale()
                ));
    }

    if (settings.get_event_time_scaler_settings().get_weight() > 0.0) {
        this->add_operator(std::make_shared<EventTimeScaler>(
                settings.get_event_time_scaler_settings().get_weight(),
                settings.get_event_time_scaler_settings().get_scale()
                ));
    }

    for (unsigned int i = 0; i < collection_settings.get_number_of_comparisons(); ++i) {
        auto comp_settings = collection_settings.get_comparison_setting(i);

        if (comp_settings.get_operator_settings().get_time_size_rate_mixer_settings().get_weight() > 0.0) {
            this->add_operator(std::make_shared<TimeMeanSizeRateMixer>(
                    i,
                    comp_settings.get_operator_settings().get_time_size_rate_mixer_settings().get_weight(),
                    comp_settings.get_operator_settings().get_time_size_rate_mixer_settings().get_scale()
                    ));
        }

        if (comp_settings.get_operator_settings().get_time_size_rate_scaler_settings().get_weight() > 0.0) {
            this->add_operator(std::make_shared<TimeMeanSizeRateScaler>(
                    i,
                    comp_settings.get_operator_settings().get_time_size_rate_scaler_settings().get_weight(),
                    comp_settings.get_operator_settings().get_time_size_rate_scaler_settings().get_scale()
                    ));
        }

        if (comp_settings.get_operator_settings().get_event_time_scaler_settings().get_weight() > 0.0) {
            this->add_operator(std::make_shared<EventTimeScaler>(
                    i,
                    comp_settings.get_operator_settings().get_event_time_scaler_settings().get_weight(),
                    comp_settings.get_operator_settings().get_event_time_scaler_settings().get_scale()
                    ));
        }

        if (comp_settings.get_operator_settings().get_mutation_rate_scaler_settings().get_weight() > 0.0) {
            this->add_operator(std::make_shared<MutationRateScaler>(
                    i,
                    comp_settings.get_operator_settings().get_mutation_rate_scaler_settings().get_weight(),
                    comp_settings.get_operator_settings().get_mutation_rate_scaler_settings().get_scale()
                    ));
        }

        if (comp_settings.get_operator_settings().get_mean_population_size_scaler_settings().get_weight() > 0.0) {
            this->add_operator(std::make_shared<MeanPopulationSizeScaler>(
                    i,
                    comp_settings.get_operator_settings().get_mean_population_size_scaler_settings().get_weight(),
                    comp_settings.get_operator_settings().get_mean_population_size_scaler_settings().get_scale()
                    ));
        }

        if (comp_settings.get_operator_settings().get_relative_population_size_mixer_settings().get_weight() > 0.0) {
            this->add_operator(std::make_shared<RelativePopulationSizeMixer>(
                    i,
                    comp_settings.get_operator_settings().get_relative_population_size_mixer_settings().get_weight(),
                    comp_settings.get_operator_settings().get_relative_population_size_mixer_settings().get_scale()
                    ));
        }

        if (comp_settings.get_operator_settings().get_root_relative_population_size_mover_settings().get_weight() > 0.0) {
            this->add_operator(std::make_shared<RootRelativePopulationSizeMover>(
                    i,
                    comp_settings.get_operator_settings().get_root_relative_population_size_mover_settings().get_weight(),
                    comp_settings.get_operator_settings().get_root_relative_population_size_mover_settings().get_window()
                    ));
        }

        if (comp_settings.get_operator_settings().get_leaf_relative_population_size_mover_settings().get_weight() > 0.0) {
            this->add_operator(std::make_shared<LeafRelativePopulationSizeMover>(
                    i,
                    comp_settings.get_operator_settings().get_leaf_relative_population_size_mover_settings().get_weight(),
                    comp_settings.get_operator_settings().get_leaf_relative_population_size_mover_settings().get_window()
                    ));
        }

        if (comp_settings.get_operator_settings().get_freq_mover_settings().get_weight() > 0.0) {
            this->add_operator(std::make_shared<FreqMover>(
                    i,
                    comp_settings.get_operator_settings().get_freq_mover_settings().get_weight(),
                    comp_settings.get_operator_settings().get_freq_mover_settings().get_window()
                    ));
        }
    }
}

void OperatorSchedule::add_operator(std::shared_ptr<OperatorInterface> o) {
    this->operators_.push_back(o);
    this->total_weight_ += o->get_weight();
    this->cumulative_probs_.push_back(0.0);
    ECOEVOLITY_ASSERT(this->operators_.size() == this->cumulative_probs_.size());
    this->cumulative_probs_.at(0) = this->operators_.at(0)->get_weight() / this->total_weight_;
    for (unsigned int i = 1; i < this->operators_.size(); ++i) {
        this->cumulative_probs_.at(i) =
                (this->operators_.at(i)->get_weight() /
                this->total_weight_) + 
                this->cumulative_probs_.at(i - 1);
    }
}

OperatorInterface& OperatorSchedule::draw_operator(RandomNumberGenerator& rng) const {
    double u = rng.uniform_real();
    for (unsigned int i = 0; i < this->cumulative_probs_.size(); ++i) {
        if (u <= this->cumulative_probs_.at(i)) {
            return *this->operators_.at(i);
        }
    }
    return *this->operators_.back();
}

OperatorInterface& OperatorSchedule::get_operator(unsigned int operator_index) const {
    return *this->operators_.at(operator_index);
}

OperatorInterface& OperatorSchedule::get_reversible_jump_operator() const {
    for (auto & op: this->operators_) {
        if (op->get_type() == OperatorInterface::OperatorTypeEnum::rj_operator) {
            return *op;
        }
    }
    ECOEVOLITY_ASSERT(true == false);
    return *this->operators_.at(0);
}

std::vector< std::shared_ptr<OperatorInterface> > OperatorSchedule::get_time_operators() const {
    std::vector< std::shared_ptr<OperatorInterface> > ops;
    for (std::shared_ptr<OperatorInterface> op: this->operators_) {
        if (op->get_type() == OperatorInterface::OperatorTypeEnum::time_operator) {
            ops.push_back(op);
        }
    }
    return ops;
}

std::vector< std::shared_ptr<OperatorInterface> > OperatorSchedule::get_time_operators(
        int tree_index) const {
    std::vector< std::shared_ptr<OperatorInterface> > ops;
    for (std::shared_ptr<OperatorInterface> op: this->operators_) {
        if ((op->get_type() == OperatorInterface::OperatorTypeEnum::time_operator) &&
                (op->get_tree_index() == tree_index)) {
            ops.push_back(op);
        }
    }
    return ops;
}

std::vector< std::shared_ptr<OperatorInterface> > OperatorSchedule::get_multivariate_time_operators() const {
    std::vector< std::shared_ptr<OperatorInterface> > ops;
    for (std::shared_ptr<OperatorInterface> op: this->operators_) {
        if (op->get_type() == OperatorInterface::OperatorTypeEnum::multivariate_time_operator) {
            ops.push_back(op);
        }
    }
    return ops;
}

std::vector< std::shared_ptr<OperatorInterface> > OperatorSchedule::get_tree_operators() const {
    std::vector< std::shared_ptr<OperatorInterface> > ops;
    for (std::shared_ptr<OperatorInterface> op: this->operators_) {
        if (op->get_type() == OperatorInterface::OperatorTypeEnum::tree_operator) {
            ops.push_back(op);
        }
    }
    return ops;
}

std::vector< std::shared_ptr<OperatorInterface> > OperatorSchedule::get_tree_operators(
        int tree_index) const {
    std::vector< std::shared_ptr<OperatorInterface> > ops;
    for (std::shared_ptr<OperatorInterface> op: this->operators_) {
        if ((op->get_type() == OperatorInterface::OperatorTypeEnum::tree_operator) &&
                (op->get_tree_index() == tree_index)) {
            ops.push_back(op);
        }
    }
    return ops;
}

double OperatorSchedule::calc_delta(const Operator& op, double log_alpha) {
    if ((this->get_auto_optimize_delay_count() < this->get_auto_optimize_delay()) ||
            (! this->auto_optimize_)) {
        ++this->auto_optimize_delay_count_;
        return 0.0;
    }
    double target = op.get_target_acceptance_probability();
    double count = (op.get_number_rejected_for_correction() +
                    op.get_number_accepted_for_correction() +
                    1.0);
    double delta_p = ((1.0 / count) * (std::exp(std::min(log_alpha, 0.0)) - target));
    double mx = std::numeric_limits<double>::max();
    if ((delta_p > -mx) && (delta_p < mx)) {
        return delta_p;
    }
    return 0.0;
}

double OperatorSchedule::get_total_weight() const {
    return this->total_weight_;
}
unsigned int OperatorSchedule::get_auto_optimize_delay_count() const {
    return this->auto_optimize_delay_count_;
}
unsigned int OperatorSchedule::get_auto_optimize_delay() const {
    return this->auto_optimize_delay_;
}
void OperatorSchedule::set_auto_optimize_delay(unsigned int delay) {
    this->auto_optimize_delay_ = delay;
}

void OperatorSchedule::write_operator_rates(std::ostream& out) const {
    OperatorInterface& op = this->get_operator(0);
    out << op.header_string();
    out << op.to_string(*this);
    for (unsigned int i = 1; i < this->operators_.size(); ++i) {
        out << this->get_operator(i).to_string(*this);
    }
    out << std::flush;
}

bool OperatorSchedule::auto_optimizing() const {
    return this->auto_optimize_;
}
void OperatorSchedule::turn_on_auto_optimize() {
    this->auto_optimize_ = true;
}
void OperatorSchedule::turn_off_auto_optimize() {
    this->auto_optimize_ = false;
}

EcoevolityOptions::ModelOperator OperatorSchedule::get_model_operator_type() const {
    for (auto op : this->operators_) {
        if (op->get_name() == "DirichletProcessGibbsSampler") {
            return EcoevolityOptions::ModelOperator::gibbs_dpp;
        }
        if (op->get_name() == "PitmanYorProcessGibbsSampler") {
            return EcoevolityOptions::ModelOperator::gibbs_pyp;
        }
        if (op->get_type() == OperatorInterface::OperatorTypeEnum::rj_operator) {
            return EcoevolityOptions::ModelOperator::rj;
        }
    }
    return EcoevolityOptions::ModelOperator::none;
}

bool OperatorSchedule::sampling_models() const {
    if (this->get_model_operator_type() == EcoevolityOptions::ModelOperator::none) {
        return false;
    }
    return true;
}
