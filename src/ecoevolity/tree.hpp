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
#include "likelihood.hpp"
#include "parameter.hpp"
#include "probability.hpp"
#include "error.hpp"
#include "assert.hpp"


class PopulationTree {
    protected:
        BiallelicData data_;
        std::shared_ptr<PopulationNode> root_;
        std::shared_ptr<ContinuousProbabilityDistribution> node_height_prior_ = std::make_shared<ExponentialDistribution>(100.0);
        std::shared_ptr<ContinuousProbabilityDistribution> population_size_prior_ = std::make_shared<GammaDistribution>(1.0, 0.001);
        double ploidy_ = 2.0;
        std::shared_ptr<PositiveRealParameter> freq_1_ = std::make_shared<PositiveRealParameter>(
                std::make_shared<BetaDistribution>(1.0, 1.0),
                0.5);
        std::shared_ptr<PositiveRealParameter> mutation_rate_ = std::make_shared<PositiveRealParameter>(
                1.0,
                true);
        LogProbabilityDensity log_likelihood_ = LogProbabilityDensity(0.0);
        LogProbabilityDensity log_likelihood_correction_ = LogProbabilityDensity(0.0);
        LogProbabilityDensity log_prior_density_ = LogProbabilityDensity(0.0);
        bool likelihood_correction_was_calculated_ = false;
        ProbabilityDensity all_green_pattern_likelihood_ = ProbabilityDensity(0.0);
        ProbabilityDensity all_red_pattern_likelihood_ = ProbabilityDensity(0.0);
        bool constant_sites_removed_ = true;
        int provided_number_of_constant_red_sites_ = -1;
        int provided_number_of_constant_green_sites_ = -1;
        // bool use_removed_constant_site_counts_ = false;
        bool population_sizes_are_constrained_ = false;
        bool state_frequencies_are_constrained_ = false;
        bool is_dirty_ = true;
        bool ignore_data_ = false;
        unsigned int number_of_likelihood_calculations_ = 0;

        // methods
        void init_tree();
        bool constant_site_counts_were_provided();
        void calculate_likelihood_correction();

        double calculate_log_binomial(
                unsigned int red_allele_count,
                unsigned int allele_count) const;


    public:
        PopulationTree() { }
        PopulationTree(
                std::string path, 
                char population_name_delimiter = ' ',
                bool population_name_is_prefix = true,
                bool genotypes_are_diploid = true,
                bool markers_are_dominant = false,
                bool constant_sites_removed = true,
                bool validate = true,
                bool strict_on_constant_sites = false,
                bool strict_on_missing_sites = false,
                double ploidy = 2.0
                );
        //~PopulationTree () { delete this->root_; }

        void init(
                std::string path, 
                char population_name_delimiter = ' ',
                bool population_name_is_prefix = true,
                bool genotypes_are_diploid = true,
                bool markers_are_dominant = false,
                bool constant_sites_removed = true,
                bool validate = true,
                bool strict_on_constant_sites = false,
                bool strict_on_missing_sites = false,
                double ploidy = 2.0
                );

        void fold_patterns();

        bool constant_sites_removed() const {
            return this->constant_sites_removed_;
        }

        int get_provided_number_of_constant_red_sites() const {
            return this->provided_number_of_constant_red_sites_;
        }
        int get_provided_number_of_constant_green_sites() const {
            return this->provided_number_of_constant_green_sites_;
        }

        bool initialized() const {return (bool)this->root_;}
        const PopulationNode& get_root() const {return *this->root_;}
        PopulationNode& get_mutable_root() const {return *this->root_;}
        void set_root_height(double height);
        void update_root_height(double height);
        double get_root_height() const;
        void store_root_height();
        void restore_root_height();

        void set_root_height_parameter(std::shared_ptr<PositiveRealParameter> h);
        std::shared_ptr<PositiveRealParameter> get_root_height_parameter() const;

        unsigned int get_degree_of_root() const {
            return this->root_->degree();
        }

        const std::vector<std::string>& get_population_labels() const {
            return this->data_.get_population_labels();
        }

        const BiallelicData& get_data() const {
            return this->data_;
        }

        void set_ploidy(double ploidy) {
            this->ploidy_ = ploidy;
        }
        double get_ploidy() const {
            return this->ploidy_;
        }

        double get_node_theta(const PopulationNode& node) const {
            return (2 * this->get_ploidy() *
                    node.get_population_size() *
                    this->get_mutation_rate());
        }
        double get_node_length_in_subs_per_site(const PopulationNode& node) const {
            return (node.get_length() * this->get_mutation_rate());
        }

        void set_freq_1(double p);
        void update_freq_1(double p);
        double get_freq_1() const;
        double get_freq_0() const;
        void store_freq_1();
        void restore_freq_1();
        double get_u() const;
        double get_v() const;

        void set_mutation_rate(double m);
        void update_mutation_rate(double m);
        double get_mutation_rate() const;
        void store_mutation_rate();
        void restore_mutation_rate();

        bool is_dirty() const;
        void make_dirty();
        void make_clean();

        void ignore_data() {
            this->ignore_data_ = true;
        }
        void use_data() {
            this->ignore_data_ = false;
        }
        bool ignoring_data() const {
            return this->ignore_data_;
        }

        void provide_number_of_constant_sites(
                unsigned int number_all_red,
                unsigned int number_all_green);

        std::shared_ptr<PositiveRealParameter> get_freq_1_parameter() const;

        void set_mutation_rate_parameter(std::shared_ptr<PositiveRealParameter> h);
        std::shared_ptr<PositiveRealParameter> get_mutation_rate_parameter() const;

        void set_root_population_size(double size);
        void set_population_size(double size);
        double get_root_population_size() const;
        std::shared_ptr<PositiveRealParameter> get_root_population_size_parameter() const;

        double get_likelihood_correction(bool force = false);

        virtual void compute_log_likelihood_and_prior(unsigned int nthreads = 1) {
            if (this->is_dirty()) {
                this->compute_log_likelihood(nthreads);
                ++this->number_of_likelihood_calculations_;
                this->compute_log_prior_density();
                this->make_clean();
            }
            return;
        }

        double compute_log_likelihood(unsigned int nthreads = 1);

        void set_log_likelihood_value(double value) {
            this->log_likelihood_.set_value(value);
        }
        double get_log_likelihood_value() const;
        double get_stored_log_likelihood_value() const;

        virtual double compute_log_prior_density();
        double compute_log_prior_density_of_state_frequencies() const;
        double compute_log_prior_density_of_mutation_rate() const;
        double compute_log_prior_density_of_node_heights() const;
        double compute_log_prior_density_of_population_sizes() const;
        double get_log_prior_density_value() const;
        double get_stored_log_prior_density_value() const;

        void store_state();
        void store_likelihood();
        void store_prior_density();
        virtual void store_parameters();
        void store_all_population_sizes();
        virtual void store_all_heights();
        void restore_state();
        void restore_likelihood();
        void restore_prior_density();
        virtual void restore_parameters();
        void restore_all_population_sizes();
        virtual void restore_all_heights();

        void set_node_height_prior(std::shared_ptr<ContinuousProbabilityDistribution> prior);
        std::shared_ptr<ContinuousProbabilityDistribution> get_node_height_prior() const {
            return this->node_height_prior_;
        }

        void set_population_size_prior(std::shared_ptr<ContinuousProbabilityDistribution> prior);
        std::shared_ptr<ContinuousProbabilityDistribution> get_population_size_prior() const {
            return this->population_size_prior_;
        }

        void set_freq_1_prior(std::shared_ptr<ContinuousProbabilityDistribution> prior);
        std::shared_ptr<ContinuousProbabilityDistribution> get_freq_1_prior() const {
            return this->freq_1_->prior;
        }

        void set_mutation_rate_prior(std::shared_ptr<ContinuousProbabilityDistribution> prior);
        std::shared_ptr<ContinuousProbabilityDistribution> get_mutation_rate_prior() const {
            return this->mutation_rate_->prior;
        }

        void fix_population_sizes() {
            this->root_->fix_all_population_sizes();
        }
        void estimate_population_sizes() {
            this->root_->estimate_all_population_sizes();
        }
        bool population_sizes_are_fixed() const {
            return this->root_->all_population_sizes_are_fixed();
        }

        void fix_state_frequencies() {
            this->freq_1_->fix();
        }
        void estimate_state_frequencies() {
            if (this->state_frequencies_are_constrained_) {
                throw EcoevolityError("Cannot estimate constrained state frequencies");
            }
            this->freq_1_->estimate();
        }
        bool state_frequencies_are_fixed() const {
            return this->freq_1_->is_fixed();
        }

        void fix_mutation_rate() {
            this->mutation_rate_->fix();
        }
        void estimate_mutation_rate() {
            this->mutation_rate_->estimate();
        }
        bool mutation_rate_is_fixed() {
            return this->mutation_rate_->is_fixed();
        }

        void constrain_population_sizes() {
            this->population_sizes_are_constrained_ = true;
            this->root_->set_all_population_size_parameters();
        }
        bool population_sizes_are_constrained() const {
            return this->population_sizes_are_constrained_;
        }
        void constrain_state_frequencies() {
            this->state_frequencies_are_constrained_ = true;
            this->freq_1_->set_value(0.5);
            this->freq_1_->fix();
            this->make_dirty();
        }
        bool state_frequencies_are_constrained() const {
            return this->state_frequencies_are_constrained_;
        }

        unsigned int get_number_of_likelihood_calculations() {
            return this->number_of_likelihood_calculations_;
        }
        unsigned int get_leaf_node_count() const {
            return this->root_->get_leaf_node_count();
        }
};

class ComparisonPopulationTree: public PopulationTree {

    friend class ComparisonPopulationTreeCollection;

    public:
        ComparisonPopulationTree() { }
        ComparisonPopulationTree(
                std::string path, 
                char population_name_delimiter = ' ',
                bool population_name_is_prefix = true,
                bool genotypes_are_diploid = true,
                bool markers_are_dominant = false,
                bool constant_sites_removed = true,
                bool validate = true,
                bool strict_on_constant_sites = false,
                bool strict_on_missing_sites = false,
                double ploidy = 2.0
                );
        ComparisonPopulationTree(
                const ComparisonSettings& settings,
                RandomNumberGenerator& rng,
                bool strict_on_constant_sites = false,
                bool strict_on_missing_sites = false
                );
        void comparison_init(
                std::string path, 
                char population_name_delimiter = ' ',
                bool population_name_is_prefix = true,
                bool genotypes_are_diploid = true,
                bool markers_are_dominant = false,
                bool constant_sites_removed = true,
                bool validate = true,
                bool strict_on_constant_sites = false,
                bool strict_on_missing_sites = false,
                double ploidy = 2.0
                );

        void set_child_population_size(unsigned int child_index, double size);
        void update_child_population_size(unsigned int child_index, double size);
        double get_child_population_size(unsigned int child_index) const;
        void store_child_population_size(unsigned int child_index);
        void restore_child_population_size(unsigned int child_index);
        std::shared_ptr<PositiveRealParameter> get_child_population_size_parameter(
                unsigned int child_index) const;

        void set_height(double height) {this->set_root_height(height);}
        void update_height(double height) {this->update_root_height(height);}
        double get_height() const {return this->get_root_height();}
        void store_height() {this->store_root_height();}
        void restore_height() {this->restore_root_height();}
        void set_height_parameter(std::shared_ptr<PositiveRealParameter> h) {
            this->set_root_height_parameter(h);
        }
        std::shared_ptr<PositiveRealParameter> get_height_parameter() const {
            return this->get_root_height_parameter();
        }

        void store_all_heights() {
            this->store_height();
        }
        void restore_all_heights() {
            this->restore_height();
        }
        void store_parameters();
        void restore_parameters();

        double compute_log_prior_density();

        void write_state_log_header(std::ostream& out,
                bool include_event_index,
                const std::string& delimiter = "\t") const;
        void log_state(std::ostream& out,
                unsigned int event_index,
                const std::string& delimiter = "\t") const;
        void log_state(std::ostream& out,
                const std::string& delimiter = "\t") const;

        std::string get_state_header_string(
                const std::string& delimiter = "\t") const;
        std::string get_state_string(
                const std::string& delimiter = "\t",
                unsigned int precision = 18) const;

        std::shared_ptr<GeneTreeSimNode> simulate_gene_tree(
                const unsigned int pattern_index,
                RandomNumberGenerator& rng) const;

        static double coalesce_in_branch(
                std::vector< std::shared_ptr<GeneTreeSimNode> >& lineages,
                double population_size,
                RandomNumberGenerator& rng,
                double bottom_of_branch_height = 0.0,
                double top_of_branch_height = std::numeric_limits<double>::infinity()
                );

        BiallelicData simulate_biallelic_data_set(
                RandomNumberGenerator& rng,
                bool validate = true) const;

        std::pair<
                std::pair<std::vector<unsigned int>, std::vector<unsigned int> >,
                std::shared_ptr<GeneTreeSimNode> >
        simulate_biallelic_site(
                const unsigned int pattern_idx,
                RandomNumberGenerator& rng) const;

        void write_data_summary(
                std::ostream& out,
                unsigned int indent_level = 0) const {
            this->data_.write_summary(out, indent_level);
        }

        void draw_from_prior(RandomNumberGenerator& rng);

};

#endif
