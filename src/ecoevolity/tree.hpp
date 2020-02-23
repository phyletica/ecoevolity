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

#include "basetree.hpp"
#include "data.hpp"
#include "node.hpp"
#include "likelihood.hpp"
#include "parameter.hpp"
#include "probability.hpp"
#include "error.hpp"
#include "assert.hpp"


class BasePopulationTree : public BaseTree<PopulationNode> {
    protected:
        BiallelicData data_;
        std::shared_ptr<ContinuousProbabilityDistribution> population_size_prior_ = std::make_shared<GammaDistribution>(1.0, 0.001);
        double ploidy_ = 2.0;
        std::shared_ptr<PositiveRealParameter> freq_1_ = std::make_shared<PositiveRealParameter>(
                std::make_shared<BetaDistribution>(1.0, 1.0),
                0.5);
        std::shared_ptr<PositiveRealParameter> mutation_rate_ = std::make_shared<PositiveRealParameter>(
                1.0,
                true);
        LogProbabilityDensity log_likelihood_correction_ = LogProbabilityDensity(0.0);
        bool likelihood_correction_was_calculated_ = false;
        bool constant_sites_removed_ = true;
        // int provided_number_of_constant_red_sites_ = -1;
        // int provided_number_of_constant_green_sites_ = -1;
        // bool use_removed_constant_site_counts_ = false;
        bool population_sizes_are_constrained_ = false;
        bool state_frequencies_are_constrained_ = false;
        bool is_dirty_ = true;

        // Vectors for storing unique allele counts and associated weights.
        // These are used for calculating the likelihood correction term for
        // constant site patterns. It is a bit weird to store data here, but
        // it's cheaper than calling 'this->data_.get_unique_allele_counts()'
        // every time the likelihood needs to be calculated.
        // These vectors are populated in 'init' method.
        std::vector< std::vector<unsigned int> > unique_allele_counts_;
        std::vector<unsigned int> unique_allele_count_weights_;

        // methods
        void process_and_vet_initialized_data(
                bool strict_on_constant_sites = false,
                bool strict_on_missing_sites = false,
                bool strict_on_triallelic_sites = true);

        void init_tree();

        // bool constant_site_counts_were_provided();
        void calculate_likelihood_correction();

        double calculate_log_binomial(
                unsigned int red_allele_count,
                unsigned int allele_count) const;

        void set_population_sizes(
                std::shared_ptr<PopulationNode> node,
                const std::vector<double> & sizes);
        
        void get_population_sizes(
                std::shared_ptr<PopulationNode> node,
                std::vector<double> & sizes) const;

        void update_unique_allele_counts();

    public:
        BasePopulationTree() { }
        BasePopulationTree(
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
                );

        BasePopulationTree(
                std::shared_ptr<PopulationNode> root,
                unsigned int number_of_loci,
                unsigned int length_of_loci,
                bool validate_data = false);

        BasePopulationTree(
                std::shared_ptr<PopulationNode> root);

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
                bool strict_on_triallelic_sites = true,
                double ploidy = 2.0,
                bool store_seq_loci_info = false
                );

        bool has_seq_loci_info() const {
            return this->data_.has_seq_loci_info();
        }

        void fold_patterns();

        bool constant_sites_removed() const {
            return this->constant_sites_removed_;
        }

        // int get_provided_number_of_constant_red_sites() const {
        //     return this->provided_number_of_constant_red_sites_;
        // }
        // int get_provided_number_of_constant_green_sites() const {
        //     return this->provided_number_of_constant_green_sites_;
        // }

        bool initialized() const {return (bool)this->root_;}

        const std::vector<std::string>& get_population_labels() const {
            return this->data_.get_population_labels();
        }

        const BiallelicData& get_data() const {
            return this->data_;
        }

        void set_data(const BiallelicData & data, bool constant_sites_removed);

        void set_ploidy(double ploidy) {
            this->ploidy_ = ploidy;
        }
        double get_ploidy() const {
            return this->ploidy_;
        }

        virtual double get_node_theta(const PopulationNode& node) const {
            return (2 * this->get_ploidy() *
                    node.get_population_size() *
                    this->get_mutation_rate());
        }
        double get_node_length_in_subs_per_site(const PopulationNode& node) const {
            return (node.get_length() * this->get_mutation_rate());
        }
        double get_node_height_in_subs_per_site(const PopulationNode& node) const {
            return (node.get_height() * this->get_mutation_rate());
        }

        void set_freq_1(double p);
        double get_freq_1() const;
        double get_freq_0() const;
        void store_freq_1();
        void restore_freq_1();
        double get_u() const;
        double get_v() const;

        void set_mutation_rate(double m);
        double get_mutation_rate() const;
        void store_mutation_rate();
        void restore_mutation_rate();

        bool is_dirty() const;
        void make_dirty();
        void make_clean();

        // void provide_number_of_constant_sites(
        //         unsigned int number_all_red,
        //         unsigned int number_all_green);

        std::shared_ptr<PositiveRealParameter> get_freq_1_parameter() const;

        void set_mutation_rate_parameter(std::shared_ptr<PositiveRealParameter> h);
        std::shared_ptr<PositiveRealParameter> get_mutation_rate_parameter() const;

        virtual void set_root_population_size(double size);
        virtual void set_all_population_sizes(double size);
        virtual unsigned int scale_all_population_sizes(double scale);
        virtual unsigned int scale_root_population_size(double scale);
        virtual double get_root_population_size() const;
        virtual std::shared_ptr<PositiveRealParameter> get_root_population_size_parameter() const;

        std::vector<double> get_population_sizes() const;
        std::vector< std::shared_ptr<PositiveRealParameter> > get_pointers_to_population_sizes() const;

        double get_mean_population_size() const;
        virtual void set_mean_population_size(double size);

        double get_leaf_mean_population_size() const;

        std::vector<double> get_population_sizes_as_proportions() const;
        std::vector<double> get_population_sizes_as_multipliers() const;
        virtual void set_population_sizes_as_proportions(const std::vector<double> & proportions);
        virtual void set_population_sizes(const std::vector<double> & sizes);

        double get_likelihood_correction(bool force = false);

        void compute_log_likelihood_and_prior(unsigned int nthreads = 1) {
            if (this->is_dirty()) {
                this->compute_log_likelihood(nthreads);
                ++this->number_of_likelihood_calculations_;
                this->compute_log_prior_density();
                this->make_clean();
            }
            return;
        }

        double compute_log_likelihood(unsigned int nthreads = 1);

        double get_derived_class_component_of_log_prior_density() const;
        double compute_log_prior_density_of_state_frequencies() const;
        double compute_log_prior_density_of_mutation_rate() const;
        virtual double compute_log_prior_density_of_population_sizes() const;

        void store_derived_class_parameters();
        virtual void store_all_population_sizes();
        void restore_derived_class_parameters();
        virtual void restore_all_population_sizes();

        virtual void set_population_size_prior(std::shared_ptr<ContinuousProbabilityDistribution> prior);
        virtual std::shared_ptr<ContinuousProbabilityDistribution> get_population_size_prior() const {
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

        virtual void fix_population_sizes() {
            this->root_->fix_all_population_sizes();
        }
        virtual void estimate_population_sizes() {
            this->root_->estimate_all_population_sizes();
        }
        virtual bool population_sizes_are_fixed() const {
            return this->root_->all_population_sizes_are_fixed();
        }

        virtual bool root_population_size_is_fixed() const {
            return this->root_->population_size_is_fixed();
        }

        virtual void constrain_population_sizes() {
            this->population_sizes_are_constrained_ = true;
            this->root_->set_all_population_size_parameters();
        }
        virtual bool population_sizes_are_constrained() const {
            return this->population_sizes_are_constrained_;
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
        bool mutation_rate_is_fixed() const {
            return this->mutation_rate_->is_fixed();
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

        void simulate_gene_tree(
                const std::shared_ptr<PopulationNode> node,
                std::unordered_map<unsigned int, std::vector< std::shared_ptr<GeneTreeSimNode> > > & branch_lineages,
                const unsigned int pattern_index,
                RandomNumberGenerator & rng,
                const bool use_max_allele_counts = false) const;

        std::shared_ptr<GeneTreeSimNode> simulate_gene_tree(
                const unsigned int pattern_index,
                RandomNumberGenerator& rng,
                const bool use_max_allele_counts = false) const;

        static double coalesce_in_branch(
                std::vector< std::shared_ptr<GeneTreeSimNode> >& lineages,
                double population_size,
                RandomNumberGenerator& rng,
                double bottom_of_branch_height = 0.0,
                double top_of_branch_height = std::numeric_limits<double>::infinity(),
                unsigned int branch_index = 0
                );

        bool sample_pattern(
                RandomNumberGenerator& rng,
                const float singleton_sample_probability,
                const std::vector<unsigned int>& red_allele_counts,
                const std::vector<unsigned int>& allele_counts) const;

        BiallelicData simulate_biallelic_data_set(
                RandomNumberGenerator& rng,
                float singleton_sample_probability = 1.0,
                bool validate = true) const;

        BiallelicData simulate_linked_biallelic_data_set(
                RandomNumberGenerator& rng,
                float singleton_sample_probability,
                bool max_one_variable_site_per_locus = false,
                bool validate = true) const;

        std::pair<BiallelicData, unsigned int>
        simulate_complete_biallelic_data_set(
                RandomNumberGenerator& rng,
                unsigned int locus_size = 1,
                float singleton_sample_probability = 1.0,
                bool validate = true) const;

        std::pair<BiallelicData, unsigned int>
        simulate_data_set_max_one_variable_site_per_locus(
                RandomNumberGenerator& rng,
                unsigned int locus_size,
                float singleton_sample_probability,
                bool validate = true) const;

        std::pair<
                std::pair<std::vector<unsigned int>, std::vector<unsigned int> >,
                std::shared_ptr<GeneTreeSimNode> >
        simulate_biallelic_site(
                const unsigned int pattern_idx,
                RandomNumberGenerator& rng,
                const bool use_max_allele_counts = false) const;

        std::pair<std::vector<unsigned int>, std::vector<unsigned int> >
        simulate_biallelic_site(
                std::shared_ptr<GeneTreeSimNode> gene_tree,
                RandomNumberGenerator& rng) const;

        std::pair<std::vector<unsigned int>, std::vector<unsigned int> >
        simulate_biallelic_site_sans_missing(
                std::shared_ptr<GeneTreeSimNode> gene_tree,
                const std::vector<unsigned int> & site_allele_counts,
                RandomNumberGenerator& rng) const;

        void write_data_summary(
                std::ostream& out,
                unsigned int indent_level = 0) const {
            this->data_.write_summary(out, indent_level);
        }

};



// TODO: PopulationTree is misnomer for this class and should be changed.
// This name made sense when this class used to be at the very base of the tree
// class hierarchy. Now it is an intermediate class in the hierarchy, an
// intermediate leading to the comparison tree classes (only one or two tips).
class PopulationTree : public BasePopulationTree {

    protected:
        std::shared_ptr<ContinuousProbabilityDistribution> root_node_height_prior_ = std::make_shared<ExponentialDistribution>(100.0);
        std::shared_ptr<DirichletDistribution> population_size_multiplier_prior_;
        bool population_size_multipliers_are_fixed_ = false;
        bool mean_population_size_is_fixed_ = false;

        void update_root_population_size() { return; }
        void update_relative_root_population_size() { return; }

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
                    store_seq_loci_info) { }

        PopulationTree(
                std::shared_ptr<PopulationNode> root,
                unsigned int number_of_loci = 10000,
                unsigned int length_of_loci = 1,
                bool validate_data = false
                ) : BasePopulationTree(
                    root,
                    number_of_loci,
                    length_of_loci,
                    validate_data) { }

        virtual bool using_population_size_multipliers() const {
            return false;
        }

        virtual bool using_relative_root_population_size() const {
            return false;
        }

        // These are overloaded by RelativeRootPopulationTree
        bool relative_root_population_size_is_fixed() const { return false; }
        void fix_relative_root_population_size() { return; }
        void estimate_relative_root_population_size() { return; }


        void set_root_node_height_prior(std::shared_ptr<ContinuousProbabilityDistribution> prior) {
            this->root_node_height_prior_ = prior;
            this->root_->set_all_node_height_priors(prior);
        }
        std::shared_ptr<ContinuousProbabilityDistribution> get_root_node_height_prior() const {
            return this->root_node_height_prior_;
        }
        void set_node_height_prior(std::shared_ptr<ContinuousProbabilityDistribution> prior) {
            this->set_root_node_height_prior(prior);
        }
        std::shared_ptr<ContinuousProbabilityDistribution> get_node_height_prior() const {
            return this->get_root_node_height_prior();
        }

        void set_root_height_parameter(std::shared_ptr<PositiveRealParameter> h);
        std::shared_ptr<PositiveRealParameter> get_root_height_parameter() const;


        void store_root_height();
        void restore_root_height();

        virtual void set_relative_root_population_size_prior(std::shared_ptr<ContinuousProbabilityDistribution> prior) { return; }

        virtual double get_relative_root_population_size() const {
            return this->get_root_population_size() / this->get_leaf_mean_population_size();
        }

        virtual void fix_population_sizes() {
            this->fix_mean_population_size();
            this->fix_population_size_multipliers();
            this->root_->fix_all_population_sizes();
        }
        virtual void estimate_population_sizes() {
            this->estimate_mean_population_size();
            this->estimate_population_size_multipliers();
            this->root_->estimate_all_population_sizes();
        }

        void fix_population_size_multipliers() {
            this->population_size_multipliers_are_fixed_ = true;
            this->make_dirty();
        }
        void estimate_population_size_multipliers() {
            this->population_size_multipliers_are_fixed_ = false;
            this->make_dirty();
        }
        bool population_size_multipliers_are_fixed() const {
            return this->population_size_multipliers_are_fixed_;
        }

        void fix_mean_population_size() {
            this->mean_population_size_is_fixed_ = true;
            this->make_dirty();
        }
        void estimate_mean_population_size() {
            this->mean_population_size_is_fixed_ = false;
            this->make_dirty();
        }
        bool mean_population_size_is_fixed() const {
            return this->mean_population_size_is_fixed_;
        }

        void set_population_size_multiplier_prior(std::shared_ptr<DirichletDistribution> prior);
        std::shared_ptr<DirichletDistribution> get_population_size_multiplier_prior() const{
            return this->population_size_multiplier_prior_;
        }

        void set_population_sizes(
                const std::vector<double> & sizes);
        
        void set_population_sizes_as_proportions(const std::vector<double> & proportions);
        
        void set_mean_population_size(double size);

        void store_parameters();
        void restore_parameters();

        // Override to old comparison prior behavior
        double compute_log_prior_density();
        // Override this method to old node height prior behavior
        double compute_log_prior_density_of_node_heights() const;

        // TODO: This PopulationTree hierarchy of classes is messy. The problem
        // is that each derived class has its own subset of methods in addition
        // to the base class methods. Thus, I can't simply use "virtual ... =
        // 0;" here, because some of the derived classes would be left with
        // invalid methods. Declaring methods here that throw errors if not
        // overloaded works for now.
        // Methods to be overloaded
        virtual void set_child_population_size(unsigned int child_index, double size) {
            throw EcoevolityError("set_child_population_size called from PopulationTree");
        }
        virtual double get_child_population_size(unsigned int child_index) const {
            throw EcoevolityError("get_child_population_size called from PopulationTree");
        }
        virtual std::shared_ptr<PositiveRealParameter> get_child_population_size_parameter(
                unsigned int child_index) const {
            throw EcoevolityError("get_child_population_size_parameter called from PopulationTree");
        }
};


class RelativeRootPopulationTree : public PopulationTree {

    protected:
        std::shared_ptr<PositiveRealParameter> relative_root_population_size_ = std::make_shared<PositiveRealParameter>(
                std::make_shared<GammaDistribution>(10.0, 0.1),
                1.0);
        bool leaf_population_sizes_are_fixed_ = false;

        void update_root_population_size();
        void update_relative_root_population_size();

    public:
        RelativeRootPopulationTree() { }
        RelativeRootPopulationTree(
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
                ) : PopulationTree(
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
                    store_seq_loci_info) { }

        bool using_relative_root_population_size() const {
            return true;
        }

        bool relative_root_population_size_is_fixed() const {
            return this->relative_root_population_size_->is_fixed();
        }
        void fix_relative_root_population_size() {
            this->relative_root_population_size_->fix();
        }
        void estimate_relative_root_population_size() {
            this->relative_root_population_size_->estimate();
            this->root_->estimate_population_size();
        }

        void set_root_population_size(double size);
        double get_root_population_size() const;
        void set_all_population_sizes(double size);
        unsigned int scale_all_population_sizes(double scale);
        unsigned int scale_root_population_size(double scale);
        void set_mean_population_size(double size);

        double get_relative_root_population_size() const {
            return this->relative_root_population_size_->get_value();
        }

        void set_population_sizes_as_proportions(
                const std::vector<double> & proportions);
        void set_population_sizes(const std::vector<double> & sizes);

        void set_relative_root_population_size_prior(
                std::shared_ptr<ContinuousProbabilityDistribution> prior);

        void fix_population_sizes() {
            PopulationTree::fix_population_sizes();
            if (! this->relative_root_population_size_is_fixed()) {
                this->root_->estimate_population_size();
            }
            this->leaf_population_sizes_are_fixed_ = true;
        }
        void estimate_population_sizes() {
            PopulationTree::estimate_population_sizes();
            this->leaf_population_sizes_are_fixed_ = false;
        }
        virtual void constrain_population_sizes() {
            PopulationTree::constrain_population_sizes();
            this->relative_root_population_size_->set_value(1.0);
            this->relative_root_population_size_->fix();
        }

        void store_all_population_sizes();
        void restore_all_population_sizes();

        double compute_log_prior_density_of_population_sizes() const;
};


class DirichletPopulationTree: public PopulationTree {

    public:
        DirichletPopulationTree() { }
        DirichletPopulationTree(
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
                );

        bool using_population_size_multipliers() const {
            return true;
        }

        void constrain_population_sizes() {
            throw EcoevolityError("This method is not supported; multipliers should be set and fixed");
        }
        bool population_sizes_are_constrained() const {
            throw EcoevolityError("This method is not supported; check if multipliers are fixed");
        }

        double compute_log_prior_density_of_population_sizes() const;
};


class ComparisonPopulationTree: public PopulationTree {

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
                bool strict_on_triallelic_sites = true,
                double ploidy = 2.0,
                bool store_seq_loci_info = false
                );
        ComparisonPopulationTree(
                const ComparisonSettings& settings,
                RandomNumberGenerator& rng,
                bool strict_on_constant_sites = false,
                bool strict_on_missing_sites = false,
                bool strict_on_triallelic_sites = true,
                bool store_seq_loci_info = false
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
                bool strict_on_triallelic_sites = true,
                double ploidy = 2.0,
                bool store_seq_loci_info = false
                );

        void set_child_population_size(unsigned int child_index, double size);
        double get_child_population_size(unsigned int child_index) const;
        std::shared_ptr<PositiveRealParameter> get_child_population_size_parameter(
                unsigned int child_index) const;

        void store_all_heights() {
            this->store_root_height();
        }
        void restore_all_heights() {
            this->restore_root_height();
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

        void draw_from_prior(RandomNumberGenerator& rng);

};

class ComparisonRelativeRootPopulationTree: public RelativeRootPopulationTree {

    public:
        ComparisonRelativeRootPopulationTree() { }
        ComparisonRelativeRootPopulationTree(
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
                );
        ComparisonRelativeRootPopulationTree(
                const RelativeRootComparisonSettings& settings,
                RandomNumberGenerator& rng,
                bool strict_on_constant_sites = false,
                bool strict_on_missing_sites = false,
                bool strict_on_triallelic_sites = true,
                bool store_seq_loci_info = false
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
                bool strict_on_triallelic_sites = true,
                double ploidy = 2.0,
                bool store_seq_loci_info = false
                );

        void set_child_population_size(unsigned int child_index, double size);
        double get_child_population_size(unsigned int child_index) const;
        std::shared_ptr<PositiveRealParameter> get_child_population_size_parameter(
                unsigned int child_index) const;

        void store_all_heights() {
            this->store_root_height();
        }
        void restore_all_heights() {
            this->restore_root_height();
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

        void draw_from_prior(RandomNumberGenerator& rng);

};

class ComparisonDirichletPopulationTree: public ComparisonPopulationTree {

    public:
        ComparisonDirichletPopulationTree() { }
        ComparisonDirichletPopulationTree(
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
                );
        ComparisonDirichletPopulationTree(
                const DirichletComparisonSettings& settings,
                RandomNumberGenerator& rng,
                bool strict_on_constant_sites = false,
                bool strict_on_missing_sites = false,
                bool strict_on_triallelic_sites = true,
                bool store_seq_loci_info = false
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
                bool strict_on_triallelic_sites = true,
                double ploidy = 2.0,
                bool store_seq_loci_info = false
                );

        bool using_population_size_multipliers() const {
            return true;
        }

        void constrain_population_sizes() {
            throw EcoevolityError("This method is not supported; multipliers should be set and fixed");
        }
        bool population_sizes_are_constrained() const {
            throw EcoevolityError("This method is not supported; check if multipliers are fixed");
        }

        double compute_log_prior_density_of_population_sizes() const;

        void write_state_log_header(std::ostream& out,
                bool include_event_index,
                const std::string& delimiter = "\t") const;
        void log_state(std::ostream& out,
                unsigned int event_index,
                const std::string& delimiter = "\t") const;
        void log_state(std::ostream& out,
                const std::string& delimiter = "\t") const;

        void draw_from_prior(RandomNumberGenerator& rng);
};

#endif
