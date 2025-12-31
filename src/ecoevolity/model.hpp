#pragma once

#include <algorithm>
#include <vector>
#include "qmatrix.hpp"
#include "asrv.hpp"
#include "parameter.hpp"
#include "partition.hpp"

namespace ecoevolity {

    class Model {

        public:
            typedef std::shared_ptr<Model>              SharedPtr;
            typedef std::vector<ASRV::SharedPtr>        asrv_vect_t;
            typedef std::vector<QMatrix::SharedPtr>     qmatrix_vect_t;
            typedef std::vector<unsigned>               subset_sizes_t;

                                            Model();
                                            ~Model();

            QMatrix::SharedPtr              get_q_matrix(unsigned subset_index) const;
            ASRV::SharedPtr                 get_asrv(unsigned subset_index) const;

            void                            init(
                                                const NucTreeAnalysisSettings::partition_settings_t & partition_settings,
                                                const PositiveRealParameterSettings & subset_rel_rate_settings,
                                                const subset_sizes_t & subset_sizes);

            void                            clear();

        private:

            qmatrix_vect_t                              _qmatrices;
            asrv_vect_t                                 _asrvs;
            ProportionParameterVector                   _subset_rel_rates;
            subset_sizes_t                              _subset_sizes;
            unsigned                                    _num_sites;

            std::vector< std::shared_ptr<ProportionParameterVector> >   _uniq_state_freq_params;
            std::vector< std::shared_ptr<ProportionParameterVector> >   _uniq_rate_matrix_params;
            std::vector< std::shared_ptr<PositiveRealParameter> >       _uniq_asrv_shape_params;
            std::vector<unsigned>                                       _uniq_asrv_num_cats;
            std::vector< std::shared_ptr<ProportionParameter> >         _uniq_asrv_pinvar_params;
        };

    inline Model::Model() {
        //std::cout << "Constructing a Model" << std::endl;
        clear();
    }

    inline Model::~Model() {
        //std::cout << "Destroying a Model" << std::endl;
    }

    inline void Model::clear() {
        this->_num_sites = 0;
        this->_subset_sizes.clear();
        this->_qmatrices->clear();
        this->_asrvs->clear();

        this->_subset_rel_rates.clear();

        std::vector<double> values = {1.0};
        bool fix = true;
        this->_subset_rel_rates = ProportionParameterVector(values, fix);

        this->_uniq_state_freq_params.clear();
        this->_uniq_rate_matrix_params.clear();
        this->_uniq_asrv_shape_params.clear();
        this->_uniq_asrv_pinvar_params.clear();

        this->_uniq_asrv_num_cats.clear();
    }

    inline void Model::init(
            // const NucTreeAnalysisSettings::param_settings_map_t & state_freq_settings, 
            // const NucTreeAnalysisSettings::param_settings_map_t & rate_matrix_settings, 
            // const NucTreeAnalysisSettings::gamma_settings_map_t & discrete_gamma_settings,
            // const NucTreeAnalysisSettings::param_settings_map_t & prop_invar_settings, 
            const NucTreeAnalysisSettings::partition_settings_t & partition_settings,
            const PositiveRealParameterSettings                 & subset_rel_rate_settings,
            const subset_sizes_t                                & subset_sizes) {

        this->clear();

        // Temp maps of subset names to parameters in case some subsets have
        // linked parameters
        std::map<std::string, std::shared_ptr<ProportionParameterVector> > subset_to_state_freqs;
        std::map<std::string, std::shared_ptr<ProportionParameterVector> > subset_to_rate_matrix;
        std::map<std::string, std::shared_ptr<ProportionParameter> > subset_to_prop_invar;
        std::map<std::string, std::shared_ptr<PositiveRealParameter> > subset_to_asrv_shape;
        std::map<std::string, std::shared_ptr<PositiveRealParameter> > subset_to_asrv_num_cats;

        for (const auto & subset_settings : partition_settings) {
            std::shared_ptr<ProportionParameterVector> state_freqs;
            if (! subset_settings.state_freq_link.empty()) {
                state_freqs = subset_to_state_freqs.at(subset_settings.state_freq_link);
            }
            else {
                state_freqs = std::make_shared<ProportionParameterVector>(
                    subset_settings.state_freq_settings,
                    rng);
                this->_uniq_state_freq_params.push_back(state_freqs);
            }
            if (! subset_settings.name.empty()) {
                subset_to_state_freqs[subset_settings.name] = state_freqs;
            }

            std::shared_ptr<ProportionParameterVector> rate_matrix;
            if (! subset_settings.rate_matrix_link.empty()) {
                rate_matrix = subset_to_rate_matrix.at(subset_settings.rate_matrix_link);
            }
            else {
                rate_matrix = std::make_shared<ProportionParameterVector>(
                    subset_settings.rate_matrix_settings,
                    rng);
                this->_uniq_rate_matrix_params.push_back(rate_matrix);
            }
            if (! subset_settings.name.empty()) {
                subset_to_rate_matrix[subset_settings.name] = rate_matrix;
            }

            std::shared_ptr<ProportionParameter> p_invar;
            if (! subset_settings.prop_invar_link.empty()) {
                p_invar = subset_to_prop_invar.at(subset_settings.prop_invar_link);
            }
            else {
                p_invar = std::make_shared<ProportionParameter>(
                    subset_settings.prop_invar_settings,
                    rng);
                this->_uniq_asrv_prop_invar_params.push_back(p_invar);
            }
            if (! subset_settings.name.empty()) {
                subset_to_prop_invar[subset_settings.name] = p_invar;
            }

            std::shared_ptr<PositiveRealParameter> asrv_shape;
            unsigned int num_cats;
            if (! subset_settings.discrete_gamma_link.empty()) {
                asrv_shape = subset_to_asrv_shape.at(subset_settings.discrete_gamma_link);
                num_cats = subset_to_asrv_num_cats.at(subset_settings.discrete_gamma_link);
            }
            else {
                asrv_shape = std::make_shared<PositiveRealParameter>(
                    subset_settings.discrete_gamma_settings.one_over_shape_settings,
                    rng);
                num_cats = subset_settings.discrete_gamma_settings.num_cats
                this->_uniq_asrv_shape_params.push_back(asrv_shape);
                this->_uniq_asrv_num_cats.push_back(num_cats);
            }
            if (! subset_settings.name.empty()) {
                subset_to_asrv_shape[subset_settings.name] = asrv_shape;
                subset_to_asrv_num_cats[subset_settings.name] = num_cats;
            }

            QMatrix::SharedPtr qmat = std::make_shared<QMatrixNucleotide>(
                    state_freqs,
                    rate_matrix);
            this->_qmatrices.push_back(qmat);

            ASRV::SharedPtr asrv_model = std::make_shared<ASRV>(
                    num_cats,
                    asrv_shape,
                    p_invar);
            this->_asrvs.push_back(asrv_model);
        }

        this->_subset_rel_rates = ProportionParameterVector(subset_rel_rate_settings, rng);

        this->_subset_sizes = subset_sizes;
        this->_num_sites = std::accumulate(this->_subset_sizes.begin(), this->_subset_sizes.end(), 0.0);
    }

    inline QMatrix::SharedPtr Model::get_q_matrix(unsigned subset_index) const {
        return _qmatrices.at(subset_index);
    }

    inline ASRV::SharedPtr Model::get_asrv(unsigned subset_index) const {
        return _asrvs.at(subset_index);
    }
}
