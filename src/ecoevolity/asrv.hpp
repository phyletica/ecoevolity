#pragma once

#include <algorithm>
#include <vector>

namespace ecoevolity {

    class Model;

    class ASRV {

        friend class Model;

        public:
            typedef std::shared_ptr<ASRV>   SharedPtr;

                                            ASRV();
                                            ASRV(
                                                unsigned num_categ,
                                                std::shared_ptr<PositiveRealParameter> one_over_shape,
                                                std::shared_ptr<ProportionParameter> prop_invar_sites);
                                            ~ASRV();

            unsigned                        get_num_categ() const;

            double                          get_shape() const;
            double                          get_one_over_shape() const;
            // void                            set_one_over_shape(double v);
            // void                            store_one_over_shape();
            // void                            restore_one_over_shape();
            // void                            fix_one_over_shape();
            // void                            estimate_one_over_shape();
            bool                            one_over_shape_is_fixed() const;

            bool                            is_invariable_sites_model() const;
            double                          get_proportion_invariable_sites() const;
            // void                            set_proportion_invariable_sites(double p);
            // void                            store_proportion_invariable_sites();
            // void                            restore_proportion_invariable_sites();
            // void                            fix_proportion_invariable_sites();
            // void                            estimate_proportion_invariable_sites();
            bool                            proportion_invariable_sites_is_fixed() const;

            void                            clear();

        private:

            std::shared_ptr<PositiveRealParameter>  _one_over_shape;
            unsigned                                _num_categ;
            std::shared_ptr<ProportionParameter>    _proportion_invariable;

            void                            set_one_over_shape_parameter(
                                                std::shared_ptr<PositiveRealParameter> p);
            void                            set_one_over_shape_prior(
                                                std::shared_ptr<ContinuousProbabilityDistribution> prior);

            void                            set_proportion_invariable_sites_parameter(
                                                std::shared_ptr<ProportionParameter> p);
            void                            set_proportion_invariable_sites_prior(
                                                std::shared_ptr<ContinuousProbabilityDistribution> prior);
        };

    inline ASRV::ASRV() {
        //std::cout << "Constructing a Model" << std::endl;
        clear();
    }
    inline void ASRV::ASRV(
            unsigned num_categ,
            std::shared_ptr<PositiveRealParameter> one_over_shape,
            std::shared_ptr<ProportionParameter> prop_invar_sites)
    : ASRV() {
        this->_num_categ = num_categ;
        this->set_one_over_shape_parameter(one_over_shape);
        this->set_proportion_invariable_sites_parameter(prop_invar_sites);
    }


    inline ASRV::~ASRV() {
        //std::cout << "Destroying a Model" << std::endl;
    }

    inline void ASRV::clear() {
        // Fix _one_over_shape to 1.0
        this->_one_over_shape = std::make_shared<PositiveRealParameter>(1.0, true);;
        this->_num_categ = 1;
        // Fix _proportion_invariable to zero
        this->_proportion_invariable = std::make_shared<ProportionParameter>(0.0, true);
    }

    inline double ASRV::get_shape() const {
        return 1.0/_one_over_shape->get_value();
    }

    inline double ASRV::get_one_over_shape() const {
        return _one_over_shape->get_value();
    }

    // inline void ASRV::set_one_over_shape(double p) {
    //     if (_one_over_shape->is_fixed()) {
    //         // throw std::exception("Tried to set fixed ASRV one over shape");
    //         return;
    //     }
    //     _one_over_shape->set_value(p);
    // }

    // inline void ASRV::store_one_over_shape() {
    //     _one_over_shape->store();
    // }

    // inline void ASRV::restore_one_over_shape() {
    //     _one_over_shape->restore();
    // }

    // inline void ASRV::fix_one_over_shape() {
    //     _one_over_shape->fix();
    // }

    // inline void ASRV::estimate_one_over_shape() {
    //     _one_over_shape->estimate();
    // }

    inline bool ASRV::one_over_shape_is_fixed() const {
        return _one_over_shape->is_fixed();
    }

    inline void ASRV::set_one_over_shape_parameter(
            std::shared_ptr<PositiveRealParameter> p) {
        _one_over_shape = p;
    }

    inline void ASRV::set_one_over_shape_prior(
            std::shared_ptr<ContinuousProbabilityDistribution> prior) {
        _one_over_shape->set_prior(prior);;
    }

    inline unsigned ASRV::get_num_categ() const {
        return _num_categ;
    }

    inline double ASRV::get_proportion_invariable_sites() const {
        return _proportion_invariable->get_value();
    }

    // inline void ASRV::set_proportion_invariable_sites(double p) {
    //     if (_proportion_invariable->is_fixed()) {
    //         // throw std::exception("Tried to set fixed proportion invariable sites");
    //         return;
    //     }
    //     _proportion_invariable->set_value(p);
    // }

    // inline void ASRV::store_proportion_invariable_sites() {
    //     _proportion_invariable->store();
    // }

    // inline void ASRV::restore_proportion_invariable_sites() {
    //     _proportion_invariable->restore();
    // }

    // inline void ASRV::fix_proportion_invariable_sites() {
    //     _proportion_invariable->fix();
    // }

    // inline void ASRV::estimate_proportion_invariable_sites() {
    //     _proportion_invariable->estimate();
    // }

    inline bool ASRV::proportion_invariable_sites_is_fixed() const {
        return _proportion_invariable->is_fixed();
    }

    inline bool ASRV::is_invariable_sites_model() const {
        if ((_proportion_invariable->get_value() == 0.0) && (_proportion_invariable->is_fixed())) {
            return false;
        }
        return true;
    }

    inline void ASRV::set_proportion_invariable_sites_parameter(
            std::shared_ptr<ProportionParameter> p) {
        _proportion_invariable = p;
    }

    inline void ASRV::set_proportion_invariable_sites_prior(
            std::shared_ptr<ContinuousProbabilityDistribution> prior) {
        _proportion_invariable->set_prior(prior);;
    }
}
