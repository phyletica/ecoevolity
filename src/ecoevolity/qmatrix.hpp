#pragma once

#include <algorithm>
#include <vector>
#include "error.hpp"

namespace ecoevolity {

    class Model;

    class QMatrix {

        friend class Model;

        public:
            ProportionParameterVector                           freq_xchg_t;
            typedef std::shared_ptr<ProportionParameterVector>  freq_xchg_ptr_t;
            typedef std::shared_ptr<QMatrix>                    SharedPtr;

                                                    QMatrix();
            virtual                                 ~QMatrix();

            virtual void                            clear() = 0;

            // virtual void                            set_equal_state_freqs(freq_xchg_ptr_t freq_ptr) = 0;
            // virtual void                            set_state_freqs(const freq_xchg_t & freq) = 0;
            virtual freq_xchg_ptr_t                 get_state_freqs_shared_ptr() = 0;
            virtual const double *                  get_state_freqs() const = 0;

            virtual bool                            is_fixed_state_freqs() const;

            // virtual void                            set_equal_exchangeabilities(freq_xchg_ptr_t xchg_ptr) = 0;
            // virtual void                            set_exchangeabilities(const freq_xchg_t & xchg) = 0;
            virtual freq_xchg_ptr_t                 get_exchangeabilities_shared_ptr() = 0;
            virtual const double *                  get_exchangeabilities() const = 0;

            virtual bool                            is_fixed_exchangeabilities() const;

        protected:

            virtual void                            set_state_freqs_shared_ptr(freq_xchg_ptr_t freq_ptr) = 0;
            virtual void                            set_exchangeabilities_shared_ptr(freq_xchg_ptr_t xchg) = 0;

            void                                    normalize_freqs_or_exchangeabilities(freq_xchg_ptr_t v);
    };

    inline QMatrix::QMatrix() {
        //std::cout << "Creating a QMatrix object" << std::endl;
    }

    inline QMatrix::~QMatrix() {
        //std::cout << "Destroying a QMatrix object" << std::endl;
    }

    inline void QMatrix::clear() { }

    inline void QMatrix::normalize_freqs_or_exchangeabilities(QMatrix::freq_xchg_ptr_t v) {
        // Be sure elements of v sum to 1.0 and assert that they are all positive
        double sum_v = std::accumulate(v->begin(), v->end(), 0.0);
        for (auto & x : *v) {
            assert(x > 0.0);
            x /= sum_v;
        }
    }

    class QMatrixNucleotide : public QMatrix {

        public:
                                        QMatrixNucleotide();
                                        ~QMatrixNucleotide();

            void                        clear();

            // void                        set_equal_state_freqs(freq_xchg_ptr_t freq_ptr);
            void                        set_state_freqs_shared_ptr(freq_xchg_ptr_t freq_ptr);
            // void                        set_state_freqs(const freq_xchg_t & freqs);
            freq_xchg_ptr_t             get_state_freqs_shared_ptr();
            const double *              get_state_freqs() const;

            bool                        is_fixed_state_freqs() const;

            // void                        set_equal_exchangeabilities(freq_xchg_ptr_t xchg_ptr);
            void                        set_exchangeabilities_shared_ptr(freq_xchg_ptr_t xchg_ptr);
            // void                        set_exchangeabilities(const freq_xchg_t & xchg);
            freq_xchg_ptr_t             get_exchangeabilities_shared_ptr();
            const double *              get_exchangeabilities() const;
            bool                        is_fixed_exchangeabilities() const;

        private:

            freq_xchg_ptr_t             _state_freqs;
            freq_xchg_ptr_t             _exchangeabilities;
    };

    inline QMatrixNucleotide::QMatrixNucleotide() {
        //std::cout << "Constructing a QMatrixNucleotide object" << std::endl;
        clear();
    }

    inline QMatrixNucleotide::QMatrixNucleotide(
            freq_xchg_ptr_t state_freqs,
            freq_xchg_ptr_t exchangeabilities)
    : QMatrixNucleotide() {
        this->set_state_freqs_shared_ptr(state_freqs);
        this->set_exchangeabilities_shared_ptr(exchangeabilities);
    }

    inline QMatrixNucleotide::~QMatrixNucleotide() {
        //std::cout << "Destroying a QMatrixNucleotide object" << std::endl;
    }

    inline void QMatrixNucleotide::clear() {
        QMatrix::clear();

        std::vector<double> values = {1/6.0, 1/6.0, 1/6.0, 1/6.0, 1/6.0, 1/6.0};
        std::vector<double> alphas = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
        bool fix = false;
        std::shared_ptr<DirichletDistribution> xchg_prior_ptr = std::make_shared<DirichletDistribution>(alphas);
        freq_xchg_ptr_t xchg_ptr = std::make_shared<QMatrix::freq_xchg_t>(xchg_prior_ptr, values, fix);
        this->set_exchangeabilities_shared_ptr(xchg_ptr);

        values = {0.25, 0.25, 0.25, 0.25};
        alphas = {1.0, 1.0, 1.0, 1.0,};
        fix = false;
        std::shared_ptr<DirichletDistribution> freq_prior_ptr = std::make_shared<DirichletDistribution>(alphas);
        freq_xchg_ptr_t freq_ptr = std::make_shared<QMatrix::freq_xchg_t>(freq_prior_ptr, values, fix);
        this->set_state_freqs_shared_ptr(freq_ptr);

    }

    inline bool QMatrixNucleotide::is_fixed_state_freqs() const {
        return _state_freqs->is_fixed();
    }

    inline bool QMatrixNucleotide::is_fixed_exchangeabilities() const {
        return _exchangeabilities->is_fixed();
    }

    inline QMatrix::freq_xchg_ptr_t QMatrixNucleotide::get_exchangeabilities_shared_ptr() {
        return _exchangeabilities;
    }

    inline QMatrix::freq_xchg_ptr_t QMatrixNucleotide::get_state_freqs_shared_ptr() {
        return _state_freqs;
    }

    inline const double * QMatrixNucleotide::get_exchangeabilities() const {
        return _exchangeabilities->get_raw_values();
    }

    inline const double * QMatrixNucleotide::get_state_freqs() const {
        return _state_freqs->get_raw_values();
    }

    // inline void QMatrixNucleotide::set_equal_exchangeabilities(QMatrix::freq_xchg_ptr_t xchg_ptr) {
    //     _exchangeabilities = xchg_ptr;
    //     _exchangeabilities->assign(6, 1.0/6.0);
    // }

    inline void QMatrixNucleotide::set_exchangeabilities_shared_ptr(QMatrix::freq_xchg_ptr_t xchg_ptr) {
        assert(xchg_ptr->get_values().size() == 6);
        _exchangeabilities = xchg_ptr;
    }

    // inline void QMatrixNucleotide::set_exchangeabilities(const QMatrix::freq_xchg_t & xchg) {
    //     assert(xchg.size() == 6);
    //     std::copy(xchg.begin(), xchg.end(), _exchangeabilities->begin());
    // }

    // inline void QMatrixNucleotide::set_equal_state_freqs(QMatrix::freq_xchg_ptr_t freq_ptr) {
    //     _state_freqs = freq_ptr;
    //     _state_freqs->assign(4, 0.25);
    // }

    inline void QMatrixNucleotide::set_state_freqs_shared_ptr(QMatrix::freq_xchg_ptr_t freq_ptr) {
        assert(freq_ptr->get_values().size() == 4);
        _state_freqs = freq_ptr;
    }

    // inline void QMatrixNucleotide::set_state_freqs(const QMatrix::freq_xchg_t & freqs) {
    //     assert(freqs.size() == 4);
    //     std::copy(freqs.begin(), freqs.end(), _state_freqs->begin());
    // }
}
