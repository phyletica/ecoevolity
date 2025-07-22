#pragma once

#include <algorithm>
#include <vector>
#include "error.hpp"

namespace ecoevolity {

    class QMatrix {

        public:
            typedef std::vector<double>             freq_xchg_t;
            typedef std::shared_ptr<freq_xchg_t>    freq_xchg_ptr_t;
            typedef std::shared_ptr<QMatrix>        SharedPtr;

                                                    QMatrix();
            virtual                                 ~QMatrix();

            virtual void                            clear() = 0;

            virtual void                            set_equal_state_freqs(freq_xchg_ptr_t freq_ptr) = 0;
            virtual void                            set_state_freqs_shared_ptr(freq_xchg_ptr_t freq_ptr) = 0;
            virtual void                            set_state_freqs(const freq_xchg_t & freq) = 0;
            virtual freq_xchg_ptr_t                 get_state_freqs_shared_ptr() = 0;
            virtual const double *                  get_state_freqs() const = 0;
            void                                    fix_state_freqs(bool is_fixed);
            bool                                    is_fixed_state_freqs() const;

            virtual void                            set_equal_exchangeabilities(freq_xchg_ptr_t xchg_ptr) = 0;
            virtual void                            set_exchangeabilities_shared_ptr(freq_xchg_ptr_t xchg) = 0;
            virtual void                            set_exchangeabilities(const freq_xchg_t & xchg) = 0;
            virtual freq_xchg_ptr_t                 get_exchangeabilities_shared_ptr() = 0;
            virtual const double *                  get_exchangeabilities() const = 0;
            void                                    fix_exchangeabilities(bool is_fixed);
            bool                                    is_fixed_exchangeabilities() const;

        protected:

            void                                    normalize_freqs_or_exchangeabilities(freq_xchg_ptr_t v);

            bool                                    _state_freqs_fixed;
            bool                                    _exchangeabilities_fixed;
    };

    inline QMatrix::QMatrix() {
        //std::cout << "Creating a QMatrix object" << std::endl;
    }

    inline QMatrix::~QMatrix() {
        //std::cout << "Destroying a QMatrix object" << std::endl;
    }

    inline void QMatrix::clear() {
        _state_freqs_fixed = false;
        _exchangeabilities_fixed = false;
    }

    inline void QMatrix::fix_state_freqs(bool is_fixed) {
        _state_freqs_fixed = is_fixed;
    }

    inline void QMatrix::fix_exchangeabilities(bool is_fixed) {
        _exchangeabilities_fixed = is_fixed;
    }

    inline bool QMatrix::is_fixed_state_freqs() const {
        return _state_freqs_fixed;
    }

    inline bool QMatrix::is_fixed_exchangeabilities() const {
        return _exchangeabilities_fixed;
    }

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

            void                        set_equal_state_freqs(freq_xchg_ptr_t freq_ptr);
            void                        set_state_freqs_shared_ptr(freq_xchg_ptr_t freq_ptr);
            void                        set_state_freqs(const freq_xchg_t & freqs);
            freq_xchg_ptr_t             get_state_freqs_shared_ptr();
            const double *              get_state_freqs() const;

            void                        set_equal_exchangeabilities(freq_xchg_ptr_t xchg_ptr);
            void                        set_exchangeabilities_shared_ptr(freq_xchg_ptr_t xchg_ptr);
            void                        set_exchangeabilities(const freq_xchg_t & xchg);
            freq_xchg_ptr_t             get_exchangeabilities_shared_ptr();
            const double *              get_exchangeabilities() const;

        private:

            freq_xchg_ptr_t             _state_freqs;
            freq_xchg_ptr_t             _exchangeabilities;
    };

    inline QMatrixNucleotide::QMatrixNucleotide() {
        //std::cout << "Constructing a QMatrixNucleotide object" << std::endl;
        clear();
    }

    inline QMatrixNucleotide::~QMatrixNucleotide() {
        //std::cout << "Destroying a QMatrixNucleotide object" << std::endl;
    }

    inline void QMatrixNucleotide::clear() {
        QMatrix::clear();

        QMatrix::freq_xchg_t xchg = {1,1,1,1,1,1};
        _exchangeabilities = std::make_shared<QMatrix::freq_xchg_t>(xchg);

        QMatrix::freq_xchg_t freq_vect = {0.25, 0.25, 0.25, 0.25};
        _state_freqs = std::make_shared<QMatrix::freq_xchg_t>(freq_vect);

    }

    inline QMatrix::freq_xchg_ptr_t QMatrixNucleotide::get_exchangeabilities_shared_ptr() {
        return _exchangeabilities;
    }

    inline QMatrix::freq_xchg_ptr_t QMatrixNucleotide::get_state_freqs_shared_ptr() {
        return _state_freqs;
    }

    inline const double * QMatrixNucleotide::get_exchangeabilities() const {
        return _exchangeabilities->data();
    }

    inline const double * QMatrixNucleotide::get_state_freqs() const {
        return _state_freqs->data();
    }

    inline void QMatrixNucleotide::set_equal_exchangeabilities(QMatrix::freq_xchg_ptr_t xchg_ptr) {
        _exchangeabilities = xchg_ptr;
        _exchangeabilities->assign(6, 1.0/6.0);
    }

    inline void QMatrixNucleotide::set_exchangeabilities_shared_ptr(QMatrix::freq_xchg_ptr_t xchg_ptr) {
        assert(xchg_ptr->size() == 6);
        _exchangeabilities = xchg_ptr;
        normalize_freqs_or_exchangeabilities(_exchangeabilities);
    }

    inline void QMatrixNucleotide::set_exchangeabilities(const QMatrix::freq_xchg_t & xchg) {
        assert(xchg.size() == 6);
        std::copy(xchg.begin(), xchg.end(), _exchangeabilities->begin());
    }

    inline void QMatrixNucleotide::set_equal_state_freqs(QMatrix::freq_xchg_ptr_t freq_ptr) {
        _state_freqs = freq_ptr;
        _state_freqs->assign(4, 0.25);
    }

    inline void QMatrixNucleotide::set_state_freqs_shared_ptr(QMatrix::freq_xchg_ptr_t freq_ptr) {
        assert(freq_ptr->size() == 4);
        double sum_of_freqs = std::accumulate(freq_ptr->begin(), freq_ptr->end(), 0.0);
        if (std::fabs(sum_of_freqs - 1.0) > 0.001) {
            std::ostringstream msg;
            msg << "Expecting sum of 4 state frequencies to be 1, but instead got " << sum_of_freqs << std::endl;
            throw EcoevolityError(msg.str());
        }
        _state_freqs = freq_ptr;
    }

    inline void QMatrixNucleotide::set_state_freqs(const QMatrix::freq_xchg_t & freqs) {
        assert(freqs.size() == 4);
        std::copy(freqs.begin(), freqs.end(), _state_freqs->begin());
    }
}
