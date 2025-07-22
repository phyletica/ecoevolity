#pragma once

#include <algorithm>
#include <vector>
#include "qmatrix.hpp"

namespace ecoevolity {

    class Model {

        public:
            typedef std::shared_ptr<Model>  SharedPtr;

                                            Model();
                                            ~Model();

            QMatrix::SharedPtr              get_q_matrix() const;
            double                          get_asrv_shape() const;
            double                          get_asrv_one_over_shape() const;
            unsigned                        get_asrv_num_categ() const;

            void                            init(
                                                QMatrix::SharedPtr qmatrix,
                                                unsigned asrv_num_categ,
                                                double asrv_one_over_shape);

            void                            clear();

        private:

            QMatrix::SharedPtr              _qmatrix;
            double                          _asrv_one_over_shape;
            unsigned                        _asrv_num_categ;
        };

    inline Model::Model() {
        //std::cout << "Constructing a Model" << std::endl;
        clear();
    }

    inline Model::~Model() {
        //std::cout << "Destroying a Model" << std::endl;
    }

    inline void Model::clear() {
        if (this->_qmatrix) {
            this->_qmatrix->clear();
        }
        this->_asrv_one_over_shape = 1.0;
    }

    inline void Model::init(
            QMatrix::SharedPtr qmatrix,
            unsigned asrv_num_categ,
            double asrv_one_over_shape) {
        this->_qmatrix = qmatrix;
        this->_asrv_num_categ = asrv_num_categ;
        this->_asrv_one_over_shape = asrv_one_over_shape;
    }

    inline QMatrix::SharedPtr Model::get_q_matrix() const {
        return _qmatrix;
    }

    inline double Model::get_asrv_shape() const {
        return 1.0/_asrv_one_over_shape;
    }

    inline double Model::get_asrv_one_over_shape() const {
        return _asrv_one_over_shape;
    }

    inline unsigned Model::get_asrv_num_categ() const {
        return _asrv_num_categ;
    }
}
