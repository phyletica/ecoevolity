#pragma once

#include <algorithm>
#include <vector>
#include "qmatrix.hpp"
#include "asrv.hpp"

namespace ecoevolity {

    class Model {

        public:
            typedef std::shared_ptr<Model>  SharedPtr;

                                            Model();
                                            ~Model();

            QMatrix::SharedPtr              get_q_matrix() const;
            ASRV::SharedPtr                 get_asrv() const;

            void                            init(
                                                QMatrix::SharedPtr qmatrix,
                                                ASRV::SharedPtr asrv);

            void                            clear();

        private:

            QMatrix::SharedPtr                      _qmatrix;
            ASRV::SharedPtr                         _asrv;
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
        if (this->_asrv) {
            this->_asrv->clear();
        }
    }

    inline void Model::init(
            QMatrix::SharedPtr qmatrix,
            ASRV::SharedPtr asrv) {
        this->_qmatrix = qmatrix;
        this->_asrv = asrv;
    }

    inline QMatrix::SharedPtr Model::get_q_matrix() const {
        return _qmatrix;
    }

    inline ASRV::SharedPtr Model::get_asrv() const {
        return _asrv;
    }
}
