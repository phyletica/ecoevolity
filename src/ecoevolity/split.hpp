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

#ifndef ECOEVOLITY_SPLIT_HPP
#define ECOEVOLITY_SPLIT_HPP

#include <vector>
#include <memory>
#include <set>
#include <map>
#include <climits>
#include <cassert>


/**
 * Object for storing set of taxa that descend from a node.
 *
 * @note    Modified from:
 *              Split class of the Strom 3 BEAGLE tutorial
 *              <https://stromtutorial.github.io/>
 *              License:    MIT
 *              Author:     Paul Lewis
 *              Citation:   D Ayers, MP Cummings, G Baele, AE Darling, PO
 *                          Lewis, DL Swofford, JP Huelsenbeck, P Lemey, A
 *                          Rambaut, and MA Suchard. 2019. BEAGLE 3: Improved
 *                          performance, scaling, and usability for a
 *                          high-performance computing library for statistical
 *                          phylogenetics. Systematic Biology 61:170-173. DOI:
 *                          10.1093/sysbio/syz020
 */
class Split {
    public:
        Split();
        Split(const Split & other);
        ~Split();

        Split & operator=(const Split & other);
        bool operator==(const Split & other) const;
        bool operator!=(const Split & other) const;
        bool operator<(const Split & other) const;

        void clear();
        unsigned int size() const;
        void resize(unsigned int number_of_leaves);

        typedef unsigned long                                   split_unit_t;
        typedef std::vector<split_unit_t>                       split_t;
        typedef std::set<Split>                                 treeid_t;
        typedef std::map< treeid_t, std::vector<unsigned int> > treemap_t;
        typedef std::tuple<unsigned int, unsigned int, unsigned int> split_metrics_t;

        typedef std::shared_ptr<Split> SharedPtr;

        split_unit_t get_bits(const unsigned int unit_index) const;
        bool get_leaf_bit(const unsigned int leaf_index) const;
        void set_leaf_bit(const unsigned int leaf_index);
        void add_split(const Split & other);

        bool is_equivalent(const Split & other,
                const bool strict_root = true) const;
        bool is_compatible(const Split & other) const;
        bool conflicts_with(const Split & other) const;

        std::string as_string(const char unset_char = '0',
                const char set_char = '1') const;
        split_metrics_t get_split_metrics() const;

    private:
        split_unit_t  mask_;
        split_t       bits_;
        unsigned int  bits_per_unit_;
        unsigned int  number_of_leaves_;
};

inline Split::Split() {
    this->mask_ = 0L;
    this->number_of_leaves_ = 0;
    this->bits_per_unit_ = (CHAR_BIT)*sizeof(Split::split_unit_t);
    this->clear();
}

inline Split::Split(const Split & other) {
    this->mask_ = other.mask_;
    this->number_of_leaves_ = other.number_of_leaves_;
    this->bits_per_unit_ = (CHAR_BIT)*sizeof(Split::split_unit_t);
    this->bits_ = other.bits_;
}

inline Split::~Split() {}

inline void Split::clear() {
    for (auto & split_u : this->bits_) {
        split_u = 0L;
    }
}

inline Split & Split::operator=(const Split & other) {
    this->mask_ = other.mask_;
    this->number_of_leaves_ = other.number_of_leaves_;
    this->bits_per_unit_ = (CHAR_BIT)*sizeof(Split::split_unit_t);
    this->bits_ = other.bits_;
    return *this;
}

inline bool Split::operator==(const Split & other) const {
    return (this->bits_ == other.bits_);
}

inline bool Split::operator!=(const Split & other) const {
    return (! (*this == other));
}

inline bool Split::operator<(const Split & other) const {
    assert(this->bits_.size() == other.bits_.size());
    return (this->bits_ < other.bits_);
}

inline unsigned int Split::size() const {
    return this->number_of_leaves_;
}

inline void Split::resize(const unsigned int nleaves) {
    this->number_of_leaves_ = nleaves;
    unsigned int nunits = 1 + ((nleaves - 1) / this->bits_per_unit_);
    this->bits_.resize(nunits);

    // create mask used to select only those bits used in final unit
    unsigned int num_unused_bits = nunits * this->bits_per_unit_ - nleaves;
    unsigned int num_used_bits = this->bits_per_unit_ - num_unused_bits;
    this->mask_ = 0L;
    split_unit_t unity = 1;
    for (unsigned int i = 0; i < num_used_bits; ++i) {
        this->mask_ |= (unity << i);
    }
    this->clear();
}

inline void Split::set_leaf_bit(const unsigned int leaf_index) {
    unsigned int unit_index = leaf_index / this->bits_per_unit_;
    unsigned int bit_index = leaf_index - unit_index * this->bits_per_unit_;
    split_unit_t unity = 1;
    split_unit_t bit_to_set = unity << bit_index;
    this->bits_.at(unit_index) |= bit_to_set;
}

inline Split::split_unit_t Split::get_bits(const unsigned int unit_index) const {
    // assert(unit_index < this->bits_.size());
    return this->bits_.at(unit_index);
}

inline bool Split::get_leaf_bit(const unsigned leaf_index) const {
    unsigned int unit_index = leaf_index / this->bits_per_unit_;
    unsigned int bit_index = leaf_index - unit_index * this->bits_per_unit_;
    split_unit_t unity = 1;
    split_unit_t bit_to_check = unity << bit_index;
    return (bool)(this->bits_.at(unit_index) & bit_to_check);
}

inline void Split::add_split(const Split & other) {
    unsigned nunits = (unsigned int)this->bits_.size();
    assert(nunits == other.bits_.size());
    for (unsigned int i = 0; i < nunits; ++i) {
        this->bits_.at(i) |= other.bits_.at(i);
    }
}

inline std::string Split::as_string(
        const char unset_char,
        const char set_char) const {
    std::ostringstream ss;
    unsigned int nleaves_added = 0;
    for (unsigned int i = 0; i < this->bits_.size(); ++i) {
        for (unsigned int j = 0; j < this->bits_per_unit_; ++j) {
            split_unit_t bitmask = ((split_unit_t)1 << j);
            bool bit_is_set = ((this->bits_.at(i) & bitmask) > (split_unit_t)0);
            if (bit_is_set) {
                ss << set_char;
            }
            else {
                ss << unset_char;
            }
            if (++nleaves_added == this->number_of_leaves_) {
                break;
            }
        }
    }
    return ss.str();
}

inline bool Split::is_equivalent(const Split & other,
        const bool strict_root) const {
    if (strict_root) {
        return (this->bits_ == other.bits_);
    }
    unsigned int nunits = (unsigned int)this->bits_.size();
    assert(nunits > 0);

    // polarity 1 means root is on the same side of both splits
    // polarity 2 means they are inverted relative to one another
    unsigned int polarity = 0;
    for (unsigned int i = 0; i < nunits; ++i) {
        split_unit_t a = this->bits_.at(i);
        split_unit_t b = other.bits_.at(i);
        bool a_equals_b = (a == b);
        bool a_equals_inverse_b = (a == ~b);
        if (i == nunits - 1) {
            a_equals_inverse_b = (a == (~b & this->mask_));
        }
        bool ok = (a_equals_b || a_equals_inverse_b);
        if (ok) {
            if (polarity == 0) {
                // First unit determines polarity
                if (a_equals_b) {
                    polarity = 1;
                }
                else {
                    polarity = 2;
                }
            }
            else {
                // Polarity determined by first unit used for all subsequent units
                if (polarity == 1 && ! a_equals_b) {
                    return false;
                }
                else if (polarity == 2 && ! a_equals_inverse_b) {
                    return false;
                }
            }
        }
        else {
            return false;
        }
    }
    return true;
}

inline bool Split::is_compatible(const Split & other) const {
    for (unsigned int i = 0; i < this->bits_.size(); ++i) {
        split_unit_t a = this->bits_.at(i);
        split_unit_t b = other.bits_.at(i);
        split_unit_t a_and_b = (a & b);
        bool equals_a = (a == a_and_b);
        bool equals_b = (b == a_and_b);
        if (a_and_b && ! (equals_a || equals_b)) {
            // A failure of any unit to be compatible makes the entire split incompatible
            return false;
        }
    }
    return true;
}

inline bool Split::conflicts_with(const Split & other) const {
    return (! this->is_compatible(other));
}

#endif
