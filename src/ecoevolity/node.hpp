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

#ifndef ECOEVOLITY_NODE_HPP
#define ECOEVOLITY_NODE_HPP

#include "basenode.hpp"

/**
 * Base class for a node of a phylogenetic tree.
 *
 * @note    Many of the class' methods modified from:
 *              BasicTNode class of Bio++ Library
 *              <http://biopp.univ-montp2.fr/wiki/index.php/Main_Page>
 *              License:    CeCILL <http://www.cecill.info>
 *              Author:     Sylvain Gaillard
 *              Copyright:  CNRS, (January 12, 2011)
 */
class Node: public BaseNode<Node>{
    private:
        typedef BaseNode<Node> BaseClass;

    public:
        Node() { }
        Node(const Node& node) : BaseClass(node) { }
        Node(std::string label) : BaseClass(label) { }
        Node(double height) : BaseClass(height) { }
        Node(std::string label, double height) : BaseClass(label, height) { }
};

class PopulationNode: public BaseNode<PopulationNode>{
    private:
        typedef BaseNode<PopulationNode> BaseClass;
        BiallelicPatternProbabilityMatrix pattern_probs_bottom_;
        BiallelicPatternProbabilityMatrix pattern_probs_top_;

    public:
        PopulationNode() { }
        PopulationNode(const Node& node) : BaseClass(node) { }
        PopulationNode(std::string label) : BaseClass(label) { }
        PopulationNode(double height) : BaseClass(height) { }
        PopulationNode(std::string label, double height) :
            BaseClass(label, height)
            { }
        PopulationNode(unsigned int allele_count) : BaseClass() {
            this->pattern_probs_bottom_.resize(allele_count);
            this->pattern_probs_top_.resize(allele_count);
        }
        PopulationNode(std::string label, unsigned int allele_count) :
            BaseClass(label)
        {
            this->pattern_probs_bottom_.resize(allele_count);
            this->pattern_probs_top_.resize(allele_count);
        }
        PopulationNode(double height, unsigned int allele_count) :
            BaseClass(height)
        {
            this->pattern_probs_bottom_.resize(allele_count);
            this->pattern_probs_top_.resize(allele_count);
        }
        PopulationNode(std::string label,
                double height,
                unsigned int allele_count) :
            BaseClass(label, height)
        {
            this->pattern_probs_bottom_.resize(allele_count);
            this->pattern_probs_top_.resize(allele_count);
        }
// 
//         // from NodeData
//         unsigned int get_allele_count() const;
//         void resize(unsigned int allele_count);
//         void reset(unsigned int allele_count);
};

#endif
