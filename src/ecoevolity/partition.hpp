#pragma once

#include <tuple>
#include <numeric>
#include <cmath>
#include "error.hpp"

namespace ecoevolity {

    class Partition {
        public:
            typedef std::match_results<std::string::const_iterator>::const_reference    regex_match_t;
            typedef std::tuple<unsigned, unsigned, unsigned>                            site_range_t;
            typedef std::vector<site_range_t>                                           subset_t;
            typedef std::vector<subset_t>                                               partition_t;
            typedef std::vector<unsigned>                                               subset_sizes_vect_t;
            typedef std::vector<std::string>                                            subset_names_vect_t;
            typedef std::shared_ptr<Partition>                                          SharedPtr;

                                                        Partition();
                                                        ~Partition();
                                                        Partition(
                                                                unsigned num_sites,
                                                                const std::vector<NucSubsetSettings> & partition_settings);
        
            unsigned                                    get_num_sites() const;
            unsigned                                    get_num_subsets() const;
            std::string                                 get_subset_name(unsigned subset) const;
        
            unsigned                                    find_subset_by_name(const std::string & subset_name) const;
            unsigned                                    find_subset_for_site(unsigned site_index) const;
            bool                                        site_in_subset(unsigned site_index, unsigned subset_index) const;
        
            unsigned                                    num_sites_in_subset(unsigned subset_index) const;

            std::pair< std::vector<unsigned>,
                       std::vector<unsigned> >          get_duplicated_and_missing_sites() const;

            void                                        clear();

        private:

            site_range_t                                parse_site_range(std::string range_definition) const;

            int                                         extract_int_from_regex_match(regex_match_t s, unsigned min_value) const;
    
            unsigned                                    _num_sites;
            subset_names_vect_t                         _subset_names;
            subset_sizes_vect_t                         _subset_sizes;
            partition_t                                 _subsets;
    };
    
    inline Partition::Partition() {
        //std::cout << "Constructing a Partition" << std::endl;
        clear();
    }  

    inline Partition::Partition(
            unsigned num_sites,
            const std::vector<NucSubsetSettings> & partition_settings) : Partition::Partition() {
        this->_num_sites = num_sites;
        for (const auto & subset_info : partition_settings) {
            this->_subset_names.push_back(subset_info.name);
            subset_t subset;
            for (const auto & sites_str : subset_info.sites) {
                subset.push_back(this->parse_site_range(sites_str));
            }
            this->_subsets.push_back(subset);
        }
        this->update_num_sites_in_subsets();

        std::pair< std::vector<unsigned>, std::vector<unsigned> > dups_missing;
        dups_missing = this->get_duplicated_and_missing_sites();
        if (dups_missing->first.size() > 0) {
            std::ostringstream msg;
            msg << dups_missing->first.size()
                << " sites were defined in multiple subsets";
            throw EcoevolityError(msg.str());
        }
        if (dups_missing->second.size() > 0) {
            std::cerr << "WARNING: "
                << dups_missing->second.size()
                << " sites were missing from partition definition." << std::endl
                << "These sites will be ignored during the analysis.";
        }
    }

    inline Partition::~Partition() {
        //std::cout << "Destroying a Partition" << std::endl;
    }

    inline unsigned Partition::get_num_sites() const {
        return this->_num_sites;
    }
    
    inline unsigned Partition::get_num_subsets() const {
        return this->_subsets.size();
    }
    
    inline std::string Partition::get_subset_name(unsigned subset) const {
        assert(subset < this->_num_subsets);
        return this->_subset_names.at(subset);
    }
    
    inline unsigned Partition::find_subset_by_name(const std::string & subset_name) const {
        auto iter = std::find(this->_subset_names.begin(), this->_subset_names.end(), subset_name);
        if (iter == this->_subset_names.end()) {
            std::ostringstream msg;
            msg << "Specified subset name \'"
                << subset_name
                << "\' not found in partition";
            throw EcoevolityError(msg.str());
        }
        return (unsigned)std::distance(this->_subset_names.begin(), iter);
    }

    inline unsigned Partition::find_subset_for_site(unsigned site_index) const {
        if (site_index >= this->_num_sites) {
            std::ostringstream msg;
            msg << "Site number \'"
                << site_index + 1
                << "\' is beyond the number of sites";
            throw EcoevolityError(msg.str());
        }
        for (unsigned subset_idx = 0; subset_idx < this->_subsets.size(); ++subset_idx) {
            for (auto & t : this->_subsets.at(subset_idx)) {
                unsigned begin_idx = std::get<0>(t);
                unsigned past_idx = std::get<1>(t);
                unsigned stride = std::get<2>(t);
                bool inside_range = site_index >= begin_idx && site_index < past_idx;
                if (inside_range && (site_index - begin_idx) % stride == 0)
                    return subset_idx;
            }
        }
        std::ostringstream msg;
        msg << "Site "
            << site_index + 1
            << " not found in any subset of partition";
        throw EcoevolityError(msg.str());
    }
    
    inline bool Partition::site_in_subset(unsigned site_index, unsigned subset_index) const {
        unsigned which_subset = this->find_subset_for_site(site_index);
        return (which_subset == subset_index ? true : false);
    }

    inline unsigned Partition::num_sites_in_subset(unsigned subset_index) const {
        return this->_subset_sizes.at(subset_index);
    }
    
    inline unsigned Partition::num_sites_in_subset(unsigned subset_index) const {
        unsigned nsites = 0;
        for (auto & t : this->_subsets.at(subset_index)) {
            unsigned begin_idx = std::get<0>(t);
            unsigned past_idx = std::get<1>(t);
            unsigned stride = std::get<2>(t);
            unsigned n = past_idx - begin_site;
            nsites += (unsigned)(floor(n/stride)) + (n % stride == 0 ? 0 : 1);
        }
        return nsites;
    }
    
    inline void Partition::clear() {
        this->_num_sites = 0;
        this->_subset_names.clear();
        this->_subsets.clear();
    }

    inline void Partition::update_num_sites_in_subsets() {
        this->_subset_sizes.clear();
        this->_subset_sizes_.reserve(this->_subsets.size());
        for (unsigned subset_index = 0; subset_index < this->_subsets.size(); ++subset_index) {
            unsigned nsites = 0;
            for (auto & t : this->_subsets.at(subset_index)) {
                unsigned begin_idx = std::get<0>(t);
                unsigned past_idx = std::get<1>(t);
                unsigned stride = std::get<2>(t);
                unsigned n = past_idx - begin_site;
                nsites += (unsigned)(floor(n/stride)) + (n % stride == 0 ? 0 : 1);
            }
            this->_subset_sizes.push_back(nsites);
        }
    }
    
    inline site_range_t Partition::parse_site_range(std::string range_definition) const {
        // match patterns like these: "1-.\3" "1-1000" "1001-."
        const char * pattern_string = R"((\d+)\s*(-\s*([0-9.]+)(\\\s*(\d+))*)*)";
        std::regex re(pattern_string);
        std::smatch match_obj;
        bool matched = std::regex_match(range_definition, match_obj, re);
        if (!matched) {
            std::ostringstream msg;
            msg << "Could not interpret \'"
                << range_definition
                << "\' as a range of sites";
            throw EcoevolityError(msg.str());
        }
        
        // match_obj always yields 6 strings that can be indexed using the operator[] function
        // match_obj[0] equals entire site_range (e.g. "1-.\3")
        // match_obj[1] equals beginning site index (e.g. "1")
        // match_obj[2] equals everything after beginning site index (e.g. "-.\3")
        // match_obj[3] equals "" or ending site index (e.g. ".")
        // match_obj[4] equals "" or everything after ending site index (e.g. "\3")
        // match_obj[5] equals "" or step value (e.g. "3")
        int ibegin = this->extract_int_from_regex_match(match_obj[1], 1);
        int iend   = this->extract_int_from_regex_match(match_obj[3], ibegin);
        int istep  = this->extract_int_from_regex_match(match_obj[5], 1);

        assert(ibegin > 0);
        if (ibegin < 1) {
            std::ostringstream msg;
            msg << "Site number \'"
                << ibegin
                << "\' in subset definition is not valid";
            throw EcoevolityError(msg.str());
        }
        if (iend < ibegin) {
            std::ostringstream msg;
            msg << "Site range definition \'"
                << range_definition
                << "\' is not valid";
            throw EcoevolityError(msg.str());
        }
        
        return std::make_tuple(ibegin - 1, iend, istep);
    }
    
    inline int Partition::extract_int_from_regex_match(regex_match_t s, unsigned min_value) const {
        assert(min_value > 0);
        int int_value = min_value;
        if (s.length() > 0) {
            std::string str_value = s.str();
            if (str_value == ".") {
                return this->_num_sites;
            }
            try {
                int_value = std::stoi(str_value);
            }
            catch(const std::invalid_argument & e) {
                std::ostringstream msg;
                msg << "Could not interpret \'"
                    << s.str()
                    << "\' as a number in partition subset definition";
                throw EcoevolityError(msg.str());
            }
            // sanity check
            if (int_value < (int)min_value) {
                std::ostringstream msg;
                msg << "Value specified in partition subset definition ("
                    << int_value
                    << ") is lower than minimum value ("
                    << min_value
                    << ")";
                throw EcoevolityError(msg.str());
            }
        }
        return int_value;
    }

    std::pair< std::vector<unsigned>, std::vector<unsigned> > Partition::get_duplicated_and_missing_sites() const {
        std::vector<unsigned> dups;
        std::vector<unsigned> missing;
        std::vector<int> sites(this->_num_sites, 1);
        for (unsigned subset_idx = 0; subset_idx < this->_subsets.size(); ++subset_idx) {
            for (auto & t : this->_subsets.at(subset_idx)) {
                unsigned begin_idx = std::get<0>(t);
                unsigned past_idx = std::get<1>(t);
                unsigned stride = std::get<2>(t);
                for (unsigned site_idx = begin_idx;
                        site_idx < past_idx;
                        site_idx += stride) {
                    if (sites.at(site_idx) == 0) {
                        dups.push_back(site_idx);
                    }
                    --sites.at(site_idx);
                }
            }
        }
        for (unsigned site_idx = 0; site_idx < sites.size(); ++site_idx) {
            if (sites.at(site_idx) == 1) {
                missing.push_back(site_idx);
            }
        }
        return std::make_pair(dups, missing);
    }
}
