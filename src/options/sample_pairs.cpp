/*
    grenedalf - Genome Analyses of Differential Allele Frequencies
    Copyright (C) 2020-2024 Lucas Czech

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Contact:
    Lucas Czech <lucas.czech@sund.ku.dk>
    University of Copenhagen, Globe Institute, Section for GeoGenetics
    Oster Voldgade 5-7, 1350 Copenhagen K, Denmark
*/

#include "options/sample_pairs.hpp"

#include "options/global.hpp"
#include "tools/cli_setup.hpp"
#include "tools/misc.hpp"

#include "genesis/utils/core/algorithm.hpp"
#include "genesis/utils/core/fs.hpp"
#include "genesis/utils/core/logging.hpp"
#include "genesis/utils/core/std.hpp"
#include "genesis/utils/text/convert.hpp"
#include "genesis/utils/text/string.hpp"

#include <algorithm>
#include <cassert>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <tuple>
#include <vector>

// =================================================================================================
//      Setup Functions
// =================================================================================================

void SamplePairsOptions::add_sample_pairs_opts_to_app(
    CLI::App* sub,
    std::string const& group
) {
    // Settings: Comparand
    comparand.option = sub->add_option(
        "--comparand",
        comparand.value,
        "By default, statistics between all pairs of samples (that are not filtered) are computed. "
        "If this option is given a sample name however, only the pairwise statistics between that "
        "sample and all others (that are not filtered) are computed."
    );
    comparand.option->group( group );

    // Settings: Second comparand
    second_comparand.option = sub->add_option(
        "--second-comparand",
        second_comparand.value,
        "If in addition to `--comparand`, this option is also given a (second) sample name, only "
        "the statistics between those two samples are computed."
    );
    second_comparand.option->group( group );
    second_comparand.option->needs( comparand.option );

    // Settings: Comparand list
    comparand_list.option = sub->add_option(
        "--comparand-list",
        comparand_list.value,
        "By default, statistics between all pairs of samples are computed. If this option is given "
        "a file containing comma- or tab-separated pairs of sample names (one pair per line) "
        "however, only these pairwise statistics are computed."
    );
    comparand_list.option->group( group );
    comparand_list.option->check( CLI::ExistingFile );
    comparand_list.option->excludes( comparand.option );
    comparand_list.option->excludes( second_comparand.option );
}

// =================================================================================================
//      Run Functions
// =================================================================================================

bool SamplePairsOptions::is_all_to_all() const
{
    return ( ! *comparand.option ) && ( ! *comparand_list.option );
}

// We want to build a hash map for pairs, which needs hashing first...
struct pair_hash {
    template <class T1, class T2>
    std::size_t operator () (const std::pair<T1,T2> &p) const {
        auto h1 = std::hash<T1>{}(p.first);
        auto h2 = std::hash<T2>{}(p.second);
        return genesis::utils::hash_combine( h1, h2 );
    }
};

std::vector<std::pair<size_t, size_t>> SamplePairsOptions::get_sample_pairs(
    std::vector<std::string> const& sample_names
) const {
    using namespace genesis::utils;

    // We here build a vector of pairs, in order to keep the order of pairs as specified.
    // This is important, as otherwise our assignmen of which samples and pairs belong
    // to which value computed will get messed up.

    // Get sample names and map to their index or throw if invalid name is given.
    // We create a map of names to their index in the list; if used, the names need to be unique.
    // We however only need it when a comparand of sorts is provided, so that we can look up the
    // names. Without comparand, we do not need the names, and hence can allow duplicates.
    // We hence use a lambda to only fill the map when needed.
    std::unordered_map<std::string, size_t> sample_map;
    auto make_sample_map = [&](){
        for( size_t i = 0; i < sample_names.size(); ++i ) {
            auto const& name = sample_names[i];
            if( sample_map.count( name ) > 0 ) {
                throw std::runtime_error(
                    "Duplicate sample name \"" + name + "\". "
                    "Cannot uniquely identify the samples used for the FST comparands."
                );
            }
            sample_map[name] = i;
        }
    };
    auto sample_index = [&]( std::string const& name ) -> size_t {
        auto it = sample_map.find( name );
        if( it == sample_map.end() ) {
            throw CLI::ValidationError(
                "Comparand sample names",
                "Invalid sample name: \"" + name  + "\" that was not found in the input, "
                "or was filtered out."
            );
        }
        return it->second;
    };

    // Get all pairs of samples for which we want to compute FST.
    std::vector<std::pair<size_t, size_t>> sample_pairs;
    if( *comparand.option ) {
        make_sample_map();
        if( *second_comparand.option ) {
            // Only exactly one pair of samples.
            auto const index_a = sample_index( comparand.value );
            auto const index_b = sample_index( second_comparand.value );
            sample_pairs.emplace_back( index_a, index_b );
        } else {
            // One sample against all others.
            auto const index_a = sample_index( comparand.value );
            for( auto const& sn : sample_names ) {
                if( sn != comparand.value ) {
                    auto const index_b = sample_index( sn );
                    sample_pairs.emplace_back( index_a, index_b );
                }
            }
        }
    } else if( *comparand_list.option ) {
        // Read list of pairs from file, and prepare the name index lookup.
        auto const lines = file_read_lines( comparand_list.value );
        make_sample_map();

        // Get all pairs from the file, and build the list.
        for( size_t i = 0; i < lines.size(); ++i ) {
            auto const& line = lines[i];
            auto const pair = split( line, ",\t", false );
            if( pair.size() != 2 ) {
                throw CLI::ValidationError(
                    comparand_list.option->get_name() + "(" +
                    comparand_list.value + ")",
                    "Invalid line that does not contain two sample names (line " +
                    std::to_string( i + 1 ) + ")."
                );
            }
            auto const index_a = sample_index( pair[0] );
            auto const index_b = sample_index( pair[1] );
            sample_pairs.emplace_back( index_a, index_b );
        }
    } else {
        // All pairs. Build upper triangle list of sample indices.
        for( size_t i = 0; i < sample_names.size(); ++i ) {
            for( size_t j = i + 1; j < sample_names.size(); ++j ) {
                sample_pairs.emplace_back( i, j );
            }
        }
    }

    // Check for duplicates. Build a hash map for fast lookup of which pairs are there already.
    // We tried with a vector first for simplicity, but that just got too expensive
    // when many samples (in the thousands) are involved, as this means millions of pairs,
    // each involving a linear search through the vector to check for duplicates...
    std::unordered_set<std::pair<size_t, size_t>, pair_hash> dups;
    for( auto const& pair : sample_pairs ) {
        auto const rev_pair = std::pair<size_t, size_t>( pair.second, pair.first );
        if( dups.count( pair ) > 0 || dups.count( rev_pair ) > 0 ) {
            LOG_WARN << "Sample pairing for sample \"" << sample_names[pair.first]
                     << "\" and sample \"" << sample_names[pair.second]
                     << "\" is provided more than once.";
        }
        dups.insert( pair );
    }

    return sample_pairs;
}
