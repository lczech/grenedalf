/*
    grenedalf - Genome Analyses of Differential Allele Frequencies
    Copyright (C) 2020-2021 Lucas Czech

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
    Lucas Czech <lczech@carnegiescience.edu>
    Department of Plant Biology, Carnegie Institution For Science
    260 Panama Street, Stanford, CA 94305, USA
*/

#include "options/poolsizes.hpp"

#include "options/global.hpp"
#include "tools/misc.hpp"

#include "genesis/utils/core/fs.hpp"
#include "genesis/utils/text/convert.hpp"
#include "genesis/utils/text/string.hpp"

#include <cassert>
#include <stdexcept>

// =================================================================================================
//      Setup Functions
// =================================================================================================

void PoolsizesOptions::add_poolsizes_opt_to_app(
    CLI::App* sub,
    bool required,
    std::string const& group
) {
    poolsizes.option = sub->add_option(
        "--poolsizes",
        poolsizes.value,
        "Poolsizes for all samples that are use (not filtered out). Either a single number that "
        "is used for all samples, or a path to a file that contains a tab-separated list of sample "
        "names and pool sizes, with one sample/size entry per line."
    );
    poolsizes.option->group( group );
    if( required ) {
        poolsizes.option->required();
    }
}

// =================================================================================================
//      Run Functions
// =================================================================================================

std::vector<size_t> PoolsizesOptions::get_pool_sizes(
    std::vector<std::string> const& sample_names,
    std::vector<bool> const& sample_filter
) const {
    using namespace genesis::utils;

    // Internal error check.
    if( ! sample_filter.empty() && sample_filter.size() != sample_names.size() ) {
        throw std::runtime_error(
            "Internal error: different sizes of sample_names and sample_filter"
        );
    }

    // Nice error message.
    if( ! *poolsizes.option ) {
        throw CLI::ValidationError(
            poolsizes.option->get_name(),
            "Pool sizes required to be provided."
        );
    }

    // Convert a pool size to a number, or throw.
    auto convert_poolsize_ = [&]( std::string const& str ){
        try {
            return convert_from_string<size_t>( str );
        } catch(...) {
            throw CLI::ValidationError(
                poolsizes.option->get_name(),
                "Invalid pool size value: " + str
            );
        }
    };

    // Get the pool sizes, depending on input type either from file, or as a single number.
    auto result = std::vector<size_t>( sample_names.size(), 0 );
    if( is_file( poolsizes.value )) {

        // Read the file line by line and make a map from sample names to pool sizes.
        auto const lines = file_read_lines( poolsizes.value );
        std::unordered_map<std::string, size_t> ps_map;
        for( size_t i = 0; i < lines.size(); ++i ) {
            auto const& line = lines[i];

            // Dissect the line and see if we got a sample name and a number.
            auto const pair = split( line, "\t", false );
            if( pair.size() != 2 ) {
                throw CLI::ValidationError(
                    poolsizes.option->get_name() + "(" +
                    poolsizes.value + ")",
                    "Invalid line that does not contain sample name and a pool size (line " +
                    std::to_string( i + 1 ) + ")."
                );
            }

            // Duplicate check.
            if( ps_map.count( pair[0] )) {
                throw CLI::ValidationError(
                    poolsizes.option->get_name() + "(" +
                    poolsizes.value + ")",
                    "Invalid line that does not contain sample name and a pool size (line " +
                    std::to_string( i + 1 ) + ")."
                );
            }

            // Add the entry.
            assert( ps_map.count( pair[0] ) == 0 );
            ps_map[ pair[0] ] = convert_poolsize_( pair[1] );
        }

        // Now, fill the vector of pool sizes with values in the correct order
        // (that is, using the sample name order).
        // Throw if the sample needs to be given (it's in the used list), but isn't in the file.
        for( size_t i = 0; i < sample_names.size(); ++i ) {
            if( ! sample_filter.empty() && ! sample_filter[i] ) {
                continue;
            }

            auto const& sn = sample_names[i];
            if( ps_map.count(sn) > 0 ) {
                result[i] = ps_map[sn];
            } else {
                throw CLI::ValidationError(
                    poolsizes.option->get_name() + "(" +
                    poolsizes.value + ")",
                    "Sample name \"" + sn +  "\" missing from pool size file."
                );
            }
        }
    } else {

        // Just give every sample the same pool size.
        auto const ps = convert_poolsize_( poolsizes.value );
        for( auto& entry : result ) {
            entry = ps;
        }
    }
    return result;
}
