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
        "--pool-sizes",
        poolsizes.value,
        "Pool sizes for all samples that are used (not filtered out). Either "
        "(1) a single pool size that is used for all samples, "
        "(2) a comma- or tab-separated list of pool sizes in the same order as the samples in the "
        "input file (including samples that will be filtered out), "
        "(3) a path to a file with one pool size per line, in the same order as the samples in the "
        "input file (that is, same as (2), but in a file instead of on the command line), or "
        "(4) a path to a file that contains a comma- or tab-separated list of sample names and "
        "pool sizes, with one name/size pair per line, in any order."
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
            "Option `--poolsizes` needs to be provided."
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

        // We cover both variants here, file with single values per line,
        // or file with name to size mapping. We check for consistency using the indicator variable.
        // 0 means: not yet set (before first line was processed), 1 and 2 mean that this
        // either single value or map from name to size.
        size_t type = 0;

        // Read the file line by line and process. We keep a map, but might not use it afterwards.
        auto const lines = file_read_lines( poolsizes.value );
        std::unordered_map<std::string, size_t> ps_map;
        for( size_t i = 0; i < lines.size(); ++i ) {
            auto const& line = lines[i];

            // Dissect the line and see if we got a sample name and a number, or a single value.
            auto const pair = split( line, ",\t", false );
            if( type != 0 && type != pair.size() ) {
                throw CLI::ValidationError(
                    poolsizes.option->get_name() + "(" +
                    poolsizes.value + ")",
                    "Invalid pool sizes file that contains an invalid line at " +
                    std::to_string( i + 1 ) + "."
                );
            }

            // Error check on first line of file.
            if( type == 0 && pair.size() == 1 && lines.size() != sample_names.size() ) {
                throw CLI::ValidationError(
                    poolsizes.option->get_name() + "(" +
                    poolsizes.value + ")",
                    "Invalid file with list of pool sizes that contains " +
                    std::to_string( lines.size() ) + " entries, which is different from number of " +
                    "samples in the input file " + std::to_string( sample_names.size() ) + "."
                );
            }

            // Now set the type for subsequent iterations, so that we can check that the file
            // is consistently either a single value, or a name value pair.
            type = pair.size();

            if( pair.size() == 1 ) {
                // For single value, set the entry in the result.
                assert( lines.size() == result.size() );
                assert( i < result.size() );
                result[i] = convert_poolsize_( line );
            } else if(  pair.size() == 2 ) {
                // For name value pairs, do a duplicate check first...
                if( ps_map.count( pair[0] ) > 0 ) {
                    throw CLI::ValidationError(
                        poolsizes.option->get_name() + "(" +
                        poolsizes.value + ")",
                        "Invalid line that contains duplicate sample names (line " +
                        std::to_string( i + 1 ) + ")."
                    );
                }

                // ... then add the entry to the map.
                assert( ps_map.count( pair[0] ) == 0 );
                ps_map[ pair[0] ] = convert_poolsize_( pair[1] );
            } else {
                // Everything above 2 values per line is an error.
                throw CLI::ValidationError(
                    poolsizes.option->get_name() + "(" +
                    poolsizes.value + ")",
                    "Invalid line that does not contain a pool size or pair of sample name and "
                    "pool size (line " + std::to_string( i + 1 ) + ")."
                );
            }
        }

        // If we filled the map, use it to set sizes.
        if( type == 2 ) {
            // Fill the vector of pool sizes with values in the correct order
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
        }
    } else {

        // Non-file case. If it is a single value, use it for all. If it's a list, split.
        auto const vals = split( poolsizes.value, ",\t" );
        if( vals.size() == 1 ) {
            // Just give every sample the same pool size.
            auto const ps = convert_poolsize_( poolsizes.value );
            for( auto& entry : result ) {
                entry = ps;
            }
        } else if( vals.size() == sample_names.size() ) {
            for( size_t i = 0; i < vals.size(); ++i ) {
                result[i] = convert_poolsize_( vals[i] );
            }
        } else {
            throw CLI::ValidationError(
                poolsizes.option->get_name() + "(" +
                poolsizes.value + ")",
                "Invalid list of pool sizes provided with " + std::to_string( vals.size() ) +
                " entries, which is different from number of samples in the input file " +
                std::to_string( sample_names.size() ) + "."
            );
        }
    }
    assert( result.size() == sample_names.size() );

    // User output for --verbose
    LOG_MSG2 << "Using pool sizes:";
    for( size_t i = 0; i < sample_names.size(); ++i ) {
        LOG_MSG2 << "  - " << sample_names[i] << ": " << result[i];
    }

    return result;
}
