/*
    grenedalf - Genome Analyses of Differential Allele Frequencies
    Copyright (C) 2020-2022 Lucas Czech

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

#include "options/variant_input_sample_names.hpp"

#include "options/global.hpp"
#include "tools/misc.hpp"

#include "genesis/utils/core/fs.hpp"
#include "genesis/utils/text/convert.hpp"
#include "genesis/utils/text/string.hpp"

#include <algorithm>
#include <cassert>
#include <stdexcept>

// =================================================================================================
//      Setup Functions
// =================================================================================================

void VariantInputSampleNamesOptions::add_sample_name_opts_to_app(
    CLI::App* sub,
    std::string const& group
) {
    // Correct setup check.
    internal_check(
        sample_name_prefix_.option == nullptr && sample_name_list_.option == nullptr,
        "Cannot use the same VariantInputSampleNamesOptions object multiple times."
    );

    // Name list option.
    sample_name_list_.option = sub->add_option(
        "--sample-name-list",
        sample_name_list_.value,
        "Some file types do not contain sample names, such as (m)pileup or sync files. For such "
        "file types, sample names can here be provided as either (1) a comma- or tab-separated "
        "list, or (2) as a file with one sample name per line, in the same order as samples are in "
        "the actual input file. We then use these names in the output and the "
        "`--filter-samples-include` and `--filter-samples-exclude` options. "
        "If not provided, we simply use numbers 1..n as sample names for these files types. "
        "Note that this option can only be used if a single file is given as input. "
        "Alternatively, use `--sample-name-prefix` to provide a prefix for this sample numbering."
    );
    sample_name_list_.option->group( group );
    // sample_name_list_.option->check( CLI::ExistingFile );

    // Prefix option.
    sample_name_prefix_.option = sub->add_option(
        "--sample-name-prefix",
        sample_name_prefix_.value,
        "Some file types do not contain sample names, such as (m)pileup or sync files. For such "
        "file types, this prefix followed by indices 1..n can be used instead to provide unique "
        "names per sample that we use in the output and the `--filter-samples-include` and "
        "`--filter-samples-exclude` options. For example, use \"Sample_\" as a prefix. "
        "If not provided, we simply use numbers 1..n as sample names for these files types. "
        "This prefix also works if multiple files are given as input. "
        "Alternatively, use `--sample-name-list` to directly provide a list of sample names."
    );
    sample_name_prefix_.option->group( group );
    // sample_name_prefix_.option->needs( pileup_file_.option );

    // The two ways of specifying sample names are mutually exclusive.
    sample_name_list_.option->excludes( sample_name_prefix_.option );
    sample_name_prefix_.option->excludes( sample_name_list_.option );

    // Add option for sample name filter.
    filter_samples_include_.option = sub->add_option(
        "--filter-samples-include",
        filter_samples_include_.value,
        "Sample names to include (all other samples are excluded); either (1) a comma- or "
        "tab-separated list, or (2) a file with one sample name per line. If no sample filter "
        "is provided, all samples in the input file are used. The option considers "
        "`--sample-name-list` or `--sample-name-prefix` for file types that do not contain sample "
        "names. Note that this option can only be used if a single file is given as input."
    );
    filter_samples_include_.option->group( group );

    // And the other way round.
    filter_samples_exclude_.option = sub->add_option(
        "--filter-samples-exclude",
        filter_samples_exclude_.value,
        "Sample names to exclude (all other samples are included); either (1) a comma- or "
        "tab-separated list, or (2) a file with one sample name per line. If no sample filter "
        "is provided, all samples in the input file are used. The option considers "
        "`--sample-name-list` or `--sample-name-prefix` for file types that do not contain sample "
        "names. Note that this option can only be used if a single file is given as input."
    );
    filter_samples_exclude_.option->group( group );

    // Include and exclude are mutually exclusive.
    filter_samples_exclude_.option->excludes( filter_samples_include_.option );
    filter_samples_include_.option->excludes( filter_samples_exclude_.option );
}

// =================================================================================================
//      Run Functions
// =================================================================================================

// -------------------------------------------------------------------------
//     process_sample_name_list_option_
// -------------------------------------------------------------------------

std::vector<std::string> VariantInputSampleNamesOptions::process_sample_name_list_option(
    std::string const& list
) const {
    using namespace genesis::utils;

    // If the input is a file, read it line by line as sample names.
    // Otherwise, split by comma or tab.
    std::vector<std::string> result;
    if( is_file( list ) ) {
        result = file_read_lines( list );
    } else {
        result = split( list, ",\t" );
    }
    if( result.empty() ) {
        throw std::runtime_error( "Invalid empty list of sample names given." );
    }

    return result;
}

// -------------------------------------------------------------------------
//     make_anonymous_sample_names_
// -------------------------------------------------------------------------

std::vector<std::string> VariantInputSampleNamesOptions::make_anonymous_sample_names(
    size_t sample_count
) const {
    // Edge case, just to make sure.
    if( sample_count == 0 ) {
        throw std::runtime_error( "Input file does not contain any samples." );
    }

    // Prepare
    std::vector<std::string> result;

    // In case we have a proper list of sample names, use that.
    if( sample_name_list_.option && *sample_name_list_.option ) {
        result = process_sample_name_list_option( sample_name_list_.value );
        if( result.size() != sample_count ) {
            throw CLI::ValidationError(
                sample_name_list_.option->get_name() + "(" + sample_name_list_.value + ")",
                "Invalid sample names list that contains " + std::to_string( result.size() ) +
                " name entries. This is incongruent with the input file, which contains " +
                std::to_string( sample_count ) +
                " samples (after filtering, if a sample name filter was given)."
            );
        }
        assert( result.size() > 0 );
        return result;
    }

    // In case we have a prefix for sample names instead, or nothing, just enumerate.
    std::string prefix;
    if( sample_name_prefix_.option && *sample_name_prefix_.option ) {
        prefix = sample_name_prefix_.value;
    }
    result.reserve( sample_count );
    for( size_t i = 0; i < sample_count; ++i ) {
        result.emplace_back( prefix + std::to_string( i + 1 ));
    }
    assert( result.size() > 0 );
    return result;
}

// -------------------------------------------------------------------------
//     find_sample_indices_from_sample_filters_
// -------------------------------------------------------------------------

std::pair<std::vector<size_t>, bool>
VariantInputSampleNamesOptions::find_sample_indices_from_sample_filters() const
{
    // Prepare result.
    std::vector<size_t> indices;

    // Get whether we want to include or exclude sample names.
    bool const is_include = ! filter_samples_include_.value.empty();
    bool const is_exclude = ! filter_samples_exclude_.value.empty();

    // Not both can be given at the same time, as we made the options mutually exclusive.
    internal_check( !( is_include && is_exclude ), "include and exclude filters are both given" );

    // If no filters are given, just return empty, indicating to the caller that we do not filter.
    // This is also the return of this function when multiple files are given, in which case
    // we cannot apply sample filtering... that would just be too tedious to provide via a CLI.
    if( ! is_include && ! is_exclude ) {
        assert( indices.empty() );
        return { indices, false };
    }

    // Get the sample names, depending on which type (inc/exc) we have.
    auto const filter_list = process_sample_name_list_option(
        is_include ? filter_samples_include_.value : filter_samples_exclude_.value
    );

    // When this function is called, we do not have opened the input file yet, so we do not
    // know how many samples there are. So here, we either need to get the list of sample names
    // (if the user provided one), or use enumeration to find the indices. The list loading is
    // duplicate work unfortunately, as it will be loaded again later. Alternatively, we could
    // split this into the case where a list of names is given, and the case where enumeration is
    // used, but that would make set_anonymous_sample_names() and its setup more complex...
    // So let's live with loading the list twice.

    if( sample_name_list_.option && *sample_name_list_.option ) {
        // In case we have a proper list of sample names, use that, and check the filters against it.
        auto const sample_names = process_sample_name_list_option( sample_name_list_.value );
        assert( sample_names.size() > 0 );

        // Go through all provided sample names for the filter, find their positon in the sample
        // name list, and fill the indices with these positions.
        for( auto const& filtered_name : filter_list ) {
            auto const it = std::find( sample_names.begin(), sample_names.end(), filtered_name );
            if( it == sample_names.end() ) {
                // Filter name not found, throw. First, helpful user output.
                LOG_MSG << "Sample names:";
                for( auto const& sn : sample_names ) {
                    LOG_MSG << " - " << sn;
                }
                throw CLI::ValidationError(
                    "Invalid sample name used for filtering: \"" + filtered_name  + "\"."
                );
            } else {
                // Add filter index to result.
                auto const index = it - sample_names.begin();
                indices.push_back( index );
            }
        }
    } else {
        // If we do not have a given list of sample names, we instead use the prefix and enumerate.
        std::string prefix;
        if( sample_name_prefix_.option && *sample_name_prefix_.option ) {
            prefix = sample_name_prefix_.value;
        }

        // Go through all filter names, remove the prefix from the filter names (might be empty),
        // and check that the remainder is a number. If so, it's the index (offset by 1).
        for( auto filtered_name : filter_list ) {
            if( ! genesis::utils::starts_with( filtered_name, prefix )) {
                throw CLI::ValidationError(
                    "Invalid sample name used for filtering: \"" + filtered_name  + "\" "
                    "that does not fit the sample prefix given by " +
                    sample_name_prefix_.option->get_name()
                );
            }

            // Here, we know that the filter name contains the prefix, so we can remove it,
            // and convert the rest to a number. We use 1 based indexing in the user provided names,
            // so zero is invalid. But then we need to offset to get our actual zero-based index.
            filtered_name = filtered_name.substr( prefix.length() );
            size_t index = 0;
            try {
                index = genesis::utils::convert_from_string<size_t>( filtered_name, true );
            } catch(...) {
                index = 0;
            }
            if( index == 0 ) {
                throw CLI::ValidationError(
                    "Invalid sample name used for filtering: \"" + filtered_name  + "\" "
                    "that does not contain a valid number after the prefix."
                );
            }
            assert( index > 0 );
            indices.push_back( index - 1 );
        }
    }

    // Return the indices, and whether they are to be interpreted as inverses.
    return { std::move( indices ), is_exclude };
}
