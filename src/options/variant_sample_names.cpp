/*
    grenedalf - Genome Analyses of Differential Allele Frequencies
    Copyright (C) 2020-2023 Lucas Czech

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

#include "options/variant_sample_names.hpp"

#include "options/global.hpp"
#include "tools/misc.hpp"

#include "genesis/population/functions/variant_input_iterator.hpp"
#include "genesis/utils/core/fs.hpp"
#include "genesis/utils/text/convert.hpp"
#include "genesis/utils/text/string.hpp"

#include <algorithm>
#include <cassert>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

// =================================================================================================
//      Setup Functions
// =================================================================================================

void VariantSampleNamesOptions::add_sample_name_opts_to_app(
    CLI::App* sub,
    std::string const& group
) {
    // Correct setup check.
    internal_check(
        rename_samples_.option == nullptr,
        "Cannot use the same VariantSampleNamesOptions object multiple times."
    );

    // Rename samples option.
    rename_samples_.option = sub->add_option(
        "--rename-samples-file",
        rename_samples_.value,
        "Allows to rename samples, by providing a file that lists the old and new sample names, "
        "one per line, separating old and new names by a tab.\n"
        "By default, we use sample names as provided in the input files. Some file types however do "
        "not contain sample names, such as (m)pileup or sync files. For such file types, sample "
        "names are automatically assigned by using their input file base name (without path and "
        "extension), followed by a dot and numbers 1..n for all samples in that file. "
        "For instance, samples in `/path/to/sample.sync` are named `sample.1`, `sample.2`, etc.\n"
        "Using this option, those names can be renamed as needed. Use verbose output (`--verbose`) "
        "to show a list of all sample names. We then use these names in the output and the "
        "`--filter-samples-include` and `--filter-samples-exclude` options. "
    );
    rename_samples_.option->group( group );
    rename_samples_.option->check( CLI::ExistingFile );

    // Add option for sample name filter.
    filter_samples_include_.option = sub->add_option(
        "--filter-samples-include",
        filter_samples_include_.value,
        "Sample names to include (all other samples are excluded); either (1) a comma- or "
        "tab-separated list given on the command line, or (2) a file with one sample name per line. "
        "If no sample filter is provided, all samples in the input file are used. "
        "The option is applied after potentially renaming the samples with `--rename-samples-file`."
    );
    filter_samples_include_.option->group( group );

    // And the other way round.
    filter_samples_exclude_.option = sub->add_option(
        "--filter-samples-exclude",
        filter_samples_exclude_.value,
        "Sample names to exclude (all other samples are included); either (1) a comma- or "
        "tab-separated list given on the command line, or (2) a file with one sample name per line. "
        "If no sample filter is provided, all samples in the input file are used. "
        "The option is applied after potentially renaming the samples with `--rename-samples-file`."
    );
    filter_samples_exclude_.option->group( group );

    // Include and exclude are mutually exclusive.
    filter_samples_include_.option->excludes( filter_samples_exclude_.option );
    filter_samples_exclude_.option->excludes( filter_samples_include_.option );
}

// =================================================================================================
//      Run Functions
// =================================================================================================

// -------------------------------------------------------------------------
//     rename_samples
// -------------------------------------------------------------------------

void VariantSampleNamesOptions::rename_samples( std::vector<std::string>& sample_names ) const
{
    using namespace genesis::utils;

    // If the option is not set or used, we do no renaming.
    if( ! rename_samples_.option || ! *rename_samples_.option ) {
        return;
    }

    // Read the rename list from the provided file.
    std::unordered_map<std::string, std::string> rename_map;
    std::unordered_set<std::string> uniq_new_names;
    auto const lines = file_read_lines( rename_samples_.value );
    for( auto const& line : lines ) {
        auto const parts = split( line, "\t" );
        if( parts.size() != 2 ) {
            throw CLI::ValidationError(
                rename_samples_.option->get_name() + "(" + rename_samples_.value + ")",
                "Sample renaming list contains invalid lines. We expect an old and a new sample "
                "name, split by a tab character, but instead found \"" + line + "\""
            );
        }

        // Check for uniqueness of old and new sample names.
        if( rename_map.count( parts[0] ) > 0 ) {
            throw CLI::ValidationError(
                rename_samples_.option->get_name() + "(" + rename_samples_.value + ")",
                "Sample renaming list contains duplicate entries for old sample name \"" +
                parts[0] + "\""
            );
        }
        if( uniq_new_names.count( parts[1] ) > 0 ) {
            throw CLI::ValidationError(
                rename_samples_.option->get_name() + "(" + rename_samples_.value + ")",
                "Sample renaming list contains duplicate entries for new sample name \"" +
                parts[1] + "\", which would result in multiple samples with the same name."
            );
        }

        // Now we can add the entry to both sets.
        rename_map[ parts[0] ] = parts[1];
        uniq_new_names.insert( parts[1] );
    }
    assert( rename_map.size() == lines.size() );
    assert( uniq_new_names.size() == lines.size() );

    // Now rename, using the map. While doing so, we do duplication check of the old names as well,
    // just to be on the save side. Also, we keep track of all possibilities here,
    // for user convenience and error checking.
    std::unordered_set<std::string> uniq_old_names;
    size_t not_in_map = 0;
    size_t properly_renamed = 0;
    size_t same_old_and_new_name = 0;
    for( auto& name : sample_names ) {
        // Old name unique check
        if( uniq_old_names.count( name ) > 0 ) {
            throw std::runtime_error(
                "Cannot rename sample names, as input sample name \"" + name + "\" is not unique "
                "and appears multiple times in the input files."
            );
        }
        uniq_old_names.insert( name );

        // See if we can do a renaming properly.
        if( rename_map.count( name ) == 0 ) {
            ++not_in_map;
            continue;
        }
        if( rename_map[name] == name ) {
            ++same_old_and_new_name;
        } else {
            ++properly_renamed;
        }
        auto const old_name = name;
        name = rename_map[old_name];
        rename_map.erase( old_name );
    }
    assert( not_in_map + properly_renamed + same_old_and_new_name == sample_names.size() );

    // Now produce user output reporting on what we did here.
    LOG_MSG << "Renamed " << ( properly_renamed + same_old_and_new_name ) << " out of "
            << sample_names.size() << " sample names.";
    if( same_old_and_new_name > 0 ) {
        LOG_WARN << "Out of those, " << same_old_and_new_name << " names contained identical "
                 << "entries for the old and new name, which hence did not change anything.";
    }
    if( not_in_map > 0 ) {
        LOG_WARN << "Did not find old and new name pair for " << not_in_map << " samples "
                 << "in the rename list, which where hence not renamed.";
    }
    if( rename_map.size() > 0 ) {
        LOG_WARN << "Rename list contains " << rename_map.size() << " entries whose old name "
                 << "is not a sample name of the input samples, and which were hence not used.";
    }
}

// -------------------------------------------------------------------------
//     add_sample_name_filter
// -------------------------------------------------------------------------

void VariantSampleNamesOptions::add_sample_name_filter(
    genesis::population::VariantInputIterator& iterator
) const {
    using namespace genesis::population;

    // Only continue if we actually have a filter.
    if( filter_samples_include_.value.empty() && filter_samples_exclude_.value.empty() ) {
        return;
    }

    // Make filters as required.
    auto const& sample_names = iterator.data().sample_names;
    std::vector<bool> sample_filter;
    if( ! filter_samples_include_.value.empty() ) {
        auto const list = process_sample_name_list_option_( filter_samples_include_.value );
        sample_filter = make_sample_filter( sample_names, list, false );
    }
    if( ! filter_samples_exclude_.value.empty() ) {
        auto const list = process_sample_name_list_option_( filter_samples_exclude_.value );
        sample_filter = make_sample_filter( sample_names, list, true );
    }
    internal_check(
        sample_filter.size() == sample_names.size(), "sample_filter.size() != sample_names.size()"
    );

    // Checks and user information.
    size_t pop_count = 0;
    for( size_t i = 0; i < sample_filter.size(); ++i ) {
        if( sample_filter[i] ) {
            ++pop_count;
        }
    }
    if( pop_count == 0 ) {
        throw std::runtime_error( "Sample filter results in no samples being used." );
    }
    LOG_MSG << "Sample filter results in " << pop_count << " samples out of "
            << sample_filter.size() << " being used.";

    // Now add the transform that does the actual filtering.
    iterator.add_transform(
        make_variant_input_iterator_sample_name_filter_transform( sample_filter )
    );
}

// -------------------------------------------------------------------------
//     process_sample_name_list_option_
// -------------------------------------------------------------------------

std::vector<std::string> VariantSampleNamesOptions::process_sample_name_list_option_(
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
