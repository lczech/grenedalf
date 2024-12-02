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

#include "options/variant_sample_names.hpp"

#include "options/global.hpp"
#include "tools/misc.hpp"

#include "genesis/population/function/variant_input_stream.hpp"
#include "genesis/population/stream/variant_input_stream_adapters.hpp"
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
        "--rename-samples-list",
        rename_samples_.value,
        "Allows to rename samples, by providing a file that lists the original and new sample names, "
        "one per line, separating original and new names by a comma or tab.\n"
        "By default, we use sample names as provided in the input files. Some file types however do "
        "not contain sample names, such as (m)pileup or sync files (unless the non-standard sync "
        "header line is provided). For such file types, sample names are automatically assigned "
        "by using their input file base name (without path and extension), followed by a dot and "
        "numbers 1..n for all samples in that file. "
        "For instance, samples in `/path/to/sample.sync` are named `sample.1`, `sample.2`, etc.\n"
        "Using this option, those names can be renamed as needed. Use verbose output (`--verbose`) "
        "to show a list of all sample names. We then use these names in the output as well as in the "
        "`--filter-samples-include` and `--filter-samples-exclude` options. "
    );
    rename_samples_.option->group( group );
    rename_samples_.option->check( CLI::ExistingFile );

    // Add option for sample name filter.
    filter_samples_include_.option = sub->add_option(
        "--filter-samples-include",
        filter_samples_include_.value,
        "Sample names to include (all other samples are excluded); either (1) a comma- or "
        "tab-separated list given on the command line (in a typical shell, this list has to be "
        "enclosed in quotation marks), or (2) a file with one sample name per line. "
        "If no sample filter is provided, all samples in the input file are used. "
        "The option is applied after potentially renaming the samples with `--rename-samples-list`."
    );
    filter_samples_include_.option->group( group );

    // And the other way round.
    filter_samples_exclude_.option = sub->add_option(
        "--filter-samples-exclude",
        filter_samples_exclude_.value,
        "Sample names to exclude (all other samples are included); either (1) a comma- or "
        "tab-separated list given on the command line (in a typical shell, this list has to be "
        "enclosed in quotation marks), or (2) a file with one sample name per line. "
        "If no sample filter is provided, all samples in the input file are used. "
        "The option is applied after potentially renaming the samples with `--rename-samples-list`."
    );
    filter_samples_exclude_.option->group( group );

    // Include and exclude are mutually exclusive.
    filter_samples_include_.option->excludes( filter_samples_exclude_.option );
    filter_samples_exclude_.option->excludes( filter_samples_include_.option );

    // We allow to in-place merge samples that the user determines belong to the same group.
    sample_group_merge_table_file_.option = sub->add_option(
        "--sample-group-merge-table",
        sample_group_merge_table_file_.value,
        "When the input contains multiple samples (either within a single input file, or by "
        "providing multiple input files), these can be merged into new samples, by summing up "
        "their nucleotide base counts at each position. This has essentially the same effect as "
        "having merged the raw fastq files or the mapped sam/bam files of the samples, that is, "
        "all reads from those samples are treated as if they were a single sample. "
        "For this grouping, the option takes a simple table file (comma- or tab-separated), "
        "with the sample names (after the above renaming, if provided) in the first column, "
        "and their assigned group names in the second column. "
        "All samples in the same group are then merged into a grouped sample, and the group names "
        "are used as the new sample names for the output. Note that the `--pool-sizes` option "
        "then need to contain the summed up pool sizes for each group, using the group names."
    );
    sample_group_merge_table_file_.option->group( group );
    sample_group_merge_table_file_.option->check( CLI::ExistingFile );
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
        auto const parts = split( line, ",\t" );
        if( parts.size() != 2 ) {
            throw CLI::ValidationError(
                rename_samples_.option->get_name() + "(" + rename_samples_.value + ")",
                "Sample renaming list contains invalid lines. We expect an original and a new sample "
                "name, split by a comma or tab character, but instead found \"" + line + "\"" +
                ( parts.size() > 2 ? ", which contains multiple separator chars." : "." )
            );
        }

        // Check for uniqueness of original and new sample names.
        if( rename_map.count( parts[0] ) > 0 ) {
            throw CLI::ValidationError(
                rename_samples_.option->get_name() + "(" + rename_samples_.value + ")",
                "Sample renaming list contains duplicate entries for original sample name \"" +
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

    // Now rename, using the map. While doing so, we do duplication check of the original names as well,
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
                 << "entries for the original and new name, which hence did not change anything.";
    }
    if( not_in_map > 0 ) {
        LOG_WARN << "Did not find original and new name pair for " << not_in_map << " samples "
                 << "in the rename list, which where hence not renamed.";
    }
    if( rename_map.size() > 0 ) {
        LOG_WARN << "Rename list contains " << rename_map.size() << " entries whose original name "
                 << "is not a sample name of the input samples, and which were hence not used.";
    }
}

// -------------------------------------------------------------------------
//     add_sample_name_filter
// -------------------------------------------------------------------------

void VariantSampleNamesOptions::add_sample_name_filter(
    genesis::population::VariantInputStream& stream
) const {
    using namespace genesis::population;

    // Only continue if we actually have a filter.
    if( filter_samples_include_.value.empty() && filter_samples_exclude_.value.empty() ) {
        return;
    }

    // Make filters as required.
    auto const& sample_names = stream.data().sample_names;
    std::vector<bool> sample_filter;
    if( ! filter_samples_include_.value.empty() ) {
        auto const list = process_sample_name_list_option_( filter_samples_include_.value );
        sample_filter = make_sample_name_filter( sample_names, list, false );
    }
    if( ! filter_samples_exclude_.value.empty() ) {
        auto const list = process_sample_name_list_option_( filter_samples_exclude_.value );
        sample_filter = make_sample_name_filter( sample_names, list, true );
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
    stream.add_transform(
        make_variant_input_stream_sample_name_filter_transform( sample_filter )
    );

    // We also need to re-set the sample names, so that they only contain the used ones.
    std::vector<std::string> new_sample_names;
    new_sample_names.reserve( pop_count );
    for( size_t i = 0; i < sample_filter.size(); ++i ) {
        if( sample_filter[i] ) {
            new_sample_names.push_back( sample_names[i] );
        }
    }
    stream.data().sample_names = new_sample_names;
}

// -------------------------------------------------------------------------
//     apply_sample_group_merging
// -------------------------------------------------------------------------

void VariantSampleNamesOptions::apply_sample_group_merging(
    genesis::population::VariantInputStream& stream
) const {
    using namespace genesis::utils;

    // If the option is not set or used, we do no grouping.
    if( ! sample_group_merge_table_file_.option || ! *sample_group_merge_table_file_.option ) {
        return;
    }

    // Read the file line by line and process.
    auto const lines = file_read_lines( sample_group_merge_table_file_.value );
    std::unordered_map<std::string, std::string> sample_name_to_group;
    for( size_t i = 0; i < lines.size(); ++i ) {
        auto const& line = lines[i];

        // Dissect the line and see if we got a sample name and a group name.
        auto const pair = split( line, ",\t", false );
        if( pair.size() != 2 ) {
            throw CLI::ValidationError(
                sample_group_merge_table_file_.option->get_name() + "(" +
                sample_group_merge_table_file_.value + ")",
                "Invalid sample grouping file that contains an invalid line at " +
                std::to_string( i + 1 ) + " not consisting of a sample name and a group name."
            );
        }

        // For name value pairs, do a duplicate check first...
        if( sample_name_to_group.count( pair[0] ) > 0 ) {
            throw CLI::ValidationError(
                sample_group_merge_table_file_.option->get_name() + "(" +
                sample_group_merge_table_file_.value + ")",
                "Invalid line that contains duplicate sample names (line " +
                std::to_string( i + 1 ) + "): \"" + pair[0] + "\""
            );
        }

        // ... then add the entry to the map.
        assert( sample_name_to_group.count( pair[0] ) == 0 );
        sample_name_to_group[ pair[0] ] = pair[1];
    }

    // We inline replace the stream with a new one, that captures the existing one
    // in its lambda capture internally, and log a bit for the user.
    auto const orig_smp_cnt = stream.data().sample_names.size();
    stream = make_variant_merging_input_stream( stream, sample_name_to_group );
    LOG_MSG1
        << "Merging " << orig_smp_cnt << " samples into "
        << stream.data().sample_names.size() << " sample groups"
    ;
    if( orig_smp_cnt != sample_name_to_group.size() ) {
        LOG_MSG1
            << "Note that a total of " << sample_name_to_group.size() << " sample name assignments "
            << "are given, but not all appear in the input and are hence ignored"
        ;
    }
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
