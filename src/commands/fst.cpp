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

#include "commands/fst.hpp"
#include "options/global.hpp"
#include "tools/cli_setup.hpp"
#include "tools/misc.hpp"

#include "genesis/population/functions/base_counts.hpp"
#include "genesis/population/functions/structure.hpp"
#include "genesis/utils/containers/transform_iterator.hpp"
#include "genesis/utils/core/fs.hpp"
#include "genesis/utils/text/convert.hpp"
#include "genesis/utils/text/string.hpp"

#include <algorithm>
#include <cassert>
#include <unordered_map>
#include <utility>
#include <vector>

#ifdef GENESIS_OPENMP
#   include <omp.h>
#endif

// =================================================================================================
//      Setup
// =================================================================================================

void setup_fst( CLI::App& app )
{
    // Create the options and subcommand objects.
    auto options = std::make_shared<FstOptions>();
    auto sub = app.add_subcommand(
        "fst",
        "Compute F_ST in windows or at individual positions along the genome."
    );

    // Required input of some frequency format (mpileup or vcf at the moment),
    // and settings for the sliding window.
    options->freq_input.add_frequency_input_opts_to_app( sub );
    options->freq_input.add_sliding_window_opts_to_app( sub );

    // Settings: Pool Sizes
    options->poolsizes.option = sub->add_option(
        "--poolsizes",
        options->poolsizes.value,
        "Poolsizes for all samples for which F_ST is to be computed. Either a single number that "
        "is used for all samples, or a path to a file that contains a tab-separated list of sample "
        "names and pool sizes, with one sample/size entry per line."
    );
    options->poolsizes.option->group( "Settings" );
    // options->poolsizes.option->required();

    // Settings: F_ST Method
    options->method.option = sub->add_option(
        "--method",
        options->method.value,
        "F_ST method to use for the computation, either the conventional F_ST statistic for "
        "pool-sequenced data, following Kofler et al, or the asymptotically unbiased F_ST "
        "estimator of Karlsson et al."
    );
    options->method.option->group( "Settings" );
    options->method.option->transform(
        CLI::IsMember({ "conventional", "karlsson" }, CLI::ignore_case )
    );

    // Settings: Omit Empty Windows
    options->omit_empty_windows.option = sub->add_flag(
        "--omit-empty-windows",
        options->omit_empty_windows.value,
        "Do not output empty windows (without any SNPs). This is particularly relevant when choosing "
        "`--window-width 1` (or other small window sizes), in order to not produce output for "
        "every position in the genome."
    );
    options->omit_empty_windows.option->group( "Settings" );

    // Settings: Comparand
    options->comparand.option = sub->add_option(
        "--comparand",
        options->comparand.value,
        "By default, F_ST between all pairs of samples (that are not filtered) is computed. "
        "If this option is given a sample name however, only the pairwise F_ST between that "
        "sample and all others (that are not filtered) is computed."
    );
    options->comparand.option->group( "Settings" );

    // Settings: Second comparand
    options->second_comparand.option = sub->add_option(
        "--second-comparand",
        options->second_comparand.value,
        "If in addition to `--comparand`, this option is also given a (second) sample name, only "
        "F_ST between those two samples is computed."
    );
    options->second_comparand.option->group( "Settings" );
    options->second_comparand.option->needs( options->comparand.option );

    // Settings: Comparand list
    options->comparand_list.option = sub->add_option(
        "--comparand-list",
        options->comparand_list.value,
        "By default, F_ST between all pairs of samples is computed. If this option is given a file "
        "containing tab-separated pairs of sample names (one pair per line) however, only these "
        "pairwise F_ST values are computed."
    );
    options->comparand_list.option->group( "Settings" );
    options->comparand_list.option->excludes( options->comparand.option );
    options->comparand_list.option->excludes( options->second_comparand.option );

    // Add an option to set the text for not-a-number.
    options->na_entry.option = sub->add_option(
        "--na-entry",
        options->na_entry.value,
        "Set the text to use for n/a entries (that is, windows that contain no values). "
        "This is useful to match formatting expectations of downstream software."
    )->group( "Settings" );

    // Output
    options->file_output.add_default_output_opts_to_app( sub );
    options->file_output.add_file_compress_opt_to_app( sub );

    // Set the run function as callback to be called when this subcommand is issued.
    // Hand over the options by copy, so that their shared ptr stays alive in the lambda.
    sub->callback( grenedalf_cli_callback(
        sub,
        {
            // Citation keys as needed
        },
        [ options ]() {
            run_fst( *options );
        }
    ));
}

// =================================================================================================
//      Helper Functions
// =================================================================================================

std::vector<std::pair<size_t, size_t>> get_sample_pairs_( FstOptions const& options )
{
    using namespace genesis::utils;

    // Get sample names and map to their index or throw if invalid name is given.
    auto const& sample_names = options.freq_input.sample_names();
    auto sample_index = [&]( std::string const& name ) -> size_t {
        auto it = std::find( sample_names.begin(), sample_names.end(), name );
        if( it == sample_names.end() ) {
            throw CLI::ValidationError(
                "Invalid sample name: \"" + name  + "\" that was not found in the input, "
                "or was filtered out."
            );
        }
        return static_cast<size_t>( it - sample_names.begin() );
    };

    // Get all pairs of samples for which we want to compute F_ST.
    std::vector<std::pair<size_t, size_t>> sample_pairs;
    if( *options.comparand.option ) {
        if( *options.second_comparand.option ) {
            // Only exactly one pair of samples.
            auto const index_a = sample_index( options.comparand.value );
            auto const index_b = sample_index( options.second_comparand.value );
            sample_pairs.emplace_back( index_a, index_b );
        } else {
            // One sample against all others.
            auto const index_a = sample_index( options.comparand.value );
            for( auto const& sn : sample_names ) {
                if( sn != options.comparand.value ) {
                    auto const index_b = sample_index( sn );
                    sample_pairs.emplace_back( index_a, index_b );
                }
            }
        }
    } else if( *options.comparand_list.option ) {
        // Read list of pairs from file.
        auto const lines = file_read_lines( options.comparand_list.value );
        for( size_t i = 0; i < lines.size(); ++i ) {
            auto const& line = lines[i];
            auto const pair = split( line, "\t", false );
            if( pair.size() != 2 ) {
                throw CLI::ValidationError(
                    options.comparand_list.option->get_name() + "(" +
                    options.comparand_list.value + ")",
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
    return sample_pairs;
}

std::vector<size_t> get_pool_sizes(
    FstOptions const& options,
    std::vector<std::pair<size_t, size_t>> const& sample_pairs
) {
    using namespace genesis::utils;

    // Convert a pool size to a number, or throw.
    auto convert_poolsize_ = [&]( std::string const& str ){
        try {
            return convert_from_string<size_t>( str );
        } catch(...) {
            throw CLI::ValidationError(
                options.poolsizes.option->get_name(),
                "Invalid pool size value: " + str
            );
        }
    };

    // Get the pool sizes, depending on input type either from file, or as a single number.
    auto const& sample_names = options.freq_input.sample_names();
    auto result = std::vector<size_t>( sample_names.size(), 0 );
    if( is_file( options.poolsizes.value )) {

        // Read the file line by line and make a map from sample names to pool sizes.
        auto const lines = file_read_lines( options.poolsizes.value );
        std::unordered_map<std::string, size_t> ps_map;
        for( size_t i = 0; i < lines.size(); ++i ) {
            auto const& line = lines[i];

            // Dissect the line and see if we got a sample name and a number.
            auto const pair = split( line, "\t", false );
            if( pair.size() != 2 ) {
                throw CLI::ValidationError(
                    options.poolsizes.option->get_name() + "(" +
                    options.poolsizes.value + ")",
                    "Invalid line that does not contain sample name and a pool size (line " +
                    std::to_string( i + 1 ) + ")."
                );
            }

            // Duplicate check.
            if( ps_map.count( pair[0] )) {
                throw CLI::ValidationError(
                    options.poolsizes.option->get_name() + "(" +
                    options.poolsizes.value + ")",
                    "Invalid line that does not contain sample name and a pool size (line " +
                    std::to_string( i + 1 ) + ")."
                );
            }

            // Add the entry.
            assert( ps_map.count( pair[0] ) == 0 );
            ps_map[ pair[0] ] = convert_poolsize_( pair[1] );
        }

        // Get all sample indices that we are actually interested in.
        // We do this so that pool sizes for samples that we ignore anyway do not need to be given.
        auto used_samples = std::vector<bool>( sample_names.size(), false );
        for( auto const& sp : sample_pairs ) {
            used_samples[ sp.first ]  = true;
            used_samples[ sp.second ] = true;
        }

        // Now, fill the vector of pool sizes with values in the correct order
        // (that is, using the sample name order).
        // Throw if the sample needs to be given (it's in the used list), but isn't in the file.
        for( size_t i = 0; i < sample_names.size(); ++i ) {
            if( ! used_samples[i] ) {
                continue;
            }

            auto const& sn = sample_names[i];
            if( ps_map.count(sn) > 0 ) {
                result[i] = ps_map[sn];
            } else {
                throw CLI::ValidationError(
                    options.poolsizes.option->get_name() + "(" +
                    options.poolsizes.value + ")",
                    "Sample name \"" + sn +  "\" missing from pool size file."
                );
            }
        }
    } else {

        // Just give every sample the same pool size.
        auto const ps = convert_poolsize_( options.poolsizes.value );
        for( auto& entry : result ) {
            entry = ps;
        }
    }
    return result;
}

// =================================================================================================
//      Run
// =================================================================================================

void run_fst( FstOptions const& options )
{
    using namespace genesis::population;
    using namespace genesis::utils;
    using BaseCountWindow = Window<std::vector<BaseCounts>>;

    // Output preparation.
    options.file_output.check_output_files_nonexistence( "fst", "csv" );

    // Use an enum for the method, which is faster to check in the main loop than doing
    // string comparisons all the time. We could use a bool here, but let's be prepared for
    // any additional future F_ST methods.
    enum class Method
    {
        kConventional,
        kKarlsson
    };
    Method method = ( options.method.value == "conventional"
        ? Method::kConventional
        : Method::kKarlsson
    );

    // Get indices of all pairs of samples for which we want to compute F_ST,
    // and their pool sizes.
    auto const sample_pairs = get_sample_pairs_( options );
    auto const pool_sizes = ( method == Method::kConventional
        ? get_pool_sizes( options, sample_pairs )
        : std::vector<size_t>{}
    );
    LOG_MSG << "Computing F_ST between " << sample_pairs.size() << " pairs of samples.";

    // Prepare output file and write header line with all pairs of samples.
    auto fst_ofs = options.file_output.get_output_target( "fst", "csv" );
    auto const& sample_names = options.freq_input.sample_names();
    (*fst_ofs) << "CHROM\tSTART\tEND\tSNPS";
    for( auto const& pair : sample_pairs ) {
        (*fst_ofs) << "\t" << sample_names[pair.first] << "." << sample_names[pair.second];
    }
    (*fst_ofs) << "\n";

    // For conventional F_ST, check that we got the right number of pool sizes.
    internal_check(
        method != Method::kConventional ||
        pool_sizes.size() == sample_names.size(),
        "Inconsistent number of samples and number of pool sizes."
    );

    // Iterate the file and compute per-window F_ST.
    auto window_it = options.freq_input.get_base_count_sliding_window_iterator();
    auto window_fst = std::vector<double>( sample_pairs.size() );
    for( ; window_it; ++window_it ) {
        auto const& window = *window_it;

        // Some user output to report progress.
        if( window_it.is_first_window() ) {
            LOG_MSG << "At chromosome " << window.chromosome();
        }

        // Skip empty windows if the user wants to.
        if( window.empty() && options.omit_empty_windows.value ) {
            continue;
        }

        // Write fixed columns.
        (*fst_ofs) << window.chromosome();
        (*fst_ofs) << "\t" << window.first_position();
        (*fst_ofs) << "\t" << window.last_position();
        (*fst_ofs) << "\t" << window.entry_count();

        // Compute F_ST in parallel over the different pairs of samples.
        #pragma omp parallel for
        for( size_t i = 0; i < sample_pairs.size(); ++i ) {
            auto const index_a = sample_pairs[i].first;
            auto const index_b = sample_pairs[i].second;

            // The window contains entries from all samples in the input file, but our fst function
            // expects iterators over two ranges of BaseCounts for which it computes fst.
            // Hence, we here select the respective entries from the window BaseCounts vector.
            // We use a range transformer that selects the respective entry from the Variant,
            // and returns by reference, so that we avoid expensive copies here.
            auto range_a = make_transform_range(
                [index_a]( BaseCountWindow::Entry const& entry ) -> BaseCounts const& {
                    internal_check(
                        index_a < entry.data.size(),
                        "Inconsistent number of samples in input file."
                    );
                    return entry.data[index_a];
                },
                window.begin(), window.end()
            );
            auto range_b = make_transform_range(
                [index_b]( BaseCountWindow::Entry const& entry ) -> BaseCounts const& {
                    internal_check(
                        index_b < entry.data.size(),
                        "Inconsistent number of samples in input file."
                    );
                    return entry.data[index_b];
                },
                window.begin(), window.end()
            );

            // Run the computation.
            if( method == Method::kConventional ) {
                window_fst[i] = f_st_conventional_pool(
                    pool_sizes[index_a], pool_sizes[index_b],
                    range_a.begin(), range_a.end(),
                    range_b.begin(), range_b.end()
                );
            } else if( method == Method::kKarlsson ) {
                window_fst[i] = f_st_asymptotically_unbiased(
                    range_a.begin(), range_a.end(),
                    range_b.begin(), range_b.end()
                );
            } else {
                throw std::domain_error( "Internal error: Invalid F_ST method." );
            }
        }

        // Write the per-pair F_ST values in the correct order.
        for( auto const& fst : window_fst ) {
            if( std::isfinite( fst ) ) {
                (*fst_ofs) << "\t" << fst;
            } else {
                (*fst_ofs) << "\t" << options.na_entry.value;
            }
        }
        (*fst_ofs) << "\n";
    }
}
