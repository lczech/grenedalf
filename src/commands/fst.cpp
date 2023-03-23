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

#include "commands/fst.hpp"
#include "options/global.hpp"
#include "tools/cli_setup.hpp"
#include "tools/misc.hpp"

#include "genesis/population/functions/functions.hpp"
#include "genesis/population/functions/fst_pool_functions.hpp"
#include "genesis/population/functions/fst_pool_karlsson.hpp"
#include "genesis/population/functions/fst_pool_kofler.hpp"
#include "genesis/population/functions/fst_pool_unbiased.hpp"
#include "genesis/population/functions/fst_pool.hpp"
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
//      Enum Mapping
// =================================================================================================

enum class FstMethod
{
    kUnbiasedNei,
    kUnbiasedHudson,
    kKofler,
    kKarlsson
};

std::vector<std::pair<std::string, FstMethod>> const fst_method_map = {
    { "unbiased-nei",    FstMethod::kUnbiasedNei },
    { "unbiased-hudson", FstMethod::kUnbiasedHudson },
    { "kofler",          FstMethod::kKofler },
    { "karlsson",        FstMethod::kKarlsson }
};

// =================================================================================================
//      Setup
// =================================================================================================

void setup_fst( CLI::App& app )
{
    // Create the options and subcommand objects.
    auto options = std::make_shared<FstOptions>();
    auto sub = app.add_subcommand(
        "fst",
        "Compute FST in windows or at individual positions along the genome."
    );

    // -------------------------------------------------------------------------
    //     Input
    // -------------------------------------------------------------------------

    // Required input of some frequency format, and settings for the sliding window.
    options->variant_input.add_variant_input_opts_to_app( sub );
    // options->variant_input.add_sample_name_opts_to_app( sub );
    // options->variant_input.add_region_filter_opts_to_app( sub );
    options->window.add_window_opts_to_app( sub, true );

    // -------------------------------------------------------------------------
    //     Settings
    // -------------------------------------------------------------------------

    // Settings: FST Method
    options->method.option = sub->add_option(
        "--method",
        options->method.value,
        "FST method to use for the computation, either the unbiased pool-sequencing statistic "
        "(in two variants, following the definition of Nei, and the definition of Hudson et al), "
        "the statistic by Kofler et al of PoPoolation2, or the asymptotically unbiased "
        "estimator of Karlsson et al (which is also implemented in PoPoolation2)."
    );
    options->method.option->group( "Settings" );
    options->method.option->required();
    options->method.option->transform(
        CLI::IsMember( enum_map_keys( fst_method_map ), CLI::ignore_case )
    );

    // Settings: Pool sizes
    options->poolsizes.add_poolsizes_opt_to_app( sub, false );

    // TODO need settings for min/max coverage etc. see prototype implementations!

    // TODO need to check if this is still needed when filtering for snps...
    // Settings: Omit Empty Windows
    options->omit_na_windows.option = sub->add_flag(
        "--omit-na-windows",
        options->omit_na_windows.value,
        "Do not output windows where all values are n/a (e.g., without any SNPs). This is "
        "particularly relevant with small window sizes (or individual positions), "
        "in order to not produce output for invariant positions in the genome."
    );
    options->omit_na_windows.option->group( "Settings" );

    // Settings: Comparand
    options->comparand.option = sub->add_option(
        "--comparand",
        options->comparand.value,
        "By default, FST between all pairs of samples (that are not filtered) is computed. "
        "If this option is given a sample name however, only the pairwise FST between that "
        "sample and all others (that are not filtered) is computed."
    );
    options->comparand.option->group( "Settings" );

    // Settings: Second comparand
    options->second_comparand.option = sub->add_option(
        "--second-comparand",
        options->second_comparand.value,
        "If in addition to `--comparand`, this option is also given a (second) sample name, only "
        "FST between those two samples is computed."
    );
    options->second_comparand.option->group( "Settings" );
    options->second_comparand.option->needs( options->comparand.option );

    // Settings: Comparand list
    options->comparand_list.option = sub->add_option(
        "--comparand-list",
        options->comparand_list.value,
        "By default, FST between all pairs of samples is computed. If this option is given a file "
        "containing comma- or tab-separated pairs of sample names (one pair per line) however, "
        "only these pairwise FST values are computed."
    );
    options->comparand_list.option->group( "Settings" );
    options->comparand_list.option->check( CLI::ExistingFile );
    options->comparand_list.option->excludes( options->comparand.option );
    options->comparand_list.option->excludes( options->second_comparand.option );

    // -------------------------------------------------------------------------
    //     Output
    // -------------------------------------------------------------------------

    // Add table output options.
    options->table_output.add_separator_char_opt_to_app( sub );
    options->table_output.add_na_entry_opt_to_app( sub );

    // Output
    options->file_output.add_default_output_opts_to_app( sub );
    options->file_output.add_file_compress_opt_to_app( sub );

    // Set the run function as callback to be called when this subcommand is issued.
    // Hand over the options by copy, so that their shared ptr stays alive in the lambda.
    sub->callback( grenedalf_cli_callback(
        sub,
        {
            // Citation keys as needed
            "Kofler2011-PoPoolation2"
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
    auto const& sample_names = options.variant_input.sample_names();
    auto sample_index = [&]( std::string const& name ) -> size_t {
        auto it = std::find( sample_names.begin(), sample_names.end(), name );
        if( it == sample_names.end() ) {
            throw CLI::ValidationError(
                "Comparand sample names",
                "Invalid sample name: \"" + name  + "\" that was not found in the input, "
                "or was filtered out."
            );
        }
        return static_cast<size_t>( it - sample_names.begin() );
    };

    // Get all pairs of samples for which we want to compute FST.
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
            auto const pair = split( line, ",\t", false );
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

/**
 * @brief Helper, copied from genesis, but amended to allow extra args for constructing the Calculator.
 *
 * Need to back-port into genesis later.
 */
template<class Calculator, typename... Args>
inline genesis::population::FstPoolProcessor make_fst_pool_processor_(
    std::vector<size_t> const& pool_sizes,
    std::vector<std::pair<size_t, size_t>> const& sample_pairs,
    Args... args
) {
    using namespace genesis::population;

    FstPoolProcessor result;
    for( auto const& p : sample_pairs ) {
        if( p.first >= pool_sizes.size() || p.second >= pool_sizes.size() ) {
            throw std::invalid_argument(
                "Invalid sample indices for computing FST Pool Kofler: " +
                std::to_string( pool_sizes.size() ) +
                " pool sizes provided, but asked to use indices " +
                std::to_string( p.first ) + " and " + std::to_string( p.second )
            );
        }
        result.add_calculator(
            p.first, p.second,
            ::genesis::utils::make_unique<Calculator>(
                pool_sizes[p.first],
                pool_sizes[p.second],
                args...
            )
        );
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
    // using VariantWindow = Window<genesis::population::Variant>;

    // Output preparation.
    options.file_output.check_output_files_nonexistence( "fst", "csv" );

    // -------------------------------------------------------------------------
    //     Preparation
    // -------------------------------------------------------------------------

    // Use an enum for the method, which is faster to check in the main loop than doing
    // string comparisons all the time.
    auto const method = get_enum_map_value( fst_method_map, options.method.value );

    // Get indices of all pairs of samples for which we want to compute FST.
    auto const sample_pairs = get_sample_pairs_( options );

    // Get all sample indices that we are actually interested in.
    // We do this so that pool sizes for samples that we ignore anyway do not need to be given.
    auto const& sample_names = options.variant_input.sample_names();
    auto used_samples = std::vector<bool>( sample_names.size(), false );
    for( auto const& sp : sample_pairs ) {
        used_samples[ sp.first ]  = true;
        used_samples[ sp.second ] = true;
    }
    assert( used_samples.size() == sample_names.size() );

    // Get the pool sizes for all samples that we are interested in.
    auto const needs_pool_sizes = (
        method == FstMethod::kUnbiasedNei ||
        method == FstMethod::kUnbiasedHudson ||
        method == FstMethod::kKofler
    );
    auto const pool_sizes = (
        needs_pool_sizes
        ? options.poolsizes.get_pool_sizes( sample_names, used_samples )
        : std::vector<size_t>( sample_names.size() )
    );
    if( sample_pairs.empty() ) {
        LOG_WARN << "No pairs of samples selected, which will produce empty output. Stopping now.";
        return;
    }
    LOG_MSG << "Computing FST between " << sample_pairs.size()
            << " pair" << ( sample_pairs.size() > 1 ? "s" : "" ) << " of samples.";

    // Get the separator char to use for table entries.
    auto const sep_char = options.table_output.get_separator_char();

    // Make the processor.
    FstPoolProcessor processor;
    switch( method ) {
        case FstMethod::kUnbiasedNei: {
            processor = make_fst_pool_processor_<FstPoolCalculatorUnbiased>(
                pool_sizes, sample_pairs, FstPoolCalculatorUnbiased::Estimator::kNei
            );
            break;
        }
        case FstMethod::kUnbiasedHudson: {
            processor = make_fst_pool_processor_<FstPoolCalculatorUnbiased>(
                pool_sizes, sample_pairs, FstPoolCalculatorUnbiased::Estimator::kHudson
            );
            break;
        }
        case FstMethod::kKofler: {
            processor = make_fst_pool_processor_<FstPoolCalculatorKofler>(
                pool_sizes, sample_pairs
            );
            break;
        }
        case FstMethod::kKarlsson: {
            processor = make_fst_pool_processor_<FstPoolCalculatorKarlsson>(
                pool_sizes, sample_pairs
            );
            break;
        }
        default: {
            throw std::domain_error( "Internal error: Invalid FST method." );
        }
    }

    // -------------------------------------------------------------------------
    //     Table Header
    // -------------------------------------------------------------------------

    // Prepare output file and write header line with all pairs of samples.
    auto fst_ofs = options.file_output.get_output_target( "fst", "csv" );
    (*fst_ofs) << "CHROM" << sep_char << "START" << sep_char << "END" << sep_char << "SNPS";
    for( auto const& pair : sample_pairs ) {
        (*fst_ofs) << sep_char << sample_names[pair.first] << "." << sample_names[pair.second];
    }
    (*fst_ofs) << "\n";

    // For our unbiased and for Kofler, check that we got the right number of pool sizes.
    internal_check(
        ! needs_pool_sizes || pool_sizes.size() == sample_names.size(),
        "Inconsistent number of samples and number of pool sizes."
    );

    // -------------------------------------------------------------------------
    //     Main Loop
    // -------------------------------------------------------------------------

    // Iterate the file and compute per-window FST.
    size_t chr_cnt = 0;
    size_t win_cnt = 0;
    size_t pos_cnt = 0;
    size_t nan_cnt = 0;

    auto window_it = options.window.get_variant_window_view_iterator(
        options.variant_input.get_iterator()
    );
    auto window_fst = std::vector<double>( sample_pairs.size() );
    for( auto cur_it = window_it->begin(); cur_it != window_it->end(); ++cur_it ) {
        auto const& window = *cur_it;

        // Reset all calculator accumulators to zero for this window.
        processor.reset();

        // Compute diversity over samples.
        // #pragma omp parallel for
        size_t entry_count = 0;
        for( auto const& variant : window ) {
            internal_check(
                variant.samples.size() == sample_names.size(),
                "Inconsistent number of samples in input file."
            );

            processor.process( variant );
            ++pos_cnt;
            ++entry_count;
        }
        auto const& window_fst = processor.get_result();

        // Write the values, unless all of them are nan, then we might skip.
        // Skip empty windows if the user wants to.
        if(
            options.omit_na_windows.value && (
                entry_count == 0 ||
                std::none_of( window_fst.begin(), window_fst.end(), []( double v ) {
                    return std::isfinite( v );
                })
            )
        ) {
            ++nan_cnt;
        } else {
            ++win_cnt;

            // Write fixed columns.
            (*fst_ofs) << window.chromosome();
            (*fst_ofs) << sep_char << window.first_position();
            (*fst_ofs) << sep_char << window.last_position();
            (*fst_ofs) << sep_char << entry_count;

            // Write the per-pair FST values in the correct order.
            for( auto const& fst : window_fst ) {
                if( std::isfinite( fst ) ) {
                    (*fst_ofs) << sep_char << fst;
                } else {
                    (*fst_ofs) << sep_char << options.table_output.get_na_entry();
                }
            }
            (*fst_ofs) << "\n";
        }
    }

    // Final user output.
    LOG_MSG << "\nProcessed " << chr_cnt << " chromosome" << ( chr_cnt != 1 ? "s" : "" )
            << " with " << pos_cnt << " total position" << ( pos_cnt != 1 ? "s" : "" )
            << " in " << win_cnt << " window" << ( win_cnt != 1 ? "s" : "" )
            << " with FST values, and skipped " << nan_cnt << " window"
            << ( nan_cnt != 1 ? "s" : "" ) << " without any FST values.";
}
