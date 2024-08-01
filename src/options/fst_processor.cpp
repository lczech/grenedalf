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

#include "options/fst_processor.hpp"

#include "options/global.hpp"
#include "tools/cli_setup.hpp"
#include "tools/misc.hpp"

#include "genesis/population/function/fst_pool_calculator.hpp"
#include "genesis/population/function/fst_pool_functions.hpp"
#include "genesis/population/function/fst_pool_karlsson.hpp"
#include "genesis/population/function/fst_pool_kofler.hpp"
#include "genesis/population/function/fst_pool_unbiased.hpp"
#include "genesis/population/function/functions.hpp"
#include "genesis/utils/containers/matrix/operators.hpp"
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
//      Enum Mapping
// =================================================================================================

std::vector<std::pair<std::string, FstProcessorOptions::FstMethod>> const fst_method_map_all = {
    { "unbiased-nei",    FstProcessorOptions::FstMethod::kUnbiasedNei },
    { "unbiased-hudson", FstProcessorOptions::FstMethod::kUnbiasedHudson },
    { "kofler",          FstProcessorOptions::FstMethod::kKofler },
    { "karlsson",        FstProcessorOptions::FstMethod::kKarlsson }
};

std::vector<std::pair<std::string, FstProcessorOptions::FstMethod>> const fst_method_map_unbiased = {
    { "unbiased-nei",    FstProcessorOptions::FstMethod::kUnbiasedNei },
    { "unbiased-hudson", FstProcessorOptions::FstMethod::kUnbiasedHudson },
};

// =================================================================================================
//      Setup Functions
// =================================================================================================

void FstProcessorOptions::add_fst_processor_opts_to_app(
    CLI::App* sub,
    VariantReferenceGenomeOptions const& ref_genome_opts,
    Params const& fst_processor_params,
    std::string const& group
) {

    // -------------------------------------------------------------------------
    //     Basic Fst Settings
    // -------------------------------------------------------------------------

    // Window averaging
    if( fst_processor_params.with_window_average_policy ) {
        window_average_policy.add_window_average_opt_to_app(
            sub, ref_genome_opts, !fst_processor_params.with_all_methods
        );
    }

    // Settings: FST Method
    if( fst_processor_params.with_all_methods ) {
        method.option = sub->add_option(
            "--method",
            method.value,
            "FST method to use for the computation.\n(1) The unbiased pool-sequencing statistic "
            "(in two variants, following the definition of Nei, and the definition of Hudson),"
            "\n(2) the statistic by Kofler et al of PoPoolation2, or \n(3) the asymptotically "
            "unbiased estimator of Karlsson et al (which is also implemented in PoPoolation2). "
            "\nAll except for the Karlsson method also require `--pool-sizes` to be provided."
        );
        method.option->transform(
            CLI::IsMember( enum_map_keys( fst_method_map_all ), CLI::ignore_case )
        );
    } else {
        method.option = sub->add_option(
            "--method",
            method.value,
            "FST method to use for the computation: The unbiased pool-sequencing statistic "
            "in two variants, following the definition of Nei, and the definition of Hudson."
        );
        method.option->transform(
            CLI::IsMember( enum_map_keys( fst_method_map_unbiased ), CLI::ignore_case )
        );
    }
    method.option->group( group );
    method.option->required();

    // Pool sizes. When not all methods are selected, we only offer the unbiased ones,
    // which do need pool sizes, so then we require this to be provided.
    poolsizes.add_poolsizes_opt_to_app( sub, !fst_processor_params.with_all_methods, group );

    // Store for later
    params = fst_processor_params;

    // -------------------------------------------------------------------------
    //     Sample Pairs
    // -------------------------------------------------------------------------

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

    // -------------------------------------------------------------------------
    //     Misc
    // -------------------------------------------------------------------------

    // Setting: threading_threshold
    threading_threshold.option = sub->add_option(
        "--threading-threshold",
        threading_threshold.value,
        "When only computing FST beween a few pairs of samples, using individual threads for each "
        "usually incurs a substantial overhead due to thread synchronozation. Hence, we only want "
        "to use threads for the computation if many pairs of samples are being processed. "
        "This setting determiens the number of sample pairs at which threads are used. "
        "(Note that we still always use threads for input file parsing.)"
    );
    threading_threshold.option->group( "" );
}

// =================================================================================================
//      Run Functions
// =================================================================================================

// -------------------------------------------------------------------------
//     Properties and Getters
// -------------------------------------------------------------------------

bool FstProcessorOptions::is_all_to_all() const
{
    return ( ! *comparand.option ) && ( ! *comparand_list.option );
}

FstProcessorOptions::FstMethod FstProcessorOptions::get_fst_method() const
{
    // Independent of whether we use all or just the unbiased methods,
    // we can resolve the string to enum using the all map.
    // No need to switch to the fst_method_map_unbiased, as it is a subset.
    return get_enum_map_value( fst_method_map_all, method.value );
}

// -------------------------------------------------------------------------
//     get_sample_pairs
// -------------------------------------------------------------------------

// We want to build a hash map for pairs, which needs hashing first...
struct pair_hash {
    template <class T1, class T2>
    std::size_t operator () (const std::pair<T1,T2> &p) const {
        auto h1 = std::hash<T1>{}(p.first);
        auto h2 = std::hash<T2>{}(p.second);
        return genesis::utils::hash_combine( h1, h2 );
    }
};

std::vector<std::pair<size_t, size_t>> FstProcessorOptions::get_sample_pairs(
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

    // User output
    if( sample_pairs.empty() ) {
        LOG_WARN << "No pairs of samples selected, which will produce empty output. Stopping now.";
        return sample_pairs;
    }
    LOG_MSG << "Computing FST between " << sample_pairs.size()
            << " pair" << ( sample_pairs.size() > 1 ? "s" : "" ) << " of samples.";

    // If we do all-to-all, we expect a triangular matrix size of entries.
    internal_check(
        ! is_all_to_all() || sample_pairs.size() == triangular_size( sample_names.size() ),
        "Expeting all-to-all FST to create a trinangular number of sample pairs."
    );

    return sample_pairs;
}

// -------------------------------------------------------------------------
//     get_fst_pool_processor
// -------------------------------------------------------------------------

genesis::population::FstPoolProcessor FstProcessorOptions::get_fst_pool_processor(
    std::vector<std::string> const& sample_names,
    std::vector<std::pair<size_t, size_t>> const& sample_pairs
) const {
    using namespace genesis::population;

    // Get the method that we want to use.
    auto const method_value = get_fst_method();

    // Get all sample indices that we are actually interested in.
    // We do this so that pool sizes for samples that we ignore anyway do not need to be given.
    auto used_samples = std::vector<bool>( sample_names.size(), false );
    for( auto const& sp : sample_pairs ) {
        used_samples[ sp.first ]  = true;
        used_samples[ sp.second ] = true;
    }
    assert( used_samples.size() == sample_names.size() );

    // Get the pool sizes for all samples that we are interested in.
    auto const needs_pool_sizes = (
        method_value == FstMethod::kUnbiasedNei ||
        method_value == FstMethod::kUnbiasedHudson ||
        method_value == FstMethod::kKofler
    );
    auto const pool_sizes = (
        needs_pool_sizes
        ? poolsizes.get_pool_sizes( sample_names, used_samples )
        : std::vector<size_t>( sample_names.size(), 0 )
    );

    // Check the window averaging setup
    auto const needs_window_average_policy = (
        method_value == FstMethod::kUnbiasedNei ||
        method_value == FstMethod::kUnbiasedHudson
    );
    switch( method_value ) {
        case FstMethod::kUnbiasedNei:
        case FstMethod::kUnbiasedHudson: {
            if(
                ! window_average_policy.get_window_average_policy_option().option ||
                !*window_average_policy.get_window_average_policy_option().option
            ) {
                throw CLI::ValidationError(
                    method.option->get_name() + ", " +
                    window_average_policy.get_window_average_policy_option().option->get_name(),
                    "Window average policy option needs to be provided with " +
                    method.option->get_name() + " " + method.value
                );
            }
            break;
        }
        case FstMethod::kKofler:
        case FstMethod::kKarlsson: {
            if(
                 window_average_policy.get_window_average_policy_option().option &&
                *window_average_policy.get_window_average_policy_option().option
            ) {
                throw CLI::ValidationError(
                    method.option->get_name() + ", " +
                    window_average_policy.get_window_average_policy_option().option->get_name(),
                    "Window average policy option cannot be used with " +
                    method.option->get_name() + " " + method.value
                );
            }
            break;
        }
        default: {
            throw std::domain_error( "Internal error: Invalid FST method." );
        }
    }

    // For our unbiased and for Kofler, check that we got the right number of pool sizes.
    internal_check(
        ! needs_pool_sizes || pool_sizes.size() == sample_names.size(),
        "Inconsistent number of samples and number of pool sizes."
    );

    // Get the window average policy to use for all processors.
    auto const win_avg_policy = (
        params.with_window_average_policy && needs_window_average_policy
        ? window_average_policy.get_window_average_policy()
        : params.fix_window_average_policy
    );

    // Make the type of processor that we need for the provided method.
    FstPoolProcessor processor;
    switch( method_value ) {
        case FstMethod::kUnbiasedNei: {
            processor = make_fst_pool_processor<FstPoolCalculatorUnbiased>(
                sample_pairs, pool_sizes, win_avg_policy,
                FstPoolCalculatorUnbiased::Estimator::kNei
            );
            break;
        }
        case FstMethod::kUnbiasedHudson: {
            processor = make_fst_pool_processor<FstPoolCalculatorUnbiased>(
                sample_pairs, pool_sizes, win_avg_policy,
                FstPoolCalculatorUnbiased::Estimator::kHudson
            );
            break;
        }
        case FstMethod::kKofler: {
            processor = make_fst_pool_processor<FstPoolCalculatorKofler>(
                sample_pairs, pool_sizes
            );
            break;
        }
        case FstMethod::kKarlsson: {
            processor = make_fst_pool_processor<FstPoolCalculatorKarlsson>(
                sample_pairs, pool_sizes
            );
            break;
        }
        default: {
            throw std::domain_error( "Internal error: Invalid FST method." );
        }
    }

    // Set the threading options as provided by the (currently hidden) setting.
    processor.thread_pool( global_options.thread_pool() );
    processor.threading_threshold( threading_threshold.value );
    return processor;
}
