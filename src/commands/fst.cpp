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

#include "genesis/population/functions/fst_pool_functions.hpp"
#include "genesis/population/functions/fst_pool_karlsson.hpp"
#include "genesis/population/functions/fst_pool_kofler.hpp"
#include "genesis/population/functions/fst_pool_unbiased.hpp"
#include "genesis/population/functions/fst_pool.hpp"
#include "genesis/population/functions/functions.hpp"
#include "genesis/utils/containers/matrix.hpp"
#include "genesis/utils/containers/matrix/operators.hpp"
#include "genesis/utils/containers/matrix/writer.hpp"
#include "genesis/utils/containers/transform_iterator.hpp"
#include "genesis/utils/core/algorithm.hpp"
#include "genesis/utils/core/fs.hpp"
#include "genesis/utils/core/std.hpp"
#include "genesis/utils/text/convert.hpp"
#include "genesis/utils/text/string.hpp"

#include <algorithm>
#include <cassert>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

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

    // Required input of some count/frequency file format.
    options->variant_input.add_variant_input_opts_to_app( sub );
    // options->variant_input.add_sample_name_opts_to_app( sub );
    // options->variant_input.add_region_filter_opts_to_app( sub );

    // Individual settings for numerical filtering.
    // We do not add the only-SNP filter here as a user option, as we always filter out invariant
    // (non-SNP) sites below anyway. They do not contribute to the FST computation.
    options->filter_numerical.add_sample_count_filter_opts_to_app( sub );
    options->filter_numerical.add_sample_coverage_filter_opts_to_app( sub );
    options->filter_numerical.add_total_coverage_filter_opts_to_app( sub );
    options->filter_numerical.add_total_snp_filter_opts_to_app( sub, false, true );
    options->filter_numerical.add_total_freq_filter_opts_to_app( sub );

    // We amend the existing biallelic SNP filter option description here.
    auto& tobs = options->filter_numerical.total_only_biallelic_snps.option;
    tobs->description(
        tobs->get_description() + "\nBy default, we already filter out invariant sites, so that "
        "only SNPs remain. With this option however, this is further reduced to only biallelic "
        "(pure) SNPs. Note that with `--method karlsson`, we implicitly also activate to filter "
        "for biallelic SNPs only, as the method requires it."
    );

    // Settings for the windowing.
    options->window.add_window_opts_to_app( sub, true );

    // -------------------------------------------------------------------------
    //     Settings
    // -------------------------------------------------------------------------

    // Settings: FST Method
    options->method.option = sub->add_option(
        "--method",
        options->method.value,
        "FST method to use for the computation.\n(1) The unbiased pool-sequencing statistic "
        "(in two variants, following the definition of Nei, and the definition of Hudson et al),"
        "\n(2) the statistic by Kofler et al of PoPoolation2, or \n(3) the asymptotically unbiased "
        "estimator of Karlsson et al (which is also implemented in PoPoolation2). "
        "\nAll except for the Karlsson method also require `--pool-sizes` to be provided."
    );
    options->method.option->group( "Settings" );
    options->method.option->required();
    options->method.option->transform(
        CLI::IsMember( enum_map_keys( fst_method_map ), CLI::ignore_case )
    );

    // Settings: Pool sizes
    options->poolsizes.add_poolsizes_opt_to_app( sub, false );

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

// -------------------------------------------------------------------------
//     get_sample_pairs_
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

std::vector<std::pair<size_t, size_t>> get_sample_pairs_( FstOptions const& options )
{
    using namespace genesis::utils;

    // We here build a vector of pairs, in order to keep the order of pairs as specified.
    // This is important, as otherwise our assignmen of which samples and pairs belong
    // to which value computed will get messed up.

    // Get sample names and map to their index or throw if invalid name is given.
    // We create a map of names to their index in the list; if used, the names need to be unique.
    // We however only need it when a comparand of sorts is provided, so that we can look up the
    // names. Without comaparnd, we do not need the names, and hence can allow duplicates.
    // We hence use a lambda to only fill the map when needed.
    auto const& sample_names = options.variant_input.sample_names();
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
    if( *options.comparand.option ) {
        make_sample_map();
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
        // Read list of pairs from file, and prepare the name index lookup.
        auto const lines = file_read_lines( options.comparand_list.value );
        make_sample_map();

        // Get all pairs from the file, and build the list.
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

// -------------------------------------------------------------------------
//     get_fst_pool_processor_
// -------------------------------------------------------------------------

genesis::population::FstPoolProcessor get_fst_pool_processor_(
    FstOptions const& options,
    FstMethod method,
    std::vector<std::pair<size_t, size_t>> const& sample_pairs
) {
    using namespace genesis::population;

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

    // For our unbiased and for Kofler, check that we got the right number of pool sizes.
    internal_check(
        ! needs_pool_sizes || pool_sizes.size() == sample_names.size(),
        "Inconsistent number of samples and number of pool sizes."
    );

    switch( method ) {
        case FstMethod::kUnbiasedNei: {
            return make_fst_pool_processor<FstPoolCalculatorUnbiased>(
                sample_pairs, pool_sizes, FstPoolCalculatorUnbiased::Estimator::kNei
            );
        }
        case FstMethod::kUnbiasedHudson: {
            return make_fst_pool_processor<FstPoolCalculatorUnbiased>(
                sample_pairs, pool_sizes, FstPoolCalculatorUnbiased::Estimator::kHudson
            );
        }
        case FstMethod::kKofler: {
            return make_fst_pool_processor<FstPoolCalculatorKofler>(
                sample_pairs, pool_sizes
            );
        }
        case FstMethod::kKarlsson: {
            return make_fst_pool_processor<FstPoolCalculatorKarlsson>(
                sample_pairs, pool_sizes
            );
        }
        default: {
            throw std::domain_error( "Internal error: Invalid FST method." );
        }
    }
}

// -------------------------------------------------------------------------
//     print_fst_list_
// -------------------------------------------------------------------------

void print_fst_list_(
    FstOptions const& options,
    std::vector<double> const& window_fst,
    std::vector<std::pair<size_t, size_t>> const& sample_pairs
) {
    using namespace genesis::utils;

    // Get basic options.
    auto const& sample_names = options.variant_input.sample_names();
    auto const sep_char = options.table_output.get_separator_char();

    // For whole genome FST, we fill a list with all pairs.
    std::shared_ptr<BaseOutputTarget> fst_list_ofs;
    fst_list_ofs = options.file_output.get_output_target( "fst-list", "csv" );
    (*fst_list_ofs) << "first" << sep_char << "second" << sep_char << "FST\n";
    assert( sample_pairs.size() == window_fst.size() );

    // Write all pairs.
    for( size_t i = 0; i < sample_pairs.size(); ++i ) {
        auto const& pair = sample_pairs[i];
        (*fst_list_ofs) << sample_names[pair.first] << sep_char << sample_names[pair.second];
        (*fst_list_ofs) << sep_char << window_fst[i] << "\n";
    }
}

// -------------------------------------------------------------------------
//     print_fst_matrix_
// -------------------------------------------------------------------------

void print_fst_matrix_(
    FstOptions const& options,
    std::vector<double> const& window_fst,
    std::vector<std::pair<size_t, size_t>> const& sample_pairs
) {
    using namespace genesis::utils;

    // Get basic options.
    auto const& sample_names = options.variant_input.sample_names();
    auto const sep_char = options.table_output.get_separator_char();

    // Make a matrix with all the entries, so that we can easily print them.
    auto mat = Matrix<double>( sample_names.size(), sample_names.size(), 0.0 );
    for( size_t i = 0; i < sample_pairs.size(); ++i ) {
        auto const& pair = sample_pairs[i];
        mat( pair.first, pair.second ) = window_fst[i];
        mat( pair.second, pair.first ) = window_fst[i];
    }

    // For whole genome all-to-all FST, we fill a matrix with all samples.
    std::shared_ptr<BaseOutputTarget> fst_matrix_ofs;
    fst_matrix_ofs = options.file_output.get_output_target( "fst-matrix", "csv" );
    auto writer = MatrixWriter<double>( std::string( 1, sep_char ));
    writer.write( mat, fst_matrix_ofs, sample_names, sample_names, "sample" );
}

// =================================================================================================
//      Run
// =================================================================================================

void run_fst( FstOptions const& options )
{
    using namespace genesis::population;
    using namespace genesis::utils;
    // using VariantWindow = Window<genesis::population::Variant>;

    // Basic output file check
    options.file_output.check_output_files_nonexistence( "fst", "csv" );

    // Also check the additional files, if we are going to produce them.
    bool const is_whole_genome = ( options.window.window_type() == WindowOptions::WindowType::kGenome );
    bool const is_all_to_all = ( ! *options.comparand.option) && ( ! *options.comparand_list.option );
    if( is_whole_genome ) {
        options.file_output.check_output_files_nonexistence( "fst-list", "csv" );

        // If we run an all-to-all FST computation, we also want to output a matrix.
        if( is_all_to_all ) {
            options.file_output.check_output_files_nonexistence( "fst-matrix", "csv" );
        }
    }

    // -------------------------------------------------------------------------
    //     Preparation
    // -------------------------------------------------------------------------

    // Use an enum for the method, which is faster to check in the main loop than doing
    // string comparisons all the time.
    auto const method = get_enum_map_value( fst_method_map, options.method.value );

    // Before accessing the variant input, we need to add the filters to it.
    // We use a variant filter that always filters out non-SNP positions here.
    // There might still be pairs of samples between which a position is invariant,
    // but that's okay, as that will be caught by the fst computation anyway.
    // We many use this pre-filter here for speed gains,
    // to get rid of invariants before they even get to the FST computation functions.
    options.variant_input.add_combined_filter_and_transforms(
        options.filter_numerical.make_sample_filter()
    );
    genesis::population::VariantFilter total_filter;
    total_filter.only_snps = true;
    if( method == FstMethod::kKarlsson ) {
        total_filter.only_biallelic_snps = true;
    }
    options.variant_input.add_combined_filter_and_transforms(
        options.filter_numerical.make_total_filter( total_filter )
    );

    // Get indices of all pairs of samples for which we want to compute FST.
    auto const& sample_names = options.variant_input.sample_names();
    auto const sample_pairs = get_sample_pairs_( options );
    if( sample_pairs.empty() ) {
        LOG_WARN << "No pairs of samples selected, which will produce empty output. Stopping now.";
        return;
    }
    LOG_MSG << "Computing FST between " << sample_pairs.size()
            << " pair" << ( sample_pairs.size() > 1 ? "s" : "" ) << " of samples.";

    // If we do all-to-all, we expect a triangular matrix size of entries.
    internal_check(
        ! is_all_to_all || sample_pairs.size() == triangular_size( sample_names.size() ),
        "Expeting all-to-all FST to create a trinangular number of sample pairs."
    );

    // Make the processor.
    FstPoolProcessor processor = get_fst_pool_processor_( options, method, sample_pairs );

    // -------------------------------------------------------------------------
    //     Prepare Tables, Print Headers
    // -------------------------------------------------------------------------

    // Get the separator char to use for table entries.
    auto const sep_char = options.table_output.get_separator_char();

    // Prepare output file and write header line with all pairs of samples.
    auto fst_ofs = options.file_output.get_output_target( "fst", "csv" );
    (*fst_ofs) << "chrom" << sep_char << "start" << sep_char << "end" << sep_char << "snps";
    for( auto const& pair : sample_pairs ) {
        (*fst_ofs) << sep_char << sample_names[pair.first] << "." << sample_names[pair.second];
    }
    (*fst_ofs) << "\n";

    // -------------------------------------------------------------------------
    //     Main Loop
    // -------------------------------------------------------------------------

    // Iterate the file and compute per-window FST.
    size_t use_cnt = 0;
    size_t nan_cnt = 0;

    auto window_it = options.window.get_variant_window_view_iterator( options.variant_input );
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

            // TODO add thread pool? might need to be implemented in the processor itself.
            processor.process( variant );
            ++entry_count;
        }
        auto const& window_fst = processor.get_result();
        internal_check(
            window_fst.size() == sample_pairs.size(),
            "Inconsistent size of window fst values and sample pairs."
        );

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
            ++use_cnt;

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

        // For whole genome, and then also for all-to-all, we produce special tables on top.
        if( is_whole_genome ) {
            internal_check(
                options.window.get_num_windows() == 1,
                "More than one window for whole genome"
            );
            print_fst_list_( options, window_fst, sample_pairs );
            if( is_all_to_all ) {
                print_fst_matrix_( options, window_fst, sample_pairs );
            }
        }
    }
    internal_check(
        use_cnt + nan_cnt == options.window.get_num_windows(),
        "Number of windows is inconsistent."
    );

    // Final user output.
    options.window.print_report();
    if( nan_cnt > 0 ) {
        LOG_MSG << "Thereof, skipped " << nan_cnt << " window"
                << ( nan_cnt != 1 ? "s" : "" ) << " without any FST values.";
    }
    options.filter_numerical.print_report();
}
