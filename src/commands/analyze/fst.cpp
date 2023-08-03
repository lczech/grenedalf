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

#include "commands/analyze/fst.hpp"
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
#include <tuple>
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
        "Compute pool-sequencing corrected measures of FST."
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

    // Settings: Write pi tables
    options->write_pi_tables.option = sub->add_flag(
        "--write-pi-tables",
        options->write_pi_tables.value,
        "When using either of the two unbiased (Nei or Hudson) estimators, also write tables "
        "of the involved pi values (within, between, and total). See our equations for details. "
        "We here scale the values by the number of variant sites (which is printed in the output "
        "table column `snps`); multiplying the pi and the snps per row gives the absolute sum of pi "
        " values, which can be used if a different scaling is required instead."
    );
    options->write_pi_tables.option->group( "Settings" );

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

    // Setting: threading_threshold
    options->threading_threshold.option = sub->add_option(
        "--threading-threshold",
        options->threading_threshold.value,
        "When only computing FST beween a few pairs of samples, doing so in individual threads "
        "usually incurs a substantial overhead due to thread synchronozation. Hence, we only want "
        "to use threads for the FST computation if many pairs of samples are being computed. "
        "This setting determiens the number of sample pairs at which threads are used for FST. "
        "(Note that we still always use threads for input file parsing.)"
    );
    options->threading_threshold.option->group( "" );

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
//      Setup Functions
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
        : std::vector<size_t>( sample_names.size(), 0 )
    );

    // For our unbiased and for Kofler, check that we got the right number of pool sizes.
    internal_check(
        ! needs_pool_sizes || pool_sizes.size() == sample_names.size(),
        "Inconsistent number of samples and number of pool sizes."
    );

    // Make the type of processor that we need for the provided method.
    FstPoolProcessor processor;
    switch( method ) {
        case FstMethod::kUnbiasedNei: {
            processor = make_fst_pool_processor<FstPoolCalculatorUnbiased>(
                sample_pairs, pool_sizes, FstPoolCalculatorUnbiased::Estimator::kNei
            );
            break;
        }
        case FstMethod::kUnbiasedHudson: {
            processor = make_fst_pool_processor<FstPoolCalculatorUnbiased>(
                sample_pairs, pool_sizes, FstPoolCalculatorUnbiased::Estimator::kHudson
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
    processor.threading_threshold( options.threading_threshold.value );
    return processor;
}

// =================================================================================================
//      Output Functions
// =================================================================================================

/**
 * @brief State of the command.
 *
 * This is a simple POD to keep all data related to the fst computation and output files
 * for per-position/window files in one place, so that we can break stuff into functions.
 */
struct FstCommandState
{
    FstMethod method;

    // The output to write to, for per-window output.
    // The pi ones are only used when using fst unbiased nei or hudson,
    // and when the respective option is set.
    std::shared_ptr<::genesis::utils::BaseOutputTarget> fst_output_target;
    std::shared_ptr<::genesis::utils::BaseOutputTarget> pi_within_output_target;
    std::shared_ptr<::genesis::utils::BaseOutputTarget> pi_between_output_target;
    std::shared_ptr<::genesis::utils::BaseOutputTarget> pi_total_output_target;

    // How many windows were actually used, and how many skipped due to only containing nan.
    size_t used_count = 0;
    size_t skipped_count = 0;

    // Are we using whole genome mode, and all to all?
    // In those cases, we need to produce some extra files.
    bool is_whole_genome;
    bool is_all_to_all;
};

/**
 * @brief Collection of vectors for pi within, between, and total, in that order.
 */
using PiVectorTuple = std::tuple<std::vector<double>, std::vector<double>, std::vector<double>>;

// -------------------------------------------------------------------------
//     check_output_files_
// -------------------------------------------------------------------------

void check_output_files_(
    FstOptions const& options,
    FstCommandState& state
) {
    // Check that the pi tables are actually computed.
    if(
        options.write_pi_tables.value && (
            state.method != FstMethod::kUnbiasedNei &&
            state.method != FstMethod::kUnbiasedHudson
        )
    ) {
        throw CLI::ValidationError(
            options.write_pi_tables.option->get_name(),
            "Can only write pi tables when the unbiased Nei or Hudson estimator is being used."
        );
    }

    // Basic output file check
    options.file_output.check_output_files_nonexistence( "fst", "csv" );
    if( options.write_pi_tables.value ) {
        options.file_output.check_output_files_nonexistence( "pi-within", "csv" );
        options.file_output.check_output_files_nonexistence( "pi-between", "csv" );
        options.file_output.check_output_files_nonexistence( "pi-total", "csv" );
    }

    // Also check the additional files, if we are going to produce them.
    if( state.is_whole_genome ) {
        options.file_output.check_output_files_nonexistence( "fst-list", "csv" );
        if( options.write_pi_tables.value ) {
            options.file_output.check_output_files_nonexistence( "pi-within-list", "csv" );
            options.file_output.check_output_files_nonexistence( "pi-between-list", "csv" );
            options.file_output.check_output_files_nonexistence( "pi-total-list", "csv" );
        }

        // If we run an all-to-all FST computation, we also want to output a matrix.
        if( state.is_all_to_all ) {
            options.file_output.check_output_files_nonexistence( "fst-matrix", "csv" );
            if( options.write_pi_tables.value ) {
                options.file_output.check_output_files_nonexistence( "pi-within-matrix", "csv" );
                options.file_output.check_output_files_nonexistence( "pi-between-matrix", "csv" );
                options.file_output.check_output_files_nonexistence( "pi-total-matrix", "csv" );
            }
        }
    }
}

// -------------------------------------------------------------------------
//     prepare_output_files_
// -------------------------------------------------------------------------

void prepare_output_files_(
    FstOptions const& options,
    std::vector<std::pair<size_t, size_t>> const& sample_pairs,
    FstCommandState& state
) {
    // Get basic options.
    auto const& sample_names = options.variant_input.sample_names();
    auto const sep_char = options.table_output.get_separator_char();

    // Helper function to write the header of our per-window file formats
    auto open_file_and_write_header_ = [&](
        std::shared_ptr<::genesis::utils::BaseOutputTarget>& target,
        std::string const& base_name
    ) {
        // Make an output target, i.e., get the file name, and open the file
        target = options.file_output.get_output_target( base_name, "csv" );

        // Write the header line to it.
        (*target) << "chrom" << sep_char << "start" << sep_char << "end" << sep_char << "snps";
        for( auto const& pair : sample_pairs ) {
            (*target) << sep_char << sample_names[pair.first] << ":" << sample_names[pair.second];
        }
        (*target) << "\n";
    };

    // Prepare output file and write header line with all pairs of samples.
    open_file_and_write_header_( state.fst_output_target, "fst" );

    // Also create the per-window output files for the pi files, if needed.
    if( options.write_pi_tables.value ) {
        open_file_and_write_header_( state.pi_within_output_target, "pi-within" );
        open_file_and_write_header_( state.pi_between_output_target, "pi-between" );
        open_file_and_write_header_( state.pi_total_output_target, "pi-total" );
    }
}

// -------------------------------------------------------------------------
//     get_pi_vectors_
// -------------------------------------------------------------------------

/**
 * @brief Get the pi values for all pairs of samples.
 *
 * We here currently need to dynamically cast the FstPoolCalculator instances, which is ugly.
 * Might want to refactor later, but good enough for now to get it working.
 *
 * Returns vectors for pi within, between, and total, in that order.
 */
PiVectorTuple get_pi_vectors_(
     ::genesis::population::FstPoolProcessor const& processor
) {
    using namespace ::genesis::population;

    // Prepare the result vectors.
    PiVectorTuple result;
    auto const result_size = processor.calculators().size();
    std::get<0>(result).resize( result_size );
    std::get<1>(result).resize( result_size );
    std::get<2>(result).resize( result_size );

    // Get the pi values from all calculators, assuming that they are of the correct type.
    for( size_t i = 0; i < result_size; ++i ) {
        auto& raw_calc = processor.calculators()[i];
        auto cast_calc = dynamic_cast<FstPoolCalculatorUnbiased const*>( raw_calc.get() );
        internal_check( cast_calc, "Unbiased Nei/Hudson calculator expected." );

        // We print the values, scaled by number of variants
        // (i.e., the number of values that were added up there in the first place).
        std::get<0>(result)[i] = cast_calc->get_pi_within()  / processor.processed_count();
        std::get<1>(result)[i] = cast_calc->get_pi_between() / processor.processed_count();
        std::get<2>(result)[i] = cast_calc->get_pi_total()   / processor.processed_count();
    }

    return result;
}

// -------------------------------------------------------------------------
//     print_pairwise_value_list_
// -------------------------------------------------------------------------

/**
 * @brief Print a pairwise list of values to a file.
 *
 * This is meant for whole genome computations, where we have a single value per pair of samples.
 * This function then prints these values, in rows of the form "Sample1 Sample2 value".
 */
void print_pairwise_value_list_(
    FstOptions const& options,
    std::string const& out_file_basename,
    std::vector<double> const& values,
    std::vector<std::pair<size_t, size_t>> const& sample_pairs
) {
    using namespace genesis::utils;

    // Get basic options.
    auto const& sample_names = options.variant_input.sample_names();
    auto const sep_char = options.table_output.get_separator_char();

    // For whole genome FST, we fill a list with all pairs.
    std::shared_ptr<BaseOutputTarget> output_target;
    output_target = options.file_output.get_output_target( out_file_basename + "-list", "csv" );
    (*output_target) << "first" << sep_char << "second" << sep_char << "FST\n";
    assert( sample_pairs.size() == values.size() );

    // Write all pairs.
    for( size_t i = 0; i < sample_pairs.size(); ++i ) {
        auto const& pair = sample_pairs[i];
        (*output_target) << sample_names[pair.first] << sep_char << sample_names[pair.second];
        (*output_target) << sep_char << values[i] << "\n";
    }
}

// -------------------------------------------------------------------------
//     print_pairwise_value_matrix_
// -------------------------------------------------------------------------

/**
 * @brief Print a pairwise matrix of values to a file.
 *
 * This is meant for whole genome computations, where we have a single value per pair of samples.
 * This function then prints these values, as a symmetrical matrix, with the samples in the rows
 * and columns.
 */
void print_pairwise_value_matrix_(
    FstOptions const& options,
    std::string const& out_file_basename,
    std::vector<double> const& values,
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
        mat( pair.first, pair.second ) = values[i];
        mat( pair.second, pair.first ) = values[i];
    }

    // For whole genome all-to-all FST, we fill a matrix with all samples.
    std::shared_ptr<BaseOutputTarget> output_target;
    output_target = options.file_output.get_output_target( out_file_basename + "-matrix", "csv" );
    auto writer = MatrixWriter<double>( std::string( 1, sep_char ));
    writer.write( mat, output_target, sample_names, sample_names, "sample" );
}

// -------------------------------------------------------------------------
//     print_output_line_
// -------------------------------------------------------------------------

void print_output_line_(
    FstOptions const& options,
    ::genesis::population::FstPoolProcessor const& processor,
    VariantWindowView const& window,
    std::vector<double> const& values,
    ::genesis::utils::BaseOutputTarget& target
) {
    // Get the separator char to use for table entries.
    auto const sep_char = options.table_output.get_separator_char();

    // Write fixed columns.
    target << window.chromosome();
    target << sep_char << window.first_position();
    target << sep_char << window.last_position();
    target << sep_char << processor.processed_count();

    // Write the per-pair FST values in the correct order.
    for( auto const& fst : values ) {
        if( std::isfinite( fst ) ) {
            target << sep_char << fst;
        } else {
            target << sep_char << options.table_output.get_na_entry();
        }
    }
    target << "\n";
}

// -------------------------------------------------------------------------
//     print_fst_files_
// -------------------------------------------------------------------------

/**
 * @brief Print the fst values for all sample pairs, for a given window.
 *
 * This function is called once per window, and writes out the per-pair results to the output file.
 */
void print_to_output_files_(
    FstOptions const& options,
    ::genesis::population::FstPoolProcessor const& processor,
    VariantWindowView const& window,
    std::vector<std::pair<size_t, size_t>> const& sample_pairs,
    FstCommandState& state
) {
    // Get the results and check them.
    auto const& window_fst = processor.get_result();
    internal_check(
        window_fst.size() == sample_pairs.size(),
        "Inconsistent size of window fst values and sample pairs."
    );

    // In case that we want to print the pi values for the unbiased estimator, get those.
    PiVectorTuple pi_vectors;
    if( options.write_pi_tables.value ) {
        // Assert that we indeed are computing one of the unbiased estimators.
        internal_check(
            state.method == FstMethod::kUnbiasedNei || state.method == FstMethod::kUnbiasedHudson,
            "Unbiased Nei/Hudson calculator expected when using --write-pi-values"
        );

        // Get the values to be printed here
        pi_vectors = get_pi_vectors_( processor );
        internal_check( std::get<0>(pi_vectors).size() == sample_pairs.size() );
        internal_check( std::get<1>(pi_vectors).size() == sample_pairs.size() );
        internal_check( std::get<2>(pi_vectors).size() == sample_pairs.size() );
    }

    // Write the values, unless all of them are nan, then we might skip.
    // Skip empty windows if the user wants to.
    if(
        options.omit_na_windows.value && (
            processor.processed_count() == 0 ||
            std::none_of( window_fst.begin(), window_fst.end(), []( double v ) {
                return std::isfinite( v );
            })
        )
    ) {
        ++state.skipped_count;
    } else {
        ++state.used_count;
        print_output_line_( options, processor, window, window_fst, *state.fst_output_target );
        if( options.write_pi_tables.value ) {
            print_output_line_(
                options, processor, window, std::get<0>(pi_vectors), *state.pi_within_output_target
            );
            print_output_line_(
                options, processor, window, std::get<1>(pi_vectors), *state.pi_between_output_target
            );
            print_output_line_(
                options, processor, window, std::get<2>(pi_vectors), *state.pi_total_output_target
            );
        }
    }

    // For whole genome, and then also for all-to-all, we produce special tables on top.
    if( state.is_whole_genome ) {
        internal_check( options.window.get_num_windows() == 1, "Whole genome with multiple windows" );
        print_pairwise_value_list_( options, "fst", window_fst, sample_pairs );
        if( options.write_pi_tables.value ) {
            print_pairwise_value_list_(
                options, "pi-within",  std::get<0>(pi_vectors), sample_pairs
            );
            print_pairwise_value_list_(
                options, "pi-between", std::get<1>(pi_vectors), sample_pairs
            );
            print_pairwise_value_list_(
                options, "pi-total",   std::get<2>(pi_vectors), sample_pairs
            );
        }

        if( state.is_all_to_all ) {
            print_pairwise_value_matrix_( options, "fst", window_fst, sample_pairs );
            if( options.write_pi_tables.value ) {
                print_pairwise_value_matrix_(
                    options, "pi-within",  std::get<0>(pi_vectors), sample_pairs
                );
                print_pairwise_value_matrix_(
                    options, "pi-between", std::get<1>(pi_vectors), sample_pairs
                );
                print_pairwise_value_matrix_(
                    options, "pi-total",   std::get<2>(pi_vectors), sample_pairs
                );
            }
        }
    }
}

// =================================================================================================
//      Run
// =================================================================================================

void run_fst( FstOptions const& options )
{
    using namespace genesis::population;
    using namespace genesis::utils;

    // The state POD that stores the run objects in one place.
    FstCommandState state;
    state.method = get_enum_map_value( fst_method_map, options.method.value );
    state.is_whole_genome = ( options.window.window_type() == WindowOptions::WindowType::kGenome );
    state.is_all_to_all = ( ! *options.comparand.option) && ( ! *options.comparand_list.option );

    // Check that none of the output files exist.
    check_output_files_( options, state );

    // -------------------------------------------------------------------------
    //     Options Preparation
    // -------------------------------------------------------------------------

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
    if( state.method == FstMethod::kKarlsson ) {
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
        ! state.is_all_to_all || sample_pairs.size() == triangular_size( sample_names.size() ),
        "Expeting all-to-all FST to create a trinangular number of sample pairs."
    );

    // Make the processor.
    FstPoolProcessor processor = get_fst_pool_processor_( options, state.method, sample_pairs );

    // -------------------------------------------------------------------------
    //     Prepare Tables, Print Headers
    // -------------------------------------------------------------------------

    // Prepare the main output file per window
    prepare_output_files_( options, sample_pairs, state );

    // -------------------------------------------------------------------------
    //     Main Loop
    // -------------------------------------------------------------------------

    auto window_it = options.window.get_variant_window_view_iterator( options.variant_input );
    for( auto cur_it = window_it->begin(); cur_it != window_it->end(); ++cur_it ) {
        auto const& window = *cur_it;

        // Reset all calculator accumulators to zero for this window.
        processor.reset();

        // Compute diversity over samples.
        // #pragma omp parallel for
        for( auto const& variant : window ) {
            internal_check(
                variant.samples.size() == sample_names.size(),
                "Inconsistent number of samples in input file."
            );
            processor.process( variant );
        }

        // Print the output to all files that we want.
        print_to_output_files_( options, processor, window, sample_pairs, state );
    }
    internal_check(
        state.used_count + state.skipped_count == options.window.get_num_windows(),
        "Number of windows is inconsistent."
    );

    // Final user output.
    options.window.print_report();
    if( state.skipped_count > 0 ) {
        LOG_MSG << "Thereof, skipped " << state.skipped_count << " window"
                << ( state.skipped_count != 1 ? "s" : "" ) << " without any FST values.";
    }
    options.filter_numerical.print_report();
}
