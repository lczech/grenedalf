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

#include "commands/analyze/fst_cathedral.hpp"
#include "options/global.hpp"
#include "tools/cli_setup.hpp"
#include "tools/misc.hpp"

#include "genesis/population/function/fst_cathedral.hpp"
#include "genesis/population/function/fst_pool_calculator.hpp"
#include "genesis/population/function/fst_pool_functions.hpp"
#include "genesis/population/function/fst_pool_karlsson.hpp"
#include "genesis/population/function/fst_pool_kofler.hpp"
#include "genesis/population/function/fst_pool_processor.hpp"
#include "genesis/population/function/fst_pool_unbiased.hpp"
#include "genesis/population/function/functions.hpp"
#include "genesis/population/plotting/cathedral_plot.hpp"
#include "genesis/utils/containers/matrix.hpp"
#include "genesis/utils/containers/matrix/operators.hpp"
#include "genesis/utils/containers/matrix/writer.hpp"
#include "genesis/utils/containers/transform_iterator.hpp"
#include "genesis/utils/core/algorithm.hpp"
#include "genesis/utils/core/fs.hpp"
#include "genesis/utils/core/logging.hpp"
#include "genesis/utils/core/std.hpp"
#include "genesis/utils/core/thread_functions.hpp"
#include "genesis/utils/core/thread_pool.hpp"
#include "genesis/utils/formats/json/document.hpp"
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

std::vector<std::pair<std::string, genesis::population::CathedralWindowWidthMethod>> const window_width_method_map = {
    { "exponential", genesis::population::CathedralWindowWidthMethod::kExponential },
    { "geometric",   genesis::population::CathedralWindowWidthMethod::kGeometric },
    { "linear",      genesis::population::CathedralWindowWidthMethod::kLinear },
};

// =================================================================================================
//      Setup
// =================================================================================================

void setup_fst_cathedral( CLI::App& app )
{
    // Create the options and subcommand objects.
    auto options = std::make_shared<FstCathedralOptions>();
    auto sub = app.add_subcommand(
        "fst-cathedral",
        "Compute the data for an FST cathedral plot."
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
    options->filter_numerical.add_total_snp_count_opts_to_app( sub );
    options->filter_numerical.add_total_freq_filter_opts_to_app( sub );

    // We amend the existing biallelic SNP filter option description here.
    auto& tobs = options->filter_numerical.total_only_biallelic_snps.option;
    tobs->description(
        tobs->get_description() + "\nBy default, we already filter out invariant sites, so that "
        "only SNPs remain. With this option however, this is further reduced to only biallelic "
        "(pure) SNPs. Note that with `--method karlsson`, we implicitly also activate to filter "
        "for biallelic SNPs only, as the method requires it."
    );

    // -------------------------------------------------------------------------
    //     Settings
    // -------------------------------------------------------------------------

    // Fst Processor
    options->fst_processor.add_fst_processor_opts_to_app( sub, false, "Settings" );

    // Width
    options->cathedral_width.option = sub->add_option(
        "--cathedral-width",
        options->cathedral_width.value,
        "Width of the plot, in pixels. In particular, this determines the resolution of the last "
        "row of the plot, where each pixel corresponds to a window of the size "
        "`genome length / plot width`."
    );
    options->cathedral_width.option->group( "Settings" );

    // Height
    options->cathedral_height.option = sub->add_option(
        "--cathedral-height",
        options->cathedral_height.value,
        "Height of the plot, in pixels. This determines the number of different window sizes "
        "displayed in the plot, from whole chromosome at the top to finest resolution at the bottom."
    );
    options->cathedral_height.option->group( "Settings" );

    // Method for scaling. Currently hidden, as exponential just makes the most sense.
    options->cathedral_method.option = sub->add_option(
        "--cathedral-scaling",
        options->cathedral_method.value,
        "Scaling of the window sizes from top to bottom of the plot. We recommend `exponential`, "
        "as this covers the orders of magnitude between the two extrems in the most intiutive way."
    );
    options->cathedral_method.option->group( "" );
    options->cathedral_method.option->transform(
        CLI::IsMember( enum_map_keys( window_width_method_map ), CLI::ignore_case )
    );

    // -------------------------------------------------------------------------
    //     Output
    // -------------------------------------------------------------------------

    // Output
    options->file_output.add_default_output_opts_to_app( sub );
    // options->file_output.add_file_compress_opt_to_app( sub );

    // Set the run function as callback to be called when this subcommand is issued.
    // Hand over the options by copy, so that their shared ptr stays alive in the lambda.
    sub->callback( grenedalf_cli_callback(
        sub,
        {
            // Citation keys as needed
        },
        [ options ]() {
            run_fst_cathedral( *options );
        }
    ));
}

// =================================================================================================
//      Run
// =================================================================================================

void run_fst_cathedral( FstCathedralOptions const& options )
{
    using namespace genesis::population;
    using namespace genesis::utils;

    // Check that none of the output files exist.
    options.file_output.check_output_files_nonexistence( "cathedral-plot-*", "csv" );
    options.file_output.check_output_files_nonexistence( "cathedral-plot-*", "json" );

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
    genesis::population::VariantFilterNumericalParams total_filter;
    total_filter.only_snps = true;
    // if( state.method == FstProcessorOptions::FstMethod::kKarlsson ) {
    //     total_filter.only_biallelic_snps = true;
    // }
    options.variant_input.add_combined_filter_and_transforms(
        options.filter_numerical.make_total_filter( total_filter )
    );

    // Get indices of all pairs of samples for which we want to compute FST.
    auto const& sample_names = options.variant_input.sample_names();
    auto const sample_pairs = options.fst_processor.get_sample_pairs( sample_names );
    if( sample_pairs.empty() ) {
        LOG_MSG << "No sample pairs selected. Nothing to do.";
        return;
    }

    // Set the basic cathedral plot params.
    CathedralPlotParameters params;
    params.width = options.cathedral_width.value;
    params.height = options.cathedral_height.value;
    params.window_width_method = get_enum_map_value(
        window_width_method_map, options.cathedral_method.value
    );

    // Make the processor.
    FstPoolProcessor processor = options.fst_processor.get_fst_pool_processor(
        sample_names, sample_pairs
    );

    // Translate from our internal enum to the one used by the unbiased estimator.
    // Not super elegant, but good enough for now.
    FstPoolCalculatorUnbiased::Estimator fst_estimator;
    switch( options.fst_processor.get_fst_method() ) {
        case FstProcessorOptions::FstMethod::kUnbiasedNei: {
            fst_estimator = FstPoolCalculatorUnbiased::Estimator::kNei;
            break;
        }
        case FstProcessorOptions::FstMethod::kUnbiasedHudson: {
            fst_estimator = FstPoolCalculatorUnbiased::Estimator::kHudson;
            break;
        }
        default: {
            throw CLI::ValidationError(
                "Invalid FST estimator, only `nei` and `hudson` can be used for cathedral plots."
            );
        }
    }

    // -------------------------------------------------------------------------
    //     Main Loop
    // -------------------------------------------------------------------------

    // We stream through the input, and compute the matrixes at each chromosome while we go.
    // That way, we do not have to keep the records in memory for more than one chromosome
    // at a time, as opposed to the all-in-one function commented out below.
    auto iterator = options.variant_input.get_stream().begin();
    while( iterator ) {
        // Get SNP data for the chromosome, for all pairs of samples.
        // The function stops at the end of the chromosome or end of the input.
        auto records = compute_fst_cathedral_records_for_chromosome(
            iterator, processor, fst_estimator, sample_names, options.variant_input.get_reference_dict()
        );
        assert( records.size() == processor.size() );

        // Compute and store all of them in files, for this chromosome.
        // This is the expensive part of the computation, which calculates the value for each
        // pixel of the plot. Hence, we parallelize this loop over sample pairs (the records).
        genesis::utils::parallel_for_each(
            records.begin(), records.end(),
            [&]( FstCathedralPlotRecord& record ){
                // Compute the matrix
                compute_fst_cathedral_matrix( params, record );

                // Get a consistent base name for the output.
                // At the moment, this is manually matched to the output of the fst-cathedral
                // command, so that all related output files share a common base name.
                // Maybe there is a better way to put this information somewhere, for consistency.
                auto const& plt_name = record.plot_name;
                auto const& chr_name = record.chromosome_name;
                auto const base_name = "cathedral-plot-" + plt_name + "-" + chr_name;

                // Store the record meta data and the matrix in files
                auto const document = fst_cathedral_plot_record_to_json_document( record );
                save_cathedral_plot_record_to_targets(
                    document,
                    record.value_matrix,
                    options.file_output.get_output_target( base_name, "json" ),
                    options.file_output.get_output_target( base_name, "csv" )
                );
            }
        );
    }
    if( iterator != options.variant_input.get_stream().end() ) {
        throw std::runtime_error( "Internal error: Iterator not at end()" );
    }

    // -------------------------------------------------------------------------
    //     Main Computation
    // -------------------------------------------------------------------------

    // Alternative version of the above, that uses the pre-coded all-chromosome function in
    // genesis, instead of looping ourselves here. It hence keeps the records for all pairs of
    // samples for _all_ chromosomes in memory, which can be a lot, and hence we do not want
    // to use it here. Still, keeping it for reference below, as an example of that function.

    // // Compute all records, for all sample pairs, and all chromosomes.
    // auto records = compute_fst_cathedral_records(
    //     options.variant_input.get_stream(), processor, fst_estimator, sample_names, options.variant_input.get_reference_dict()
    // );
    //
    // // Compute and store all of them in files
    // for( auto& record : records ) {
    //     // Compute the matrix
    //     compute_fst_cathedral_matrix( params, record );
    //
    //     // Store the record meta data and the matrix in files
    //     auto const base_name = "cathedral-plot-" + record.plot_name + "-" + record.chromosome_name;
    //     auto const document = fst_cathedral_plot_record_to_json_document( record );
    //     save_cathedral_plot_record_to_targets(
    //         document,
    //         record.value_matrix,
    //         options.file_output.get_output_target( base_name, "json" ),
    //         options.file_output.get_output_target( base_name, "csv" )
    //     );
    // }
}
