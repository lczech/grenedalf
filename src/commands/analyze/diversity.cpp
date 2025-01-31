/*
    grenedalf - Genome Analyses of Differential Allele Frequencies
    Copyright (C) 2020-2025 Lucas Czech

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

#include "commands/analyze/diversity.hpp"
#include "options/global.hpp"
#include "tools/cli_setup.hpp"
#include "tools/misc.hpp"

#include "genesis/population/function/diversity_pool_calculator.hpp"
#include "genesis/population/function/diversity_pool_functions.hpp"
#include "genesis/population/filter/sample_counts_filter_numerical.hpp"
#include "genesis/population/filter/sample_counts_filter.hpp"
#include "genesis/population/filter/variant_filter_numerical.hpp"
#include "genesis/population/filter/variant_filter.hpp"
#include "genesis/population/function/functions.hpp"
#include "genesis/population/window/functions.hpp"
#include "genesis/utils/text/string.hpp"

#include <cassert>

using namespace genesis::population;
using namespace genesis::utils;

// =================================================================================================
//      Setup
// =================================================================================================

void setup_diversity( CLI::App& app )
{
    // Create the options and subcommand objects.
    auto options = std::make_shared<DiversityOptions>();
    auto sub = app.add_subcommand(
        "diversity",
        "Compute pool-sequencing corrected diversity measures "
        "Theta Pi, Theta Watterson, and Tajima's D."
    );

    // -------------------------------------------------------------------------
    //     Input
    // -------------------------------------------------------------------------

    // Required input of some count/frequency file format.
    options->variant_input.add_variant_input_opts_to_app( sub );
    // options->variant_input.add_sample_name_opts_to_app( sub );
    // options->variant_input.add_region_filter_opts_to_app( sub );

    // We do not add the only-SNP filter here as a user option, as we always filter out invariant
    // (non-SNP) sites below anyway. They do not contribute to the diversity.
    options->filter_numerical.add_sample_count_filter_opts_to_app( sub );
    options->filter_numerical.add_sample_read_depth_filter_opts_to_app( sub );
    options->filter_numerical.add_total_read_depth_filter_opts_to_app( sub );
    options->filter_numerical.add_total_snp_filter_opts_to_app( sub, false, true );
    options->filter_numerical.add_total_snp_count_opts_to_app( sub );
    options->filter_numerical.add_total_freq_filter_opts_to_app( sub );

    // We always need the min count option, and it cannot be 0.
    options->filter_numerical.sample_min_count.option->required();
    options->filter_numerical.sample_min_count.option->check( CLI::PositiveNumber );

    // Also offer subsampling options for this command.
    options->transform_subsample.add_subsample_opts_to_app( sub );

    // Settings for the windowing.
    options->window.add_window_opts_to_app( sub, true );

    // -------------------------------------------------------------------------
    //     Settings
    // -------------------------------------------------------------------------

    // Diversity Processor
    options->diversity_processor.add_diversity_processor_opts_to_app(
        sub, options->variant_input.get_variant_reference_genome_options()
    );

    // Settings: Omit Empty Windows
    options->no_extra_columns.option = sub->add_flag(
        "--no-extra-columns",
        options->no_extra_columns.value,
        "Do not output the extra columns containing counts for each position and sample pair "
        "that summarize the effects of the filtering. Only the window coordinates and the fst values "
        "are printed in that case."
    );
    options->no_extra_columns.option->group( "Settings" );

    // // Settings: Omit Empty Windows
    // options->no_nan_windows.option = sub->add_flag(
    //     "--no-nan-windows",
    //     options->no_nan_windows.value,
    //     "Do not output windows where all values are n/a. This is can be relevant "
    //     "with small window sizes (or individual positions), to reduce output clutter."
    // );
    // options->no_nan_windows.option->group( "Settings" );

    // -------------------------------------------------------------------------
    //     Output
    // -------------------------------------------------------------------------

    // Add table output options.
    auto sep_opt = options->table_output.add_separator_char_opt_to_app( sub );
    auto nan_opt = options->table_output.add_na_entry_opt_to_app( sub );

    // PoPoolation output format
    options->popoolation_format.option = sub->add_flag(
        "--popoolation-format",
        options->popoolation_format.value,
        "If set, instead of writing one output table for all measures and all samples, "
        "write the results in separate files for each sample and for each measure of "
        "Theta Pi, Theta Watterson, and Tajima's D, following the format of PoPoolation."
    )->group( "Formatting" );

    // Exclude separator char option and na entry in PoPoolation compatibility mode,
    // as we have to use their values in that case.
    sep_opt->excludes( options->popoolation_format.option );
    nan_opt->excludes( options->popoolation_format.option );

    // Output
    options->file_output.add_default_output_opts_to_app( sub );
    options->file_output.add_file_compress_opt_to_app( sub );

    // Set the run function as callback to be called when this subcommand is issued.
    // Hand over the options by copy, so that their shared ptr stays alive in the lambda.
    sub->callback( grenedalf_cli_callback(
        sub,
        {
            // Citation keys as needed
            "Kofler2011-PoPoolation"
        },
        [ options ]() {
            run_diversity( *options );
        }
    ));
}

// =================================================================================================
//      Helper Functions
// =================================================================================================

// -------------------------------------------------------------------------
//     Diversity Command State
// -------------------------------------------------------------------------

/**
 * @brief Collect all output_data that we might need in one place.
 */
struct DiversityCommandState
{
    // Provided loci, if needed for window averaging
    std::shared_ptr<genesis::population::GenomeLocusSet> provided_loci;

    // Which values are actually being computed?
    bool compute_theta_pi;
    bool compute_theta_wa;
    bool compute_tajima_d;

    // Store the column separator and nan entry, so that we do not need to look them up every time.
    char sep_char;
    std::string na_entry;
    bool write_extra_columns;

    // We need all pointers, for our format and for the PoPoolation format,
    // in order to avoid code duplication in the output.
    std::shared_ptr<BaseOutputTarget> table_ofs;
    std::vector<std::shared_ptr<BaseOutputTarget>> popoolation_theta_pi_ofss;
    std::vector<std::shared_ptr<BaseOutputTarget>> popoolation_theta_wa_ofss;
    std::vector<std::shared_ptr<BaseOutputTarget>> popoolation_tajima_d_ofss;
};

// -------------------------------------------------------------------------
//     prepare_output_data_
// -------------------------------------------------------------------------

/**
 * @brief Prepare output file and write fixed header fields.
 */
DiversityCommandState prepare_output_data_(
    DiversityOptions const& options,
    std::vector<std::string> const& sample_names
) {
    DiversityCommandState output;

    // Provided loci, if needed for window averaging
    output.provided_loci = options.diversity_processor.get_provided_loci();

    // Get the measures to compute. At least one of them will be active.
    output.compute_theta_pi = options.diversity_processor.compute_theta_pi();
    output.compute_theta_wa = options.diversity_processor.compute_theta_watterson();
    output.compute_tajima_d = options.diversity_processor.compute_tajima_d();

    // Get the separator char to use for table entries.
    output.sep_char = options.table_output.get_separator_char();
    output.na_entry = options.table_output.get_na_entry();
    output.write_extra_columns = !options.no_extra_columns.value;

    // -----------------------------------
    //     Settings checks
    // -----------------------------------

    if( output.compute_tajima_d ) {
        if( options.filter_numerical.sample_min_count.value != 2 ) {
            throw CLI::ValidationError(
                options.filter_numerical.sample_min_count.option->get_name(),
                "The pool-seq corrected computation of Tajima's D requires the minimum allele count "
                "to be exactly 2, i.e., exclusing singletons, according to Kofler et al. "
                "In case 2 is insufficient, we recommend to subsample the reads to a smaller "
                "read depth. Alternatively, deactivate the compuation of Tajima's D."
            );
        }
        if( options.filter_numerical.sample_min_read_depth.value == 0 ) {
            throw CLI::ValidationError(
                options.filter_numerical.sample_min_read_depth.option->get_name(),
                "The pool-seq corrected computation of Tajima's D "
                "requires the minimum read depth to be set to a value greater than 0. "
                "Alternatively, deactivate the compuation of Tajima's D."
            );
        }
    }

    // -----------------------------------
    //     Output file checks
    // -----------------------------------

    if( options.popoolation_format.value ) {
        for( auto const& sample_name : sample_names ) {
            if( output.compute_theta_pi ) {
                options.file_output.check_output_files_nonexistence(
                    "diversity-" + sample_name + "-theta-pi", "csv"
                );
            }
            if( output.compute_theta_wa ) {
                options.file_output.check_output_files_nonexistence(
                    "diversity-" + sample_name + "-theta-watterson", "csv"
                );
            }
            if( output.compute_tajima_d ) {
                options.file_output.check_output_files_nonexistence(
                    "diversity-" + sample_name + "-tajimas-d", "csv"
                );
            }
        }
    } else {
        options.file_output.check_output_files_nonexistence( "diversity", "csv" );
    }

    // -----------------------------------
    //     PoPoolation format
    // -----------------------------------

    if( options.popoolation_format.value ) {

        // Open the three PoPoolation file formats for each sample.
        for( auto const& sample_name : sample_names ) {
            if( output.compute_theta_pi ) {
                output.popoolation_theta_pi_ofss.emplace_back(
                    options.file_output.get_output_target(
                        "diversity-" + sample_name + "-theta-pi", "csv"
                    )
                );
            }
            if( output.compute_theta_wa ) {
                output.popoolation_theta_wa_ofss.emplace_back(
                    options.file_output.get_output_target(
                        "diversity-" + sample_name + "-theta-watterson", "csv"
                    )
                );
            }
            if( output.compute_tajima_d ) {
                output.popoolation_tajima_d_ofss.emplace_back(
                    options.file_output.get_output_target(
                        "diversity-" + sample_name + "-tajimas-d", "csv"
                    )
                );
            }
        }

    // -----------------------------------
    //     Table Format
    // -----------------------------------

    } else {

        // Open and prepare our table format
        output.table_ofs = options.file_output.get_output_target( "diversity", "csv" );
        (*output.table_ofs) << "chrom" << output.sep_char << "start" << output.sep_char << "end";

        // Print the extra column headers if needed.
        if( output.write_extra_columns ) {
            (*output.table_ofs) << output.sep_char << "total.masked";
            (*output.table_ofs) << output.sep_char << "total.missing";
            (*output.table_ofs) << output.sep_char << "total.empty";
            (*output.table_ofs) << output.sep_char << "total.numeric";
            (*output.table_ofs) << output.sep_char << "total.invariant";
            (*output.table_ofs) << output.sep_char << "total.passed";
        }

        // Make the header per-sample fields.
        std::vector<std::string> fields;
        if( output.write_extra_columns ) {
            fields.push_back( "missing" );
            fields.push_back( "numeric" );
            fields.push_back( "passed" );
        }
        if( output.compute_theta_pi ) {
            fields.push_back( "theta_pi" );
        }
        if( output.compute_theta_wa ) {
            fields.push_back( "theta_watterson" );
        }
        if( output.compute_tajima_d ) {
            fields.push_back( "tajimas_d" );
        }

        // Write all fields for all samples.
        for( auto const& sample : sample_names ) {
            for( auto const& field : fields ) {
                (*output.table_ofs) << output.sep_char << sample << "." << field;
            }
        }
        (*output.table_ofs) << "\n";
    }

    return output;
}

// -------------------------------------------------------------------------
//     write_output_popoolation_line_
// -------------------------------------------------------------------------

/**
 * @brief Write output for PoPoolation style tables
 */
void write_output_popoolation_line_(
    VariantWindowView const& window,
    genesis::population::DiversityPoolProcessor const& processor,
    DiversityCommandState& output
) {
    using namespace genesis::population;

    // The popoolation table formats use the same format, so let's make one function to rule them all!
    // We take the DiversityPoolCalculator::Result here for the general values, and then the actual
    // data value again, so that we don't have to switch to get it.
    // Format: "2R	19500	0	0.000	na" or "A	1500	101	1.000	1.920886709" for example.
    auto write_popoolation_line_ = [](
        std::shared_ptr<genesis::utils::BaseOutputTarget>& ofs,
        VariantWindowView const& window,
        VariantFilterCategoryStats const& variant_stats,
        SampleCountsFilterCategoryStats const& sample_stats,
        double value
    ) {
        // TODO Check if that is actually the correct way to do it.
        // This is a bit messy, as our approach of handling sample vs total filters differs
        // from PoPoolation, but this is as close as we can get, I think.
        auto const read_depth = static_cast<double>(
            variant_stats[VariantFilterTagCategory::kInvariant] +
            sample_stats[SampleCountsFilterTagCategory::kPassed]
        );
        auto const window_width = static_cast<double>( window.width() );
        auto const read_depth_fraction = read_depth / window_width;

        // Write fixed columns.
        (*ofs) << window.chromosome();
        (*ofs) << "\t" << anchor_position( window, WindowAnchorType::kIntervalMidpoint );
        (*ofs) << "\t" << sample_stats[SampleCountsFilterTagCategory::kPassed];
        (*ofs) << "\t" << std::fixed << std::setprecision( 3 ) << read_depth_fraction;
        if( std::isfinite( value ) ) {
            (*ofs) << "\t" << std::fixed << std::setprecision( 9 ) << value;
        } else {
            (*ofs) << "\tna";
        }
        (*ofs) << "\n";
    };

    // Get the results for the whole set of calculators.
    auto const results = processor.get_result( window, output.provided_loci );
    auto const total_stats = variant_filter_stats_category_counts(
        processor.get_filter_stats()
    );

    // Write to all individual files for each sample and each value.
    for( size_t i = 0; i < processor.calculators().size(); ++i ) {
        auto const& calculator = processor.calculators()[i];
        auto const sample_stats = sample_counts_filter_stats_category_counts(
            calculator.get_filter_stats()
        );

        // Theta Pi
        if( output.compute_theta_pi ) {
            write_popoolation_line_(
                output.popoolation_theta_pi_ofss[i],
                window,
                total_stats,
                sample_stats,
                results[i].theta_pi
            );
        }

        // Theta Watterson
        if( output.compute_theta_wa ) {
            write_popoolation_line_(
                output.popoolation_theta_wa_ofss[i],
                window,
                total_stats,
                sample_stats,
                results[i].theta_watterson
            );
        }

        // Tajima's D
        if( output.compute_tajima_d ) {
            write_popoolation_line_(
                output.popoolation_tajima_d_ofss[i],
                window,
                total_stats,
                sample_stats,
                results[i].tajima_d
            );
        }
    }
}

// -------------------------------------------------------------------------
//     write_output_table_line_
// -------------------------------------------------------------------------

/**
 * @brief Write output for our table format
 */
void write_output_table_line_(
    VariantWindowView const& window,
    genesis::population::DiversityPoolProcessor const& processor,
    DiversityCommandState& output
) {
    // Shorthand
    using namespace genesis::population;
    auto& table_ofs = output.table_ofs;
    auto const sep_char = output.sep_char;

    // Helper function to write a field value to one of the tables.
    auto write_table_field_ = [&]( double value ){
        if( std::isfinite( value ) ) {
            (*table_ofs) << sep_char;
            (*table_ofs) << std::defaultfloat << std::setprecision( 9 ) << value;
        } else {
            (*table_ofs) << sep_char << output.na_entry;
        }
    };

    // Write fixed columns.
    (*table_ofs) << window.chromosome();
    (*table_ofs) << sep_char << window.first_position();
    (*table_ofs) << sep_char << window.last_position();
    if( output.write_extra_columns ) {
        // Summarize the exact statistics into their category summaries.
        // Good enough for the user output here. We do not need to know which exact filter
        // failed each time; it is good enough that it was a numerical one etc.
        auto const total_stats = variant_filter_stats_category_counts(
            processor.get_filter_stats()
        );
        // Print the counters.
        // This obviously requires the column header line to exactly match this order.
        // Change both if needed, throughout this file.
        (*table_ofs) << sep_char << total_stats[ VariantFilterTagCategory::kMasked ];
        (*table_ofs) << sep_char << total_stats[ VariantFilterTagCategory::kMissingInvalid ];
        (*table_ofs) << sep_char << total_stats[ VariantFilterTagCategory::kSamplesFailed ];
        (*table_ofs) << sep_char << total_stats[ VariantFilterTagCategory::kNumeric ];
        (*table_ofs) << sep_char << total_stats[ VariantFilterTagCategory::kInvariant ];
        (*table_ofs) << sep_char << total_stats[ VariantFilterTagCategory::kPassed ];
    }

    // Get the results for the whole set of calculators.
    auto const results = processor.get_result( window, output.provided_loci );

    // Write to all individual files for each sample and each value.
    for( size_t i = 0; i < processor.calculators().size(); ++i ) {
    }

    // Write the per-pair diversity values in the correct order.
    for( size_t i = 0; i < processor.calculators().size(); ++i ) {
        auto const& calculator = processor.calculators()[i];
        auto const sample_stats = sample_counts_filter_stats_category_counts(
            calculator.get_filter_stats()
        );

        // Write all columns. The extra ones need to match the ones written in the header,
        // so if either needs updating, both need to be changed.
        if( output.write_extra_columns ) {
            (*table_ofs) << sep_char << sample_stats[ SampleCountsFilterTagCategory::kMissingInvalid ];
            (*table_ofs) << sep_char << sample_stats[ SampleCountsFilterTagCategory::kNumeric ];
            (*table_ofs) << sep_char << sample_stats[ SampleCountsFilterTagCategory::kPassed ];
        }
        if( output.compute_theta_pi ) {
            write_table_field_( results[i].theta_pi );
        }
        if( output.compute_theta_wa ) {
            write_table_field_( results[i].theta_watterson );
        }
        if( output.compute_tajima_d ) {
            write_table_field_( results[i].tajima_d );
        }
    }
    (*table_ofs) << "\n";
}

// -------------------------------------------------------------------------
//     write_output_line_
// -------------------------------------------------------------------------

/**
 * @brief Write output to output
 */
void write_output_line_(
    DiversityOptions const& options,
    VariantWindowView const& window,
    genesis::population::DiversityPoolProcessor const& processor,
    DiversityCommandState& output
) {
    // Write the data, depending on the format.
    if( options.popoolation_format.value ) {
        write_output_popoolation_line_(
            window, processor, output
        );
    } else {
        write_output_table_line_(
            window, processor, output
        );
    }
}

// =================================================================================================
//      Run
// =================================================================================================

void run_diversity( DiversityOptions const& options )
{
    // -------------------------------------------------------------------------
    //     Options Preparation
    // -------------------------------------------------------------------------

    // Before accessing the variant input, we need to add the filters to it.
    // We use a variant filter that always filters out non-SNP positions here.
    // There might still be samples where a position is invariant,
    // but that's okay, as that will be caught by the diversity computation anyway.
    options.variant_input.add_combined_filter_and_transforms(
        options.filter_numerical.make_sample_filter()
    );
    genesis::population::VariantFilterNumericalParams total_filter;
    total_filter.only_snps = true;
    options.variant_input.add_combined_filter_and_transforms(
        options.filter_numerical.make_total_filter( total_filter )
    );

    // TODO Previously, we used SampleCountsFilterNumericalParams::only_snps here instead of the
    // only_snps filter for the total Variant. This is what we documented now, and the old
    // behaviour can be achieved by simply running each sample indivudally. But maybe there is a
    // better way, where this can be made a user option instead?

    // Lastly, apply the subsampling. It is important that this happens after the above numercial
    // filters above, as otherwise we might subsample to a lower read depth, and then want to apply
    // a read depth filter, which would not work any more.
    // However, in this case, we might subsample to lower read counts that are below the threshold
    // of, e.g., the Tajuma D minimum count. So, we then apply the filter again, so that the lower
    // count filters get applied again. This is super ugly, but should work.
    // See also https://github.com/lczech/grenedalf/issues/38
    if( options.transform_subsample.add_subsample_transformation( options.variant_input )) {
        options.variant_input.add_combined_filter_and_transforms(
            options.filter_numerical.make_sample_filter()
        );
        options.variant_input.add_combined_filter_and_transforms(
            options.filter_numerical.make_total_filter( total_filter )
        );
    }

    // Get all samples names from the input file, and make a processor for them.
    auto const& sample_names = options.variant_input.sample_names();
    auto processor = options.diversity_processor.get_diversity_pool_processor(
        sample_names, options.filter_numerical.get_sample_filter_params().first
    );

    // Prepare file output_data and print headers.
    auto output_data = prepare_output_data_( options, sample_names );

    // -------------------------------------------------------------------------
    //     Main Loop
    // -------------------------------------------------------------------------

    // Iterate the file and compute per-window diversity measures.
    // We run the samples in parallel, storing their results before writing to the output file.
    auto window_it = options.window.get_variant_window_view_stream( options.variant_input );
    for( auto cur_it = window_it->begin(); cur_it != window_it->end(); ++cur_it ) {
        auto const& window = *cur_it;

        // Reset all calculator accumulators to zero for this window.
        processor.reset();

        // Compute diversity over samples.
        size_t processed_cnt = 0;
        for( auto const& variant : window ) {
            internal_check(
                variant.samples.size() == sample_names.size(),
                "Inconsistent number of samples in input file."
            );
            processor.process( variant );
            ++processed_cnt;
        }
        internal_check(
            processed_cnt == processor.get_filter_stats().sum(),
            "Number of processed positions does not equal filter sum"
        );

        // Write the output to files.
        write_output_line_( options, window, processor, output_data );
    }

    // Final user output.
    options.window.print_report();
    options.variant_input.print_report();
    // options.filter_numerical.print_report();
}
