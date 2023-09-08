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

#include "commands/analyze/diversity.hpp"
#include "options/global.hpp"
#include "tools/cli_setup.hpp"
#include "tools/misc.hpp"

#include "genesis/population/functions/diversity_pool_calculator.hpp"
#include "genesis/population/functions/diversity_pool_functions.hpp"
#include "genesis/population/functions/filter_transform.hpp"
#include "genesis/population/functions/functions.hpp"
#include "genesis/population/window/functions.hpp"
#include "genesis/utils/text/string.hpp"

#include <cassert>

// Keep it simple...
using namespace genesis::population;
using namespace genesis::utils;
using VariantWindowView = WindowView<genesis::population::Variant>;

// =================================================================================================
//      Setup
// =================================================================================================

void setup_diversity( CLI::App& app )
{
    // Create the options and subcommand objects.
    auto options = std::make_shared<DiversityOptions>();
    auto sub = app.add_subcommand(
        "diversity",
        "Compute pool-sequencing corrected diversity measures Theta Pi, "
        "Theta Watterson, and Tajima's D."
    );

    // -------------------------------------------------------------------------
    //     Input
    // -------------------------------------------------------------------------

    // Required input of some count/frequency file format.
    options->variant_input.add_variant_input_opts_to_app( sub );
    // options->variant_input.add_sample_name_opts_to_app( sub );
    // options->variant_input.add_region_filter_opts_to_app( sub );

    // Individual settings for numerical filtering.
    // We do not add the SNP filters here as user options, as we want to do our own filtering later
    // on, so that we can keep track of covered invariant sites for the relative theta computation.
    // We also do not add total filters here, as diversity is computed per-sample.
    options->filter_numerical.add_sample_count_filter_opts_to_app( sub );
    options->filter_numerical.add_sample_coverage_filter_opts_to_app( sub );

    // We always need the min count option, and it cannot be 0.
    options->filter_numerical.sample_min_count.option->required();
    options->filter_numerical.sample_min_count.option->check( CLI::PositiveNumber );

    // Settings for the windowing.
    options->window.add_window_opts_to_app( sub, true );

    // -------------------------------------------------------------------------
    //     Settings
    // -------------------------------------------------------------------------

    // Settings: Pool sizes
    options->poolsizes.add_poolsizes_opt_to_app( sub );

    // Measure
    options->measure.option = sub->add_option(
        "--measure",
        options->measure.value,
        "Diversity measure to compute."
    );
    options->measure.option->group( "Settings" );
    options->measure.option->transform(
        CLI::IsMember({ "all", "theta-pi", "theta-watterson", "tajimas-d" }, CLI::ignore_case )
    );

    // Below options are deactivated, as we use the more standardized numerical filters instead now.
    // However, we currently to not apply them in the input stream, as we want to keep track
    // of invariant positions, so that for file formats that contain data for all variants,
    // we can accurately compute the relative Theta values.

    // // Minimum allele count
    // options->min_count.option = sub->add_option(
    //     "--min-allele-count",
    //     options->min_count.value,
    //     "Minimum allele count needed for a base count (`ACGT`) to be considered. Bases below this "
    //     "count are filtered out. Used for the identification of SNPs."
    // )->group( "Settings" );
    //
    // // Minimum coverage
    // options->min_coverage.option = sub->add_option(
    //     "--min-coverage",
    //     options->min_coverage.value,
    //     "Minimum coverage of a site. Sites with a lower coverage will not be considered "
    //     "for SNP identification and coverage estimation."
    // )->group( "Settings" );
    //
    // // Maximum coverage
    // options->max_coverage.option = sub->add_option(
    //     "--max-coverage",
    //     options->max_coverage.value,
    //     "Maximum coverage used for SNP identification. Coverage in ALL populations has to be lower "
    //     "or equal to this threshold, otherwise no SNP will be called."
    // )->group( "Settings" );

    // TODO reactivate, but per sample, instead of for all populations

    // Minimum coverage fraction
    // options->min_coverage_fraction.option = sub->add_option(
    //     "--min-coverage-fraction",
    //     options->min_coverage_fraction.value,
    //     "Minimum coverage fraction of a window being between `--min-coverage` and `--max-coverage` "
    //     "in ALL populations that needs to be reached in order to compute the diversity measures."
    // )->group( "Settings" );

    // -------------------------------------------------------------------------
    //     Formatting and Compatibility
    // -------------------------------------------------------------------------

    // Add table output options.
    auto sep_opt = options->table_output.add_separator_char_opt_to_app( sub );
    auto nan_opt = options->table_output.add_na_entry_opt_to_app( sub );

    // Add an option to purposely activate the PoPoolation bugs that we discovered.
    options->with_popoolation_bugs.option = sub->add_flag(
        "--popoolation-corrected-tajimas-d",
        options->with_popoolation_bugs.value,
        "If set, use the pool-sequencing correction term for Tajima's D as implemented in "
        "PoPoolation <=1.2.2, which contains two major bugs that alter numerical results. "
        "We here offer to imitate these bugs for comparability with PoPoolation."
    )->group( "Compatibility" );

    // PoPoolation output format
    options->popoolation_format.option = sub->add_flag(
        "--popoolation-format",
        options->popoolation_format.value,
        "If set, instead of writing one output table for all measures and all samples, "
        "write the results in separate files for each sample and for each measure of "
        "Theta Pi, Theta Watterson, and Tajima's D, following the format of PoPoolation."
    )->group( "Compatibility" );

    // Exclude separator char option and na entry in PoPoolation compatibility mode,
    // as we have to use their values in that case.
    sep_opt->excludes( options->popoolation_format.option );
    nan_opt->excludes( options->popoolation_format.option );

    // -------------------------------------------------------------------------
    //     Output
    // -------------------------------------------------------------------------

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
//     DiversityOutputData
// -------------------------------------------------------------------------

/**
 * @brief Collect all output_data that we might need in one place.
 */
struct DiversityOutputData
{
    // Which values are actually being computed?
    bool compute_theta_pi;
    bool compute_theta_wa;
    bool compute_tajima_d;

    // Store the column separator and nan entry, so that we do not need to look them up every time.
    char sep_char;
    std::string na_entry;

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
DiversityOutputData prepare_output_data_(
    DiversityOptions const& options,
    std::vector<std::string> const& sample_names
) {
    DiversityOutputData output_data;

    // Get the measures to compute. At least one of them will be active.
    auto const measure = to_lower( options.measure.value );
    output_data.compute_theta_pi = ( measure == "all" || measure == "theta-pi" );
    output_data.compute_theta_wa = ( measure == "all" || measure == "theta-watterson" );
    output_data.compute_tajima_d = ( measure == "all" || measure == "tajimas-d" );
    assert(
        output_data.compute_theta_pi ||
        output_data.compute_theta_wa ||
        output_data.compute_tajima_d
    );

    // Get the separator char to use for table entries.
    output_data.sep_char = options.table_output.get_separator_char();
    output_data.na_entry = options.table_output.get_na_entry();

    // -----------------------------------
    //     Settings checks
    // -----------------------------------

    if( output_data.compute_tajima_d ) {
        if( options.filter_numerical.sample_min_count.value != 2 ) {
            throw CLI::ValidationError(
                options.filter_numerical.sample_min_count.option->get_name(),
                "The pool-seq corrected computation of Tajima's D requires the minimum allele count "
                "to be exactly 2, according to Kofler et al. In case 2 is insufficient, "
                "we recommend to subsample the reads to a smaller coverage. "
                "Alternatively, deactivate the compuation of Tajima's D."
            );
        }
        if( options.filter_numerical.sample_min_coverage.value == 0 ) {
            throw CLI::ValidationError(
                options.filter_numerical.sample_min_coverage.option->get_name(),
                "The pool-seq corrected computation of Tajima's D "
                "requires the minimum coverage to be set to a value greater than 0. "
                "Alternatively, deactivate the compuation of Tajima's D."
            );
        }
    }

    // -----------------------------------
    //     Output file checks
    // -----------------------------------

    if( options.popoolation_format.value ) {
        for( auto const& sample_name : sample_names ) {
            if( output_data.compute_theta_pi ) {
                options.file_output.check_output_files_nonexistence(
                    "diversity-" + sample_name + "-theta-pi", "csv"
                );
            }
            if( output_data.compute_theta_wa ) {
                options.file_output.check_output_files_nonexistence(
                    "diversity-" + sample_name + "-theta-watterson", "csv"
                );
            }
            if( output_data.compute_tajima_d ) {
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
            if( output_data.compute_theta_pi ) {
                output_data.popoolation_theta_pi_ofss.emplace_back(
                    options.file_output.get_output_target(
                        "diversity-" + sample_name + "-theta-pi", "csv"
                    )
                );
            }
            if( output_data.compute_theta_wa ) {
                output_data.popoolation_theta_wa_ofss.emplace_back(
                    options.file_output.get_output_target(
                        "diversity-" + sample_name + "-theta-watterson", "csv"
                    )
                );
            }
            if( output_data.compute_tajima_d ) {
                output_data.popoolation_tajima_d_ofss.emplace_back(
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
        output_data.table_ofs = options.file_output.get_output_target( "diversity", "csv" );
        (*output_data.table_ofs) << "chrom" << output_data.sep_char << "start" << output_data.sep_char << "end";

        // Make the header per-sample fields.
        std::vector<std::string> fields;

        // fields.push_back( "variant_count" );
        // fields.push_back( "coverage_count" );
        fields.push_back( "snp_count" );
        fields.push_back( "coverage_fraction" );

        if( output_data.compute_theta_pi ) {
            fields.push_back( "theta_pi_abs" );
            fields.push_back( "theta_pi_rel" );
        }
        if( output_data.compute_theta_wa ) {
            fields.push_back( "theta_watterson_abs" );
            fields.push_back( "theta_watterson_rel" );
        }
        if( output_data.compute_tajima_d ) {
            fields.push_back( "tajimas_d" );
        }

        // Write all fields for all samples.
        for( auto const& sample : sample_names ) {
            for( auto const& field : fields ) {
                (*output_data.table_ofs) << output_data.sep_char << sample << "." << field;
            }
        }
        (*output_data.table_ofs) << "\n";
    }

    return output_data;
}

// -------------------------------------------------------------------------
//     get_diversity_calculators_
// -------------------------------------------------------------------------

/**
 * @brief Get the calculators for all samples.
 */
std::vector<genesis::population::DiversityPoolCalculator> get_diversity_calculators_(
    DiversityOptions const& options,
    BaseCountsFilter const& filter,
    std::vector<std::string> const& sample_names,
    std::vector<size_t> const& pool_sizes,
    DiversityOutputData const& output_data
) {
    // Prepare pool settings for each sample.
    // We here re-use the numerical filter settings as provided by the user.
    DiversityPoolSettings pool_settings;
    pool_settings.min_count    = filter.min_count;
    pool_settings.min_coverage = filter.min_coverage;
    pool_settings.max_coverage = filter.max_coverage;
    pool_settings.tajima_denominator_policy
        = options.with_popoolation_bugs.value
        ? TajimaDenominatorPolicy::kWithPopoolatioBugs
        : TajimaDenominatorPolicy::kWithoutPopoolatioBugs
    ;
    std::vector<DiversityPoolCalculator> sample_diversity_calculators;
    for( size_t i = 0; i < sample_names.size(); ++i ) {
        sample_diversity_calculators.emplace_back(
            pool_settings, pool_sizes[i]
        );
        sample_diversity_calculators.back().enable_theta_pi( output_data.compute_theta_pi );
        sample_diversity_calculators.back().enable_theta_watterson( output_data.compute_theta_wa );
        sample_diversity_calculators.back().enable_tajima_d( output_data.compute_tajima_d );
    }
    return sample_diversity_calculators;
}

// -------------------------------------------------------------------------
//     write_output_popoolation_
// -------------------------------------------------------------------------

/**
 * @brief Write output for PoPoolation style tables
 */
void write_output_popoolation_(
    std::vector<std::string> const& sample_names,
    VariantWindowView const& window,
    std::vector<DiversityPoolCalculator> const& sample_diversity_calculators,
    std::vector<BaseCountsFilterStats> const& sample_filter_stats,
    DiversityOutputData& output_data
) {
    // The popoolation table formats use the same format, so let's make one function to rule them all!
    // We take the DiversityPoolCalculator::Result here for the general values, and then the actual
    // data value again, so that we don't have to switch to get it.
    // Format: "2R	19500	0	0.000	na" or "A	1500	101	1.000	1.920886709" for example.
    auto write_popoolation_line_ = [](
        std::shared_ptr<genesis::utils::BaseOutputTarget>& ofs,
        VariantWindowView const& window,
        BaseCountsFilterStats const& stats,
        DiversityPoolCalculator::Result const& results,
        double value
    ) {
        internal_check(
            stats.passed == results.processed_count,
            "stats.passed != results.processed_count"
        );
        auto const coverage = static_cast<double>( stats.passed + stats.not_snp );
        auto const window_width = static_cast<double>( window.width() );
        auto const coverage_fraction = coverage / window_width;

        // Write fixed columns.
        (*ofs) << window.chromosome();
        (*ofs) << "\t" << anchor_position( window, WindowAnchorType::kIntervalMidpoint );
        (*ofs) << "\t" << results.processed_count;
        (*ofs) << "\t" << std::fixed << std::setprecision( 3 ) << coverage_fraction;
        if( std::isfinite( value ) ) {
            (*ofs) << "\t" << std::fixed << std::setprecision( 9 ) << value;
        } else {
            (*ofs) << "\tna";
        }
        (*ofs) << "\n";
    };

    // Write to all individual files for each sample and each value.
    for( size_t i = 0; i < sample_names.size(); ++i ) {
        auto const& div_calc = sample_diversity_calculators[i];
        auto const div_calc_result = div_calc.get_result( window.width() );

        // Theta Pi
        if( output_data.compute_theta_pi ) {
            write_popoolation_line_(
                output_data.popoolation_theta_pi_ofss[i],
                window,
                sample_filter_stats[i],
                div_calc_result,
                div_calc_result.theta_pi_relative
            );
        }

        // Theta Watterson
        if( output_data.compute_theta_wa ) {
            write_popoolation_line_(
                output_data.popoolation_theta_wa_ofss[i],
                window,
                sample_filter_stats[i],
                div_calc_result,
                div_calc_result.theta_watterson_relative
            );
        }

        // Tajima's D
        if( output_data.compute_tajima_d ) {
            write_popoolation_line_(
                output_data.popoolation_tajima_d_ofss[i],
                window,
                sample_filter_stats[i],
                div_calc_result,
                div_calc_result.tajima_d
            );
        }
    }
}

// -------------------------------------------------------------------------
//     write_output_table_
// -------------------------------------------------------------------------

/**
 * @brief Write output for our table format
 */
void write_output_table_(
    std::vector<std::string> const& sample_names,
    VariantWindowView const& window,
    std::vector<DiversityPoolCalculator> const& sample_diversity_calculators,
    std::vector<BaseCountsFilterStats> const& sample_filter_stats,
    DiversityOutputData& output_data
) {
    // Shorthand
    auto& table_ofs = output_data.table_ofs;

    // Helper function to write a field value to one of the tables.
    auto write_table_field_ = [&]( double value ){
        if( std::isfinite( value ) ) {
            (*table_ofs) << output_data.sep_char;
            (*table_ofs) << std::defaultfloat << std::setprecision( 9 ) << value;
        } else {
            (*table_ofs) << output_data.sep_char << output_data.na_entry;
        }
    };

    // Write fixed columns.
    (*table_ofs) << window.chromosome();
    (*table_ofs) << output_data.sep_char << window.first_position();
    (*table_ofs) << output_data.sep_char << window.last_position();

    // Write the per-pair diversity values in the correct order.
    for( size_t i = 0; i < sample_names.size(); ++i ) {
        auto const coverage = sample_filter_stats[i].passed + sample_filter_stats[i].not_snp;
        auto const window_width = static_cast<double>( window.width() );
        auto const coverage_fraction = static_cast<double>( coverage ) / window_width;

        auto const& div_calc = sample_diversity_calculators[i];
        auto const div_calc_result = div_calc.get_result( coverage );
        // auto const div_calc_result = div_calc.get_result( window.width() );

        // Meta info per sample
        // (*table_ofs) << output_data.sep_char << div_calc.variant_count;
        // (*table_ofs) << output_data.sep_char << div_calc.coverage_count;
        (*table_ofs) << output_data.sep_char << div_calc_result.processed_count;
        (*table_ofs) << output_data.sep_char << std::fixed << std::setprecision( 3 )
                     << coverage_fraction;

        // Values
        if( output_data.compute_theta_pi ) {
            write_table_field_( div_calc_result.theta_pi_absolute );
            write_table_field_( div_calc_result.theta_pi_relative );
        }
        if( output_data.compute_theta_wa ) {
            write_table_field_( div_calc_result.theta_watterson_absolute );
            write_table_field_( div_calc_result.theta_watterson_relative );
        }
        if( output_data.compute_tajima_d ) {
            write_table_field_( div_calc_result.tajima_d );
        }
    }
    (*table_ofs) << "\n";
}

// -------------------------------------------------------------------------
//     write_output_
// -------------------------------------------------------------------------

/**
 * @brief Write output to output_data
 */
void write_output_(
    DiversityOptions const& options,
    std::vector<std::string> const& sample_names,
    VariantWindowView const& window,
    std::vector<DiversityPoolCalculator> const& sample_diversity_calculators,
    std::vector<BaseCountsFilterStats> const& sample_filter_stats,
    DiversityOutputData& output_data
) {
    // Write the data, depending on the format.
    if( options.popoolation_format.value ) {
        write_output_popoolation_(
            sample_names, window, sample_diversity_calculators, sample_filter_stats, output_data
        );
    } else {
        write_output_table_(
            sample_names, window, sample_diversity_calculators, sample_filter_stats, output_data
        );
    }
}

// =================================================================================================
//      Run
// =================================================================================================

void run_diversity( DiversityOptions const& options )
{
    // -------------------------------------------------------------------------
    //     Settings
    // -------------------------------------------------------------------------

    // Warn the user about the PoPoolation bugs.
    if( options.with_popoolation_bugs.value ) {
        LOG_WARN << "Option `" << options.with_popoolation_bugs.option->get_name()
                 << "` was activated. Note that this yields numerically different results "
                 << "from the intended measure as described by the publication of Kofler et al, "
                 << "and that this option is provided solely for comparability with results "
                 << "obtained from PoPoolation.";
    }

    // Add a filter for just SNPs, per sample. Invariant sites don't change the absolute diversity.
    // Right now, we however use the invariant sites to compute relative Theta, and so we do not
    // want to filter them out beforehand. This is instead done in the actual processing,
    // sample by sample. Bit slower, but not by much.
    // BaseCountsFilter snp_filter;
    // snp_filter.only_snps = true;
    // options.variant_input.add_combined_filter_and_transforms(
    //     [snp_filter]( genesis::population::Variant& variant ){
    //         return genesis::population::filter_base_counts( variant, snp_filter );
    //     }
    // );

    // Prepare the base counts filter. We default the snps filter,
    // and then override everything else with the user provided values.
    BaseCountsFilter filter;
    filter.only_snps = true;
    filter = options.filter_numerical.get_sample_filter( filter ).first;

    // Get all samples names from the input file.
    auto const& sample_names = options.variant_input.sample_names();

    // Get the pool sizes.
    auto const pool_sizes = options.poolsizes.get_pool_sizes( sample_names );
    internal_check(
        pool_sizes.size() == sample_names.size(),
        "Inconsistent number of pool sizes and samples."
    );

    // Prepare file output_data and print headers.
    auto output_data = prepare_output_data_( options, sample_names );

    // Prepare pool settings for each sample.
    auto sample_diversity_calculators = get_diversity_calculators_(
        options, filter, sample_names, pool_sizes, output_data
    );
    std::vector<BaseCountsFilterStats> sample_filter_stats{ sample_names.size() };

    // -------------------------------------------------------------------------
    //     Main Loop
    // -------------------------------------------------------------------------

    // Iterate the file and compute per-window diversitye measures.
    // We run the samples in parallel, storing their results before writing to the output file.
    // For now, we compute all of them, in not the very most efficient way, but the easiest.
    auto window_it = options.window.get_variant_window_view_iterator( options.variant_input );
    for( auto cur_it = window_it->begin(); cur_it != window_it->end(); ++cur_it ) {
        auto& window = *cur_it;

        // Reset all calculator accumulators to zero for this window.
        for( size_t i = 0; i < sample_diversity_calculators.size(); ++i ) {
            sample_diversity_calculators[i].reset();
            reset( sample_filter_stats[i] );
        }

        // Compute diversity over samples.
        // #pragma omp parallel for
        for( auto& variant : window ) {
            internal_check(
                variant.samples.size() == sample_names.size(),
                "Inconsistent number of samples in input file."
            );

            // TODO refine the below to apply the filtering beforehand?!
            // that would make computing relative theta more tricky though.

            // Compute diversity for each sample.
            for( size_t i = 0; i < sample_names.size(); ++i ) {
                if( filter_base_counts( variant.samples[i], filter, sample_filter_stats[i] )) {
                    sample_diversity_calculators[i].process( variant.samples[i] );
                }
            }

            // It does not help here to parallelize, at least with the current layout...
            // We lose too much to the synchronization overhead.
            // if( sample_names.size() > 1 ) {
            //     global_options.thread_pool()->parallel_for(
            //         0, sample_names.size(),
            //         [&]( size_t i ){
            //             // Currently, we need to do some filtering here...
            //             if( filter_base_counts( variant.samples[i], filter, sample_filter_stats[i] )) {
            //                 sample_diversity_calculators[i].process( variant.samples[i] );
            //             }
            //         }
            //     ).wait();
            // } else if( sample_names.size() == 1 ) {
            //     if( filter_base_counts( variant.samples[0], filter, sample_filter_stats[0] )) {
            //         sample_diversity_calculators[0].process( variant.samples[0] );
            //     }
            // }
        }

        // Write the output to files.
        write_output_(
            options, sample_names, window,
            sample_diversity_calculators, sample_filter_stats, output_data
        );
    }

    // Final user output.
    options.window.print_report();
    options.filter_numerical.print_report();
}
