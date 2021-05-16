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

#include "commands/diversity.hpp"
#include "options/global.hpp"
#include "tools/cli_setup.hpp"
#include "tools/misc.hpp"

#include "genesis/population/functions/base_counts.hpp"
#include "genesis/population/functions/diversity.hpp"
#include "genesis/population/functions/variant.hpp"

#include <cassert>

// =================================================================================================
//      Setup
// =================================================================================================

void setup_diversity( CLI::App& app )
{
    // Create the options and subcommand objects.
    auto options = std::make_shared<DiversityOptions>();
    auto sub = app.add_subcommand(
        "diversity",
        "Compute the pool-sequencing corrected diversity measures Theta Pi, "
        "Theta Watterson, and Tajima's D, following the PoPoolation approach, "
        "for each sample, in windows along the genome."
    );

    // -------------------------------------------------------------------------
    //     Input
    // -------------------------------------------------------------------------

    // Required input of some frequency format, and settings for the sliding window.
    options->freq_input.add_frequency_input_opts_to_app( sub );
    options->freq_input.add_sample_name_opts_to_app( sub );
    options->freq_input.add_filter_opts_to_app( sub );
    options->freq_input.add_sliding_window_opts_to_app( sub );

    // -------------------------------------------------------------------------
    //     Settings
    // -------------------------------------------------------------------------

    // Settings: Pool sizes
    options->poolsizes.add_poolsizes_opt_to_app( sub );

    // Minimum allele count
    options->min_allele_count.option = sub->add_option(
        "--min-allele-count",
        options->min_allele_count.value,
        "Minimum allele count"
    )->group( "Settings" );

    // Minimum coverage
    options->min_coverage.option = sub->add_option(
        "--min-coverage",
        options->min_coverage.value,
        "Minimum coverage"
    )->group( "Settings" );

    // Maximum coverage
    options->max_coverage.option = sub->add_option(
        "--max-coverage",
        options->max_coverage.value,
        "Maximum coverage"
    )->group( "Settings" );

    // Minimum coverage fraction
    options->min_coverage_fraction.option = sub->add_option(
        "--min-coverage-fraction",
        options->min_coverage_fraction.value,
        "Minimum coverage fraction"
    )->group( "Settings" );

    // Add an option to purposely activate the PoPollation bugs that we discovered.
    // The option is hidden for now, to not raise too much awareness before the issue is solved
    // in a way that works well for everyone.
    options->with_popoolation_bugs.option = sub->add_flag(
        "--with-popoolation-bugs",
        options->with_popoolation_bugs.value,
        "If set, imitate the two major bugs of PoPoolation that we (likely) found, which yield "
        "numerically different results. We offer this option for comparability with PoPoolation."
    )->group("");

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
            "Kofler2011-PoPoolation"
        },
        [ options ]() {
            run_diversity( *options );
        }
    ));
}

// =================================================================================================
//      Run
// =================================================================================================

void run_diversity( DiversityOptions const& options )
{
    using namespace genesis::population;
    using namespace genesis::utils;
    using BaseCountWindow = Window<std::vector<BaseCounts>>;

    // Output preparation.
    options.file_output.check_output_files_nonexistence( "diversity", "csv" );

    // Get the separator char to use for table entries.
    auto const sep_char = options.table_output.get_separator_char();

    // -------------------------------------------------------------------------
    //     Settings
    // -------------------------------------------------------------------------

    // Get the pool sizes and sample names.
    auto const& sample_names = options.freq_input.sample_names();
    auto const pool_sizes = options.poolsizes.get_pool_sizes( sample_names );

    if( options.with_popoolation_bugs.value ) {
        LOG_WARN << "Option `--with-popoolation-bugs` was activated. "
                 << "Note that this yields numerically different results and that this option is "
                 << "provided solely for comparability with results obtained with PoPoolation.";
    }

    // Prepare pool settings for each sample. We need to copy the options over to our settings class.
    // Bit cumbersome, but makes the internal handling easier...
    // Some are commented out, as they are currently not used by any diversity computation.
    auto pool_settings = std::vector<PoolDiversitySettings>( sample_names.size() );
    auto const window_width_and_stride = options.freq_input.get_window_width_and_stride();
    for( size_t i = 0; i < sample_names.size(); ++i ) {
        pool_settings[i].window_width          = window_width_and_stride.first;
        pool_settings[i].window_stride         = window_width_and_stride.second;
        // pool_settings[i].min_phred_score       = options.freq_input...;
        pool_settings[i].poolsize              = pool_sizes[i];
        pool_settings[i].min_allele_count      = options.min_allele_count.value;
        pool_settings[i].min_coverage          = options.min_coverage.value;
        pool_settings[i].max_coverage          = options.max_coverage.value;
        pool_settings[i].min_coverage_fraction = options.min_coverage_fraction.value;
        pool_settings[i].with_popoolation_bugs = options.with_popoolation_bugs.value;
    }

    // -------------------------------------------------------------------------
    //     Table Header
    // -------------------------------------------------------------------------
    // Prepare output file and write fixed header fields.
    auto div_ofs = options.file_output.get_output_target( "diversity", "csv" );
    (*div_ofs) << "CHROM" << sep_char << "START" << sep_char << "END" << sep_char << "SNPS";

    // Make the header per-sample fields.
    std::vector<std::string> fields;
    fields.push_back( "theta_pi_abs" );
    fields.push_back( "theta_pi_rel" );
    fields.push_back( "theta_watterson_abs" );
    fields.push_back( "theta_watterson_rel" );
    fields.push_back( "tajimas_d" );
    fields.push_back( "variant_count" );
    fields.push_back( "snp_count" );
    fields.push_back( "coverage_count" );
    fields.push_back( "coverage_fraction" );

    // Write all fields for all samples.
    for( auto const& sample : sample_names ) {
        for( auto const& field : fields ) {
            (*div_ofs) << sep_char << sample << "." << field;
        }
    }
    (*div_ofs) << "\n";

    // -------------------------------------------------------------------------
    //     Main Loop
    // -------------------------------------------------------------------------

    // Helper function to write a field value to the table.
    auto write_field = [&]( double value ){
        if( std::isfinite( value ) ) {
            (*div_ofs) << sep_char << value;
        } else {
            (*div_ofs) << sep_char << options.table_output.get_na_entry();
        }
    };

    // Iterate the file and compute per-window diversitye measures.
    // We run the samples in parallel, storing their results before writing to the output file.
    // For now, we compute all of them, in not the very most efficient way, but the easiest.
    auto window_it = options.freq_input.get_base_count_sliding_window_iterator();
    auto window_div = std::vector<PoolDiversityResults>( sample_names.size() );
    for( ; window_it; ++window_it ) {
        auto const& window = *window_it;

        // Some user output to report progress.
        if( window_it.is_first_window() ) {
            LOG_MSG << "At chromosome " << window.chromosome();
        }

        // Skip empty windows if the user wants to.
        // if( window.empty() && options.omit_empty_windows.value ) {
        //     continue;
        // }

        // Write fixed columns.
        (*div_ofs) << window.chromosome();
        (*div_ofs) << sep_char << window.first_position();
        (*div_ofs) << sep_char << window.last_position();
        (*div_ofs) << sep_char << window.entry_count();

        // Compute diversity in parallel over samples.
        #pragma omp parallel for
        for( size_t i = 0; i < sample_names.size(); ++i ) {

            // Select sample i within the current window.
            auto range = make_transform_range(
                [i]( BaseCountWindow::Entry const& entry ) -> BaseCounts const& {
                    internal_check(
                        i < entry.data.size(),
                        "Inconsistent number of samples in input file."
                    );
                    return entry.data[i];
                },
                window.begin(), window.end()
            );

            // Compute diversity measures for the sample.
            window_div[i] = pool_diversity_measures( pool_settings[i], range.begin(), range.end() );
        }

        // Write the per-pair diversity values in the correct order.
        for( auto const& div : window_div ) {
            write_field( div.theta_pi_absolute );
            write_field( div.theta_pi_relative );
            write_field( div.theta_watterson_absolute );
            write_field( div.theta_watterson_relative );
            write_field( div.tajima_d );
            write_field( div.variant_count );
            write_field( div.snp_count );
            write_field( div.coverage_count );
            write_field( div.coverage_fraction );
        }
        (*div_ofs) << "\n";
    }
}
