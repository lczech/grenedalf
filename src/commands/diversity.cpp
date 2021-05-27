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
#include "genesis/utils/text/string.hpp"

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

    // Minimum allele count
    options->min_allele_count.option = sub->add_option(
        "--min-allele-count",
        options->min_allele_count.value,
        "Minimum allele count of the minor allele. Used for the identification of SNPs."
    )->group( "Settings" );

    // Minimum coverage
    options->min_coverage.option = sub->add_option(
        "--min-coverage",
        options->min_coverage.value,
        "Minimum coverage of a site. Sites with a lower coverage will not be considered "
        "for SNP identification and coverage estimation."
    )->group( "Settings" );

    // Maximum coverage
    options->max_coverage.option = sub->add_option(
        "--max-coverage",
        options->max_coverage.value,
        "Maximum coverage used for SNP identification. Coverage in ALL populations has to be lower "
        "or equal to this threshold, otherwise no SNP will be called."
    )->group( "Settings" );

    // Minimum coverage fraction
    options->min_coverage_fraction.option = sub->add_option(
        "--min-coverage-fraction",
        options->min_coverage_fraction.value,
        "Minimum coverage fraction of a window being between `--min-coverage` and `--max-coverage` "
        "in ALL populations that needs to be reached in order to compute the diversity measures."
    )->group( "Settings" );

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
//      Run
// =================================================================================================

void run_diversity( DiversityOptions const& options )
{
    using namespace genesis::population;
    using namespace genesis::utils;
    using BaseCountWindow = Window<std::vector<BaseCounts>>;

    // Get all samples names from the input file.
    auto const& sample_names = options.freq_input.sample_names();

    // Get the measures to compute. At least one of them will be active.
    auto const measure = to_lower( options.measure.value );
    bool const compute_theta_pi = ( measure == "all" || measure == "theta-pi" );
    bool const compute_theta_wa = ( measure == "all" || measure == "theta-watterson" );
    bool const compute_tajima_d = ( measure == "all" || measure == "tajimas-d" );
    assert( compute_theta_pi || compute_theta_wa || compute_tajima_d );

    // Output file checks.
    if( options.popoolation_format.value ) {
        for( auto const& sample_name : sample_names ) {
            if( compute_theta_pi ) {
                options.file_output.check_output_files_nonexistence(
                    "diversity-" + sample_name + "-theta-pi", "csv"
                );
            }
            if( compute_theta_wa ) {
                options.file_output.check_output_files_nonexistence(
                    "diversity-" + sample_name + "-theta-watterson", "csv"
                );
            }
            if( compute_tajima_d ) {
                options.file_output.check_output_files_nonexistence(
                    "diversity-" + sample_name + "-tajimas-d", "csv"
                );
            }
        }
    } else {
        options.file_output.check_output_files_nonexistence( "diversity", "csv" );
    }

    // -------------------------------------------------------------------------
    //     Settings
    // -------------------------------------------------------------------------

    // Warn the user about the PoPoolation bugs.
    if( options.with_popoolation_bugs.value ) {
        LOG_WARN << "Option `" << options.with_popoolation_bugs.option->get_name() << "` was activated. "
                 << "Note that this yields numerically different results and that this option is "
                 << "provided solely for comparability with results obtained from PoPoolation.";
    }

    // Get the pool sizes.
    auto const pool_sizes = options.poolsizes.get_pool_sizes( sample_names );

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

    // Get the separator char to use for table entries.
    auto const sep_char = options.table_output.get_separator_char();

    // -------------------------------------------------------------------------
    //     File Output and Table Header
    // -------------------------------------------------------------------------

    // Prepare output file and write fixed header fields. We need all pointers, for our
    // and for the PoPoolation format, in order to avoid code duplication in the computation.
    std::shared_ptr<genesis::utils::BaseOutputTarget> table_ofs;
    std::vector<std::shared_ptr<genesis::utils::BaseOutputTarget>> popoolation_theta_pi_ofss;
    std::vector<std::shared_ptr<genesis::utils::BaseOutputTarget>> popoolation_theta_wa_ofss;
    std::vector<std::shared_ptr<genesis::utils::BaseOutputTarget>> popoolation_tajima_d_ofss;

    if( options.popoolation_format.value ) {

        // Open the three PoPoolation file formats for each sample.
        for( auto const& sample_name : sample_names ) {
            if( compute_theta_pi ) {
                popoolation_theta_pi_ofss.emplace_back(
                    options.file_output.get_output_target(
                        "diversity-" + sample_name + "-theta-pi", "csv"
                    )
                );
            }
            if( compute_theta_wa ) {
                popoolation_theta_wa_ofss.emplace_back(
                    options.file_output.get_output_target(
                        "diversity-" + sample_name + "-theta-watterson", "csv"
                    )
                );
            }
            if( compute_tajima_d ) {
                popoolation_tajima_d_ofss.emplace_back(
                    options.file_output.get_output_target(
                        "diversity-" + sample_name + "-tajimas-d", "csv"
                    )
                );
            }
        }

    } else {

        // Open and prepare our table format
        table_ofs = options.file_output.get_output_target( "diversity", "csv" );
        (*table_ofs) << "chrom" << sep_char << "start" << sep_char << "end";

        // Make the header per-sample fields.
        std::vector<std::string> fields;

        // fields.push_back( "variant_count" );
        // fields.push_back( "coverage_count" );
        fields.push_back( "snp_count" );
        fields.push_back( "coverage_fraction" );

        if( compute_theta_pi ) {
            fields.push_back( "theta_pi_abs" );
            fields.push_back( "theta_pi_rel" );
        }
        if( compute_theta_wa ) {
            fields.push_back( "theta_watterson_abs" );
            fields.push_back( "theta_watterson_rel" );
        }
        if( compute_tajima_d ) {
            fields.push_back( "tajimas_d" );
        }

        // Write all fields for all samples.
        for( auto const& sample : sample_names ) {
            for( auto const& field : fields ) {
                (*table_ofs) << sep_char << sample << "." << field;
            }
        }
        (*table_ofs) << "\n";
    }

    // -------------------------------------------------------------------------
    //     Write Helper Functions
    // -------------------------------------------------------------------------

    // Helper function to write a field value to one of the tables.
    // Only used for our table format, as PoPoolation needs a bit of a different formatting.
    auto write_table_field_ = [&]( double value ){
        if( std::isfinite( value ) ) {
            (*table_ofs) << sep_char << std::defaultfloat << std::setprecision( 9 ) << value;
        } else {
            (*table_ofs) << sep_char << options.table_output.get_na_entry();
        }
    };

    // The popoolation table formats use the same format, so let's make one function to rule them all!
    // We take the PoolDiversityResults here for the general values, and then the actual data
    // value again, so that we don't have to switch to get it.
    // Format: "2R	19500	0	0.000	na" or "A	1500	101	1.000	1.920886709" for example.
    auto write_popoolation_line_ = [](
        std::shared_ptr<genesis::utils::BaseOutputTarget>& ofs,
        BaseCountWindow const& window,
        PoolDiversityResults const& results,
        double value
    ){
        // Write fixed columns.
        (*ofs) << window.chromosome();
        (*ofs) << "\t" << window.anchor_position( WindowAnchorType::kIntervalMidpoint );
        (*ofs) << "\t" << results.snp_count;
        (*ofs) << "\t" << std::fixed << std::setprecision( 3 ) << results.coverage_fraction;
        if( std::isfinite( value ) ) {
            (*ofs) << "\t" << std::fixed << std::setprecision( 9 ) << value;
        } else {
            (*ofs) << "\tna";
        }
        (*ofs) << "\n";
    };

    // -------------------------------------------------------------------------
    //     Main Loop
    // -------------------------------------------------------------------------

    // A bit of user output to keep 'em happy. Chromsomes, windows, positions.
    size_t chr_cnt = 0;
    size_t win_cnt = 0;
    size_t pos_cnt = 0;

    // Iterate the file and compute per-window diversitye measures.
    // We run the samples in parallel, storing their results before writing to the output file.
    // For now, we compute all of them, in not the very most efficient way, but the easiest.
    auto window_it = options.freq_input.get_base_count_sliding_window_iterator();
    auto sample_divs = std::vector<PoolDiversityResults>( sample_names.size() );
    for( ; window_it; ++window_it ) {
        auto const& window = *window_it;
        ++win_cnt;
        pos_cnt += window.size();

        // Some user output to report progress.
        if( window_it.is_first_window() ) {
            LOG_MSG << "At chromosome " << window.chromosome();
            ++chr_cnt;
        }
        LOG_MSG2 << "    At window "
                 << window.chromosome() << ":"
                 << window.first_position() << "-"
                 <<  window.last_position();

        // Skip empty windows if the user wants to.
        // if( window.empty() && options.omit_empty_windows.value ) {
        //     continue;
        // }

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

            // Compute diversity measures for the sample. We always compute all measures,
            // even if not all of them will be written afterwards. It's fast enough anyway,
            // and most of the compute time is spent in parsing, so that's okay and easier.
            sample_divs[i] = pool_diversity_measures( pool_settings[i], range.begin(), range.end() );
        }

        // Write the data, depending on the format.
        if( options.popoolation_format.value ) {

            // Write to all individual files for each sample and each value.
            for( size_t i = 0; i < sample_divs.size(); ++i ) {
                // Theta Pi
                if( compute_theta_pi ) {
                    write_popoolation_line_(
                        popoolation_theta_pi_ofss[i],
                        window,
                        sample_divs[i],
                        sample_divs[i].theta_pi_relative
                    );
                }

                // Theta Watterson
                if( compute_theta_wa ) {
                    write_popoolation_line_(
                        popoolation_theta_wa_ofss[i],
                        window,
                        sample_divs[i],
                        sample_divs[i].theta_watterson_relative
                    );
                }

                // Tajima's D
                if( compute_tajima_d ) {
                    write_popoolation_line_(
                        popoolation_tajima_d_ofss[i],
                        window,
                        sample_divs[i],
                        sample_divs[i].tajima_d
                    );
                }
            }

        } else {

            // Write fixed columns.
            (*table_ofs) << window.chromosome();
            (*table_ofs) << sep_char << window.first_position();
            (*table_ofs) << sep_char << window.last_position();

            // Write the per-pair diversity values in the correct order.
            for( auto const& sample_div : sample_divs ) {
                // Meta info per sample
                // (*table_ofs) << sep_char << sample_div.variant_count;
                // (*table_ofs) << sep_char << sample_div.coverage_count;
                (*table_ofs) << sep_char << sample_div.snp_count;
                (*table_ofs) << sep_char << std::fixed << std::setprecision( 3 )
                             << sample_div.coverage_fraction;

                // Values
                if( compute_theta_pi ) {
                    write_table_field_( sample_div.theta_pi_absolute );
                    write_table_field_( sample_div.theta_pi_relative );
                }
                if( compute_theta_wa ) {
                    write_table_field_( sample_div.theta_watterson_absolute );
                    write_table_field_( sample_div.theta_watterson_relative );
                }
                if( compute_tajima_d ) {
                    write_table_field_( sample_div.tajima_d );
                }
            }
            (*table_ofs) << "\n";
        }
    }

    LOG_MSG << "\nProcessed " << chr_cnt << " chromosome" << ( chr_cnt != 1 ? "s" : "" )
            << " with " << pos_cnt << " total position" << ( pos_cnt != 1 ? "s" : "" )
            << " in " << win_cnt << " window" << ( win_cnt != 1 ? "s" : "" );
}
