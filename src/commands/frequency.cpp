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

#include "commands/frequency.hpp"
#include "options/global.hpp"
#include "tools/cli_setup.hpp"

#include "genesis/population/functions/base_counts.hpp"

// =================================================================================================
//      Setup
// =================================================================================================

void setup_frequency( CLI::App& app )
{
    // Create the options and subcommand objects.
    auto options = std::make_shared<FrequencyOptions>();
    auto sub = app.add_subcommand(
        "frequency",
        "Create a table with per-sample and/or total base counts and frequencies "
        "at each position in the genome."
    );

    // Required input of some frequency format (mpileup or vcf at the moment).
    options->freq_input.add_frequency_input_opts_to_app( sub );
    options->freq_input.add_sample_name_opts_to_app( sub );
    options->freq_input.add_filter_opts_to_app( sub );

    // What type of columns to produce.
    options->type.option = sub->add_option(
        "--type",
        options->type.value,
        "Select which type of columns to output: Either per-sample counts and frequencies, "
        "or the total counts and frequencies across all samples, "
        "using the sum of the per-sample REF and ALT counts, or both (default)."
    );
    options->type.option->group( "Settings" );
    options->type.option->transform(
        CLI::IsMember({ "samples", "total", "both" }, CLI::ignore_case )
    );

    // Which columns to write and not to write
    options->no_coverage.option = sub->add_flag(
        "--no-coverage",
        options->no_coverage.value,
        "If set, do not write the 'COV' columns (per sample and/or total), "
        "which contain the coverage (sum of REF and ALT) counts."
    )->group( "Settings" );
    options->no_frequency.option = sub->add_flag(
        "--no-frequency",
        options->no_frequency.value,
        "If set, do not write write the 'FREQ' columns (per sample and/or total), "
        "which contain the frequency, computed as REF/(REF+ALT) counts."
    )->group( "Settings" );
    options->no_counts.option = sub->add_flag(
        "--no-counts",
        options->no_counts.value,
        "If set, do not write the 'REF_CNT' and 'ALT_CNT' columns (per sample and/or total), "
        "which contain the REF and ALT counts."
    )->group( "Settings" );

    // Invariant site handling
    options->omit_invariants.option = sub->add_flag(
        "--omit-invariants",
        options->omit_invariants.value,
        "If set, do not write rows that have no alternative counts. This is particularly useful "
        "for (m)pileup input, which otherwise produces an output row for each position."
    )->group( "Settings" );

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
        },
        [ options ]() {
            run_frequency( *options );
        }
    ));
}

// =================================================================================================
//      Run
// =================================================================================================

void run_frequency( FrequencyOptions const& options )
{
    using namespace genesis::population;

    options.file_output.check_output_files_nonexistence( "frequency", "csv" );
    auto freq_ofs = options.file_output.get_output_target( "frequency", "csv" );

    // Make the header fields.
    std::vector<std::string> fields;
    if( ! options.no_coverage.value ) {
        fields.emplace_back( "COV" );
    }
    if( ! options.no_frequency.value ) {
        fields.emplace_back( "FREQ" );
    }
    if( ! options.no_counts.value ) {
        fields.emplace_back( "REF_CNT" );
        fields.emplace_back( "ALT_CNT" );
    }
    if( fields.empty() ) {
        LOG_WARN << "Warning: All count and frequency output columns are deselected; "
                 << "the output will hence only contain the columns CHROM, POS, REF, ALT.";
    }

    // Get the separator char to use for table entries.
    auto const sep_char = options.table_output.get_separator_char();

    // Which types of columns to output: per sample and/or totals
    bool const write_columns = options.type.value == "samples" || options.type.value == "both";
    bool const write_total   = options.type.value == "total"   || options.type.value == "both";

    // Write the csv header line.
    (*freq_ofs) << "CHROM" << sep_char << "POS" << sep_char << "REF" << sep_char << "ALT";
    if( write_columns ) {
        for( auto const& sample : options.freq_input.sample_names() ) {
            for( auto const& field : fields ) {
                (*freq_ofs) << sep_char << sample << "." << field;
            }
        }
    }
    if( write_total ) {
        for( auto const& field : fields ) {
            (*freq_ofs) << sep_char << "TOTAL." << field;
        }
    }
    (*freq_ofs) << "\n";

    // Process the input file line by line and write the table data
    size_t line_cnt = 0;
    size_t skip_cnt = 0;
    for( auto const& freq_it : options.freq_input.get_iterator() ) {
        ++line_cnt;

        // If we want to omit invariant sites from the output, we need to do a prior check
        // whether the position is invariant or not. We only do this if needed, in order to
        // not slow down the case when we do not want to omit.
        if( options.omit_invariants.value ) {
            size_t ref_cnt = 0;
            size_t alt_cnt = 0;
            for( auto const& sample : freq_it.samples ) {
                ref_cnt += get_base_count( sample, freq_it.reference_base );
                alt_cnt += get_base_count( sample, freq_it.alternative_base );
            }
            if( ref_cnt == 0 || alt_cnt == 0 ) {
                ++skip_cnt;
                continue;
            }
        }

        // Write fixed columns
        (*freq_ofs) << freq_it.chromosome;
        (*freq_ofs) << sep_char << freq_it.position;
        (*freq_ofs) << sep_char << freq_it.reference_base;
        (*freq_ofs) << sep_char << freq_it.alternative_base;

        // Keep track of totals
        size_t total_ref_cnt = 0;
        size_t total_alt_cnt = 0;
        size_t total_cnt_sum = 0;

        // Per sample processing
        for( auto const& sample : freq_it.samples ) {

            // Get raw counts
            auto const ref_cnt = get_base_count( sample, freq_it.reference_base );
            auto const alt_cnt = get_base_count( sample, freq_it.alternative_base );
            auto const cnt_sum = ref_cnt + alt_cnt;
            total_ref_cnt += ref_cnt;
            total_alt_cnt += alt_cnt;
            total_cnt_sum += cnt_sum;

            // Write per sample columns
            if( write_columns ) {
                if( ! options.no_coverage.value ) {
                    (*freq_ofs) << sep_char << cnt_sum;
                }
                if( ! options.no_frequency.value ) {
                    if( cnt_sum > 0 ) {
                        auto const freq = static_cast<double>( ref_cnt ) / static_cast<double>( cnt_sum );
                        (*freq_ofs) << sep_char << freq;
                    } else {
                        (*freq_ofs) << sep_char << options.table_output.get_na_entry();
                    }
                }
                if( ! options.no_counts.value ) {
                    (*freq_ofs) << sep_char << ( ref_cnt );
                    (*freq_ofs) << sep_char << ( alt_cnt );
                }
            }
        }

        // Total columns
        if( write_total ) {
            if( ! options.no_coverage.value ) {
                (*freq_ofs) << sep_char << total_cnt_sum;
            }
            if( ! options.no_frequency.value ) {
                if( total_cnt_sum > 0 ) {
                    auto const freq
                        = static_cast<double>( total_ref_cnt )
                        / static_cast<double>( total_cnt_sum )
                    ;
                    (*freq_ofs) << sep_char << freq;
                } else {
                    (*freq_ofs) << sep_char << options.table_output.get_na_entry();
                }
            }
            if( ! options.no_counts.value ) {
                (*freq_ofs) << sep_char << ( total_ref_cnt );
                (*freq_ofs) << sep_char << ( total_alt_cnt );
            }
        }

        (*freq_ofs) << "\n";
    }

    // Output, depending on the setting, to keep it clean.
    if( options.omit_invariants.value ) {
        LOG_MSG << "Processed " << line_cnt << " genome positions of the input file, "
                << "and thereof skipped " << skip_cnt << " due to being invariant sites.";
    } else {
        LOG_MSG << "Processed " << line_cnt << " genome positions of the input file.";
    }
}
