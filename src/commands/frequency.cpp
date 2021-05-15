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
        "Create a table with per-sample base counts and frequencies at each position in the genome."
    );

    // Required input of some frequency format (mpileup or vcf at the moment).
    options->freq_input.add_frequency_input_opts_to_app( sub );

    // Which columns to write
    auto write_coverage_opt = sub->add_flag(
        "--write-coverage",
        options->write_coverage,
        "If set, write a column 'COV' per sample containing the coverage (sum of REF and ALT) counts."
    )->group( "Settings" );
    auto write_freq_opt = sub->add_flag(
        "--write-frequency",
        options->write_frequency,
        "If set, write a column 'FREQ' per sample containing the frequency, "
        "computed as REF/(REF+ALT) counts."
    )->group( "Settings" );
    auto write_counts_opt = sub->add_flag(
        "--write-counts",
        options->write_counts,
        "If set, write columns 'REF_CNT' and 'ALT_CNT' per sample containing the REF and ALT counts."
    )->group( "Settings" );

    // Shortcut for all.
    auto write_all_opt = sub->add_flag(
        "--write-all",
        options->write_all,
        "If set, write all the above columns (COV, FREQ, REF_CNT, ALT_CNT)."
    )->group( "Settings" );
    write_all_opt->excludes( write_coverage_opt );
    write_all_opt->excludes( write_freq_opt );
    write_all_opt->excludes( write_counts_opt );

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
    if( options.write_coverage || options.write_all ) {
        fields.emplace_back( "COV" );
    }
    if( options.write_frequency || options.write_all ) {
        fields.emplace_back( "FREQ" );
    }
    if( options.write_counts || options.write_all ) {
        fields.emplace_back( "REF_CNT" );
        fields.emplace_back( "ALT_CNT" );
    }
    if( fields.empty() ) {
        LOG_WARN << "Warning: No output columns are selected; the output will hence only contain "
                 << "the columns CHROM, POS, REF, ALT. Use the --write-... options to select "
                 << "which additional columns to write.";
    }

    // Get the separator char to use for table entries.
    auto const sep_char = options.table_output.get_separator_char();

    // Write the csv header line.
    (*freq_ofs) << "CHROM" << sep_char << "POS" << sep_char << "REF" << sep_char << "ALT";
    for( auto const& sample : options.freq_input.sample_names() ) {
        for( auto const& field : fields ) {
            (*freq_ofs) << sep_char << sample << "." << field;
        }
    }
    (*freq_ofs) << "\n";

    // Write the table data
    for( auto const& freq_it : options.freq_input.get_iterator() ) {
        (*freq_ofs) << freq_it.chromosome;
        (*freq_ofs) << sep_char << freq_it.position;
        (*freq_ofs) << sep_char << freq_it.reference_base;
        (*freq_ofs) << sep_char << freq_it.alternative_base;

        for( auto const& sample : freq_it.samples ) {
            auto const ref_cnt = get_base_count( sample, freq_it.reference_base );
            auto const alt_cnt = get_base_count( sample, freq_it.alternative_base );
            auto const cnt_sum = ref_cnt + alt_cnt;

            if( options.write_coverage || options.write_all ) {
                (*freq_ofs) << sep_char << cnt_sum;
            }
            if( options.write_frequency || options.write_all ) {
                if( cnt_sum > 0 ) {
                    auto const freq = static_cast<double>( ref_cnt ) / static_cast<double>( cnt_sum );
                    (*freq_ofs) << sep_char << freq;
                } else {
                    (*freq_ofs) << sep_char << options.table_output.get_na_entry();
                }
            }
            if( options.write_counts || options.write_all ) {
                (*freq_ofs) << sep_char << ( ref_cnt );
                (*freq_ofs) << sep_char << ( alt_cnt );
            }
        }

        (*freq_ofs) << "\n";
    }
}
