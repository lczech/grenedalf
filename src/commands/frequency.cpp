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

#include "genesis/population/functions/functions.hpp"

#include <unordered_set>

// =================================================================================================
//      Setup
// =================================================================================================

void setup_frequency( CLI::App& app )
{
    // Create the options and subcommand objects.
    auto options = std::make_shared<FrequencyOptions>();
    auto sub = app.add_subcommand(
        "frequency",
        "Create a table with per-sample and/or total base counts and/or frequencies "
        "at positions in the genome."
    );

    // Required input of some frequency format (mpileup or vcf at the moment).
    options->variant_input.add_variant_input_opts_to_app( sub );
    // options->variant_input.add_sample_name_opts_to_app( sub );
    // options->variant_input.add_region_filter_opts_to_app( sub );

    // Add per-sample output selection.
    options->write_sample_counts.option = sub->add_flag(
        "--write-sample-counts",
        options->write_sample_counts.value,
        "If set, write 'REF_CNT' and 'ALT_CNT' columns per sample, "
        "containing the REF and ALT base counts at the position for each sample."
    )->group( "Settings" );
    options->write_sample_coverage.option = sub->add_flag(
        "--write-sample-coverage",
        options->write_sample_coverage.value,
        "If set, write a 'COV' column per sample, "
        "containing the coverage (sum of REF and ALT) counts of each sample."
    )->group( "Settings" );
    options->write_sample_ref_freq.option = sub->add_flag(
        "--write-sample-ref-freq",
        options->write_sample_ref_freq.value,
        "If set, write a 'FREQ' column per sample, containing the reference allele frequency, "
        "computed as REF/(REF+ALT) of the counts of each sample."
    )->group( "Settings" );
    options->write_sample_alt_freq.option = sub->add_flag(
        "--write-sample-alt-freq",
        options->write_sample_alt_freq.value,
        "If set, write a 'FREQ' column per sample, containing the alternative allele frequency, "
        "computed as ALT/(REF+ALT) of the counts of each sample."
    )->group( "Settings" );
    options->write_sample_ref_freq.option->excludes( options->write_sample_alt_freq.option );
    options->write_sample_alt_freq.option->excludes( options->write_sample_ref_freq.option );

    // Add total output selection.
    options->write_total_counts.option = sub->add_flag(
        "--write-total-counts",
        options->write_total_counts.value,
        "If set, write the 'REF_CNT' and 'ALT_CNT' columns for the total, "
        "which contain the REF and ALT base counts at the position across all samples."
    )->group( "Settings" );
    options->write_total_coverage.option = sub->add_flag(
        "--write-total-coverage",
        options->write_total_coverage.value,
        "If set, write the 'COV' column for the total, "
        "containing the coverage (sum of REF and ALT) counts across all samples."
    )->group( "Settings" );
    options->write_total_frequency.option = sub->add_flag(
        "--write-total-frequency",
        options->write_total_frequency.value,
        "If set, write the 'FREQ' column for the total, "
        "containing the frequency, computed as REF/(REF+ALT) of the counts across all samples."
    )->group( "Settings" );

    // Special cases handling
    options->write_invariants.option = sub->add_flag(
        "--write-invariants",
        options->write_invariants.value,
        "If set, write rows that have no alternative (ALT) counts. "
        "By default, we omit these positions, which is useful to keep the output small. "
        "For example sam/bam/cream or (m)pileup files otherwise produce an output row "
        "for each position in the input."
    )->group( "Settings" );
    options->omit_ref_alt_bases.option = sub->add_flag(
        "--omit-ref-and-alt-bases",
        options->omit_ref_alt_bases.value,
        "If set, do not write the columns containing the reference and alternative bases. "
        "This can be useful when the input is obtained from a source that does not contain "
        "those anyway. In that case, we internally assign them to 'A' and 'G', respectively, "
        "which usually is not correct, and hence should be omitted from the output."
    )->group( "Settings" );
    options->omit_alt_bases.option = sub->add_flag(
        "--omit-alt-bases",
        options->omit_alt_bases.value,
        "If set, do not write the column containing the alternative bases. "
        "This can be useful when the input is obtained from a source that does not contain "
        "them anyway. In that case, we internally assign the alternative base to be the "
        "transition base of the reference ('A' <-> 'G' and 'C' <-> 'T'), "
        "which usually is not correct, and hence should be omitted from the output. "
        "Note: To at least set the reference bases, "
        "consider providing the `--reference-genome-file` option."
    )->group( "Settings" );
    options->omit_ref_alt_bases.option->excludes( options->omit_alt_bases.option );
    options->omit_alt_bases.option->excludes( options->omit_ref_alt_bases.option );

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

    // Check the output file, and get write access to it.
    // We directly obtain the ostream here, to avoid having to dereference the target all the time.
    options.file_output.check_output_files_nonexistence( "frequency", "csv" );
    auto freq_oft = options.file_output.get_output_target( "frequency", "csv" );
    auto& freq_ofs = freq_oft->ostream();

    // If we have a reference genome, use it to correctly get the bases for the output here.
    if( options.variant_input.get_reference_genome() ) {
        auto const ref_gen = options.variant_input.get_reference_genome();
        options.variant_input.add_combined_filter_and_transforms(
            [ref_gen]( genesis::population::Variant& var ){
                guess_and_set_ref_and_alt_bases( var, *ref_gen );
                return true;
            }
        );
    } else {
        // Most of our input sources do not provide ref, and almost non provide alt bases.
        // So we use our guess function to augment the data. The function is idempotent
        // (unless we set the `force` parameter, which we do not do here), so for sources that do
        // contain ref and/or alt bases, nothing changes.
        options.variant_input.add_combined_filter_and_transforms(
            []( genesis::population::Variant& var ){
                guess_and_set_ref_and_alt_bases( var );
                return true;
            }
        );
    }

    // -------------------------------------------------------------------------
    //     Output preparation
    // -------------------------------------------------------------------------

    // Get the separator char to use for table entries.
    auto const sep_char = options.table_output.get_separator_char();
    auto const write_ref_bases = ! options.omit_ref_alt_bases.value;
    auto const write_alt_bases = ! options.omit_ref_alt_bases.value && ! options.omit_alt_bases.value;

    // Write the fixed part of the csv header line.
    freq_ofs << "CHROM" << sep_char << "POS";
    if( write_ref_bases ) {
        freq_ofs << sep_char << "REF";
    }
    if( write_alt_bases ) {
        freq_ofs << sep_char << "ALT";
    }
    bool write_anything = false;

    // Make the header fields for the samples.
    std::vector<std::string> fields;
    if( options.write_sample_counts.value ) {
        fields.emplace_back( "REF_CNT" );
        fields.emplace_back( "ALT_CNT" );
        write_anything = true;
    }
    if( options.write_sample_coverage.value ) {
        fields.emplace_back( "COV" );
        write_anything = true;
    }
    if( options.write_sample_ref_freq.value || options.write_sample_alt_freq.value ) {
        fields.emplace_back( "FREQ" );
        write_anything = true;
    }
    for( auto const& sample : options.variant_input.sample_names() ) {
        for( auto const& field : fields ) {
            freq_ofs << sep_char << sample << "." << field;
        }
    }

    // Write the header fields for the total.
    if( options.write_total_counts.value ) {
        freq_ofs << sep_char << "TOTAL.REF_CNT";
        freq_ofs << sep_char << "TOTAL.ALT_CNT";
        write_anything = true;
    }
    if( options.write_total_coverage.value ) {
        freq_ofs << sep_char << "TOTAL.COV";
        write_anything = true;
    }
    if( options.write_total_frequency.value ) {
        freq_ofs << sep_char << "TOTAL.FREQ";
        write_anything = true;
    }
    freq_ofs << "\n";

    // Let's help the user a bit.
    if( ! write_anything ) {
        throw CLI::ValidationError(
            "--write-...",
            "Error: No count or frequency output columns are selected; "
            "the output will hence not contain any actual frequency or count data."
        );
    }

    // -------------------------------------------------------------------------
    //     Loop over input data
    // -------------------------------------------------------------------------

    // We keep a list of chromosome names, for user output.
    std::unordered_set<std::string> chr_names;

    // Process the input file line by line and write the table data
    size_t line_cnt = 0;
    size_t skip_cnt = 0;
    size_t tot_off_cnt = 0;
    for( auto const& variant : options.variant_input.get_iterator() ) {
        ++line_cnt;

        // User output, so that we can see that something is happing.
        // It only outputs each chromosome once, even if they are not sorted.
        // That's probably okay for now. Might add code to issue a warning later when not sorted.
        if( chr_names.count( variant.chromosome ) == 0 ) {
            chr_names.insert( variant.chromosome );
            // LOG_MSG << "At chromosome " << variant.chromosome;
        }

        // If we want to omit invariant sites from the output, we need to do a prior check
        // whether the position is invariant or not. We only do this if needed, in order to
        // not slow down the case when we do not want to omit.
        if( ! options.write_invariants.value ) {
            size_t ref_cnt = 0;
            size_t alt_cnt = 0;
            size_t tot_cnt = 0;
            for( auto const& sample : variant.samples ) {
                ref_cnt += get_base_count( sample, variant.reference_base );
                alt_cnt += get_base_count( sample, variant.alternative_base );
                tot_cnt += total_base_count_sum( sample );
            }
            if( 2 * ( ref_cnt + alt_cnt ) < tot_cnt ) {
                // If the ref and alt counts make up less than half of the total counts,
                // we likely have some issue in the data, mismatching the bases. Report that.
                ++tot_off_cnt;
            }
            if( ref_cnt == 0 || alt_cnt == 0 ) {
                ++skip_cnt;
                continue;
            }
        }

        // Write fixed columns
        freq_ofs << variant.chromosome;
        freq_ofs << sep_char << variant.position;
        if( write_ref_bases ) {
            freq_ofs << sep_char << variant.reference_base;
        }
        if( write_alt_bases ) {
            freq_ofs << sep_char << variant.alternative_base;
        }

        // Keep track of totals
        size_t total_ref_cnt = 0;
        size_t total_alt_cnt = 0;
        size_t total_cnt_sum = 0;

        // Per sample columns processing
        for( auto const& sample : variant.samples ) {

            // Get raw counts
            auto const ref_cnt = get_base_count( sample, variant.reference_base );
            auto const alt_cnt = get_base_count( sample, variant.alternative_base );
            auto const cnt_sum = ref_cnt + alt_cnt;
            total_ref_cnt += ref_cnt;
            total_alt_cnt += alt_cnt;
            total_cnt_sum += cnt_sum;

            // Write per sample columns
            if( options.write_sample_counts.value ) {
                freq_ofs << sep_char << ref_cnt;
                freq_ofs << sep_char << alt_cnt;
            }
            if( options.write_sample_coverage.value ) {
                freq_ofs << sep_char << cnt_sum;
            }
            if( options.write_sample_ref_freq.value || options.write_sample_alt_freq.value ) {
                if( cnt_sum > 0 ) {
                    double freq = 0.0;
                    if( options.write_sample_ref_freq.value ) {
                        freq = static_cast<double>( ref_cnt ) / static_cast<double>( cnt_sum );
                    } else if( options.write_sample_alt_freq.value ) {
                        freq = static_cast<double>( alt_cnt ) / static_cast<double>( cnt_sum );
                    }
                    freq_ofs << sep_char << freq;
                } else {
                    freq_ofs << sep_char << options.table_output.get_na_entry();
                }
            }
        }

        // Total columns
        if( options.write_total_counts.value ) {
            freq_ofs << sep_char << total_ref_cnt;
            freq_ofs << sep_char << total_alt_cnt;
        }
        if( options.write_total_coverage.value ) {
            freq_ofs << sep_char << total_cnt_sum;
        }
        if( options.write_total_frequency.value ) {
            if( total_cnt_sum > 0 ) {
                auto const freq
                    = static_cast<double>( total_ref_cnt )
                    / static_cast<double>( total_cnt_sum )
                ;
                freq_ofs << sep_char << freq;
            } else {
                freq_ofs << sep_char << options.table_output.get_na_entry();
            }
        }

        freq_ofs << "\n";
    }

    // -------------------------------------------------------------------------
    //     Final user output
    // -------------------------------------------------------------------------

    // Output, depending on the setting, to keep it clean.
    if( options.write_invariants.value ) {
        LOG_MSG << "Processed " << line_cnt << " genome positions of the input file, "
                << "from " << chr_names.size() << " chromosomes.";
    } else {
        LOG_MSG << "Processed " << line_cnt << " genome positions of the input file, "
                << "from " << chr_names.size() << " chromosomes, "
                << "and thereof skipped " << skip_cnt << " due to being invariant sites.";
    }
    if( tot_off_cnt > line_cnt / 100 ) {
        // Report this only if it happened in more thatn 1% of the cases.
        LOG_WARN << "There were " << tot_off_cnt << " positions where the sum of the reference "
                 << "and alternative base counts were less than have of the total counts at "
                 << "the position. That is, the two remaining nucleotide counts were more in "
                 << "sum than the two that we used here. This likely indicates some data error, "
                 << "such as a mismatch of the reference base. Please check your data.";
    }
}
