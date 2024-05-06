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
    Lucas Czech <lczech@carnegiescience.edu>
    Department of Plant Biology, Carnegie Institution For Science
    260 Panama Street, Stanford, CA 94305, USA
*/

#include "options/variant_file_frequency_table.hpp"

#include "options/global.hpp"
#include "tools/misc.hpp"

#include "genesis/population/formats/frequency_table_input_stream.hpp"

#include <cassert>
#include <stdexcept>

// =================================================================================================
//      Setup Functions
// =================================================================================================

CLI::Option* VariantFileFrequencyTableOptions::add_file_input_opt_to_app_(
    CLI::App* sub,
    bool required,
    std::string const& group
) {
    // Correct setup check.
    internal_check(
        file_input_.option() == nullptr,
        "Cannot use the same VariantFileFrequencyTableOptions object multiple times."
    );

    // Add the base file input option
    file_input_.add_multi_file_input_opt_to_app(
        sub, "frequency-table", "frequency table",
        "(csv|tsv)(\\.gz)?",
        "(csv|tsv)[.gz]",
        required, group
    );

    // Separator char
    separator_char_.option = sub->add_option(
        "--frequency-table-separator-char",
        separator_char_.value,
        "Separator char between fields of the frequency table input."
    );
    separator_char_.option->group( group );
    separator_char_.option->transform(
        CLI::IsMember({ "comma", "tab", "space", "semicolon" }, CLI::ignore_case )
    );
    separator_char_.option->needs( file_input_.option() );

    // Int factor
    int_factor_.option = sub->add_option(
        "--frequency-table-cov-factor",
        int_factor_.value,
        "For frequency table input that only contains allele frequencies, without any information "
        "on coverage, we need to transform those frequencies into counts for our internal processing. "
        "This number is multiplied by the frequency to obtain these pseudo-counts. By default, "
        "we use 1000000, to get a reasonable interger approximation of the floating point frequency. "
        "This is of course above any typical coverage, but allows for more accurate counts when "
        "using for instance haplotype-corrected frequencies such as those from HAF-pipe."
    );
    int_factor_.option->group( group );
    int_factor_.option->transform( CLI::PositiveNumber );
    int_factor_.option->needs( file_input_.option() );

    // Flags frequency is alt
    frequency_is_ref_.option = sub->add_flag(
        "--frequency-table-freq-is-ref",
        frequency_is_ref_.value,
        "For frequency table input that contains allele frequencies, we need to decide whether those "
        "frequencies represent the reference or the alternative allele. By default, we assume the "
        "latter, i.e., values are interpreted as alternative allele frequencies. "
        "Use this flag to instead interpret them as reference allele frequencies."
    );
    frequency_is_ref_.option->group( group );
    frequency_is_ref_.option->needs( file_input_.option() );

    // Set the extra options
    if( add_extra_opts_ ) {
        add_extra_file_input_opts_to_app_( sub, group );
    }

    return file_input_.option();
}

void VariantFileFrequencyTableOptions::add_extra_file_input_opts_to_app_(
    CLI::App* sub,
    std::string const& group
) {
    // Chromosome column
    usr_chr_name_.option = sub->add_option(
        "--frequency-table-chr-column",
        usr_chr_name_.value,
        "Specify the name of the chromosome column in the header, case sensitive. "
        "By default, we look for columns named \"chromosome\", \"chrom\", \"chr\", or \"contig\", "
        "case insensitive."
        // "With this option however, the given string is searched instead in the header, and "
        // "the respective column is used for the chromosome information when parsing the table."
    );
    usr_chr_name_.option->group( group );
    usr_chr_name_.option->needs( file_input_.option() );

    // Position column
    usr_pos_name_.option = sub->add_option(
        "--frequency-table-pos-column",
        usr_pos_name_.value,
        "Specify the name of the position column in the header, case sensitive. "
        "By default, we look for columns named \"position\" or \"pos\", case insensitive."
    );
    usr_pos_name_.option->group( group );
    usr_pos_name_.option->needs( file_input_.option() );

    // Reference base column
    usr_ref_name_.option = sub->add_option(
        "--frequency-table-ref-base-column",
        usr_ref_name_.value,
        "Specify the name of the reference base column in the header, case sensitive. "
        "By default, we look for columns named \"reference\", \"referencebase\", \"ref\", or "
        "\"refbase\", case insensitive, and ignoring any extra punctuation marks."
    );
    usr_ref_name_.option->group( group );
    usr_ref_name_.option->needs( file_input_.option() );

    // Alternative base column
    usr_alt_name_.option = sub->add_option(
        "--frequency-table-alt-base-column",
        usr_alt_name_.value,
        "Specify the name of the alternative base column in the header, case sensitive. "
        "By default, we look for columns named \"alternative\", \"alternativebase\", \"alt\", or "
        "\"altbase\", case insensitive, and ignoring any extra punctuation marks."
    );
    usr_alt_name_.option->group( group );
    usr_alt_name_.option->needs( file_input_.option() );

    // Sample reference count column
    usr_smp_ref_name_.option = sub->add_option(
        "--frequency-table-sample-ref-count-column",
        usr_smp_ref_name_.value,
        "Specify the exact prefix or suffix of the per-sample reference count columns in the "
        "header, case sensitive. "
        "By default, we look leniently for column names that combine any of \"reference\", "
        "\"referencebase\", \"ref\", or \"refbase\" with any of \"counts\", \"count\", \"cnt\", or "
        "\"ct\", case insensitive, and ignoring any extra punctuation marks, as a prefix or suffix, "
        "with the remainder of the column name used as the sample name. For example, \"S1.ref_cnt\" "
        "indicates the reference count column for sample \"S1\"."
    );
    usr_smp_ref_name_.option->group( group );
    usr_smp_ref_name_.option->needs( file_input_.option() );

    // Sample alternative count column
    usr_smp_alt_name_.option = sub->add_option(
        "--frequency-table-sample-alt-count-column",
        usr_smp_alt_name_.value,
        "Specify the exact prefix or suffix of the per-sample alternative count columns in the "
        "header, case sensitive. "
        "By default, we look leniently for column names that combine any of \"alternative\", "
        "\"alternativebase\", \"alt\", or \"altbase\" with any of \"counts\", \"count\", \"cnt\", or "
        "\"ct\", case insensitive, and ignoring any extra punctuation marks, as a prefix or suffix, "
        "with the remainder of the column name used as the sample name. For example, \"S1.alt_cnt\" "
        "indicates the alternative count column for sample \"S1\"."
    );
    usr_smp_alt_name_.option->group( group );
    usr_smp_alt_name_.option->needs( file_input_.option() );

    // Sample frequency column
    usr_smp_frq_name_.option = sub->add_option(
        "--frequency-table-sample-freq-column",
        usr_smp_frq_name_.value,
        "Specify the exact prefix or suffix of the per-sample frequency columns in the "
        "header, case sensitive. "
        "By default, we look for column names having \"frequency\", \"freq\", \"maf\", \"af\", "
        "or \"allelefrequency\", case insensitive, and ignoring any extra punctuation marks, "
        "as a prefix or suffix, with the remainder of the column name used as the sample name. "
        "For example, \"S1.freq\" indicates the frequency column for sample \"S1\". "
        "Note that when the input data contains frequencies, but no reference or alternative base "
        "columns, such as HAF-pipe output tables, we cannot know the bases, and will hence guess. "
        "To properly set the reference bases, consider providing the `--reference-genome-fasta-file` option."
    );
    usr_smp_frq_name_.option->group( group );
    usr_smp_frq_name_.option->needs( file_input_.option() );

    // Sample coverage column
    usr_smp_cov_name_.option = sub->add_option(
        "--frequency-table-sample-cov-column",
        usr_smp_cov_name_.value,
        "Specify the exact prefix or suffix of the per-sample coverage (i.e., depth) columns "
        "in the header, case sensitive. "
        "By default, we look for column names having \"coverage\", \"cov\", \"depth\", or \"ad\", "
        "case insensitive, and ignoring any extra punctuation marks, as a prefix or suffix, "
        "with the remainder of the column name used as the sample name. "
        "For example, \"S1.cov\" indicates the coverage column for sample \"S1\"."
    );
    usr_smp_cov_name_.option->group( group );
    usr_smp_cov_name_.option->needs( file_input_.option() );
}

// =================================================================================================
//      Run Functions
// =================================================================================================

char VariantFileFrequencyTableOptions::get_separator_char_() const
{
    return translate_separator_char( separator_char_ );
}

VariantFileFrequencyTableOptions::VariantInputStream VariantFileFrequencyTableOptions::get_stream_(
    std::string const& filename
) const {
    using namespace genesis::population;

    // Prepare the reader with our settings. Sep char is set twice - not needed, but okay.
    FrequencyTableInputStream reader;
    reader.separator_char( get_separator_char_() );
    reader.frequency_is_ref( frequency_is_ref_.value );
    if( *int_factor_.option ) {
        reader.int_factor( int_factor_.value );
    }
    if( reference_genome_ ) {
        reader.reference_genome( reference_genome_ );
    }

    // If provided, set the user specified header name fragments.
    if( *usr_chr_name_.option ) {
        reader.header_chromosome_string( usr_chr_name_.value );
    }
    if( *usr_pos_name_.option ) {
        reader.header_position_string( usr_pos_name_.value );
    }
    if( *usr_ref_name_.option ) {
        reader.header_reference_base_string( usr_ref_name_.value );
    }
    if( *usr_alt_name_.option ) {
        reader.header_alternative_base_string( usr_alt_name_.value );
    }
    if( *usr_smp_ref_name_.option ) {
        reader.header_sample_reference_count_substring( usr_smp_ref_name_.value );
    }
    if( *usr_smp_alt_name_.option ) {
        reader.header_sample_alternative_count_substring( usr_smp_alt_name_.value );
    }
    if( *usr_smp_frq_name_.option ) {
        reader.header_sample_frequency_substring( usr_smp_frq_name_.value );
    }
    if( *usr_smp_cov_name_.option ) {
        reader.header_sample_coverage_substring( usr_smp_cov_name_.value );
    }

    // Prepare the stream.
    return make_variant_input_stream_from_frequency_table_file(
        filename, get_separator_char_(), reader
    );
}
