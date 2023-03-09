/*
    grenedalf - Genome Analyses of Differential Allele Frequencies
    Copyright (C) 2020-2022 Lucas Czech

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

#include "options/variant_input_frequency_table.hpp"

#include "options/global.hpp"
#include "options/variant_input_sample_names.hpp"
#include "tools/misc.hpp"

#include "genesis/population/formats/frequency_table_input_iterator.hpp"

#include <cassert>
#include <stdexcept>

// =================================================================================================
//      Setup Functions
// =================================================================================================

CLI::Option* VariantInputFrequencyTableOptions::add_file_input_opt_to_app_(
    CLI::App* sub,
    bool required,
    std::string const& group
) {
    // Correct setup check.
    internal_check(
        file_input_.option() == nullptr,
        "Cannot use the same VariantInputFrequencyTableOptions object multiple times."
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
        "--frequency-table-int-factor",
        int_factor_.value,
        "For frequency table input that only contains allele frequencies, without any information "
        "on coverage, we need to transform those frequencies into counts for our internal processing. "
        "This number is multiplied by the frequency to obtain these pseudo-counts. By default, "
        "we use the largest number that is possible with our numerical types."
    );
    int_factor_.option->group( group );
    int_factor_.option->transform( CLI::PositiveNumber );
    int_factor_.option->needs( file_input_.option() );

    // Flags frequency is alt
    frequency_is_alt_.option = sub->add_flag(
        "--frequency-table-frequency-is-alt",
        frequency_is_alt_.value,
        "For frequency table input that contains allele frequencies, we need to decide whether those "
        "frequencies represent the reference or the alternative allele. By default, we assume the "
        "former; use this flag to instead interpret them as alternative allele frequencies."
    );
    frequency_is_alt_.option->group( group );
    frequency_is_alt_.option->needs( file_input_.option() );

    return file_input_.option();
}

// =================================================================================================
//      Run Functions
// =================================================================================================

char VariantInputFrequencyTableOptions::get_separator_char_() const
{
    return translate_separator_char( separator_char_ );
}

VariantInputFrequencyTableOptions::VariantInputIterator VariantInputFrequencyTableOptions::get_iterator_(
    std::string const& filename,
    VariantInputSampleNamesOptions const& sample_names_options
) const {
    using namespace genesis::population;

    // Prepare the reader with our settings. Sep char is set twice - not needed, but okay.
    FrequencyTableInputIterator reader;
    reader.separator_char( get_separator_char_() );
    reader.int_factor( int_factor_.value );
    reader.frequency_is_ref( ! frequency_is_alt_.value );

    // Prepare the iterator.
    // See if we want to filter by sample name, and if so, resolve the name list.
    VariantInputIterator iterator;
    if( ! sample_names_options.get_filter_samples_include().value.empty() ) {
        auto const list = sample_names_options.process_sample_name_list_option(
            sample_names_options.get_filter_samples_include().value
        );
        iterator = make_variant_input_iterator_from_frequency_table_file(
            filename, list, false, get_separator_char_(), reader
        );
    } else if( ! sample_names_options.get_filter_samples_exclude().value.empty() ) {
        auto const list = sample_names_options.process_sample_name_list_option(
            sample_names_options.get_filter_samples_exclude().value
        );
        iterator = make_variant_input_iterator_from_frequency_table_file(
            filename, list, true, get_separator_char_(), reader
        );
    } else {
        iterator = make_variant_input_iterator_from_frequency_table_file(
            filename, get_separator_char_(), reader
        );
    }

    // As opposed to most other file formats, this contains sample names (only the filtered ones).
    // So here we do not need to set them, and can directly return.
    return iterator;
}
