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

#include "options/table_output.hpp"

#include "options/global.hpp"
#include "tools/misc.hpp"

#include "genesis/utils/core/fs.hpp"
#include "genesis/utils/text/convert.hpp"
#include "genesis/utils/text/string.hpp"

#include <cassert>
#include <stdexcept>

// =================================================================================================
//      Setup Functions
// =================================================================================================

CLI::Option* TableOutputOptions::add_separator_char_opt_to_app(
    CLI::App* sub,
    std::string const& group
) {
    separator_char_.option = sub->add_option(
        "--separator-char",
        separator_char_.value,
        "Separator char between fields of output tabular data."
    )->transform(
        CLI::IsMember({ "comma", "tab", "space", "semicolon" }, CLI::ignore_case )
    );
    separator_char_.option->group( group );
    return separator_char_.option;
}

CLI::Option* TableOutputOptions::add_na_entry_opt_to_app(
    CLI::App* sub,
    std::string const& group
) {
    // Add an option to set the text for not-a-number.
    na_entry_.option = sub->add_option(
        "--na-entry",
        na_entry_.value,
        "Set the text to use in the output for n/a and NaN entries "
        "(e.g., resulting from positions with no counts, or windows with no variants). "
        "This is useful to match formatting expectations of downstream software."
    );
    na_entry_.option->group( group );
    return na_entry_.option;
}

// =================================================================================================
//      Run Functions
// =================================================================================================

char TableOutputOptions::get_separator_char() const
{
    // Get sep char.
    if(
        separator_char_.value == "comma" ||
        separator_char_.value == ","
    ) {
        return ',';
    } else if(
        separator_char_.value == "tab" ||
        separator_char_.value == "tabulator" ||
        separator_char_.value == "\t"
    ) {
        return '\t';
    } else if(
        separator_char_.value == "space" ||
        separator_char_.value == " "
    ) {
        return ' ';
    } else if(
        separator_char_.value == "semicolon" ||
        separator_char_.value == ";"
    ) {
        return ';';
    } else {
        throw CLI::ValidationError(
            "--separator-char",
            "Invalid separator char '" + separator_char_.value + "'."
        );
    }
}
