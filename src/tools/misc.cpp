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

#include "tools/misc.hpp"

#include "genesis/utils/core/info.hpp"
#include "genesis/utils/core/options.hpp"
#include "genesis/utils/text/string.hpp"

#include <cmath>
#include <cstdlib>
#include <stdexcept>

// =================================================================================================
//      Legacy Commands
// =================================================================================================

void add_legacy_command(
    CLI::App& app,
    std::string const& old_name,
    std::string const& new_path
) {
    auto sub = app.add_subcommand(
        old_name,
        "Command has been renamed to `grenedalf " + new_path + "`"
    );

    sub->group("");
    sub->callback( [ new_path ]() {
        throw RenamedCommandError( "Command has been renamed to `grenedalf " + new_path + "`" );
    });
}

// =================================================================================================
//      Formatting
// =================================================================================================

std::string format_columns(
    std::string const& left,
    std::string const& right,
    size_t left_w
) {
    // If std out is a terminal, we use its width for the maximal line length.
    unsigned long twidth = 0;
    if( genesis::utils::info_stdout_is_terminal() ) {
        twidth = static_cast<unsigned long>(
            genesis::utils::info_terminal_size().first
        );
    }

    // Get the widths of the columns. If there is not room for the second one,
    // make it 0 length, meaning all is written in one line.
    auto const lwidth = left_w;
    auto const rwidth = ( twidth > lwidth ? twidth - lwidth : 0 );

    // Write.
    std::stringstream out;
    write_columns( out, left, right, lwidth, rwidth );
    return out.str();
}

void write_columns(
    std::ostream& out,
    std::string const& left,
    std::string const& right,
    size_t left_w,
    size_t right_w
) {
    // Write left column.
    out << std::setw(static_cast<int>( left_w )) << std::left << left;

    // Write right column.
    if( ! right.empty() ) {
        auto rcpy = right;

        // If the left is already longer than the column allows, start a new line.
        if( left.length() >= left_w ) {
            out << "\n" << std::setw(static_cast<int>( left_w )) << "";
        }

        // If we have an actual useful width for the right column, wrap it.
        // Otherwise, we just put everything in one long line.
        if( right_w > 0 ) {
            rcpy = genesis::utils::wrap( rcpy, right_w );
        }

        // Indent and then trim again. The trimming removes the leading whitespace,
        // which we do not want, as we already inserted enough, and it removes
        // the trailing new line from the wrapping, which we do not want, as we output
        // one later anyway.
        rcpy = genesis::utils::indent( rcpy, std::string( left_w, ' ' ));
        rcpy = genesis::utils::trim( rcpy );
        out << rcpy;
    }
    out << "\n";
}

char translate_separator_char( CliOption<std::string> const& option )
{
    if(
        option.value == "comma" ||
        option.value == ","
    ) {
        return ',';
    } else if(
        option.value == "tab" ||
        option.value == "tabulator" ||
        option.value == "\t"
    ) {
        return '\t';
    } else if(
        option.value == "space" ||
        option.value == " "
    ) {
        return ' ';
    } else if(
        option.value == "semicolon" ||
        option.value == ";"
    ) {
        return ';';
    } else {
        throw CLI::ValidationError(
            option.option->get_name(),
            "Invalid separator char '" + option.value + "'."
        );
    }
}
