#ifndef GRENEDALF_OPTIONS_TABLE_OUTPUT_H_
#define GRENEDALF_OPTIONS_TABLE_OUTPUT_H_

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

#include "CLI/CLI.hpp"

#include "options/file_output.hpp"
#include "tools/cli_option.hpp"

#include <string>
#include <vector>

// =================================================================================================
//      Poolsizes Options
// =================================================================================================

/**
 * @brief Output data in table form.
 *
 * We have several commands that need to output data in form of a table, and in particular,
 * a table that has several fields (types of values) for each input sample. This is a bit of
 * a tricky setup, because the user might want to name these fields in different ways, or even
 * have separate output files for each. This class helps to organize this.
 */
class TableOutputOptions
{
public:

    // -------------------------------------------------------------------------
    //     Constructor and Rule of Five
    // -------------------------------------------------------------------------

    TableOutputOptions()  = default;
    virtual ~TableOutputOptions() = default;

    TableOutputOptions( TableOutputOptions const& other ) = default;
    TableOutputOptions( TableOutputOptions&& )            = default;

    TableOutputOptions& operator= ( TableOutputOptions const& other ) = default;
    TableOutputOptions& operator= ( TableOutputOptions&& )            = default;

    // -------------------------------------------------------------------------
    //     Setup Functions
    // -------------------------------------------------------------------------

    void add_separator_char_opt_to_app(
        CLI::App* sub,
        std::string const& group = "Formatting"
    );

    void add_na_entry_opt_to_app(
        CLI::App* sub,
        std::string const& group = "Formatting"
    );

    // -------------------------------------------------------------------------
    //     Run Functions
    // -------------------------------------------------------------------------

    char get_separator_char() const;

    std::string const& get_na_entry() const
    {
        return na_entry_.value;
    }

    // -------------------------------------------------------------------------
    //     Option Members
    // -------------------------------------------------------------------------

private:

    // FileOutputOptions  file_output_;

    CliOption<std::string> separator_char_ = "tab";
    CliOption<std::string> na_entry_ = "NA";

};

#endif // include guard
