#ifndef GRENEDALF_COMMANDS_FST_H_
#define GRENEDALF_COMMANDS_FST_H_

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
#include "options/frequency_input.hpp"
#include "options/poolsizes.hpp"
#include "tools/cli_option.hpp"

#include <string>
#include <vector>

// =================================================================================================
//      Options
// =================================================================================================

class FstOptions
{
public:

    FrequencyInputOptions freq_input;
    PoolsizesOptions poolsizes;

    CliOption<std::string> method = "conventional";
    CliOption<bool>        omit_empty_windows = false;
    CliOption<std::string> comparand = "";
    CliOption<std::string> second_comparand = "";
    CliOption<std::string> comparand_list = "";
    CliOption<std::string> na_entry = "NA";

    FileOutputOptions  file_output;

};

// =================================================================================================
//      Functions
// =================================================================================================

void setup_fst( CLI::App& app );
void run_fst( FstOptions const& options );

#endif // include guard
