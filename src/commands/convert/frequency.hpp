#ifndef GRENEDALF_COMMANDS_CONVERT_FREQUENCY_H_
#define GRENEDALF_COMMANDS_CONVERT_FREQUENCY_H_

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

#include "CLI/CLI.hpp"

#include "options/file_output.hpp"
#include "options/variant_input.hpp"
#include "options/table_output.hpp"
#include "tools/cli_option.hpp"

#include <string>
#include <vector>

// =================================================================================================
//      Options
// =================================================================================================

class FrequencyOptions
{
public:

    CliOption<bool> write_sample_counts     = false;
    CliOption<bool> write_sample_read_depth = false;
    CliOption<bool> write_sample_ref_freq   = false;
    CliOption<bool> write_sample_alt_freq   = false;

    CliOption<bool> write_total_counts      = false;
    CliOption<bool> write_total_read_depth  = false;
    CliOption<bool> write_total_frequency   = false;

    CliOption<bool> write_invariants        = false;
    CliOption<bool> omit_ref_alt_bases      = false;
    CliOption<bool> omit_alt_bases          = false;

    VariantInputOptions variant_input;
    TableOutputOptions  table_output;
    FileOutputOptions   file_output;

};

// =================================================================================================
//      Functions
// =================================================================================================

void setup_frequency( CLI::App& app );
void run_frequency( FrequencyOptions const& options );

#endif // include guard
