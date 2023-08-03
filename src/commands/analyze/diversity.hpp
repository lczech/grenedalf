#ifndef GRENEDALF_COMMANDS_ANALYZE_DIVERSITY_H_
#define GRENEDALF_COMMANDS_ANALYZE_DIVERSITY_H_

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
#include "options/poolsizes.hpp"
#include "options/table_output.hpp"
#include "options/variant_filter_numerical.hpp"
#include "options/variant_input.hpp"
#include "options/window.hpp"
#include "tools/cli_option.hpp"

#include <string>
#include <vector>

// =================================================================================================
//      Options
// =================================================================================================

class DiversityOptions
{
public:

    // Input options
    VariantInputOptions           variant_input;
    VariantFilterNumericalOptions filter_numerical;
    WindowOptions                 window;
    PoolsizesOptions              poolsizes;

    // Using defaults from PoPoolation, see Variance-sliding.pl
    CliOption<std::string> measure = "all";
    // CliOption<size_t> min_count = 2;
    // CliOption<size_t> min_coverage = 4;
    // CliOption<size_t> max_coverage = 1000000;
    // CliOption<double> min_coverage_fraction = 0.6;

    // Compatibility
    CliOption<bool> with_popoolation_bugs = false;
    CliOption<bool> popoolation_format = false;

    // Output options
    TableOutputOptions table_output;
    FileOutputOptions  file_output;

};

// =================================================================================================
//      Functions
// =================================================================================================

void setup_diversity( CLI::App& app );
void run_diversity( DiversityOptions const& options );

#endif // include guard
