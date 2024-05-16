#ifndef GRENEDALF_COMMANDS_ANALYZE_FST_CATHEDRAL_H_
#define GRENEDALF_COMMANDS_ANALYZE_FST_CATHEDRAL_H_

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
    Lucas Czech <lucas.czech@sund.ku.dk>
    University of Copenhagen, Globe Institute, Section for GeoGenetics
    Oster Voldgade 5-7, 1350 Copenhagen K, Denmark
*/

#include "CLI/CLI.hpp"

#include "options/file_output.hpp"
#include "options/poolsizes.hpp"
#include "options/fst_processor.hpp"
#include "options/variant_filter_numerical.hpp"
#include "options/variant_input.hpp"
#include "options/window.hpp"
#include "tools/cli_option.hpp"

#include <string>
#include <vector>

// =================================================================================================
//      Options
// =================================================================================================

class FstCathedralOptions
{
public:

    // Input options
    VariantInputOptions           variant_input;
    VariantFilterNumericalOptions filter_numerical;

    // Specific settings
    FstProcessorOptions    fst_processor;
    CliOption<size_t>      cathedral_width  = 1500;
    CliOption<size_t>      cathedral_height = 500;
    CliOption<std::string> cathedral_method = "exponential";

    // Output options
    FileOutputOptions  file_output;

};

// =================================================================================================
//      Functions
// =================================================================================================

void setup_fst_cathedral( CLI::App& app );
void run_fst_cathedral( FstCathedralOptions const& options );

#endif // include guard
