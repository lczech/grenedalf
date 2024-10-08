#ifndef GRENEDALF_COMMANDS_ANALYZE_DIVERSITY_H_
#define GRENEDALF_COMMANDS_ANALYZE_DIVERSITY_H_

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

#include "options/diversity_processor.hpp"
#include "options/file_output.hpp"
#include "options/poolsizes.hpp"
#include "options/table_output.hpp"
#include "options/variant_filter_numerical.hpp"
#include "options/variant_input.hpp"
#include "options/variant_transform_subsample.hpp"
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
    VariantInputOptions              variant_input;
    VariantFilterNumericalOptions    filter_numerical;
    VariantTransformSubsampleOptions transform_subsample;
    WindowOptions                    window;

    // Specific settings
    DiversityProcessorOptions diversity_processor;
    CliOption<bool> no_extra_columns   = false;
    CliOption<bool> no_nan_windows     = false;
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
