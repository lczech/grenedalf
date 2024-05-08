#ifndef GRENEDALF_COMMANDS_ANALYZE_CATHEDRAL_PLOT_H_
#define GRENEDALF_COMMANDS_ANALYZE_CATHEDRAL_PLOT_H_

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

#include "options/file_input.hpp"
#include "options/file_output.hpp"
#include "options/color_map.hpp"
#include "tools/cli_option.hpp"

#include "genesis/population/plotting/cathedral_plot.hpp"

#include <limits>
#include <string>
#include <vector>

// =================================================================================================
//      Options
// =================================================================================================

class CathedralPlotOptions
{
public:

    // Input options
    FileInputOptions json_input;
    FileInputOptions csv_input;

    // Plot settings
    ColorMapOptions color_map;
    CliOption<std::string> color_norm = "linear";
    CliOption<std::string> color_norm_range = "all";
    CliOption<double> min_value = std::numeric_limits<double>::quiet_NaN();
    CliOption<double> max_value = std::numeric_limits<double>::quiet_NaN();

    // Output options
    FileOutputOptions  file_output;

};

// =================================================================================================
//      Functions
// =================================================================================================

void setup_cathedral_plot( CLI::App& app );
void run_cathedral_plot( CathedralPlotOptions const& options );

#endif // include guard
