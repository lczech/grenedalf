#ifndef GRENEDALF_OPTIONS_WINDOW_AVERAGE_H_
#define GRENEDALF_OPTIONS_WINDOW_AVERAGE_H_

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

#include "tools/cli_option.hpp"

#include "genesis/population/function/window_average.hpp"

#include <string>
#include <vector>

// =================================================================================================
//      Window Average Options
// =================================================================================================

/**
 * @brief Window averaging options.
 *
 * We have several commands that want to compute averages of statistics across windows. There
 * are different ways of doing that, so we offer this option here to select.
 */
class WindowAverageOptions
{
public:

    // -------------------------------------------------------------------------
    //     Constructor and Rule of Five
    // -------------------------------------------------------------------------

    WindowAverageOptions()  = default;
    ~WindowAverageOptions() = default;

    WindowAverageOptions( WindowAverageOptions const& other ) = default;
    WindowAverageOptions( WindowAverageOptions&& )            = default;

    WindowAverageOptions& operator= ( WindowAverageOptions const& other ) = default;
    WindowAverageOptions& operator= ( WindowAverageOptions&& )            = default;

    // -------------------------------------------------------------------------
    //     Setup Functions
    // -------------------------------------------------------------------------

    CLI::Option* add_window_average_opt_to_app(
        CLI::App* sub,
        std::string const& group = "Settings"
    );

    // -------------------------------------------------------------------------
    //     Run Functions
    // -------------------------------------------------------------------------

    genesis::population::WindowAveragePolicy get_window_average_policy() const;

    // -------------------------------------------------------------------------
    //     Option Members
    // -------------------------------------------------------------------------

private:

    CliOption<std::string> window_average_;

};

#endif // include guard
