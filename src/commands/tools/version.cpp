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

#include "commands/tools/version.hpp"
#include "options/global.hpp"
#include "tools/references.hpp"
#include "tools/version.hpp"

#include "genesis/utils/core/options.hpp"

// =================================================================================================
//      Setup
// =================================================================================================

void setup_version( CLI::App& app )
{
    // Good place to check citations. This is executed every time with grenedalf,
    // so we never miss the check when editing the citation list.
    check_all_citations();

    // Create the options and subcommand objects.
    auto options = std::make_shared<VersionOptions>();
    auto sub = app.add_subcommand(
        "version",
        "Extended version information about grenedalf."
    );

    // Set the run function as callback to be called when this subcommand is issued.
    // Hand over the options by copy, so that their shared ptr stays alive in the lambda.
    sub->callback( [ options ]() {
        run_version( *options );
    });
}

// =================================================================================================
//      Run
// =================================================================================================

void run_version( VersionOptions const& options )
{
    (void) options;

    LOG_BOLD << grenedalf_header();
    LOG_BOLD;
    LOG_BOLD << "grenedalf version: " << grenedalf_version();
    LOG_BOLD;
    LOG_BOLD << genesis::utils::Options::get().info_compile_time();
    LOG_BOLD;
    LOG_BOLD << "For citation information, call  `grenedalf tools citation`";
    LOG_BOLD << "For license information, call  `grenedalf tools license`";
    LOG_BOLD;
    LOG_BOLD << grenedalf_title();
}
