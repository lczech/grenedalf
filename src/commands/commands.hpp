#ifndef GRENEDALF_COMMANDS_COMMANDS_H_
#define GRENEDALF_COMMANDS_COMMANDS_H_

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

#include "commands/afs_heatmap.hpp"
#include "commands/frequency.hpp"
#include "commands/fst.hpp"
#include "commands/sync_file.hpp"

#include "options/global.hpp"
#include "tools/cli_setup.hpp"

#include <string>
#include <vector>

// =================================================================================================
//      Functions
// =================================================================================================

inline void setup_commands( CLI::App& app )
{
    // Create the module subcommand objects.
    // auto sub = app.add_subcommand(
    //     "tools",
    //     "Auxiliary commands of grenedalf."
    // );
    // sub->require_subcommand( 1 );

    // Add module subcommands.
    setup_afs_heatmap( app );
    setup_frequency( app );
    setup_fst( app );
    setup_sync_file( app );

    // Add the global options to each of the above subcommands.
    global_options.add_to_module( app );
    set_module_help_group( app );
}

#endif // include guard
