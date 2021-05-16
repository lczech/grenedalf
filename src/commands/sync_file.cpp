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

#include "commands/sync_file.hpp"
#include "options/global.hpp"
#include "tools/cli_setup.hpp"

#include "genesis/population/functions/base_counts.hpp"
#include "genesis/population/functions/variant.hpp"

// =================================================================================================
//      Setup
// =================================================================================================

void setup_sync_file( CLI::App& app )
{
    // Create the options and subcommand objects.
    auto options = std::make_shared<SyncFileOptions>();
    auto sub = app.add_subcommand(
        "sync-file",
        "Create a PoPoolation2 sync file that lists per-sample base counts at each position "
        "in the genome."
    );

    // Required input of some frequency format (mpileup or vcf at the moment).
    options->freq_input.add_frequency_input_opts_to_app( sub );
    options->freq_input.add_sample_name_opts_to_app( sub );
    options->freq_input.add_filter_opts_to_app( sub );

    // Output
    options->file_output.add_default_output_opts_to_app( sub );
    options->file_output.add_file_compress_opt_to_app( sub );

    // Set the run function as callback to be called when this subcommand is issued.
    // Hand over the options by copy, so that their shared ptr stays alive in the lambda.
    sub->callback( grenedalf_cli_callback(
        sub,
        {
            // Citation keys as needed
            "Kofler2011-PoPoolation2"
        },
        [ options ]() {
            run_sync_file( *options );
        }
    ));
}

// =================================================================================================
//      Run
// =================================================================================================

void run_sync_file( SyncFileOptions const& options )
{
    using namespace genesis::population;

    options.file_output.check_output_files_nonexistence( "counts", "sync" );
    auto sync_ofs = options.file_output.get_output_target( "counts", "sync" );

    // Write the sync data
    for( auto const& freq_it : options.freq_input.get_iterator() ) {
        to_sync( freq_it, sync_ofs->ostream() );
        (*sync_ofs) << "\n";
    }
}
