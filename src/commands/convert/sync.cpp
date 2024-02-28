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
    Lucas Czech <lczech@carnegiescience.edu>
    Department of Plant Biology, Carnegie Institution For Science
    260 Panama Street, Stanford, CA 94305, USA
*/

#include "commands/convert/sync.hpp"
#include "options/global.hpp"
#include "tools/cli_setup.hpp"

#include "genesis/population/functions/functions.hpp"
#include "genesis/population/formats/sync_common.hpp"

// =================================================================================================
//      Setup
// =================================================================================================

void setup_sync( CLI::App& app )
{
    // Create the options and subcommand objects.
    auto options = std::make_shared<SyncOptions>();
    auto sub = app.add_subcommand(
        "sync",
        "Create a PoPoolation2 sync file that lists per-sample base counts at each position "
        "in the genome."
    );

    // Required input of some frequency format (mpileup or vcf at the moment).
    options->variant_input.add_variant_input_opts_to_app( sub );
    // options->variant_input.add_sample_name_opts_to_app( sub );
    // options->variant_input.add_region_filter_opts_to_app( sub );

    // With Header
    options->with_header.option = sub->add_flag(
        "--with-header",
        options->with_header.value,
        "We provide an ad-hoc extension of the sync format that allows to store sample names in "
        "sync files. When using this flag, a header line is added to the output file of the form: "
        "`#chr pos ref S1...`, where `S1...` is the list of sample names. Not all tools that read "
        "sync files will be able to parse this though."
    );
    options->with_header.option->group( "Settings" );

    // Guess Reference Base
    options->guess_reference_base.option = sub->add_flag(
        "--guess-reference-base",
        options->guess_reference_base.value,
        "By default, when reading from input file formats that do not store the reference base, "
        "we do not attempt to guess it. When set however, we use the base with the highest count "
        "as the reference base for the output. Alternatively, when a reference genome is provided, "
        "we use that to correctly set the reference bases, independently of whether this flag is set."
    );
    options->guess_reference_base.option->group( "Settings" );

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
            run_sync( *options );
        }
    ));
}

// =================================================================================================
//      Run
// =================================================================================================

void run_sync( SyncOptions const& options )
{
    using namespace genesis::population;

    options.file_output.check_output_files_nonexistence( "sync", "sync" );
    auto sync_ofs = options.file_output.get_output_target( "sync", "sync" );

    // TODO The below ref genome code is almost identical to the one in frequency.
    // Refactor to avoid code duplication.

    // If we have a reference genome, use it to correctly get the bases for the output here.
    if( options.variant_input.get_reference_genome() ) {
        auto const ref_gen = options.variant_input.get_reference_genome();
        options.variant_input.add_combined_filter_and_transforms(
            [ref_gen]( genesis::population::Variant& var ){
                guess_and_set_ref_and_alt_bases( var, *ref_gen );
                return true;
            }
        );
    } else if( options.guess_reference_base.value ) {
        // Most of our input sources do not provide ref, and almost non provide alt bases.
        // So we use our guess function to augment the data. The function is idempotent
        // (unless we set the `force` parameter, which we do not do here), so for sources that do
        // contain ref and/or alt bases, nothing changes.
        options.variant_input.add_combined_filter_and_transforms(
            []( genesis::population::Variant& var ){
                guess_and_set_ref_and_alt_bases( var );
                return true;
            }
        );
    }

    // Add the header line
    if( options.with_header.value ) {
        (*sync_ofs) << "#chr\tpos\tref";
        for( auto const& sample_name :  options.variant_input.sample_names() ) {
            (*sync_ofs) << "\t" << sample_name;
        }
        (*sync_ofs) << "\n";
    }

    // Write the sync data
    for( auto const& freq_it : options.variant_input.get_stream() ) {
        to_sync( freq_it, sync_ofs->ostream() );
    }
}
