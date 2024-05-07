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

#include "commands/convert/sync.hpp"
#include "options/global.hpp"
#include "tools/cli_setup.hpp"

#include "genesis/population/function/functions.hpp"
#include "genesis/population/format/sync_common.hpp"
#include "genesis/population/stream/variant_input_stream_adapters.hpp"

// =================================================================================================
//      Setup
// =================================================================================================

void setup_sync( CLI::App& app )
{
    // Create the options and subcommand objects.
    auto options = std::make_shared<SyncOptions>();
    auto sub = app.add_subcommand(
        "sync",
        "Create a sync file that lists per-sample base counts at each position in the genome."
    );

    // -------------------------------------------------------------------------
    //     Input
    // -------------------------------------------------------------------------

    // Required input of some frequency format (mpileup or vcf at the moment).
    options->variant_input.add_variant_input_opts_to_app( sub );

    // -------------------------------------------------------------------------
    //     Settings
    // -------------------------------------------------------------------------

    // Without header
    options->without_header.option = sub->add_flag(
        "--no-header",
        options->without_header.value,
        "We provide an extension of the sync format that allows to store sample names in "
        "sync files, where a header line is added to the output file of the form: "
        "`#chr pos ref S1...`, where `S1...` is the list of sample names. "
        "Not all other tools that read sync files will be able to parse this, "
        "and it hence can be deactivated with this option."
    );
    options->without_header.option->group( "Settings" );

    // Without missing
    options->without_missing.option = sub->add_flag(
        "--no-missing-marker",
        options->without_missing.value,
        "We provide an extension of the sync format that allows to mark positions as missing, "
        "by using the format `.:.:.:.:.:.`, instead of setting zero counts for these positions. "
        "This is particularly useful when storing data from multiple samples, to for instance "
        "distinguish true missing data from positions that have been filtered out. "
        "Not all other tools that read sync files will be able to parse this, and it hence can be "
        "deactivated here, in which case zero counts are written at these positions instead."
    );
    options->without_missing.option->group( "Settings" );

    // Gapless
    options->gapless.option = sub->add_flag(
        "--gapless-gsync",
        options->gapless.value,
        "By default, only the positions for which there is data are printed in the output. "
        "However, it might make processing with other tools easier if all files contain all "
        "positions, which one might call a gsync file (following gvcf). With this option, "
        "all missing positions are filled with the missing data indicator or with zero counts, "
        "depending on the `--no-missing-marker` option."
    );
    options->gapless.option->group( "Settings" );

    // Guess Reference Base
    options->guess_ref_base.option = sub->add_flag(
        "--guess-reference-base",
        options->guess_ref_base.value,
        "By default, when reading from input file formats that do not store the reference base, "
        "we do not attempt to guess it. When set however, we use the base with the highest count "
        "as the reference base for the output. Alternatively, when a reference genome is provided, "
        "we use that to correctly set the reference bases, independently of whether this flag is set."
    );
    options->guess_ref_base.option->group( "Settings" );

    // // max_concurrent_files
    // options->max_concurrent_files.option = sub->add_option(
    //     "--max-concurrent-files",
    //     options->max_concurrent_files.value,
    //     "Typical Unix-like systems have a maximum number of files that can be opened simultaneously "
    //     "by a process (usually, ~1024, hence the default of 1000 for this setting). "
    //     "For cases where more files need to be processed than that, this setting can be used "
    //     "to avoid errors given by the operating system. This works by processing files in batches "
    //     "of the given size into temporary files, and then combining them in a second step."
    // );
    // options->max_concurrent_files.option->group( "Settings" );
    // options->max_concurrent_files.option->transform( CLI::Range( 2, 1000 ));
    //
    // // temp_dir
    // options->temp_dir.option = sub->add_option(
    //     "--temp-dir",
    //     options->temp_dir.value,
    //     "Temporary directory to write intermediate files to. This is only used if more input "
    //     "files are provided than the above `" + options->max_concurrent_files.option->get_name() +
    //     "` setting, so that the processing has to be done in batches. Files are deleted after "
    //     "finishing successfully. We use a hash of the current time and the process ID of the program "
    //     "as a prefix to avoid name collisions for multiple simultaneous invocations of this command."
    // );
    // options->temp_dir.option->group( "Settings" );

    // -------------------------------------------------------------------------
    //     Output
    // -------------------------------------------------------------------------

    // Output
    options->file_output.add_default_output_opts_to_app( sub );
    options->file_output.add_file_compress_opt_to_app( sub );

    // -------------------------------------------------------------------------
    //     Callback
    // -------------------------------------------------------------------------

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

    // Pre-check for file existence
    options.file_output.check_output_files_nonexistence( "sync", "sync" );

    // If we have a reference genome, use it to correctly get the bases for the output here.
    if( options.variant_input.get_reference_genome() ) {
        auto const ref_gen = options.variant_input.get_reference_genome();
        options.variant_input.add_combined_filter_and_transforms(
            [ref_gen]( genesis::population::Variant& var ){
                guess_and_set_ref_and_alt_bases( var, *ref_gen );
                return true;
            }
        );
    } else if( options.guess_ref_base.value ) {
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

    // If we want to create a gsync file, activate this in the input stream.
    // Needs to be set before any reading access to the stream.
    if( options.gapless.value ) {
        options.variant_input.gapless_stream( true );
    }

    // Open the file.
    LOG_MSG << "Writing to file " << options.file_output.get_output_filename( "sync", "sync", true );
    auto sync_ofs = options.file_output.get_output_target( "sync", "sync" );

    // Add the header line if needed.
    if( ! options.without_header.value ) {
        (*sync_ofs) << "#chr\tpos\tref";
        for( auto const& sample_name :  options.variant_input.sample_names() ) {
            (*sync_ofs) << "\t" << sample_name;
        }
        (*sync_ofs) << "\n";
    }

    // TODO implement merging here as well?! --> options are there already, just in case ;-)

    // Write the sync data
    for( auto const& freq_it : options.variant_input.get_stream() ) {
        to_sync( freq_it, sync_ofs->ostream(), ! options.without_missing.value );
    }
}
