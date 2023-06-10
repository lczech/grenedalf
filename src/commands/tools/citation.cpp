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

#include "commands/tools/citation.hpp"
#include "options/global.hpp"
#include "tools/references.hpp"

// =================================================================================================
//      Setup
// =================================================================================================

void setup_citation( CLI::App& app )
{
    // Good place to check citations. This is executed every time with grenedalf,
    // so we never miss the check when editing the citation list.
    check_all_citations();

    // Create the options and subcommand objects.
    auto options = std::make_shared<CitationOptions>();
    auto sub = app.add_subcommand(
        "citation",
        "Print references to be cited when using grenedalf."
    );

    // Which keys to cite.
    auto keys_opt = sub->add_option(
        "keys",
        options->keys,
        "Only print the citations for the given keys."
    );
    keys_opt->transform( CLI::IsMember( get_all_citation_keys() ));

    // Format.
    auto format_opt = sub->add_option(
        "--format",
        options->format,
        "Output format for citations.",
        true
    );
    format_opt->transform( CLI::IsMember({ "bibtex", "markdown", "both" }, CLI::ignore_case ));

    // Flags
    sub->add_flag(
        "--all",
        options->all,
        "Print all relevant citations used by commands in grenedalf."
    );
    sub->add_flag(
        "--list",
        options->list,
        "List all available citation keys."
    );

    // Set the run function as callback to be called when this subcommand is issued.
    // Hand over the options by copy, so that their shared ptr stays alive in the lambda.
    sub->callback( [ options ]() {
        run_citation( *options );
    });
}

// =================================================================================================
//      Run
// =================================================================================================

void run_citation( CitationOptions const& options )
{
    // If the --list flag is given, simply list all citation keys.
    if( options.list ) {
        LOG_BOLD << "Available citation keys:";
        for( auto const& key : get_all_citation_keys() ) {
            LOG_BOLD << " - " << key;
        }
        return;
    }

    // If --all is set, output all references that we have.
    // If not, by default, just print the grenedalf reference itself.
    // If keys are given, use those.
    std::vector<std::string> list;
    if( options.keys.empty() ) {
        if( options.all ) {
            list = get_all_citation_keys();
        } else {
            list = { "Czech2023-grenedalf" };
        }
    } else {
        list = options.keys;
    }

    // Do the printing in all desired formats.
    if( options.format == "bibtex" ) {
        LOG_BOLD << cite_bibtex( list );
    } else if( options.format == "markdown" ) {
        LOG_BOLD << cite_markdown( list );
    } else if( options.format == "both" ) {
        LOG_BOLD << cite_bibtex( list );
        LOG_BOLD;
        LOG_BOLD << cite_markdown( list );
    } else {
        throw CLI::ValidationError( "--format", "Invalid citation format: " + options.format );
    }
}
