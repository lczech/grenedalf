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

#include "commands/tools/wiki.hpp"
#include "options/global.hpp"
#include "tools/cli_setup.hpp"
#include "tools/references.hpp"
#include "tools/version.hpp"

#include "genesis/utils/core/fs.hpp"
#include "genesis/utils/text/char.hpp"
#include "genesis/utils/text/string.hpp"

#include <algorithm>
#include <iostream>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

// =================================================================================================
//      Setup
// =================================================================================================

void setup_wiki( CLI::App& app )
{
    // Create the options and subcommand objects.
    auto options = std::make_shared<WikiOptions>();
    auto sub = app.add_subcommand(
        "wiki",
        "Create wiki help pages for grenedalf."
    )->group( "" );

    // Need to capture the main app, as the wiki needs this to run.
    options->app = &app;
    while( options->app->get_parent() ) {
        options->app = options->app->get_parent();
    }

    // Markdown dir option.
    auto md_dir_opt = sub->add_option(
        "--md-dir",
        options->md_dir,
        "Directory with the Markdown files documenting the grenedalf commands."
    );
    md_dir_opt->group( "Settings" );
    md_dir_opt->check( CLI::ExistingDirectory );
    // md_dir_opt->required();

    // Out dir option.
    auto out_dir_opt = sub->add_option(
        "--out-dir",
        options->out_dir,
        "Directory to write Wiki files to. Should be a git clone of the wiki repository."
    );
    out_dir_opt->group( "Settings" );
    out_dir_opt->check( CLI::ExistingDirectory );
    // out_dir_opt->required();

    // Set the run function as callback to be called when this subcommand is issued.
    // Hand over the options by copy, so that their shared ptr stays alive in the lambda.
    sub->callback( [ options ]() {
        run_wiki( *options );
    });
}

// =================================================================================================
//      Helper Functions
// =================================================================================================

// -------------------------------------------------------------------------
//     Markdown Helpers
// -------------------------------------------------------------------------

/**
 * @brief Take a markdown text and replace all pairs of backticks by html code environments.
 */
std::string codify_markdown( std::string const& text )
{
    // Replace pairs of backticks.
    auto result = text;
    size_t pos = 0;
    size_t cnt = 0;
    while( pos < result.size() ) {
        if( result[pos] == '`' ) {
            result.replace( pos, 1, cnt % 2 == 0 ? "<code>" : "</code>" );
            ++cnt;
        }
        ++pos;
    }

    // We are calling this function only internally, and independent of user input,
    // so we have control over what arguments we call it with. Let's hence be on the safe side,
    // and only replace stuff that we know is good.
    if( cnt % 2 != 0 ) {
        throw std::runtime_error( "Invalid markdown with uneven number of backticks." );
    }

    return result;
}

// -------------------------------------------------------------------------
//     App Subcommand Helpers
// -------------------------------------------------------------------------

/**
 * @brief Get the immediate subcommands of an App, sorted by name.
 */
std::vector<CLI::App const*> get_sorted_subcommands( CLI::App const* app )
{
    std::vector<CLI::App const*> subcomms;
    for( auto const& subcom : app->get_subcommands({}) ) {
        if( subcom->get_group() != "" ) {
            subcomms.push_back( subcom );
        }
    }

    std::sort(
        subcomms.begin(), subcomms.end(),
        []( CLI::App const* lhs, CLI::App const* rhs ){
            return lhs->get_name() < rhs->get_name();
        }
    );
    return subcomms;
}

/**
 * @brief Get all subcommands, recursively, of an App, sorted by name.
 */
std::vector<CLI::App const*> get_all_subcommands( CLI::App const* app )
{
    std::vector<CLI::App const*> result;

    // Fill with subcommands.
    auto subcomms = get_sorted_subcommands( app );
    for( auto const& sc : subcomms ) {
        result.push_back( sc );

        // Recurse. The copying is wasteful, but it's not that much, and I am lazy today.
        auto subsubcomms = get_all_subcommands( sc );
        for( auto const& ssc : subsubcomms ) {
            result.push_back( ssc );
        }
    }

    return result;
}

/**
 * @brief Add the contents of a file to a stream.
 */
bool add_markdown_content(
    WikiOptions const& options, std::string const& md_file, std::ostream& os, bool add_header = false
) {
    using namespace genesis::utils;

    // Add markdown file content.
    std::string const fn = dir_normalize_path( options.md_dir ) + md_file + ".md";
    if( file_is_readable( fn ) ) {
        std::ifstream mds( fn );
        if( add_header ) {
            os << "# Description\n\n";
        }
        os << mds.rdbuf();
        return true;
    } else {
        LOG_MSG << " - No documentation markdown found: " << md_file;
    }
    return false;
}

/**
 * @brief Check if a command has any options at all, and any required options.
 */
std::pair<bool, bool> command_has_options( CLI::App const& command )
{
    // We do not count the help option, so we need to manually check if there are any others.
    bool has_options = false;
    bool has_required_options = false;
    for( auto const& opt : command.get_options() ) {
        if( opt->get_name() != "-h,--help" && opt->get_name() != "--help" ) {
            has_options = true;
        }
        if( opt->get_required() ) {
            has_required_options = true;
        }
    }
    return std::make_pair( has_options, has_required_options );
}

// -------------------------------------------------------------------------
//     Add Command Header Synopsis
// -------------------------------------------------------------------------

void add_wiki_command_header_synopsis(
    CLI::App const& command, std::ostream& os
) {
    // Get the usage line.
    std::string usage = command.get_name();
    auto parent = const_cast< CLI::App& >( command ).get_parent();
    while( parent ) {
        usage  = parent->get_name() + " " + usage;
        parent = parent->get_parent();
    }

    // Write command header.
    os << "**Synopsis:** " << command.get_description() << "\n\n";
    os << "**Usage:** `" << usage;
    if( command_has_options( command ).first ) {
        os << " [options]";
    }
    if( ! command.get_subcommands({}).empty() ) {
        if( command.get_require_subcommand_min() > 0 ) {
            os << " subcommand";
        } else {
            os << " [subcommand]";
        }
    }
    os << "`\n\n";
    if( command_has_options( command ).second ) {
        os << "**Required options:**\n\n";
        for( auto const& opt : command.get_options() ) {
            if( opt->get_required() ) {
                os << "  * `" << opt->get_name() << "`\n";
            }
        }
        os << "\n";
    }
    os << "Documentation for grenedalf " << grenedalf_version() << "\n\n";
}

// -------------------------------------------------------------------------
//     Add Command Header TOC
// -------------------------------------------------------------------------

void add_wiki_command_header_toc(
    std::string const& md_file
) {
    // Get the markdown for the command, and find all its headings.
    using namespace genesis::utils;

    std::stringstream os;
    os << "**Table of contents:**\n\n";

    // Helper to get the link for a given markdown header
    auto get_link_from_header_ = []( std::string const& header )
    {
        auto const link = to_lower( replace_all( header, " ", "-" ));
        return remove_all_chars_pred(
            link,
            []( char c ){
                return !(is_alpha(c) || is_digit(c) || c == '-');
            }
        );
    };

    // Add markdown file content headers.
    auto const md_lines = file_read_lines( md_file );
    for( auto const& line : md_lines ) {
        std::string suffix;

        // Heading levels
        if( starts_with( line, "# ", suffix )) {
            os << "  * [" << suffix << "](#" << get_link_from_header_( suffix ) << ")\n";
        }
        if( starts_with( line, "## ", suffix )) {
            os << "      * [" << suffix << "](#" << get_link_from_header_( suffix ) << ")\n";
        }
        if( starts_with( line, "### ", suffix )) {
            os << "          * [" << suffix << "](#" << get_link_from_header_( suffix ) << ")\n";
        }
    }

    // Now read the file as a string again, and replace the TOC.
    auto md_content = file_read( md_file );
    md_content = replace_all( md_content, "<!--TOC-->", os.str() );

    // We do not want to use the normal file write function here,
    // as it throws when the file already exists, which it always does here.
    // file_write( md_content, md_file );
    std::ofstream of( md_file );
    of << md_content;
    of.close();
}

// -------------------------------------------------------------------------
//     Make Command Header
// -------------------------------------------------------------------------

void make_wiki_command_header(
    CLI::App const& command, std::ostream& os
) {
    // Add the table header for the two-column layout. Best we can do in GitHub flavored markdown.
    os << "<table border=\"0\">\n";
    os << "<tr>\n";
    os << "<td width=\"441\" valign=\"top\">\n\n";

    add_wiki_command_header_synopsis( command, os );

    // Now add the html for starting the second column.
    os << "</td>\n";
    os << "<td width=\"441\" valign=\"top\">\n\n";

    // add_wiki_command_header_toc( wiki_options, command, os );
    os << "<!--TOC-->\n";

    // Finally, wrap up the two column table.
    os << "</td>\n";
    os << "</tr>\n";
    os << "</table>\n\n";
}

// -------------------------------------------------------------------------
//     Make Options Table
// -------------------------------------------------------------------------

void make_options_table( WikiOptions const& wiki_options, CLI::App const& command, std::ostream& os )
{
    // Get the options that are part of this command.
    auto const options = command.get_options();

    // map from group name to table contents, with the tuple
    // (1) group name, (2) contents, (3) is any option required
    using OptGroup = std::tuple<std::string, std::string, bool>;
    // we use a vec to keep order.
    // std::map<std::string, std::string> opt_helps;
    std::vector<OptGroup> opt_helps;

    // Add lines for each group.
    for( auto const& opt : options ) {

        // Do not add help option.
        if( opt->get_name() == "-h,--help" || opt->get_name() == "--help" ) {
            continue;
        }

        // Do not add hidden options.
        if( opt->get_group() == "" ) {
            continue;
        }

        // Write to temporary stream.
        std::stringstream tmp_os;

        // Simple version that makes a markdown table.
        // tmp_os << "| <nobr>`" << opt->get_name() << "`</nobr> ";
        // tmp_os << "|";
        // if( opt->get_required() ) {
        //     tmp_os << " **Required.**";
        // }
        // if( ! opt->help_aftername().empty() ) {
        //     // print stuff without leading space.
        //     tmp_os << " `" << opt->help_aftername().substr( 1 ) << "`<br />";
        // }
        //
        // auto descr = opt->get_description();
        // tmp_os << " " << descr << " |\n";
        // // tmp_os << " " << opt->get_description() << " |\n";
        // // tmp_os << "| " << opt->get_description() << " |\n";

        // Add content to the group help.
        // tmp_os << "<tr><td><code>" << opt->get_name() << "</code></td>";
        // tmp_os << "<td>";
        tmp_os << "<dt><code>" << opt->get_name() << "</code></dt>";

        // Get option type, special flags and validator strings.
        auto formatter = dynamic_cast<CLI::Formatter const*>( command.get_formatter().get() );
        auto opt_str = formatter->make_option_opts( opt );
        opt_str = genesis::utils::replace_all(
            opt_str,
            command.get_formatter()->get_label("REQUIRED"),
            ""
        );
        if( opt->get_type_name().empty() ) {
            opt_str = "FLAG " + opt_str;
        }

        // Little special case: --threads defaults to the number of cores on the current system
        // where this wiki command is being run. Make this nicer.
        if( opt->get_name() == "--threads" ) {
            std::string tn;
            if( !opt->get_type_name().empty() ) {
                tn = formatter->get_label( opt->get_type_name() );
            }
            std::string search;
            if( !opt->get_default_str().empty() ) {
                search = tn + "=" + opt->get_default_str();
            }
            if( !search.empty() ) {
                opt_str = genesis::utils::replace_all( opt_str, search, tn );
            }
        }

        // Now print to the output.
        tmp_os << "<dd>";
        if( ! opt_str.empty() ) {
            tmp_os << "<code>" << genesis::utils::trim( opt_str ) << "</code><br />";
        }
        // tmp_os << "</dt>\n";

        // Add description. Replace backticks by html code elements here, as markdown
        // does not replace them automatically within html environments.
        // tmp_os << "<dd>";
        if( opt->get_required() ) {
            tmp_os << "<strong>Required.</strong> ";
        }
        auto descr = opt->get_description();
        if( descr.substr( 0, 10 ) == "Required. " ) {
            descr = descr.substr( 10 );
        }
        descr = genesis::utils::replace_all( descr, "\n", "<br />" );
        tmp_os << codify_markdown( descr ) << "</dd>\n";
        // tmp_os << " " << codify_markdown( descr ) << "</td></tr>\n";
        // tmp_os << " " << opt->get_description() << " |\n";
        // tmp_os << "| " << opt->get_description() << " |\n";

        // Add content to the group help.
        // first check if the group was already used, and if not add it.
        auto get_group = [&]( std::string const& name ) -> OptGroup& {
            for( auto& elem : opt_helps ) {
                if( std::get<0>(elem) == name ) {
                    return elem;
                }
            }
            opt_helps.push_back({ name, "", false });
            return opt_helps.back();
        };
        auto& opt_group = get_group( opt->get_group() );
        std::get<1>( opt_group ) += tmp_os.str();
        if( opt->get_required() ) {
            std::get<2>( opt_group ) = true;
        }
        // opt_helps[ opt->get_group() ] += tmp_os.str();
    }

    // Simple markdown verison to print the groups and their tables.
    // for( auto const& gr : opt_helps ) {
    //     os << "**" << gr.first << ":**\n\n";
    //     os << "| Option  | Description |\n";
    //     os << "| ------- | ----------- |\n";
    //     // os << "| Option  | Type | Description |\n";
    //     // os << "| ------- | ---- | ----------- |\n";
    //     os << gr.second << "\n";
    // }

    // Print the groups and their tables
    // os << "<table>\n";
    // bool done_first_group = false;
    // for( auto const& gr : opt_helps ) {
    //     if( done_first_group ) {
    //         os << "<tr height=30px></tr>\n";
    //     }
    //     os << "<thead><tr><th colspan=\"2\" align=\"left\">" << gr.first << "</th></tr></thead>\n";
    //     os << "<tbody>\n";
    //     os << gr.second;
    //     os << "</tbody>\n";
    //     done_first_group = true;
    // }
    // os << "</table>\n\n";

    // Print description lists
    for( auto const& gr : opt_helps ) {
        // os << "## " << gr.first << "\n\n";
        // os << "<dl>\n";
        // os << gr.second;
        // os << "</dl>\n\n";

        if( wiki_options.use_details ) {
            os << "<details" << ( std::get<2>(gr) ? " open" : "" ) << ">\n";
            os << "<summary>" << std::get<0>(gr) << "</summary>\n";
            os << "<dl>\n";
            os << std::get<1>(gr);
            os << "</dl>\n";
            os << "</details>\n\n";
        } else {
            os << "## " << std::get<0>(gr) << "\n\n";
            os << "<dl>\n";
            os << std::get<1>(gr);
            os << "</dl>\n\n";
        }
    }
}

// -------------------------------------------------------------------------
//     Make Subcommands Table
// -------------------------------------------------------------------------

void make_subcommands_table(
    std::vector<CLI::App const*> subcomms,
    std::ostream& os,
    std::string const& heading = "Subcommand"
) {
    os << "| " << heading << "  | Description |\n";
    os << "| ----------- | ----------- |\n";

    for( auto const& subcomm : subcomms ) {
        os << "| [" << subcomm->get_name() << "](../wiki/Subcommand:-" << subcomm->get_name() << ") ";
        os << "| " << subcomm->get_description() << " |\n";
    }
    os << "\n";
}

// -------------------------------------------------------------------------
//     Make Wiki Page
// -------------------------------------------------------------------------

void make_wiki_command_page( WikiOptions const& wiki_options, CLI::App const& command )
{
    using namespace genesis::utils;

    // User output.
    LOG_MSG << "Subcommand: " << command.get_name();

    // Open out file stream.
    std::string const out_file
        = dir_normalize_path( wiki_options.out_dir )
        + "Subcommand:-" + command.get_name() + ".md"
    ;
    if( ! file_is_readable( out_file )) {
        LOG_MSG << " - No existing wiki file!";
    }
    std::ofstream os( out_file );

    // Write styles. --> Nope, not accepted by GitHub wiki...
    // os << "<style>\n";
    // os << "code { white-space: nowrap; }\n";
    // os << "dt {\n";
    // os << "    font-weight: normal;\n";
    // os << "    font-style: normal;\n";
    // // os << "    margin-top: 8px;\n";
    // // os << "    margin-bottom: 8px;\n";
    // os << "}\n";
    // // os << "dd {\n";
    // // os << "    margin-left: 18px;\n";
    // // os << "}\n";
    // if( wiki_options.use_details ) {
    //     os << "details > summary {\n";
    //     os << "    font-size: 1.5em;\n";
    //     os << "    font-weight: bold;\n";
    //     os << "    line-height: 1.2;\n";
    //     os << "    margin-top: 1em;\n";
    //     os << "    margin-bottom: 16px;\n";
    //     // os << "    border-bottom: 1px solid #aaa;\n";
    //     os << "}\n";
    // }
    // os << "</style>\n\n";

    // Add a header, with synopsis, required options, and a table of contents.
    make_wiki_command_header( command, os );

    // Add markdown file content that contains the actual documentation.
    add_markdown_content( wiki_options, command.get_name(), os );

    // Print the options of the command.
    if( command_has_options( command ).first ) {
        os << "# Options\n\n";
        make_options_table( wiki_options, command, os );
    }

    // Print the subcommands of this command.
    auto const subcomms = command.get_subcommands({});
    if( ! subcomms.empty() ) {
        os << "# Subommands\n\n";
        make_subcommands_table( subcomms, os );
    }

    // If there is a citation list for this command, add it in a nice format.
    if( citation_list.count( &command ) > 0 ) {
        os << "\n";
        os << "# Citation\n\n";
        os << "When using this method, please do not forget to cite\n\n";
        os << cite_markdown( citation_list[ &command ], true, false );
    }

    // Now close the file, and parse it again to add the TOC.
    os.close();
    add_wiki_command_header_toc( out_file );
}

// -------------------------------------------------------------------------
//     Make Wiki Home Page
// -------------------------------------------------------------------------

void make_wiki_home_page( WikiOptions const& options )
{
    using namespace genesis::utils;

    // Make Home page.
    LOG_MSG << "Home";

    // Open stream
    std::string const out_file = dir_normalize_path( options.out_dir ) + "Home.md";
    if( ! file_is_readable( out_file )) {
        LOG_MSG << " - No existing wiki file!";
    }
    std::ofstream os( out_file );

    // Add home header.
    add_markdown_content( options, "Home_header", os );
    os << "\n";

    // Get all top level commands.
    auto subcomms = get_sorted_subcommands( options.app );

    // Make a list of all top level commands that do not have sub commands, that is,
    // which are not modules. We make this first, so that we get a nice order.
    std::vector<CLI::App const*> direct_subcommands;
    for( auto const& subcomm : subcomms ) {
        if( subcomm->get_group() == "" ) {
            continue;
        }

        // Get all subcommands of the comamnd, and if there are none, this is a top level command,
        // and not a module, so then we list it here.
        auto subsubcomms = get_sorted_subcommands( subcomm );
        if( subsubcomms.empty() ) {
            direct_subcommands.push_back( subcomm );
        }
    }

    // If we found direct subcommands of the main, make a commands table of them.
    if( ! direct_subcommands.empty() ) {
        os << "# Commands\n\n";
        make_subcommands_table( direct_subcommands, os, "Command" );
    }

    // Go through the list again, and add submodule lists for all subcommands
    // that have child commands of their own.
    for( auto const& subcomm : subcomms ) {
        if( subcomm->get_group() == "" ) {
            continue;
        }

        // Get all subcommands of the subcommand, and if there are any, make a list of them.
        // If not, this is a top level command, and not a module, so skip it.
        auto subsubcomms = get_sorted_subcommands( subcomm );
        if( !subsubcomms.empty() ) {
            os << "<br />\n";
            os << "### Module `" << subcomm->get_name() << "`\n\n";
            os << subcomm->get_description() << "\n\n";
            make_subcommands_table( subsubcomms, os );
        }
    }

    // Add home footer.
    add_markdown_content( options, "Home_footer", os );
    os.close();
}

// -------------------------------------------------------------------------
//     Side Bar
// -------------------------------------------------------------------------

void make_wiki_sidebar( WikiOptions const& options )
{
    using namespace genesis::utils;

    // Make Sidebar page.
    LOG_MSG << "Sidebar";

    // Open stream
    std::string const out_file = dir_normalize_path( options.out_dir ) + "_Sidebar.md";
    if( ! file_is_readable( out_file )) {
        LOG_MSG << " - No existing wiki file!";
    }
    std::ofstream os( out_file );

    // Add standard entries
    // os << "[Home](../wiki)\n\n";
    // os << "[Build](../wiki/Build)\n\n";

    // Add sidebar header, which includes the standard entries.
    bool const found_sidebar_header = add_markdown_content( options, "_Sidebar_header", os );
    if( found_sidebar_header ) {
        os << "\n";
    }

    // Get all top level commands.
    auto subcomms = get_sorted_subcommands( options.app );

    // Make a list of all top level commands that do not have sub commands, that is,
    // which are not modules. We make this first, so that we get a nice order.
    bool has_direct_subcommands = false;
    for( auto const& subcomm : subcomms ) {
        if( subcomm->get_group() == "" ) {
            continue;
        }

        // Get all subcommands of the comamnd, and if there are none, this is a top level command,
        // and not a module, so then we list it here.
        auto subsubcomms = get_sorted_subcommands( subcomm );
        if( subsubcomms.empty() ) {

            // First time we find a command that we want to print, we print a heading as well.
            if( ! has_direct_subcommands ) {
                os << "Commands\n\n";
                has_direct_subcommands = true;
            }

            os << " * [" << subcomm->get_name() << "](../wiki/Subcommand:-";
            os << subcomm->get_name() << ")\n";
        }
    }

    // After finishing, if we printed any direct commands,
    // we need some space before moving on to the modules.
    if( has_direct_subcommands ) {
        os << "\n\n";
    }

    // Go through the list again, and add submodule lists for all subcommands
    // that have child commands of their own.
    for( auto const& subcomm : subcomms ) {
        if( subcomm->get_group() == "" ) {
            continue;
        }

        // Get all subcommands of the subcommand, and if there are any, make a list of them.
        // If not, this is a top level command, and not a module, so skip it.
        auto subsubcomms = get_sorted_subcommands( subcomm );
        if( !subsubcomms.empty() ) {
            os << "Module `" << subcomm->get_name() << "`\n\n";
            for( auto const& subcomm : subsubcomms ) {
                os << " * [" << subcomm->get_name() << "](../wiki/Subcommand:-";
                os << subcomm->get_name() << ")\n";
            }
            os << "\n";
        }
    }

    os.close();
}

// =================================================================================================
//      Run
// =================================================================================================

void run_wiki( WikiOptions const& options )
{
    // Home page.
    make_wiki_home_page( options );
    make_wiki_sidebar( options );

    // Now, make pages for the commands of the main app and the modules.
    // This hence works both if we have sub-commands in modules, and if all commands are directly
    // commands below main.
    auto subcomms = get_sorted_subcommands( options.app );
    for( auto const& sc : subcomms ) {
        make_wiki_command_page( options, *sc );
        for( auto const& ssc : get_sorted_subcommands( sc )) {
            make_wiki_command_page( options, *ssc );
        }
    }
}
