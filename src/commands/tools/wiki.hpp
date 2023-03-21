#ifndef GRENEDALF_COMMANDS_TOOLS_WIKI_H_
#define GRENEDALF_COMMANDS_TOOLS_WIKI_H_

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

#include <string>

// =================================================================================================
//      Options
// =================================================================================================

class WikiOptions
{
public:

    std::string md_dir  = "doc/md/";
    std::string out_dir = "../grenedalf.wiki/";
    bool use_details = false;

    CLI::App*   app;

};

// =================================================================================================
//      Functions
// =================================================================================================

void setup_wiki( CLI::App& app );
void run_wiki( WikiOptions const& options );

#endif // include guard
