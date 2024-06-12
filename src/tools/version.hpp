#ifndef GRENEDALF_TOOLS_VERSION_H_
#define GRENEDALF_TOOLS_VERSION_H_

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

#include <string>

// =================================================================================================
//      Grenedalf Version
// =================================================================================================

inline std::string grenedalf_version()
{
    return "v0.5.1"; // #GRENEDALF_VERSION#
}

inline std::string grenedalf_header()
{
    return "\
                                                 ,--;                           \n\
                                                 \\  \\          ----.            \n\
    ,-----.           ,----.;---.  ,--.  ,----.   \\  \\   ---.   |  |    ;-.----.\n\
   /  ---./ ;-,----. / .---'|  , `.|  | / .---' ,--'  \\   \\  \\  |  |    | .---' \n\
  |  ( .---.|  ´` .'|  `--  |  |`  `  ||  `--  / .--.  )/ .-. \\ |  |    | `--,  \n\
   \\  --'  ||  |\\  \\ \\ `---.|  |  \\   | \\ `---.\\ `--' //  `-'  \\|  '---.| |`    \n\
    `----´ '`--' '--' `----'`--'   `--:  `----' `----' `------'.`------'| |     \n\
                                                                        | |     \n\
       ========================================================////-    ,-'     \n\
       " + grenedalf_version() + " (c) 2020-2024 by Lucas Czech\n";
}

inline std::string grenedalf_title()
{
    return "grenedalf: population genetic statistics for the next generation of pool sequencing";
}

#endif // include guard
