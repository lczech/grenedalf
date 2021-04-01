#ifndef GRENEDALF_TOOLS_CLI_OPTION_H_
#define GRENEDALF_TOOLS_CLI_OPTION_H_

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
//      CLI11 Option Helper
// =================================================================================================

/**
 * @brief Helper that encapsulates an option for the command line interface,
 * storing its value and the CLI11 object used in the interface to change that value.
 */
template<typename T>
struct CliOption
{
    CliOption() = default;

    CliOption( T const& val )
        : value( val )
    {}

    CliOption& operator =( CLI::Option* opt )
    {
        option = opt;
        return *this;
    }

    T            value  = {};
    CLI::Option* option = nullptr;
};

/**
 * @brief Specialized version that allows construction from char arrays,
 * so that we can easility initialize the class in standard use cases.
 */
template<>
struct CliOption<std::string>
{
    CliOption() = default;

    CliOption( std::string const& val )
        : value( val )
    {}

    CliOption( char const* val )
        : value( val )
    {}

    CliOption& operator =( CLI::Option* opt )
    {
        option = opt;
        return *this;
    }

    // CliOption& operator =( std::string const& val )
    // {
    //     value = val;
    //     return *this;
    // }
    //
    // CliOption& operator =( char const* val )
    // {
    //     value = val;
    //     return *this;
    // }

    std::string  value  = {};
    CLI::Option* option = nullptr;
};

#endif // include guard
