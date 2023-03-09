#ifndef GRENEDALF_TOOLS_MISC_H_
#define GRENEDALF_TOOLS_MISC_H_

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

#include "tools/cli_option.hpp"

#include "CLI/CLI.hpp"

#include "genesis/utils/text/string.hpp"

#include <iosfwd>
#include <string>
#include <stdexcept>
#include <vector>

// =================================================================================================
//      Legacy Commands
// =================================================================================================

class RenamedCommandError : public std::runtime_error {

public:

    RenamedCommandError( std::string message )
        : std::runtime_error( message )
    {}
};

void add_legacy_command(
    CLI::App& app,
    std::string const& old_name,
    std::string const& new_path
);

// =================================================================================================
//      Formatting
// =================================================================================================

std::string format_columns(
    std::string const& left,
    std::string const& right,
    size_t left_width
);

void write_columns(
    std::ostream& out,
    std::string const& left,
    std::string const& right,
    size_t left_width,
    size_t right_width
);

// =================================================================================================
//      Misc
// =================================================================================================

/**
 * @brief Alternative for normal `assert()` that allows to specify an error message,
 * throws an exception instead of terminating, and is always used, also in release mode.
 */
inline void internal_check(
    bool condition,
    std::string const& error_message
) {
    if( ! condition ) {
        throw std::domain_error(
            "Internal error: " + error_message
        );
    }
}

/**
 * @brief Helper function to get the keys of a map (we use a vector of pairs to keep the order).
 * This allows to easilty set up CLI options with enums and translate to the underlying value
 * from a string, while also giving nicer user output and help than the CLI11 default for enums.
 */
template<class T>
std::vector<std::string> enum_map_keys( std::vector<std::pair<std::string, T>> const& map )
{
    std::vector<std::string> result;
    result.reserve( map.size() );
    for( auto const& kv : map ) {
        result.emplace_back( kv.first );
    }
    return result;
}

/**
 * @brief Helper function to translate the keys of a map (we use a vector of pairs to keep the order).
 * This allows to easilty set up CLI options with enums and translate to the underlying value
 * from a string, while also giving nicer user output and help than the CLI11 default for enums.
 */
template<class T>
T get_enum_map_value( std::vector<std::pair<std::string, T>> const& map, std::string const& key )
{
    // Case insensitive comparison.
    auto const key_lower = genesis::utils::to_lower( key );
    for( auto const& kv : map ) {
        if( genesis::utils::to_lower( kv.first ) == key_lower ) {
            return kv.second;
        }
    }

    // This case should not happen if we used enum_map_keys before to restrict the CLI11 options,
    // as CLI11 will to the check already for us.
    throw std::domain_error(
        "Internal error: Key \"" + key + "\" not found in list of possible values."
    );
}

/**
 * @brief Helper function to get the char representation for table separator chars,
 * given its textual description from the option.
 */
char translate_separator_char( CliOption<std::string> const& option );

#endif // include guard
