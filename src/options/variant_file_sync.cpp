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

#include "options/variant_file_sync.hpp"

#include "options/global.hpp"
#include "tools/misc.hpp"

#include "genesis/population/format/sync_input_stream.hpp"
#include "genesis/population/format/sync_reader.hpp"
#include "genesis/population/stream/variant_input_stream_sources.hpp"

#include <cassert>
#include <stdexcept>

// =================================================================================================
//      Setup Functions
// =================================================================================================

CLI::Option* VariantFileSyncOptions::add_file_input_opt_to_app_(
    CLI::App* sub,
    bool required,
    std::string const& group
) {
    // Correct setup check.
    internal_check(
        file_input_.option() == nullptr,
        "Cannot use the same VariantFileSyncOptions object multiple times."
    );

    // Add the option
    file_input_.add_multi_file_input_opt_to_app(
        sub, "sync", "sync (as specified by PoPoolation2)",
        "sync(\\.gz)?", "sync[.gz]",
        required, group
    );

    return file_input_.option();
}

// =================================================================================================
//      Run Functions
// =================================================================================================

VariantFileSyncOptions::VariantInputStream VariantFileSyncOptions::get_stream_(
    std::string const& filename
) const {
    using namespace genesis::population;

    // Make an stream.
    return make_variant_input_stream_from_sync_file( filename );
}
