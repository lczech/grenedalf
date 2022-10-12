/*
    grenedalf - Genome Analyses of Differential Allele Frequencies
    Copyright (C) 2020-2022 Lucas Czech

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

#include "options/variant_input_sync.hpp"

#include "options/global.hpp"
#include "options/variant_input_sample_names.hpp"
#include "tools/misc.hpp"

#include "genesis/population/formats/sync_input_iterator.hpp"
#include "genesis/population/formats/sync_reader.hpp"

#include <cassert>
#include <stdexcept>

// =================================================================================================
//      Setup Functions
// =================================================================================================

CLI::Option* VariantInputSyncOptions::add_file_input_opt_to_app_(
    CLI::App* sub,
    bool required,
    std::string const& group
) {
    // Correct setup check.
    internal_check(
        file_input_.option() == nullptr,
        "Cannot use the same VariantInputSyncOptions object multiple times."
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

VariantInputSyncOptions::VariantInputIterator VariantInputSyncOptions::get_iterator_(
    std::string const& filename,
    VariantInputSampleNamesOptions const& sample_names_options
) const {
    using namespace genesis::population;

    // We can use the sample filter settings, see prepare_pileup_iterator_() for details.
    auto const sample_filter = sample_names_options.find_sample_indices_from_sample_filters();

    // Make an iterator.
    auto iterator = make_variant_input_iterator_from_sync_file(
        filename, sample_filter.first, sample_filter.second
    );

    // Sync does not have sample names, so set them based on user input or simple enumeration.
    // make_variant_input_iterator_from_sync_file() returns a list with as many
    // empty strings as the file has samples (exactly one in that case).
    iterator.data().sample_names = sample_names_options.make_anonymous_sample_names(
        iterator.data().sample_names.size()
    );

    return iterator;
}
