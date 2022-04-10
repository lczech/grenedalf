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

#include "options/window.hpp"

#include "options/global.hpp"
#include "tools/misc.hpp"

#include <algorithm>
#include <cassert>
#include <memory>
#include <stdexcept>

// =================================================================================================
//      Setup Functions
// =================================================================================================

void WindowOptions::add_window_opts_to_app(
    CLI::App* sub,
    std::string const& group
) {
    // Width
    window_width_.option = sub->add_option(
        "--window-width",
        window_width_.value,
        "Width of each window along the chromosome."
    );
    window_width_.option->group( group );
    window_width_.option->required();

    // Stride
    window_stride_.option = sub->add_option(
        "--window-stride",
        window_stride_.value,
        "Stride between windows along the chromosome, that is how far to move to get to the next "
        "window. If set to 0 (default), this is set to the same value as the `--window-width`."
    );
    window_stride_.option->group( group );
}

// =================================================================================================
//      Run Functions
// =================================================================================================

std::pair<size_t, size_t> WindowOptions::get_window_width_and_stride() const
{
    // We need to check the stride here ourselves. For the actual window iterator, this is done
    // in the iterator constructor, but when this function here is called, that might not have
    // happened yet, and so we need to do the check ourselves.
    if( window_stride_.value == 0 ) {
        // stride == 0 --> use stride == width
        return { window_width_.value, window_width_.value };
    } else {
        // stride != 0 --> just use it as is
        return { window_width_.value, window_stride_.value };
    }
}

// -------------------------------------------------------------------------
//     get_variant_window_iterator
// -------------------------------------------------------------------------

VariantWindowIterator WindowOptions::get_variant_window_iterator(
    genesis::population::VariantInputIterator& input
) const {
    using namespace genesis;
    using namespace genesis::population;
    return make_default_sliding_interval_window_iterator(
        input.begin(), input.end(),  window_width_.value, window_stride_.value
    );
}
