#ifndef GRENEDALF_OPTIONS_WINDOW_H_
#define GRENEDALF_OPTIONS_WINDOW_H_

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

#include "CLI/CLI.hpp"

#include "tools/cli_option.hpp"
#include "options/variant_input.hpp"

#include "genesis/population/formats/variant_input_iterator.hpp"
#include "genesis/population/variant.hpp"
#include "genesis/population/window/base_window_iterator.hpp"
#include "genesis/population/window/region_window_iterator.hpp"
#include "genesis/population/window/sliding_interval_window_iterator.hpp"
#include "genesis/population/window/window.hpp"

#include <functional>
#include <memory>
#include <string>
#include <utility>
#include <vector>

// =================================================================================================
//      Typedefs and Enums
// =================================================================================================

/**
 * @brief Generic window Iterator type that we are using throughout grenedalf.
 *
 * We use the base class for types of window iterators here,
 * meaning that this type has to be used through pointers.
 */
using VariantWindowIterator = genesis::population::BaseWindowIterator<
    genesis::population::VariantInputIterator::Iterator
>;

/**
 * @brief Typedef for sliding interval window iterator, for convenience.
 */
using VariantSlidingIntervalWindowIterator = genesis::population::SlidingIntervalWindowIterator<
    genesis::population::VariantInputIterator::Iterator
>;

/**
 * @brief Typedef for region window iterator, for convenience.
 */
using VariantRegionWindowIterator = genesis::population::RegionWindowIterator<
    genesis::population::VariantInputIterator::Iterator
>;

// =================================================================================================
//      Window Options
// =================================================================================================

/**
 * @brief
 */
class WindowOptions
{
public:

    // -------------------------------------------------------------------------
    //     Constructor and Rule of Five
    // -------------------------------------------------------------------------

    WindowOptions()  = default;
    ~WindowOptions() = default;

    WindowOptions( WindowOptions const& other ) = default;
    WindowOptions( WindowOptions&& )            = default;

    WindowOptions& operator= ( WindowOptions const& other ) = default;
    WindowOptions& operator= ( WindowOptions&& )            = default;

    // -------------------------------------------------------------------------
    //     Setup Functions
    // -------------------------------------------------------------------------

    void add_window_opts_to_app(
        CLI::App* sub,
        std::string const& group = "Window"
    );

    // -------------------------------------------------------------------------
    //     Run Functions
    // -------------------------------------------------------------------------
    //
    // /**
    //  * @brief Get the specified window width and stride.
    //  */
    // std::pair<size_t, size_t> get_window_width_and_stride() const;

    /**
     * @brief Get a window iterator that dereferences to Variant%s,
     * using the @p input to get its data from.
     *
     * This is meant to be called with the get_iterator() function of the grenedalf
     * VariantInputOptions class.
     */
    std::unique_ptr<VariantWindowIterator> get_variant_window_iterator(
        genesis::population::VariantInputIterator& input
    ) const;

    // -------------------------------------------------------------------------
    //     Internal Members
    // -------------------------------------------------------------------------

private:

    std::unique_ptr<VariantWindowIterator> get_variant_window_iterator_sliding_interval_(
        genesis::population::VariantInputIterator& input
    ) const;

    std::unique_ptr<VariantWindowIterator> get_variant_window_iterator_regions_(
        genesis::population::VariantInputIterator& input
    ) const;

    // -------------------------------------------------------------------------
    //     Option Members
    // -------------------------------------------------------------------------

private:

    // Window type selection
    CliOption<std::string> window_type_ = "sliding";

    // Sliding window settings
    CliOption<size_t> window_width_  = 0;
    CliOption<size_t> window_stride_ = 0;

    // Region window settings
    CliOption<std::vector<std::string>> region_ ;
    CliOption<std::vector<std::string>> region_list_ ;
    CliOption<std::vector<std::string>> region_bed_ ;
    CliOption<std::vector<std::string>> region_gff_ ;
    CliOption<std::vector<std::string>> region_bim_ ;
    CliOption<std::vector<std::string>> region_vcf_ ;
    CliOption<bool> skip_empty_regions_ = false;

};

#endif // include guard
