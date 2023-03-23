#ifndef GRENEDALF_OPTIONS_WINDOW_H_
#define GRENEDALF_OPTIONS_WINDOW_H_

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

#include "CLI/CLI.hpp"

#include "tools/cli_option.hpp"
#include "options/variant_input.hpp"

#include "genesis/population/formats/variant_input_iterator.hpp"
#include "genesis/population/variant.hpp"
#include "genesis/population/window/base_window_iterator.hpp"
#include "genesis/population/window/base_window.hpp"
#include "genesis/population/window/chromosome_iterator.hpp"
#include "genesis/population/window/region_window_iterator.hpp"
#include "genesis/population/window/sliding_entries_window_iterator.hpp"
#include "genesis/population/window/sliding_interval_window_iterator.hpp"
#include "genesis/population/window/variant_window_iterator.hpp"
#include "genesis/population/window/window_view_iterator.hpp"
#include "genesis/population/window/window_view.hpp"
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
 * @brief Generic Window Iterator type that we are using throughout grenedalf.
 *
 * We use the base class for types of window iterators here,
 * meaning that this type has to be used through pointers.
 */
using VariantWindowIterator = genesis::population::VariantWindowIterator;

/**
 * @brief Generic Window View Iterator type that we are using throughout grenedalf.
 *
 * We use the base class for types of window iterators here,
 * meaning that this type has to be used through pointers.
 */
using VariantWindowViewIterator = genesis::population::VariantWindowViewIterator;

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
    //     Typedefs and Enums
    // -------------------------------------------------------------------------

    // Same as the above global typedefs, just to have them here as well.
    using VariantWindowIterator = genesis::population::VariantWindowIterator;
    using VariantWindowViewIterator = genesis::population::VariantWindowViewIterator;

    // Typedefs for the Window-based iterators.
    using VariantSlidingIntervalWindowIterator = genesis::population::SlidingIntervalWindowIterator<
        genesis::population::VariantInputIterator::Iterator
    >;
    using VariantSlidingEntriesWindowIterator = genesis::population::SlidingEntriesWindowIterator<
        genesis::population::VariantInputIterator::Iterator
    >;
    using VariantRegionWindowIterator = genesis::population::RegionWindowIterator<
        genesis::population::VariantInputIterator::Iterator
    >;

    // Typedefs for the WindowView-based iterators.
    using ChromosomeIterator = genesis::population::ChromosomeIterator<
        genesis::population::VariantInputIterator::Iterator
    >;

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

    /**
     * @brief Add options for selecting the windowing approach to a command.
     *
     * For efficiency reasons, genesis offers two types of iterators over windows: WindowIterator
     * and WindowViewIterator. The first one yields Windows that keep all their data in memory,
     * while the second one only points to existing data withing storing it, and is hence more
     * suitable for example when iterating a whole chromosome as a window.
     *
     * This means, we here also need to make this distinction. Some algorithms might need to iterate
     * a window multiple times in order to compute their thing, which means that all the data has
     * to be in memory (or we'd have to read the file multiple times... but that's not supported
     * at the moment). Hence, for these types of algorithms, we can only use the WindowIterator.
     *
     * For other algorithms that only need to stream through the data once, we can also offer the
     * WindowViewIterator types (at the time of writing: whole chromosomes, and the whole genome).
     *
     * Depending on the choice here, WindowViewIterator types are added and available, or not,
     * for the command where this class is used.
     */
    void add_window_opts_to_app(
        CLI::App* sub,
        bool include_window_view_types,
        std::string const& group = "Window"
    );

    // -------------------------------------------------------------------------
    //     Run Functions
    // -------------------------------------------------------------------------

public:

    /**
     * @brief Get a Window iterator over Variants, using the @p input to get its data from.
     *
     * This is meant to be called with the get_iterator() function of the VariantInputOptions class.
     *
     * It is only suppored if add_window_opts_to_app() above had been called with
     * `include_window_view_types = false`, as otherwise, Window View based iterators might be
     * requested by the user, which however cannot be packed into a VariantWindowIterator.
     * For this, use get_variant_window_view_iterator() instead, which wraps both types of
     * iterators into the same, namely into WindowViewIterator.
     */
    std::unique_ptr<VariantWindowIterator> get_variant_window_iterator(
        genesis::population::VariantInputIterator& input
    ) const;

    /**
     * @brief Get a Window View iterator over Variants, using the @p input to get its data from.
     *
     * This is meant to be called with the get_iterator() function of the VariantInputOptions class.
     */
    std::unique_ptr<VariantWindowViewIterator> get_variant_window_view_iterator(
        genesis::population::VariantInputIterator& input
    ) const;

    // -------------------------------------------------------------------------
    //     Internal Members
    // -------------------------------------------------------------------------

private:

    void check_options_() const;

    VariantSlidingIntervalWindowIterator get_variant_window_iterator_sliding_(
        genesis::population::VariantInputIterator& input
    ) const;

    VariantSlidingEntriesWindowIterator get_variant_window_iterator_queue_(
        genesis::population::VariantInputIterator& input
    ) const;

    VariantSlidingIntervalWindowIterator get_variant_window_iterator_single_(
        genesis::population::VariantInputIterator& input
    ) const;

    VariantRegionWindowIterator get_variant_window_iterator_regions_(
        genesis::population::VariantInputIterator& input
    ) const;

    ChromosomeIterator get_variant_window_view_iterator_chromosomes_(
        genesis::population::VariantInputIterator& input
    ) const;

    ChromosomeIterator get_variant_window_view_iterator_genome_(
        genesis::population::VariantInputIterator& input
    ) const;

    // -------------------------------------------------------------------------
    //     Option Members
    // -------------------------------------------------------------------------

private:

    // Window type selection
    bool include_window_view_types_ = false;
    CliOption<std::string> window_type_ = "sliding";

    // Sliding interval window settings
    CliOption<size_t> sliding_width_  = 0;
    CliOption<size_t> sliding_stride_ = 0;

    // Sliding entries window settings
    CliOption<size_t> queue_count_  = 0;
    CliOption<size_t> queue_stride_ = 0;

    // Region window settings
    CliOption<std::vector<std::string>> region_ ;
    CliOption<std::vector<std::string>> region_list_ ;
    CliOption<std::vector<std::string>> region_bed_ ;
    CliOption<std::vector<std::string>> region_gff_ ;
    CliOption<std::vector<std::string>> region_bim_ ;
    CliOption<std::vector<std::string>> region_vcf_ ;
    CliOption<bool> region_skip_empty_ = false;

};

#endif // include guard
