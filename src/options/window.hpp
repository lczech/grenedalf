#ifndef GRENEDALF_OPTIONS_WINDOW_H_
#define GRENEDALF_OPTIONS_WINDOW_H_

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

#include "CLI/CLI.hpp"

#include "tools/cli_option.hpp"
#include "options/variant_input.hpp"

#include "genesis/population/stream/variant_input_stream.hpp"
#include "genesis/population/variant.hpp"
#include "genesis/population/window/base_window_stream.hpp"
#include "genesis/population/window/base_window.hpp"
#include "genesis/population/window/chromosome_window_stream.hpp"
#include "genesis/population/window/genome_window_stream.hpp"
#include "genesis/population/window/interval_window_stream.hpp"
#include "genesis/population/window/position_window_stream.hpp"
#include "genesis/population/window/queue_window_stream.hpp"
#include "genesis/population/window/region_window_stream.hpp"
#include "genesis/population/window/variant_window_stream.hpp"
#include "genesis/population/window/window_view_stream.hpp"
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
 * @brief Window type that we use throughout, where each position in the Window is a Variant.
 */
using VariantWindow = ::genesis::population::Window<::genesis::population::Variant>;

/**
 * @brief Window type that we use throughout, where each position in the Window is a Variant.
 */
using VariantWindowView = ::genesis::population::WindowView<::genesis::population::Variant>;

/**
 * @brief Generic Window Stream type that we are using throughout grenedalf.
 *
 * We use the base class for types of window streams here,
 * meaning that this type has to be used through pointers.
 */
using VariantWindowStream = ::genesis::population::VariantWindowStream;

/**
 * @brief Generic Window View Stream type that we are using throughout grenedalf.
 *
 * We use the base class for types of window streams here,
 * meaning that this type has to be used through pointers.
 */
using VariantWindowViewStream = ::genesis::population::VariantWindowViewStream;

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

    enum class WindowType
    {
        kInterval,
        kQueue,
        kSingle,
        kRegions,
        kChromosomes,
        kGenome
    };

    // Same as the above global typedefs, just to have them here as well.
    using VariantWindowStream = genesis::population::VariantWindowStream;
    using VariantWindowViewStream = genesis::population::VariantWindowViewStream;

    // Typedefs for the Window-based streams.
    using VariantIntervalWindowStream = genesis::population::IntervalWindowStream<
        genesis::population::VariantInputStream::Iterator
    >;
    using VariantQueueWindowStream = genesis::population::QueueWindowStream<
        genesis::population::VariantInputStream::Iterator
    >;
    using VariantPositionWindowStream = genesis::population::PositionWindowStream<
        genesis::population::VariantInputStream::Iterator
    >;
    using VariantRegionWindowStream = genesis::population::RegionWindowStream<
        genesis::population::VariantInputStream::Iterator
    >;

    // Typedefs for the WindowView-based streams.
    using ChromosomeWindowStream = genesis::population::ChromosomeWindowStream<
        genesis::population::VariantInputStream::Iterator
    >;
    using GenomeWindowStream = genesis::population::GenomeWindowStream<
        genesis::population::VariantInputStream::Iterator
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
     * For efficiency reasons, genesis offers two types of streams over windows: WindowStream
     * and WindowViewStream. The first one yields Windows that keep all their data in memory,
     * while the second one only points to existing data withing storing it, and is hence more
     * suitable for example when iterating a whole chromosome as a window.
     *
     * This means, we here also need to make this distinction. Some algorithms might need to iterate
     * a window multiple times in order to compute their thing, which means that all the data has
     * to be in memory (or we'd have to read the file multiple times... but that's not supported
     * at the moment). Hence, for these types of algorithms, we can only use the WindowStream.
     *
     * For other algorithms that only need to stream through the data once, we can also offer the
     * WindowViewStream types (at the time of writing: whole chromosomes, and the whole genome).
     *
     * Depending on the choice here, WindowViewStream types are added and available, or not,
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

    /**
     * @brief Get the type of window that the user selected.
     */
    WindowType window_type() const;

    /**
     * @brief Get a Window stream over Variants, using the @p input to get its data from.
     *
     * This is meant to be called with the get_stream() function of the VariantInputOptions class.
     *
     * It is only suppored if add_window_opts_to_app() above had been called with
     * `include_window_view_types = false`, as otherwise, Window View based streams might be
     * requested by the user, which however cannot be packed into a VariantWindowStream.
     * For this, use get_variant_window_view_stream() instead, which wraps both types of
     * streams into the same, namely into WindowViewStream.
     */
    std::unique_ptr<VariantWindowStream> get_variant_window_stream(
        VariantInputOptions const& variant_input
    ) const;

    /**
     * @brief Get a Window View stream over Variants, using the @p input to get its data from.
     *
     * This is meant to be called with the get_stream() function of the VariantInputOptions class.
     */
    std::unique_ptr<VariantWindowViewStream> get_variant_window_view_stream(
        VariantInputOptions const& variant_input
    ) const;

    // -------------------------------------------------------------------------
    //     Reporting Functions
    // -------------------------------------------------------------------------

    /**
     * @brief Get the number of windows that have been processed in total.
     */
    size_t get_num_windows() const
    {
        return num_windows_;
    }

    /**
     * @brief Print a report of the number of processed windows.
     *
     * This includes the information offered by VariantInputOptions::print_report(),
     * so that if a command uses windows, only this one here needs to be called.
     */
    void print_report() const;

    // -------------------------------------------------------------------------
    //     Internal Members
    // -------------------------------------------------------------------------

private:

    void check_options_() const;

    VariantIntervalWindowStream get_variant_window_stream_interval_(
        genesis::population::VariantInputStream& input
    ) const;

    VariantQueueWindowStream get_variant_window_stream_queue_(
        genesis::population::VariantInputStream& input
    ) const;

    VariantPositionWindowStream get_variant_window_stream_single_(
        genesis::population::VariantInputStream& input,
        bool is_gapless_stream
    ) const;

    VariantRegionWindowStream get_variant_window_stream_regions_(
        genesis::population::VariantInputStream& input
    ) const;

    ChromosomeWindowStream get_variant_window_view_stream_chromosomes_(
        genesis::population::VariantInputStream& input
    ) const;

    GenomeWindowStream get_variant_window_view_stream_genome_(
        genesis::population::VariantInputStream& input
    ) const;

    // -------------------------------------------------------------------------
    //     Option Members
    // -------------------------------------------------------------------------

private:

    // -------------------------------------
    //     CLI Options
    // -------------------------------------

    // Window type selection
    bool include_window_view_types_ = false;
    CliOption<std::string> window_type_ = "interval";

    // Interval window settings
    CliOption<size_t> interval_width_  = 0;
    CliOption<size_t> interval_stride_ = 0;

    // Queue window settings
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

    // -------------------------------------
    //     Run Data
    // -------------------------------------

    mutable size_t num_windows_ = 0;
    mutable VariantInputOptions const* variant_input_ = nullptr;

};

#endif // include guard
