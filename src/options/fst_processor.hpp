#ifndef GRENEDALF_OPTIONS_FST_PROCESSOR_H_
#define GRENEDALF_OPTIONS_FST_PROCESSOR_H_

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

#include "options/poolsizes.hpp"
#include "options/window_average.hpp"
#include "tools/cli_option.hpp"

#include "genesis/population/function/fst_pool_processor.hpp"
#include "genesis/population/function/window_average.hpp"

#include <string>
#include <utility>
#include <tuple>
#include <vector>

// =================================================================================================
//      Fst Processor Options
// =================================================================================================

/**
 * @brief Create an Fst Processor, including the selection of pairs of samples.
 *
 * This set of options allows to select a set of sample pairs, for instance to compute pairwise
 * measures such as FST. Using these options, the subset of samples can be selected by the user.
 */
class FstProcessorOptions
{
public:

    // -------------------------------------------------------------------------
    //     Typedefs and Enums
    // -------------------------------------------------------------------------

    enum class FstMethod
    {
        kUnbiasedNei,
        kUnbiasedHudson,
        kKofler,
        kKarlsson
    };

    struct Params
    {
        // Offer all FST methods (default), or just the
        // unbiased Nei/Hudson variants (when set to false).
        bool with_all_methods = true;

        // Offer to chose the window average policy (default),
        // or use the fixed one instead.
        bool with_window_average_policy = true;
        genesis::population::WindowAveragePolicy fix_window_average_policy;
    };

    // -------------------------------------------------------------------------
    //     Constructor and Rule of Five
    // -------------------------------------------------------------------------

    FstProcessorOptions()  = default;
    virtual ~FstProcessorOptions() = default;

    FstProcessorOptions( FstProcessorOptions const& other ) = default;
    FstProcessorOptions( FstProcessorOptions&& )            = default;

    FstProcessorOptions& operator= ( FstProcessorOptions const& other ) = default;
    FstProcessorOptions& operator= ( FstProcessorOptions&& )            = default;

    // -------------------------------------------------------------------------
    //     Setup Functions
    // -------------------------------------------------------------------------

    void add_fst_processor_opts_to_app(
        CLI::App* sub,
        std::string const& group = "Settings"
    ) {
        // Yet again we need an overload due to the stupid compiler bug about
        // default constructed arguments, see https://stackoverflow.com/q/43819314/4184258
        Params params;
        return add_fst_processor_opts_to_app( sub, params, group );
    }

    void add_fst_processor_opts_to_app(
        CLI::App* sub,
        Params const& fst_processor_params,
        std::string const& group = "Settings"
    );

    // -------------------------------------------------------------------------
    //     Run Functions
    // -------------------------------------------------------------------------

    bool is_all_to_all() const;

    FstMethod get_fst_method() const;

    std::vector<std::pair<size_t, size_t>> get_sample_pairs(
        std::vector<std::string> const& sample_names
    ) const;

    genesis::population::FstPoolProcessor get_fst_pool_processor(
        std::vector<std::string> const& sample_names,
        std::vector<std::pair<size_t, size_t>> const& sample_pairs
    ) const;

    // -------------------------------------------------------------------------
    //     Option Members
    // -------------------------------------------------------------------------

private:

    CliOption<std::string> method = "unbiased-nei";
    PoolsizesOptions       poolsizes;
    WindowAverageOptions   window_average_policy;
    Params                 params;

    CliOption<std::string> comparand = "";
    CliOption<std::string> second_comparand = "";
    CliOption<std::string> comparand_list = "";

    CliOption<size_t> threading_threshold = 4096;

};

#endif // include guard
