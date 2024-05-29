#ifndef GRENEDALF_OPTIONS_DIVERSITY_PROCESSOR_H_
#define GRENEDALF_OPTIONS_DIVERSITY_PROCESSOR_H_

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

#include "genesis/population/filter/sample_counts_filter_numerical.hpp"
#include "genesis/population/filter/sample_counts_filter.hpp"
#include "genesis/population/function/diversity_pool_processor.hpp"

#include <string>
#include <utility>
#include <tuple>
#include <vector>

// =================================================================================================
//      Diversity Processor Options
// =================================================================================================

/**
 * @brief Create a Diversity Processor, including all policy settings.
 */
class DiversityProcessorOptions
{
public:

    // -------------------------------------------------------------------------
    //     Typedefs and Enums
    // -------------------------------------------------------------------------

    // -------------------------------------------------------------------------
    //     Constructor and Rule of Five
    // -------------------------------------------------------------------------

    DiversityProcessorOptions()  = default;
    virtual ~DiversityProcessorOptions() = default;

    DiversityProcessorOptions( DiversityProcessorOptions const& other ) = default;
    DiversityProcessorOptions( DiversityProcessorOptions&& )            = default;

    DiversityProcessorOptions& operator= ( DiversityProcessorOptions const& other ) = default;
    DiversityProcessorOptions& operator= ( DiversityProcessorOptions&& )            = default;

    // -------------------------------------------------------------------------
    //     Setup Functions
    // -------------------------------------------------------------------------

    void add_diversity_processor_opts_to_app(
        CLI::App* sub,
        std::string const& group = "Settings"
    );

    // -------------------------------------------------------------------------
    //     Run Functions
    // -------------------------------------------------------------------------

    genesis::population::DiversityPoolProcessor get_diversity_pool_processor(
        std::vector<std::string> const& sample_names,
        genesis::population::SampleCountsFilterNumericalParams filter_params
    ) const;

    // -------------------------------------------------------------------------
    //     Option Members
    // -------------------------------------------------------------------------

private:

    // Diversity specific settings
    CliOption<std::string> tajima_d_denominator_policy = "empirical-min-read-depth";
    CliOption<bool> no_theta_pi        = false;
    CliOption<bool> no_theta_watterson = false;
    CliOption<bool> no_tajima_d        = false;

    // General settings
    PoolsizesOptions       poolsizes;
    WindowAverageOptions   window_average_policy;

    // Remnant settings. These are simply re-used from the numerical filters now.
    // Kept here for reference, just in case they are needed for later.
    // CliOption<size_t> min_count = 2;
    // CliOption<size_t> min_read_depth = 4;
    // CliOption<size_t> max_read_depth = 1000000;
    // CliOption<double> min_read_depth_fraction = 0.6;

    CliOption<size_t> threading_threshold = 4096;

};

#endif // include guard
