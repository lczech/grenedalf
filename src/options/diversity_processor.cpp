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

#include "options/diversity_processor.hpp"

#include "options/global.hpp"
#include "options/variant_filter_numerical.hpp"
#include "tools/cli_setup.hpp"
#include "tools/misc.hpp"

#include "genesis/population/filter/sample_counts_filter_numerical.hpp"
#include "genesis/population/filter/sample_counts_filter.hpp"
#include "genesis/population/filter/variant_filter_numerical.hpp"
#include "genesis/population/filter/variant_filter.hpp"
#include "genesis/population/function/diversity_pool_calculator.hpp"
#include "genesis/population/function/diversity_pool_functions.hpp"
#include "genesis/population/function/functions.hpp"
#include "genesis/population/window/functions.hpp"
#include "genesis/utils/core/algorithm.hpp"
#include "genesis/utils/core/fs.hpp"
#include "genesis/utils/core/logging.hpp"
#include "genesis/utils/core/std.hpp"
#include "genesis/utils/text/convert.hpp"
#include "genesis/utils/text/string.hpp"

#include <algorithm>
#include <cassert>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <tuple>
#include <vector>

// =================================================================================================
//      Enum Mapping
// =================================================================================================

using namespace genesis::population;

std::vector<std::pair<std::string, TajimaDenominatorPolicy>> const
tajima_d_denominator_policy_map = {
    { "empirical-min-read-depth", TajimaDenominatorPolicy::kEmpiricalMinReadDepth },
    { "provided-min-read-depth",  TajimaDenominatorPolicy::kProvidedMinReadDepth },
    { "popoolation-bugs",         TajimaDenominatorPolicy::kWithPopoolationBugs },
    { "pool-size",                TajimaDenominatorPolicy::kPoolsize },
    { "uncorrected",              TajimaDenominatorPolicy::kUncorrected }
};

// =================================================================================================
//      Setup Functions
// =================================================================================================

void DiversityProcessorOptions::add_diversity_processor_opts_to_app(
    CLI::App* sub,
    std::string const& group
) {

    // -------------------------------------------------------------------------
    //     Basic Diversity Settings
    // -------------------------------------------------------------------------

    // Pool sizes and window averaging
    poolsizes.add_poolsizes_opt_to_app( sub, true, group );
    window_average_policy.add_window_average_opt_to_app( sub, group );

    // Settings: Tajima Denominator Policy
    tajima_d_denominator_policy.option = sub->add_option(
        "--tajima-d-denominator-policy",
        tajima_d_denominator_policy.value,
        "Estimator for the expected number of distinct individuals sequenced in the denominator of "
        "the Achaz (2008) correction of Tajima's D, following its adaptation by PoPoolation. "
        "With pool seq data, there is no simple way to obtain a statistic that is numerically "
        "comparable to the classic Tajima's D with individual data. Hence, all of the below are "
        "simplicications that introduce some bias."
        "\n(1) `empirical-min-read-depth`: Use the lowest empirical read depth found in each window, "
        "and the pool size, to compute the expected number of individuals sequenced. "
        "This is a conservative estimator that we recommend by default."
        "\n(2) `provided-min-read-depth`: Same as (1), but use the user-provided "
        "`--filter-sample-min-read-depth` instead of the empirical minum read depth. "
        "This is what PoPoolation uses."
        "\n(3) `popoolation-bugs`: Same as (2), but additionally re-introduce their bugs. "
        "We offer this for comparability with PoPoolation."
        "\n(4) `pool-size`: Directly use the pool size as an estimate of the number of individuals, "
        "instead of computing the expected value. This assumes the number of individuals sequenced "
        "to be equal to the pool size, and is good under high read depths."
        "\n(5) `uncorrected`: The Achaz correction is not applied, so that the result is simply "
        "Theta Pi minus Theta Watterson. Hence, magnitudes of values are not comparable to classic "
        "Tajima's D. Still, using their sign, and comparing them across windows can be useful."
    );
    tajima_d_denominator_policy.option->transform(
        CLI::IsMember( enum_map_keys( tajima_d_denominator_policy_map ), CLI::ignore_case )
    );
    tajima_d_denominator_policy.option->group( group );
    // tajima_d_denominator_policy.option->required();

    // Also add some switches to turn off specific parts of the computation.

    // No Theta Pi
    no_theta_pi.option = sub->add_flag(
        "--no-theta-pi",
        no_theta_pi.value,
        "Do not compute or output Theta Pi. Note that if Tajima's D is computed, "
        "we also need to compute the two thetas, but will then not print them in the output."
    );
    no_theta_pi.option->group( group );

    // No Theta Watterson
    no_theta_watterson.option = sub->add_flag(
        "--no-theta-watterson",
        no_theta_watterson.value,
        "Do not compute or output Theta Watterson. Note that if Tajima's D is computed, "
        "we also need to compute the two thetas, but will then not print them in the output."
    );
    no_theta_watterson.option->group( group );

    // No Tajima's D
    no_tajima_d.option = sub->add_flag(
        "--no-tajima-d",
        no_tajima_d.value,
        "Do not compute Tajmias' D."
    );
    no_tajima_d.option->group( group );

    // -------------------------------------------------------------------------
    //     Misc
    // -------------------------------------------------------------------------

    // Setting: threading_threshold
    threading_threshold.option = sub->add_option(
        "--threading-threshold",
        threading_threshold.value,
        "When computing the diversity of a few samples, using individual threads for each "
        "usually incurs a substantial overhead due to thread synchronozation. Hence, we only want "
        "to use threads for the computation if many samples are being processed. "
        "This setting determiens the number of samples at which threads are used. "
        "(Note that we still always use threads for input file parsing.)"
    );
    threading_threshold.option->group( "" );
}

// =================================================================================================
//      Run Functions
// =================================================================================================

// -------------------------------------------------------------------------
//     get_diversity_pool_processor
// -------------------------------------------------------------------------

genesis::population::DiversityPoolProcessor DiversityProcessorOptions::get_diversity_pool_processor(
    std::vector<std::string> const& sample_names,
    genesis::population::SampleCountsFilterNumericalParams filter_params
) const {
    using namespace genesis::population;

    // User input check
    if( no_theta_pi.value && no_theta_watterson.value && no_tajima_d.value ) {
        throw CLI::ValidationError(
            no_theta_pi.option->get_name() + ", " +
            no_theta_watterson.option->get_name() + ", " +
            no_tajima_d.option->get_name(),
            "All diversity metrics have been deactivate. Nothing is being computed."
        );
    }

    // Get and check the pool sizes.
    auto const pool_sizes = poolsizes.get_pool_sizes( sample_names );
    internal_check(
        pool_sizes.size() == sample_names.size(),
        "Inconsistent number of samples and number of pool sizes."
    );

    // Get the window average policy to use for all processors.
    auto const win_avg_policy = window_average_policy.get_window_average_policy();

    // Prepare the diversity settings we want to use.
    DiversityPoolSettings settings;
    settings.min_count      = filter_params.min_count;
    settings.min_read_depth = filter_params.min_read_depth;
    settings.max_read_depth = filter_params.max_read_depth;
    settings.tajima_denominator_policy = get_enum_map_value(
        tajima_d_denominator_policy_map, tajima_d_denominator_policy.value
    );

    // Make the type of processor that we need for the provided method.
    DiversityPoolProcessor processor = make_diversity_pool_processor(
        win_avg_policy, settings, pool_sizes
    );

    // Set the per-calculator settings
    for( auto& calculator : processor ) {
        calculator.enable_theta_pi(        !no_theta_pi.value );
        calculator.enable_theta_watterson( !no_theta_watterson.value );
        calculator.enable_tajima_d(        !no_tajima_d.value );
    }

    // Set the threading options as provided by the (currently hidden) setting.
    processor.thread_pool( global_options.thread_pool() );
    processor.threading_threshold( threading_threshold.value );
    return processor;
}
