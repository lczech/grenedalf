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

#include "options/window_average.hpp"

#include "options/global.hpp"
#include "tools/misc.hpp"

#include "genesis/utils/text/string.hpp"

#include <cassert>
#include <stdexcept>

// =================================================================================================
//      Enum Mapping
// =================================================================================================

using namespace genesis::population;

std::vector<std::pair<std::string, WindowAveragePolicy>> const window_average_policy_map = {
    { "window-length",  WindowAveragePolicy::kWindowLength },
    { "available-loci", WindowAveragePolicy::kAvailableLoci },
    { "valid-loci",     WindowAveragePolicy::kValidLoci },
    { "valid-snps",     WindowAveragePolicy::kValidSnps },
    { "sum",            WindowAveragePolicy::kAbsoluteSum }
};

// =================================================================================================
//      Setup Functions
// =================================================================================================

CLI::Option* WindowAverageOptions::add_window_average_opt_to_app(
    CLI::App* sub,
    std::string const& group
) {
    // Correct setup check.
    internal_check(
        window_average_policy_.option == nullptr,
        "Cannot use the same WindowAverageOptions object multiple times."
    );

    // Add the option. We need quite the bit of documentation here...
    window_average_policy_.option = sub->add_option(
        "--window-average-policy",
        window_average_policy_.value,
        "Denominator to use when computing the average of a metric in a window: "
        "\n(1) `window-length`: Simply use the window length, which likely underestimates the metric, "
        "in particular in regions with low coverage and high missing data."
        "\n(2) `available-loci`: Use the number of positions for which there was data at all, "
        "independent of all filter settings."
        "\n(3) `valid-loci`: Use the number of positions that passed all quality and numerical "
        "filters (that is, excluding the SNP-related filters). This uses all positions of high "
        "quality, and is the recommended policy when the input contains data for all positions."
        "\n(4) `valid-snps`: Use the number of SNPs only. This might overestimate the metric, "
        "but can be useful when the data only consists of SNPs."
        "\n(5) `sum`: Simply report the sum of the per-site values, with no averaging "
        "applied to it. This can be used to apply custom averaging later."
    );
    window_average_policy_.option->transform(
        CLI::IsMember( enum_map_keys( window_average_policy_map ), CLI::ignore_case )
    );
    window_average_policy_.option->group( group );
    window_average_policy_.option->required();
    return window_average_policy_.option;
}

// =================================================================================================
//      Run Functions
// =================================================================================================

genesis::population::WindowAveragePolicy WindowAverageOptions::get_window_average_policy() const
{
    return get_enum_map_value(
        window_average_policy_map,
        window_average_policy_.value
    );
}
