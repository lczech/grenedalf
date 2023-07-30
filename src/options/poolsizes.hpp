#ifndef GRENEDALF_OPTIONS_POOLSIZES_H_
#define GRENEDALF_OPTIONS_POOLSIZES_H_

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

#include <functional>
#include <string>
#include <utility>
#include <vector>

// =================================================================================================
//      Poolsizes Options
// =================================================================================================

/**
 * @brief
 */
class PoolsizesOptions
{
public:

    // -------------------------------------------------------------------------
    //     Constructor and Rule of Five
    // -------------------------------------------------------------------------

    PoolsizesOptions()  = default;
    virtual ~PoolsizesOptions() = default;

    PoolsizesOptions( PoolsizesOptions const& other ) = default;
    PoolsizesOptions( PoolsizesOptions&& )            = default;

    PoolsizesOptions& operator= ( PoolsizesOptions const& other ) = default;
    PoolsizesOptions& operator= ( PoolsizesOptions&& )            = default;

    // -------------------------------------------------------------------------
    //     Setup Functions
    // -------------------------------------------------------------------------

    /**
     * @brief
     */
    void add_poolsizes_opt_to_app(
        CLI::App* sub,
        bool required = true,
        std::string const& group = "Settings"
    );

    // -------------------------------------------------------------------------
    //     Run Functions
    // -------------------------------------------------------------------------

    /**
     * @brief Get the pool sizes for the given samples, optionally filtered.
     *
     * If the @p sample_filter is provided, it needs to have the same length as the @p sample_names.
     * Then, only pool sizes for sample names for which the filter is true are required to be given.
     */
    std::vector<size_t> get_pool_sizes(
        std::vector<std::string> const& sample_names,
        std::vector<bool> const& sample_filter = std::vector<bool>{}
    ) const;

    // -------------------------------------------------------------------------
    //     Option Members
    // -------------------------------------------------------------------------

private:

    CliOption<std::string> poolsizes = "";

};

#endif // include guard
