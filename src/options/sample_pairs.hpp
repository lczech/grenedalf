#ifndef GRENEDALF_OPTIONS_SAMPLE_PAIRS_H_
#define GRENEDALF_OPTIONS_SAMPLE_PAIRS_H_

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

#include <string>
#include <utility>
#include <tuple>
#include <vector>

// =================================================================================================
//      Sample Pairs Options
// =================================================================================================

/**
 * @brief Select pairs of samples.
 *
 * This set of options allows to select a set of sample pairs, for instance to compute pairwise
 * measures such as FST. Using these options, the subset of samples can be selected by the user.
 */
class SamplePairsOptions
{
public:

    // -------------------------------------------------------------------------
    //     Constructor and Rule of Five
    // -------------------------------------------------------------------------

    SamplePairsOptions()  = default;
    virtual ~SamplePairsOptions() = default;

    SamplePairsOptions( SamplePairsOptions const& other ) = default;
    SamplePairsOptions( SamplePairsOptions&& )            = default;

    SamplePairsOptions& operator= ( SamplePairsOptions const& other ) = default;
    SamplePairsOptions& operator= ( SamplePairsOptions&& )            = default;

    // -------------------------------------------------------------------------
    //     Setup Functions
    // -------------------------------------------------------------------------

    void add_sample_pairs_opts_to_app(
        CLI::App* sub,
        std::string const& group = "Settings"
    );

    // -------------------------------------------------------------------------
    //     Run Functions
    // -------------------------------------------------------------------------

    bool is_all_to_all() const;

    std::vector<std::pair<size_t, size_t>> get_sample_pairs(
        std::vector<std::string> const& sample_names
    ) const;

    // -------------------------------------------------------------------------
    //     Option Members
    // -------------------------------------------------------------------------

private:

    CliOption<std::string> comparand = "";
    CliOption<std::string> second_comparand = "";
    CliOption<std::string> comparand_list = "";

};

#endif // include guard
