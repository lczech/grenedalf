#ifndef GRENEDALF_OPTIONS_VARIANT_SAMPLE_NAMES_H_
#define GRENEDALF_OPTIONS_VARIANT_SAMPLE_NAMES_H_

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

#include "genesis/population/streams/variant_input_stream.hpp"
#include "genesis/population/variant.hpp"

#include <string>
#include <vector>

// =================================================================================================
//      Variant Sample Names Options
// =================================================================================================

/**
 * @brief
 */
class VariantSampleNamesOptions
{
public:

    // -------------------------------------------------------------------------
    //     Constructor and Rule of Five
    // -------------------------------------------------------------------------

    VariantSampleNamesOptions()  = default;
    ~VariantSampleNamesOptions() = default;

    VariantSampleNamesOptions( VariantSampleNamesOptions const& other ) = default;
    VariantSampleNamesOptions( VariantSampleNamesOptions&& )            = default;

    VariantSampleNamesOptions& operator= ( VariantSampleNamesOptions const& other ) = default;
    VariantSampleNamesOptions& operator= ( VariantSampleNamesOptions&& )            = default;

    // -------------------------------------------------------------------------
    //     Setup Functions
    // -------------------------------------------------------------------------

    void add_sample_name_opts_to_app(
        CLI::App* sub,
        std::string const& group = "Sample Names and Filters"
    );

    // -------------------------------------------------------------------------
    //     Run Functions
    // -------------------------------------------------------------------------

    void rename_samples( std::vector<std::string>& sample_names ) const;

    void add_sample_name_filter( genesis::population::VariantInputStream& iterator ) const;

private:

    std::vector<std::string> process_sample_name_list_option_( std::string const& list ) const;

    // -------------------------------------------------------------------------
    //     Option Members
    // -------------------------------------------------------------------------

private:

    // Sample names and filters
    CliOption<std::string> rename_samples_ = "";
    CliOption<std::string> filter_samples_include_ = "";
    CliOption<std::string> filter_samples_exclude_ = "";

};

#endif // include guard
