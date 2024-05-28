#ifndef GRENEDALF_OPTIONS_VARIANT_TRANSFORM_SUBSAMPLE_H_
#define GRENEDALF_OPTIONS_VARIANT_TRANSFORM_SUBSAMPLE_H_

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

#include <string>
#include <vector>

// =================================================================================================
//      Variant Transform Subsample Options
// =================================================================================================

/**
 * @brief
 */
class VariantTransformSubsampleOptions
{
public:

    // -------------------------------------------------------------------------
    //     Constructor and Rule of Five
    // -------------------------------------------------------------------------

    VariantTransformSubsampleOptions()  = default;
    ~VariantTransformSubsampleOptions() = default;

    VariantTransformSubsampleOptions( VariantTransformSubsampleOptions const& ) = default;
    VariantTransformSubsampleOptions( VariantTransformSubsampleOptions&& )      = default;

    VariantTransformSubsampleOptions& operator= ( VariantTransformSubsampleOptions const& ) = default;
    VariantTransformSubsampleOptions& operator= ( VariantTransformSubsampleOptions&& )      = default;

    // -------------------------------------------------------------------------
    //     Setup Functions
    // -------------------------------------------------------------------------

    void add_subsample_opts_to_app(
        CLI::App* sub,
        std::string const& group = "Sample Subsampling"
    );

    // -------------------------------------------------------------------------
    //     Run Functions
    // -------------------------------------------------------------------------

    void add_subsample_transformation( VariantInputOptions const& variant_input ) const;

    // TODO add rescaling as well?

    // -------------------------------------------------------------------------
    //     Option Members
    // -------------------------------------------------------------------------

private:

    // Sample subsample transform
    CliOption<size_t> max_read_depth_ = 0;
    CliOption<std::string> method_ = "subscale";

};

#endif // include guard
