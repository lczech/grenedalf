#ifndef GRENEDALF_OPTIONS_VARIANT_FILTER_REGION_H_
#define GRENEDALF_OPTIONS_VARIANT_FILTER_REGION_H_

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

#include "genesis/population/stream/variant_input_stream.hpp"
#include "genesis/population/function/variant_input_stream.hpp"
#include "genesis/population/genome_locus_set.hpp"
#include "genesis/population/variant.hpp"

#include <functional>
#include <string>
#include <utility>
#include <vector>

// =================================================================================================
//      Variant Filter Region Options
// =================================================================================================

/**
 * @brief
 */
class VariantFilterRegionOptions
{
public:

    // -------------------------------------------------------------------------
    //     Typedefs and Enums
    // -------------------------------------------------------------------------

    using Variant = genesis::population::Variant;
    using GenomeLocusSet = genesis::population::GenomeLocusSet;
    using VariantInputStream = genesis::population::VariantInputStream;

    // -------------------------------------------------------------------------
    //     Constructor and Rule of Five
    // -------------------------------------------------------------------------

    VariantFilterRegionOptions()  = default;
    ~VariantFilterRegionOptions() = default;

    VariantFilterRegionOptions( VariantFilterRegionOptions const& other ) = default;
    VariantFilterRegionOptions( VariantFilterRegionOptions&& )            = default;

    VariantFilterRegionOptions& operator= ( VariantFilterRegionOptions const& other ) = default;
    VariantFilterRegionOptions& operator= ( VariantFilterRegionOptions&& )            = default;

    // -------------------------------------------------------------------------
    //     Setup Functions
    // -------------------------------------------------------------------------

    void add_region_filter_opts_to_app(
        CLI::App* sub,
        std::string const& group = "Region Filters"
    );

    // -------------------------------------------------------------------------
    //     Run Functions
    // -------------------------------------------------------------------------

    /**
     * @brief Parse the region filter files, e.g., BED or GFF, and make a filter from them.
     */
    void prepare_region_filters() const;

    std::shared_ptr<GenomeLocusSet> get_region_filter() const
    {
        prepare_region_filters();
        return region_filter_;
    }

    // -------------------------------------------------------------------------
    //     Option Members
    // -------------------------------------------------------------------------

private:

    // Filters for rows and columns
    CliOption<std::vector<std::string>> filter_region_;
    CliOption<std::vector<std::string>> filter_region_list_;
    CliOption<std::vector<std::string>> filter_region_bed_;
    CliOption<std::vector<std::string>> filter_region_gff_;
    CliOption<std::vector<std::string>> filter_region_bim_;
    CliOption<std::vector<std::string>> filter_region_vcf_;
    CliOption<std::string> filter_region_set_ = "union";

    // We keep the region filter here, so that it can be re-used for all inputs.
    // This filter is created by combining (union or intersection) all input filter files.
    mutable std::shared_ptr<GenomeLocusSet> region_filter_;

};

#endif // include guard
