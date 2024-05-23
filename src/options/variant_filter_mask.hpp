#ifndef GRENEDALF_OPTIONS_VARIANT_FILTER_MASK_H_
#define GRENEDALF_OPTIONS_VARIANT_FILTER_MASK_H_

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
#include "genesis/sequence/sequence_dict.hpp"

#include <functional>
#include <string>
#include <utility>
#include <vector>

// =================================================================================================
//      Variant Filter Mask Options
// =================================================================================================

/**
 * @brief
 */
class VariantFilterMaskOptions
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

    VariantFilterMaskOptions()  = default;
    ~VariantFilterMaskOptions() = default;

    VariantFilterMaskOptions( VariantFilterMaskOptions const& other ) = default;
    VariantFilterMaskOptions( VariantFilterMaskOptions&& )            = default;

    VariantFilterMaskOptions& operator= ( VariantFilterMaskOptions const& other ) = default;
    VariantFilterMaskOptions& operator= ( VariantFilterMaskOptions&& )            = default;

    // -------------------------------------------------------------------------
    //     Setup Functions
    // -------------------------------------------------------------------------

    void add_mask_filter_opts_to_app(
        CLI::App* sub,
        std::string const& group = "Masking Filters"
    );

    // -------------------------------------------------------------------------
    //     Run Functions
    // -------------------------------------------------------------------------

    /**
     * @brief Parse the region filter files, e.g., BED or GFF, and make a filter from them.
     */
    void prepare_mask() const;

    /**
     * @brief Get the mask, where set bits indicate positions that are masked out.
     */
    std::shared_ptr<GenomeLocusSet> get_mask() const
    {
        prepare_mask();
        return mask_;
    }

    /**
     * @brief Check that the mask fits with a given reference genome, if given.
     */
    void check_mask_against_reference(
        std::shared_ptr<genesis::sequence::SequenceDict> ref_dict
    ) const;

    /**
     * @brief Create the tranform function to be applied to the stream for masking.
     */
    std::function<void( genesis::population::Variant& )> make_mask_transform() const;

    // -------------------------------------------------------------------------
    //     Option Members
    // -------------------------------------------------------------------------

private:

    // Filters for rows and columns
    CliOption<std::string> filter_mask_bed_;
    CliOption<std::string> filter_mask_fasta_;
    CliOption<size_t>      filter_mask_fasta_min_ = 0;
    CliOption<bool>        filter_mask_fasta_inv_ = false;

    // We keep the mask here. Bits that are set here are masked!
    mutable std::shared_ptr<GenomeLocusSet> mask_;

};

#endif // include guard
