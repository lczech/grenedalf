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

#include "options/variant_reference_genome.hpp"
#include "tools/cli_option.hpp"

#include "genesis/population/stream/variant_input_stream.hpp"
#include "genesis/population/function/variant_input_stream.hpp"
#include "genesis/population/genome_locus_set.hpp"
#include "genesis/population/variant.hpp"
#include "genesis/sequence/sequence_dict.hpp"

#include <functional>
#include <string>
#include <unordered_map>
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

    using Variant            = genesis::population::Variant;
    using GenomeLocusSet     = genesis::population::GenomeLocusSet;
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
        VariantReferenceGenomeOptions const& ref_genome_opts,
        std::string const& group = "Masking Filters"
    );

    void add_mask_filter_sample_opts_to_app(
        CLI::App* sub,
        VariantReferenceGenomeOptions const& ref_genome_opts,
        std::string const& group = "Masking Filters"
    );

    void add_mask_filter_total_opts_to_app(
        CLI::App* sub,
        VariantReferenceGenomeOptions const& ref_genome_opts,
        std::string const& group = "Masking Filters"
    );

    // -------------------------------------------------------------------------
    //     Run Functions
    // -------------------------------------------------------------------------

    /**
     * @brief Parse the region filter files, e.g., BED or GFF, and make a filter from them.
     */
    void prepare_masks() const;

    /**
     * @brief Get the sample masks, where set bits indicated positions that are masked out,
     * mapped from sample names to their masks.
     */
    std::unordered_map<std::string, std::shared_ptr<GenomeLocusSet>> const& get_sample_masks() const
    {
        prepare_sample_masks_();
        return sample_masks_;
    }

    /**
     * @brief Get the mask, where set bits indicate positions that are masked out.
     */
    std::shared_ptr<GenomeLocusSet> get_total_mask() const
    {
        prepare_total_mask_();
        return total_mask_;
    }

    /**
     * @brief Check that the mask fits with a given reference genome, if given.
     */
    void check_masks_against_reference() const;

    /**
     * @brief Check that the provided per-sample masks and the input sample names
     * share their names, as otherwise, there is something off.
     */
    void check_sample_masks_name_list(
        std::vector<std::string> const& sample_names
    ) const;

    /**
     * @brief Create the tranform function to be applied to the stream for masking samples.
     */
    void add_sample_mask_transform_to_stream(
        genesis::population::VariantInputStream& stream
    ) const;

    /**
     * @brief Create the tranform function to be applied to the stream for masking the total.
     */
    void add_total_mask_transform_to_stream(
        genesis::population::VariantInputStream& stream
    ) const;

    // -------------------------------------------------------------------------
    //     Internal Run Functions
    // -------------------------------------------------------------------------

private:

    void prepare_sample_masks_() const;
    void prepare_total_mask_() const;

    void check_sample_masks_name_list_(
        std::vector<std::string> const& sample_names
    ) const;
    void check_reference_and_masks_compatibility_() const;
    void check_inter_masks_compatibility_() const;

    // -------------------------------------------------------------------------
    //     Option Members
    // -------------------------------------------------------------------------

private:

    // Filters for each sample individually
    CliOption<std::string> filter_mask_sample_bed_list_;
    CliOption<bool>        filter_mask_sample_bed_inv_ = false;
    CliOption<std::string> filter_mask_sample_fasta_list_;
    CliOption<size_t>      filter_mask_sample_fasta_min_ = 0;
    CliOption<bool>        filter_mask_sample_fasta_inv_ = false;

    // Filters for the whole variant
    CliOption<std::string> filter_mask_total_bed_;
    CliOption<bool>        filter_mask_total_bed_inv_ = false;
    CliOption<std::string> filter_mask_total_fasta_;
    CliOption<size_t>      filter_mask_total_fasta_min_ = 0;
    CliOption<bool>        filter_mask_total_fasta_inv_ = false;

    // We need a reference dict, for verification, and for inverting bed files
    VariantReferenceGenomeOptions const* ref_genome_opts_ = nullptr;

    // We keep the masks here. Bits that are set here are masked!
    mutable std::unordered_map<std::string, std::shared_ptr<GenomeLocusSet>> sample_masks_;
    mutable std::shared_ptr<GenomeLocusSet> total_mask_;

};

#endif // include guard
