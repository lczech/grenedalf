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

#include "options/variant_filter_mask.hpp"

#include "options/global.hpp"
#include "tools/misc.hpp"

#include "genesis/population/filter/variant_filter_positional.hpp"
#include "genesis/population/format/bed_reader.hpp"
#include "genesis/population/format/vcf_common.hpp"
#include "genesis/population/format/vcf_input_stream.hpp"
#include "genesis/population/function/functions.hpp"
#include "genesis/population/function/genome_locus_set.hpp"
#include "genesis/population/function/genome_locus_set.hpp"
#include "genesis/population/function/genome_region.hpp"
#include "genesis/population/genome_region.hpp"
#include "genesis/sequence/functions/dict.hpp"
#include "genesis/utils/io/input_source.hpp"

#include <algorithm>
#include <cassert>
#include <memory>
#include <stdexcept>

// =================================================================================================
//      Setup Functions
// =================================================================================================

void VariantFilterMaskOptions::add_mask_filter_opts_to_app(
    CLI::App* sub,
    std::string const& group
) {
    // Correct setup check.
    internal_check(
        filter_mask_bed_.option == nullptr,
        "Cannot use the same VariantFilterMaskOptions object multiple times."
    );

    // Add option for genomic mask filter by BED file.
    filter_mask_bed_.option = sub->add_option(
        "--filter-mask-bed",
        filter_mask_bed_.value,
        "Genomic positions to mask (mark as missing), as a BED file."
        "\nThe regions listed in the BED file are masked; this is in line with, e.g., smcpp, but "
        "is the inverse of the above usage of a BED file for selection regions, where instead "
        "the listed regions are kept. "
        "Note that this also conceptually differs from the region BED above. We here do not remove "
        "the masked positions, but instead just mark them as masked, so that they can still "
        "contribute to, e.g., denominators in the statistics for certain settings."
        "\nThis only uses the chromosome, as well as start and end information per line, and "
        "ignores everything else in the file. Note that BED uses 0-based positions, and a half-open "
        "`[)` interval for the end position; simply using columns extracted from other file formats "
        "(such as vcf or gff) will not work."

    );
    filter_mask_bed_.option->check( CLI::ExistingFile );
    filter_mask_bed_.option->group( group );

    // Add option for genomic mask filter by fasta-style mask file.
    filter_mask_fasta_.option = sub->add_option(
        "--filter-mask-fasta",
        filter_mask_fasta_.value,
        "Genomic positions to mask, as a FASTA-like mask file (such as used by vcftools). "
        "The file contains a sequence of integer digits `[0-9]`, one for each position on the "
        "chromosomes, which specify if the position should be masked or not. Any positions "
        "with digits above the `--filter-mask-fasta-min` value are tagged as being masked. "
        "Note that this conceptually differs from the region fasta above. We here do not remove the "
        "the masked positions, but instead just mark them as masked, so that they can still "
        "contribute to, e.g., denominators in the statistics for certain settings."
    );
    filter_mask_fasta_.option->check( CLI::ExistingFile );
    filter_mask_fasta_.option->group( group );

    // Add min threshold option for above mask option.
    filter_mask_fasta_min_.option = sub->add_option(
        "--filter-mask-fasta-min",
        filter_mask_fasta_min_.value,
        "When using `--filter-mask-fasta`, set the cutoff threshold for the masked digits. "
        "All positions above that value are masked. The default is 0, meaning that only exactly "
        "the positons with value 0 will not be masked."
    );
    filter_mask_fasta_min_.option->check(CLI::Range(0,9));
    filter_mask_fasta_min_.option->group( group );
    filter_mask_fasta_min_.option->needs( filter_mask_fasta_.option );

    // Add inversion option for above mask option.
    filter_mask_fasta_inv_.option = sub->add_flag(
        "--filter-mask-fasta-invert",
        filter_mask_fasta_inv_.value,
        "When using `--filter-mask-fasta`, invert the mask. This option has the same effect "
        "as the equivalent in vcftools, but instead of specifying the file, this here is a flag. "
        "When it is set, the mask specified above is inverted."
    );
    filter_mask_fasta_inv_.option->check(CLI::Range(0,9));
    filter_mask_fasta_inv_.option->group( group );
    filter_mask_fasta_inv_.option->needs( filter_mask_fasta_.option );

    // We only want to be able to specify a single mask file, as everything else would lead to chaos.
    filter_mask_bed_.option->excludes( filter_mask_fasta_.option );
    filter_mask_fasta_.option->excludes( filter_mask_bed_.option );
}

// =================================================================================================
//      Run Functions
// =================================================================================================

// -------------------------------------------------------------------------
//     prepare_mask
// -------------------------------------------------------------------------

void VariantFilterMaskOptions::prepare_mask() const
{
    using namespace genesis;
    using namespace genesis::population;
    using namespace genesis::utils;

    // Not running again if we already have set up a filter (i.e., if the shared pointer has data).
    if( mask_ ) {
        return;
    }

    // Add the mask from bed files.
    if( *filter_mask_bed_.option ) {
        LOG_MSG2 << "Reading mask BED file " << filter_mask_bed_.value;
        mask_ = std::make_shared<GenomeLocusSet>(
            BedReader().read_as_genome_locus_set( from_file( filter_mask_bed_.value ))
        );
    }

    // Add the mask from mask files.
    if( *filter_mask_fasta_.option ) {
        LOG_MSG2 << "Reading mask FASTA file " << filter_mask_fasta_.value;
        mask_ = std::make_shared<GenomeLocusSet>(
            read_mask_fasta(
                from_file( filter_mask_fasta_.value ),
                filter_mask_fasta_min_.value,
                filter_mask_fasta_inv_.value
            )
        );
    }
}

void VariantFilterMaskOptions::check_mask_against_reference(
    std::shared_ptr<genesis::sequence::SequenceDict> ref_dict
) const {
    // This only makes sense for the fixed length mask provided by a fasta-like format.
    if( ! *filter_mask_fasta_.option || ! ref_dict ) {
        return;
    }

    // Turn the mask into a dict for the comparison.
    auto const mask_dict = genesis::population::reference_locus_set_to_dict( *get_mask() );

    // Strict comparison: ref genome and mask need to have the same sizes exactly.
    if( ! genesis::sequence::compatible_references( *ref_dict, mask_dict )) {
        throw std::invalid_argument(
            "Provided reference genome and mask filter do not contain the exact same "
            "chromosomes or lengths of chromosomes. We fail here out of caution, "
            "as this typically indicates some inconsistency with the data."
        );
    }
}

std::function<void( genesis::population::Variant& )> VariantFilterMaskOptions::make_mask_transform() const
{
    // Edge case if no mask is given. Should not occur, as we usually do check anyway
    // before calling this function, but doesn't hurt to have the check just in case for later.
    if( ! get_mask() ) {
        return {};
    }

    // We make a tagging positional filter, using the mask tag.
    // We need to use the complement here, as the filter expects the "good" positions to be set.
    return make_variant_filter_by_region_tagging(
        get_mask(), genesis::population::VariantFilterTag::kMaskedPosition, true
    );
}
