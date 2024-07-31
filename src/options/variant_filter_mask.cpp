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
#include "genesis/population/filter/sample_counts_filter_positional.hpp"
#include "genesis/population/format/bed_reader.hpp"
#include "genesis/population/format/vcf_common.hpp"
#include "genesis/population/format/vcf_input_stream.hpp"
#include "genesis/population/function/functions.hpp"
#include "genesis/population/function/genome_locus_set.hpp"
#include "genesis/population/function/genome_locus_set.hpp"
#include "genesis/population/function/genome_region.hpp"
#include "genesis/population/genome_region.hpp"
#include "genesis/sequence/functions/dict.hpp"
#include "genesis/utils/core/fs.hpp"
#include "genesis/utils/io/input_source.hpp"
#include "genesis/utils/text/convert.hpp"
#include "genesis/utils/text/string.hpp"

#include <algorithm>
#include <cassert>
#include <memory>
#include <stdexcept>
#include <unordered_set>

// =================================================================================================
//      Setup Functions
// =================================================================================================

void VariantFilterMaskOptions::add_mask_filter_opts_to_app(
    CLI::App* sub,
    std::string const& group
) {
    add_mask_filter_sample_opts_to_app( sub, group );
    add_mask_filter_total_opts_to_app( sub, group );
}

void VariantFilterMaskOptions::add_mask_filter_sample_opts_to_app(
    CLI::App* sub,
    std::string const& group
) {
    // Correct setup check.
    internal_check(
        filter_mask_sample_bed_list_.option == nullptr,
        "Cannot use the same VariantFilterMaskOptions object multiple times."
    );

    // Add option for genomic mask filter by BED file.
    filter_mask_sample_bed_list_.option = sub->add_option(
        "--filter-mask-samples-bed-list",
        filter_mask_sample_bed_list_.value,
        "For each sample, genomic positions to mask (mark as missing), as a set of BED files."
        "\nSee the below `--filter-mask-total-bed` for details. Here, individual BED files can "
        "be provided for each sample, for fine-grained control over the masking. The option takes "
        "a path to a file that contains a comma- or tab-separated list of sample names and "
        "BED file paths, with one name/path pair per line, in any order of lines."

    );
    filter_mask_sample_bed_list_.option->check( CLI::ExistingFile );
    filter_mask_sample_bed_list_.option->group( group );

    // Add option for genomic mask filter by fasta-style mask file.
    filter_mask_sample_fasta_list_.option = sub->add_option(
        "--filter-mask-samples-fasta-list",
        filter_mask_sample_fasta_list_.value,
        "For each sample, genomic positions to mask, as a FASTA-like mask file.\n"
        "\nSee the below `--filter-mask-total-fasta` for details. Here, individual FASTA files can "
        "be provided for each sample, for fine-grained control over the masking. The option takes "
        "a path to a file that contains a comma- or tab-separated list of sample names and "
        "FASTA file paths, with one name/path pair per line, in any order of lines."
    );
    filter_mask_sample_fasta_list_.option->check( CLI::ExistingFile );
    filter_mask_sample_fasta_list_.option->group( group );

    // Add min threshold option for above mask option.
    filter_mask_sample_fasta_min_.option = sub->add_option(
        "--filter-mask-samples-fasta-min",
        filter_mask_sample_fasta_min_.value,
        "When using `--filter-mask-samples-fasta-list`, set the cutoff threshold for the masked "
        "digits. All positions above that value are masked. The default is 0, meaning that only "
        "exactly the positons with value 0 will not be masked."
    );
    filter_mask_sample_fasta_min_.option->check(CLI::Range(0,9));
    filter_mask_sample_fasta_min_.option->group( group );
    filter_mask_sample_fasta_min_.option->needs( filter_mask_sample_fasta_list_.option );

    // Add inversion option for above mask option.
    filter_mask_sample_fasta_inv_.option = sub->add_flag(
        "--filter-mask-samples-fasta-invert",
        filter_mask_sample_fasta_inv_.value,
        "When using `--filter-mask-samples-fasta-list`, invert the mask."
        "When this flag is set, the mask specified above is inverted."
    );
    filter_mask_sample_fasta_inv_.option->check(CLI::Range(0,9));
    filter_mask_sample_fasta_inv_.option->group( group );
    filter_mask_sample_fasta_inv_.option->needs( filter_mask_sample_fasta_list_.option );

    // We only want to be able to specify a single mask file, as everything else would lead to chaos.
    filter_mask_sample_bed_list_.option->excludes( filter_mask_sample_fasta_list_.option );
    filter_mask_sample_fasta_list_.option->excludes( filter_mask_sample_bed_list_.option );
}

void VariantFilterMaskOptions::add_mask_filter_total_opts_to_app(
    CLI::App* sub,
    std::string const& group
) {
    // Correct setup check.
    internal_check(
        filter_mask_total_bed_.option == nullptr,
        "Cannot use the same VariantFilterMaskOptions object multiple times."
    );

    // Add option for genomic mask filter by BED file.
    filter_mask_total_bed_.option = sub->add_option(
        "--filter-mask-total-bed",
        filter_mask_total_bed_.value,
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
    filter_mask_total_bed_.option->check( CLI::ExistingFile );
    filter_mask_total_bed_.option->group( group );

    // Add option for genomic mask filter by fasta-style mask file.
    filter_mask_total_fasta_.option = sub->add_option(
        "--filter-mask-total-fasta",
        filter_mask_total_fasta_.value,
        "Genomic positions to mask, as a FASTA-like mask file (such as used by vcftools).\n"
        "The file contains a sequence of integer digits `[0-9]`, one for each position on the "
        "chromosomes, which specify if the position should be masked or not. Any positions "
        "with digits above the `--filter-mask-total-fasta-min` value are tagged as being masked. "
        "Note that this conceptually differs from the region fasta above. We here do not remove the "
        "the masked positions, but instead just mark them as masked, so that they can still "
        "contribute to, e.g., denominators in the statistics for certain settings."
    );
    filter_mask_total_fasta_.option->check( CLI::ExistingFile );
    filter_mask_total_fasta_.option->group( group );

    // Add min threshold option for above mask option.
    filter_mask_total_fasta_min_.option = sub->add_option(
        "--filter-mask-total-fasta-min",
        filter_mask_total_fasta_min_.value,
        "When using `--filter-mask-total-fasta`, set the cutoff threshold for the masked digits. "
        "All positions above that value are masked. The default is 0, meaning that only exactly "
        "the positons with value 0 will not be masked."
    );
    filter_mask_total_fasta_min_.option->check(CLI::Range(0,9));
    filter_mask_total_fasta_min_.option->group( group );
    filter_mask_total_fasta_min_.option->needs( filter_mask_total_fasta_.option );

    // Add inversion option for above mask option.
    filter_mask_total_fasta_inv_.option = sub->add_flag(
        "--filter-mask-total-fasta-invert",
        filter_mask_total_fasta_inv_.value,
        "When using `--filter-mask-total-fasta`, invert the mask. This option has the same effect "
        "as the equivalent in vcftools, but instead of specifying the file, this here is a flag. "
        "When it is set, the mask specified above is inverted."
    );
    filter_mask_total_fasta_inv_.option->check(CLI::Range(0,9));
    filter_mask_total_fasta_inv_.option->group( group );
    filter_mask_total_fasta_inv_.option->needs( filter_mask_total_fasta_.option );

    // We only want to be able to specify a single mask file, as everything else would lead to chaos.
    filter_mask_total_bed_.option->excludes( filter_mask_total_fasta_.option );
    filter_mask_total_fasta_.option->excludes( filter_mask_total_bed_.option );
}

// =================================================================================================
//      Run Functions
// =================================================================================================

// -------------------------------------------------------------------------
//     prepare_masks
// -------------------------------------------------------------------------

void VariantFilterMaskOptions::prepare_masks() const
{
    prepare_sample_masks_();
    prepare_total_mask_();
}

// -------------------------------------------------------------------------
//     check_masks_against_reference
// -------------------------------------------------------------------------

void VariantFilterMaskOptions::check_masks_against_reference(
    std::shared_ptr<genesis::sequence::SequenceDict> ref_dict
) const {
    check_reference_and_masks_compatibility_( ref_dict );
    check_inter_masks_compatibility_();
}

// -------------------------------------------------------------------------
//     check_sample_masks_name_list
// -------------------------------------------------------------------------

void VariantFilterMaskOptions::check_sample_masks_name_list(
    std::vector<std::string> const& sample_names
) const {
    check_sample_masks_name_list_( sample_names );
}

// -------------------------------------------------------------------------
//     add_sample_mask_transform_to_stream
// -------------------------------------------------------------------------

void VariantFilterMaskOptions::add_sample_mask_transform_to_stream(
    genesis::population::VariantInputStream& stream
) const {
    // Edge case if no sample masks are given. Should not occur, as we usually do check anyway
    // before calling this function, but doesn't hurt to have the check just in case for later.
    prepare_sample_masks_();
    if( sample_masks_.empty() ) {
        return;
    }

    // Get the sample masks for the samples of the stream, in the correct order.
    // If there is no mask for a given sample, we just add a default constructed shared ptr,
    // which is recognized by the transformation to not do anything.
    // We warn about this in check_sample_masks_name_list_() already, so no need to warn here.
    auto sample_masks_subset_ = std::vector<std::shared_ptr<GenomeLocusSet>>{};
    sample_masks_subset_.reserve( stream.data().sample_names.size() );
    for( auto const& sample_name : stream.data().sample_names ) {
        if( sample_masks_.count( sample_name ) > 0 ) {
            sample_masks_subset_.push_back( sample_masks_[sample_name] );
        } else {
            sample_masks_subset_.emplace_back();
        }
    }
    assert( sample_masks_subset_.size() == stream.data().sample_names.size() );

    // We make a tagging positional filter, using the mask tag.
    // We need to use the complement here, as the filter expects the "good" positions to be set.
    stream.add_transform( make_sample_counts_filter_by_region_tagging(
        sample_masks_subset_, genesis::population::SampleCountsFilterTag::kMaskedPosition, true
    ));
}

// -------------------------------------------------------------------------
//     add_total_mask_transform_to_stream
// -------------------------------------------------------------------------

void VariantFilterMaskOptions::add_total_mask_transform_to_stream(
    genesis::population::VariantInputStream& stream
) const {
    // Edge case if no mask is given. Should not occur, as we usually do check anyway
    // before calling this function, but doesn't hurt to have the check just in case for later.
    prepare_total_mask_();
    if( ! total_mask_ ) {
        return;
    }

    // We make a tagging positional filter, using the mask tag.
    // We need to use the complement here, as the filter expects the "good" positions to be set.
    stream.add_transform( make_variant_filter_by_region_tagging(
        total_mask_, genesis::population::VariantFilterTag::kMaskedPosition, true
    ));
}

// =================================================================================================
//      Internal Run Functions
// =================================================================================================

// -------------------------------------------------------------------------
//     prepare_sample_masks_
// -------------------------------------------------------------------------

void VariantFilterMaskOptions::prepare_sample_masks_() const
{
    using namespace genesis;
    using namespace genesis::population;
    using namespace genesis::utils;

    // Not running again if we already have set up a filter (i.e., if the shared pointer has data).
    if( ! sample_masks_.empty() ) {
        return;
    }

    auto get_list_file_pairs_ = []( CliOption<std::string> const& opt )
    {
        // Read the file line by line and process.
        auto const lines = file_read_lines( opt.value );
        std::unordered_map<std::string, std::string> sample_to_file;
        for( size_t i = 0; i < lines.size(); ++i ) {
            auto const& line = lines[i];

            // Dissect the line and see if we got a sample name and a number.
            auto const pair = split( line, ",\t", false );
            if( pair.size() != 2 ) {
                throw CLI::ValidationError(
                    opt.option->get_name() + "(" +
                    opt.value + ")",
                    "Invalid samples maks file that contains an invalid line at " +
                    std::to_string( i + 1 ) + " not consisting of a sample name and a mask file."
                );
            }
            if( ! file_exists( pair[1] )) {
                throw CLI::ValidationError(
                    opt.option->get_name() + "(" +
                    opt.value + ")",
                    "Invalid samples maks file that contains an invalid line at " +
                    std::to_string( i + 1 ) + " with a non-existing mask file `" + pair[1] + "`"
                );
            }

            // For name value pairs, do a duplicate check first...
            if( sample_to_file.count( pair[0] ) > 0 ) {
                throw CLI::ValidationError(
                    opt.option->get_name() + "(" +
                    opt.value + ")",
                    "Invalid line that contains duplicate sample names (line " +
                    std::to_string( i + 1 ) + "): \"" + pair[0] + "\""
                );
            }

            // ... then add the entry to the map.
            assert( sample_to_file.count( pair[0] ) == 0 );
            sample_to_file[ pair[0] ] = pair[1];
        }
        return sample_to_file;
    };

    // Add the masks from bed files.
    if( *filter_mask_sample_bed_list_.option ) {
        LOG_MSG2 << "Reading sample masks BED files " << filter_mask_sample_bed_list_.value;
        auto const sample_to_file = get_list_file_pairs_( filter_mask_sample_bed_list_ );
        for( auto const& sample_file_pair : sample_to_file ) {
            LOG_MSG2 << "  - " << sample_file_pair.first << ": " << sample_file_pair.second;
            assert( sample_masks_( sample_file_pair.first ).count() == 0 );
            sample_masks_[ sample_file_pair.first ] = std::make_shared<GenomeLocusSet>(
                BedReader().read_as_genome_locus_set( from_file( sample_file_pair.second ))
            );
        }
    }

    // Add the mask from mask files.
    if( *filter_mask_sample_fasta_list_.option ) {
        LOG_MSG2 << "Reading sample masks FASTA files " << filter_mask_sample_fasta_list_.value;
        auto const sample_to_file = get_list_file_pairs_( filter_mask_sample_fasta_list_ );
        for( auto const& sample_file_pair : sample_to_file ) {
            LOG_MSG2 << "  - " << sample_file_pair.first << ": " << sample_file_pair.second;
            assert( sample_masks_( sample_file_pair.first ).count() == 0 );
            sample_masks_[ sample_file_pair.first ] = std::make_shared<GenomeLocusSet>(
                read_mask_fasta(
                    from_file( sample_file_pair.second ),
                    filter_mask_sample_fasta_min_.value,
                    filter_mask_sample_fasta_inv_.value
                )
            );
        }
    }
}

// -------------------------------------------------------------------------
//     prepare_total_mask_
// -------------------------------------------------------------------------

void VariantFilterMaskOptions::prepare_total_mask_() const
{
    using namespace genesis;
    using namespace genesis::population;
    using namespace genesis::utils;

    // Not running again if we already have set up a filter (i.e., if the shared pointer has data).
    if( total_mask_ ) {
        return;
    }

    // Add the mask from a bed file.
    if( *filter_mask_total_bed_.option ) {
        LOG_MSG2 << "Reading mask BED file " << filter_mask_total_bed_.value;
        total_mask_ = std::make_shared<GenomeLocusSet>(
            BedReader().read_as_genome_locus_set( from_file( filter_mask_total_bed_.value ))
        );
    }

    // Add the mask from a fasta file.
    if( *filter_mask_total_fasta_.option ) {
        LOG_MSG2 << "Reading mask FASTA file " << filter_mask_total_fasta_.value;
        total_mask_ = std::make_shared<GenomeLocusSet>(
            read_mask_fasta(
                from_file( filter_mask_total_fasta_.value ),
                filter_mask_total_fasta_min_.value,
                filter_mask_total_fasta_inv_.value
            )
        );
    }
}

// -------------------------------------------------------------------------
//     check_sample_masks_name_list_
// -------------------------------------------------------------------------

void VariantFilterMaskOptions::check_sample_masks_name_list_(
    std::vector<std::string> const& sample_names
) const {
    // If there are no samples masks, we have nothing to do here.
    prepare_sample_masks_();
    if( sample_masks_.empty() ) {
        return;
    }

    // Make a set of the input sample names, as that's faster to work with here.
    auto const sample_name_set = std::unordered_set<std::string>{
        sample_names.begin(), sample_names.end()
    };

    // We compare both ways, and warn about missing names in either list in different ways.
    // We could use some algorithms here such as std::set_intersection(), but it seems
    // that would not make it much easier in this case, as we'd have to sort the elements first.
    // So instead we do it our way.

    // First, get all sample names that are in the masks, but not in the input file samples.
    std::unordered_set<std::string> mask_extra_samples;
    for( auto const& sample_mask : sample_masks_ ) {
        if( sample_name_set.count( sample_mask.first ) == 0 ) {
            mask_extra_samples.insert( sample_mask.first );
        }
    }
    if( mask_extra_samples.size() > 0 ) {
        LOG_WARN << "Provided per-sample masks contain masks for sample names "
                 << "that are not present in the input:";
        for( auto const& extra_sample : mask_extra_samples ) {
            LOG_WARN << " - " << extra_sample;
        }
        LOG_WARN << "We are continuing now, but those masks will not be used. "
                 << "If this is intentional, this warning can be ignored.";
    }

    // Second, the other way round, get all sample names that are in the input, but not in the masks.
    std::unordered_set<std::string> input_extra_samples;
    for( auto const& input_sample : sample_name_set ) {
        if( sample_masks_.count( input_sample ) == 0 ) {
            input_extra_samples.insert( input_sample );
        }
    }
    if( input_extra_samples.size() > 0 ) {
        LOG_WARN << "Provided per-sample masks do not contain masks for sample names "
                 << "that are present in the input:";
        for( auto const& extra_sample : input_extra_samples ) {
            LOG_WARN << " - " << extra_sample;
        }
        LOG_WARN << "We are continuing now, but those samples will not be masked. "
                 << "If this is intentional, this warning can be ignored.";
    }
}

// -------------------------------------------------------------------------
//     check_reference_and_masks_compatibility_
// -------------------------------------------------------------------------

void VariantFilterMaskOptions::check_reference_and_masks_compatibility_(
    std::shared_ptr<genesis::sequence::SequenceDict> ref_dict
) const {
    using namespace genesis::sequence;

    // See if there is a ref dict at all to compare to.
    if( !ref_dict ) {
        return;
    }

    // The below checks only make sense for fixed length masks provided by a fasta-like format.
    // For the bed format masks, we do not execute these checks, as those are not fixed length.

    // Check the total mask against the reference dict
    if( *filter_mask_total_fasta_.option ) {
        // Turn the mask into a dict for the comparison.
        prepare_total_mask_();
        assert( total_mask_ );
        auto const mask_dict = genesis::population::reference_locus_set_to_dict( *total_mask_ );

        // Strict comparison: ref genome and mask need to have the same sizes exactly.
        if( ! compatible_references( *ref_dict, mask_dict, ReferenceComparisonMode::kStrict )) {
            throw std::invalid_argument(
                "Provided reference genome and mask file (" + filter_mask_total_fasta_.value +
                ") do not contain the exact same chromosomes or lengths of chromosomes. "
                "We fail here out of caution, as this typically indicates some inconsistency "
                "with the data."
            );
        }
    }

    // Check the samples masks against the reference dict
    if( *filter_mask_sample_fasta_list_.option ) {
        prepare_sample_masks_();
        for( auto const& sample_mask : sample_masks_ ) {
            // Turn the mask into a dict for the comparison.
            auto const mask_dict = genesis::population::reference_locus_set_to_dict(
                *sample_mask.second
            );

            // Strict comparison: ref genome and mask need to have the same sizes exactly.
            if( ! compatible_references( *ref_dict, mask_dict, ReferenceComparisonMode::kStrict )) {
                throw std::invalid_argument(
                    "Provided reference genome and sample mask file for sample " +
                    sample_mask.first + " do not contain the exact same chromosomes or lengths "
                    "of chromosomes. We fail here out of caution, as this typically indicates "
                    "some inconsistency with the data."
                );
            }
        }
    }
}

// -------------------------------------------------------------------------
//     check_inter_masks_compatibility_
// -------------------------------------------------------------------------

void VariantFilterMaskOptions::check_inter_masks_compatibility_() const
{
    using namespace genesis::population;
    using namespace genesis::sequence;

    // Make sure we have read the samples masks.
    prepare_sample_masks_();
    prepare_total_mask_();

    // If both a total and sample masks are given, we check those.
    if( total_mask_ ) {
        // Turn the total mask into a dict for the comparison.
        auto const mask_dict = reference_locus_set_to_dict( *total_mask_ );

        // Same for each sample masks, and then run the comparison.
        for( auto const& sample_mask : sample_masks_ ) {
            auto const sample_dict = reference_locus_set_to_dict( *sample_mask.second );

            // Strict comparison: masks need to have the same sizes exactly.
            auto const compatible = compatible_references(
                sample_dict, mask_dict, ReferenceComparisonMode::kStrict
            );
            if( !compatible ) {
                throw std::invalid_argument(
                    "Provided total mask file and sample mask for sample " +
                    sample_mask.first + " do not contain the exact same chromosomes or lengths "
                    "of chromosomes. We fail here out of caution, as this typically indicates "
                    "some inconsistency with the data."
                );
            }
        }

        // We have checked the total against the sample masks. If that worked, all the sample masks
        // are also compatible with each other, so we do not need to check that again below.
        return;
    }

    // Alternatively, we check that the given sample masks are compatible with each other,
    // to avoid weird bugs down the line. Only needed if there is more than one.
    if( sample_masks_.size() < 2 ) {
        return;
    }

    // Turn the first sample mask into a dict for the comparison.
    // We just need any of them, and as long as they all are compatible with that one,
    // they are with each other. So we just use the first here.
    auto const first_mask_dict = reference_locus_set_to_dict( *sample_masks_.begin()->second );

    // Same for each sample masks, and then run the comparison.
    for( auto const& sample_mask : sample_masks_ ) {
        auto const sample_dict = reference_locus_set_to_dict( *sample_mask.second );

        // Strict comparison: masks need to have the same sizes exactly.
        auto const compatible = compatible_references(
            sample_dict, first_mask_dict, ReferenceComparisonMode::kStrict
        );
        if( !compatible ) {
            throw std::invalid_argument(
                "Provided sample mask files for samples " + sample_masks_.begin()->first +
                " and " + sample_mask.first + " do not contain the exact same chromosomes or "
                "lengths of chromosomes. We fail here out of caution, as this typically indicates "
                "some inconsistency with the data."
            );
        }
    }
}
