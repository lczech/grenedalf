/*
    grenedalf - Genome Analyses of Differential Allele Frequencies
    Copyright (C) 2020-2022 Lucas Czech

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

#include "options/window.hpp"

#include "options/global.hpp"
#include "tools/misc.hpp"

#include "genesis/population/formats/bed_reader.hpp"
#include "genesis/population/formats/genome_region_reader.hpp"
#include "genesis/population/formats/gff_reader.hpp"
#include "genesis/population/formats/map_bim_reader.hpp"
#include "genesis/population/functions/functions.hpp"
#include "genesis/population/functions/genome_region.hpp"
#include "genesis/population/genome_region.hpp"
#include "genesis/utils/core/std.hpp"

#include <cassert>
#include <memory>
#include <stdexcept>

// =================================================================================================
//      Setup Functions
// =================================================================================================

void WindowOptions::add_window_opts_to_app(
    CLI::App* sub,
    std::string const& group
) {
    // Add the set combination of the genom regions.
    window_type_.option = sub->add_option(
        "--window-type",
        window_type_.value,
        "Type of window to use: sliding windows over intervals along the genome, or "
        "windows corresponding to some regions of interest."
    );
    window_type_.option->transform(
        CLI::IsMember(
            { "sliding", "regions" },
            CLI::ignore_case
        )
    );
    window_type_.option->group( group );
    window_type_.option->required();

    // -------------------------------------------------------------------------
    //     Sliding interval window
    // -------------------------------------------------------------------------

    // Width
    window_width_.option = sub->add_option(
        "--window-sliding-width",
        window_width_.value,
        "When using `" + window_type_.option->get_name() + " sliding`: "
        "Width of each window along the chromosome. Has to be provided when using sliding window."
    );
    window_width_.option->group( group );
    // window_width_.option->required();

    // Stride
    window_stride_.option = sub->add_option(
        "--window-sliding-stride",
        window_stride_.value,
        "When using `" + window_type_.option->get_name() + " sliding`: "
        "Stride between windows along the chromosome, that is how far to move to get to the next "
        "window. If set to 0 (default), this is set to the same value as the `--window-width`."
    );
    window_stride_.option->group( group );

    // -------------------------------------------------------------------------
    //     Region window
    // -------------------------------------------------------------------------

    // Add option for genomic region list.
    region_.option = sub->add_option(
        "--window-region",
        region_.value,
        "When using `" + window_type_.option->get_name() + " regions`: "
        "Genomic region to process as windows, in the format \"chr\" (for whole chromosomes), "
        "\"chr:position\", \"chr:start-end\", or \"chr:start..end\". "
        "Positions are 1-based and inclusive (closed intervals). "
        "The option can be provided multiple times to add region windows to be processed."
    );
    region_.option->group( group );

    // Add option for genomic region window.
    region_list_.option = sub->add_option(
        "--window-region-list",
        region_list_.value,
        "When using `" + window_type_.option->get_name() + " regions`: "
        "Genomic regions to process as windows, as a file with one region per line, "
        "either in the format \"chr\" (for whole chromosomes), \"chr:position\", \"chr:start-end\", "
        "\"chr:start..end\", or tab- or space-delimited \"chr position\" or \"chr start end\". "
        "Positions are 1-based and inclusive (closed intervals). "
        "The option can be provided multiple times to add region windows to be processed."
    );
    region_list_.option->check( CLI::ExistingFile );
    region_list_.option->group( group );

    // Add option for genomic region window by BED file.
    region_bed_.option = sub->add_option(
        "--window-region-bed",
        region_bed_.value,
        "When using `" + window_type_.option->get_name() + " regions`: "
        "Genomic regions to process as windows, as a BED file. This only uses the chromosome, "
        "as well as start and end information per line, and ignores everything else in the file. "
        "Note that BED uses 0-based positions, and a half-open `[)` interval for the end position; "
        "simply using columns extracted from other file formats (such as vcf or gff) will not work. "
        "The option can be provided multiple times to add region windows to be processed."
    );
    region_bed_.option->check( CLI::ExistingFile );
    region_bed_.option->group( group );

    // Add option for genomic region window by GFF/GTF file.
    region_gff_.option = sub->add_option(
        "--window-region-gff",
        region_gff_.value,
        "When using `" + window_type_.option->get_name() + " regions`: "
        "Genomic regions to process as windows, as a GFF2/GFF3/GTF file. This only uses the chromosome, "
        "as well as start and end information per line, and ignores everything else in the file. "
        "The option can be provided multiple times to add region windows to be processed."
    );
    region_gff_.option->check( CLI::ExistingFile );
    region_gff_.option->group( group );

    // // MAP and BIM files are specifying coordinates, and not regions, so we leave them out for now.
    // // Add option for genomic region window by PLINK MAP/BIM file.
    // region_bim_.option = sub->add_option(
    //     "--window-region-map-bim",
    //     region_bim_.value,
    //     "When using `" + window_type_.option->get_name() + " regions`: "
    //     "Genomic regions to process as windows, as a MAP or BIM file as used in PLINK. This only "
    //     "uses the chromosome and coordinate per line, and ignores everything else in the file. "
    //     "Consecutive positions are merged into one interval. "
    //     "The option can be provided multiple times to add region windows to be processed."
    // );
    // region_bim_.option->check( CLI::ExistingFile );
    // region_bim_.option->group( group );

    // // VCF files are specifying coordinates, and not regions, so we leave them out for now.
    // // Add option for genomic region window by VCF file.
    // region_vcf_.option = sub->add_option(
    //     "--window-region-vcf",
    //     region_vcf_.value,
    //     "When using `" + window_type_.option->get_name() + " regions`: "
    //     "Genomic regions to process as windows, as a VCF/BCF file (such as a known-variants file). "
    //     "This only uses the chromosome and coordinate per line, and ignores everything else in the file. "
    //     "Consecutive positions are merged into one interval. "
    //     "The option can be provided multiple times to add region windows to be processed."
    // );
    // region_vcf_.option->check( CLI::ExistingFile );
    // region_vcf_.option->group( group );

    // Skip empty regions
    skip_empty_regions_.option = sub->add_flag(
        "--window-region-skip-empty",
        skip_empty_regions_.value,
        "When using `" + window_type_.option->get_name() + " regions`: "
        "In cases where there is no data in the input files for a region window, by default, "
        "we produce some \"empty\" or NaN output. "
        "With this option however, regions without data are skipped in the output."
    );
    skip_empty_regions_.option->group( group );
}

// =================================================================================================
//      Run Functions
// =================================================================================================

// std::pair<size_t, size_t> WindowOptions::get_window_width_and_stride() const
// {
//     // We need to check the stride here ourselves. For the actual window iterator, this is done
//     // in the iterator constructor, but when this function here is called, that might not have
//     // happened yet, and so we need to do the check ourselves.
//     if( window_stride_.value == 0 ) {
//         // stride == 0 --> use stride == width
//         return { window_width_.value, window_width_.value };
//     } else {
//         // stride != 0 --> just use it as is
//         return { window_width_.value, window_stride_.value };
//     }
// }

// -------------------------------------------------------------------------
//     get_variant_window_iterator
// -------------------------------------------------------------------------

std::unique_ptr<VariantWindowIterator> WindowOptions::get_variant_window_iterator(
    genesis::population::VariantInputIterator& input
) const {

    if( window_type_.value == "sliding" ) {
        return get_variant_window_iterator_sliding_interval_( input );
    } else if( window_type_.value == "regions" ) {
        return get_variant_window_iterator_regions_( input );
    } else {
        throw CLI::ValidationError(
            window_type_.option->get_name(),
            "Invalid window type '" + window_type_.value + "'."
        );
    }

    // Cannot happen - just here for compiler satisfaction.
    assert( false );
    return nullptr;
}

// -------------------------------------------------------------------------
//     get_variant_window_iterator_sliding_interval_
// -------------------------------------------------------------------------

std::unique_ptr<VariantWindowIterator> WindowOptions::get_variant_window_iterator_sliding_interval_(
    genesis::population::VariantInputIterator& input
) const {
    if( window_width_.value == 0 ) {
        throw CLI::ValidationError(
            window_width_.option->get_name(),
            "Invalid window width for sliding window, has to be greater than 0."
        );
    }

    return genesis::utils::make_unique<VariantSlidingIntervalWindowIterator>(
        genesis::population::make_default_sliding_interval_window_iterator(
            input.begin(), input.end(), window_width_.value, window_stride_.value
        )
    );
}

// -------------------------------------------------------------------------
//     get_variant_window_iterator_regions_
// -------------------------------------------------------------------------

std::unique_ptr<VariantWindowIterator> WindowOptions::get_variant_window_iterator_regions_(
    genesis::population::VariantInputIterator& input
) const {
    using namespace genesis::population;
    using namespace genesis::utils;

    // Create a region list to which we add all provided regions.
    auto region_list = std::make_shared<GenomeRegionList>();
    size_t in_regions = 0;

    // Add the regions from strings.
    for( auto const& value : region_.value ) {
        region_list->add( parse_genome_region( value ));
        ++in_regions;
    }

    // Add the region list files.
    for( auto const& list_file : region_list_.value ) {
        LOG_MSG2 << "Reading regions list file " << list_file;
        GenomeRegionReader().read_as_genome_region_list( from_file( list_file ), *region_list );
        ++in_regions;
    }

    // Add the regions from bed files.
    for( auto const& file : region_bed_.value ) {
        LOG_MSG2 << "Reading regions BED file " << file;
        BedReader().read_as_genome_region_list( from_file( file ), *region_list );
        ++in_regions;
    }

    // Add the regions from gff files.
    for( auto const& file : region_gff_.value ) {
        LOG_MSG2 << "Reading regions GFF2/GFF3/GTF file " << file;
        GffReader().read_as_genome_region_list( from_file( file ), *region_list );
        ++in_regions;
    }

    // Add the regions from map/bim files.
    for( auto const& file : region_bim_.value ) {
        LOG_MSG2 << "Reading regions MAP/BIM file " << file;
        MapBimReader().read_as_genome_region_list( from_file( file ), *region_list );
        ++in_regions;
    }

    // Add the regions from vcf files.
    for( auto const& file : region_vcf_.value ) {
        LOG_MSG2 << "Reading regions VCF/BCF file " << file;
        genome_region_list_from_vcf_file( file, *region_list );
        ++in_regions;
    }

    // User error. We only check that some input regions were provided, but not whether
    // they actually contain anything... That would be on the user.
    if( in_regions == 0) {
        throw CLI::ValidationError(
            window_type_.option->get_name(),
            "Window type '" + window_type_.value + "' specified, but no region inputs are given. "
            "At least one type of input containing regions to process has to be provided."
        );
    }

    // Count how many regions we have, in total, for user convenience.
    LOG_MSG << "Iterating over " << region_list->total_region_count() << " total regions across "
            << region_list->chromosome_count() << " chromosomes.";
    for( auto const& chr : region_list->chromosome_names() ) {
        LOG_MSG2 << " - Chromosome \"" << chr << "\" with "
                 << region_list->region_count( chr ) << " regions";
    }

    // Return an iterator over the regions that we just read.
    // It copies the shared pointer to the regions, keeping it alive.
    return genesis::utils::make_unique<VariantRegionWindowIterator>(
        genesis::population::make_default_region_window_iterator(
            input.begin(), input.end(), region_list
        )
    );
}
