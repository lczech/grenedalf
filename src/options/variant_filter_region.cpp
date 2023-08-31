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

#include "options/variant_filter_region.hpp"

#include "options/global.hpp"
#include "tools/misc.hpp"

#include "genesis/population/formats/bed_reader.hpp"
#include "genesis/population/formats/genome_region_reader.hpp"
#include "genesis/population/formats/gff_reader.hpp"
#include "genesis/population/formats/map_bim_reader.hpp"
#include "genesis/population/formats/vcf_input_iterator.hpp"
#include "genesis/population/functions/filter_transform.hpp"
#include "genesis/population/functions/functions.hpp"
#include "genesis/population/functions/genome_region.hpp"
#include "genesis/population/genome_region.hpp"

#include <algorithm>
#include <cassert>
#include <memory>
#include <stdexcept>

// =================================================================================================
//      Setup Functions
// =================================================================================================

void VariantFilterRegionOptions::add_region_filter_opts_to_app(
    CLI::App* sub,
    std::string const& group
) {
    // Correct setup check.
    internal_check(
        filter_region_.option == nullptr,
        "Cannot use the same VariantFilterRegionOptions object multiple times."
    );

    // Add option for genomic region filter.
    filter_region_.option = sub->add_option(
        "--filter-region",
        filter_region_.value,
        "Genomic region to filter for, in the format \"chr\" (for whole chromosomes), "
        "\"chr:position\", \"chr:start-end\", or \"chr:start..end\". "
        "Positions are 1-based and inclusive (closed intervals). "
        "Multiple region options can be provided, see also `--filter-region-set`."
    );
    filter_region_.option->group( group );

    // Add option for genomic region filter.
    filter_region_list_.option = sub->add_option(
        "--filter-region-list",
        filter_region_list_.value,
        "Genomic regions to filter for, as a file with one region per line, "
        "either in the format \"chr\" (for whole chromosomes), \"chr:position\", \"chr:start-end\", "
        "\"chr:start..end\", or tab- or space-delimited \"chr position\" or \"chr start end\". "
        "Positions are 1-based and inclusive (closed intervals). "
        "Multiple region options can be provided, see also `--filter-region-set`."
    );
    filter_region_list_.option->check( CLI::ExistingFile );
    filter_region_list_.option->group( group );

    // Add option for genomic region filter by BED file.
    filter_region_bed_.option = sub->add_option(
        "--filter-region-bed",
        filter_region_bed_.value,
        "Genomic regions to filter for, as a BED file. This only uses the chromosome, "
        "as well as start and end information per line, and ignores everything else in the file. "
        "Note that BED uses 0-based positions, and a half-open `[)` interval for the end position; "
        "simply using columns extracted from other file formats (such as vcf or gff) will not work. "
        "Multiple region options can be provided, see also `--filter-region-set`."
    );
    filter_region_bed_.option->check( CLI::ExistingFile );
    filter_region_bed_.option->group( group );

    // Add option for genomic region filter by GFF/GTF file.
    filter_region_gff_.option = sub->add_option(
        "--filter-region-gff",
        filter_region_gff_.value,
        "Genomic regions to filter for, as a GFF2/GFF3/GTF file. This only uses the chromosome, "
        "as well as start and end information per line, and ignores everything else in the file. "
        "Multiple region options can be provided, see also `--filter-region-set`."
    );
    filter_region_gff_.option->check( CLI::ExistingFile );
    filter_region_gff_.option->group( group );

    // Add option for genomic region filter by PLINK MAP/BIM file.
    filter_region_bim_.option = sub->add_option(
        "--filter-region-map-bim",
        filter_region_bim_.value,
        "Genomic positions to filter for, as a MAP or BIM file as used in PLINK. This only "
        "uses the chromosome and coordinate per line, and ignores everything else in the file. "
        "Multiple region options can be provided, see also `--filter-region-set`."
    );
    filter_region_bim_.option->check( CLI::ExistingFile );
    filter_region_bim_.option->group( group );

    // Add option for genomic region filter by VCF file.
    filter_region_vcf_.option = sub->add_option(
        "--filter-region-vcf",
        filter_region_vcf_.value,
        "Genomic positions to filter for, as a VCF/BCF file (such as a known-variants file). This "
        "only uses the chromosome and position per line, and ignores everything else in the file. "
        "Multiple region options can be provided, see also `--filter-region-set`."
    );
    filter_region_vcf_.option->check( CLI::ExistingFile );
    filter_region_vcf_.option->group( group );

    // Add the set combination of the genom regions.
    filter_region_set_.option = sub->add_option(
        // If flag name is changed in the future, change it in the above options as well.
        "--filter-region-set",
        filter_region_set_.value,
        "If multiple genomic region filters are set, "
        "decide on how to combine the loci of these filters."
    );
    filter_region_set_.option->transform(
        CLI::IsMember(
            { "union", "intersection" },
            CLI::ignore_case
        )
    );
    filter_region_set_.option->group( group );
}

// =================================================================================================
//      Run Functions
// =================================================================================================

// -------------------------------------------------------------------------
//     prepare_region_filters
// -------------------------------------------------------------------------

void VariantFilterRegionOptions::prepare_region_filters() const
{
    using namespace genesis;
    using namespace genesis::population;
    using namespace genesis::utils;

    // Not running again if we already have set up a filter (i.e., if the shared pointer has data).
    if( region_filter_ ) {
        return;
    }

    // Keep track of how many filters we have added, for user output.
    size_t filter_cnt = 0;

    // Helper function to add filters according to the set union/intersection setting.
    auto add_filter_ = [&]( GenomeLocusSet&& filter )
    {
        ++filter_cnt;

        // If this is the first filter that we add, just copy it over to the shared pointer.
        if( ! region_filter_ ) {
            region_filter_ = std::make_shared<GenomeLocusSet>( std::move( filter ));
            return;
        }

        // If this is not the first filter, combine the new one with the existing one.
        if( to_lower( filter_region_set_.value ) == "union" ) {
            region_filter_->set_union( filter );
        } else if( to_lower( filter_region_set_.value ) == "intersection" ) {
            region_filter_->set_intersect( filter );
        } else {
            throw CLI::ValidationError(
                filter_region_set_.option->get_name() + "(" + filter_region_set_.value + ")",
                "Invalid value."
            );
        }
    };

    // Add the region string filters.
    for( auto const& value : filter_region_.value ) {
        GenomeLocusSet loci;
        loci.add( parse_genome_region( value ));
        add_filter_( std::move( loci ));
    }

    // Add the region list files.
    for( auto const& list_file : filter_region_list_.value ) {
        LOG_MSG2 << "Reading regions list file " << list_file;
        add_filter_( GenomeRegionReader().read_as_genome_locus_set( from_file( list_file )));
    }

    // Add the regions from bed files.
    for( auto const& file : filter_region_bed_.value ) {
        LOG_MSG2 << "Reading regions BED file " << file;
        add_filter_( BedReader().read_as_genome_locus_set( from_file( file )));
    }

    // Add the regions from gff files.
    for( auto const& file : filter_region_gff_.value ) {
        LOG_MSG2 << "Reading regions GFF2/GFF3/GTF file " << file;
        add_filter_( GffReader().read_as_genome_locus_set( from_file( file )));
    }

    // Add the regions from map/bim files.
    for( auto const& file : filter_region_bim_.value ) {
        LOG_MSG2 << "Reading regions MAP/BIM file " << file;
        add_filter_( MapBimReader().read_as_genome_locus_set( from_file( file )));
    }

    // Add the regions from vcf files.
    for( auto const& file : filter_region_vcf_.value ) {
        LOG_MSG2 << "Reading regions VCF/BCF file " << file;
        add_filter_( genome_locus_set_from_vcf_file( file ));
    }

    // User output.
    if( filter_cnt > 0 ) {
        assert( region_filter_ );

        // Get counts of positions per chromosome.
        size_t full_chr = 0;
        size_t spec_chr = 0;
        size_t spec_pos = 0;
        for( auto const& chr_name : region_filter_->chromosome_names() ) {
            auto const& bv = region_filter_->chromosome_positions( chr_name );
            if( !bv.empty() && bv.get(0) ) {
                ++full_chr;
            } else {
                ++spec_chr;
                spec_pos += bv.count();
            }
        }
        auto const chr_cnt = region_filter_->chromosome_count();
        assert( chr_cnt == full_chr + spec_chr );

        // Nice log message telling us what the filters will do.
        LOG_MSG << "Applying " << to_lower( filter_region_set_.value ) << " of " << filter_cnt
                << " region filters, which in total leads to filtering for " << chr_cnt
                << " chromosome" << ( chr_cnt != 1 ? "s" : "" ) << ", with "
                << (( full_chr > 0 ) ? std::to_string( full_chr ) + " full chromosome(s)" : "" )
                << (( full_chr > 0 && spec_pos > 0 ) ? " and " : "" )
                << (
                    ( spec_pos > 0 ) ? std::to_string( spec_chr ) +
                    " chromosome(s) containing a total of " +
                    std::to_string( spec_pos ) + " individually specified positions" : ""
                );

        for( auto const& chr_name : region_filter_->chromosome_names() ) {
            auto const& bv = region_filter_->chromosome_positions( chr_name );
            if( !bv.empty() && bv.get(0) ) {
                LOG_MSG2 << " - Chromosome \"" << chr_name << "\" fully included";
            } else {
                LOG_MSG2 << " - Chromosome \"" << chr_name << "\" with "
                         << bv.count() << " specified positions";
            }
        }
    }
}
