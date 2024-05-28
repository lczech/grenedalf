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

#include "options/window.hpp"

#include "options/global.hpp"
#include "tools/misc.hpp"

#include "genesis/population/format/bed_reader.hpp"
#include "genesis/population/format/genome_region_reader.hpp"
#include "genesis/population/format/gff_reader.hpp"
#include "genesis/population/format/map_bim_reader.hpp"
#include "genesis/population/function/functions.hpp"
#include "genesis/population/function/genome_region.hpp"
#include "genesis/population/genome_region.hpp"
#include "genesis/utils/core/logging.hpp"
#include "genesis/utils/core/std.hpp"
#include "genesis/utils/text/string.hpp"

#include <cassert>
#include <memory>
#include <stdexcept>

// =================================================================================================
//      Setup Functions
// =================================================================================================

void WindowOptions::add_window_opts_to_app(
    CLI::App* sub,
    bool include_window_view_types,
    std::string const& group
) {
    // Store what type of windows we are offering.
    include_window_view_types_ = include_window_view_types;

    // Add the set combination of the genom regions.
    window_type_.option = sub->add_option(
        "--window-type",
        window_type_.value,
        // We need a bit of string trickery to use the ternary operator here...
        std::string() +
        "Type of window to use. Depending on the type, additional options might need to be provided. "
        "\n(1) `interval`: Typical sliding window over intervals of fixed length (in bases) "
        "along the genome. "
        "\n(2) `queue`: Sliding window, but instead of using a fixed length of bases along the genome, "
        "it uses a fixed number of positions of the input data. Typically used for windowing over "
        "variant positions such as (biallelic) SNPs, and useful for example when SNPs are sparse "
        "in the genome. "
        "\n(3) `single`: Treat each position of the input data as an individual window of size 1. "
        "This is typically used to process single SNPs, and equivalent to `interval` or `queue` "
        "with a width/count of 1. "
        "\n(4) `regions`: Windows corresponding to some regions of interest, such as genes. "
        "The regions need to be provided, and can be overlapping or nested as needed. " +
        (
            include_window_view_types_
            ? std::string() +
            "\n(5) `chromosomes`: Each window covers a whole chromosome. "
            "\n(6) `genome`: The whole genome is treated as a single window."
            : std::string()
        )
    );
    window_type_.option->group( group );
    window_type_.option->required();

    // The distinction between the two window stream types boils down to just allowing the
    // WindowView streams to be specified as a type here or not. At the moment, none of them
    // (i.e., ChromosomeStream) offer any extra options anyway. If that changes, we will need to
    // also conditionally add these extra options, as done for example for the interval windows below.
    if( include_window_view_types_ ) {
        window_type_.option->transform(
            CLI::IsMember(
                { "interval", "queue", "single", "regions", "chromosomes", "genome" },
                CLI::ignore_case
            )
        );
    } else {
        window_type_.option->transform(
            CLI::IsMember(
                { "interval", "queue", "single", "regions" },
                CLI::ignore_case
            )
        );
    }

    // -------------------------------------------------------------------------
    //     Interval window
    // -------------------------------------------------------------------------

    // Width
    interval_width_.option = sub->add_option(
        "--window-interval-width",
        interval_width_.value,
        "Required when using `" + window_type_.option->get_name() + " interval`: "
        "Width of each window along the chromosome, in bases."
    );
    interval_width_.option->group( group );
    // interval_width_.option->required();

    // Stride
    interval_stride_.option = sub->add_option(
        "--window-interval-stride",
        interval_stride_.value,
        "When using `" + window_type_.option->get_name() + " interval`: "
        "Stride between windows along the chromosome, that is how far to move to get to the next "
        "window. If set to 0 (default), this is set to the same value as the `" +
        interval_width_.option->get_name() + "`."
    );
    interval_stride_.option->group( group );

    // -------------------------------------------------------------------------
    //     Queue window
    // -------------------------------------------------------------------------

    // Count
    queue_count_.option = sub->add_option(
        "--window-queue-count",
        queue_count_.value,
        "Required when using `" + window_type_.option->get_name() + " queue`: "
        "Number of positions in the genome in each window. This is most commonly used when also "
        "filtering for variant positions such as (biallelic) SNPs (which most commands do "
        "implicitly), so that each window of the analysis conists of a fixed number of SNPs, "
        "instead of a fixed length alogn the genome."
    );
    queue_count_.option->group( group );
    // queue_count_.option->required();

    // Stride
    queue_stride_.option = sub->add_option(
        "--window-queue-stride",
        queue_stride_.value,
        "When using `" + window_type_.option->get_name() + " queue`: "
        "Stride of positions between windows along the chromosome, that is how many positions "
        "does the window move forward each time. If set to 0 (default), this is set to the same "
        "value as the `" + queue_count_.option->get_name() + "`, meaning that each new window "
        "consists of new positions."
    );
    queue_stride_.option->group( group );

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
        "Multiple region options can be provided to add region windows to be processed."
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
        "Multiple region options can be provided to add region windows to be processed."
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
        "Multiple region options can be provided to add region windows to be processed."
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
        "Multiple region options can be provided to add region windows to be processed."
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
    //     "Multiple region options can be provided to add region windows to be processed."
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
    //     "Multiple region options can be provided to add region windows to be processed."
    // );
    // region_vcf_.option->check( CLI::ExistingFile );
    // region_vcf_.option->group( group );

    // Skip empty regions
    region_skip_empty_.option = sub->add_flag(
        "--window-region-skip-empty",
        region_skip_empty_.value,
        "When using `" + window_type_.option->get_name() + " regions`: "
        "In cases where there is no data in the input files for a region window, by default, "
        "we produce some \"empty\" or NaN output. "
        "With this option however, regions without data are skipped in the output."
    );
    region_skip_empty_.option->group( group );
}

// =================================================================================================
//      Run Functions
// =================================================================================================

// -------------------------------------------------------------------------
//     window_type
// -------------------------------------------------------------------------

WindowOptions::WindowType WindowOptions::window_type() const
{
    auto const win_type = genesis::utils::to_lower( window_type_.value );
    if( win_type == "interval" ) {
        return WindowType::kInterval;
    } else if( win_type == "queue" ) {
        return WindowType::kQueue;
    } else if( win_type == "single" ) {
        return WindowType::kSingle;
    } else if( win_type == "regions" ) {
        return WindowType::kRegions;
    } else if( win_type == "chromosomes" ) {
        return WindowType::kChromosomes;
    } else if( win_type == "genome" ) {
        return WindowType::kGenome;
    } else {
        throw CLI::ValidationError(
            window_type_.option->get_name(),
            "Invalid window type '" + window_type_.value + "'"
        );
    }
}

// -------------------------------------------------------------------------
//     get_variant_window_stream
// -------------------------------------------------------------------------

std::unique_ptr<VariantWindowStream> WindowOptions::get_variant_window_stream(
    VariantInputOptions const& variant_input
) const {

    // Safety check. If this is set, we've made a mistake in a command setup,
    // by adding both the Window and Window View streams, but requesting the Window stream here,
    // which is not available in this situation.
    if( include_window_view_types_ ) {
        throw std::domain_error(
            "Internal error: Cannot use Window Stream when Window View Streams are available."
        );
    }

    // Check that no extra options were provided.
    check_options_();

    // Get the input stream, and store the options for the report later.
    auto& input_stream = variant_input.get_stream();
    variant_input_ = &variant_input;

    // Get the window types that are available as Window streams.
    std::unique_ptr<VariantWindowStream> result;
    switch( window_type() ) {
        case WindowType::kInterval: {
            result = genesis::utils::make_unique<VariantIntervalWindowStream>(
                get_variant_window_stream_interval_( input_stream )
            );
            break;
        }
        case WindowType::kQueue: {
            result = genesis::utils::make_unique<VariantQueueWindowStream>(
                get_variant_window_stream_queue_( input_stream )
            );
            break;
        }
        case WindowType::kSingle: {
            result = genesis::utils::make_unique<VariantIntervalWindowStream>(
                get_variant_window_stream_single_( input_stream )
            );
            break;
        }
        case WindowType::kRegions: {
            result = genesis::utils::make_unique<VariantRegionWindowStream>(
                get_variant_window_stream_regions_( input_stream )
            );
            break;
        }
        default: {
            // This would also catch the chromosome and genome types,
            // but they should already be caught by the CLI11 validation IsMember check anyway.
            throw CLI::ValidationError(
                window_type_.option->get_name(),
                "Invalid window type '" + window_type_.value + "'."
            );
        }
    }
    assert( result );

    // Add logging to the stream, and return it.
    using VariantWindowType = genesis::population::Window<genesis::population::Variant>;
    result->add_on_enter_observer([ this ]( VariantWindowType const& window ){
        ++num_windows_;
        LOG_MSG2 << "    At window "
                 << window.chromosome() << ":"
                 << window.first_position() << "-"
                 << window.last_position();
    });
    return result;
}

// -------------------------------------------------------------------------
//     get_variant_window_view_stream
// -------------------------------------------------------------------------

std::unique_ptr<VariantWindowViewStream> WindowOptions::get_variant_window_view_stream(
    VariantInputOptions const& variant_input
) const {
    // Derived type of the wrapper class that we need for the Window Streams
    using WindowViewStream = genesis::population::WindowViewStream<
        genesis::population::VariantInputStream::Iterator
    >;

    // Safety check. If this is set, we've made a mistake in a command setup,
    // by not adding both the Window and Window View streams, but requesting them here.
    if( ! include_window_view_types_ ) {
        throw std::domain_error(
            "Internal error: Window View Streams are not available."
        );
    }

    // Check that no extra options were provided.
    check_options_();

    // Get the input stream, and store the options for the report later.
    auto& input_stream = variant_input.get_stream();
    variant_input_ = &variant_input;

    // Longer switch between all supported window types.
    // For the ones that yield Window Streams, we additionally need to wrap them,
    // so that they become Window View streams instead.
    std::unique_ptr<VariantWindowViewStream> result;
    switch( window_type() ) {
        case WindowType::kInterval: {
            result = genesis::utils::make_unique<WindowViewStream>(
                make_window_view_stream(
                    get_variant_window_stream_interval_( input_stream )
                )
            );
            break;
        }
        case WindowType::kQueue: {
            result = genesis::utils::make_unique<WindowViewStream>(
                make_window_view_stream(
                    get_variant_window_stream_queue_( input_stream )
                )
            );
            break;
        }
        case WindowType::kSingle: {
            result = genesis::utils::make_unique<WindowViewStream>(
                make_window_view_stream(
                    get_variant_window_stream_single_( input_stream )
                )
            );
            break;
        }
        case WindowType::kRegions: {
            result = genesis::utils::make_unique<WindowViewStream>(
                make_window_view_stream(
                    get_variant_window_stream_regions_( input_stream )
                )
            );
            break;
        }
        case WindowType::kChromosomes: {
            result = genesis::utils::make_unique<ChromosomeWindowStream>(
                get_variant_window_view_stream_chromosomes_( input_stream )
            );
            break;
        }
        case WindowType::kGenome: {
            result = genesis::utils::make_unique<GenomeWindowStream>(
                get_variant_window_view_stream_genome_( input_stream )
            );
            break;
        }
        default: {
            throw CLI::ValidationError(
                window_type_.option->get_name(),
                "Invalid window type '" + window_type_.value + "'."
            );
        }
    }
    assert( result );

    // Add logging to the stream, and return it.
    using VariantWindowViewType = genesis::population::WindowView<genesis::population::Variant>;
    result->add_on_enter_observer([ this ]( VariantWindowViewType const& window ){
        ++num_windows_;
        LOG_MSG2 << "    At window "
                 << window.chromosome() << ":"
                 << window.first_position() << "-"
                 << window.last_position();
    });
    return result;
}

// =================================================================================================
//      Reporting Functions
// =================================================================================================

// -------------------------------------------------------------------------
//     print_report
// -------------------------------------------------------------------------

void WindowOptions::print_report() const
{
    if( !variant_input_ ) {
        throw std::domain_error(
            "Internal error: Window report called without actually using a Window."
        );
    }

    // If phrasing here is changed, it should also be changed in VariantInputOptions::print_report()
    auto const chr_cnt = variant_input_->get_num_chromosomes();
    auto const pos_cnt = variant_input_->get_num_positions();
    LOG_MSG << "\nProcessed " << chr_cnt << " chromosome" << ( chr_cnt != 1 ? "s" : "" )
            << " with " << pos_cnt << " (non-filtered) position" << ( pos_cnt != 1 ? "s" : "" )
            << " in " << num_windows_ << " window" << ( num_windows_ != 1 ? "s" : "" ) << ".";
}

// =================================================================================================
//      Internal Functions
// =================================================================================================

// -------------------------------------------------------------------------
//     check_options_
// -------------------------------------------------------------------------

void WindowOptions::check_options_() const
{
    // Check that no interval window opts are provided unless interval window was selected.
    bool const has_interval_opt = ( *interval_width_.option || *interval_stride_.option );
    if( has_interval_opt && window_type_.value != "interval" ) {
        throw CLI::ValidationError(
            window_type_.option->get_name(),
            "Window type \"" + window_type_.value +
            "\" specified, but options for type \"interval\" were provided."
        );
    }

    // Check that no queue window opts are provided unless queue window was selected.
    bool const has_queue_opt   = ( *queue_count_.option || *queue_stride_.option );
    if( has_queue_opt && window_type_.value != "queue" ) {
        throw CLI::ValidationError(
            window_type_.option->get_name(),
            "Window type \"" + window_type_.value +
            "\" specified, but options for type \"queue\" were provided."
        );
    }

    // Check that no region window opts are provided unless region window was selected.
    // The CLI options might not all be set, as we currently leave out map/bim and vcf regions.
    // We hence here check if those are present, and rely on short circuit eval of it to not check
    // the pointers unless they are valid. That makes it easy in the future to activate them.
    bool const has_region_opt  = (
        *region_.option || *region_list_.option ||
        *region_bed_.option || *region_gff_.option ||
        ( region_bim_.option && *region_bim_.option) ||
        ( region_vcf_.option && *region_vcf_.option ) ||
        *region_skip_empty_.option
    );
    if( has_region_opt && window_type_.value != "regions" ) {
        throw CLI::ValidationError(
            window_type_.option->get_name(),
            "Window type \"" + window_type_.value +
            "\" specified, but options for type \"regions\" were provided."
        );
    }
}

// -------------------------------------------------------------------------
//     get_variant_window_stream_interval_
// -------------------------------------------------------------------------

WindowOptions::VariantIntervalWindowStream
WindowOptions::get_variant_window_stream_interval_(
    genesis::population::VariantInputStream& input
) const {
    if( ! *interval_width_.option ) {
        throw CLI::ValidationError(
            interval_width_.option->get_name(),
            "Window width has to be provided when using `" +
            window_type_.option->get_name() + " interval`."
        );
    }
    if( interval_width_.value == 0 ) {
        throw CLI::ValidationError(
            interval_width_.option->get_name(),
            "Invalid window width for interval window, has to be greater than 0."
        );
    }

    return genesis::population::make_default_interval_window_stream(
        input.begin(), input.end(), interval_width_.value, interval_stride_.value
    );
}

// -------------------------------------------------------------------------
//     get_variant_window_stream_queue_
// -------------------------------------------------------------------------

WindowOptions::VariantQueueWindowStream
WindowOptions::get_variant_window_stream_queue_(
    genesis::population::VariantInputStream& input
) const {
    if( ! *queue_count_.option ) {
        throw CLI::ValidationError(
            queue_count_.option->get_name(),
            "Position count has to be provided when using `" +
            window_type_.option->get_name() + " queue`."
        );
    }
    if( queue_count_.value == 0 ) {
        throw CLI::ValidationError(
            queue_count_.option->get_name(),
            "Invalid position count for queue window, has to be greater than 0."
        );
    }

    return genesis::population::make_passing_variant_queue_window_stream(
        input.begin(), input.end(), queue_count_.value, queue_stride_.value
    );
}

// -------------------------------------------------------------------------
//     get_variant_window_stream_single_
// -------------------------------------------------------------------------

WindowOptions::VariantIntervalWindowStream
WindowOptions::get_variant_window_stream_single_(
    genesis::population::VariantInputStream& input
) const {
    // Always return a interval interval stream with window width 1.
    return genesis::population::make_default_interval_window_stream(
        input.begin(), input.end(), 1, 1
    );

    // TODO change this to not use a interval window, but just a view based on the stream!
}

// -------------------------------------------------------------------------
//     get_variant_window_stream_regions_
// -------------------------------------------------------------------------

WindowOptions::VariantRegionWindowStream
WindowOptions::get_variant_window_stream_regions_(
    genesis::population::VariantInputStream& input
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
            "Window type \"" + window_type_.value + "\" specified, but no region inputs are given. "
            "At least one type of input containing regions to process has to be provided."
        );
    }

    // Count how many regions we have, in total, for user convenience.
    LOG_MSG << "Streaming over " << region_list->total_region_count() << " total regions across "
            << region_list->chromosome_count() << " chromosomes.";
    for( auto const& chr : region_list->chromosome_names() ) {
        LOG_MSG2 << " - Chromosome \"" << chr << "\" with "
                 << region_list->region_count( chr ) << " regions";
    }

    // Return a stream over the regions that we just read.
    // It copies the shared pointer to the regions, keeping it alive.
    auto result = genesis::population::make_default_region_window_stream(
        input.begin(), input.end(), region_list
    );
    result.skip_empty_regions( region_skip_empty_.value );

    return result;
}

// -------------------------------------------------------------------------
//     get_variant_window_view_stream_chromosomes_
// -------------------------------------------------------------------------

WindowOptions::ChromosomeWindowStream
WindowOptions::get_variant_window_view_stream_chromosomes_(
    genesis::population::VariantInputStream& input
) const {
    assert( variant_input_ );
    auto it = genesis::population::make_default_chromosome_window_stream(
        input.begin(), input.end()
    );
    it.sequence_dict( variant_input_->get_reference_dict() );
    return it;
}

// -------------------------------------------------------------------------
//     get_variant_window_view_stream_genome_
// -------------------------------------------------------------------------

WindowOptions::GenomeWindowStream
WindowOptions::get_variant_window_view_stream_genome_(
    genesis::population::VariantInputStream& input
) const {
    return genesis::population::make_default_genome_window_stream(
        input.begin(), input.end()
    );
}
