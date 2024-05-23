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

#include "options/variant_input.hpp"

#include "options/global.hpp"
#include "options/variant_file_frequency_table.hpp"
#include "options/variant_file_pileup.hpp"
#include "options/variant_file_sam.hpp"
#include "options/variant_file_sync.hpp"
#include "options/variant_file_vcf.hpp"
#include "tools/misc.hpp"

#include "genesis/population/stream/variant_parallel_input_stream.hpp"
#include "genesis/population/filter/sample_counts_filter.hpp"
#include "genesis/population/filter/variant_filter_positional.hpp"
#include "genesis/population/filter/variant_filter.hpp"
#include "genesis/population/function/functions.hpp"
#include "genesis/population/function/genome_locus_set.hpp"
#include "genesis/population/function/variant_input_stream.hpp"
#include "genesis/sequence/functions/dict.hpp"
#include "genesis/utils/core/algorithm.hpp"
#include "genesis/utils/core/fs.hpp"
#include "genesis/utils/core/info.hpp"
#include "genesis/utils/core/logging.hpp"
#include "genesis/utils/core/options.hpp"
#include "genesis/utils/core/std.hpp"
#include "genesis/utils/text/char.hpp"
#include "genesis/utils/text/convert.hpp"
#include "genesis/utils/text/string.hpp"

#include <algorithm>
#include <cassert>
#include <memory>
#include <stdexcept>

// =================================================================================================
//      Enum Mapping
// =================================================================================================

// Translate strings to enums for the VariantParallelInputStream.
// As described below, we simplify this from the more complex and more powerful approach that is
// offered by VariantParallelInputStream to just the choice of union or intersection of all
// input loci when having multiple input files, in order to keep the command line interface simple.
std::vector<
    std::pair<std::string, genesis::population::VariantParallelInputStream::ContributionType>
> const multi_file_contribution_type_map_ = {
    { "union",        genesis::population::VariantParallelInputStream::ContributionType::kCarrying },
    { "intersection", genesis::population::VariantParallelInputStream::ContributionType::kFollowing }
};

// =================================================================================================
//      CLI Setup Functions
// =================================================================================================

// -------------------------------------------------------------------------
//     All Input File Types
// -------------------------------------------------------------------------

void VariantInputOptions::add_variant_input_opts_to_app(
    CLI::App* sub,
    Suboptions const& suboptions,
    std::string const& group
) {

    // Add the basic options first.
    add_input_files_opts_to_app( sub );

    // Additional options next.
    if( suboptions.with_reference_genome_opts ) {
        reference_genome_options_.add_reference_genome_opts_to_app( sub, group );
    }
    if( suboptions.with_sample_name_opts ) {
        sample_name_options_.add_sample_name_opts_to_app( sub );
    }
    if( suboptions.with_region_filter_opts ) {
        region_filter_options_.add_region_filter_opts_to_app( sub );
    }
    if( suboptions.with_mask_filter_opts ) {
        mask_filter_options_.add_mask_filter_opts_to_app( sub );
    }
}

void VariantInputOptions::add_variant_input_opts_to_app(
    CLI::App* sub,
    std::string const& group
) {
    auto subopts = Suboptions();
    add_variant_input_opts_to_app( sub, subopts, group );
}

// -------------------------------------------------------------------------
//     add_input_files_opts_to_app
// -------------------------------------------------------------------------

void VariantInputOptions::add_input_files_opts_to_app(
    CLI::App* sub,
    std::string const& group
) {

    // Add input file type options. This is the only point where we explicitly state
    // which file types we want to add. If we want to make some of them optional later for
    // certain commands - here is the place to do so.
    input_files_.emplace_back( genesis::utils::make_unique<VariantFileSamOptions>() );
    input_files_.emplace_back( genesis::utils::make_unique<VariantFilePileupOptions>() );
    input_files_.emplace_back( genesis::utils::make_unique<VariantFileSyncOptions>() );
    input_files_.emplace_back( genesis::utils::make_unique<VariantFileVcfOptions>() );
    input_files_.emplace_back( genesis::utils::make_unique<VariantFileFrequencyTableOptions>() );

    // Now add all command line arguments of these file types to the CLI app,
    // in the order in which we added them above.
    for( auto const& input_file : input_files_ ) {
        input_file->add_file_input_opt_to_app( sub );
    }

    // Multi file set operation. In the VariantParallelInputStream, we have a different way
    // of expressing which loci of which input source to visit, but that would be way too complex
    // to specify via a command line interface. So, at least for now, we simply boil it down
    // to either the union of all loci, or their intersection.
    multi_file_loci_set_.option = sub->add_option(
        "--multi-file-locus-set",
        multi_file_loci_set_.value,
        "When multiple input files are provided, select whether the union of all their loci is "
        "used (outer join), or their intersection (inner join). For their union, input files that "
        "do not have data at a particular locus are considered as missing at that locus. "
        "Note that we allow to use multiple input files even with different file types."
    );
    multi_file_loci_set_.option->group( group );
    multi_file_loci_set_.option->transform(
        CLI::IsMember( enum_map_keys( multi_file_contribution_type_map_ ), CLI::ignore_case )
    );

    // Hidden options to set the Generic Input Stream block sizes for speed.

    // First for the main block size of the stream that is collecing all Variants,
    // which is either a single file, or the parallel input stream over multiple files.
    iterator_block_size_.option = sub->add_option(
        "--block-size",
        iterator_block_size_.value,
        "Size of the buffer block that is used for the multithreaded file parsing in the background. "
        "This option is an optimization feature that can be experiment with for larger datasets."
    );
    iterator_block_size_.option->group( "" );

    // Second for the inner streams of _all_ input files when using the parallel input
    // stream over multiple files. This will spawn a thread for each input file if set to > 0.
    parallel_block_size_.option = sub->add_option(
        "--parallel-block-size",
        parallel_block_size_.value,
        "Size of the buffer blocks that is used for the multithreaded file parsing in the background, "
        "which is used for each input file when iterating over multiple files in parallel. "
        "This will spawn a separate reading thread for each input file when set to a value > 0."
        "This option is an optimization feature that can be experiment with for larger datasets."
    );
    parallel_block_size_.option->group( "" );
}

// =================================================================================================
//      Command Setup Functions
// =================================================================================================

void VariantInputOptions::add_individual_filter_and_transforms(
    std::function<bool(genesis::population::Variant&)> const& func
) const {
    // Checks for internal correct setup
    if( static_cast<bool>( stream_ )) {
        throw std::domain_error(
            "Internal error: Calling add_individual_filter_and_transforms() "
            "after input source iteration has already been started."
        );
    }

    // Do not add a filter if the function is empty.
    // This can happen with the numerical filters - if the user did not provide any filter options,
    // we do not need to add it to the input stream here.
    if( !func ) {
        return;
    }

    individual_filters_and_transforms_.push_back( func );
}

void VariantInputOptions::add_combined_filter_and_transforms(
    std::function<bool(genesis::population::Variant&)> const& func
) const {
    // Checks for internal correct setup
    if( static_cast<bool>( stream_ )) {
        throw std::domain_error(
            "Internal error: Calling add_combined_filter_and_transforms() "
            "after input source iteration has already been started."
        );
    }

    // Do not add a filter if the function is empty.
    // This can happen with the numerical filters - if the user did not provide any filter options,
    // we do not need to add it to the input stream here.
    if( !func ) {
        return;
    }

    combined_filters_and_transforms_.push_back( func );
}

// =================================================================================================
//      Reporting Functions
// =================================================================================================

size_t VariantInputOptions::get_input_file_count() const
{
    size_t file_count = 0;
    for( auto const& input_file : input_files_ ) {
        file_count += input_file->get_file_input_options().file_count();
    }
    return file_count;
}

void VariantInputOptions::print_report() const
{
    // If phrasing here is changed, it should also be changed in WindowOptions::print_report()
    auto const chr_cnt = num_chromosomes_;
    auto const pos_cnt = num_positions_;
    LOG_MSG << "\nProcessed " << chr_cnt << " chromosome" << ( chr_cnt != 1 ? "s" : "" )
            << " with " << pos_cnt << " (non-filtered) position" << ( pos_cnt != 1 ? "s" : "" ) << ".";
    if( num_mismatch_bases_ > 0 ) {
        LOG_WARN << "Out of these, " << num_mismatch_bases_ << " positions had a mismatch between "
                 << "the bases as determined from the input file(s) compared to the provided "
                 << "reference genome. Please check that you provided a fitting reference genome!";
    }
}

// =================================================================================================
//      Stream Setup
// =================================================================================================

// -------------------------------------------------------------------------
//     prepare_inputs_
// -------------------------------------------------------------------------

void VariantInputOptions::prepare_inputs_() const
{
    // Checks for internal correct setup
    if( static_cast<bool>( stream_ )) {
        // Nothing to be done. We already prepared the data.
        return;
    }

    // First prepare the reference genome and dict, as this might be needed by the steam.
    // Very hacky... at the moment, only the Frequency Table Input actually makes use of this,
    // but we somehow need to get the ref genome to it before we initialize its stream...
    // So we do this here, and in order to avoid ugly casting of our different input file types,
    // we just "set" the ref genome for all, which is a dummy function for all other types.
    for( auto const& input_file : input_files_ ) {
        input_file->add_reference_genome( get_reference_genome() );
    }

    // We also do a quick check if the different types of references are compatible with each other.
    // We can only check those whose formats cover the full length of the genome that we expect.
    // At the moment, for instance, we do not check that the region filters cover every position,
    // as that would only make sense for the fasta mask type region filter, which is kind of an
    // outlyer anyway.
    // We do test the reference genome and mask file here though, if both are given.
    mask_filter_options_.check_mask_against_reference( get_reference_dict() );
}

// -------------------------------------------------------------------------
//     prepare_stream_
// -------------------------------------------------------------------------

void VariantInputOptions::prepare_stream_() const
{
    using namespace genesis;
    using namespace genesis::population;
    using namespace genesis::utils;

    // Checks for internal correct setup
    if( static_cast<bool>( stream_ )) {
        // Nothing to be done. We already prepared the data.
        return;
    }

    // Get and check the number of input files provided.
    size_t const file_count = get_input_file_count();
    size_t const max_file_count = genesis::utils::info_max_file_count();
    size_t const file_count_margin = 10;
    if( max_file_count > 0 && file_count > max_file_count - file_count_margin ) {
        LOG_WARN << "In total, " <<  file_count << " input files are provided. However, the system "
                 << "limit for the number of concurrently opened files is " << max_file_count
                 << ( file_count <= max_file_count
                      ? ", leaving not much margin for other needed files such as outputs. "
                      : ". "
                  ) << "The command might hence fail. If that happens, process your files in "
                  << "batches, or merge them into fewer files first.";
    }

    // Check how many files are given. If it is a single one, we just use that for the input
    // stream, which is faster than piping it through a parallel stream. For multiple files,
    // we create a parallel stream.
    // Set up stream depending on how many files we found in total across all inputs.
    if( file_count == 0 ) {
        throw CLI::ValidationError(
            "Input sources", "At least one input file has to be provided."
        );
    } else if( file_count == 1 ) {
        prepare_stream_single_file_();
    } else {
        prepare_stream_multiple_files_();
    }
    assert( stream_ );

    // Some user output.
    LOG_MSG1 << "Processing " << sample_names().size()
             << " sample" << ( sample_names().size() == 1 ? "" : "s" );
    for( auto const& sn : sample_names() ) {
        LOG_MSG2 << " - " << sn;
    }
    LOG_MSG1;
    if( contains_duplicates( sample_names() )) {
        LOG_WARN << "The input contains duplicate sample names. Some grenedalf commands can work "
                 << "with that internally, other not (and will then fail). Either way, this will "
                 << "make working with the output more difficult downstream, and we recommend "
                 << "to clean up the names first. "
                 << "Use verbose output (`--verbose`) to show a list of all sample names.";
        LOG_MSG1;
    }
}

// -------------------------------------------------------------------------
//     prepare_stream_
// -------------------------------------------------------------------------

void VariantInputOptions::prepare_stream_( size_t first, size_t last ) const
{
    using namespace genesis;
    using namespace genesis::population;
    using namespace genesis::utils;

    // Boundary checks, just in case. Should be done correctly internally already.
    size_t const file_count = get_input_file_count();
    if( first >= file_count || last > file_count || first >= last ) {
        throw std::logic_error( "Internal error: Requesting invalid range of file streams." );
    }
    assert( last - first > 0 );

    // Check how many files are given. Other than in the base class, we here always create a parallel
    // stream. We do this also in the edge case that we only want a range of files which happens
    // to be of size one (last-first==1), as that should be rare, and avoids duplicating special
    // range code for the single file function.
    if( file_count == 0 ) {
        throw CLI::ValidationError(
            "Input sources", "At least one input file has to be provided."
        );
    } else {
        prepare_stream_multiple_files_( first, last );
    }
    // assert( stream_ );

    // Some user output.
    LOG_MSG1 << "Processing batch of " << ( last - first ) << " input file"
             << ( last - first == 1 ? "" : "s" ) << " [ " << ( first + 1 ) << ", " << last << " ]"
             << " with " << sample_names().size()
             << " sample" << ( sample_names().size() == 1 ? "" : "s" );
    for( auto const& sn : sample_names() ) {
        LOG_MSG2 << " - " << sn;
    }
    LOG_MSG1;
    if( contains_duplicates( sample_names() )) {
        LOG_WARN << "The input contains duplicate sample names. Some grenedalf commands can work "
                 << "with that internally, other not (and will then fail). Either way, this will "
                 << "make working with the output more difficult downstream, and we recommend "
                 << "to clean up the names first. "
                 << "Use verbose output (`--verbose`) to show a list of all sample names.";
        LOG_MSG1;
    }
}

// -------------------------------------------------------------------------
//     prepare_stream_single_file_
// -------------------------------------------------------------------------

void VariantInputOptions::prepare_stream_single_file_() const
{
    // Assert that this function is only called in a context where the data is not yet prepared.
    internal_check(
        ! static_cast<bool>( stream_ ),
        "prepare_stream_single_file_() called in an invalid context."
    );

    // Get the input file as a pointer, using the one input that the user provided,
    // and asserting that it only was one that was provided.
    VariantFileOptions const* provided_input_file = nullptr;
    for( auto const& input_file : input_files_ ) {
        if( input_file->get_file_input_options().provided() ) {
            if( provided_input_file == nullptr ) {
                provided_input_file = input_file.get();
            } else {
                internal_check(
                    false, "prepare_stream_single_file_() called with more than one file provided"
                );
            }
        }
    }
    internal_check(
        provided_input_file, "prepare_stream_single_file_() called with no file provided"
    );

    // Prepare the stream depending on the input file format, using the pointer.
    assert( provided_input_file );
    assert( provided_input_file->get_file_input_options().file_count() == 1 );
    stream_ = provided_input_file->get_stream(
        provided_input_file->get_file_input_options().file_paths()[0]
    );

    // Apply sample renaming and filtering. We need to do this here first, before processing
    // any other filters and transforms, to make sure that they all see the correct
    // sample names and numbers.
    sample_name_options_.rename_samples( stream_.data().sample_names );
    sample_name_options_.add_sample_name_filter( stream_ );
    sample_name_options_.apply_sample_group_merging( stream_ );
    conditionally_make_gapless_stream_();

    // Copy over the sample names from the stream, so that they are accessible.
    if( sample_names().empty() ) {
        throw std::runtime_error( "Invalid input that does not contain any samples." );
    }
    internal_check(
        static_cast<bool>( stream_ ),
        "prepare_stream_single_file_() call to file prepare function did not succeed."
    );

    // Add an observer that checks chromosome order and length. We need to check order here,
    // as the single stream does not do this already (as opposed to the parallel one below).
    // Without a dict (when sequence_dict_ is nullptr), this checks lexicographically.
    stream_.add_observer(
        genesis::population::make_variant_input_stream_sequence_order_observer(
            get_reference_dict()
        )
    );

    // Add the region filters.
    // need to refactor, rename, and add a filter for each sample.
    // we want:
    //  - filter for each input source (to skip regions etc), taking a variant
    //  - filter for combined variant of the final stream, eg for ref base guessing
    //  - filter for each BaseCounts, eg for min coverage per sample
    add_individual_filters_and_transforms_to_stream_( stream_ );
    add_combined_filters_and_transforms_to_stream_( stream_ );

    // Set the thread pool of the Generic Input Stream to be used. We want to use the global thread pool
    // here, as otherwise, each input file would spawn it's own thread, which could easily go into
    // the hundreds or more, and hence lead to slowdown due to blocking issues, or even problems
    // on clusters that monitor CPU usage. Well, not for a single file, but still, better to
    // use the global pool. This is more relevant for multiple files, see below.
    stream_.thread_pool( global_options.thread_pool() );

    // Set the buffer block size of the stream, for multi-threaded processing speed.
    // Needs to be tested - might give more or less advantage depending on setting.
    stream_.block_size( iterator_block_size_.value );

    // TODO also add filters depending on file type. sam pileup and sync might need
    // biallelic snp filters (using the merged variants filter type setting), while vcf might not?!
    // or has it already set below or in the variant input stream - need to check.
}

// -------------------------------------------------------------------------
//     prepare_stream_multiple_files_
// -------------------------------------------------------------------------

void VariantInputOptions::prepare_stream_multiple_files_() const
{
    using namespace genesis;
    using namespace genesis::population;
    using namespace genesis::utils;

    // Assert that this function is only called in a context where the data is not yet prepared.
    internal_check(
        ! static_cast<bool>( stream_ ),
        "prepare_stream_multiple_files_() called in an invalid context."
    );

    // Get whether the user wants the union or intersection of all parallel input loci.
    auto const contribution = get_enum_map_value(
        multi_file_contribution_type_map_, multi_file_loci_set_.value
    );

    // Make a parallel input stream, and add all input from all file formats to it,
    // using the same contribution type for all of them, which either results in the union
    // or the intersection of all input loci. See VariantParallelInputStream for details.
    size_t file_count = 0;
    VariantParallelInputStream parallel_stream;
    for( auto const& input_file : input_files_ ) {
        for( auto const& file : input_file->get_file_input_options().file_paths() ) {
            parallel_stream.add_variant_input_stream(
                input_file->get_stream( file ), contribution
            );
            ++file_count;
        }
    }

    // Check that we have multiple input files.
    internal_check(
        file_count > 1, "prepare_stream_multiple_files_() called with just one file provided"
    );

    // We here split the second step of turning the parallel iteratr into an actual stream
    // into a separate function, solely so that the VariantInputOptionsRange class that derives
    // from this for processing ranges of input files can re-use that function.
    prepare_stream_from_parallel_stream_( std::move( parallel_stream ));
}

// -------------------------------------------------------------------------
//     prepare_stream_multiple_files_
// -------------------------------------------------------------------------

void VariantInputOptions::prepare_stream_multiple_files_( size_t first, size_t last ) const
{
    using namespace genesis;
    using namespace genesis::population;
    using namespace genesis::utils;

    // Assert correct usage.
    size_t const file_count = get_input_file_count();
    internal_check(
        first < last && first < file_count && last <= file_count,
        "Invalid subset for multiple file range."
    );

    // Get whether the user wants the union or intersection of all parallel input loci.
    auto const contribution = get_enum_map_value(
        multi_file_contribution_type_map_, multi_file_loci_set_.value
    );

    // Make a parallel input stream, and add all input from all file formats to it,
    // using the same contribution type for all of them, which either results in the union
    // or the intersection of all input loci. See VariantParallelInputStream for details.
    // As we only want to use a range of the input files, subset accordingly.
    // We cannot simply select the range of files directly, as they are spread across diffrent
    // input file types, which each could have a varying number of files provided. It's easier
    // to just loop through all of them, than to try to select sub-ranges on that.
    size_t input_count = 0;
    VariantParallelInputStream parallel_stream;
    for( auto const& input_file : input_files_ ) {
        for( auto const& file : input_file->get_file_input_options().file_paths() ) {
            if( input_count >= first && input_count < last ) {
                parallel_stream.add_variant_input_stream(
                    input_file->get_stream( file ), contribution
                );
            }
            ++input_count;
        }
    }

    // Check that we have the correct range.
    internal_check( input_count == file_count, "Invalid input file count." );
    internal_check( parallel_stream.input_size() == last - first, "Invalid input file subset." );

    // Now turn the parallel stream into the type that we want, with all additional settings
    // applied as well. This function is provided by the base class, so that we do not have to
    // duplicate it here.
    prepare_stream_from_parallel_stream_( std::move( parallel_stream ));
}

// -------------------------------------------------------------------------
//     prepare_stream_from_parallel_stream_
// -------------------------------------------------------------------------

void VariantInputOptions::prepare_stream_from_parallel_stream_(
     genesis::population::VariantParallelInputStream&& parallel_stream
) const {
    using namespace genesis::population;

    // If a sequence dict in some format was provided, we use it for the stream,
    // so that the chromosome order can be used correctly. If not provided, nullptr is also okay.
    // We check chromosome length below as well.
    parallel_stream.sequence_dict( get_reference_dict() );

    // Go through all sources again, build the sample names list from them,
    // and add the individual samples filters to all of them, e.g., so that regions are filtered
    // before the data even reaches the parallel stream. Faster!
    for( auto& input : parallel_stream.inputs() ) {
        add_individual_filters_and_transforms_to_stream_( input );

        // Following the reasoning from above: we need to use the global thread pool,
        // as otherwise reading hundreds of files will spawn too many threads, which is not good.
        input.thread_pool( global_options.thread_pool() );

        // Also set the buffer block size of the stream, using the second buffer size option
        // that is used for the inner streams of parallel input.
        // Needs to be tested - might give more or less advantage depending on setting.
        input.block_size( parallel_block_size_.value );
    }

    // Finally, create the parallel stream.
    stream_ = make_variant_input_stream_from_variant_parallel_input_stream( parallel_stream );

    // Apply sample renaming and filtering. We need to do this here first, before processing
    // any other filters and transforms, to make sure that they all see the correct
    // sample names and numbers.
    sample_name_options_.rename_samples( stream_.data().sample_names );
    sample_name_options_.add_sample_name_filter( stream_ );
    sample_name_options_.apply_sample_group_merging( stream_ );
    conditionally_make_gapless_stream_();

    // Add an observer that checks chromosome length. We do not need to check order,
    // as this is done with the above sequence dict already internally.
    if( get_reference_dict() ) {
        stream_.add_observer(
            make_variant_input_stream_sequence_length_observer( get_reference_dict() )
        );
    }

    // Add the filters and transformations that are to be applied to all samples combined.
    add_combined_filters_and_transforms_to_stream_( stream_ );

    // We also need to set the thread pool of the parallel input stream itself,
    // for the same reasons as described above.
    stream_.thread_pool( global_options.thread_pool() );

    // Set the buffer block size of the stream, for multi-threaded speed.
    // We only buffer the final parallel stream, and not its individual sources,
    // in order to not keep too many threads from spawning... Might need testing and refinement.
    stream_.block_size( iterator_block_size_.value );
}

// -------------------------------------------------------------------------
//     conditionally_make_gapless_stream_
// -------------------------------------------------------------------------

void VariantInputOptions::conditionally_make_gapless_stream_() const
{
    using namespace genesis::population;
    using namespace genesis::sequence;

    // The function checks all cases in which we want to make the stream gapless.
    // This is the case when either the command set the gapless_stream_ flag,
    // or if we have a mask file provided. In the latter case, we also set the flag first,
    // for completeness.

    // If we have a mask file provided, we make this a gapless stream either way.
    if( mask_filter_options_.get_mask() ) {
        gapless_stream_ = true;
        auto const mask_dict = std::make_shared<SequenceDict>(
            reference_locus_set_to_dict( *mask_filter_options_.get_mask() )
        );
        stream_ = make_variant_gapless_input_stream( stream_, mask_dict );
        return;
    }

    // Now we check if we want to make a gapless stream, e.g., by demand from a command.
    // If so, we further check if any references are given, which we then want to use
    // for the stream in order to properly get the chromosome lengths.
    // If none are given, the best we can do is to just fill in the gaps in the data, but we will
    // stop once the end of the data is reached, even if that is not the end of the chromosome.
    if( ! gapless_stream_ ) {
        return;
    }
    if( get_reference_genome() ) {
        stream_ = make_variant_gapless_input_stream( stream_, get_reference_genome() );
        return;
    }
    if( get_reference_dict() ) {
        stream_ = make_variant_gapless_input_stream( stream_, get_reference_dict() );
        return;
    }
    stream_ = make_variant_gapless_input_stream( stream_ );
}

// =================================================================================================
//     Filters and Transformations
// =================================================================================================

// -------------------------------------------------------------------------
//     add_individual_filters_and_transforms_to_stream_
// -------------------------------------------------------------------------

void VariantInputOptions::add_individual_filters_and_transforms_to_stream_(
    genesis::population::VariantInputStream& stream
) const {
    using namespace genesis::population;

    // These filters and transformations are applied to each input source individually.
    // For example, filtering out regions there means that everything downstream does not have
    // to deal with positions that we are about to discard anyway.

    // Add the genome region list as a filter. We do this first, to remove everything that we do
    // not want to keep anyway from having to be processed downstream as early as we can.
    // This is also the reason to apply this on each input stream individually: This way, the
    // filtering can be done in the generic input stream buffering already, and hence use
    // threading and buffering to our advantage here.
    if( region_filter_options_.get_region_filter() ) {
        stream.add_filter( make_variant_filter_by_region_excluding(
            region_filter_options_.get_region_filter()
        ));
    }

    // Add the addtional filters and transformations that might have been set by the commands.
    for( auto const& func : individual_filters_and_transforms_ ) {
        stream.add_transform_filter( func );
    }
}

// -------------------------------------------------------------------------
//     add_combined_filters_and_transforms_to_stream_
// -------------------------------------------------------------------------

void VariantInputOptions::add_combined_filters_and_transforms_to_stream_(
    genesis::population::VariantInputStream& stream
) const {
    using namespace genesis::population;

    // Filters and transformations here are applied to the combiend Variant sample,
    // which is important in the case of parallel input sources: For example, we need to take
    // all of them into account at the same time to get a proper ref and alt base guess,
    // otherwise we might end up with contradicting bases.

    // First we add the processing of the mask file, if one is given. This needs to happen on
    // the combined stream, after all region filtering, but before any numeric filters.
    if( mask_filter_options_.get_mask() ) {
        internal_check( gapless_stream_, "Stream not set to be gapless" );
        stream.add_transform( mask_filter_options_.make_mask_transform() );
    }

    // Add the addtional filters and transformations that might have been set by the commands.
    for( auto const& func : combined_filters_and_transforms_ ) {
        stream.add_transform_filter( func );
    }

    // In addition to the transforms and filters, we here also add an observer function.
    // We always print out where the input is at, at the moment. That makes sure that we always
    // get some progress update, which is probably more useful for the user than waiting too long
    // without any.
    // Note though that region filters are applied per input file, and so might have already removed
    // chromosomes, so that they won't be printed here. We use lambda capture by value to create
    // a copy of current_chr that is kept in the lambda, and updated there. We are not in C++14 yet.
    std::string current_chr;
    stream.add_observer(
        [ current_chr, this ]( Variant const& variant ) mutable {
            ++num_positions_;
            if( current_chr != variant.chromosome ) {
                LOG_MSG << "At chromosome " << variant.chromosome;
                current_chr = variant.chromosome;
                ++num_chromosomes_;
            }
        }
    );

    // If we have a reference genome, we also check that its bases match what we find in the
    // input files, and report if that's not fitting.
    if( get_reference_genome() ) {
        stream.add_observer(
            [ this ]( Variant const& variant ) mutable {
                auto const var_base = genesis::utils::to_upper( variant.reference_base );
                auto const ref_base = get_reference_genome()->get_base(
                    variant.chromosome, variant.position
                );
                if( var_base != 'N' && ref_base != 'N' && var_base != ref_base ) {
                    ++num_mismatch_bases_;
                }
            }
        );
    }
}
