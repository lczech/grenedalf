/*
    grenedalf - Genome Analyses of Differential Allele Frequencies
    Copyright (C) 2020-2023 Lucas Czech

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

#include "options/variant_input.hpp"

#include "options/global.hpp"
#include "options/variant_file_frequency_table.hpp"
#include "options/variant_file_pileup.hpp"
#include "options/variant_file_sam.hpp"
#include "options/variant_file_sync.hpp"
#include "options/variant_file_vcf.hpp"
#include "tools/misc.hpp"

#include "genesis/population/formats/variant_parallel_input_iterator.hpp"
#include "genesis/population/functions/filter_transform.hpp"
#include "genesis/population/functions/functions.hpp"
#include "genesis/population/functions/variant_input_iterator.hpp"
#include "genesis/utils/core/algorithm.hpp"
#include "genesis/utils/core/fs.hpp"
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

// Translate strings to enums for the VariantParallelInputIterator.
// As described below, we simplify this from the more complex and more powerful approach that is
// offered by VariantParallelInputIterator to just the choice of union or intersection of all
// input loci when having multiple input files, in order to keep the command line interface simple.
std::vector<
    std::pair<std::string, genesis::population::VariantParallelInputIterator::ContributionType>
> const multi_file_contribution_type_map_ = {
    { "union",        genesis::population::VariantParallelInputIterator::ContributionType::kCarrying },
    { "intersection", genesis::population::VariantParallelInputIterator::ContributionType::kFollowing }
};

// =================================================================================================
//      CLI Setup Functions
// =================================================================================================

// -------------------------------------------------------------------------
//     All Input File Types
// -------------------------------------------------------------------------

void VariantInputOptions::add_variant_input_opts_to_app(
    CLI::App* sub,
    bool with_reference_genome_opts,
    bool with_sample_name_opts,
    bool with_region_filter_opts,
    std::string const& group
) {

    // Add the basic options first.
    add_input_files_opts_to_app( sub );

    // Additional options next.
    if( with_reference_genome_opts ) {
        reference_genome_options_.add_reference_genome_opts_to_app( sub, group );
    }
    if( with_sample_name_opts ) {
        // sample_name_options_.add_sample_name_opts_to_app( sub, group );
        sample_name_options_.add_sample_name_opts_to_app( sub );
    }
    if( with_region_filter_opts ) {
        // region_filter_options_.add_region_filter_opts_to_app( sub, group );
        region_filter_options_.add_region_filter_opts_to_app( sub );
    }
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

    // Multi file set operation. In the VariantParallelInputIterator, we have a different way
    // of expressing which loci of which input source to visit, but that would be way too complex
    // to specify via a command line interface. So, at least for now, we simply boil it down
    // to either the union of all loci, or their intersection.
    multi_file_loci_set_.option = sub->add_option(
        "--multi-file-locus-set",
        multi_file_loci_set_.value,
        "When multiple input files are provided, select whether the union of all their loci is "
        "used, or their intersection. For their union, input files that do not have data at a "
        "particular locus are considered to have zero coverage at that locus. "
        "Note that we allow to use multiple input files even with different file types."
    );
    multi_file_loci_set_.option->group( group );
    multi_file_loci_set_.option->transform(
        CLI::IsMember( enum_map_keys( multi_file_contribution_type_map_ ), CLI::ignore_case )
    );

    // Hidden options to set the LambdaIterator block sizes for speed.

    // First for the main block size of the iterator that is collecing all Variants,
    // which is either a single file, or the parallel input iterator over multiple files.
    iterator_block_size_.option = sub->add_option(
        "--block-size",
        iterator_block_size_.value,
        "Size of the buffer block that is used for the multithreaded file parsing in the background. "
        "This option is an optimization feature that can be experiment with for larger datasets."
    );
    iterator_block_size_.option->group( "" );

    // Second for the inner iterators of _all_ input files when using the parallel input
    // iterator over multiple files. This will spawn a thread for each input file if set to > 0.
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
    if( static_cast<bool>( iterator_ )) {
        throw std::domain_error(
            "Internal error: Calling add_individual_filter_and_transforms() "
            "after input source iteration has already been started."
        );
    }

    // Do not add a filter if the function is empty.
    // This can happen with the numerical filters - if the user did not provide any filter options,
    // we do not need to add it to the input iterator here.
    if( !func ) {
        return;
    }

    individual_filters_and_transforms_.push_back( func );
}

void VariantInputOptions::add_combined_filter_and_transforms(
    std::function<bool(genesis::population::Variant&)> const& func
) const {
    // Checks for internal correct setup
    if( static_cast<bool>( iterator_ )) {
        throw std::domain_error(
            "Internal error: Calling add_combined_filter_and_transforms() "
            "after input source iteration has already been started."
        );
    }

    // Do not add a filter if the function is empty.
    // This can happen with the numerical filters - if the user did not provide any filter options,
    // we do not need to add it to the input iterator here.
    if( !func ) {
        return;
    }

    combined_filters_and_transforms_.push_back( func );
}

// =================================================================================================
//      Reporting Functions
// =================================================================================================

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
//      Iterator Setup
// =================================================================================================

// -------------------------------------------------------------------------
//     prepare_
// -------------------------------------------------------------------------

void VariantInputOptions::prepare_() const
{
    // Checks for internal correct setup
    if( static_cast<bool>( iterator_ )) {
        // Nothing to be done. We already prepared the data.
        return;
    }

    // First prepare the reference genome and dict, as this might be needed by the iterator.
    // Very hacky... at the moment, only the Frequency Table Input actually makes use of this,
    // but we somehow need to get the ref genome to it before we initialize its iterator...
    // So we do this here, and in order to avoid ugly casting of our different input file types,
    // we just "set" the ref genome for all, which is a dummy function for all other types.
    for( auto const& input_file : input_files_ ) {
        input_file->add_reference_genome( reference_genome_options_.get_reference_genome() );
    }

    // Then, prepare the iterator, either for a single, or for multiple input files.
    prepare_iterator_();
}

// -------------------------------------------------------------------------
//     prepare_iterator_
// -------------------------------------------------------------------------

void VariantInputOptions::prepare_iterator_() const
{
    using namespace genesis;
    using namespace genesis::population;
    using namespace genesis::utils;

    // Check how many files are given. If it is a single one, we just use that for the input
    // iterator, which is faster than piping it through a parallel iterator. For multiple files,
    // we create a parallel iterator.
    size_t file_count = 0;
    for( auto const& input_file : input_files_ ) {
        file_count += input_file->get_file_input_options().file_count();
    }

    // Set up iterator depending on how many files we found in total across all inputs.
    if( file_count == 0 ) {
        throw CLI::ValidationError(
            "Input sources", "At least one input file has to be provided."
        );
    } else if( file_count == 1 ) {
        prepare_iterator_single_file_();
    } else {
        prepare_iterator_multiple_files_();
    }
    assert( iterator_ );

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
//     prepare_iterator_single_file_
// -------------------------------------------------------------------------

void VariantInputOptions::prepare_iterator_single_file_() const
{
    // Assert that this function is only called in a context where the data is not yet prepared.
    internal_check(
        ! static_cast<bool>( iterator_ ),
        "prepare_iterator_single_file_() called in an invalid context."
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
                    false, "prepare_iterator_single_file_() called with more than one file provided"
                );
            }
        }
    }
    internal_check(
        provided_input_file, "prepare_iterator_single_file_() called with no file provided"
    );

    // Prepare the iterator depending on the input file format, using the pointer.
    assert( provided_input_file );
    assert( provided_input_file->get_file_input_options().file_count() == 1 );
    iterator_ = provided_input_file->get_iterator(
        provided_input_file->get_file_input_options().file_paths()[0]
    );

    // Apply sample renaming and filtering. We need to do this here first, before processing
    // any other filters and transforms, to make sure that they all see the correct
    // sample names and numbers.
    sample_name_options_.rename_samples( iterator_.data().sample_names );
    sample_name_options_.add_sample_name_filter( iterator_ );

    // Copy over the sample names from the iterator, so that they are accessible.
    if( sample_names().empty() ) {
        throw std::runtime_error( "Invalid input file that does not contain any samples." );
    }
    internal_check(
        static_cast<bool>( iterator_ ),
        "prepare_iterator_single_file_() call to file prepare function did not succeed."
    );

    // Add a visitor that checks chromosome order and length. We need to check order here,
    // as the single iterator does not do this already (as opposed to the parallel one below).
    // Without a dict (when sequence_dict_ is nullptr), this checks lexicographically.
    iterator_.add_visitor(
        genesis::population::make_variant_input_iterator_sequence_order_visitor(
            get_reference_dict()
        )
    );

    // Add the region filters.
    // need to refactor, rename, and add a filter for each sample.
    // we want:
    //  - filter for each input source (to skip regions etc), taking a variant
    //  - filter for combined variant of the final iterator, eg for ref base guessing
    //  - filter for each BaseCounts, eg for min coverage per sample
    add_individual_filters_and_transforms_to_iterator_( iterator_ );
    add_combined_filters_and_transforms_to_iterator_( iterator_ );

    // Set the thread pool of the LambdaIterator to be used. We want to use the global thread pool
    // here, as otherwise, each input file would spawn it's own thread, which could easily go into
    // the hundreds or more, and hence lead to slowdown due to blocking issues, or even problems
    // on clusters that monitor CPU usage. Well, not for a single file, but still, better to
    // use the global pool. This is more relevant for multiple files, see below.
    iterator_.thread_pool( global_options.thread_pool() );

    // Set the buffer block size of the iterator, for multi-threaded processing speed.
    // Needs to be tested - might give more or less advantage depending on setting.
    iterator_.block_size( iterator_block_size_.value );

    // TODO also add filters depending on file type. sam pileup and sync might need
    // biallelic snp filters (using the merged variants filter type setting), while vcf might not?!
    // or has it already set below or in the variant input iterator - need to check.
}

// -------------------------------------------------------------------------
//     prepare_iterator_multiple_files_
// -------------------------------------------------------------------------

void VariantInputOptions::prepare_iterator_multiple_files_() const
{
    using namespace genesis;
    using namespace genesis::population;
    using namespace genesis::utils;

    // Assert that this function is only called in a context where the data is not yet prepared.
    internal_check(
        ! static_cast<bool>( iterator_ ),
        "prepare_iterator_multiple_files_() called in an invalid context."
    );

    // Get whether the user wants the union or intersection of all parallel input loci.
    auto const contribution = get_enum_map_value(
        multi_file_contribution_type_map_, multi_file_loci_set_.value
    );

    // Make a parallel input iterator, and add all input from all file formats to it,
    // using the same contribution type for all of them, which either results in the union
    // or the intersection of all input loci. See VariantParallelInputIterator for details.
    size_t file_count = 0;
    VariantParallelInputIterator parallel_it;
    for( auto const& input_file : input_files_ ) {
        for( auto const& file : input_file->get_file_input_options().file_paths() ) {
            parallel_it.add_variant_input_iterator(
                input_file->get_iterator( file ), contribution
            );
            ++file_count;
        }
    }

    // Check that we have multiple input files.
    internal_check(
        file_count > 1, "prepare_iterator_multiple_files_() called with just one file provided"
    );

    // If a sequence dict in some format was provided, we use it for the iterator,
    // so that the chromosome order can be used correctly. If not provided, nullptr is also okay.
    // We check chromosome length below as well.
    parallel_it.sequence_dict( get_reference_dict() );

    // Go through all sources again, build the sample names list from them,
    // and add the individual samples filters to all of them, e.g., so that regions are filtered
    // before the data even reaches the parallel iterator. Faster!
    for( auto& input : parallel_it.inputs() ) {
        add_individual_filters_and_transforms_to_iterator_( input );

        // Following the reasoning from above: we need to use the global thread pool,
        // as otherwise reading hundreds of files will spawn too many threads, which is not good.
        input.thread_pool( global_options.thread_pool() );

        // Also set the buffer block size of the iterator, using the second buffer size option
        // that is used for the inner iterators of parallel input.
        // Needs to be tested - might give more or less advantage depending on setting.
        input.block_size( parallel_block_size_.value );
    }

    // Finally, create the parallel iterator.
    iterator_ = make_variant_input_iterator_from_variant_parallel_input_iterator( parallel_it );

    // Apply sample renaming and filtering. We need to do this here first, before processing
    // any other filters and transforms, to make sure that they all see the correct
    // sample names and numbers.
    sample_name_options_.rename_samples( iterator_.data().sample_names );
    sample_name_options_.add_sample_name_filter( iterator_ );

    // Add a visitor that checks chromosome length. We do not need to check order,
    // as this is done with the above sequence dict already internally.
    if( get_reference_dict() ) {
        iterator_.add_visitor(
            make_variant_input_iterator_sequence_length_visitor( get_reference_dict() )
        );
    }

    // Add the filters and transformations that are to be applied to all samples combined.
    add_combined_filters_and_transforms_to_iterator_( iterator_ );

    // We also need to set the thread pool of the parallel input iterator itself,
    // for the same reasons as described above.
    iterator_.thread_pool( global_options.thread_pool() );

    // Set the buffer block size of the iterator, for multi-threaded speed.
    // We only buffer the final parallel iterator, and not its individual sources,
    // in order to not keep too many threads from spawning... Might need testing and refinement.
    iterator_.block_size( iterator_block_size_.value );
}

// =================================================================================================
//     Filters and Transformations
// =================================================================================================

// -------------------------------------------------------------------------
//     add_individual_filters_and_transforms_to_iterator_
// -------------------------------------------------------------------------

void VariantInputOptions::add_individual_filters_and_transforms_to_iterator_(
    genesis::population::VariantInputIterator& iterator
) const {
    using namespace genesis::population;

    // These filters and transformations are applied to each input source individually.
    // For example, filtering out regions there means that everything downstream does not have
    // to deal with positions that we are about to discard anyway.

    // Add the genome region list as a filter.
    if( region_filter_options_.get_region_filter() ) {
        iterator.add_filter( make_filter_by_region( region_filter_options_.get_region_filter() ));
    }

    // Add the addtional filters and transformations that might have been set by the commands.
    for( auto const& func : individual_filters_and_transforms_ ) {
        iterator.add_transform_filter( func );
    }
}

// -------------------------------------------------------------------------
//     add_combined_filters_and_transforms_to_iterator_
// -------------------------------------------------------------------------

void VariantInputOptions::add_combined_filters_and_transforms_to_iterator_(
    genesis::population::VariantInputIterator& iterator
) const {
    using namespace genesis::population;

    // Filters and transformations here are applied to the combiend Variant sample,
    // which is important in the case of parallel input sources: For example, we need to take
    // all of them into account at the same time to get a proper ref and alt base guess,
    // otherwise we might end up with contradicting bases.

    // Add the addtional filters and transformations that might have been set by the commands.
    for( auto const& func : combined_filters_and_transforms_ ) {
        iterator.add_transform_filter( func );
    }

    // In addition to the transforms and filters, we here also add a visitor function.
    // We always print out where the input is at, at the moment. That makes sure that we always
    // get some progress update, which is probably more useful for the user than waiting too long
    // without any.
    // Note though that region filters are applied per input file, and so might have already removed
    // chromosomes, so that they won't be printed here. We use lambda capture by value to create
    // a copy of current_chr that is kept in the lambda, and updated there. We are not in C++14 yet.
    std::string current_chr;
    iterator.add_visitor(
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
        iterator.add_visitor(
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
