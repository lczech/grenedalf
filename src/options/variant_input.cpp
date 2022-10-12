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

#include "options/variant_input.hpp"

#include "options/global.hpp"
#include "tools/misc.hpp"

#include "genesis/population/formats/variant_parallel_input_iterator.hpp"
#include "genesis/population/functions/filter_transform.hpp"
#include "genesis/population/functions/functions.hpp"
#include "genesis/utils/core/algorithm.hpp"
#include "genesis/utils/core/fs.hpp"
#include "genesis/utils/core/options.hpp"
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
    bool with_sample_name_opts,
    bool with_region_filter_opts,
    std::string const& group
) {

    // Add input file type options.
    input_sam_.add_sam_input_opt_to_app( sub );
    input_pileup_.add_pileup_input_opt_to_app( sub );
    input_sync_.add_sync_input_opt_to_app( sub );
    input_vcf_.add_vcf_input_opt_to_app( sub );

    // Multi file set operation. In the VariantParallelInputIterator, we have a different way
    // of expressing which loci of which input source to visit, but that would be way too complex
    // to specify via a command line interface. So, at least for now, we simply boil it down
    // to either the union of all loci, or their intersection.
    multi_file_loci_set_.option = sub->add_option(
        "--multi-file-locus-set",
        multi_file_loci_set_.value,
        "When multiple input files are provided, select whether the union of all their loci is "
        "used, or their intersection. For their union, input files that do not have data at a "
        "particular locus are considered to have zero counts at every base at that locus. "
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

    // Additional options.
    if( with_sample_name_opts ) {
        input_sample_names_.add_sample_name_opts_to_app( sub, group );
    }
    if( with_region_filter_opts ) {
        region_filter_.add_region_filter_opts_to_app( sub, group );
    }
}

// =================================================================================================
//      Command Setup Functions
// =================================================================================================

void VariantInputOptions::add_individual_filter_and_transforms(
    std::function<bool(genesis::population::Variant&)> const& func
) const {
    // Checks for internal correct setup
    if( static_cast<bool>( iterator_ ) || ! sample_names_.empty() ) {
        throw std::domain_error(
            "Internal error: Calling add_individual_filter_and_transforms() "
            "after input source iteration has already been started."
        );
    }
    individual_filters_and_transforms_.push_back( func );
}

void VariantInputOptions::add_combined_filter_and_transforms(
    std::function<bool(genesis::population::Variant&)> const& func
) const {
    // Checks for internal correct setup
    if( static_cast<bool>( iterator_ ) || ! sample_names_.empty() ) {
        throw std::domain_error(
            "Internal error: Calling add_combined_filter_and_transforms() "
            "after input source iteration has already been started."
        );
    }
    combined_filters_and_transforms_.push_back( func );
}

// =================================================================================================
//      Run Functions
// =================================================================================================

// -------------------------------------------------------------------------
//     sample_names
// -------------------------------------------------------------------------

std::vector<std::string> const& VariantInputOptions::sample_names() const
{
    prepare_data_();
    return sample_names_;
}

// -------------------------------------------------------------------------
//     get_iterator
// -------------------------------------------------------------------------

genesis::population::VariantInputIterator& VariantInputOptions::get_iterator() const
{
    prepare_data_();
    return iterator_;
}

// =================================================================================================
//      Iterator Setup
// =================================================================================================

// -------------------------------------------------------------------------
//     prepare_data_
// -------------------------------------------------------------------------

void VariantInputOptions::prepare_data_() const
{
    using namespace genesis;
    using namespace genesis::population;
    using namespace genesis::utils;

    // Checks for internal correct setup
    if( static_cast<bool>( iterator_ ) || ! sample_names_.empty() ) {
        // Nothing to be done. We already prepared the data.
        return;
    }

    // We first read all region filters, so that they can be re-used for all inputs.
    region_filter_.prepare_region_filters();

    // Check how many files are given. If it is a single one, we just use that for the input
    // iterator, which is faster than piping it through a parallel iterator. For multiple files,
    // we create a parallel iterator.
    size_t file_count = 0;
    if( input_sam_.get_file_input_options().provided() ) {
        file_count += input_sam_.get_file_input_options().file_count();
    }
    if( input_pileup_.get_file_input_options().provided() ) {
        file_count += input_pileup_.get_file_input_options().file_count();
    }
    if( input_sync_.get_file_input_options().provided() ) {
        file_count += input_sync_.get_file_input_options().file_count();
    }
    if( input_vcf_.get_file_input_options().provided() ) {
        file_count += input_vcf_.get_file_input_options().file_count();
    }

    // Set up iterator depending on how many files we found in total across all inputs.
    if( file_count == 0 ) {
        throw CLI::ValidationError(
            "At least one input file has to be provided."
        );
    } else if( file_count == 1 ) {
        prepare_data_single_file_();
    } else {
        prepare_data_multiple_files_();
    }
    assert( iterator_ );
    assert( ! sample_names_.empty() );

    // Some user output.
    LOG_MSG1 << "Processing " << sample_names_.size()
             << " sample" << ( sample_names_.size() == 1 ? "" : "s" );
    for( auto const& sn : sample_names_ ) {
        LOG_MSG2 << " - " << sn;
    }
    LOG_MSG1;
    if( contains_duplicates( sample_names_ )) {
        LOG_WARN << "The input contains duplicate sample names. We can work with that internally, "
                 << "but it will make working with the output more difficult. Use verbose output "
                 << "(--verbose) to show a list of all sample names.";
        LOG_MSG1;
    }
}

// -------------------------------------------------------------------------
//     prepare_data_single_file_
// -------------------------------------------------------------------------

void VariantInputOptions::prepare_data_single_file_() const
{
    // Assert that this function is only called in a context where the data is not yet prepared.
    internal_check(
        ! static_cast<bool>( iterator_ ) && sample_names_.empty(),
        "prepare_data_single_file_() called in an invalid context."
    );

    // Check that we have exactly one input files.
    auto const is_sam    = input_sam_.get_file_input_options().provided();
    auto const is_pileup = input_pileup_.get_file_input_options().provided();
    auto const is_sync   = input_sync_.get_file_input_options().provided();
    auto const is_vcf    = input_vcf_.get_file_input_options().provided();
    internal_check(
        ( is_sam + is_pileup + is_sync + is_vcf == 1 ) &&
        (
            input_sam_.get_file_input_options().file_count()    + input_pileup_.get_file_input_options().file_count() +
            input_sync_.get_file_input_options().file_count()   + input_vcf_.get_file_input_options().file_count()    == 1
        ),
        "prepare_data_single_file_() called with more than one file provided"
    );

    // If a sample name prefix or a list is given,
    // we check that this is only for the allowed file types.
    auto const is_sample_name_list = (
        input_sample_names_.sample_name_list_.option &&
        *input_sample_names_.sample_name_list_.option
    );
    auto const is_sample_name_pref = (
        input_sample_names_.sample_name_prefix_.option &&
        *input_sample_names_.sample_name_prefix_.option
    );
    if( is_sample_name_list || is_sample_name_pref ) {
        // Not both can be given, as the options are mutually exclusive.
        assert( is_sample_name_list ^ is_sample_name_pref );
        if( !( is_sam && ! input_sam_.sam_split_by_rg() ) && ! is_pileup && ! is_sync ) {
            throw CLI::ValidationError(
                "Can only use " + input_sample_names_.sample_name_list_.option->get_name() + " or " +
                input_sample_names_.sample_name_prefix_.option->get_name() + " for input file "
                "formats that do not already have sample sames, such as (m)pileup or sync files, "
                "or sam/bam/cram files without splitting by @RG read group tags (which then is "
                "considered as a single sample that contains all reads independent of their @RG)."
            );
        }
    }

    // Prepare the iterator depending on the input file format.
    if( input_sam_.get_file_input_options().provided() ) {
        assert( input_sam_.get_file_input_options().file_count() == 1 );
        iterator_ = input_sam_.prepare_sam_iterator(
            input_sam_.get_file_input_options().file_paths()[0],
            input_sample_names_
        );
    }
    if( input_pileup_.get_file_input_options().provided() ) {
        assert( input_pileup_.get_file_input_options().file_count() == 1 );
        iterator_ = input_pileup_.prepare_pileup_iterator(
            input_pileup_.get_file_input_options().file_paths()[0],
            input_sample_names_
        );
    }
    if( input_sync_.get_file_input_options().provided() ) {
        assert( input_sync_.get_file_input_options().file_count() == 1 );
        iterator_ = input_sync_.prepare_sync_iterator(
            input_sync_.get_file_input_options().file_paths()[0],
            input_sample_names_
        );
    }
    if( input_vcf_.get_file_input_options().provided() ) {
        assert( input_vcf_.get_file_input_options().file_count() == 1 );
        iterator_ = input_vcf_.prepare_vcf_iterator(
            input_vcf_.get_file_input_options().file_paths()[0],
            input_sample_names_
        );
    }

    // Copy over the sample names from the iterator, so that they are accessible.
    sample_names_ = iterator_.data().sample_names;
    if( sample_names_.empty() ) {
        throw std::runtime_error( "Invalid input file that does not contain any samples." );
    }
    internal_check(
        static_cast<bool>( iterator_ ) && ! sample_names_.empty(),
        "prepare_data_single_file_() call to file prepare function did not succeed."
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
    iterator_.thread_pool( genesis::utils::Options::get().global_thread_pool() );

    // Set the buffer block size of the iterator, for multi-threaded processing speed.
    // Needs to be tested - might give more or less advantage depending on setting.
    iterator_.block_size( iterator_block_size_.value );

    // TODO also add filters depending on file type. sam pileup and sync might need
    // biallelic snp filters (using the merged variants filter type setting), while vcf might not?!
    // or has it already set below or in the variant input iterator - need to check.
}

// -------------------------------------------------------------------------
//     prepare_data_multiple_files_
// -------------------------------------------------------------------------

void VariantInputOptions::prepare_data_multiple_files_() const
{
    using namespace genesis;
    using namespace genesis::population;
    using namespace genesis::utils;

    // Assert that this function is only called in a context where the data is not yet prepared.
    internal_check(
        ! static_cast<bool>( iterator_ ) && sample_names_.empty(),
        "prepare_data_multiple_files_() called in an invalid context."
    );

    // Check that we have multiple input files.
    internal_check(
        input_sam_.get_file_input_options().file_count()    +
        input_pileup_.get_file_input_options().file_count() +
        input_sync_.get_file_input_options().file_count()   +
        input_vcf_.get_file_input_options().file_count()    > 1,
        "prepare_data_multiple_files_() called with just one file provided"
    );

    // Sample name list cannot be used with multiple files.
    if( input_sample_names_.sample_name_list_.option && *input_sample_names_.sample_name_list_.option ) {
        throw CLI::ValidationError(
            "Can only use " + input_sample_names_.sample_name_list_.option->get_name() +
            " for single input files, but not when multiple input files are given."
        );
    }

    // No sample name filters can be given, as that would just be too tedious to specify via a CLI.
    if(
        ! input_sample_names_.filter_samples_include_.value.empty() ||
        ! input_sample_names_.filter_samples_exclude_.value.empty()
    ) {
        throw CLI::ValidationError(
            input_sample_names_.filter_samples_include_.option->get_name() + "(" +
            input_sample_names_.filter_samples_include_.value + "), " +
            input_sample_names_.filter_samples_exclude_.option->get_name() + "(" +
            input_sample_names_.filter_samples_exclude_.value + ")",
            "Can only use sample name filters for single input files, "
            "but not when multiple input files are given, "
            "as specifying filters per file via a command line interface is just too tedious."
        );
    }

    // Get whether the user wants the union or intersection of all parallel input loci.
    auto const contribution = get_enum_map_value(
        multi_file_contribution_type_map_, multi_file_loci_set_.value
    );

    // Make a parallel input iterator, and add all input from all file formats to it,
    // using the same contribution type for all of them, which either results in the union
    // or the intersection of all input loci. See VariantParallelInputIterator for details.
    VariantParallelInputIterator parallel_it;
    for( auto const& file : input_sam_.get_file_input_options().file_paths() ) {
        parallel_it.add_variant_input_iterator(
            input_sam_.prepare_sam_iterator( file, input_sample_names_ ), contribution
        );
    }
    for( auto const& file : input_pileup_.get_file_input_options().file_paths() ) {
        parallel_it.add_variant_input_iterator(
            input_pileup_.prepare_pileup_iterator( file, input_sample_names_ ), contribution
        );
    }
    for( auto const& file : input_sync_.get_file_input_options().file_paths() ) {
        parallel_it.add_variant_input_iterator(
            input_sync_.prepare_sync_iterator( file, input_sample_names_ ), contribution
        );
    }
    for( auto const& file : input_vcf_.get_file_input_options().file_paths() ) {
        parallel_it.add_variant_input_iterator(
            input_vcf_.prepare_vcf_iterator( file, input_sample_names_ ), contribution
        );
    }

    // Go through all sources again, build the sample names list from them,
    // and add the individual samples filters to all of them, e.g., so that regions are filtered
    // before the data even reaches the parallel iterator. Faster!
    assert( sample_names_.empty() );
    for( auto& input : parallel_it.inputs() ) {
        for( auto const& sample_name : input.data().sample_names ) {
            sample_names_.push_back( input.data().source_name + ":" + sample_name );
        }
        add_individual_filters_and_transforms_to_iterator_( input );

        // Following the reasoning from above: we need to use the global thread pool,
        // as otherwise reading hundreds of files will spawn too many threads, which is not good.
        input.thread_pool( genesis::utils::Options::get().global_thread_pool() );

        // Also set the buffer block size of the iterator, using the second buffer size option
        // that is used for the inner iterators of parallel input.
        // Needs to be tested - might give more or less advantage depending on setting.
        input.block_size( parallel_block_size_.value );
    }

    // Finally, create the parallel iterator.
    iterator_ = make_variant_input_iterator_from_variant_parallel_input_iterator( parallel_it );
    internal_check(
        iterator_.data().sample_names.size() == sample_names_.size(),
        "prepare_data_multiple_files_() with different number of samples in iterator and list"
    );

    // Add the filters and transformations that are to be applied to all samples combined.
    add_combined_filters_and_transforms_to_iterator_( iterator_ );

    // We also need to set the thread pool of the parallel input iterator itself,
    // for the same reasons as described above.
    iterator_.thread_pool( genesis::utils::Options::get().global_thread_pool() );

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
    if( region_filter_.get_region_filter() ) {
        iterator.add_filter( filter_by_region( region_filter_.get_region_filter() ));
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
}
