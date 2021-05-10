/*
    grenedalf - Genome Analyses of Differential Allele Frequencies
    Copyright (C) 2020-2021 Lucas Czech

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

#include "options/frequency_input.hpp"

#include "options/global.hpp"
#include "tools/misc.hpp"

#include "genesis/population/functions/genome_region.hpp"
#include "genesis/population/functions/variant.hpp"
#include "genesis/population/genome_region.hpp"
#include "genesis/utils/containers/filter_iterator.hpp"
#include "genesis/utils/core/fs.hpp"
#include "genesis/utils/text/string.hpp"

#include <algorithm>
#include <cassert>
#include <memory>
#include <stdexcept>

// =================================================================================================
//      Setup Functions
// =================================================================================================

void FrequencyInputOptions::add_frequency_input_opts_to_app(
    CLI::App* sub,
    // bool required,
    std::string const& group
) {
    // (void) required;
    add_vcf_input_opt_to_app( sub, false, group );
    add_pileup_input_opt_to_app( sub, false, group );
    vcf_file_.option->excludes( pileup_file_.option );
    pileup_file_.option->excludes( vcf_file_.option );

    // Correct setup check.
    internal_check(
        filter_region_.option == nullptr,
        "Cannot use the same FrequencyInputOptions object multiple times."
    );

    // Add option for genomic region filter.
    filter_region_.option = sub->add_option(
        "--filter-region",
        filter_region_.value,
        "Genomic region to filter for, in the format \"chr\", \"chr:position\", \"chr:start-end\", "
        "or \"chr:start..end\". If not provided, the whole input file is used."
    );
    filter_region_.option->group( group );

    // Add option for sample name filter.
    filter_samples_include_.option = sub->add_option(
        "--filter-samples-include",
        filter_samples_include_.value,
        "Sample names to include (all other samples are excluded); either a comma- or tab-separated "
        " list, or a file with one sample name per line. If no sample filter is provided, all "
        "samples in the input file are used."
    );
    filter_samples_include_.option->group( group );

    // And the other way round.
    filter_samples_exclude_.option = sub->add_option(
        "--filter-samples-exclude",
        filter_samples_exclude_.value,
        "Sample names to exclude (all other samples are included); either a comma- or tab-separated "
        " list, or a file with one sample name per line. If no sample filter is provided, all "
        "samples in the input file are used."
    );
    filter_samples_exclude_.option->group( group );
    filter_samples_exclude_.option->excludes( filter_samples_include_.option );
}

CLI::Option* FrequencyInputOptions::add_pileup_input_opt_to_app(
    CLI::App* sub,
    bool required,
    std::string const& group
) {
    // Correct setup check.
    internal_check(
        pileup_file_.option == nullptr,
        "Cannot use the same FrequencyInputOptions object multiple times."
    );

    // Add the option
    pileup_file_.option = sub->add_option(
        "--pileup-file",
        pileup_file_.value,
        "Path to an (m)pileup file."
    );
    pileup_file_.option->check( CLI::ExistingFile );
    pileup_file_.option->group( group );
    if( required ) {
        pileup_file_.option->required();
    }

    // Pileup does not have sample names, so offer a prefix option for this.
    pileup_sample_prefix_.option = sub->add_option(
        "--pileup-sample-prefix",
        pileup_sample_prefix_.value,
        "Pileup files do not contain sample names. This prefix followed by indices 1..n "
        "is used instead to provide unique names per sample that we use in the output and the "
        "`--filter-samples-include` and `--filter-samples-exclude` options. You can for example "
        "use \"Sample_\" as a prefix. If not provided, we simply use numbers 1..n as sample names "
        "in pileup files."
    );
    pileup_sample_prefix_.option->group( group );
    pileup_sample_prefix_.option->needs( pileup_file_.option );

    return pileup_sample_prefix_.option;
}

CLI::Option* FrequencyInputOptions::add_vcf_input_opt_to_app(
    CLI::App* sub,
    bool required,
    std::string const& group
) {
    // Correct setup check.
    internal_check(
        vcf_file_.option == nullptr,
        "Cannot use the same FrequencyInputOptions object multiple times."
    );

    // Add the option
    vcf_file_.option = sub->add_option(
        "--vcf-file",
        vcf_file_.value,
        "Path to a VCF file."
    );
    vcf_file_.option->check( CLI::ExistingFile );
    vcf_file_.option->group( group );
    if( required ) {
        vcf_file_.option->required();
    }

    return vcf_file_.option;
}

void FrequencyInputOptions::add_sliding_window_opts_to_app(
    CLI::App* sub,
    std::string const& group
) {
    // Width
    window_width_.option = sub->add_option(
        "--window-width",
        window_width_.value,
        "Width of each window along the chromosome."
    );
    window_width_.option->group( group );

    // Stride
    window_stride_.option = sub->add_option(
        "--window-stride",
        window_stride_.value,
        "Stride between windows along the chromosome, that is how far to move to get to the next "
        "window. If set to 0 (default), this is set to the same value as the `--window-width`."
    );
    window_stride_.option->group( group );
}

// =================================================================================================
//      Run Functions
// =================================================================================================

// -------------------------------------------------------------------------
//     sample_names
// -------------------------------------------------------------------------

std::vector<std::string> const& FrequencyInputOptions::sample_names() const
{
    prepare_data_();
    return sample_names_;
}

// -------------------------------------------------------------------------
//     get_iterator
// -------------------------------------------------------------------------

genesis::utils::Range<genesis::utils::LambdaIterator<genesis::population::Variant>>
FrequencyInputOptions::get_iterator() const
{
    // Return the range of the data that we just captured and converted.
    prepare_data_();
    return { generator_.begin(), generator_.end() };
}

// -------------------------------------------------------------------------
//     get_base_count_sliding_window_iterator
// -------------------------------------------------------------------------

genesis::population::SlidingWindowIterator<
    genesis::utils::LambdaIterator<genesis::population::Variant>,
    genesis::population::Variant,
    std::vector<genesis::population::BaseCounts>
>
FrequencyInputOptions::get_base_count_sliding_window_iterator() const
{
    using namespace genesis;
    using namespace genesis::population;

    // User-provided sliding window settings
    SlidingWindowIteratorSettings<Variant, std::vector<BaseCounts>> settings;
    settings.width = window_width_.value;
    settings.stride = window_stride_.value;

    // Conversion functions for the sliding window iterator.
    settings.entry_input_function = []( Variant const& variant ) -> std::vector<BaseCounts> const& {
        return variant.samples;
    };
    settings.chromosome_function = []( Variant const& variant ) -> std::string const& {
        return variant.chromosome;
    };
    settings.position_function = []( Variant const& variant ){
        return variant.position;
    };

    // Make sure that we have the iterator over the input file set up, and then return the
    // window iterator.
    prepare_data_();
    return make_sliding_window_iterator( settings, generator_.begin(), generator_.end() );
}

// -------------------------------------------------------------------------
//     get_variant_sliding_window_iterator
// -------------------------------------------------------------------------

genesis::population::SlidingWindowIterator<
    genesis::utils::LambdaIterator<genesis::population::Variant>,
    genesis::population::Variant
>
FrequencyInputOptions::get_variant_sliding_window_iterator() const
{
    using namespace genesis;
    using namespace genesis::population;

    // User-provided sliding window settings
    SlidingWindowIteratorSettings<Variant> settings;
    settings.width = window_width_.value;
    settings.stride = window_stride_.value;

    // Conversion functions for the sliding window iterator.
    settings.entry_input_function = []( Variant const& variant ) -> Variant const& {
        return variant;
    };
    settings.chromosome_function = []( Variant const& variant ) -> std::string const& {
        return variant.chromosome;
    };
    settings.position_function = []( Variant const& variant ){
        return variant.position;
    };

    // Make sure that we have the iterator over the input file set up, and then return the
    // window iterator.
    prepare_data_();
    return make_sliding_window_iterator( settings, generator_.begin(), generator_.end() );
}

// =================================================================================================
//      Internal Helpers
// =================================================================================================

// -------------------------------------------------------------------------
//     prepare_data_
// -------------------------------------------------------------------------

void FrequencyInputOptions::prepare_data_() const
{
    using namespace genesis;
    using namespace genesis::population;
    using namespace genesis::utils;

    // Checks for internal correct setup
    if( static_cast<bool>( generator_ ) || !sample_names_.empty() ) {
        // Nothing to be done.
        return;

        // throw std::runtime_error(
        //     "Internal error: Generator setup called multiple times."
        // );
    }

    // Check for correct input
    if( pileup_file_.value == "" && vcf_file_.value == "" ) {
        throw CLI::ValidationError(
            "Have to use either " + pileup_file_.option->get_name()  + " or " +
            vcf_file_.option->get_name() + " as input."
        );
    }
    if( pileup_file_.value != "" && vcf_file_.value != "" ) {
        throw CLI::ValidationError(
            "Cannot use both " + pileup_file_.option->get_name()  + " and " +
            vcf_file_.option->get_name() + " as input at the same time."
        );
    }
    assert(( pileup_file_.value == "" ) ^ ( vcf_file_.value == "" ));

    // Here, we need to select the different input sources and transform them into a uniform
    // iterator, using lambdas with std::function for type erasure. Additionally, we want region
    // filtering. Ideally, we would want to apply that afterwards, but that would just introduce
    // yet another level of indirection, as we'd need another lambda iterator to achieve that.
    // Sometimes, static types suck... Well, so instead, we apply the filter beforehand (which has
    // the additonal advantage that we do not need to convert irrelevant positions), but which
    // comes with the need to have individual filter iterators for each input data source...

    if( pileup_file_.value != "" ) {
        prepare_data_pileup_();
    }
    if( vcf_file_.value != "" ) {
        prepare_data_vcf_();
    }
}

// -------------------------------------------------------------------------
//     prepare_data_pileup_
// -------------------------------------------------------------------------

void FrequencyInputOptions::prepare_data_pileup_() const
{
    using namespace genesis;
    using namespace genesis::population;
    using namespace genesis::utils;

    // Assert that this function is only called in a context where the data is not yet prepared.
    assert( ! static_cast<bool>( generator_ ) && sample_names_.empty() );

    // Prepare the base Reader.
    auto reader = SimplePileupReader();
    // reader.quality_encoding( genesis::sequence::QualityEncoding::kIllumina13 );

    // Open the file, which aleady reads the first line. We use this to get the number of
    // samples in the pileup, and create dummy names for them.
    // We mgiht later open it again to incorporate the sample name filtering, because to do that,
    // we first need to know how many samples there are in total...
    // Maybe there is a smart way to avoid that, but for now, that seems easiest.
    auto it = SimplePileupInputIterator( utils::from_file( pileup_file_.value ), reader );
    auto const smp_cnt = it->samples.size();
    for( size_t i = 0; i < smp_cnt; ++i ) {
        sample_names_.push_back( pileup_sample_prefix_.value + std::to_string(i+1) );
    }

    // Filter sample names as needed. This is a bit cumbersome, but gets the job done.
    if( ! filter_samples_include_.value.empty() || ! filter_samples_exclude_.value.empty() ) {
        // Not both can be given at the same time, as we made the options mutually exclusive.
        assert( filter_samples_include_.value.empty() != filter_samples_exclude_.value.empty() );

        // Get the sample names and prepare a filter vector.
        // We init the filter so that in the include case, all are false first, and in the exclude
        // case, all are true first. Then, we can set the listed indices accordingly.
        // We also need to keep track of which indices we find, so that we can build the proper
        // name list afterwards. How cumbersome... maybe fix in the future.
        auto const list = get_sample_name_list( filter_samples_include_.value );
        auto sample_filter = std::vector<bool>( smp_cnt, filter_samples_exclude_.value.empty() );
        std::vector<size_t> sample_indices;
        sample_indices.reserve( list.size() );
        for( auto const& sn : list ) {
            auto it = std::find( sample_names_.begin(), sample_names_.end(), sn );
            if( it == sample_names_.end() ) {
                throw CLI::ValidationError(
                    "Invalid sample name used for filtering: \"" + sn  + "\"."
                );
            } else {
                auto const index = it - sample_names_.begin();
                sample_filter[ index ] = filter_samples_include_.value.empty();
                sample_indices.push_back( static_cast<size_t>( index ));
            }
        }

        // Now restart the iteration, this time with the filtering, and renew the sample names.
        it = SimplePileupInputIterator(
            utils::from_file( pileup_file_.value ), sample_filter, reader
        );
        sample_names_.clear();
        for( auto idx : sample_indices ) {
            sample_names_.push_back( pileup_sample_prefix_.value + std::to_string(idx+1) );
        }
    }

    // Apply region filter if necessary.
    if( filter_region_.value.empty() ) {
        // Create a generator that reads pileup.
        generator_ = LambdaIteratorGenerator<Variant>(
            [it]() mutable -> std::shared_ptr<Variant>{
                if( it ) {
                    auto res = std::make_shared<Variant>( convert_to_variant(*it) );
                    ++it;
                    return res;
                } else {
                    return nullptr;
                }
            }
        );
    } else {
        auto const region = parse_genome_region( filter_region_.value );
        auto region_filtered_range = make_filter_range(
            [region]( SimplePileupReader::Record const& record ){
                return genesis::population::is_covered( region, record );
            },
            // Use the iterator and a default constructed dummy as begin and end.
            it, SimplePileupInputIterator()
        );

        // Create a generator that reads pileup and filters by region.
        auto beg = region_filtered_range.begin();
        auto end = region_filtered_range.end();
        generator_ = LambdaIteratorGenerator<Variant>(
            [beg, end]() mutable -> std::shared_ptr<Variant>{
                if( beg != end ) {
                    // TODO this should be simply *beg here, but somewhere along the way,
                    // this would drop a const qualifier. need to figure out a neat way to fix.
                    auto res = std::make_shared<Variant>( convert_to_variant(*beg.base()) );
                    ++beg;
                    return res;
                } else {
                    return nullptr;
                }
            }
        );
    }
}

// -------------------------------------------------------------------------
//     prepare_data_vcf_
// -------------------------------------------------------------------------

void FrequencyInputOptions::prepare_data_vcf_() const
{
    using namespace genesis;
    using namespace genesis::population;
    using namespace genesis::utils;

    // Assert that this function is only called in a context where the data is not yet prepared.
    assert( ! static_cast<bool>( generator_ ) && sample_names_.empty() );

    // TODO buffer size

    // Prepare the base iterator.
    // See if we want to filter by sample name, and if so, resolve the name list.
    auto vcf_in = VcfInputIterator();
    if( ! filter_samples_include_.value.empty() ) {
        auto const list = get_sample_name_list( filter_samples_include_.value );
        vcf_in = VcfInputIterator( vcf_file_.value, list );
    } else if( ! filter_samples_exclude_.value.empty() ) {
        auto const list = get_sample_name_list( filter_samples_exclude_.value );
        vcf_in = VcfInputIterator( vcf_file_.value, list, true );
    } else {
        vcf_in = VcfInputIterator( vcf_file_.value );
    }
    if( !vcf_in.header().has_format("AD") ) {
        throw std::runtime_error(
            "Cannot use VCF input file that does not have the `AD` format field."
        );
    }

    // Get the sample names. This will only contain the filtered names.
    // Then, create a filter over the input that only allows biallelic SNPs with the AD format field.
    // Everything else we cannot use anyway for our subsequent conversion steps.
    sample_names_ = vcf_in.header().get_sample_names();
    auto vcf_range = utils::make_filter_range([]( VcfRecord const& record ){
        bool const is_snp = record.is_snp() && record.get_alternatives_count() == 1;
        bool const has_ad = record.has_format( "AD" );
        return is_snp && has_ad;
    }, vcf_in, {} );

    // Apply region filter if necessary.
    if( filter_region_.value.empty() ) {
        // Need variables that can be captures by copy...
        auto beg = vcf_range.begin();
        auto end = vcf_range.end();

        // Create an iterator that erases the type (which here is a complex templated type with
        // the VcfInputIterator and FilterIterator and all that). We use a lambda capture by
        // mutable value, so that the actual iterator is stored (as a copy, but that works) in the
        // lambda capture even after this whole function here is executed. That is, the generator_
        // keeps the underlying VCF iterator alive via the lambda capture. Wow.
        // Addendum: There was a nasty bug hidden in the implementation of VcfInputIterator,
        // as it uses a thread pool to pre-parse Vcf Records, but that thread pool was storing
        // invalid pointers after a race condition where the above setup deletes the original
        // copy of the iterator (`vcf_in`), but the thread would still want to use it...
        // That took a while to figure out, and is fixed now by having the thread pool keep copies
        // of the internal members of VcfFormatIterator of its own.
        generator_ = LambdaIteratorGenerator<Variant>(
            [beg, end]() mutable -> std::shared_ptr<Variant>{
                if( beg != end ) {
                    auto res = std::make_shared<Variant>( convert_to_variant(*beg) );
                    ++beg;
                    return res;
                } else {
                    return nullptr;
                }
            }
        );
    } else {
        auto const region = parse_genome_region( filter_region_.value );
        auto region_filtered_range = make_filter_range(
            [region]( VcfRecord const& record ){
                return genesis::population::is_covered( region, record );
            },
            vcf_range.begin(), vcf_range.end()
        );

        // Create a generator that reads pileup and filters by region.
        auto beg = region_filtered_range.begin();
        auto end = region_filtered_range.end();
        generator_ = LambdaIteratorGenerator<Variant>(
            [beg, end]() mutable -> std::shared_ptr<Variant>{
                if( beg != end ) {
                    auto res = std::make_shared<Variant>( convert_to_variant(*beg) );
                    ++beg;
                    return res;
                } else {
                    return nullptr;
                }
            }
        );
    }
}

// -------------------------------------------------------------------------
//     get_sample_name_list
// -------------------------------------------------------------------------

std::vector<std::string> FrequencyInputOptions::get_sample_name_list( std::string const& value ) const
{
    using namespace genesis::utils;

    // If the input is a file, read it line by line as sample names. Otherwise, split by comma.
    if( is_file( value ) ) {
        return file_read_lines( value );
    } else {
        return split( value, ",\t" );
    }
}
