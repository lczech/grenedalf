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

#include "genesis/population/formats/bed_reader.hpp"
#include "genesis/population/formats/gff_reader.hpp"
#include "genesis/population/formats/simple_pileup_input_iterator.hpp"
#include "genesis/population/formats/simple_pileup_reader.hpp"
#include "genesis/population/formats/sync_input_iterator.hpp"
#include "genesis/population/formats/sync_reader.hpp"
#include "genesis/population/formats/variant_parallel_input_iterator.hpp"
#include "genesis/population/formats/vcf_input_iterator.hpp"
#include "genesis/population/functions/filter_transform.hpp"
#include "genesis/population/functions/functions.hpp"
#include "genesis/population/functions/genome_region.hpp"
#include "genesis/population/genome_region.hpp"
#include "genesis/sequence/functions/quality.hpp"
#include "genesis/utils/containers/filter_iterator.hpp"
#include "genesis/utils/core/algorithm.hpp"
#include "genesis/utils/core/fs.hpp"
#include "genesis/utils/text/convert.hpp"
#include "genesis/utils/text/string.hpp"

#include <algorithm>
#include <cassert>
#include <memory>
#include <stdexcept>

// =================================================================================================
//      Setup Functions
// =================================================================================================

// -------------------------------------------------------------------------
//     All Input File Types
// -------------------------------------------------------------------------

void VariantInputOptions::add_frequency_input_opts_to_app(
    CLI::App* sub,
    // bool required,
    // bool with_sample_name_opts,
    // bool with_filter_opts,
    std::string const& group
) {
    // (void) required;
    (void) group;

    // Add input file type options.
    add_sam_input_opt_to_app( sub );
    add_pileup_input_opt_to_app( sub );
    add_sync_input_opt_to_app( sub );
    add_vcf_input_opt_to_app( sub );

    // Ha! Deactivated the below, as now we are allowing multiple input files!
    // Only one input file format allowed at a time, at least for now.
    // We could just list all pairs and call these options on all of them.
    // Or we can be fancy and do that in a nested loop over all pairs.
    // std::vector<CliOption<std::string>*> file_opt_ptrs {
    //     &sam_file_, &pileup_file_, &sync_file_, &vcf_file_
    // };
    // for( size_t i = 0; i < file_opt_ptrs.size(); ++i ) {
    //     for( size_t j = 0; j < file_opt_ptrs.size(); ++j ) {
    //         if( i != j ) {
    //             file_opt_ptrs[i]->option->excludes( file_opt_ptrs[j]->option );
    //         }
    //     }
    // }

    // Hidden option to set the LambdaIterator block size for speed.
    block_size_.option = sub->add_option(
        "--block-size",
        block_size_.value,
        "Size of the buffer block that is used for the multithreaded file parsing in the background. "
        "This option is an optimization feature that can be experiment with for larger datasets."
    );
    block_size_.option->group( "" );

    // // Additional options.
    // if( with_sample_name_opts ) {
    //     add_sample_name_opts_to_app( sub, group );
    // }
    // if( with_filter_opts ) {
    //     add_filter_opts_to_app( sub, group );
    // }
}

// -------------------------------------------------------------------------
//     SAM/BAM/CRAN Input
// -------------------------------------------------------------------------

CLI::Option* VariantInputOptions::add_sam_input_opt_to_app(
    CLI::App* sub,
    bool required,
    std::string const& group
) {
    // Correct setup check.
    internal_check(
        sam_file_.option() == nullptr,
        "Cannot use the same VariantInputOptions object multiple times."
    );

    // TODO add min_depth max_depth and max_accumulation_depth, and add FLAG to reaader and here.
    // but at least the first two should actually be general filter settings.

    // Add the option
    sam_file_.add_multi_file_input_opt_to_app(
        sub, "sam", "sam/bam/cram", "(sam(\\.gz)?|bam|cram)", "sam[.gz]|bam|cram", required, group
    );

    // Min mapping quality
    sam_min_map_qual_.option = sub->add_option(
        "--sam-min-map-qual",
        sam_min_map_qual_.value,
        "Minimum phred-scaled mapping quality score [0-90] for a read in sam/bam/cram files to be "
        "considered. Any read that is below the given value of mapping quality will be completely "
        "discarded, and its bases not taken into account. "
        "Default is 0, meaning no filtering by base quality qual."
    );
    sam_min_map_qual_.option->group( group );
    sam_min_map_qual_.option->check( CLI::Range( static_cast<size_t>(0), static_cast<size_t>(90) ));
    sam_min_map_qual_.option->needs( sam_file_.option() );

    // Min base qual
    sam_min_base_qual_.option = sub->add_option(
        "--sam-min-base-qual",
        sam_min_base_qual_.value,
        "Minimum phred-scaled quality score [0-90] for a base in sam/bam/cram files to be "
        "considered. Bases below this are ignored when computing allele frequencies. "
        "Default is 0, meaning no filtering by base quality qual."
    );
    sam_min_base_qual_.option->group( group );
    sam_min_base_qual_.option->check( CLI::Range( static_cast<size_t>(0), static_cast<size_t>(90) ));
    sam_min_base_qual_.option->needs( sam_file_.option() );

    // Split by RG read group tag.
    sam_split_by_rg_.option = sub->add_flag(
        "--sam-split-by-rg",
        sam_split_by_rg_.value,
        "Instead of considering the whole sam/bam/cram file as one large colletion of reads, "
        "use the `@RG` read group tag to split reads. Each read group is then considered a sample. "
        "Reads with an invalid (not in the header) read group tag or without a tag are ignored."
    );
    sam_split_by_rg_.option->group( group );
    sam_split_by_rg_.option->needs( sam_file_.option() );

    return sam_file_.option();
}

// -------------------------------------------------------------------------
//     Pileup Input
// -------------------------------------------------------------------------

CLI::Option* VariantInputOptions::add_pileup_input_opt_to_app(
    CLI::App* sub,
    bool required,
    std::string const& group
) {
    // Correct setup check.
    internal_check(
        pileup_file_.option() == nullptr,
        "Cannot use the same VariantInputOptions object multiple times."
    );

    // TODO add options for reading: with quality, with ancestral base

    // Add the option
    pileup_file_.add_multi_file_input_opt_to_app(
        sub, "pileup", "(m)pileup",
        "(plp|mplp|pileup|mpileup)(\\.gz)?",
        "(plp|mplp|pileup|mpileup)[.gz]",
        required, group
    );

    // Quality encoding.
    pileup_quality_encoding_.option = sub->add_option(
        "--pileup-quality-encoding",
        pileup_quality_encoding_.value,
        "Encoding of the quality scores of the bases in (m)pileup files. "
        "Default is `\"sanger\"`, which seems to be the most common these days. "
        "Both `\"sanger\"` and `\"illumina-1.8\"` are identical and use an ASCII offset of 33, "
        "while `\"illumina-1.3\"` and `\"illumina-1.5\"` are identical with an ASCII offset of 64 "
        "(we provide different names for completeness). Lastly, `\"solexa\"` has an offset of 64, "
        "but uses a different equation (not phred score) for the encoding."
    );
    pileup_quality_encoding_.option->group( group );
    pileup_quality_encoding_.option->transform(
        CLI::IsMember(
            { "sanger", "illumina-1.3", "illumina-1.5", "illumina-1.8", "solexa" },
            CLI::ignore_case
        )
    );
    pileup_quality_encoding_.option->needs( pileup_file_.option() );

    // Min phred score
    pileup_min_base_qual_.option = sub->add_option(
        "--pileup-min-base-qual",
        pileup_min_base_qual_.value,
        "Minimum phred quality score [0-90] for a base in (m)pileup files to be considered. "
        "Bases below this are ignored when computing allele frequencies. "
        "Default is 0, meaning no filtering by phred quality score."
    );
    pileup_min_base_qual_.option->group( group );
    pileup_min_base_qual_.option->check( CLI::Range( static_cast<size_t>(0), static_cast<size_t>(90) ));
    pileup_min_base_qual_.option->needs( pileup_file_.option() );

    return pileup_file_.option();
}

// -------------------------------------------------------------------------
//     Sync Input
// -------------------------------------------------------------------------

CLI::Option* VariantInputOptions::add_sync_input_opt_to_app(
    CLI::App* sub,
    bool required,
    std::string const& group
) {
    // Correct setup check.
    internal_check(
        sync_file_.option() == nullptr,
        "Cannot use the same VariantInputOptions object multiple times."
    );

    // Add the option
    sync_file_.add_multi_file_input_opt_to_app(
        sub, "sync", "sync (as specified by PoPoolation2)",
        "sync(\\.gz)?", "sync[.gz]",
        required, group
    );

    return sync_file_.option();
}

// -------------------------------------------------------------------------
//     VCF Input
// -------------------------------------------------------------------------

CLI::Option* VariantInputOptions::add_vcf_input_opt_to_app(
    CLI::App* sub,
    bool required,
    std::string const& group
) {
    // Correct setup check.
    internal_check(
        vcf_file_.option() == nullptr,
        "Cannot use the same VariantInputOptions object multiple times."
    );

    // Add the option
    vcf_file_.add_multi_file_input_opt_to_app(
        sub, "vcf", "vcf/bcf", "(vcf(\\.gz)?|bcf)", "vcf[.gz]|bcf", required, group,
        "This expects that the input file has the per-sample VCF FORMAT field `AD` (alleleic depth) "
        "given, containing the counts of the reference and alternative base. "
        "This assumes that the data that was used to create the VCF file was actually a pool of "
        "individuals (e.g., from pool sequencing) for each sample (column) of the VCF file. "
        "In this context, the `AD` field can then be interpreted as describing the allele "
        "frequencines of each pool of individuals."
    );

    return vcf_file_.option();
}

// -------------------------------------------------------------------------
//     Sample Naming
// -------------------------------------------------------------------------

void VariantInputOptions::add_sample_name_opts_to_app(
    CLI::App* sub,
    std::string const& group
) {
    // Correct setup check.
    internal_check(
        sample_name_prefix_.option == nullptr && sample_name_list_.option == nullptr,
        "Cannot use the same VariantInputOptions object multiple times."
    );

    // Name list option.
    sample_name_list_.option = sub->add_option(
        "--sample-name-list",
        sample_name_list_.value,
        "Some file types do not contain sample names, such as (m)pileup or sync files. For such "
        "file types, sample names can here be provided as either (1) a comma- or tab-separated "
        "list, or (2) as a file with one sample name per line, in the same order as samples are in "
        "the actual input file. We then use these names in the output and the "
        "`--filter-samples-include` and `--filter-samples-exclude` options. "
        "If not provided, we simply use numbers 1..n as sample names for these files types. "
        "Note that this option can only be used if a single file is given as input."
        "Alternatively, use `--sample-name-prefix` to provide a prefix for this sample numbering."
    );
    sample_name_list_.option->group( group );
    // sample_name_list_.option->check( CLI::ExistingFile );

    // Prefix option.
    sample_name_prefix_.option = sub->add_option(
        "--sample-name-prefix",
        sample_name_prefix_.value,
        "Some file types do not contain sample names, such as (m)pileup or sync files. For such "
        "file types, this prefix followed by indices 1..n can be used instead to provide unique "
        "names per sample that we use in the output and the `--filter-samples-include` and "
        "`--filter-samples-exclude` options. For example, use \"Sample_\" as a prefix. "
        "If not provided, we simply use numbers 1..n as sample names for these files types. "
        "This prefix also works if multiple files are given as input."
        "Alternatively, use `--sample-name-list` to directly provide a list of sample names."
    );
    sample_name_prefix_.option->group( group );
    // sample_name_prefix_.option->needs( pileup_file_.option );

    // The two ways of specifying sample names are mutually exclusive.
    sample_name_list_.option->excludes( sample_name_prefix_.option );
    sample_name_prefix_.option->excludes( sample_name_list_.option );
}

// -------------------------------------------------------------------------
//     Filtering
// -------------------------------------------------------------------------

void VariantInputOptions::add_filter_opts_to_app(
    CLI::App* sub,
    std::string const& group
) {
    // Correct setup check.
    internal_check(
        filter_region_.option == nullptr,
        "Cannot use the same VariantInputOptions object multiple times."
    );

    // Add option for genomic region filter.
    filter_region_.option = sub->add_option(
        "--filter-region",
        filter_region_.value,
        "Genomic region to filter for, in the format \"chr\", \"chr:position\", \"chr:start-end\", "
        "or \"chr:start..end\"."
    );
    filter_region_.option->group( group );

    // Add option for genomic region filter by BED file
    filter_region_bed_.option = sub->add_option(
        "--filter-region-bed",
        filter_region_bed_.value,
        "Genomic regions to filter for, as a BED file. In its simplest form, this file may contain "
        "a list of regions, one per line, with tab-separated chromosome and start and end positions. "
        // "Note that BED uses 0-based positions, and a half-open `[)` interval for the end position."
    );
    filter_region_bed_.option->check( CLI::ExistingFile );
    filter_region_bed_.option->group( group );

    // Add option for genomic region filter by BED file
    filter_region_gff_.option = sub->add_option(
        "--filter-region-gff",
        filter_region_gff_.value,
        "Genomic regions to filter for, as a GFF2/GFF3/GTF file. This only uses the chromosome, "
        "as well as start and end information per line, and ignores everything else in the file."
    );
    filter_region_gff_.option->check( CLI::ExistingFile );
    filter_region_gff_.option->group( group );

    // Region filter could be mutually exclusive... but why not allow all of them.
    // filter_region_.option->excludes( filter_region_bed_.option );
    // filter_region_.option->excludes( filter_region_gff_.option );
    // filter_region_bed_.option->excludes( filter_region_.option );
    // filter_region_bed_.option->excludes( filter_region_gff_.option );
    // filter_region_gff_.option->excludes( filter_region_.option );
    // filter_region_gff_.option->excludes( filter_region_bed_.option );

    // Add option for sample name filter.
    filter_samples_include_.option = sub->add_option(
        "--filter-samples-include",
        filter_samples_include_.value,
        "Sample names to include (all other samples are excluded); either (1) a comma- or "
        "tab-separated list, or (2) a file with one sample name per line. If no sample filter "
        "is provided, all samples in the input file are used. The option considers "
        "`--sample-name-list` or `--sample-name-prefix` for file types that do not contain sample "
        "names. Note that this option can only be used if a single file is given as input."
    );
    filter_samples_include_.option->group( group );

    // And the other way round.
    filter_samples_exclude_.option = sub->add_option(
        "--filter-samples-exclude",
        filter_samples_exclude_.value,
        "Sample names to exclude (all other samples are included); either (1) a comma- or "
        "tab-separated list, or (2) a file with one sample name per line. If no sample filter "
        "is provided, all samples in the input file are used. The option considers "
        "`--sample-name-list` or `--sample-name-prefix` for file types that do not contain sample "
        "names. Note that this option can only be used if a single file is given as input."
    );
    filter_samples_exclude_.option->group( group );

    // Include and exclude are mutually exclusive.
    filter_samples_exclude_.option->excludes( filter_samples_include_.option );
    filter_samples_include_.option->excludes( filter_samples_exclude_.option );
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

    // Check how many files are given. If it is a single one, we just use that for the input
    // iterator, which is faster than piping it through a parallel iterator. For multiple files,
    // we create a parallel iterator.
    size_t file_count = 0;
    if( sam_file_.provided() ) {
        file_count += sam_file_.file_count();
    }
    if( pileup_file_.provided() ) {
        file_count += pileup_file_.file_count();
    }
    if( sync_file_.provided() ) {
        file_count += sync_file_.file_count();
    }
    if( vcf_file_.provided() ) {
        file_count += vcf_file_.file_count();
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
        LOG_WARN << "The input contains duplicate sample names. ""We can work with that internally, "
                 << "but it will make working with the output more difficult.";
        LOG_MSG1;
    }
}

// -------------------------------------------------------------------------
//     prepare_data_single_file_
// -------------------------------------------------------------------------

void VariantInputOptions::prepare_data_single_file_() const
{
    // Check that we have exactly one input files.
    auto const is_sam    = sam_file_.provided();
    auto const is_pileup = pileup_file_.provided();
    auto const is_sync   = sync_file_.provided();
    auto const is_vcf    = vcf_file_.provided();
    internal_check(
        ( is_sam + is_pileup + is_sync + is_vcf == 1 ) &&
        (
            sam_file_.file_count()  + pileup_file_.file_count() +
            sync_file_.file_count() + vcf_file_.file_count()    == 1
        ),
        "prepare_data_single_file_() called with more than one file provided"
    );

    // If a sample name prefix or a list is given,
    // we check that this is only for the allowed file types.
    auto const is_sample_name_list = ( sample_name_list_.option   && *sample_name_list_.option );
    auto const is_sample_name_pref = ( sample_name_prefix_.option && *sample_name_prefix_.option );
    if( is_sample_name_list || is_sample_name_pref ) {
        // Not both can be given, as the options are mutually exclusive.
        assert( is_sample_name_list ^ is_sample_name_pref );
        if( !( is_sam && ! sam_split_by_rg_.value ) && ! is_pileup && ! is_sync ) {
            throw CLI::ValidationError(
                "Can only use " + sample_name_list_.option->get_name() + " or " +
                sample_name_prefix_.option->get_name() + " for input file "
                "formats that do not already have sample sames, such as (m)pileup or sync files, "
                "or sam/bam/cram files without splitting by @RG read group tags (which then is "
                "considered as a single sample that contains all reads independent of their @RG)."
            );
        }
    }

    // Prepare the iterator depending on the input file format.
    if( sam_file_.provided() ) {
        iterator_ = prepare_sam_iterator_( sam_file_.file_paths()[0] );
    }
    if( pileup_file_.provided() ) {
        iterator_ = prepare_pileup_iterator_( pileup_file_.file_paths()[0] );
    }
    if( sync_file_.provided() ) {
        iterator_ = prepare_sync_iterator_( sync_file_.file_paths()[0] );
    }
    if( vcf_file_.provided() ) {
        iterator_ = prepare_vcf_iterator_( vcf_file_.file_paths()[0] );
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
    add_region_filters_to_iterator_( iterator_ );

    // Set the buffer block size of the iterator, for multi-threaded processing speed.
    // Needs to be tested - might give more or less advantage depending on setting.
    iterator_.block_size( block_size_.value );

    // TODO also add filters depending on file type.. sam pileup and sync might need
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

    // Check that we have multiple input files.
    internal_check(
        sam_file_.file_count()  + pileup_file_.file_count() +
        sync_file_.file_count() + vcf_file_.file_count()   > 1,
        "prepare_data_multiple_files_() called with just one file provided"
    );

    // Sample name list cannot be used with multiple files.
    if( sample_name_list_.option && *sample_name_list_.option ) {
        throw CLI::ValidationError(
            "Can only use " + sample_name_list_.option->get_name() + " for single input files, "
            "but not when multiple input files are given."
        );
    }

    // No sample name filters can be given, as that would just be too tedious to specify via a CLI.
    if( ! filter_samples_include_.value.empty() || ! filter_samples_exclude_.value.empty() ) {
        throw CLI::ValidationError(
            filter_samples_include_.option->get_name() + "(" + filter_samples_include_.value + "), " +
            filter_samples_exclude_.option->get_name() + "(" + filter_samples_exclude_.value + ")",
            "Can only use sample name filters for single input files, "
            "but not when multiple input files are given, "
            "as specifying filters per file via a command line interface is just too tedious."
        );
    }

    // Make a parallel input iterator, and add all input from all file formats to it.
    // We currently add all inputs as "carrying", meaning that _all_ their loci are considered,
    // and not only those that also occurr in other input files.
    VariantParallelInputIterator parallel_it;
    auto const contribution = VariantParallelInputIterator::ContributionType::kCarrying;
    for( auto const& file : sam_file_.file_paths() ) {
        parallel_it.add_variant_input_iterator( prepare_sam_iterator_( file ), contribution );
    }
    for( auto const& file : pileup_file_.file_paths() ) {
        parallel_it.add_variant_input_iterator( prepare_pileup_iterator_( file ), contribution );
    }
    for( auto const& file : sync_file_.file_paths() ) {
        parallel_it.add_variant_input_iterator( prepare_sync_iterator_( file ), contribution );
    }
    for( auto const& file : vcf_file_.file_paths() ) {
        parallel_it.add_variant_input_iterator( prepare_vcf_iterator_( file ), contribution );
    }

    // Go through all sources again, build the sample names list from them,
    // and add the region filters to all of them, so that regions are filtered before
    // the data even reaches the parallel iterator. Faster!
    assert( sample_names_.empty() );
    for( auto& input : parallel_it.inputs() ) {
        for( auto const& sample_name : input.data().sample_names ) {
            sample_names_.push_back( input.data().source_name + ":" + sample_name );
        }
        add_region_filters_to_iterator_( input );
    }

    // Finally, create the parallel iterator.
    iterator_ = make_variant_input_iterator_from_variant_parallel_input_iterator( parallel_it );
    internal_check(
        iterator_.data().sample_names.size() == sample_names_.size(),
        "prepare_data_multiple_files_() with different number of samples in iterator and list"
    );

    // Set the buffer block size of the iterator, for multi-threaded speed.
    // We only buffer the final parallel iterator, and not its individual sources,
    // in order to not keep too many threads from spawning... Might need testing and refinement.
    iterator_.block_size( block_size_.value );
}

// -------------------------------------------------------------------------
//     prepare_sam_iterator_
// -------------------------------------------------------------------------

VariantInputOptions::VariantInputIterator VariantInputOptions::prepare_sam_iterator_(
    std::string const& filename
) const {
    using namespace genesis;
    using namespace genesis::population;
    using namespace genesis::utils;

    // Assert that this function is only called in a context where the data is not yet prepared.
    internal_check(
        ! static_cast<bool>( iterator_ ) && sample_names_.empty(),
        "prepare_sam_iterator_() called in an invalid context."
    );

    // Prepare the reader with all its settings.
    SamVariantInputIterator reader;
    reader.min_map_qual( sam_min_map_qual_.value );
    reader.min_base_qual( sam_min_base_qual_.value );
    reader.split_by_rg( sam_split_by_rg_.value );

    if( sam_split_by_rg_.value ) {

        // If we split by RG tag, we can use filters. Check which ones are given by the user.
        if( ! filter_samples_include_.value.empty() ) {
            auto const list = process_sample_name_list_option_( filter_samples_include_.value );
            reader.rg_tag_filter( std::unordered_set<std::string>{ list.begin(), list.end() });
        } else if( ! filter_samples_exclude_.value.empty() ) {
            auto const list = process_sample_name_list_option_( filter_samples_exclude_.value );
            reader.rg_tag_filter( std::unordered_set<std::string>{ list.begin(), list.end() });
            reader.inverse_rg_tag_filter( true );
        }

    } else {

        // Assert that no filters are set, as this does not make sense with just one sample.
        // This cannot really be done via CLI11 settings directly, as we'd need to excludes()
        // the case that the --split-by-rg was _not_ set... So let's do this here.
        auto const sample_filter = find_sample_indices_from_sample_filters_();
        if( ! sample_filter.first.empty() ) {
            throw CLI::ValidationError(
                sam_split_by_rg_.option->get_name() + ", " +
                filter_samples_include_.option->get_name() +
                "(" + filter_samples_include_.value + "), " +
                filter_samples_exclude_.option->get_name() +
                "(" + filter_samples_exclude_.value + ")",
                "Cannot use sample name filters for sam/bam/cram files without splitting by "
                "@RG read group tags. Without splitting, the file is considered as a single sample, "
                "in which case filtering does not make sense."
            );
        }
    }

    // Make an iterator.
    auto iterator = make_variant_input_iterator_from_sam_file( filename, reader );

    // If we do not split by RG tag, there is only a single sample. But still let's name it.
    // make_variant_input_iterator_from_sam_file() with rg splitting returns a list with as many
    // empty strings as the file has samples (exactly one in that case).
    if( ! sam_split_by_rg_.value ) {
        assert( iterator.data().sample_names.size() == 1 );
        iterator.data().sample_names = make_anonymous_sample_names_(
            iterator.data().sample_names.size()
        );
    }

    return iterator;
}

// -------------------------------------------------------------------------
//     prepare_pileup_iterator_
// -------------------------------------------------------------------------

VariantInputOptions::VariantInputIterator VariantInputOptions::prepare_pileup_iterator_(
    std::string const& filename
) const {
    using namespace genesis;
    using namespace genesis::population;
    using namespace genesis::sequence;
    using namespace genesis::utils;

    // Assert that this function is only called in a context where the data is not yet prepared.
    internal_check(
        ! static_cast<bool>( iterator_ ) && sample_names_.empty(),
        "prepare_pileup_iterator_() called in an invalid context."
    );

    // TODO min_phred_score is the old name. we currently use it here to distinguish it from
    // the setting of sam, but a better way should be found in the future.

    // We can use the sample filter settings to obtain a list of indices of samples
    // that we want to restrict the reading to. If no filter is given, that list is empty.
    // The second value of the returned pair indicates whether the list in inversed.
    auto const sample_filter = find_sample_indices_from_sample_filters_();

    // Prepare the base Reader with settings as needed.
    auto reader = SimplePileupReader();
    reader.quality_encoding( guess_quality_encoding_from_name( pileup_quality_encoding_.value ));
    reader.min_base_quality( pileup_min_base_qual_.value );

    // Make an iterator.
    auto iterator = make_variant_input_iterator_from_pileup_file(
        filename, sample_filter.first, sample_filter.second, reader
    );

    // Pileup does not have sample names, so set them based on user input or simple enumeration.
    // make_variant_input_iterator_from_pileup_file() returns a list with as many
    // empty strings as the file has samples (exactly one in that case).
    iterator.data().sample_names = make_anonymous_sample_names_(
        iterator.data().sample_names.size()
    );

    return iterator;
}

// -------------------------------------------------------------------------
//     prepare_sync_iterator_
// -------------------------------------------------------------------------

VariantInputOptions::VariantInputIterator VariantInputOptions::prepare_sync_iterator_(
    std::string const& filename
) const {
    using namespace genesis;
    using namespace genesis::population;
    using namespace genesis::utils;

    // Assert that this function is only called in a context where the data is not yet prepared.
    internal_check(
        ! static_cast<bool>( iterator_ ) && sample_names_.empty(),
        "prepare_sync_iterator_() called in an invalid context."
    );

    // We can use the sample filter settings, see prepare_pileup_iterator_() for details.
    auto const sample_filter = find_sample_indices_from_sample_filters_();

    // Make an iterator.
    auto iterator = make_variant_input_iterator_from_sync_file(
        filename, sample_filter.first, sample_filter.second
    );

    // Sync does not have sample names, so set them based on user input or simple enumeration.
    // make_variant_input_iterator_from_sync_file() returns a list with as many
    // empty strings as the file has samples (exactly one in that case).
    iterator.data().sample_names = make_anonymous_sample_names_(
        iterator.data().sample_names.size()
    );

    return iterator;
}

// -------------------------------------------------------------------------
//     prepare_vcf_iterator_
// -------------------------------------------------------------------------

VariantInputOptions::VariantInputIterator VariantInputOptions::prepare_vcf_iterator_(
    std::string const& filename
) const {
    using namespace genesis;
    using namespace genesis::population;
    using namespace genesis::utils;

    // Assert that this function is only called in a context where the data is not yet prepared.
    internal_check(
        ! static_cast<bool>( iterator_ ) && sample_names_.empty(),
        "prepare_vcf_iterator_() called in an invalid context."
    );

    // Prepare the iterator.
    // See if we want to filter by sample name, and if so, resolve the name list.
    // By default, this also already filters for biallelic SNPs.
    VariantInputIterator iterator;
    bool const only_biallelic = true;
    if( ! filter_samples_include_.value.empty() ) {
        auto const list = process_sample_name_list_option_( filter_samples_include_.value );
        iterator = make_variant_input_iterator_from_pool_vcf_file(
            filename, list, false, only_biallelic
        );
    } else if( ! filter_samples_exclude_.value.empty() ) {
        auto const list = process_sample_name_list_option_( filter_samples_exclude_.value );
        iterator = make_variant_input_iterator_from_pool_vcf_file(
            filename, list, true, only_biallelic
        );
    } else {
        iterator = make_variant_input_iterator_from_pool_vcf_file(
            filename, only_biallelic
        );
    }

    // As opposed to the above file formats, VCF contains sample names (only the filtered ones).
    // So here we do not need to set them, and can directly return.
    return iterator;
}

// =================================================================================================
//     Region Filters
// =================================================================================================

// -------------------------------------------------------------------------
//     add_region_filters_to_iterator_
// -------------------------------------------------------------------------

void VariantInputOptions::add_region_filters_to_iterator_(
    genesis::population::VariantInputIterator& iterator
) const {
    using namespace genesis;
    using namespace genesis::population;
    using namespace genesis::utils;

    // Add the region filter. This is the first one we do, so that all others do not need to be
    // applied to positions that are going to be filtered out anyway.
    if( ! filter_region_.value.empty() ) {
        auto const region = parse_genome_region( filter_region_.value );
        iterator.add_filter( filter_by_region( region ));
    }

    // Apply the region filter by bed file.
    // We use a shared ptr that stays alive in the lambda produced by filter_by_region().
    if( ! filter_region_bed_.value.empty() ) {
        auto const region_list = std::make_shared<GenomeRegionList>(
            BedReader().read_as_genome_region_list( from_file( filter_region_bed_.value ))
        );
        iterator.add_filter( filter_by_region( region_list ));
    }

    // Apply the region filter by gff file. Same as above, just different reader.
    if( ! filter_region_gff_.value.empty() ) {
        auto const region_list = std::make_shared<GenomeRegionList>(
            GffReader().read_as_genome_region_list( from_file( filter_region_gff_.value ))
        );
        iterator.add_filter( filter_by_region( region_list ));
    }
}

// =================================================================================================
//     Sample Name Processing
// =================================================================================================

// -------------------------------------------------------------------------
//     process_sample_name_list_option_
// -------------------------------------------------------------------------

std::vector<std::string> VariantInputOptions::process_sample_name_list_option_(
    std::string const& list
) const {
    using namespace genesis::utils;

    // If the input is a file, read it line by line as sample names.
    // Otherwise, split by comma or tab.
    std::vector<std::string> result;
    if( is_file( list ) ) {
        result = file_read_lines( list );
    } else {
        result = split( list, ",\t" );
    }
    if( result.empty() ) {
        throw std::runtime_error( "Invalid empty list of sample names given." );
    }

    return result;
}

// -------------------------------------------------------------------------
//     make_anonymous_sample_names_
// -------------------------------------------------------------------------

std::vector<std::string> VariantInputOptions::make_anonymous_sample_names_(
    size_t sample_count
) const {
    // Edge case, just to make sure.
    if( sample_count == 0 ) {
        throw std::runtime_error( "Input file does not contain any samples." );
    }

    // Prepare
    std::vector<std::string> result;

    // In case we have a proper list of sample names, use that.
    if( sample_name_list_.option && *sample_name_list_.option ) {
        result = process_sample_name_list_option_( sample_name_list_.value );
        if( result.size() != sample_count ) {
            throw CLI::ValidationError(
                sample_name_list_.option->get_name() + "(" + sample_name_list_.value + ")",
                "Invalid sample names list that contains " + std::to_string( result.size() ) +
                " name entries. This is incongruent with the input file, which contains " +
                std::to_string( sample_count ) + " samples."
            );
        }
        assert( result.size() > 0 );
        return result;
    }

    // In case we have a prefix for sample names instead, or nothing, just enumerate.
    std::string prefix;
    if( sample_name_prefix_.option && *sample_name_prefix_.option ) {
        prefix = sample_name_prefix_.value;
    }
    result.reserve( sample_count );
    for( size_t i = 0; i < sample_count; ++i ) {
        result.emplace_back( prefix + std::to_string( i + 1 ));
    }
    assert( result.size() > 0 );
    return result;
}

// -------------------------------------------------------------------------
//     find_sample_indices_from_sample_filters_
// -------------------------------------------------------------------------

std::pair<std::vector<size_t>, bool>
VariantInputOptions::find_sample_indices_from_sample_filters_() const
{
    // Prepare result.
    std::vector<size_t> indices;

    // Get whether we want to include or exclude sample names.
    bool const is_include = ! filter_samples_include_.value.empty();
    bool const is_exclude = ! filter_samples_exclude_.value.empty();

    // Not both can be given at the same time, as we made the options mutually exclusive.
    internal_check( !( is_include && is_exclude ), "include and exclude filters are both given" );

    // If no filters are given, just return empty, indicating to the caller that we do not filter.
    // This is also the return of this function when multiple files are given, in which case
    // we cannot apply sample filtering... that would just be too tedious to provide via a CLI.
    if( ! is_include && ! is_exclude ) {
        assert( indices.empty() );
        return { indices, false };
    }

    // Get the sample names, depending on which type (inc/exc) we have.
    auto const filter_list = process_sample_name_list_option_(
        is_include ? filter_samples_include_.value : filter_samples_exclude_.value
    );

    // When this function is called, we do not have opened the input file yet, so we do not
    // know how many samples there are. So here, we either need to get the list of sample names
    // (if the user provided one), or use enumeration to find the indices. The list loading is
    // duplicate work unfortunately, as it will be loaded again later. Alternatively, we could
    // split this into the case where a list of names is given, and the case where enumeration is
    // used, but that would make set_anonymous_sample_names_() and its setup more complex...
    // So let's live with loading the list twice.

    if( sample_name_list_.option && *sample_name_list_.option ) {
        // In case we have a proper list of sample names, use that, and check the filters against it.
        auto const sample_names = process_sample_name_list_option_( sample_name_list_.value );
        assert( sample_names.size() > 0 );

        // Go through all provided sample names for the filter, find their positon in the sample
        // name list, and fill the indices with these positions.
        for( auto const& filtered_name : filter_list ) {
            auto const it = std::find( sample_names.begin(), sample_names.end(), filtered_name );
            if( it == sample_names.end() ) {
                // Filter name not found, throw. First, helpful user output.
                LOG_MSG << "Sample names:";
                for( auto const& sn : sample_names ) {
                    LOG_MSG << " - " << sn;
                }
                throw CLI::ValidationError(
                    "Invalid sample name used for filtering: \"" + filtered_name  + "\"."
                );
            } else {
                // Add filter index to result.
                auto const index = it - sample_names.begin();
                indices.push_back( index );
            }
        }
    } else {
        // If we do not have a given list of sample names, we instead use the prefix and enumerate.
        std::string prefix;
        if( sample_name_prefix_.option && *sample_name_prefix_.option ) {
            prefix = sample_name_prefix_.value;
        }

        // Go through all filter names, remove the prefix from the filter names (might be empty),
        // and check that the remainder is a number. If so, it's the index (offset by 1).
        for( auto filtered_name : filter_list ) {
            if( ! genesis::utils::starts_with( filtered_name, prefix )) {
                throw CLI::ValidationError(
                    "Invalid sample name used for filtering: \"" + filtered_name  + "\" "
                    "that does not fit the sample prefix given by " +
                    sample_name_prefix_.option->get_name()
                );
            }

            // Here, we know that the filter name contains the prefix, so we can remove it,
            // and convert the rest to a number. We use 1 based indexing in the user provided names,
            // so zero is invalid. But then we need to offset to get our actual zero-based index.
            filtered_name = filtered_name.substr( prefix.length() );
            size_t index = 0;
            try {
                index = genesis::utils::convert_from_string<size_t>( filtered_name, true );
            } catch(...) {
                index = 0;
            }
            if( index == 0 ) {
                throw CLI::ValidationError(
                    "Invalid sample name used for filtering: \"" + filtered_name  + "\" "
                    "that does not contain a valid number after the prefix."
                );
            }
            assert( index > 0 );
            indices.push_back( index - 1 );
        }
    }

    // Return the indices, and whether they are to be interpreted as inverses.
    return { std::move( indices ), is_exclude };
}
