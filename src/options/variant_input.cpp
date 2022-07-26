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
#include "genesis/population/formats/genome_region_reader.hpp"
#include "genesis/population/formats/gff_reader.hpp"
#include "genesis/population/formats/map_bim_reader.hpp"
#include "genesis/population/formats/sam_flags.hpp"
#include "genesis/population/formats/sam_variant_input_iterator.hpp"
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
    multi_file_loci_set_.option->group( "Input Settings" );
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
        sub, "sam", "sam/bam/cram", "(sam(\\.gz)?|bam|cram)", ".sam[.gz]|.bam|.cram", required, group
    );

    // Min mapping quality
    sam_min_map_qual_.option = sub->add_option(
        "--sam-min-map-qual",
        sam_min_map_qual_.value,
        "Minimum phred-scaled mapping quality score [0-90] for a read in sam/bam/cram files to be "
        "considered. Any read that is below the given value of mapping quality will be completely "
        "discarded, and its bases not taken into account. "
        "Default is 0, meaning no filtering by base quality."
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
        "Default is 0, meaning no filtering by base quality."
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

    // Flags include all
    sam_flags_include_all_.option = sub->add_flag(
        "--sam-flags-include-all",
        sam_flags_include_all_.value,
        "Only use reads with all bits in the given value present in the FLAG field of the read. "
        "This is equivalent to the `-f` / `--require-flags` setting in `samtools view`. "
        "The value can be specified in hex by beginning with `0x` (i.e., `/^0x[0-9A-F]+/`), "
        "in octal by beginning with `0` (i.e., `/^0[0-7]+/`), as a decimal number not beginning "
        "with '0', or as a comma-, plus-, space-, or vertiacal-bar-separated list of flag names. "
        "We are more lenient in parsing flag names then `samtools`, and allow different "
        "capitalization and delimiteres such as dashes and underscores in the flag names as well."
    );
    sam_flags_include_all_.option->group( group );
    sam_flags_include_all_.option->needs( sam_file_.option() );

    // Flags include any
    sam_flags_include_any_.option = sub->add_flag(
        "--sam-flags-include-any",
        sam_flags_include_any_.value,
        "Only use reads with any bits set in the given value present in the FLAG field of the read. "
        "This is equivalent to the `--rf` / `--incl-flags` / `--include-flags` setting in "
        "`samtools view`. See `--sam-flags-include-all` above for how to specify the value."
    );
    sam_flags_include_any_.option->group( group );
    sam_flags_include_any_.option->needs( sam_file_.option() );

    // Flags exclude all
    sam_flags_exclude_all_.option = sub->add_flag(
        "--sam-flags-exclude-all",
        sam_flags_exclude_all_.value,
        "Do not use reads with all bits set in the given value present in the FLAG field of the read. "
        "This is equivalent to the `-G` setting in `samtools view`. "
        "See `--sam-flags-include-all` above for how to specify the value."
    );
    sam_flags_exclude_all_.option->group( group );
    sam_flags_exclude_all_.option->needs( sam_file_.option() );

    // Flags exclude any
    sam_flags_exclude_any_.option = sub->add_flag(
        "--sam-flags-exclude-any",
        sam_flags_exclude_any_.value,
        "Do not use reads with any bits set in the given value present in the FLAG field of the read. "
        "This is equivalent to the `-F` / `--excl-flags` / `--exclude-flags` setting in "
        "`samtools view`. See `--sam-flags-include-all` above for how to specify the value."
    );
    sam_flags_exclude_any_.option->group( group );
    sam_flags_exclude_any_.option->needs( sam_file_.option() );

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

    // Quality encoding.
    pileup_quality_encoding_.option = sub->add_option(
        "--pileup-quality-encoding",
        pileup_quality_encoding_.value,
        "Encoding of the quality scores of the bases in (m)pileup files, when using "
        "--pileup-min-base-qual. "
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
        sub, "vcf", "vcf/bcf", "(vcf(\\.gz)?|bcf)", ".vcf[.gz]|.bcf", required, group,
        "This expects that the input file has the per-sample VCF FORMAT field `AD` (alleleic depth) "
        "given, containing the counts of the reference and alternative base. "
        "This assumes that the data that was used to create the VCF file was actually a pool of "
        "individuals (e.g., from pool sequencing) for each sample (column) of the VCF file. "
        "We then interpret the `AD` field as the allele counts of each pool of individuals."
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
        "Note that this option can only be used if a single file is given as input. "
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
        "This prefix also works if multiple files are given as input. "
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
        "Genomic region to filter for, in the format \"chr\" (for whole chromosomes), "
        "\"chr:position\", \"chr:start-end\", or \"chr:start..end\". "
        "Positions are 1-based and inclusive (closed intervals). "
        "The option can be provided multiple times, see also `--filter-region-set`."
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
        "The option can be provided multiple times, see also `--filter-region-set`."
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
        "The option can be provided multiple times, see also `--filter-region-set`."
    );
    filter_region_bed_.option->check( CLI::ExistingFile );
    filter_region_bed_.option->group( group );

    // Add option for genomic region filter by GFF/GTF file.
    filter_region_gff_.option = sub->add_option(
        "--filter-region-gff",
        filter_region_gff_.value,
        "Genomic regions to filter for, as a GFF2/GFF3/GTF file. This only uses the chromosome, "
        "as well as start and end information per line, and ignores everything else in the file. "
        "The option can be provided multiple times, see also `--filter-region-set`."
    );
    filter_region_gff_.option->check( CLI::ExistingFile );
    filter_region_gff_.option->group( group );

    // Add option for genomic region filter by PLINK MAP/BIM file.
    filter_region_bim_.option = sub->add_option(
        "--filter-region-map-bim",
        filter_region_bim_.value,
        "Genomic regions to filter for, as a MAP or BIM file as used in PLINK. This only "
        "uses the chromosome and coordinate per line, and ignores everything else in the file. "
        "The option can be provided multiple times, see also `--filter-region-set`."
    );
    filter_region_bim_.option->check( CLI::ExistingFile );
    filter_region_bim_.option->group( group );

    // Add option for genomic region filter by VCF file.
    filter_region_vcf_.option = sub->add_option(
        "--filter-region-vcf",
        filter_region_vcf_.value,
        "Genomic regions to filter for, as a VCF/BCF file (such as a known-variants file). This "
        "only uses the chromosome and position per line, and ignores everything else in the file. "
        "The option can be provided multiple times, see also `--filter-region-set`."
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

    // We first read all region filters, so that they can be re-used for all inputs.
    prepare_region_filters_();

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
    add_individual_filters_and_transforms_to_iterator_( iterator_ );
    add_combined_filters_and_transforms_to_iterator_( iterator_ );

    // Set the buffer block size of the iterator, for multi-threaded processing speed.
    // Needs to be tested - might give more or less advantage depending on setting.
    iterator_.block_size( iterator_block_size_.value );

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

    // Get whether the user wants the union or intersection of all parallel input loci.
    auto const contribution = get_enum_map_value(
        multi_file_contribution_type_map_, multi_file_loci_set_.value
    );

    // Make a parallel input iterator, and add all input from all file formats to it,
    // using the same contribution type for all of them, which either results in the union
    // or the intersection of all input loci. See VariantParallelInputIterator for details.
    VariantParallelInputIterator parallel_it;
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
    // and add the individual samples filters to all of them, e.g., so that regions are filtered
    // before the data even reaches the parallel iterator. Faster!
    assert( sample_names_.empty() );
    for( auto& input : parallel_it.inputs() ) {
        for( auto const& sample_name : input.data().sample_names ) {
            sample_names_.push_back( input.data().source_name + ":" + sample_name );
        }
        add_individual_filters_and_transforms_to_iterator_( input );

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

    // Set the buffer block size of the iterator, for multi-threaded speed.
    // We only buffer the final parallel iterator, and not its individual sources,
    // in order to not keep too many threads from spawning... Might need testing and refinement.
    iterator_.block_size( iterator_block_size_.value );
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
    reader.flags_include_all( string_to_sam_flag( sam_flags_include_all_.value ));
    reader.flags_include_any( string_to_sam_flag( sam_flags_include_any_.value ));
    reader.flags_exclude_all( string_to_sam_flag( sam_flags_exclude_all_.value ));
    reader.flags_exclude_any( string_to_sam_flag( sam_flags_exclude_any_.value ));

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
//     Filters and Transformations
// =================================================================================================

// -------------------------------------------------------------------------
//     prepare_region_filters_
// -------------------------------------------------------------------------

void VariantInputOptions::prepare_region_filters_() const
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
        assert( chr_cnt = full_chr + spec_chr );

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
    if( region_filter_ ) {
        iterator.add_filter( filter_by_region( region_filter_ ));
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

    // Most of our input sources do not provide ref, and almost non provide alt bases.
    // So we use our guess function to augment the data. The function is idempotent
    // (unless we set the `force` parameter, which we do not do here), so for sources that do
    // contain ref and/or alt bases, nothing changes.
    iterator.add_transform( []( Variant& var ){
        guess_and_set_ref_and_alt_bases( var );
    });
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
                std::to_string( sample_count ) +
                " samples (after filtering, if a sample name filter was given)."
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
