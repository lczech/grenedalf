#ifndef GRENEDALF_OPTIONS_VARIANT_INPUT_H_
#define GRENEDALF_OPTIONS_VARIANT_INPUT_H_

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

#include "CLI/CLI.hpp"

#include "options/file_input.hpp"
#include "tools/cli_option.hpp"

#include "genesis/population/variant.hpp"
#include "genesis/population/formats/variant_input_iterator.hpp"
#include "genesis/utils/containers/lambda_iterator.hpp"
#include "genesis/utils/containers/range.hpp"

#include <functional>
#include <memory>
#include <string>
#include <utility>
#include <unordered_set>
#include <vector>

// =================================================================================================
//      Variant Input Options
// =================================================================================================

/**
 * @brief
 */
class VariantInputOptions
{
public:

    // -------------------------------------------------------------------------
    //     Typedefs and Enums
    // -------------------------------------------------------------------------

    using VariantInputIterator = genesis::population::VariantInputIterator;

    // -------------------------------------------------------------------------
    //     Constructor and Rule of Five
    // -------------------------------------------------------------------------

    VariantInputOptions()  = default;
    ~VariantInputOptions() = default;

    VariantInputOptions( VariantInputOptions const& other ) = default;
    VariantInputOptions( VariantInputOptions&& )            = default;

    VariantInputOptions& operator= ( VariantInputOptions const& other ) = default;
    VariantInputOptions& operator= ( VariantInputOptions&& )            = default;

    // -------------------------------------------------------------------------
    //     Setup Functions
    // -------------------------------------------------------------------------

    void add_frequency_input_opts_to_app(
        CLI::App* sub,
        // bool required = true,
        // bool with_sample_name_opts = true,
        // bool with_filter_opts = true,
        std::string const& group = "Input"
    );

private:

    CLI::Option* add_sam_input_opt_to_app(
        CLI::App* sub,
        bool required = false,
        std::string const& group = "Input SAM/BAM/CRAM"
    );

    CLI::Option* add_pileup_input_opt_to_app(
        CLI::App* sub,
        bool required = false,
        std::string const& group = "Input (m)pileup"
    );

    CLI::Option* add_sync_input_opt_to_app(
        CLI::App* sub,
        bool required = false,
        std::string const& group = "Input sync"
    );

    CLI::Option* add_vcf_input_opt_to_app(
        CLI::App* sub,
        bool required = false,
        std::string const& group = "Input VCF"
    );

public:

    void add_sample_name_opts_to_app(
        CLI::App* sub,
        std::string const& group = "Sample Names"
    );

    void add_filter_opts_to_app(
        CLI::App* sub,
        std::string const& group = "Filtering"
    );

    // -------------------------------------------------------------------------
    //     Run Functions
    // -------------------------------------------------------------------------

public:

    // -------------------------------------
    //     Getters
    // -------------------------------------

    /**
     * @brief Get all sample names given in the input file that are not filtered out.
     */
    std::vector<std::string> const& sample_names() const;

    // -------------------------------------
    //     Iteration
    // -------------------------------------

    /**
     * @brief Get an iterator over the positions in the input file.
     *
     * This takes care of any filtering of samples, chromosomes, and positions.
     */
    genesis::population::VariantInputIterator& get_iterator() const;

    // -------------------------------------------------------------------------
    //     Internal Helpers
    // -------------------------------------------------------------------------

private:

    void prepare_data_() const;
    void prepare_data_single_file_() const;
    void prepare_data_multiple_files_() const;

    VariantInputIterator prepare_sam_iterator_( std::string const& filename ) const;
    VariantInputIterator prepare_pileup_iterator_( std::string const& filename ) const;
    VariantInputIterator prepare_sync_iterator_( std::string const& filename ) const;
    VariantInputIterator prepare_vcf_iterator_( std::string const& filename ) const;

    void prepare_region_filters_() const;

    void add_region_filters_to_iterator_(
        genesis::population::VariantInputIterator& iterator
    ) const;

    /**
     * @brief Get a list of sample names, for example for pileup or sync files that do not have
     * sample names in the file, or to filter by sample name. The given @p list is interpreted
     * either as a file with one sample name per line, or as a list of sample names, tab or comma
     * separated.
     */
    std::vector<std::string> process_sample_name_list_option_(
        std::string const& list
    ) const;

    /**
     * @brief For file formats that do not have sample names, use this function to get
     * their sample names, taking the sample naming options into account.
     */
    std::vector<std::string> make_anonymous_sample_names_( size_t sample_count ) const;

    /**
     * @brief Return a list of the sample indices, which is needed for some of the readers
     * that only take indices for filtering.
     *
     * Get a list of sample indices that remain after filtering, based on sample names
     * and filtering options, and whether that list is inversed (use all _but_ the given indices).
     *
     * This looks at the filter_samples_include_ and filter_samples_exclude_ list, and determines
     * which one to use as a list of input filters. If neither is given by the user, we return
     * an empty vector. If the exclude list is given, the bool of the returned pair is true,
     * indicating that the indices are to be interpreted as inverses.
     *
     * As this function is only needed for readers that process file types without sample names,
     * we here translate to indices. If the user provided a list of sample names to be used,
     * we use that to find the indices. If not, we use the sample name prefix and simply count
     * sample numbers.
     */
    std::pair<std::vector<size_t>, bool> find_sample_indices_from_sample_filters_() const;

    // -------------------------------------------------------------------------
    //     Option Members
    // -------------------------------------------------------------------------

private:

    // Input file types

    // SAM/BAM/CRAM
    FileInputOptions       sam_file_;
    CliOption<size_t>      sam_min_map_qual_      = 0;
    CliOption<size_t>      sam_min_base_qual_     = 0;
    CliOption<bool>        sam_split_by_rg_       = false;
    CliOption<std::string> sam_flags_include_all_;
    CliOption<std::string> sam_flags_include_any_;
    CliOption<std::string> sam_flags_exclude_all_;
    CliOption<std::string> sam_flags_exclude_any_;

    // Pileup
    FileInputOptions       pileup_file_;
    CliOption<size_t>      pileup_min_base_qual_    = 0;
    CliOption<std::string> pileup_quality_encoding_ = "sanger";
    // CliOption<bool> with_quality_string_ = true;
    // CliOption<bool> with_ancestral_base_ = false;

    // Sync, VCF
    FileInputOptions       sync_file_;
    FileInputOptions       vcf_file_;

    // General input settings
    CliOption<std::string> multi_file_loci_set_ = "union";

    // Sample names
    CliOption<std::string> sample_name_list_ = "";
    CliOption<std::string> sample_name_prefix_ = ""; // "Sample_"

    // Filters for rows and columns
    CliOption<std::vector<std::string>> filter_region_ ;
    CliOption<std::vector<std::string>> filter_region_list_ ;
    CliOption<std::vector<std::string>> filter_region_bed_ ;
    CliOption<std::vector<std::string>> filter_region_gff_ ;
    CliOption<std::vector<std::string>> filter_region_bim_ ;
    CliOption<std::vector<std::string>> filter_region_vcf_ ;
    CliOption<std::string> filter_region_set_ = "union";
    CliOption<std::string> filter_samples_include_ = "";
    CliOption<std::string> filter_samples_exclude_ = "";

    // Hidden options to set the LambdaIterator block size for speed.
    CliOption<size_t> iterator_block_size_ = 4096;
    CliOption<size_t> parallel_block_size_ = 0;

    // We have different input data formats, but want to convert all of them to Variant.
    // Working with a statically typed language, this is a bit tricky, so let's introduce a level of
    // abstraction that gives us an iterator over Variants that type-erases the input data format
    // by using a std::function with a lambda that simply returns Variant objects.
    mutable VariantInputIterator iterator_;

    // Not all formats have sample names, so we need to cache those.
    mutable std::vector<std::string> sample_names_;

    // We keep the region filter here, so that they can be re-used for all inputs.
    // This filter is created by combining (union or intersection) all input filter files.
    mutable std::shared_ptr<genesis::population::GenomeLocusSet> region_filter_;

};

#endif // include guard
