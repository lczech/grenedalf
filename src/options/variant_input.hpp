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
#include "options/variant_filter_region.hpp"
#include "options/variant_input_file.hpp"
#include "options/variant_input_sample_names.hpp"
#include "tools/cli_option.hpp"

#include "genesis/population/formats/variant_input_iterator.hpp"
#include "genesis/population/genome_locus_set.hpp"
#include "genesis/population/variant.hpp"
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

    using Variant = genesis::population::Variant;
    using GenomeLocusSet = genesis::population::GenomeLocusSet;
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
    //     CLI Setup Functions
    // -------------------------------------------------------------------------

    void add_variant_input_opts_to_app(
        CLI::App* sub,
        bool with_sample_name_opts = true,
        bool with_region_filter_opts = true,
        std::string const& group = "Input Settings"
    );

    // -------------------------------------------------------------------------
    //     Command Setup Functions
    // -------------------------------------------------------------------------

    /**
     * @brief Transformations and filters for individual input sources.
     *
     * These can be added by the commands as needed for their processing, and are executed
     * on each input source. The function has to be called before calling any of the run functions,
     * i.e., sample_names() or get_iterator().
     *
     * We always add the regions filter interally here alreay as the first filter per source,
     * as this is one that is provided in this class here already, and hence offered to all
     * commands that make use of this VariantInputOptions class.
     */
    void add_individual_filter_and_transforms( std::function<bool(Variant&)> const& func ) const;

    /**
     * @brief Transformations and filters for the compbined input VariantInputIterator.
     *
     * These can be added by the commands as needed for their processing, and are executed
     * on the combined input; if only a single file is provided as input by the user,
     * both the individual and these combined functions are applied to that.
     * The function has to be called before calling any of the run functions,
     * i.e., sample_names() or get_iterator().
     */
    void add_combined_filter_and_transforms( std::function<bool(Variant&)> const& func ) const;

    // -------------------------------------------------------------------------
    //     Run Functions
    // -------------------------------------------------------------------------

    /**
     * @brief Get all sample names given in the input file that are not filtered out.
     */
    std::vector<std::string> const& sample_names() const
    {
        prepare_data_();
        return sample_names_;
    }

    /**
     * @brief Get an iterator over the positions in the input file.
     *
     * This takes care of any filtering of samples, chromosomes, and positions.
     */
    VariantInputIterator& get_iterator() const
    {
        prepare_data_();
        return iterator_;
    }

    // -------------------------------------------------------------------------
    //     Internal Helpers
    // -------------------------------------------------------------------------

private:

    void prepare_data_() const;
    void prepare_data_single_file_() const;
    void prepare_data_multiple_files_() const;

    /**
     * @brief Add filters and transformations that are to be applied to each input individually.
     */
    void add_individual_filters_and_transforms_to_iterator_(
        VariantInputIterator& iterator
    ) const;

    /**
     * @brief Add filters and transformations that are to be applied to the combined sample.
     */
    void add_combined_filters_and_transforms_to_iterator_(
        VariantInputIterator& iterator
    ) const;

private:

    // -------------------------------------------------------------------------
    //     Option Members
    // -------------------------------------------------------------------------

    // -------------------------------------
    //     CLI Options
    // -------------------------------------

    // Input file types
    // We use our abstract base class to model different types of input files,
    // to avoid code repetition when adding and processing them here.
    std::vector<std::unique_ptr<VariantInputFileOptions>> input_files_;

    // Sample names, for file types without them, and filters for sample names
    VariantInputSampleNamesOptions input_sample_names_;

    // Filters
    VariantFilterRegionOptions region_filter_;

    // General input settings
    CliOption<std::string> multi_file_loci_set_ = "union";

    // Hidden options to set the LambdaIterator block size for speed.
    CliOption<size_t> iterator_block_size_ = 8192;
    CliOption<size_t> parallel_block_size_ = 4096;

    // Transformations and filters, can be added by the commands as needed for their processing.
    // The individual is applied per input source, and the compbined is applied on the whole
    // VariantInputIterator at the end; if only a single file is provided as input by the user,
    // both are applied to that.
    // We make them mutable, so that these filters/transformations can be added at runtime
    // of the command. Doesn't really matter - we could add them in the setup of the command
    // as well, but somehow, that part of the command seems that it should be more about the CLI.
    mutable std::vector<std::function<bool(Variant&)>> individual_filters_and_transforms_;
    mutable std::vector<std::function<bool(Variant&)>> combined_filters_and_transforms_;

    // -------------------------------------
    //     Run Data
    // -------------------------------------

    // We have different input data formats, but want to convert all of them to our `Variant` type.
    // Working with a statically typed language, this is a bit tricky, so let's introduce a level of
    // abstraction that gives us an iterator over Variants that type-erases the input data format
    // by using a std::function with a lambda that simply returns Variant objects.
    // The `VariantInputIterator` takes care of all of this, and just gives us an iterable,
    // yielding a Variant for each position.
    mutable VariantInputIterator iterator_;

    // Not all formats have sample names, so we need to cache those.
    mutable std::vector<std::string> sample_names_;

};

#endif // include guard
