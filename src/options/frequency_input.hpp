#ifndef GRENEDALF_OPTIONS_FREQUENCY_INPUT_H_
#define GRENEDALF_OPTIONS_FREQUENCY_INPUT_H_

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

#include "CLI/CLI.hpp"

#include "tools/cli_option.hpp"

#include "genesis/population/variant.hpp"
#include "genesis/population/window/sliding_window_iterator.hpp"
#include "genesis/population/window/window.hpp"
#include "genesis/utils/containers/lambda_iterator.hpp"
#include "genesis/utils/containers/range.hpp"

#include <functional>
#include <string>
#include <utility>
#include <vector>

// =================================================================================================
//      Frequency Input Options
// =================================================================================================

/**
 * @brief
 */
class FrequencyInputOptions
{
public:

    // -------------------------------------------------------------------------
    //     Constructor and Rule of Five
    // -------------------------------------------------------------------------

    FrequencyInputOptions()  = default;
    virtual ~FrequencyInputOptions() = default;

    FrequencyInputOptions( FrequencyInputOptions const& other ) = default;
    FrequencyInputOptions( FrequencyInputOptions&& )            = default;

    FrequencyInputOptions& operator= ( FrequencyInputOptions const& other ) = default;
    FrequencyInputOptions& operator= ( FrequencyInputOptions&& )            = default;

    // -------------------------------------------------------------------------
    //     Setup Functions
    // -------------------------------------------------------------------------

    /**
     * @brief
     */
    void add_frequency_input_opts_to_app(
        CLI::App* sub,
        // bool required = true,
        std::string const& group = "Input",
        bool with_filters = true
    );

private:

    CLI::Option* add_pileup_input_opt_to_app(
        CLI::App* sub,
        bool required = true,
        std::string const& group = "Input"
    );

    CLI::Option* add_sync_input_opt_to_app(
        CLI::App* sub,
        bool required = true,
        std::string const& group = "Input"
    );

    CLI::Option* add_vcf_input_opt_to_app(
        CLI::App* sub,
        bool required = true,
        std::string const& group = "Input"
    );

    CLI::Option* add_sample_name_prefix_opt_to_app(
        CLI::App* sub,
        std::string const& group = "Input"
    );

    void add_filter_opts_to_app(
        CLI::App* sub,
        std::string const& group = "Input"
    );

public:

    void add_sliding_window_opts_to_app(
        CLI::App* sub,
        std::string const& group = "Sliding Window"
    );

    // -------------------------------------------------------------------------
    //     Run Functions
    // -------------------------------------------------------------------------

public:

    /**
     * @brief Get all sample names given in the input file that are not filtered out.
     */
    std::vector<std::string> const& sample_names() const;

    /**
     * @brief Get an iterator over the positions in the input file.
     *
     * This takes care of any filtering of samples, chromosomes, and positions.
     */
    genesis::utils::Range<genesis::utils::LambdaIterator<genesis::population::Variant>>
    get_iterator() const;

    /**
     * @brief Get a sliding window iterator that dereferences to a vector of BaseCounts.
     *
     * This is slightly faster than iterating over Variant%s, because we can save the extra
     * effor of copying chromosome names for each entry in the window.
     */
    genesis::population::SlidingWindowIterator<
        genesis::utils::LambdaIterator<genesis::population::Variant>,
        genesis::population::Variant,
        std::vector<genesis::population::BaseCounts>
    >
    get_base_count_sliding_window_iterator() const;

    /**
     * @brief Get a sliding window iterator that dereferences to Variant%s.
     *
     * This is useful if the per-position reference and alternative base are also needed.
     * We could also refactor a bit and include those in a new type of class if needed,
     * but for now, this approach works well enough.
     */
    genesis::population::SlidingWindowIterator<
        genesis::utils::LambdaIterator<genesis::population::Variant>,
        genesis::population::Variant
    >
    get_variant_sliding_window_iterator() const;

    // -------------------------------------------------------------------------
    //     Internal Helpers
    // -------------------------------------------------------------------------

private:

    void prepare_data_() const;
    void prepare_data_pileup_() const;
    void prepare_data_sync_() const;
    void prepare_data_vcf_() const;

    /**
     * @brief Get the list of sample names to filter, interpreting the given @filter_samples_value
     * either as a file or as a list of sample names.
     */
    std::vector<std::string> get_sample_name_list( std::string const& filter_samples_value ) const;

    /**
     * @brief Use the sample name filter options to get a vector of bools determining which samples
     * we want.
     *
     * This list of bools is needed by the pileup and sync readers for example, as these file
     * formats do not have sample names that we can use for filtering.
     */
    std::vector<bool> get_sample_filter( std::vector<std::string> const& sample_names ) const;

    /**
     * @brief Get a list of sample indices that remain after filtering, that is, the indices of
     * the positions where the provided vector is true.
     */
    std::vector<size_t> get_sample_filter_indices( std::vector<bool> const& sample_filter ) const;

    // -------------------------------------------------------------------------
    //     Option Members
    // -------------------------------------------------------------------------

private:

    // Input file types
    CliOption<std::string> pileup_file_ = "";
    CliOption<std::string> sync_file_   = "";
    CliOption<std::string> vcf_file_    = "";
    CliOption<std::string> sample_name_prefix_ = ""; // "Sample_"

    // Filters for rows and columns
    CliOption<std::string> filter_region_ = "";
    CliOption<std::string> filter_samples_include_ = "";
    CliOption<std::string> filter_samples_exclude_ = "";

    // Window settings
    CliOption<size_t> window_width_  = 1000;
    CliOption<size_t> window_stride_ = 0;

    // We have different input data formats, but want to convert all of them to Variant.
    // This is a bit tricky, as we are working with templates for things such as SlidingWindowIterator,
    // but don't want all our grenedalf code to be templated. So, let's introduce a level of
    // abstraction that gives us an iterator over Variants that type erases the input data format
    // by using a std::function with a lambda that simply returns Variant objects.
    mutable genesis::utils::LambdaIteratorGenerator<genesis::population::Variant> generator_;

    // Not all formats have sample names, so we need to cache those.
    mutable std::vector<std::string> sample_names_;

};

#endif // include guard
