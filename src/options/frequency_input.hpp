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

#include "genesis/population/formats/simple_pileup_input_iterator.hpp"
#include "genesis/population/formats/simple_pileup_reader.hpp"
#include "genesis/population/formats/vcf_input_iterator.hpp"
#include "genesis/population/window/sliding_window_iterator.hpp"
#include "genesis/population/window/window.hpp"
#include "genesis/population/variant.hpp"
#include "genesis/utils/containers/lambda_iterator.hpp"
#include "genesis/utils/containers/range.hpp"

#include <functional>
#include <string>
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
        std::string const& group = "Input"
    );

private:

    CLI::Option* add_pileup_input_opt_to_app(
        CLI::App* sub,
        bool required = true,
        std::string const& group = "Input"
    );

    CLI::Option* add_vcf_input_opt_to_app(
        CLI::App* sub,
        bool required = true,
        std::string const& group = "Input"
    );

    // -------------------------------------------------------------------------
    //     Run Functions
    // -------------------------------------------------------------------------

public:

    std::vector<std::string> const& sample_names() const;

    genesis::utils::Range<genesis::utils::LambdaIterator<genesis::population::Variant>>
    get_iterator() const;

    genesis::population::SlidingWindowIterator<
        genesis::utils::LambdaIterator<genesis::population::Variant>,
        genesis::population::Variant
    >
    get_sliding_window_iterator() const;

    // -------------------------------------------------------------------------
    //     Internal Helpers
    // -------------------------------------------------------------------------

    void prepare_data_() const;
    void prepare_data_pileup_() const;
    void prepare_data_vcf_() const;

    // -------------------------------------------------------------------------
    //     Option Members
    // -------------------------------------------------------------------------

private:

    CliOption<std::string> pileup_file_ = "";
    CliOption<std::string> pileup_sample_prefix_ = "Sample_";
    CliOption<std::string> vcf_file_ = "";
    CliOption<std::string> region_ = "";

    CliOption<size_t> window_width_  = 1000;
    CliOption<size_t> window_stride_ = 1000;

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
