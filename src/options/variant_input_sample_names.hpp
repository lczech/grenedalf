#ifndef GRENEDALF_OPTIONS_VARIANT_INPUT_SAMPLE_NAMES_H_
#define GRENEDALF_OPTIONS_VARIANT_INPUT_SAMPLE_NAMES_H_

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

#include "tools/cli_option.hpp"

#include <string>
#include <vector>

// Forward Declaration
class VariantInputOptions;

// =================================================================================================
//      VariantInputSampleNames Options
// =================================================================================================

/**
 * @brief
 */
class VariantInputSampleNamesOptions
{
public:

    // -------------------------------------------------------------------------
    //     Constructor and Rule of Five
    // -------------------------------------------------------------------------

    VariantInputSampleNamesOptions()  = default;
    ~VariantInputSampleNamesOptions() = default;

    VariantInputSampleNamesOptions( VariantInputSampleNamesOptions const& other ) = default;
    VariantInputSampleNamesOptions( VariantInputSampleNamesOptions&& )            = default;

    VariantInputSampleNamesOptions& operator= ( VariantInputSampleNamesOptions const& other ) = default;
    VariantInputSampleNamesOptions& operator= ( VariantInputSampleNamesOptions&& )            = default;

    // -------------------------------------------------------------------------
    //     Setup Functions
    // -------------------------------------------------------------------------

    void add_sample_name_opts_to_app(
        CLI::App* sub,
        std::string const& group = "Sample Names and Filters"
    );

    // -------------------------------------------------------------------------
    //     Access Functions
    // -------------------------------------------------------------------------

    /**
     * @brief Get the option for setting the sample names.
     */
    CliOption<std::string> const& get_sample_name_list() const
    {
        return sample_name_list_;
    }

    /**
     * @brief Get the option for setting the sample prefix.
     */
    CliOption<std::string> const& get_sample_name_prefix() const
    {
        return sample_name_prefix_;
    }

    /**
     * @brief Get the option for including sample names.
     */
    CliOption<std::string> const& get_filter_samples_include() const
    {
        return filter_samples_include_;
    }

    /**
     * @brief Get the option for excluding sample names.
     */
    CliOption<std::string> const& get_filter_samples_exclude() const
    {
        return filter_samples_exclude_;
    }

    // -------------------------------------------------------------------------
    //     Run Functions
    // -------------------------------------------------------------------------

    /**
     * @brief Get a list of sample names, for example for pileup or sync files that do not have
     * sample names in the file, or to filter by sample name. The given @p list is interpreted
     * either as a file with one sample name per line, or as a list of sample names, tab or comma
     * separated.
     */
    std::vector<std::string> process_sample_name_list_option(
        std::string const& list
    ) const;

    /**
     * @brief For file formats that do not have sample names, use this function to get
     * their sample names, taking the sample naming options into account.
     */
    std::vector<std::string> make_anonymous_sample_names( size_t sample_count ) const;

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
    std::pair<std::vector<size_t>, bool> find_sample_indices_from_sample_filters() const;

    // -------------------------------------------------------------------------
    //     Option Members
    // -------------------------------------------------------------------------

private:

    // Sample names and filters
    CliOption<std::string> sample_name_list_ = "";
    CliOption<std::string> sample_name_prefix_ = ""; // "Sample_"
    CliOption<std::string> filter_samples_include_ = "";
    CliOption<std::string> filter_samples_exclude_ = "";

};

#endif // include guard
