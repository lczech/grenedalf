#ifndef GRENEDALF_OPTIONS_VARIANT_FILTER_NUMERICAL_H_
#define GRENEDALF_OPTIONS_VARIANT_FILTER_NUMERICAL_H_

/*
    grenedalf - Genome Analyses of Differential Allele Frequencies
    Copyright (C) 2020-2024 Lucas Czech

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
    Lucas Czech <lucas.czech@sund.ku.dk>
    University of Copenhagen, Globe Institute, Section for GeoGenetics
    Oster Voldgade 5-7, 1350 Copenhagen K, Denmark
*/

#include "CLI/CLI.hpp"

#include "tools/cli_option.hpp"

#include "genesis/population/filter/sample_counts_filter.hpp"
#include "genesis/population/filter/sample_counts_filter_numerical.hpp"
#include "genesis/population/filter/variant_filter.hpp"
#include "genesis/population/filter/variant_filter_numerical.hpp"
#include "genesis/population/function/functions.hpp"

#include <functional>
#include <string>
#include <utility>
#include <vector>

// =================================================================================================
//      Variant Filter Numerical Options
// =================================================================================================

/**
 * @brief
 */
class VariantFilterNumericalOptions
{
public:

    // -------------------------------------------------------------------------
    //     Constructor and Rule of Five
    // -------------------------------------------------------------------------

    VariantFilterNumericalOptions()  = default;
    ~VariantFilterNumericalOptions() = default;

    VariantFilterNumericalOptions( VariantFilterNumericalOptions const& other ) = default;
    VariantFilterNumericalOptions( VariantFilterNumericalOptions&& )            = default;

    VariantFilterNumericalOptions& operator= ( VariantFilterNumericalOptions const& other ) = default;
    VariantFilterNumericalOptions& operator= ( VariantFilterNumericalOptions&& )            = default;

    // -------------------------------------------------------------------------
    //     Setup Functions
    // -------------------------------------------------------------------------

    void add_numerical_filter_opts_to_app(
        CLI::App* sub,
        std::string const& group = "Numerical Filters"
    );

    // --------------------------------------------
    //     Sample Filters
    // --------------------------------------------

    void add_sample_filter_opts_to_app(
        CLI::App* sub,
        std::string const& group = "Numerical Filters"
    );

    void add_sample_count_filter_opts_to_app(
        CLI::App* sub,
        bool add_sample_min_count = true,
        bool add_sample_max_count = true,
        std::string const& group = "Numerical Filters"
    );

    void add_sample_read_depth_filter_opts_to_app(
        CLI::App* sub,
        bool add_sample_min_read_depth = true,
        bool add_sample_max_read_depth = true,
        std::string const& group = "Numerical Filters"
    );

    void add_sample_snp_filter_opts_to_app(
        CLI::App* sub,
        bool add_sample_only_snps = true,
        bool add_sample_only_biallelic_snps = true,
        std::string const& group = "Numerical Filters"
    );

    // --------------------------------------------
    //     Total Filters
    // --------------------------------------------

    void add_total_filter_opts_to_app(
        CLI::App* sub,
        std::string const& group = "Numerical Filters"
    );

    void add_total_read_depth_filter_opts_to_app(
        CLI::App* sub,
        bool add_total_min_read_depth = true,
        bool add_total_max_read_depth = true,
        std::string const& group = "Numerical Filters"
    );

    void add_total_snp_filter_opts_to_app(
        CLI::App* sub,
        bool add_total_only_snps = true,
        bool add_total_only_biallelic_snps = true,
        std::string const& group = "Numerical Filters"
    );

    void add_total_snp_count_opts_to_app(
        CLI::App* sub,
        bool add_total_min_count = true,
        bool add_total_max_count = true,
        std::string const& group = "Numerical Filters"
    );

    void add_total_freq_filter_opts_to_app(
        CLI::App* sub,
        std::string const& group = "Numerical Filters"
    );

    // -------------------------------------------------------------------------
    //     Run Functions
    // -------------------------------------------------------------------------

    /**
     * @brief Get the samples filter, along with a bool indicatin whether any options
     * was provided/changed by the user.
     */
    std::pair<genesis::population::SampleCountsFilterNumericalParams, bool> get_sample_filter_params(
        genesis::population::SampleCountsFilterNumericalParams filter_params =
        genesis::population::SampleCountsFilterNumericalParams{}
    ) const;

    /**
     * @brief Make the sample filter.
     *
     * If the user did not provide any option that changed the filter from the default values,
     * an empty filter function is returned, so that we do not need to add it to the input filters.
     *
     * The return value of the filter function is `true` if any of the samples pass the filters.
     * Only if all of them fail can we savely remove this position from the input.
     */
    std::function<bool( genesis::population::Variant& )> make_sample_filter() const;

    /**
     * @brief Overload that takes default filter values.
     *
     * The values provided by the user will overwrite the ones in the filter.
     * This always returns a valid filter function, as we assume that if default filter settings
     * are provided, those shall be used.
     */
    std::function<bool( genesis::population::Variant& )> make_sample_filter(
        genesis::population::SampleCountsFilterNumericalParams filter_params
    ) const;

    /**
     * @brief Get the total filter, along with a bool indicatin whether any options
     * was provided/changed by the user.
     */
    std::pair<genesis::population::VariantFilterNumericalParams, bool> get_total_filter_params(
        genesis::population::VariantFilterNumericalParams filter_params =
        genesis::population::VariantFilterNumericalParams{}
    ) const;

    /**
     * @brief Make the total filter.
     *
     * If the user did not provide any option that changed the filter from the default values,
     * an empty filter function is returned, so that we do not need to add it to the input filters.
     */
    std::function<bool( genesis::population::Variant&)> make_total_filter() const;

    /**
     * @brief Overload that takes default filter values.
     *
     * The values provided by the user will overwrite the ones in the filter.
     * This always returns a valid filter function, as we assume that if default filter settings
     * are provided, those shall be used.
     */
    std::function<bool( genesis::population::Variant&)> make_total_filter(
        genesis::population::VariantFilterNumericalParams filter_params
    ) const;

    /**
     * @brief Print a report of the statistics of filtering.
     */
    void print_report() const;

    // -------------------------------------------------------------------------
    //     Option Members
    // -------------------------------------------------------------------------

    // Public, so that the defaults can be changed by the commands.

    // Filters for BaseCounts, modelled after SampleCountsFilterNumericalParams
    CliOption<size_t> sample_min_count           = 0;
    CliOption<size_t> sample_max_count           = 0;
    CliOption<size_t> sample_del_count           = 0;
    CliOption<size_t> sample_min_read_depth      = 0;
    CliOption<size_t> sample_max_read_depth      = 0;
    CliOption<bool>   sample_only_snps           = false;
    CliOption<bool>   sample_only_biallelic_snps = false;

    // Filters for Variant, modelled after VariantFilter
    CliOption<size_t> total_min_read_depth       = 0;
    CliOption<size_t> total_max_read_depth       = 0;
    CliOption<size_t> total_del_count            = 0;
    CliOption<bool>   total_only_snps            = false;
    CliOption<bool>   total_only_biallelic_snps  = false;
    CliOption<size_t> total_min_count            = 0;
    CliOption<size_t> total_max_count            = 0;
    CliOption<double> total_min_frequency        = 0.0;

private:

    // Run variables.
    mutable genesis::population::SampleCountsFilterStats sample_stats_;
    mutable genesis::population::VariantFilterStats      total_stats_;

};

#endif // include guard
