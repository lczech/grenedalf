#ifndef GRENEDALF_OPTIONS_VARIANT_INPUT_FREQUENCY_TABLE_H_
#define GRENEDALF_OPTIONS_VARIANT_INPUT_FREQUENCY_TABLE_H_

/*
    grenedalf - Genome Analyses of Differential Allele Frequencies
    Copyright (C) 2020-2023 Lucas Czech

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
#include "options/variant_input_file.hpp"
#include "tools/cli_option.hpp"

#include <string>
#include <vector>

// =================================================================================================
//      VariantInputFrequencyTable Options
// =================================================================================================

/**
 * @brief
 */
class VariantInputFrequencyTableOptions final : public VariantInputFileOptions
{
public:

    // -------------------------------------------------------------------------
    //     Constructor and Rule of Five
    // -------------------------------------------------------------------------

    /**
     * @brief Create an instance, with additional settings.
     *
     * We cheat here a bit and have the extra settings provided in the constructor,
     * instead of normal getters/setters, as this allows us to set them in the VariantInputOptions
     * when creating the instance, instead of having to cast later etc... hacky, but good enough.
     */
    VariantInputFrequencyTableOptions(
        bool add_extra_opts = true
    )
        : add_extra_opts_( add_extra_opts )
    {}

    ~VariantInputFrequencyTableOptions() = default;

    VariantInputFrequencyTableOptions( VariantInputFrequencyTableOptions const& other ) = default;
    VariantInputFrequencyTableOptions( VariantInputFrequencyTableOptions&& )            = default;

    VariantInputFrequencyTableOptions& operator= (
        VariantInputFrequencyTableOptions const& other
    ) = default;
    VariantInputFrequencyTableOptions& operator= (
        VariantInputFrequencyTableOptions&&
    )= default;

    // -------------------------------------------------------------------------
    //     Helper Functions
    // -------------------------------------------------------------------------

private:

    char get_separator_char_() const;

    // -------------------------------------------------------------------------
    //     Virtual Functions
    // -------------------------------------------------------------------------

private:

    void add_reference_genome_(
        std::shared_ptr<genesis::sequence::ReferenceGenome> reference_genome
    ) const override {
        reference_genome_ = reference_genome;
    }

    CLI::Option* add_file_input_opt_to_app_(
        CLI::App* sub,
        bool required,
        std::string const& group
    ) override;

    void add_extra_file_input_opts_to_app_(
        CLI::App* sub,
        std::string const& group
    );

    std::string get_default_group_name_() const override
    {
        return "Input frequency table";
    }

    VariantInputIterator get_iterator_(
        std::string const& filename
    ) const override;

    // -------------------------------------------------------------------------
    //     Option Members
    // -------------------------------------------------------------------------

private:

    // Basic frequency table options
    CliOption<std::string> separator_char_ = "comma";
    CliOption<double>      int_factor_;
    CliOption<bool>        frequency_is_ref_ = false;

    // Extra options
    bool add_extra_opts_ = true;
    CliOption<std::string> usr_chr_name_;
    CliOption<std::string> usr_pos_name_;
    CliOption<std::string> usr_ref_name_;
    CliOption<std::string> usr_alt_name_;
    CliOption<std::string> usr_smp_ref_name_;
    CliOption<std::string> usr_smp_alt_name_;
    CliOption<std::string> usr_smp_frq_name_;
    CliOption<std::string> usr_smp_cov_name_;

    // Store a pointer to the ref genome, provided by the main VariantInputOptions,
    // so that we can use it to correctly phase input frequencies, e.g., from HAF-pipe.
    mutable std::shared_ptr<genesis::sequence::ReferenceGenome> reference_genome_;

};

#endif // include guard
