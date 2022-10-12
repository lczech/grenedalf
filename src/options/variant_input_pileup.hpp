#ifndef GRENEDALF_OPTIONS_VARIANT_INPUT_PILEUP_H_
#define GRENEDALF_OPTIONS_VARIANT_INPUT_PILEUP_H_

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
#include "options/variant_input_file.hpp"
#include "tools/cli_option.hpp"

#include <string>
#include <vector>

// =================================================================================================
//      VariantInputPileup Options
// =================================================================================================

/**
 * @brief
 */
class VariantInputPileupOptions final : public VariantInputFileOptions
{
public:

    // -------------------------------------------------------------------------
    //     Constructor and Rule of Five
    // -------------------------------------------------------------------------

    VariantInputPileupOptions()  = default;
    ~VariantInputPileupOptions() = default;

    VariantInputPileupOptions( VariantInputPileupOptions const& other ) = default;
    VariantInputPileupOptions( VariantInputPileupOptions&& )            = default;

    VariantInputPileupOptions& operator= ( VariantInputPileupOptions const& other ) = default;
    VariantInputPileupOptions& operator= ( VariantInputPileupOptions&& )            = default;

    // -------------------------------------------------------------------------
    //     Virtual Functions
    // -------------------------------------------------------------------------

private:

    CLI::Option* add_file_input_opt_to_app_(
        CLI::App* sub,
        bool required,
        std::string const& group
    );

    std::string get_default_group_name_() const
    {
        return "Input (m)pileup";
    }

    bool has_sample_names_() const
    {
        return false;
    }

    VariantInputIterator get_iterator_(
        std::string const& filename,
        VariantInputSampleNamesOptions const& sample_names_options
    ) const;

    // -------------------------------------------------------------------------
    //     Option Members
    // -------------------------------------------------------------------------

private:

    CliOption<size_t>      pileup_min_base_qual_    = 0;
    CliOption<std::string> pileup_quality_encoding_ = "sanger";
    // CliOption<bool> with_quality_string_ = true;
    // CliOption<bool> with_ancestral_base_ = false;


};

#endif // include guard
