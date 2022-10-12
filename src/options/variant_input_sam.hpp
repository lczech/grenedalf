#ifndef GRENEDALF_OPTIONS_VARIANT_INPUT_SAM_H_
#define GRENEDALF_OPTIONS_VARIANT_INPUT_SAM_H_

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
//      Variant Input Sam Options
// =================================================================================================

/**
 * @brief
 */
class VariantInputSamOptions final : public VariantInputFileOptions
{
public:

    // -------------------------------------------------------------------------
    //     Constructor and Rule of Five
    // -------------------------------------------------------------------------

    VariantInputSamOptions()  = default;
    ~VariantInputSamOptions() = default;

    VariantInputSamOptions( VariantInputSamOptions const& other ) = default;
    VariantInputSamOptions( VariantInputSamOptions&& )            = default;

    VariantInputSamOptions& operator= ( VariantInputSamOptions const& other ) = default;
    VariantInputSamOptions& operator= ( VariantInputSamOptions&& )            = default;

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
        return "Input SAM/BAM/CRAM";
    }

    bool has_sample_names_() const
    {
        return sam_split_by_rg_.value;
    }

    VariantInputIterator get_iterator_(
        std::string const& filename,
        VariantInputSampleNamesOptions const& sample_names_options
    ) const;

    // -------------------------------------------------------------------------
    //     Option Members
    // -------------------------------------------------------------------------

private:

    // SAM/BAM/CRAM
    CliOption<size_t>      sam_min_map_qual_      = 0;
    CliOption<size_t>      sam_min_base_qual_     = 0;
    CliOption<bool>        sam_split_by_rg_       = false;
    CliOption<std::string> sam_flags_include_all_;
    CliOption<std::string> sam_flags_include_any_;
    CliOption<std::string> sam_flags_exclude_all_;
    CliOption<std::string> sam_flags_exclude_any_;


};

#endif // include guard
