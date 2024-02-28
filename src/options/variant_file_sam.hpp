#ifndef GRENEDALF_OPTIONS_VARIANT_FILE_SAM_H_
#define GRENEDALF_OPTIONS_VARIANT_FILE_SAM_H_

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
    Lucas Czech <lczech@carnegiescience.edu>
    Department of Plant Biology, Carnegie Institution For Science
    260 Panama Street, Stanford, CA 94305, USA
*/

#include "CLI/CLI.hpp"

#include "options/file_input.hpp"
#include "options/variant_file.hpp"
#include "tools/cli_option.hpp"

#include <string>
#include <vector>

// =================================================================================================
//      Variant File Sam Options
// =================================================================================================

/**
 * @brief
 */
class VariantFileSamOptions final : public VariantFileOptions
{
public:

    // -------------------------------------------------------------------------
    //     Constructor and Rule of Five
    // -------------------------------------------------------------------------

    VariantFileSamOptions()  = default;
    ~VariantFileSamOptions() = default;

    VariantFileSamOptions( VariantFileSamOptions const& other ) = default;
    VariantFileSamOptions( VariantFileSamOptions&& )            = default;

    VariantFileSamOptions& operator= ( VariantFileSamOptions const& other ) = default;
    VariantFileSamOptions& operator= ( VariantFileSamOptions&& )            = default;

    // -------------------------------------------------------------------------
    //     Virtual Functions
    // -------------------------------------------------------------------------

private:

    CLI::Option* add_file_input_opt_to_app_(
        CLI::App* sub,
        bool required,
        std::string const& group
    ) override;

    std::string get_default_group_name_() const override
    {
        return "Input SAM/BAM/CRAM";
    }

    VariantInputStream get_stream_(
        std::string const& filename
    ) const override ;

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
