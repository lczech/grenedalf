#ifndef GRENEDALF_OPTIONS_VARIANT_FILE_VCF_H_
#define GRENEDALF_OPTIONS_VARIANT_FILE_VCF_H_

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
#include "options/variant_file.hpp"
#include "tools/cli_option.hpp"

#include <string>
#include <vector>

// =================================================================================================
//      Variant File Vcf Options
// =================================================================================================

/**
 * @brief
 */
class VariantFileVcfOptions final : public VariantFileOptions
{
public:

    // -------------------------------------------------------------------------
    //     Constructor and Rule of Five
    // -------------------------------------------------------------------------

    VariantFileVcfOptions()  = default;
    ~VariantFileVcfOptions() = default;

    VariantFileVcfOptions( VariantFileVcfOptions const& other ) = default;
    VariantFileVcfOptions( VariantFileVcfOptions&& )            = default;

    VariantFileVcfOptions& operator= ( VariantFileVcfOptions const& other ) = default;
    VariantFileVcfOptions& operator= ( VariantFileVcfOptions&& )            = default;

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
        return "Input VCF/BCF";
    }

    VariantInputIterator get_iterator_(
        std::string const& filename
    ) const override;

    // -------------------------------------------------------------------------
    //     Option Members
    // -------------------------------------------------------------------------

private:

    // Vcf
    // no extra options for the user here at the moment

};

#endif // include guard
