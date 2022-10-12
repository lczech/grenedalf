#ifndef GRENEDALF_OPTIONS_VARIANT_INPUT_VCF_H_
#define GRENEDALF_OPTIONS_VARIANT_INPUT_VCF_H_

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
#include "tools/cli_option.hpp"

#include "genesis/population/formats/variant_input_iterator.hpp"
#include "genesis/population/genome_locus_set.hpp"
#include "genesis/population/variant.hpp"

#include <string>
#include <vector>

// Forward Declaration
class VariantInputOptions;
class VariantInputSampleNamesOptions;

// =================================================================================================
//      VariantInputVcf Options
// =================================================================================================

/**
 * @brief
 */
class VariantInputVcfOptions
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

    VariantInputVcfOptions()  = default;
    ~VariantInputVcfOptions() = default;

    VariantInputVcfOptions( VariantInputVcfOptions const& other ) = default;
    VariantInputVcfOptions( VariantInputVcfOptions&& )            = default;

    VariantInputVcfOptions& operator= ( VariantInputVcfOptions const& other ) = default;
    VariantInputVcfOptions& operator= ( VariantInputVcfOptions&& )            = default;

    // -------------------------------------------------------------------------
    //     Setup Functions
    // -------------------------------------------------------------------------

    CLI::Option* add_vcf_input_opt_to_app(
        CLI::App* sub,
        bool required = false,
        std::string const& group = "Input VCF/BCF"
    );

    // -------------------------------------------------------------------------
    //     Run Functions
    // -------------------------------------------------------------------------

    FileInputOptions const& get_file_input_options() const
    {
        return vcf_file_;
    }

    VariantInputIterator prepare_vcf_iterator(
        std::string const& filename,
        VariantInputSampleNamesOptions const& sample_names_options
    ) const;

    // -------------------------------------------------------------------------
    //     Option Members
    // -------------------------------------------------------------------------

private:

    // Vcf
    FileInputOptions vcf_file_;

};

#endif // include guard
