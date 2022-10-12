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
//      Variant Input Sam Options
// =================================================================================================

/**
 * @brief
 */
class VariantInputSamOptions
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

    VariantInputSamOptions()  = default;
    ~VariantInputSamOptions() = default;

    VariantInputSamOptions( VariantInputSamOptions const& other ) = default;
    VariantInputSamOptions( VariantInputSamOptions&& )            = default;

    VariantInputSamOptions& operator= ( VariantInputSamOptions const& other ) = default;
    VariantInputSamOptions& operator= ( VariantInputSamOptions&& )            = default;

    // -------------------------------------------------------------------------
    //     Setup Functions
    // -------------------------------------------------------------------------

    CLI::Option* add_sam_input_opt_to_app(
        CLI::App* sub,
        bool required = false,
        std::string const& group = "Input SAM/BAM/CRAM"
    );

    // -------------------------------------------------------------------------
    //     Run Functions
    // -------------------------------------------------------------------------

    FileInputOptions const& get_file_input_options() const
    {
        return sam_file_;
    }

    bool sam_split_by_rg() const
    {
        return sam_split_by_rg_.value;
    }

    VariantInputIterator prepare_sam_iterator(
        std::string const& filename,
        VariantInputSampleNamesOptions const& sample_names_options
    ) const;

    // -------------------------------------------------------------------------
    //     Option Members
    // -------------------------------------------------------------------------

private:

    // SAM/BAM/CRAM
    FileInputOptions       sam_file_;
    CliOption<size_t>      sam_min_map_qual_      = 0;
    CliOption<size_t>      sam_min_base_qual_     = 0;
    CliOption<bool>        sam_split_by_rg_       = false;
    CliOption<std::string> sam_flags_include_all_;
    CliOption<std::string> sam_flags_include_any_;
    CliOption<std::string> sam_flags_exclude_all_;
    CliOption<std::string> sam_flags_exclude_any_;


};

#endif // include guard
