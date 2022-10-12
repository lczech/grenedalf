#ifndef GRENEDALF_OPTIONS_VARIANT_INPUT_SYNC_H_
#define GRENEDALF_OPTIONS_VARIANT_INPUT_SYNC_H_

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
//      VariantInputSync Options
// =================================================================================================

/**
 * @brief
 */
class VariantInputSyncOptions
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

    VariantInputSyncOptions()  = default;
    ~VariantInputSyncOptions() = default;

    VariantInputSyncOptions( VariantInputSyncOptions const& other ) = default;
    VariantInputSyncOptions( VariantInputSyncOptions&& )            = default;

    VariantInputSyncOptions& operator= ( VariantInputSyncOptions const& other ) = default;
    VariantInputSyncOptions& operator= ( VariantInputSyncOptions&& )            = default;

    // -------------------------------------------------------------------------
    //     Setup Functions
    // -------------------------------------------------------------------------

    CLI::Option* add_sync_input_opt_to_app(
        CLI::App* sub,
        bool required = false,
        std::string const& group = "Input sync"
    );

    // -------------------------------------------------------------------------
    //     Run Functions
    // -------------------------------------------------------------------------

    FileInputOptions const& get_file_input_options() const
    {
        return sync_file_;
    }

    VariantInputIterator prepare_sync_iterator(
        std::string const& filename,
        VariantInputSampleNamesOptions const& sample_names_options
    ) const;

    // -------------------------------------------------------------------------
    //     Option Members
    // -------------------------------------------------------------------------

private:

    // Sync
    FileInputOptions sync_file_;

};

#endif // include guard
