#ifndef GRENEDALF_OPTIONS_VARIANT_INPUT_FILE_H_
#define GRENEDALF_OPTIONS_VARIANT_INPUT_FILE_H_

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
#include "tools/cli_option.hpp"

#include "genesis/population/formats/variant_input_iterator.hpp"
#include "genesis/population/genome_locus_set.hpp"
#include "genesis/population/variant.hpp"
#include "genesis/sequence/reference_genome.hpp"

#include <string>
#include <vector>

// Forward Declarations
class VariantInputOptions;
class VariantInputSampleNamesOptions;

// =================================================================================================
//      Variant Input Sam Options
// =================================================================================================

/**
 * @brief
 */
class VariantInputFileOptions
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

    VariantInputFileOptions()  = default;
    virtual ~VariantInputFileOptions() = default;

    VariantInputFileOptions( VariantInputFileOptions const& other ) = default;
    VariantInputFileOptions( VariantInputFileOptions&& )            = default;

    VariantInputFileOptions& operator= ( VariantInputFileOptions const& other ) = default;
    VariantInputFileOptions& operator= ( VariantInputFileOptions&& )            = default;

    // -------------------------------------------------------------------------
    //     Setup Functions
    // -------------------------------------------------------------------------

    CLI::Option* add_file_input_opt_to_app(
        CLI::App* sub,
        bool required = false,
        std::string const& group = ""
    ) {
        return add_file_input_opt_to_app_(
            sub, required, group.empty() ? get_default_group_name_() : group
        );
    }

    // -------------------------------------------------------------------------
    //     Run Functions
    // -------------------------------------------------------------------------

    void add_reference_genome(
        std::shared_ptr<genesis::sequence::ReferenceGenome> reference_genome
    ) const {
        add_reference_genome_( reference_genome );
    }

    /**
     * @brief Get the file input options object that is used to store the input file paths.
     */
    FileInputOptions const& get_file_input_options() const
    {
        return file_input_;
    }

    /**
     * @brief Get the iterator for one of the files, parsing it as the file type that
     * the derived class is for.
     */
    VariantInputIterator get_iterator( std::string const& filename ) const {
        return get_iterator_( filename );
    }

    // -------------------------------------------------------------------------
    //     Non-virtual interface
    // -------------------------------------------------------------------------

protected:

    virtual void add_reference_genome_(
        std::shared_ptr<genesis::sequence::ReferenceGenome> reference_genome
    ) const {
        (void) reference_genome;
    }

    virtual CLI::Option* add_file_input_opt_to_app_(
        CLI::App* sub,
        bool required,
        std::string const& group
    ) = 0;

    virtual std::string get_default_group_name_() const = 0;

    virtual VariantInputIterator get_iterator_( std::string const& filename ) const = 0;

    // -------------------------------------------------------------------------
    //     Option Members
    // -------------------------------------------------------------------------

protected:

    FileInputOptions file_input_;

};

#endif // include guard
