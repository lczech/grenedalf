#ifndef GRENEDALF_OPTIONS_VARIANT_REFERENCE_GENOME_H_
#define GRENEDALF_OPTIONS_VARIANT_REFERENCE_GENOME_H_

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

#include "tools/cli_option.hpp"

#include "genesis/sequence/reference_genome.hpp"
#include "genesis/sequence/sequence_dict.hpp"

#include <string>
#include <vector>

// =================================================================================================
//      Variant Reference Genome Options
// =================================================================================================

/**
 * @brief
 */
class VariantReferenceGenomeOptions
{
public:

    // -------------------------------------------------------------------------
    //     Constructor and Rule of Five
    // -------------------------------------------------------------------------

    VariantReferenceGenomeOptions()  = default;
    ~VariantReferenceGenomeOptions() = default;

    VariantReferenceGenomeOptions( VariantReferenceGenomeOptions const& other ) = default;
    VariantReferenceGenomeOptions( VariantReferenceGenomeOptions&& )            = default;

    VariantReferenceGenomeOptions& operator= ( VariantReferenceGenomeOptions const& other ) = default;
    VariantReferenceGenomeOptions& operator= ( VariantReferenceGenomeOptions&& )            = default;

    // -------------------------------------------------------------------------
    //     Setup Functions
    // -------------------------------------------------------------------------

    void add_reference_genome_opts_to_app(
        CLI::App* sub,
        std::string const& group = "Input Settings"
    );

    // -------------------------------------------------------------------------
    //     Run Functions
    // -------------------------------------------------------------------------

    /**
     * @brief Return the pointer to the reference genome, if provided, or nullptr.
     */
    std::shared_ptr<genesis::sequence::ReferenceGenome> get_reference_genome() const
    {
        prepare_reference_();
        return reference_genome_;
    }

    /**
     * @brief Return the pointer to the reference sequence dictionary, if provided, or nullptr.
     */
    std::shared_ptr<genesis::sequence::SequenceDict> get_reference_dict() const
    {
        prepare_reference_();
        return sequence_dict_;
    }

private:

    void prepare_reference_() const;

    // -------------------------------------------------------------------------
    //     Option Members
    // -------------------------------------------------------------------------

private:

    // Options
    CliOption<std::string> reference_genome_fasta_file_;
    CliOption<std::string> reference_genome_dict_file_;
    CliOption<std::string> reference_genome_fai_file_;

    // Run variables
    mutable std::shared_ptr<genesis::sequence::ReferenceGenome> reference_genome_;
    mutable std::shared_ptr<genesis::sequence::SequenceDict> sequence_dict_;

};

#endif // include guard
