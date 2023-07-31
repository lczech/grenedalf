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

#include "options/variant_reference_genome.hpp"

#include "options/global.hpp"
#include "tools/misc.hpp"

#include "genesis/population/functions/variant_input_iterator.hpp"
#include "genesis/sequence/formats/fasta_reader.hpp"
#include "genesis/sequence/functions/dict.hpp"
#include "genesis/utils/core/fs.hpp"
#include "genesis/utils/text/convert.hpp"
#include "genesis/utils/text/string.hpp"

#include <cassert>
#include <stdexcept>
#include <string>
#include <vector>

// =================================================================================================
//      Setup Functions
// =================================================================================================

void VariantReferenceGenomeOptions::add_reference_genome_opts_to_app(
    CLI::App* sub,
    std::string const& group
) {
    // Correct setup check.
    internal_check(
        reference_genome_fasta_file_.option == nullptr,
        "Cannot use the same VariantReferenceGenomeOptions object multiple times."
    );

    // Add option for reading the reference genome.
    reference_genome_fasta_file_.option = sub->add_option(
        "--reference-genome-fasta-file",
        reference_genome_fasta_file_.value,
        "Provide a reference genome in `.fasta[.gz]` format. This allows to correctly assign the "
        "reference bases in file formats that do not store them, and serves as an integrity check "
        "in those that do. It further is used as a sequence dictionary to determine the chromosome "
        "order and length, on behalf of a dict or fai file."
    );
    reference_genome_fasta_file_.option->group( group );
    reference_genome_fasta_file_.option->check( CLI::ExistingFile );

    // Add option for reading a sequence dict.
    reference_genome_dict_file_.option = sub->add_option(
        "--reference-genome-dict-file",
        reference_genome_dict_file_.value,
        "Provide a reference genome sequence dictionary in `.dict` format. It is used to determine "
        "the chromosome order and length, without having to provide the full reference genome."
    );
    reference_genome_dict_file_.option->group( group );
    reference_genome_dict_file_.option->check( CLI::ExistingFile );

    // Add option for reading a sequence fai.
    reference_genome_fai_file_.option = sub->add_option(
        "--reference-genome-fai-file",
        reference_genome_fai_file_.value,
        "Provide a reference genome sequence dictionary in `.fai` format. It is used to determine "
        "the chromosome order and length, without having to provide the full reference genome."
    );
    reference_genome_fai_file_.option->group( group );
    reference_genome_fai_file_.option->check( CLI::ExistingFile );

    // Mutally exclusive, to keep it simple. Otherwise, we'd have to check
    // that those files agree with each other. Might add later.
    reference_genome_fasta_file_.option->excludes( reference_genome_dict_file_.option );
    reference_genome_fasta_file_.option->excludes( reference_genome_fai_file_.option );
    reference_genome_dict_file_.option->excludes( reference_genome_fasta_file_.option );
    reference_genome_dict_file_.option->excludes( reference_genome_fai_file_.option );
    reference_genome_fai_file_.option->excludes( reference_genome_fasta_file_.option );
    reference_genome_fai_file_.option->excludes( reference_genome_dict_file_.option );
}

// =================================================================================================
//      Run Functions
// =================================================================================================

void VariantReferenceGenomeOptions::prepare_reference_() const
{
    using namespace genesis::sequence;
    using namespace genesis::utils;

    // Check if we already ran this function (which can happen depending on the order in which
    // we downstream request this). As we always set the sequence_dict_ if there is any
    // reference-related input option, we just use this as a check whether we already ran.
    if( sequence_dict_ ) {
        return;
    }

    // Internally check that at most one was set.
    size_t used_ref_opts = 0;

    // Read the reference genome first, if provided, as some formats might want to use it.
    if( *reference_genome_fasta_file_.option ) {
        ++used_ref_opts;

        LOG_MSG << "Reading reference genome fasta";
        auto reader = FastaReader();
        reference_genome_ = std::make_shared<ReferenceGenome>(
            reader.read_reference_genome( from_file( reference_genome_fasta_file_.value ))
        );

        // Some user output. Nope, we print the sequence dict instead, for consistency.
        // LOG_MSG1 << "Reference genome contains " << reference_genome_->size() << " chromosome"
        //          << ( reference_genome_->size() != 1 ? "s" : "" );
        // for( auto const& chr : *reference_genome_ ) {
        //     LOG_MSG2 << " - " << chr.label();
        // }

        // Also use the ref genome to build a dict for the parallel iterator and the seq check.
        sequence_dict_ = std::make_shared<SequenceDict>(
            reference_genome_to_dict( *reference_genome_ )
        );
    }

    // Read the dict or fai as an alternative.
    if( *reference_genome_dict_file_.option ) {
        ++used_ref_opts;

        LOG_MSG << "Reading reference genome dict";
        sequence_dict_ = std::make_shared<SequenceDict>(
            read_sequence_dict( from_file( reference_genome_dict_file_.value ))
        );
    }
    if( *reference_genome_fai_file_.option ) {
        ++used_ref_opts;

        LOG_MSG << "Reading reference genome fai";
        sequence_dict_ = std::make_shared<SequenceDict>(
            read_sequence_fai( from_file( reference_genome_fai_file_.value ))
        );
    }

    // User output of what we just read.
    if( sequence_dict_ ) {
        LOG_MSG1 << "Reference genome contains " << sequence_dict_->size() << " chromosome"
                 << ( sequence_dict_->size() != 1 ? "s" : "" );
        for( auto const& chr : *sequence_dict_ ) {
            LOG_MSG2 << " - " << chr.name;
        }
    }

    // Final check.
    internal_check(
        used_ref_opts < 2,
        "More than one option related to reference genome and sequence dictionary is set."
    );
}
