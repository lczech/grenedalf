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

#include "options/variant_input_pileup.hpp"

#include "options/global.hpp"
#include "options/variant_input_sample_names.hpp"
#include "tools/misc.hpp"

#include "genesis/population/formats/simple_pileup_input_iterator.hpp"
#include "genesis/population/formats/simple_pileup_reader.hpp"
#include "genesis/sequence/functions/quality.hpp"

#include <algorithm>
#include <cassert>
#include <memory>
#include <stdexcept>

// =================================================================================================
//      Setup Functions
// =================================================================================================

CLI::Option* VariantInputPileupOptions::add_file_input_opt_to_app_(
    CLI::App* sub,
    bool required,
    std::string const& group
) {
    // Correct setup check.
    internal_check(
        file_input_.option() == nullptr,
        "Cannot use the same VariantInputPileupOptions object multiple times."
    );

    // TODO add options for reading: with quality, with ancestral base

    // Add the option
    file_input_.add_multi_file_input_opt_to_app(
        sub, "pileup", "(m)pileup",
        "(plp|mplp|pileup|mpileup)(\\.gz)?",
        "(plp|mplp|pileup|mpileup)[.gz]",
        required, group
    );

    // Min phred score
    pileup_min_base_qual_.option = sub->add_option(
        "--pileup-min-base-qual",
        pileup_min_base_qual_.value,
        "Minimum phred quality score [0-90] for a base in (m)pileup files to be considered. "
        "Bases below this are ignored when computing allele frequencies. "
        "Default is 0, meaning no filtering by phred quality score."
    );
    pileup_min_base_qual_.option->group( group );
    pileup_min_base_qual_.option->check( CLI::Range( static_cast<size_t>(0), static_cast<size_t>(90) ));
    pileup_min_base_qual_.option->needs( file_input_.option() );

    // Quality encoding.
    pileup_quality_encoding_.option = sub->add_option(
        "--pileup-quality-encoding",
        pileup_quality_encoding_.value,
        "Encoding of the quality scores of the bases in (m)pileup files, when using "
        "`--pileup-min-base-qual`. "
        "Default is `\"sanger\"`, which seems to be the most common these days. "
        "Both `\"sanger\"` and `\"illumina-1.8\"` are identical and use an ASCII offset of 33, "
        "while `\"illumina-1.3\"` and `\"illumina-1.5\"` are identical with an ASCII offset of 64 "
        "(we provide different names for completeness). Lastly, `\"solexa\"` has an offset of 64, "
        "but uses a different equation (not phred score) for the encoding."
    );
    pileup_quality_encoding_.option->group( group );
    pileup_quality_encoding_.option->transform(
        CLI::IsMember(
            { "sanger", "illumina-1.3", "illumina-1.5", "illumina-1.8", "solexa" },
            CLI::ignore_case
        )
    );
    pileup_quality_encoding_.option->needs( file_input_.option() );

    return file_input_.option();
}

// =================================================================================================
//      Run Functions
// =================================================================================================

VariantInputPileupOptions::VariantInputIterator VariantInputPileupOptions::get_iterator_(
    std::string const& filename,
    VariantInputSampleNamesOptions const& sample_names_options
) const {
    using namespace genesis::population;

    // We can use the sample filter settings to obtain a list of indices of samples
    // that we want to restrict the reading to. If no filter is given, that list is empty.
    // The second value of the returned pair indicates whether the list in inversed.
    auto const sample_filter = sample_names_options.find_sample_indices_from_sample_filters();

    // Prepare the base Reader with settings as needed.
    auto reader = SimplePileupReader();
    reader.quality_encoding(
        genesis::sequence::guess_quality_encoding_from_name( pileup_quality_encoding_.value )
    );
    reader.min_base_quality( pileup_min_base_qual_.value );

    // Make an iterator.
    auto iterator = make_variant_input_iterator_from_pileup_file(
        filename, sample_filter.first, sample_filter.second, reader
    );

    // Pileup does not have sample names, so set them based on user input or simple enumeration.
    // make_variant_input_iterator_from_pileup_file() returns a list with as many
    // empty strings as the file has samples (exactly one in that case).
    iterator.data().sample_names = sample_names_options.make_anonymous_sample_names(
        iterator.data().sample_names.size()
    );

    return iterator;
}