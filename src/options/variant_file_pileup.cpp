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

#include "options/variant_file_pileup.hpp"

#include "options/global.hpp"
#include "tools/misc.hpp"

#include "genesis/population/formats/simple_pileup_common.hpp"
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

CLI::Option* VariantFilePileupOptions::add_file_input_opt_to_app_(
    CLI::App* sub,
    bool required,
    std::string const& group
) {
    // Correct setup check.
    internal_check(
        file_input_.option() == nullptr,
        "Cannot use the same VariantFilePileupOptions object multiple times."
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

VariantFilePileupOptions::VariantInputIterator VariantFilePileupOptions::get_iterator_(
    std::string const& filename
) const {
    using namespace genesis::population;
    using namespace genesis::sequence;

    // Get the encoding set by the user.
    auto const user_enc = guess_quality_encoding_from_name(
        pileup_quality_encoding_.value
    );

    // Check that the quality encoding is correct.
    // We run in a try block, as this might fail if the first lines do not contain any data...
    try {
        size_t const max_guess_lines = 100;
        auto const guess_enc = guess_pileup_quality_encoding(
            genesis::utils::from_file( filename ), max_guess_lines
        );
        if( ! compatible_quality_encodings( guess_enc, user_enc )) {
            LOG_WARN << pileup_quality_encoding_.option->get_name() << " set to "
                     << quality_encoding_name( user_enc, true ) << ", but input file \""
                     << filename << "\" seems to use " << quality_encoding_name( guess_enc, true )
                     << " instead, based on the quality codes found in the first " << max_guess_lines
                     << " lines of the file. Will continue now, but you should verify that the "
                     << "correct quality encoding is being used.";
        }
    } catch(...) {}

    // Prepare the base Reader with settings as needed.
    auto reader = SimplePileupReader();
    reader.quality_encoding( user_enc );
    reader.min_base_quality( pileup_min_base_qual_.value );

    // Make an iterator.
    return make_variant_input_iterator_from_pileup_file( filename, reader );
}
