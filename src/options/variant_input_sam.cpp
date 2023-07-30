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

#include "options/variant_input_sam.hpp"

#include "options/global.hpp"
#include "tools/misc.hpp"

#include "genesis/population/formats/sam_flags.hpp"
#include "genesis/population/formats/sam_variant_input_iterator.hpp"

#include <cassert>
#include <unordered_set>
#include <stdexcept>

// =================================================================================================
//      Setup Functions
// =================================================================================================

CLI::Option* VariantInputSamOptions::add_file_input_opt_to_app_(
    CLI::App* sub,
    bool required,
    std::string const& group
) {
    // Correct setup check.
    internal_check(
        file_input_.option() == nullptr,
        "Cannot use the same VariantInputSamOptions object multiple times."
    );

    // TODO add min_depth max_depth and max_accumulation_depth, and add FLAG to reader and here.
    // but at least the first two should actually be general filter settings.

    // Add the option
    file_input_.add_multi_file_input_opt_to_app(
        sub, "sam", "sam/bam/cram", "(sam(\\.gz)?|bam|cram)", ".sam[.gz]|.bam|.cram", required, group
    );

    // Min mapping quality
    sam_min_map_qual_.option = sub->add_option(
        "--sam-min-map-qual",
        sam_min_map_qual_.value,
        "Minimum phred-scaled mapping quality score [0-90] for a read in sam/bam/cram files to be "
        "considered. Any read that is below the given value of mapping quality will be completely "
        "discarded, and its bases not taken into account. "
        "Default is 0, meaning no filtering by base quality."
    );
    sam_min_map_qual_.option->group( group );
    sam_min_map_qual_.option->check( CLI::Range( static_cast<size_t>(0), static_cast<size_t>(90) ));
    sam_min_map_qual_.option->needs( file_input_.option() );

    // Min base qual
    sam_min_base_qual_.option = sub->add_option(
        "--sam-min-base-qual",
        sam_min_base_qual_.value,
        "Minimum phred-scaled quality score [0-90] for a base in sam/bam/cram files to be "
        "considered. Bases below this are ignored when computing allele frequencies. "
        "Default is 0, meaning no filtering by base quality."
    );
    sam_min_base_qual_.option->group( group );
    sam_min_base_qual_.option->check( CLI::Range( static_cast<size_t>(0), static_cast<size_t>(90) ));
    sam_min_base_qual_.option->needs( file_input_.option() );

    // Split by RG read group tag.
    sam_split_by_rg_.option = sub->add_flag(
        "--sam-split-by-rg",
        sam_split_by_rg_.value,
        "Instead of considering the whole sam/bam/cram file as one large colletion of reads, "
        "use the `@RG` read group tag to split reads. Each read group is then considered a sample. "
        "Reads with an invalid (not in the header) read group tag or without a tag are ignored."
    );
    sam_split_by_rg_.option->group( group );
    sam_split_by_rg_.option->needs( file_input_.option() );

    // Flags include all
    sam_flags_include_all_.option = sub->add_option(
        "--sam-flags-include-all",
        sam_flags_include_all_.value,
        "Only use reads with all bits in the given value present in the FLAG field of the read. "
        "This is equivalent to the `-f` / `--require-flags` setting in `samtools view`, "
        "and uses the same flag names and their corresponding binary values. "
        "The value can be specified in hex by beginning with `0x` (i.e., `/^0x[0-9A-F]+/`), "
        "in octal by beginning with `0` (i.e., `/^0[0-7]+/`), as a decimal number not beginning "
        "with '0', or as a comma-, plus-, space-, or vertiacal-bar-separated list of flag names "
        "as specified by samtools. "
        "We are more lenient in parsing flag names than `samtools`, and allow different "
        "capitalization and delimiteres such as dashes and underscores in the flag names as well."
    );
    sam_flags_include_all_.option->group( group );
    sam_flags_include_all_.option->needs( file_input_.option() );

    // Flags include any
    sam_flags_include_any_.option = sub->add_option(
        "--sam-flags-include-any",
        sam_flags_include_any_.value,
        "Only use reads with any bits set in the given value present in the FLAG field of the read. "
        "This is equivalent to the `--rf` / `--incl-flags` / `--include-flags` setting in "
        "`samtools view`. See `--sam-flags-include-all` above for how to specify the value."
    );
    sam_flags_include_any_.option->group( group );
    sam_flags_include_any_.option->needs( file_input_.option() );

    // Flags exclude all
    sam_flags_exclude_all_.option = sub->add_option(
        "--sam-flags-exclude-all",
        sam_flags_exclude_all_.value,
        "Do not use reads with all bits set in the given value present in the FLAG field of the read. "
        "This is equivalent to the `-G` setting in `samtools view`. "
        "See `--sam-flags-include-all` above for how to specify the value."
    );
    sam_flags_exclude_all_.option->group( group );
    sam_flags_exclude_all_.option->needs( file_input_.option() );

    // Flags exclude any
    sam_flags_exclude_any_.option = sub->add_option(
        "--sam-flags-exclude-any",
        sam_flags_exclude_any_.value,
        "Do not use reads with any bits set in the given value present in the FLAG field of the read. "
        "This is equivalent to the `-F` / `--excl-flags` / `--exclude-flags` setting in "
        "`samtools view`. See `--sam-flags-include-all` above for how to specify the value."
    );
    sam_flags_exclude_any_.option->group( group );
    sam_flags_exclude_any_.option->needs( file_input_.option() );

    return file_input_.option();
}

// =================================================================================================
//      Run Functions
// =================================================================================================

VariantInputSamOptions::VariantInputIterator VariantInputSamOptions::get_iterator_(
    std::string const& filename
) const {
    using namespace genesis::population;

    // Prepare the reader with all its settings.
    SamVariantInputIterator reader;
    reader.min_map_qual( sam_min_map_qual_.value );
    reader.min_base_qual( sam_min_base_qual_.value );
    reader.split_by_rg( sam_split_by_rg_.value );
    reader.flags_include_all( string_to_sam_flag( sam_flags_include_all_.value ));
    reader.flags_include_any( string_to_sam_flag( sam_flags_include_any_.value ));
    reader.flags_exclude_all( string_to_sam_flag( sam_flags_exclude_all_.value ));
    reader.flags_exclude_any( string_to_sam_flag( sam_flags_exclude_any_.value ));

    // Make an iterator.
    return make_variant_input_iterator_from_sam_file( filename, reader );
}
