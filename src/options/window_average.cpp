/*
    grenedalf - Genome Analyses of Differential Allele Frequencies
    Copyright (C) 2020-2024 Lucas Czech

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
    Lucas Czech <lucas.czech@sund.ku.dk>
    University of Copenhagen, Globe Institute, Section for GeoGenetics
    Oster Voldgade 5-7, 1350 Copenhagen K, Denmark
*/

#include "options/window_average.hpp"

#include "options/global.hpp"
#include "tools/misc.hpp"

#include "genesis/utils/text/string.hpp"
#include "genesis/population/format/bed_reader.hpp"
#include "genesis/population/function/functions.hpp"
#include "genesis/population/function/genome_locus_set.hpp"
#include "genesis/sequence/functions/dict.hpp"
#include "genesis/utils/core/fs.hpp"
#include "genesis/utils/io/input_source.hpp"
#include "genesis/utils/text/convert.hpp"
#include "genesis/utils/text/string.hpp"

#include <cassert>
#include <stdexcept>

// =================================================================================================
//      Enum Mapping
// =================================================================================================

using namespace genesis::population;

std::vector<std::pair<std::string, WindowAveragePolicy>> const window_average_policy_map = {
    { "window-length",  WindowAveragePolicy::kWindowLength },
    { "available-loci", WindowAveragePolicy::kAvailableLoci },
    { "valid-loci",     WindowAveragePolicy::kValidLoci },
    { "valid-snps",     WindowAveragePolicy::kValidSnps },
    { "sum",            WindowAveragePolicy::kSum },
    { "provided-loci",  WindowAveragePolicy::kProvidedLoci }
};

// =================================================================================================
//      Setup Functions
// =================================================================================================

CLI::Option* WindowAverageOptions::add_window_average_opt_to_app(
    CLI::App* sub,
    VariantReferenceGenomeOptions const& ref_genome_opts,
    bool required,
    std::string const& group
) {
    // Correct setup check.
    internal_check(
        window_average_policy_.option == nullptr,
        "Cannot use the same WindowAverageOptions object multiple times."
    );
    internal_check(
        !ref_genome_opts_ || ref_genome_opts_ == &ref_genome_opts,
        "Cannot use the same VariantFilterMaskOptions with different VariantReferenceGenomeOptions"
    );
    ref_genome_opts_ = &ref_genome_opts;

    // -------------------------------------------------------------------------
    //     Policy
    // -------------------------------------------------------------------------

    // Add the option. We need quite the bit of documentation here...
    window_average_policy_.option = sub->add_option(
        "--window-average-policy",
        window_average_policy_.value,
        "Denominator to use when computing the average of a metric in a window: "
        "\n(1) `window-length`: Simply use the window length, which likely underestimates the metric, "
        "in particular in regions with low coverage and high missing data."
        "\n(2) `available-loci`: Use the number of positions for which there was data at all "
        "(that is, absent or missing data is excluded), independent of all other filter settings."
        "\n(3) `valid-loci`: Use the number of positions that passed all quality and numerical "
        "filters (that is, excluding the SNP-related filters). This uses all positions of high "
        "quality, and is the recommended policy when the input contains data for all positions."
        "\n(4) `valid-snps`: Use the number of SNPs only. This might overestimate the metric, "
        "but can be useful when the data only consists of SNPs."
        "\n(5) `sum`: Simply report the sum of the per-site values, with no averaging "
        "applied to it. This can be used to apply custom averaging later."
        "\n(6) `provided-loci`: Use the exact loci provided via `--window-average-loci-bed` or "
        "`--window-average-loci-fasta` to determine the denominator for the window averaging, "
        "by counting all positions set in this mask in the given window."
    );
    window_average_policy_.option->transform(
        CLI::IsMember( enum_map_keys( window_average_policy_map ), CLI::ignore_case )
    );
    window_average_policy_.option->group( group );
    if( required ) {
        window_average_policy_.option->required();
    }

    // -------------------------------------------------------------------------
    //     Provided Loci BED
    // -------------------------------------------------------------------------

    // Add option for provided loci by BED file.
    window_average_loci_bed_.option = sub->add_option(
        "--window-average-loci-bed",
        window_average_loci_bed_.value,
        "Genomic positions to use for `--window-average-policy provided-loci`, as a BED file. "
        "The regions listed in the BED file are counted towards the window average denominator."
        "\nThis only uses the chromosome, as well as start and end information per line, and "
        "ignores everything else in the file. Note that BED uses 0-based positions, and a half-open "
        "`[)` interval for the end position; simply using columns extracted from other file formats "
        "(such as vcf or gff) will not work."

    );
    window_average_loci_bed_.option->check( CLI::ExistingFile );
    window_average_loci_bed_.option->group( group );
    window_average_loci_bed_.option->needs( window_average_policy_.option );

    // Add inversion option for above mask option.
    window_average_loci_bed_inv_.option = sub->add_flag(
        "--window-average-loci-bed-invert",
        window_average_loci_bed_inv_.value,
        "When using `--window-average-loci-bed`, invert the set of loci. "
        "When this flag is set, the loci that are not set are used for the denominator. "
        "Needs one of " + ref_genome_opts_->get_reference_option_names() +
        " to determine chromosome lengths."
    );
    window_average_loci_bed_inv_.option->group( group );
    window_average_loci_bed_inv_.option->needs( window_average_loci_bed_.option );

    // -------------------------------------------------------------------------
    //     Provided Loci FASTA
    // -------------------------------------------------------------------------

    // Add option for provided loci by fasta-style mask file.
    window_average_loci_fasta_.option = sub->add_option(
        "--window-average-loci-fasta",
        window_average_loci_fasta_.value,
        "Genomic positions to use for `--window-average-policy provided-loci`, as a FASTA-like mask "
        "file (such as used by vcftools).\n"
        "The file contains a sequence of integer digits `[0-9]`, one for each position on the "
        "chromosomes, which specify if the position should be counted towards the window denominator. "
        "Any positions with digits above the `--window-average-loci-fasta-min` value are used. "
    );
    window_average_loci_fasta_.option->check( CLI::ExistingFile );
    window_average_loci_fasta_.option->group( group );
    window_average_loci_fasta_.option->needs( window_average_policy_.option );

    // Add min threshold option for above mask option.
    window_average_loci_fasta_min_.option = sub->add_option(
        "--window-average-loci-fasta-min",
        window_average_loci_fasta_min_.value,
        "When using `--window-average-loci-fasta`, set the cutoff threshold for the counted digits. "
        "All positions above that value are counted. The default is 0, meaning that only exactly "
        "the positons with value 0 will not be counted."
    );
    window_average_loci_fasta_min_.option->check(CLI::Range(0,9));
    window_average_loci_fasta_min_.option->group( group );
    window_average_loci_fasta_min_.option->needs( window_average_loci_fasta_.option );

    // Add inversion option for above mask option.
    window_average_loci_fasta_inv_.option = sub->add_flag(
        "--window-average-loci-fasta-invert",
        window_average_loci_fasta_inv_.value,
        "When using `--window-average-loci-fasta`, invert the set of loci. "
        "When it is set, all positions in the FASTA-like file below or equal to the threshold "
        "are counted towards the window average denominator."
    );
    window_average_loci_fasta_inv_.option->group( group );
    window_average_loci_fasta_inv_.option->needs( window_average_loci_fasta_.option );

    // We only want to be able to specify a single mask file, as everything else would lead to chaos.
    window_average_loci_bed_.option->excludes( window_average_loci_fasta_.option );
    window_average_loci_fasta_.option->excludes( window_average_loci_bed_.option );

    return window_average_policy_.option;
}

// =================================================================================================
//      Run Functions
// =================================================================================================

genesis::population::WindowAveragePolicy WindowAverageOptions::get_window_average_policy() const
{
    // Use the enum map to get the value we want
    auto const result = get_enum_map_value(
        window_average_policy_map,
        window_average_policy_.value
    );

    // User error checks when option not provided.
    // We do not always need it, but when requested here, we need to get it.
    if( !*window_average_policy_.option ) {
        throw CLI::ValidationError(
            window_average_policy_.option->get_name(),
            "Window average policy is required to be set for this combination of options"
        );
    }

    // User error checks for provided loci
    if( result == WindowAveragePolicy::kProvidedLoci ) {
        if( !*window_average_loci_bed_.option && !*window_average_loci_fasta_.option ) {
            throw CLI::ValidationError(
                window_average_policy_.option->get_name(),
                "Window average policy `provided-loci` is set, but no loci file was provided"
            );
        }
    } else {
        if( *window_average_loci_bed_.option ) {
            throw CLI::ValidationError(
                window_average_policy_.option->get_name() + ", " +
                window_average_loci_bed_.option->get_name(),
                "Provided loci BED file is given, but window average policy `provided-loci` is not set"
            );
        }
        if( *window_average_loci_fasta_.option ) {
            throw CLI::ValidationError(
                window_average_policy_.option->get_name() + ", " +
                window_average_loci_fasta_.option->get_name(),
                "Provided loci FASTA file is given, but window average policy `provided-loci` is not set"
            );
        }
    }

    return result;
}

void WindowAverageOptions::prepare_provided_loci_() const
{
    using namespace genesis;
    using namespace genesis::population;
    using namespace genesis::sequence;
    using namespace genesis::utils;
    internal_check(
        ref_genome_opts_, "WindowAverageOptions needs VariantReferenceGenomeOptions"
    );

    // Not running again if we already have set up a filter (i.e., if the shared pointer has data).
    if( provided_loci_ ) {
        return;
    }

    // Add the provided loci from a bed file.
    if( *window_average_loci_bed_.option ) {
        LOG_MSG2 << "Reading provided loci BED file " << window_average_loci_bed_.value;
        provided_loci_ = std::make_shared<GenomeLocusSet>(
            BedReader().read_as_genome_locus_set( from_file( window_average_loci_bed_.value ))
        );
        if( window_average_loci_bed_inv_.value ) {
            auto ref_dict = ref_genome_opts_->get_reference_dict();
            if( !ref_dict ) {
                throw CLI::ValidationError(
                    window_average_loci_bed_inv_.option->get_name(),
                    "Cannot invert BED provided loci without one of " +
                    ref_genome_opts_->get_reference_option_names() +
                    " being provided to determine the chromosome lengths"
                );
            }
            provided_loci_->invert( *ref_dict );
        }
    }

    // Add the provided loci from a fasta file.
    if( *window_average_loci_fasta_.option ) {
        LOG_MSG2 << "Reading provided loci FASTA file " << window_average_loci_fasta_.value;
        provided_loci_ = std::make_shared<GenomeLocusSet>(
            read_mask_fasta(
                from_file( window_average_loci_fasta_.value ),
                window_average_loci_fasta_min_.value,
                window_average_loci_fasta_inv_.value
            )
        );
    }
}
