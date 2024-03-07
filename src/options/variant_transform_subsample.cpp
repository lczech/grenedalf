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
    Lucas Czech <lczech@carnegiescience.edu>
    Department of Plant Biology, Carnegie Institution For Science
    260 Panama Street, Stanford, CA 94305, USA
*/

#include "options/variant_transform_subsample.hpp"

#include "options/global.hpp"
#include "tools/misc.hpp"

#include "genesis/population/functions/subsample.hpp"
#include "genesis/population/functions/variant_input_stream.hpp"
#include "genesis/utils/text/string.hpp"

#include <stdexcept>
#include <string>

// =================================================================================================
//      Setup Functions
// =================================================================================================

void VariantTransformSubsampleOptions::add_subsample_opts_to_app(
    CLI::App* sub,
    std::string const& group
) {
    // Correct setup check.
    internal_check(
        max_coverage_.option == nullptr,
        "Cannot use the same VariantTransformSubsampleOptions object multiple times."
    );

    // Rename samples option.
    max_coverage_.option = sub->add_option(
        "--subsample-max-coverage",
        max_coverage_.value,
        "If provided, the nucleotide counts of each sample are subsampled so that they do not "
        "exceed this given maximum total coverage (sum of the four nucleotide counts). "
        "If they are below this value anyway, they are not changed. "
        "This transformation is useful to limit the maximum coverage. For instance, the diversity "
        "estimators for Theta Pi and Theta Watterson have terms that depend on coverage. "
        "In particular when merging samples such as with `--sample-group-merge-table-file`, "
        "having an upper limit can hence avoid long compute times."
        // "This transformation is applied after any filters, so that, e.g., filters high coverage "
        // "remove any unwanted positions first. See `--subsample-method` for the subsampling method."
    );
    max_coverage_.option->group( group );

    // Add option for sample name filter.
    method_.option = sub->add_option(
        "--subsample-method",
        method_.value,
        "When using `--subsample-max-coverage`, decide which method to use. The default `subscale` "
        "simply re-scales the base counts to the given max coverage, and hence maintains the "
        "allele frequencies (within integer precision). We recommend to use this to subsample to, "
        "e.g., a max coverage of 10,000, which is a good compromise in most cases."
        "The two alternative options re-sample instead, with and without replacement, by drawing "
        "from a multinomial or multivariate hypergeometric distribution, respectively, based on "
        "the original counts of the sample."
    );
    method_.option->group( group );
    method_.option->needs( max_coverage_.option );
    method_.option->transform(
        CLI::IsMember({
            "subscale", "subsample-with-replacement", "subsample-without-replacement"
        }, CLI::ignore_case )
    );
}

// =================================================================================================
//      Run Functions
// =================================================================================================

// -------------------------------------------------------------------------
//     add_transform_subsample_filter
// -------------------------------------------------------------------------

void VariantTransformSubsampleOptions::add_subsample_transformation(
    VariantInputOptions const& variant_input
) const {
    using namespace genesis::population;

    // Nothing to do.
    if( ! max_coverage_.option || ! *max_coverage_.option || max_coverage_.value == 0 ) {
        return;
    }

    // Add a transformation to the variant stream.
    // We need to use a bool return value in the lambdas,
    // as that's expected by the VariantInputOptions.
    auto const method = genesis::utils::to_lower( method_.value );
    auto const max_coverage = max_coverage_.value;
    if( method == "subscale" ) {
        variant_input.add_combined_filter_and_transforms(
            [ max_coverage ]( Variant& variant ){
                transform_subscale( variant, max_coverage );
                return true;
            }
        );
    } else if( method == "subsample-with-replacement" ) {
        variant_input.add_combined_filter_and_transforms(
            [ max_coverage ]( Variant& variant ){
                transform_subsample_with_replacement( variant, max_coverage );
                return true;
            }
        );
    } else if( method == "subsample-without-replacement" ) {
        variant_input.add_combined_filter_and_transforms(
            [ max_coverage ]( Variant& variant ){
                transform_subsample_without_replacement( variant, max_coverage );
                return true;
            }
        );
    } else {
        throw CLI::ValidationError(
            method_.option->get_name() + "(" + method_.value + ")",
            "Invalid value."
        );
    }
}
