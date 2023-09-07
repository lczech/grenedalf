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

#include "options/variant_file_vcf.hpp"

#include "options/global.hpp"
#include "tools/misc.hpp"

#include "genesis/population/formats/vcf_input_iterator.hpp"

#include <cassert>
#include <stdexcept>

// =================================================================================================
//      Setup Functions
// =================================================================================================

CLI::Option* VariantFileVcfOptions::add_file_input_opt_to_app_(
    CLI::App* sub,
    bool required,
    std::string const& group
) {
    // Correct setup check.
    internal_check(
        file_input_.option() == nullptr,
        "Cannot use the same VariantFileVcfOptions object multiple times."
    );

    // Add the option
    file_input_.add_multi_file_input_opt_to_app(
        sub, "vcf", "vcf/bcf", "(vcf(\\.gz)?|bcf)", ".vcf[.gz]|.bcf", required, group,
        "This expects that the input file has the per-sample VCF FORMAT field `AD` (alleleic depth) "
        "given, containing the counts of the reference and alternative base. "
        "This assumes that the data that was used to create the VCF file was actually a pool of "
        "individuals (e.g., from pool sequencing) for each sample (column) of the VCF file. "
        "We then interpret the `AD` field as the allele counts of each pool of individuals. "
        "Note that only SNP positions are used; positions that contain indels and other non-SNP "
        "variants are skipped."
    );

    return file_input_.option();
}

// =================================================================================================
//      Run Functions
// =================================================================================================

VariantFileVcfOptions::VariantInputIterator VariantFileVcfOptions::get_iterator_(
    std::string const& filename
) const {
    using namespace genesis::population;

    // Prepare the iterator. Simple.
    return make_variant_input_iterator_from_pool_vcf_file( filename );
}
