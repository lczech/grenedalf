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

#include "options/variant_input_vcf.hpp"

#include "options/global.hpp"
#include "options/variant_input_sample_names.hpp"
#include "tools/misc.hpp"

#include "genesis/population/formats/vcf_input_iterator.hpp"

#include <cassert>
#include <stdexcept>

// =================================================================================================
//      Setup Functions
// =================================================================================================

CLI::Option* VariantInputVcfOptions::add_vcf_input_opt_to_app(
    CLI::App* sub,
    bool required,
    std::string const& group
) {
    // Correct setup check.
    internal_check(
        vcf_file_.option() == nullptr,
        "Cannot use the same VariantInputVcfOptions object multiple times."
    );

    // Add the option
    vcf_file_.add_multi_file_input_opt_to_app(
        sub, "vcf", "vcf/bcf", "(vcf(\\.gz)?|bcf)", ".vcf[.gz]|.bcf", required, group,
        "This expects that the input file has the per-sample VCF FORMAT field `AD` (alleleic depth) "
        "given, containing the counts of the reference and alternative base. "
        "This assumes that the data that was used to create the VCF file was actually a pool of "
        "individuals (e.g., from pool sequencing) for each sample (column) of the VCF file. "
        "We then interpret the `AD` field as the allele counts of each pool of individuals."
    );

    return vcf_file_.option();
}

// =================================================================================================
//      Run Functions
// =================================================================================================

VariantInputVcfOptions::VariantInputIterator VariantInputVcfOptions::prepare_vcf_iterator(
    std::string const& filename,
    VariantInputSampleNamesOptions const& sample_names_options
) const {
    using namespace genesis::population;

    // Prepare the iterator.
    // See if we want to filter by sample name, and if so, resolve the name list.
    // By default, this also already filters for biallelic SNPs.
    VariantInputIterator iterator;
    bool const only_biallelic = true;
    if( ! sample_names_options.get_filter_samples_include().value.empty() ) {
        auto const list = sample_names_options.process_sample_name_list_option(
            sample_names_options.get_filter_samples_include().value
        );
        iterator = make_variant_input_iterator_from_pool_vcf_file(
            filename, list, false, only_biallelic
        );
    } else if( ! sample_names_options.get_filter_samples_exclude().value.empty() ) {
        auto const list = sample_names_options.process_sample_name_list_option(
            sample_names_options.get_filter_samples_exclude().value
        );
        iterator = make_variant_input_iterator_from_pool_vcf_file(
            filename, list, true, only_biallelic
        );
    } else {
        iterator = make_variant_input_iterator_from_pool_vcf_file(
            filename, only_biallelic
        );
    }

    // As opposed to most other file formats, VCF contains sample names (only the filtered ones).
    // So here we do not need to set them, and can directly return.
    return iterator;
}
