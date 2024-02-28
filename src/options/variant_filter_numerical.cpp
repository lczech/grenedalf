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

#include "options/variant_filter_numerical.hpp"

#include "options/global.hpp"
#include "tools/misc.hpp"

#include <cassert>
#include <stdexcept>

// =================================================================================================
//      Setup Functions
// =================================================================================================

void VariantFilterNumericalOptions::add_numerical_filter_opts_to_app(
    CLI::App* sub,
    std::string const& group
) {
    // Add both types of filters.
    add_sample_filter_opts_to_app( sub, group );
    add_total_filter_opts_to_app( sub, group );
}

// ========================================================
//     Sample Filters
// ========================================================

void VariantFilterNumericalOptions::add_sample_filter_opts_to_app(
    CLI::App* sub,
    std::string const& group
) {
    // Add all three types of sample filters.
    add_sample_count_filter_opts_to_app(    sub, true, true, group );
    add_sample_coverage_filter_opts_to_app( sub, true, true, group );
    add_sample_snp_filter_opts_to_app(      sub, true, true, group );
}


void VariantFilterNumericalOptions::add_sample_count_filter_opts_to_app(
    CLI::App* sub,
    bool add_sample_min_count,
    bool add_sample_max_count,
    std::string const& group
) {
    // Add option for base/allele min count.
    if( add_sample_min_count ) {
        sample_min_count.option = sub->add_option(
            "--filter-sample-min-count",
            sample_min_count.value,
            "Minimum base count for a nucleotide (in `ACGT`) to be considered as an allele. "
            "Counts below that are set to zero, and hence ignored as an allele/variant. "
            "For example, singleton read sequencing errors can be filtered out this way. "
            "This min count is also used for filtering base counts with "
            "too many deletions, unless `--filter-sample-tolerate-deletions` is set."
        );
        sample_min_count.option->group( group );
    }

    // Add option for base/allele max count.
    if( add_sample_max_count ) {
        sample_max_count.option = sub->add_option(
            "--filter-sample-max-count",
            sample_max_count.value,
            "Maximum base count for a nucleotide (in `ACGT`) to be considered as an allele. "
            "Counts above that are set to zero, and hence ignored as an allele/variant. "
            "For example, spuriously high read counts can be filtered out this way."
        );
        sample_max_count.option->group( group );
    }

    // Add tolerate deletions flag
    sample_tolerate_deletions.option = sub->add_flag(
        "--filter-sample-tolerate-deletions",
        sample_tolerate_deletions.value,
        "By default, we filter out any samples that have too many deletions (based on "
        "`--filter-sample-min-count`). With this flag however, such samples are retained."
    );
    sample_tolerate_deletions.option->group( group );
}

void VariantFilterNumericalOptions::add_sample_coverage_filter_opts_to_app(
    CLI::App* sub,
    bool add_sample_min_coverage,
    bool add_sample_max_coverage,
    std::string const& group
) {
    // Add option for min coverage.
    if( add_sample_min_coverage ) {
        sample_min_coverage.option = sub->add_option(
            "--filter-sample-min-coverage",
            sample_min_coverage.value,
            "Minimum coverage expected for a position in a sample to be considered covered. "
            "If the sum of nucleotide counts (in `ACGT`) at a given position in a sample "
            "is less than the provided value, the sample is ignored at this position."
        );
        sample_min_coverage.option->group( group );
    }

    // Add option for max coverage.
    if( add_sample_max_coverage ) {
        sample_max_coverage.option = sub->add_option(
            "--filter-sample-max-coverage",
            sample_max_coverage.value,
            "Maximum coverage expected for a position in a sample to be considered covered. "
            "If the sum of nucleotide counts (in `ACGT`) at a given position in a sample "
            "is greater than the provided value, the sample is ignored at this position. "
            "This can for example be used to filter spuriously high coverage positions."
        );
        sample_max_coverage.option->group( group );
    }
}

void VariantFilterNumericalOptions::add_sample_snp_filter_opts_to_app(
    CLI::App* sub,
    bool add_sample_only_snps,
    bool add_sample_only_biallelic_snps,
    std::string const& group
) {
    // Add flag for SNPs only.
    if( add_sample_only_snps ) {
        sample_only_snps.option = sub->add_flag(
            "--filter-sample-only-snps",
            sample_only_snps.value,
            "Filter out any positions in a sample that do not have two or more alleles "
            "(i.e., that are invariant). "
            "That is, after applying `" + sample_min_count.option->get_name() + "` and `" +
            sample_max_count.option->get_name() + "`, if less than two counts (in `ACGT`) are "
            "non-zero, the position is not considered a SNP for the sample, and ignored."
        );
        sample_only_snps.option->group( group );
    }

    // Add flag for biallelic SNPs only.
    if( add_sample_only_biallelic_snps ) {
        sample_only_biallelic_snps.option = sub->add_flag(
            "--filter-sample-only-biallelic-snps",
            sample_only_biallelic_snps.value,
            "Filter out any positions in a sample that do not have exactly two alleles. "
            "That is, after applying `" + sample_min_count.option->get_name() + "` and `" +
            sample_max_count.option->get_name() + "`, if not exactly two counts (in `ACGT`) are "
            "non-zero, the position is not considered a biallelic SNP for the sample, and ignored."
        );
        sample_only_biallelic_snps.option->group( group );
    }

    // Can only want SNPs or biallelic SNPs.
    if( add_sample_only_snps && add_sample_only_biallelic_snps ) {
        sample_only_snps.option->excludes( sample_only_biallelic_snps.option );
        sample_only_biallelic_snps.option->excludes( sample_only_snps.option );
    }
}

// ========================================================
//     Total Filters
// ========================================================

void VariantFilterNumericalOptions::add_total_filter_opts_to_app(
    CLI::App* sub,
    std::string const& group
) {
    // Add all three types of sample filters.
    add_total_coverage_filter_opts_to_app( sub, true, true, group );
    add_total_snp_filter_opts_to_app(      sub, true, true, group );
    add_total_snp_count_opts_to_app(       sub, true, true, group );
    add_total_freq_filter_opts_to_app(     sub, group );
}

void VariantFilterNumericalOptions::add_total_coverage_filter_opts_to_app(
    CLI::App* sub,
    bool add_total_min_coverage,
    bool add_total_max_coverage,
    std::string const& group
) {
    // Add option for min coverage.
    if( add_total_min_coverage ) {
        total_min_coverage.option = sub->add_option(
            "--filter-total-min-coverage",
            total_min_coverage.value,
            "Minimum coverage expected for a position in total to be considered covered. "
            "If the sum of nucleotide counts (in `ACGT`) at a given position in total "
            "(across all samples) is less than the provided value, the position is ignored."
        );
        total_min_coverage.option->group( group );
    }

    // Add option for max coverage.
    if( add_total_max_coverage ) {
        total_max_coverage.option = sub->add_option(
            "--filter-total-max-coverage",
            total_max_coverage.value,
            "Maximum coverage expected for a position in total to be considered covered. "
            "If the sum of nucleotide counts (in `ACGT`) at a given position in total "
            "(across all samples) is greater than the provided value, the position is ignored. "
            "This can for example be used to filter spuriously high coverage positions."
        );
        total_max_coverage.option->group( group );
    }
}

void VariantFilterNumericalOptions::add_total_snp_filter_opts_to_app(
    CLI::App* sub,
    bool add_total_only_snps,
    bool add_total_only_biallelic_snps,
    std::string const& group
) {
    // Add flag for SNPs only.
    if( add_total_only_snps ) {
        total_only_snps.option = sub->add_flag(
            "--filter-total-only-snps",
            total_only_snps.value,
            "Filter out any positions that do not have two or more alleles (i.e., that are invariant) "
            "across all samples. "
            "That is, after applying all previous filters, if less than two counts (in `ACGT`) are "
            "non-zero in total across all samples, the position is not considered a SNP, and ignored."
        );
        total_only_snps.option->group( group );
    }

    // Add flag for biallelic SNPs only.
    if( add_total_only_biallelic_snps ) {
        total_only_biallelic_snps.option = sub->add_flag(
            "--filter-total-only-biallelic-snps",
            total_only_biallelic_snps.value,
            "Filter out any positions that do not have exactly two alleles across all samples. "
            "That is, after applying all previous filters, if not exactly two counts (in `ACGT`) are "
            "non-zero in total across all samples, the position is not considered a biallelic SNP, "
            "and ignored."
        );
        total_only_biallelic_snps.option->group( group );
    }

    // Can only want SNPs or biallelic SNPs.
    if( add_total_only_snps && add_total_only_biallelic_snps ) {
        total_only_snps.option->excludes( total_only_biallelic_snps.option );
        total_only_biallelic_snps.option->excludes( total_only_snps.option );
    }
}

void VariantFilterNumericalOptions::add_total_snp_count_opts_to_app(
    CLI::App* sub,
    bool add_total_min_count,
    bool add_total_max_count,
    std::string const& group
) {
    // Add min count for snps filter
    if( add_total_min_count ) {
        total_min_count.option = sub->add_option(
            "--filter-total-min-count",
            total_min_count.value,
            "When filtering for positions that are SNPs, use this minimum count (summed across all "
            "samples) to identify what is considered a SNP. Positions where the counts are below "
            "this are filtered out. This min count is also used for filtering positions with "
            "too many deletions, unless `--filter-total-tolerate-deletions` is set."
        );
        total_min_count.option->group( group );
    }

    // Add max count for snps filter
    if( add_total_max_count ) {
        total_max_count.option = sub->add_option(
            "--filter-total-max-count",
            total_max_count.value,
            "When filtering for positions that are SNPs, use this maximum count (summed across all "
            "samples) to identify what is considered a SNP. Positions where the counts are above "
            "this are filtered out; probably not relevant in practice, but offered for completeness."
        );
        total_max_count.option->group( group );
    }

    // Add tolerate deletions flag
    total_tolerate_deletions.option = sub->add_flag(
        "--filter-total-tolerate-deletions",
        total_tolerate_deletions.value,
        "By default, we filter out any positions that have too many deletions (based on "
        "`--filter-total-min-count`). With this flag however, such positions are retained."
    );
    total_tolerate_deletions.option->group( group );
}

void VariantFilterNumericalOptions::add_total_freq_filter_opts_to_app(
    CLI::App* sub,
    std::string const& group
) {
    // Add option for min frequency.
    total_min_frequency.option = sub->add_option(
        "--filter-total-min-frequency",
        total_min_frequency.value,
        "Minimum allele frequency that needs to be reached for a position to be used. "
        "Positions where the allele frequency `af` across all samples, or `1 - af`, "
        "is below this value, are ignored. We compute allele frequency based on the counts "
        "of the reference and alternative base."
    );
    total_min_frequency.option->group( group );

    // TODO add max freq instead as well, for correctly phased input?!
}

// =================================================================================================
//      Run Functions
// =================================================================================================

// -------------------------------------------------------------------------
//     Sample Filter
// -------------------------------------------------------------------------

std::pair<genesis::population::BaseCountsFilter, bool>
VariantFilterNumericalOptions::get_sample_filter(
    genesis::population::BaseCountsFilter filter
) const {
    bool any_provided = false;

    // We here want to check if the option was actually provided by the user,
    // to make it possible that an unspecified option can still be overriden.
    if( sample_min_count.option && *sample_min_count.option ) {
        filter.min_count = sample_min_count.value;
        any_provided = true;
    }
    if( sample_max_count.option && *sample_max_count.option ) {
        filter.max_count = sample_max_count.value;
        any_provided = true;
    }
    if( sample_min_coverage.option && *sample_min_coverage.option ) {
        filter.min_coverage = sample_min_coverage.value;
        any_provided = true;
    }
    if( sample_max_coverage.option && *sample_max_coverage.option ) {
        filter.max_coverage = sample_max_coverage.value;
        any_provided = true;
    }
    if( sample_only_snps.option && *sample_only_snps.option ) {
        filter.only_snps = sample_only_snps.value;
        any_provided = true;
    }
    if( sample_only_biallelic_snps.option && *sample_only_biallelic_snps.option ) {
        filter.only_biallelic_snps = sample_only_biallelic_snps.value;
        any_provided = true;
    }
    if( sample_tolerate_deletions.option && *sample_tolerate_deletions.option ) {
        filter.tolerate_deletions = sample_tolerate_deletions.value;
        any_provided = true;
    }

    return std::make_pair( filter, any_provided );
}

std::function<bool( genesis::population::Variant& )>
VariantFilterNumericalOptions::make_sample_filter() const
{
    // Get the filter settings as provided by the user.
    auto const filter_pair = get_sample_filter();
    auto const filter = filter_pair.first;

    // If none were provided, we can return an empty function.
    if( !filter_pair.second ) {
        return {};
    }

    // If provided, make a filter function.
    return [ filter, this ]( genesis::population::Variant& variant ){
        return genesis::population::filter_base_counts( variant, filter, sample_stats_ );
    };
}

std::function<bool( genesis::population::Variant& )>
VariantFilterNumericalOptions::make_sample_filter(
    genesis::population::BaseCountsFilter filter
) const {
    // Get the filter settings as provided by the user.
    auto const filter_pair = get_sample_filter( filter );
    filter = filter_pair.first;

    // Always make a filter function.
    return [ filter, this ]( genesis::population::Variant& variant ){
        return genesis::population::filter_base_counts( variant, filter, sample_stats_ );
    };
}

// -------------------------------------------------------------------------
//     Total Filter
// -------------------------------------------------------------------------

std::pair<genesis::population::VariantFilter, bool>
VariantFilterNumericalOptions::get_total_filter(
    genesis::population::VariantFilter filter
) const {
    bool any_provided = false;

    // Same as above.
    if( total_min_coverage.option && *total_min_coverage.option ) {
        filter.min_coverage = total_min_coverage.value;
        any_provided = true;
    }
    if( total_max_coverage.option && *total_max_coverage.option ) {
        filter.max_coverage = total_max_coverage.value;
        any_provided = true;
    }
    if( total_only_snps.option && *total_only_snps.option ) {
        filter.only_snps = total_only_snps.value;
        any_provided = true;
    }
    if( total_only_biallelic_snps.option && *total_only_biallelic_snps.option ) {
        filter.only_biallelic_snps = total_only_biallelic_snps.value;
        any_provided = true;
    }
    if( total_min_count.option && *total_min_count.option ) {
        filter.min_count = total_min_count.value;
        any_provided = true;
    }
    if( total_max_count.option && *total_max_count.option ) {
        filter.max_count = total_max_count.value;
        any_provided = true;
    }
    if( total_min_frequency.option && *total_min_frequency.option ) {
        filter.min_frequency = total_min_frequency.value;
        any_provided = true;
    }
    if( total_tolerate_deletions.option && *total_tolerate_deletions.option ) {
        filter.tolerate_deletions = total_tolerate_deletions.value;
        any_provided = true;
    }

    return std::make_pair( filter, any_provided );
}

std::function<bool( genesis::population::Variant const& )>
VariantFilterNumericalOptions::make_total_filter() const
{
    // Get the filter settings as provided by the user.
    auto const filter_pair = get_total_filter();
    auto const filter = filter_pair.first;

    // If none were provided, we can return an empty function.
    if( !filter_pair.second ) {
        return {};
    }

    // If provided, make a filter function.
    return [ filter, this ]( genesis::population::Variant const& variant ){
        return genesis::population::filter_variant( variant, filter, total_stats_ );
    };
}

std::function<bool( genesis::population::Variant const& )>
VariantFilterNumericalOptions::make_total_filter(
    genesis::population::VariantFilter filter
) const {
    // Get the filter settings as provided by the user.
    auto const filter_pair = get_total_filter( filter );
    filter = filter_pair.first;

    // Always make a filter function.
    return [ filter, this ]( genesis::population::Variant const& variant ){
        return genesis::population::filter_variant( variant, filter, total_stats_ );
    };
}

// -------------------------------------------------------------------------
//     print_report
// -------------------------------------------------------------------------

void VariantFilterNumericalOptions::print_report() const
{
    using namespace genesis::population;

    auto const sample_text = print_base_counts_filter_stats( sample_stats_ );
    if( ! sample_text.empty() ) {
        LOG_MSG1 << "Sample filter summary (summed up across all samples):\n"
                 << sample_text << "\n";
    }
    auto const total_text = print_variant_filter_stats( total_stats_ );
    if( ! total_text.empty() ) {
        LOG_MSG1 << "Total filter summary (after applying all sample filters):\n"
                 << total_text << "\n";
    }
}
