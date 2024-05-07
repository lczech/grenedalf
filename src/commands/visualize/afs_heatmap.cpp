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

#include "commands/visualize/afs_heatmap.hpp"
#include "options/global.hpp"
#include "tools/cli_setup.hpp"
#include "tools/misc.hpp"

#include "genesis/population/function/functions.hpp"
#include "genesis/population/plotting/genome_heatmap.hpp"
#include "genesis/population/plotting/heatmap_colorization.hpp"
#include "genesis/utils/math/statistics.hpp"
#include "genesis/utils/color/list_diverging.hpp"
#include "genesis/utils/color/list_sequential.hpp"

#include <algorithm>
#include <cassert>
#include <limits>
#include <stdexcept>
#include <utility>

// #ifdef GENESIS_OPENMP
// #   include <omp.h>
// #endif

// =================================================================================================
//      Enum Mapping
// =================================================================================================

std::vector<std::pair<std::string, FrequencySpectrumType>> const spectrum_type_map = {
    { "folded",   FrequencySpectrumType::kFolded },
    { "unfolded", FrequencySpectrumType::kUnfolded }
};

std::vector<std::pair<std::string, FrequencyAverageMethod>> const average_method_map = {
    { "arithmetic", FrequencyAverageMethod::kArithmetic },
    { "geometric",  FrequencyAverageMethod::kGeometric },
    { "harmonic",   FrequencyAverageMethod::kHarmonic },
    { "counts",     FrequencyAverageMethod::kCounts }
};

// =================================================================================================
//      Setup
// =================================================================================================

void setup_afs_heatmap( CLI::App& app )
{
    // Create the options and subcommand objects.
    auto options = std::make_shared<AfsHeatmapOptions>();
    auto sub = app.add_subcommand(
        "afs-heatmap",
        "Create a per-window heatmap of the allele frequency spectrum along the genome."
    );

    // Required input of some frequency format, and settings for the sliding window.
    options->variant_input.add_variant_input_opts_to_app( sub );
    // options->variant_input.add_sample_name_opts_to_app( sub );
    // options->variant_input.add_region_filter_opts_to_app( sub );
    options->window.add_window_opts_to_app( sub, false );

    // -------------------------------------------------------------------------
    //     Settings
    // -------------------------------------------------------------------------

    // Resolution / number of bins to use
    options->resolution.option = sub->add_option(
        "--resolution",
        options->resolution.value,
        "Resolution of the spectrum histogram, that is, the number of bins of frequencies. "
        "This is hence also the vertical resolution (in pixels) of the resulting images."
    );
    options->resolution.option->group( "Settings" );
    options->resolution.option->check( CLI::PositiveNumber );

    // Max frequency to plot on the y-axis
    options->max_frequency.option = sub->add_option(
        "--max-frequency",
        options->max_frequency.value,
        "Maximum frequency, that is, the y-axis cutoff; default is 1.0. Frequencies above the "
        "maximum will be assinged to the highest bin. When using `--spectrum-type folded`, "
        "consider setting this option to 0.5, as that is the maximum of the folded spectrum."
        // "If not provided, defaults to 1.0 when using `--spectrum-type unfolded` (default), "
        // "and to 0.5 when using `--spectrum-type folded`."
    );
    options->max_frequency.option->group( "Settings" );
    options->max_frequency.option->check( CLI::Range( 0.0, 1.0 ) & CLI::PositiveNumber );

    // TODO color, invert vert, log scale, max per column, empty color

    // Which type of allele to use
    options->spectrum_type.option = sub->add_option(
        "--spectrum-type",
        options->spectrum_type.value,
        "Type of spectrum to compute. That is, which frequencies do the values in the heat map "
        "represent: Either the frequency of the alternative allele (unfolded, default; uses the "
        "highest count/frequency allele that is not the reference allele as given in the input "
        "file, with frequencies in range [0, 1]), "
        "or the frequency of the minor allele (folded; uses the allele with the second highest "
        "count/frequency after the major (most common) allele, with frequencies in range [0, 0.5])."
    );
    options->spectrum_type.option->group( "Settings" );
    options->spectrum_type.option->transform(
        CLI::IsMember( enum_map_keys( spectrum_type_map ), CLI::ignore_case )
    );
    // options->spectrum_type.option->transform(
    //     CLI::CheckedTransformer( spectrum_type_map, CLI::ignore_case )
    // );

    // How to compute the average of the per sample frequencies
    options->average_method.option = sub->add_option(
        "--average-method",
        options->average_method.value,
        "If multiple samples are used (present in the file, and not filtered out), either compute "
        "the average of the frequencies per sample (using either arithmetic, geometric, or harmonic "
        "mean), or sum up the allele counts per sample, and then compute the frequency from that. "
        "The former gives each sample the same weight, while in the latter case, samples with "
        "more counts (more sequence reads) have a higher influence. "
        "If harmonic mean is used, we apply a correction for frequencies of zero, which happens "
        "in samples that do not have any alternative/minor allele counts."
    );
    options->average_method.option->group( "Settings" );
    options->average_method.option->transform(
        CLI::IsMember( enum_map_keys( average_method_map ), CLI::ignore_case )
    );
    // options->average_method.option->transform(
    //     CLI::CheckedTransformer( average_method_map, CLI::ignore_case )
    // );

    // What to do with undetermined ("N") reference allele positions
    options->fold_undetermined.option = sub->add_flag(
        "--fold-undetermined-positions",
        options->fold_undetermined.value,
        "Only relevant if `--spectrum-type unfolded` is set. By default, positions with undetermined "
        "reference alleles ('N') are skipped in the unfolded spectrum. If however this option is "
        "set, these positions are instead folded; that is, we then assume the major allele with the "
        "highest count/frequency to be the reference allele."
    );
    options->fold_undetermined.option->group( "Settings" );

    // Which type of allele to use
    options->individual_bmps.option = sub->add_flag(
        "--write-individual-bmps",
        options->individual_bmps.value,
        "If set, write individual heatmap files for each chromosome, in bitmap format, "
        "in addition to the svg heatmap that contains all (not filtered) chromosomes."
    );
    options->individual_bmps.option->group( "Settings" );

    // -------------------------------------------------------------------------
    //     Output
    // -------------------------------------------------------------------------

    // Output
    options->file_output.add_default_output_opts_to_app( sub );
    options->file_output.add_file_compress_opt_to_app( sub );

    // Set the run function as callback to be called when this subcommand is issued.
    // Hand over the options by copy, so that their shared ptr stays alive in the lambda.
    sub->callback( grenedalf_cli_callback(
        sub,
        {
            // Citation keys as needed
        },
        [ options ]() {
            run_afs_heatmap( *options );
        }
    ));
}

// =================================================================================================
//      Helper Functions
// =================================================================================================

/**
 * @bief Use an enum for the reference/major allele. That is hopefully a tad quicker that
 * char comparisons, because it should allow the compiler to use a jump table.
 */
enum class Nucleotide
{
    kA,
    kC,
    kG,
    kT,
    kN
};

/**
 * @brief Get the main allele of a variant, that is, either the reference base, or the major allele,
 * depending on the setting.
 */
Nucleotide get_main_allele_(
    AfsHeatmapOptions const& options,
    genesis::population::Variant const& variant
) {
    // Helper function to get the enum value from a char.
    // At some point, we could refactor the underyling code to use the enum in the first place...
    // Lots of unnecessary back and forth between data representations here... But good enough for now.
    auto get_nucleotide_ = []( char c ){
        // Test the actual good ones.
        switch( c ) {
            case 'a':
            case 'A': {
                return Nucleotide::kA;
            }
            case 'c':
            case 'C': {
                return Nucleotide::kC;
            }
            case 'g':
            case 'G': {
                return Nucleotide::kG;
            }
            case 't':
            case 'T': {
                return Nucleotide::kT;
            }
        }

        // If none of the above, return a default. We also turn deletions or any others
        // into N here, for simplicity, as they are all treated the same here anyway.
        return Nucleotide::kN;
    };

    // First, do the unfolded case, that is, get the reference base from the variant.
    if( options.spectrum_type_enum == FrequencySpectrumType::kUnfolded ) {
        auto const nc = get_nucleotide_( variant.reference_base );
        if( nc != Nucleotide::kN ) {
            // If we found a proper reference base, we are good, and can return it.
            return nc;
        } else if( ! options.fold_undetermined.value ) {
            // If not, and we also do not want to fold undetermined reference bases,
            // we simply return that we did not find a reference base.
            return Nucleotide::kN;
        }

        // If we are here then there was no good reference base in the variant,
        // but we do want to fold instead and use the major one as main allele.
        // So just continue, and do the major allele finding below.
    }

    // Second, do the folded case, that is, look for the highest count base.
    // Also, if the above did not yield a valid reference base, e.g. the reference base of the
    // variant is 'N' or '.', we instead use the folded approach.
    // First, get the sum of all counts of the variant sorted, most common first.
    auto const order = sorted_sample_counts(
        variant, false, genesis::population::SampleCountsFilterPolicy::kOnlyPassing
    );

    // If the most common one has a count of 0, this is not really a variant, so return N.
    // If not, return the symbol of the most common one.
    if( order.data[0].count == 0 ) {
        assert( order.data[1].count == 0 );
        assert( order.data[2].count == 0 );
        assert( order.data[3].count == 0 );
        return Nucleotide::kN;
    }
    return get_nucleotide_( order.data[0].base );
}

/**
 * @brief Compute the average frequency from a given Variant object,
 * using the folded or unfolded approach.
 */
double compute_frequency_(
    AfsHeatmapOptions const& options,
    genesis::population::Variant const& variant
) {
    using namespace genesis::population;
    using namespace genesis::utils;

    // Get the main allele of the variant, that is, the reference base, or the major allele,
    // depending on the setting. If this is N, there is no frequency to compute.
    auto const main_allele = get_main_allele_( options, variant );
    if( main_allele == Nucleotide::kN ) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    // Helper function to extract the counts of the main (reference or major) allele,
    // and the sum of all alleles, from a given SampleCounts object.
    // Returns a pair with first = main, second = sum of all.
    auto get_bc_counts_ = [main_allele]( SampleCounts const& bc ) -> std::pair<size_t, size_t> {
        size_t all_cnt = bc.a_count + bc.c_count + bc.g_count + bc.t_count;
        switch( main_allele ) {
            case Nucleotide::kA: {
                return std::make_pair( bc.a_count, all_cnt );
            }
            case Nucleotide::kC: {
                return std::make_pair( bc.c_count, all_cnt );
            }
            case Nucleotide::kG: {
                return std::make_pair( bc.g_count, all_cnt );
            }
            case Nucleotide::kT: {
                return std::make_pair( bc.t_count, all_cnt );
            }
            default: {}
        }
        throw std::domain_error( "Internal error: Invalid reference/major allele." );
    };

    if( options.average_method_enum == FrequencyAverageMethod::kCounts ) {

        // Sum up all counts of the samples. We simply sum up all four counts for the "alternative",
        // with no checks, as this is way quicker than first finding the second most common count,
        // and also allows to compute the spectrum for non biallelic variants, if we ever want this.
        // At the moment though, the iterator over the input file already filters for biallelic
        // snps only, so two of the counts will always be zero anyway.
        size_t main_cnt = 0;
        size_t all_cnt = 0;
        for( auto const& bc : variant.samples ) {
            auto const bc_counts = get_bc_counts_( bc );
            main_cnt += bc_counts.first;
            all_cnt  += bc_counts.second;
        }

        // Compute the frequency of the alternative/minor.
        assert( main_cnt <= all_cnt );
        auto const freq = static_cast<double>( all_cnt - main_cnt ) / static_cast<double>( all_cnt );
        assert( 0.0 <= freq && freq <= 1.0 );
        return freq;

    } else {

        // Prepare a buffer for computing the per sample frequencies, to compute their averages.
        // We want that buffer in order to simply use our existing mean functions... we could
        // instead also use transform iterators for this, but that somehow seems more convolulted,
        // so let's keep it simple.
        auto frequencies = std::vector<double>( variant.samples.size() );
        assert( variant.samples.size() == options.variant_input.sample_names().size() );

        // Go through all samples and get their alternative/minor frequencies.
        // Not parallelizing here, as we already parallelized the outer loop over the whole window.
        for( size_t i = 0; i < variant.samples.size(); ++i ) {
            auto const bc_counts = get_bc_counts_( variant.samples[i] );
            assert( bc_counts.first <= bc_counts.second );
            auto const numerator = static_cast<double>( bc_counts.second - bc_counts.first );
            frequencies[i] = numerator / static_cast<double>( bc_counts.second );

            // LOG_MSG << variant.chromosome << ":" << variant.position << " " << i << " --> " << frequencies[i];
        }

        // Average the frequencies.
        switch( options.average_method_enum ) {
            case FrequencyAverageMethod::kArithmetic: {
                return arithmetic_mean( frequencies );
            }
            case FrequencyAverageMethod::kGeometric: {
                return geometric_mean( frequencies );
            }
            case FrequencyAverageMethod::kHarmonic: {
                return harmonic_mean( frequencies, HarmonicMeanZeroPolicy::kCorrection );
            }
            default: {}
        }
    }

    // Make the compiler happy.
    throw std::domain_error( "Internal error: Invalid average frequency method." );
}

// =================================================================================================
//      Run
// =================================================================================================

void run_afs_heatmap( AfsHeatmapOptions const& options )
{
    using namespace genesis::population;
    using namespace genesis::utils;

    options.file_output.check_output_files_nonexistence( "afs-heatmap", "svg" );
    if( options.individual_bmps.value ) {
        options.file_output.check_output_files_nonexistence( "afs-heatmap-*", "bmp" );
    }

    // -------------------------------------------------------------------------
    //     Settings
    // -------------------------------------------------------------------------

    // Check resolution
    if( options.resolution.value <= 1 ) {
        throw CLI::ValidationError(
            options.resolution.option->get_name(),
            "Invalid resolution value: " + std::to_string( options.resolution.value ) + " <= 1"
        );
    }

    // Set the enum values.
    options.spectrum_type_enum  = get_enum_map_value(
        spectrum_type_map, options.spectrum_type.value
    );
    options.average_method_enum = get_enum_map_value(
        average_method_map, options.average_method.value
    );

    // Set max frequency to 0.5 if it is not given and we are using the folded spectrum.
    // double max_frequency = options->max_frequency.value;
    // if(
    //     options.spectrum_type.value == FrequencySpectrumType::kFolded &&
    //     ! *options->max_frequency.option
    // ) {
    //     max_frequency = 0.5;
    // }
    // if( max_frequency <= 0.0 ) {
    //     throw CLI::ValidationError(
    //         options.max_frequency.option->get_name(),
    //         "Invalid maximum frequency value: " + std::to_string( max_frequency )
    //     );
    // }

    // -------------------------------------------------------------------------
    //     Preparation
    // -------------------------------------------------------------------------

    // TODO
    auto colorization = HeatmapColorization( color_list_blues() );
    colorization.log_scale( true );
    colorization.color_map().mask_color( Color(1,1,1) );
    colorization.color_map().under_color( Color(1,1,1) );
    colorization.max_per_column( true );

    GenomeHeatmap heatmap;
    HeatmapColorization::Spectrum spectrum;

    size_t skipped_nonfinites = 0;

    // -------------------------------------------------------------------------
    //     Run
    // -------------------------------------------------------------------------

    // Iterate the file in windows. We then go through all entries in the window again,
    // which makes it a bit wasteful to copy everything into a window first, but makes handling
    // of the surrounding code (of keeping track of positions etc) so much easier.
    // Might optimize in the future.
    auto window_it = options.window.get_variant_window_stream( options.variant_input );
    for( auto cur_it = window_it->begin(); cur_it != window_it->end(); ++cur_it ) {
        auto const& window = *cur_it;

        // Things to do when we start with a new chromosome.
        if( cur_it.is_first_window() ) {
            // Some user output to report progress.
            LOG_MSG << "At chromosome " << window.chromosome();

            // Prepare a new spectrum for this chromosome.
            spectrum = HeatmapColorization::Spectrum();
            spectrum.chromosome = window.chromosome();
        }

        // Beginning of a new window. Add a column to the spectrum.
        spectrum.values.emplace_back( options.resolution.value, 0.0 );
        auto& spectrum_column = spectrum.values.back();

        // #pragma omp parallel for
        for( size_t i = 0; i < window.size(); ++i ) {

            // Compute the frequency at this position, taking all settings into account.
            auto const freq = compute_frequency_( options, window[i].data );

            // Some final checks.
            if( ! std::isfinite(freq) ) {
                // That means we cannot work with that frequency. Low counts here. Llet's just skip.
                ++skipped_nonfinites;
                continue;
            }
            assert( std::isfinite( freq ));
            assert( 0.0 <= freq && freq <= 1.0 );

            // Add the frequency to the respective bin, that is, increase the count in the bin.
            // Bring the value into the bins and store it. We need a special case for the exact value
            // of 1.0 here and all values above that (which can happen if the max frequency set by
            // the user is less than a computed frequency), so that we don't get an overflow
            // of the bin index.
            size_t index = std::floor(
                freq / options.max_frequency.value * static_cast<double>( options.resolution.value )
            );
            index = std::min( index, options.resolution.value - 1 );
            assert( spectrum_column.size() == options.resolution.value );
            assert( index < spectrum_column.size() );

            // Atomic increment for double values.
            // #pragma omp atomic
            spectrum_column[ index ] += 1.0;
        }

        // Things to do when we finish with a chromosome.
        if( cur_it.is_last_window() ) {
            // Add the spectrum to the genome heatmap.
            heatmap.add(
                spectrum.chromosome, colorization.spectrum_to_image( spectrum ).first
            );

            // Write out individual bitmaps for the chromosome.
            if( options.individual_bmps.value ) {
                auto const max_val = colorization.spectrum_to_bmp_file(
                    spectrum,
                    options.file_output.get_output_target(
                        "afs-heatmap-" + spectrum.chromosome, "bmp"
                    )
                );
                if( max_val != 1.0 ) {
                    LOG_MSG << "Max y-axis value for " << spectrum.chromosome
                            << " heatmap bmp: " << max_val;
                }
            }
        }
    }

    if( skipped_nonfinites > 0 ) {
        LOG_MSG << "Skipped " << skipped_nonfinites << " positions due to low counts.";
    }
    heatmap.write( options.file_output.get_output_target( "afs-heatmap", "svg" ));
}
