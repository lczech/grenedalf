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

#include "commands/simulate.hpp"
#include "options/global.hpp"
#include "tools/cli_setup.hpp"

#include "genesis/population/functions/functions.hpp"
#include "genesis/sequence/functions/quality.hpp"
#include "genesis/utils/text/convert.hpp"
#include "genesis/utils/text/string.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#include <random>
#include <stdexcept>
#include <unordered_set>
#include <utility>

// =================================================================================================
//      Setup
// =================================================================================================

void setup_simulate( CLI::App& app )
{
    // Create the options and subcommand objects.
    auto options = std::make_shared<SimulateOptions>();
    auto sub = app.add_subcommand(
        "simulate",
        "Create a file with simulated random frequency data."
    );

    // Format
    options->format.option = sub->add_option(
        "--format",
        options->format.value,
        "Select the output file format, either (m)pileup, or PoPoolation2 sync."
    );
    options->format.option->group( "Settings" );
    options->format.option->transform(
        CLI::IsMember({ "pileup", "sync" }, CLI::ignore_case )
    );

    // Random Seed
    options->random_seed.option = sub->add_option(
        "--random-seed",
        options->random_seed.value,
        "Set the random seed for generating values, which allows reproducible results. "
        "If not provided, the system clock is used to obtain a random seed."
    );
    options->random_seed.option->group( "Settings" );

    // Pool Sizes
    options->coverages.option = sub->add_option(
        "--coverages",
        options->coverages.value,
        "Coverages of the samples to simulate, as a comma- or tab-separated list. "
        "The coverage of each sample is used at the total count per position to randomly distribute "
        "across nucleotides. Per sample, the list can either contain a single number, which will "
        "be used as the coverage for that sample at each position, or it can be two numbers "
        "separated by a slash, which will be used as min/max to generate random coverage at "
        "each position. "
        "The length of this list is also used to determine the number of samples to simulate."
    );
    options->coverages.option->group( "Samples" );
    options->coverages.option->required();

    // Chromosome
    options->chromosome.option = sub->add_option(
        "--chromosome",
        options->chromosome.value,
        "Name of the chromosome. This is simply used as the first column in the output file. "
        "At the moment, only one chromosome is supported."
    );
    options->chromosome.option->group( "Genome" );

    // Mutation Rate
    options->mutation_rate.option = sub->add_option(
        "--mutation-rate",
        options->mutation_rate.value,
        "Mutation rate to simulate. This rate times the `--length` is used as the number of "
        "mutations to generate in total (which can alternatively be directly provided via "
        "`--mutation-count`)."
    );
    options->mutation_rate.option->group( "Genome" );
    options->mutation_rate.option->check( CLI::Range( 0.0, 1.0 ) & CLI::PositiveNumber );

    // Mutation Count
    options->mutation_count.option = sub->add_option(
        "--mutation-count",
        options->mutation_count.value,
        "Number of mutations to simulate in total across the chromosome, "
        "spread across the `--length`."
    );
    options->mutation_count.option->group( "Genome" );

    // Only one way to provide mutations
    options->mutation_rate.option->excludes( options->mutation_count.option );
    options->mutation_count.option->excludes( options->mutation_rate.option );

    // Length
    options->length.option = sub->add_option(
        "--length",
        options->length.value,
        "Total length of the chromosome to simulate. Mutations are spread across this length."
    );
    options->length.option->group( "Genome" );
    options->length.option->required();

    // Omit non mutated positions
    options->omit_invariant_positions.option = sub->add_flag(
        "--omit-invariant-positions",
        options->omit_invariant_positions.value,
        "If set, only write the mutated positions in the output file. Note that these are not "
        "standard (m)pileup or sync files any more; still this option might be useful."
    );
    options->omit_invariant_positions.option->group( "Genome" );

    // With Quality Scores
    options->with_quality_scores.option = sub->add_flag(
        "--with-quality-scores",
        options->with_quality_scores.value,
        "If set, phred-scaled quality scores are written when simulating an (m)pileup file, "
        "using the `--min-phred-score` and `--max-phred-score` settings. Ignored otherwise."
    );
    options->with_quality_scores.option->group( "Pileup" );

    // Min Phred Score
    options->min_phred_score.option = sub->add_option(
        "--min-phred-score",
        options->min_phred_score.value,
        "Minimum phred score to use when simulating an (m)pileup file. Ignored otherwise."
    );
    options->min_phred_score.option->group( "Pileup" );
    options->min_phred_score.option->check(
        CLI::Range( static_cast<size_t>(0), static_cast<size_t>(90) )
    );

    // Max Phred Score
    options->max_phred_score.option = sub->add_option(
        "--max-phred-score",
        options->max_phred_score.value,
        "Maximum phred score to use when simulating an (m)pileup file. Ignored otherwise."
    );
    options->max_phred_score.option->group( "Pileup" );
    options->max_phred_score.option->check(
        CLI::Range( static_cast<size_t>(0), static_cast<size_t>(90) )
    );

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
            run_simulate( *options );
        }
    ));
}

// =================================================================================================
//      Helper Functions
// =================================================================================================

/**
 * @brief List of possible output formats, for speedy access.
 */
enum class SimulateFormat
{
    kPileup,
    kSync
};

SimulateFormat get_format( SimulateOptions const& options )
{
    using namespace genesis::utils;

    if( to_lower( options.format.value ) == "pileup" ) {
        return SimulateFormat::kPileup;
    } else if( to_lower( options.format.value ) == "sync" ) {
        return SimulateFormat::kSync;
    } else {
        throw std::runtime_error( "Internal error: Invalid format " + options.format.value );
    }
}

/**
 * @brief Parse the coverages option, which gives the number of samples, and their coverages.
 */
std::vector<std::pair<size_t,size_t>> get_coverages_( SimulateOptions const& options )
{
    using namespace genesis::utils;

    auto const vals = split( options.coverages.value, ",\t" );
    auto result = std::vector<std::pair<size_t,size_t>>( vals.size() );
    for( size_t i = 0; i < vals.size(); ++i ) {
        try {
            auto const minmax = split( trim( vals[i] ), "/" );
            if( minmax.size() == 1 ) {
                auto const num = convert_from_string<size_t>( trim( minmax[0] ));
                result[i] = { num, num };
            } else if( minmax.size() == 2 ) {
                auto const min = convert_from_string<size_t>( trim( minmax[0] ));
                auto const max = convert_from_string<size_t>( trim( minmax[1] ));
                if( min > max ) {
                    throw std::runtime_error("Invalid coverage value with min > max");
                }
                result[i] = { min, max };
            } else {
                throw std::runtime_error("Invalid coverage value.");
            }
        } catch(...) {
            throw CLI::ValidationError(
                options.coverages.option->get_name(),
                "Invalid coverage value: " + vals[i]
            );
        }
    }
    return result;
}

/**
 * @brief Turn a number from the distributions into a nucleotide.
 */
char allele_to_char_( size_t a )
{
    switch( a ) {
        // Use the order of sync files. If we add more file types later,
        // we might need to adjust the approach.
        case 0: return 'A';
        case 1: return 'T';
        case 2: return 'C';
        case 3: return 'G';
        default: {
            throw std::runtime_error( "Internal error: Invalid allele " + std::to_string(a) );
        }
    }
}

// =================================================================================================
//      Run
// =================================================================================================

void run_simulate( SimulateOptions const& options )
{
    using namespace genesis::sequence;
    using namespace genesis::utils;

    // -------------------------------------------------------------------------
    //     Set up and checks
    // -------------------------------------------------------------------------

    // User friendly check.
    options.file_output.check_output_files_nonexistence( "simulate", options.format.value );

    // Get the format, and store it in a an enum for access speed.
    auto const format = get_format( options );

    // Get a random seed, either from the option if used, or using the current time.
    auto const seed
        = *options.random_seed.option
        ? options.random_seed.value
        : std::chrono::system_clock::now().time_since_epoch().count()
    ;
    std::default_random_engine engine(seed);

    // Check phred score validity.
    if( options.min_phred_score.value > options.max_phred_score.value ) {
        throw CLI::ValidationError(
            options.min_phred_score.option->get_name() + " " +
            options.max_phred_score.option->get_name(),
            "Invalid phred score values with min (" + std::to_string( options.min_phred_score.value ) +
            ") > max (" + std::to_string( options.max_phred_score.value ) + ")."
        );
    }

    // Get coverages and sample count (length of coverage list).
    auto const sample_coverages = get_coverages_( options );

    // Get and check mutation count to avoid infinite loop when picking mutation positions.
    if( ! *options.mutation_count.option && ! *options.mutation_rate.option ) {
        throw CLI::ValidationError(
            "Either " + options.mutation_count.option->get_name() + " or " +
            options.mutation_rate.option->get_name() + " have to be provided."
        );
    }
    size_t const mutation_count
        = *options.mutation_count.option
        ? options.mutation_count.value
        : static_cast<size_t>( options.mutation_rate.value * options.length.value )
    ;
    if( options.mutation_count.value >= options.length.value ) {
        throw CLI::ValidationError(
            options.mutation_count.option->get_name() + "(" +
            std::to_string( options.mutation_count.value ) + ")",
            "Cannot create more mutation positions than genome length."
        );
    }
    LOG_MSG << "Generating genome of length " << options.length.value << " with a total of "
            << mutation_count << " mutated positions, for " << sample_coverages.size()
            << " sample" << ( sample_coverages.size() == 1 ? "" : "s" ) << ".";
    if( mutation_count == 0 ) {
        LOG_WARN << "Zero mutations will be produced. All positions in the file will be invariant. "
                 << "Consider increasing " << options.mutation_rate.option->get_name();
    }

    // Use an unordered_set to collect all positions that we simulate as mutations beforehand.
    // Generating random positions a priori allows us to ensure that we hit the exact genome length.
    // Alternatively, we could use a geometric distribution to draw "skip" distances between
    // mutations, but this might by chance lead to different lengths of the resulting genome.
    // We use a set here, so that we can draw positions without duplicates until we have enough.
    std::unordered_set<size_t> mutation_positions;
    std::uniform_int_distribution<size_t> mutation_positions_distrib( 1, options.length.value );
    while( mutation_positions.size() < mutation_count ) {
        mutation_positions.emplace( mutation_positions_distrib( engine ));
    }

    // Prepare distributions for allele identities, allele frequencies, and phred scores.
    std::uniform_int_distribution<size_t> first_allele_distrib( 0, 3 );
    std::uniform_int_distribution<size_t> second_allele_distrib( 1, 3 );
    std::uniform_real_distribution<double> allele_freq_distrib( 0.0, 1.0 );
    std::uniform_int_distribution<unsigned char> phred_scores_distrib(
        options.min_phred_score.value,
        options.max_phred_score.value
    );

    // Prepare distributions for per sample coverages.
    std::vector<std::uniform_int_distribution<size_t>> sample_coverage_distribs;
    for( auto const& coverage : sample_coverages ) {
        sample_coverage_distribs.emplace_back( coverage.first, coverage.second );
    }

    // -------------------------------------------------------------------------
    //     Run the generation
    // -------------------------------------------------------------------------

    // We iterate all positions, but might skip invariant as needed.
    auto sim_ofs = options.file_output.get_output_target( "simulate", options.format.value );
    for( size_t position = 1; position <= options.length.value; ++position ) {

        // Might not want to write to the file if this is an invariant position.
        if( options.omit_invariant_positions.value && mutation_positions.count( position ) == 0 ) {
            continue;
        }

        // Draw two alleles (we are only doing biallelic for now).
        // The second is drawn so that it differes from the first one.
        // It might not be used if this is an invariant position. We draw it anyway.
        auto const a1 = first_allele_distrib( engine );
        auto const a2 = ( a1 + second_allele_distrib( engine ) ) % 4;

        // Write fixed columns: chromosome, position, referece base.
        // Same for pileup and sync. If more formats are added, this might need to be changed.
        (*sim_ofs) << options.chromosome.value;
        (*sim_ofs) << "\t" << position;
        (*sim_ofs) << "\t" << allele_to_char_( a1 );

        // Go through all samples.
        for( size_t s = 0; s < sample_coverages.size(); ++s ) {
            // Simulate a coverage for the sample.
            auto const coverage = sample_coverage_distribs[s]( engine );

            // Draw major allele frequency for the sample. If this is an invariant position,
            // we use a frequency of 1, and with that can simply use the same code for
            // mutated and invariant positions for the writing below.
            auto const fraction
                = mutation_positions.count( position ) == 0
                ? 1.0
                : allele_freq_distrib( engine )
            ;

            // Distribute the coverage to the two alleles. For simplicity, we assume
            // read count = coverage at every position. Might elaborate on this in the future.
            // We use an array of counts to store them, in the order of sync files, ATCG N D.
            // This works for now. For pileup, we just access the two alleles that we need,
            // but ignore the rest.
            std::array<size_t, 6> counts = {{ 0, 0, 0, 0, 0, 0 }};
            counts[ a1 ] = coverage * fraction;
            counts[ a2 ] = coverage - counts[ a1 ];

            // Write sample.
            switch( format ) {
                case SimulateFormat::kPileup: {
                    // Pileup needs repeated chars of each of the two alleles.
                    // For simplicity, we just output both of them in order.
                    (*sim_ofs) << "\t" << coverage << "\t";
                    (*sim_ofs) << std::string( counts[ a1 ], allele_to_char_( a1 ) );
                    (*sim_ofs) << std::string( counts[ a2 ], allele_to_char_( a2 ) );

                    // Draw and write some random quality scores if needed.
                    if( options.with_quality_scores.value ) {
                        (*sim_ofs) << "\t";
                        for( size_t i = 0; i < coverage; ++i ) {
                            auto const score = phred_scores_distrib( engine );
                            (*sim_ofs) << quality_encode_from_phred_score( score );
                        }
                    }
                    break;
                }
                case SimulateFormat::kSync: {
                    // Sync is simpler, and just needs the counts of each of the 6 different values
                    // (ACGT, as well as N and D for deletions).
                    (*sim_ofs) << "\t" << counts[0];
                    for( size_t i = 1; i < 6; ++i ) {
                        (*sim_ofs) << ":" << counts[i];
                    }
                    break;
                }
            }
        }
        (*sim_ofs) << "\n";
    }
}
