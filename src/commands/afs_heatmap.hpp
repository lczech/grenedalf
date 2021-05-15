#ifndef GRENEDALF_COMMANDS_AFS_HEATMAP_H_
#define GRENEDALF_COMMANDS_AFS_HEATMAP_H_

/*
    grenedalf - Genome Analyses of Differential Allele Frequencies
    Copyright (C) 2020-2021 Lucas Czech

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

#include "CLI/CLI.hpp"

#include "options/file_output.hpp"
#include "options/frequency_input.hpp"
#include "tools/cli_option.hpp"

#include <string>
#include <vector>

// =================================================================================================
//      Options
// =================================================================================================

/**
 * @brief Helper enum for which kind of allele frequency spectrum to compute,
 * in order to allow direct CLI conversion, and to avoid string comparisons in the inner loop.
 */
enum class FrequencySpectrumType : int
{
    kFolded,
    kUnfolded
};

/**
 * @brief Helper enum for how to compute average allele frequencies,
 * in order to allow direct CLI conversion, and to avoid string comparisons in the inner loop.
 */
enum class FrequencyAverageMethod : int
{
    kArithmetic,
    kGeometric,
    kHarmonic,
    kCounts
};

/**
 * @brief Options class for the allele frequency spectrum heatmap.
 */
class AfsHeatmapOptions
{
public:

    CliOption<size_t> resolution = 100;
    CliOption<double> max_frequency = 1.0;

    // For the enum types, we use a string for the CLI option for nice user output and help,
    // but internally convert to the enum for speed reasons.
    CliOption<std::string> spectrum_type = "unfolded";
    CliOption<std::string> average_method = "harmonic";
    mutable FrequencySpectrumType spectrum_type_enum;
    mutable FrequencyAverageMethod average_method_enum;

    CliOption<bool> fold_undetermined = false;
    CliOption<bool> individual_bmps = false;

    FrequencyInputOptions freq_input;
    FileOutputOptions  file_output;

};

// =================================================================================================
//      Functions
// =================================================================================================

void setup_afs_heatmap( CLI::App& app );
void run_afs_heatmap( AfsHeatmapOptions const& options );

#endif // include guard
