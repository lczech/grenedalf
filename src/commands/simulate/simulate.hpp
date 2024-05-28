#ifndef GRENEDALF_COMMANDS_SIMULATE_SIMULATE_H_
#define GRENEDALF_COMMANDS_SIMULATE_SIMULATE_H_

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

#include "CLI/CLI.hpp"

#include "options/file_output.hpp"
#include "tools/cli_option.hpp"

#include <cstdint>
#include <string>
#include <vector>

// =================================================================================================
//      Options
// =================================================================================================

class SimulateOptions
{
public:

    // General options
    CliOption<std::string> format = "pileup";
    CliOption<std::uint32_t> random_seed;

    // Sample read depths, as single numbers, or as min/max entries per sample.
    CliOption<std::string> read_depths;

    // Chromosome and positions
    CliOption<std::string> chromosome = "A";
    CliOption<double> mutation_rate = 1e-8;
    CliOption<size_t> mutation_count;
    CliOption<size_t> length;
    CliOption<bool> omit_invariant_positions = false;

    // Pileup quality scores
    CliOption<bool> with_quality_scores = false;
    CliOption<size_t> min_phred_score = 10;
    CliOption<size_t> max_phred_score = 40;

    FileOutputOptions  file_output;

};

// =================================================================================================
//      Functions
// =================================================================================================

void setup_simulate( CLI::App& app );
void run_simulate( SimulateOptions const& options );

#endif // include guard
