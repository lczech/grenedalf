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

#include "options/global.hpp"

#include "tools/version.hpp"

#include "genesis/utils/core/info.hpp"

#include <thread>

// =================================================================================================
//      Setup Functions
// =================================================================================================

void GlobalOptions::initialize( int const argc, char const* const* argv )
{
    // By default, use the hardware threads, taking hypterthreding into account
    opt_threads.value = genesis::utils::guess_number_of_threads();

    // If hardware value is not available, just use 1 thread.
    // This is executed if the call to the above function fails.
    if( opt_threads.value == 0 ) {
        opt_threads.value = 1;
    }

    // Set verbosity to max, just in case.
    genesis::utils::Logging::max_level( genesis::utils::Logging::LoggingLevel::kDebug4 );

    // Store all arguments in the array.
    command_line_.clear();
    for (int i = 0; i < argc; i++) {
        command_line_.push_back(argv[i]);
    }
}

void GlobalOptions::add_to_module( CLI::App& module )
{
    for( auto subcomm : module.get_subcommands({}) ) {
        add_to_subcommand( *subcomm );
    }
}

void GlobalOptions::add_to_subcommand( CLI::App& subcommand )
{
    // Allow to overwrite files.
    opt_allow_file_overwriting = subcommand.add_flag(
        // "--allow-file-overwriting",
        allow_file_overwriting_flag,
        opt_allow_file_overwriting.value,
        "Allow to overwrite existing output files instead of aborting the command."
    );
    opt_allow_file_overwriting.option->group( "Global Options" );

    // Verbosity
    opt_verbose = subcommand.add_flag(
        "--verbose",
        opt_verbose.value,
        "Produce more verbose output."
    );
    opt_verbose.option->group( "Global Options" );

    // Threads
    opt_threads = subcommand.add_option(
        "--threads",
        opt_threads.value,
        "Number of threads to use for calculations. If not set, we guess a reasonable number of "
        "threads, by looking at the environmental variables (1) `OMP_NUM_THREADS` (OpenMP) and "
        "(2) `SLURM_CPUS_PER_TASK` (slurm), as well as (3) the hardware concurrency, taking "
        "hyperthreads into account, in the given order of precedence."
    );
    opt_threads.option->group( "Global Options" );

    // Log File
    opt_log_file = subcommand.add_option(
        "--log-file",
        opt_log_file.value,
        "Write all output to a log file, in addition to standard output to the terminal."
    );
    opt_log_file.option->group( "Global Options" );

    // TODO add random seed option
}

// =================================================================================================
//      Run Functions
// =================================================================================================

void GlobalOptions::run_global()
{
    // If user did not provide number, use hardware value (taking care of hyperthreads as well).
    if( opt_threads.value == 0 ) {
        opt_threads.value = genesis::utils::guess_number_of_threads();
    }

    // If hardware value is not available, just use 1 thread.
    // This is executed if the call above fails.
    if( opt_threads.value == 0 ) {
        opt_threads.value = 1;
    }

    // Initialize the global thread pool. We use one fewer than the number of specified thread
    // here, as we need to acount for the main thread doing work as well. We have implemented
    // a type of proactive future that also does work from the pool when waiting for results,
    // meaning that the main thread will also participate in the pool.
    genesis::utils::Options::get().init_global_thread_pool( opt_threads.value - 1 );

    // Allow to overwrite files. Has to be done before adding the log file (coming below),
    // as this might already fail if the log file exists.
    if( opt_allow_file_overwriting.value ) {
        genesis::utils::Options::get().allow_file_overwriting( true );
    }

    // Set log file.
    if( ! opt_log_file.value.empty() ) {
        genesis::utils::Logging::log_to_file( opt_log_file.value );
    }

    // Set verbosity level for logging output.
    if( opt_verbose.value ) {
        genesis::utils::Logging::max_level( genesis::utils::Logging::LoggingLevel::kMessage2 );
    } else {
        genesis::utils::Logging::max_level( genesis::utils::Logging::LoggingLevel::kMessage1 );
    }
}

// =================================================================================================
//      Getters
// =================================================================================================

std::string GlobalOptions::command_line() const
{
    std::string ret = "";
    for (size_t i = 0; i < command_line_.size(); ++i) {
        ret += ( i==0 ? "" : " " ) + command_line_[i];
    }
    return ret;
}

// =================================================================================================
//      Global Instance
// =================================================================================================

/**
 * @brief Instanciation of the global options object. This is alive during the whole program run.
 */
GlobalOptions global_options;

std::string const allow_file_overwriting_flag = "--allow-file-overwriting";
