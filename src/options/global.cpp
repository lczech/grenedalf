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
    opt_threads.value = genesis::utils::guess_number_of_threads(false);

    // If hardware value is not available, just use 1 thread.
    // This is executed if the call to the above function fails.
    if( opt_threads.value == 0 ) {
        opt_threads.value = 1;
    }

    // Set number of threads for genesis.
    genesis::utils::Options::get().number_of_threads( opt_threads.value );

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
        opt_threads.value = genesis::utils::guess_number_of_threads(false);
    }

    // If hardware value is not available, just use 1 thread.
    // This is executed if the call above fails.
    if( opt_threads.value == 0 ) {
        opt_threads.value = 1;
    }

    // Set number of threads for genesis. At the moment, this only sets the number of threads
    // to use for OpenMP. We might refactor this in the future, to not use OpenMP any more,
    // because users had trouble getting it to compile, and instead use our own multi-threading
    // solution using a global thread pool for all CPU-heavy computation. See below for that pool;
    // once we have refactored everything to use that pool, this comment here will be outdated...
    genesis::utils::Options::get().number_of_threads( opt_threads.value );

    // At the moment, we use a thread pool for the LambdaIterator of the VariantInputOptions
    // reading (but only for that, see comment above as well), because otherwise each input file
    // spawns a separate thread, which could be hundreds  or a thousand threads all trying to
    // get CPU time to parse their file at the same time... which would invalidate all thread
    // settings, and might cause trouble on clusters.
    // Note that this technique is independent of OpenMP, and so we will end of having two times
    // the number of threads spawned as the moment... Hence the need for the refactor later.
    genesis::utils::Options::get().init_global_thread_pool( opt_threads.value );

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
