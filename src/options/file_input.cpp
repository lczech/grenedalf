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

#include "options/file_input.hpp"

#include "options/global.hpp"

#include "genesis/utils/core/fs.hpp"
#include "genesis/utils/text/string.hpp"

#include <algorithm>
#include <iostream>
#include <stdexcept>

#ifdef GENESIS_OPENMP
#   include <omp.h>
#endif

// =================================================================================================
//      Setup Functions
// =================================================================================================

CLI::Option* FileInputOptions::add_multi_file_input_opt_to_app(
    CLI::App* sub,
    std::string const& option_name_infix,
    std::string const& option_name_nice,
    std::string const& extension_regex,
    std::string const& extension_nice,
    bool               required,
    std::string const& group,
    std::string const& extra_help
){
    // Correct setup check.
    if( option_ != nullptr ) {
        throw std::domain_error( "Cannot use the same FileInputOptions object multiple times." );
    }

    // Store file type info.
    file_type_ = option_name_infix;
    file_ext_  = extension_regex;

    // Input files.
    option_ = sub->add_option(
        "--" + option_name_infix + "-path",
        raw_paths_,
        "List of " + option_name_nice + " files or directories to process. For directories, " +
        "only files with the extension `." + extension_nice + "` are processed." +
        ( extra_help.empty() ? "" : " " + extra_help )
    );
    if( required ) {
        option_->required();
    }

    // Check if it is a path.
    option_->check( CLI::ExistingPath );
    option_->group( group );

    return option_;
}

// =================================================================================================
//      Run Functions
// =================================================================================================

bool FileInputOptions::provided() const
{
    return option_ && *option_;
}

size_t FileInputOptions::file_count() const
{
    return file_paths().size();
}

std::vector<std::string> const& FileInputOptions::file_paths() const
{
    #pragma omp critical(GRENEDALF_FILE_INPUT_PATHS)
    {
        if( resolved_paths_.empty() ) {
            using namespace genesis::utils;
            for( auto const& path : raw_paths_ ) {
                if( is_file( path ) ) {

                    resolved_paths_.push_back( path );

                } else if( is_dir( path ) ) {

                    // Get all files in dir.
                    auto list = dir_list_files( path, true, ".*\\." + file_ext_ + "$" );
                    for( auto const& jplace : list ) {
                        resolved_paths_.push_back( jplace );
                    }

                } else {
                    // throw std::runtime_error( "Not a valid file or directory: " + path );
                    throw CLI::ValidationError(
                        "--" + file_type_ + "-path", "Not a valid file or directory: " + path
                    );
                }
            }

            // If required, we actually need files!
            if( resolved_paths_.empty() && option_->get_required() ) {
                throw CLI::ValidationError(
                    "--" + file_type_ + "-path", "No files found."
                );
            }

            // We sort them to get reproducible order.
            std::sort( resolved_paths_.begin(), resolved_paths_.end() );
        }
    }

    return resolved_paths_;
}

std::string const& FileInputOptions::file_path( size_t index ) const
{
    auto const& files = file_paths();
    if( index >= files.size() ) {
        throw std::runtime_error( "Invalid file index." );
    }
    return files[ index ];
}

std::vector<std::string> const& FileInputOptions::raw_file_paths() const
{
    return raw_paths_;
}

std::vector<std::string> FileInputOptions::base_file_names() const
{
    using namespace genesis::utils;

    auto paths = file_paths();
    for( auto& path : paths ) {
        auto fn = file_basename( path );
        if( ends_with( fn, ".gz" ) ) {
            fn.erase( fn.size() - 3 );
        }
        path = file_filename( fn );
    }
    return paths;
}

std::string FileInputOptions::base_file_name( size_t index ) const
{
    using namespace genesis::utils;
    auto fn = file_basename( file_path( index ));
    if( ends_with( fn, ".gz" ) ) {
        fn.erase( fn.size() - 3 );
    }
    return file_filename( fn );
}

void FileInputOptions::print() const
{
    std::string type = file_type_;
    if( ! type.empty() ) {
        type = " " + type;
    }

    // Print list of files, depending on verbosity.
    auto const& files = file_paths();
    LOG_MSG1 << "Found " << files.size() << type << " file"
             << ( files.size() > 1 ? "s" : "" );
    // LOG_MSG2 << genesis::utils::join( base_file_names(), ", " );

    // auto const& files = file_paths();
    // if( global_options.verbosity() == 0 ) {
    //     return;
    // } else if( global_options.verbosity() == 1 ) {
    //     std::cout << "Found " << files.size() << type << " file";
    //     std::cout << ( files.size() > 1 ? "s" : "" ) << ".\n";
    // } else if( global_options.verbosity() == 2 ) {
    //     std::cout << "Found " << files.size() << type << " file";
    //     std::cout << ( files.size() > 1 ? "s" : "" ) << ": ";
    //     for( auto const& file : files ) {
    //         if( &file != &files[0] ) {
    //             std::cout << ", ";
    //         }
    //         std::cout << genesis::utils::file_basename( file );
    //     }
    //     std::cout << "\n";
    // } else {
    //     std::cout << "Found " << files.size() << type << " file";
    //     std::cout << ( files.size() > 1 ? "s" : "" ) << ":\n";
    //
    //     for( auto const& file : files ) {
    //         std::string rp;
    //         try{
    //             rp = genesis::utils::real_path( file );
    //         } catch(...) {
    //             rp = file;
    //         }
    //         std::cout << "  - " << rp << "\n";
    //     }
    // }
}
