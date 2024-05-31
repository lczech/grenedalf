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

#include "tools/references.hpp"

#include "genesis/utils/text/string.hpp"

#include <algorithm>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

// =================================================================================================
//      Citations
// =================================================================================================

// -------------------------------------------------------------------------
//     List of References
// -------------------------------------------------------------------------

struct Citation
{
    struct Author
    {
        std::string first;
        std::string last;
    };

    std::vector<Author> authors;
    std::string title;
    std::string journal;
    std::string volume;
    std::string issue;
    std::string year;
    std::string doi;
};

static std::unordered_map<std::string, Citation> citations_ = {
    { "Czech2023-grenedalf", {
        {{ "Lucas", "Czech" }, { "Jeffrey", "Spence" }, { "Moises", "Exposito-Alonso" }},
        "grenedalf: population genetic statistics for the next generation of pool sequencing",
        "arXiv",
        "",
        "",
        "2023",
        "10.48550/arXiv.2306.11622"
    }},
    { "Kofler2011-PoPoolation", {
        {
            { "Robert", "Kofler" }, { "Pablo", "Orozco-terWengel" }, { "Nicola", "De Maio" },
            { "Ram Vinay", "Pandey" }, { "Viola", "Nolte" }, { "Andreas", "Futschik" },
            { "Carolin", "Kosiol" }, { "Christian", "Schlötterer" }
        },
        "PoPoolation: A Toolbox for Population Genetic Analysis of Next Generation Sequencing Data from Pooled Individuals",
        "PLoS ONE",
        "6",
        "1",
        "2011",
        "10.1371/journal.pone.0015925"
    }},
    { "Kofler2011-PoPoolation2", {
        {{ "Robert", "Kofler" }, { "Ram Vinay", "Pandey" }, { "Christian", "Schlötterer" }},
        "PoPoolation2: identifying differentiation between populations using sequencing of pooled DNA samples (Pool-Seq)",
        "Bioinformatics",
        "27",
        "24",
        "2011",
        "10.1093/bioinformatics/btr589"
    }},
};

// -------------------------------------------------------------------------
//     Helper Functions
// -------------------------------------------------------------------------

void check_citation_duplicates( std::vector<std::string> keys )
{
    std::sort( keys.begin(), keys.end() );
    auto uniq = std::unique( keys.begin(), keys.end() );
    if( uniq != keys.end() ) {
        throw std::runtime_error( "Duplicate citation keys: " + (*uniq) );
    }
}

Citation const& get_citation( std::string const& key )
{
    if( citations_.count(key) == 0 ) {
        throw std::runtime_error( "Invalid citation key: " + key );
    }
    auto const& entry = citations_.at(key);

    if(
        entry.authors.empty() || entry.title.empty() ||
        entry.journal.empty() || entry.year.empty()  || entry.doi.empty()
    ) {
        throw std::runtime_error( "Citation is missing some information: " + key );
    }
    for( auto const& author : entry.authors ) {
        if( author.first.empty() || author.last.empty() ) {
            throw std::runtime_error( "Citation is missing author information: " + key );
        }
    }

    return entry;
}

std::string cite_authors( Citation const& entry, bool first_last, std::string const& delim )
{
    std::vector<std::string> authors;
    for( auto const& author : entry.authors ) {
        if( first_last ) {
            authors.push_back( author.first + " " + author.last );
        } else {
            authors.push_back( author.last + ", " + author.first );
        }
    }
    return genesis::utils::join( authors, delim );
}

void check_all_citations()
{
    // Simply retrieve all citations, which triggers their validity checks.
    for( auto const& entry : citations_ ) {
        (void) get_citation( entry.first );
    }
}

void check_citation( std::string const& key )
{
    // Simply retrieve the citation, which triggers its validity checks.
    (void) get_citation( key );
}

void check_citations( std::vector<std::string> const& keys )
{
    // Simply retrieve all citations, which triggers their validity checks.
    for( auto const& entry : keys ) {
        (void) get_citation( entry );
    }

    check_citation_duplicates( keys );
}

std::vector<std::string> get_all_citation_keys()
{
    std::vector<std::string> result;
    for( auto const& entry : citations_ ) {
        result.push_back( entry.first );
    }
    std::sort( result.begin(), result.end() );
    return result;
}

// -------------------------------------------------------------------------
//     Run Functions
// -------------------------------------------------------------------------

std::string cite_bibtex( std::string const& key )
{
    auto const& entry = get_citation( key );

    std::stringstream ss;
    ss << "@article{" << key << ",\n";
    ss << "    author = {" << cite_authors( entry, false, " and " ) << "},\n";
    ss << "    title = {{" << entry.title << "}},\n";
    ss << "    journal = {" << entry.journal << "},\n";
    ss << "    year = {" << entry.year << "},\n";
    if( ! entry.volume.empty() ) {
        ss << "    volume = {" << entry.volume << "},\n";
    }
    if( ! entry.issue.empty() ) {
        ss << "    number = {" << entry.issue << "},\n";
    }
    ss << "    doi = {" << entry.doi << "}\n";
    ss << "}\n";
    return ss.str();
}

std::string cite_markdown( std::string const& key, bool with_quote_block, bool with_key )
{
    auto const& entry = get_citation( key );
    std::string const in = ( with_quote_block ? "> " : "" );

    std::stringstream ss;
    if( with_key ) {
        ss << key << ":\n";
    }
    ss << in << cite_authors( entry, true, ", " ) << ".\n";
    ss << in << "**" << entry.title << ".**\n";
    ss << in << "*" << entry.journal << "*";
    if( ! entry.volume.empty() ) {
        ss << ", vol. " << entry.volume;
    }
    if( ! entry.issue.empty() ) {
        ss << ", no. " << entry.issue;
    }
    ss << ", " << entry.year << ".\n";
    ss << in << "doi:[" << entry.doi << "](https://doi.org/" << entry.doi << ")\n";
    return ss.str();
}

std::string cite_bibtex( std::vector<std::string> const& keys )
{
    check_citation_duplicates( keys );

    std::string result;
    for( auto const& key : keys ) {
        if( &key != &keys[0] ) {
            result += "\n";
        }
        result += cite_bibtex( key );
    }
    return result;
}

std::string cite_markdown( std::vector<std::string> const& keys, bool with_quote_block, bool with_key )
{
    check_citation_duplicates( keys );

    std::string result;
    for( auto const& key : keys ) {
        if( &key != &keys[0] ) {
            result += "\n";
        }
        result += cite_markdown( key, with_quote_block, with_key );
    }
    return result;
}
