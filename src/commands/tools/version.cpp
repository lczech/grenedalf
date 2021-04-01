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

#include "commands/tools/version.hpp"
#include "options/global.hpp"
#include "tools/references.hpp"
#include "tools/version.hpp"

#include "genesis/utils/core/options.hpp"
#include "genesis/utils/io/base64.hpp"

#include <algorithm>
#include <chrono>
#include <random>
#include <vector>

// =================================================================================================
//      Setup
// =================================================================================================

// Local forward declaraction.
void run_ee();

void setup_version( CLI::App& app )
{
    // Create the options and subcommand objects.
    auto options = std::make_shared<VersionOptions>();
    auto sub_version = app.add_subcommand(
        "version",
        "Extended version information about grenedalf."
    );

    // Set the run function as callback to be called when this subcommand is issued.
    // Hand over the options by copy, so that their shared ptr stays alive in the lambda.
    sub_version->callback( [ options ]() {
        run_version( *options );
    });

    // Same for the ee command, but this time without options at all.
    auto sub_ee = app.add_subcommand(
        "quote",
        "Quotes from Gandalf the Grey and Ralph Wiggum. Perfectly balanced, as all things should be."
    )->group( "" );
    sub_ee->callback( []() {
        run_ee();
    });
}

// =================================================================================================
//      Run
// =================================================================================================

void run_version( VersionOptions const& options )
{
    (void) options;

    LOG_BOLD << grenedalf_header();
    LOG_BOLD;
    LOG_BOLD << "grenedalf version: " << grenedalf_version();
    LOG_BOLD;
    LOG_BOLD << genesis::utils::Options::get().info_compile_time();
    LOG_BOLD;
    LOG_BOLD << "For citation information, call  `grenedalf tools citation`";
    LOG_BOLD << "For license information, call  `grenedalf tools license`";
    LOG_BOLD;
    LOG_BOLD << grenedalf_title();
}

void run_ee()
{
    std::vector<std::string> qs = {
        "R2FuZGFsZj8gWWVzLCB0aGF0IHdhcyB3aGF0IHRoZXkgdXNlZCB0byBjYWxsIG1lLi4uIEdhbmRhbGYgdGhlIEdyZX"
        "kuIFRoYXQgd2FzIG15IG5hbWUu", "SXQgaXMgbm90IGRlc3BhaXIsIGZvciBkZXNwYWlyIGlzIG9ubHkgZm9yIHRo"
        "b3NlIHdobyBzZWUgdGhlIGVuZCBiZXlvbmQgYWxsIGRvdWJ0LiBXZSBkbyBub3Qu", "SSB3YXMgdGFsa2luZyBhbG"
        "91ZCB0byBteXNlbGYuIEEgaGFiaXQgb2YgdGhlIG9sZDogdGhleSBjaG9vc2UgdGhlIHdpc2VzdCBwZXJzb24gcHJl"
        "c2VudCB0byBzcGVhayB0by4=", "SSBhbSBhIHNlcnZhbnQgb2YgdGhlIFNlY3JldCBGaXJlLCB3aWVsZGVyIG9mIH"
        "RoZSBmbGFtZSBvZiBBbm9yLiBZb3UgY2Fubm90IHBhc3MuVGhlIGRhcmsgZmlyZSB3aWxsIG5vdCBhdmFpbCB5b3Us"
        "IGZsYW1lIG9mIFVkdW4uIEdvIGJhY2sgdG8gdGhlIFNoYWRvdyEgWW91IGNhbm5vdCBwYXNzLg==", "U29tZXRoaW"
        "5nIGhhcyBjcmVwdCwgb3IgaGFzIGJlZW4gZHJpdmVuIG91dCBvZiBkYXJrIHdhdGVycyB1bmRlciB0aGUgbW91bnRh"
        "aW5zLiBUaGVyZSBhcmUgb2xkZXIgYW5kIGZvdWxlciB0aGluZ3MgdGhhbiBPcmNzIGluIHRoZSBkZWVwIHBsYWNlcy"
        "BvZiB0aGUgd29ybGQu", "SSBtdXN0IHJlc3QgaGVyZSBhIG1vbWVudCwgZXZlbiBpZiBhbGwgdGhlIG9yY3MgZXZl"
        "ciBzcGF3bmVkIGFyZSBhZnRlciB1cy4=", "SXQgaXMgYSBjb21mb3J0IG5vdCB0byBiZSBtaXN0YWtlbiBhdCBhbG"
        "wgcG9pbnRzLiBEbyBJIG5vdCBrbm93IGl0IG9ubHkgdG9vIHdlbGwh", "SGUgdGhhdCBicmVha3MgYSB0aGluZyB0"
        "byBmaW5kIG91dCB3aGF0IGl0IGlzIGhhcyBsZWZ0IHRoZSBwYXRoIG9mIHdpc2RvbS4=", "SXQgaXMgbm90IG91ci"
        "BwYXJ0IHRvIG1hc3RlciBhbGwgdGhlIHRpZGVzIG9mIHRoZSB3b3JsZCwgYnV0IHRvIGRvIHdoYXQgaXMgaW4gdXMg"
        "Zm9yIHRoZSBzdWNjb3VyIG9mIHRob3NlIHllYXJzIHdoZXJlaW4gd2UgYXJlIHNldCwgdXByb290aW5nIHRoZSBldm"
        "lsIGluIHRoZSBmaWVsZHMgdGhhdCB3ZSBrbm93LCBzbyB0aGF0IHRob3NlIHdobyBsaXZlIGFmdGVyIG1heSBoYXZl"
        "IGNsZWFuIGVhcnRoIHRvIHRpbGwuIFdoYXQgd2VhdGhlciB0aGV5IHNoYWxsIGhhdmUgaXMgbm90IG91cnMgdG8gcn"
        "VsZS4=", "VGhlbiBkYXJrbmVzcyB0b29rIG1lLCBhbmQgSSBzdHJheWVkIG91dCBvZiB0aG91Z2h0IGFuZCB0aW1l"
        "LCBhbmQgSSB3YW5kZXJlZCBmYXIgb24gcm9hZHMgdGhhdCBJIHdpbGwgbm90IHRlbGwuIE5ha2VkIEkgd2FzIHNlbn"
        "QgYmFjayDigJMgZm9yIGEgYnJpZWYgdGltZSwgdW50aWwgbXkgdGFzayBpcyBkb25lLiA=", "SXQgaXMgbm90IG91"
        "ciBwYXJ0IGhlcmUgdG8gdGFrZSB0aG91Z2h0IG9ubHkgZm9yIGEgc2Vhc29uLCBvciBmb3IgYSBmZXcgbGl2ZXMgb2"
        "YgTWVuLCBvciBmb3IgYSBwYXNzaW5nIGFnZSBvZiB0aGUgd29ybGQuIFdlIHNob3VsZCBzZWVrIGEgZmluYWwgZW5k"
        "IG9mIHRoaXMgbWVuYWNlLCBldmVuIGlmIHdlIGRvIG5vdCBob3BlIHRvIG1ha2Ugb25lLg==", "VGhlIHRyZWFjaG"
        "Vyb3VzIGFyZSBldmVyIGRpc3RydXN0ZnVsLg==", "T25seSBhIHNtYWxsIHBhcnQgaXMgcGxheWVkIGluIGdyZWF0"
        "IGRlZWRzIGJ5IGFueSBoZXJvLg==", "WWV0IGEgdHJlYWNoZXJvdXMgd2VhcG9uIGlzIGV2ZXIgYSBkYW5nZXIgdG"
        "8gdGhlIGhhbmQu", "SSB3aWxsIG5vdCBzYXk6IGRvIG5vdCB3ZWVwOyBmb3Igbm90IGFsbCB0ZWFycyBhcmUgYW4g"
        "ZXZpbC4=", "SSdtIGEgU3RhciBXYXJzLg==", "VGhpcyBpcyBteSBzd2luZyBzZXQuIFRoaXMgaXMgbXkgc2FuZG"
        "JveCwgSSdtIG5vdCBhbGxvd2VkIHRvIGdvIGluIHRoZSBkZWVwIGVuZC4=", "U2xvdyBkb3duLCBCYXJ0ISBNeSBs"
        "ZWdzIGRvbuKAmXQga25vdyBob3cgdG8gYmUgYXMgbG9uZyBhcyB5b3Vycy4=", "SWYgTW9tbXkncyBwdXJzZSBkaW"
        "RuJ3QgYmVsb25nIGluIHRoZSBtaWNyb3dhdmUsIHdoeSBkaWQgaXQgZml0Pw==", "SGkgU3VwZXIgTmludGVuZG8g"
        "Q2hhbG1lcnMsIEknbSBsZWFybmVkaW5nLg==", "SSdtIGEgdW5pdGFyZC4=", "QW5kIHdoZW4gdGhlIGRvY3Rvci"
        "BzYWlkIEkgZGlkbid0IGhhdmUgd29ybXMgYW55bW9yZSwgdGhhdCB3YXMgdGhlIGhhcHBpZXN0IGRheSBvZiBteSBs"
        "aWZlLg==", "UGxhbnQgaXQgYW5kIHlvdSdsbCBncm93IGEgbmV3IFJhbHBoLg==", "WW91ciBleWVzIG5lZWQgZG"
        "lhcGVycy4=", "SSdtIGluIGRhbmdlci4=", "TGllcyBhcmUgbGlrZSBzdGFycywgdGhleSBhbHdheXMgY29tZSBv"
        "dXQuLi4gSSBoYXZlIGZpdmUgZmFjZSBob2xlcy4=", "Q2FuIHlvdSBvcGVuIG15IG1pbGssIE1vbW15Lg==", "SS"
        "B3YXMgZG9uZSBiZWZvcmUgd2UgY2FtZSBpbi4=", "SSB3YW50IHRvIGJlIGEgdHJpYW5nbGUu", "TWUgZmFpbCBF"
        "bmdsaXNoLi4udGhhdCdzIHVucG9zc2libGUu"
    };

    auto const seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine engine( seed );
    std::shuffle( qs.begin(), qs.end(), engine );
    LOG_BOLD << genesis::utils::base64_decode_string( qs[0] );
}
