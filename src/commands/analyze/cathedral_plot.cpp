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

#include "commands/analyze/cathedral_plot.hpp"
#include "options/global.hpp"
#include "tools/cli_setup.hpp"
#include "tools/misc.hpp"

#include "genesis/population/plotting/cathedral_plot.hpp"
#include "genesis/population/plotting/genome_heatmap.hpp"
#include "genesis/utils/color/heat_map.hpp"
#include "genesis/utils/formats/bmp/writer.hpp"
#include "genesis/utils/formats/json/document.hpp"
#include "genesis/utils/formats/svg/svg.hpp"
#include "genesis/utils/text/string.hpp"

// =================================================================================================
//      Enum Mapping
// =================================================================================================

std::vector<std::pair<std::string, genesis::utils::HeatmapParameters::ColorNorm>> const color_norm_map = {
    { "linear",      genesis::utils::HeatmapParameters::ColorNorm::kLinear },
    { "logarithmic", genesis::utils::HeatmapParameters::ColorNorm::kLogarithmic },
    // { "diverging",   genesis::utils::HeatmapParameters::ColorNorm::kDiverging },
};

std::vector<std::pair<std::string, genesis::utils::HeatmapParameters::NormalizationRange>> const color_norm_range_map = {
    { "all", genesis::utils::HeatmapParameters::NormalizationRange::kAll },
    { "row", genesis::utils::HeatmapParameters::NormalizationRange::kRow },
    { "col", genesis::utils::HeatmapParameters::NormalizationRange::kCol },
};

// =================================================================================================
//      Setup
// =================================================================================================

void setup_cathedral_plot( CLI::App& app )
{
    // Create the options and subcommand objects.
    auto options = std::make_shared<CathedralPlotOptions>();
    auto sub = app.add_subcommand(
        "cathedral-plot",
        "Create a cathedral plot, using the pre-computated cathedral data."
    );

    // -------------------------------------------------------------------------
    //     Input
    // -------------------------------------------------------------------------

    // Add the file input options, separate for both types
    options->json_input.add_multi_file_input_opt_to_app( sub, "json", "json", "json", "json", false );
    options->csv_input.add_multi_file_input_opt_to_app( sub, "csv", "csv", "csv", "csv", false );
    options->json_input.option()->excludes( options->csv_input.option() );
    options->csv_input.option()->excludes( options->json_input.option() );

    // -------------------------------------------------------------------------
    //     Settings
    // -------------------------------------------------------------------------

    // Add the color map
    options->color_map.add_color_list_opt_to_app( sub, "inferno" );
    options->color_map.add_under_opt_to_app( sub );
    options->color_map.add_over_opt_to_app( sub );
    options->color_map.add_clip_opt_to_app( sub );

    // Add color norm
    options->color_norm.option = sub->add_option(
        "--color-normalization",
        options->color_norm.value,
        "To create the cathedral plot, the value of each pixel needs to be translated into a color, "
        "by mapping from the range of values into the range of the color map. This translation "
        "can be done as a simple linear transform, or logarithmic, so that low values can be "
        "distinguished with more detail."
    );
    options->color_norm.option->group( "Color" );
    options->color_norm.option->transform(
        CLI::IsMember( enum_map_keys( color_norm_map ), CLI::ignore_case )
    );

    // Add color normalization range option.
    // For now, we make this a hidden option, as it does not seem utterly useful
    // for the types of plots we want to make here.
    options->color_norm_range.option = sub->add_option(
        "--color-normalization-range",
        options->color_norm_range.value,
        "When applying the above normalization of the colors, we need to figure out the range of "
        "values that exist in the matrix, so that this range can be transformed onto the color map. "
        "This means, we need to know the min and max values in the data. This can be done for the "
        "whole matrix, in which case the color values will be comparable to each other. "
        "Alternatively, the range can be determined per row or column individually, in order to "
        "reveal more detail in rows or columns that only contain a small sub-range of the overall "
        "range of values across the matrix. Generally, the default `all` to normalize across "
        "the whole matrix makes the most sense though."
    );
    // options->color_norm_range.option->group( "Color" );
    options->color_norm_range.option->group("");
    options->color_norm_range.option->transform(
        CLI::IsMember( enum_map_keys( color_norm_range_map ), CLI::ignore_case )
    );

    // Add min max value for the heat map
    options->min_value.option = sub->add_option(
        "--min-value",
        options->min_value.value,
        "As an alternative to determining the range of values automatically, the range limits "
        "can be set explicitly. This allows for instance to cap the visualization in cases of "
        "outliers that would otherwise hide detail in the lower values. "
        "Any value that is below the min specified here will then be mapped to the `under` color, "
        "or clipped to the lowest value in the color map."
    );
    options->min_value.option->group( "Color" );
    options->max_value.option = sub->add_option(
        "--max-value",
        options->max_value.value,
        "See `--min-value`; this is the equivalent upper limit of values."
        "Any value that is above the max specified here will then be mapped to the `over` color, "
        "or be clipped to the highest value in the color map."
    );
    options->max_value.option->group( "Color" );

    // -------------------------------------------------------------------------
    //     Output
    // -------------------------------------------------------------------------

    // Output. Here, compression does not really make sense, as we want to be able to see the
    // picture files in the end anyway, so we do not add compression here.
    options->file_output.add_default_output_opts_to_app( sub );
    // options->file_output.add_file_compress_opt_to_app( sub );

    // -------------------------------------------------------------------------
    //     Callback
    // -------------------------------------------------------------------------

    // Set the run function as callback to be called when this subcommand is issued.
    // Hand over the options by copy, so that their shared ptr stays alive in the lambda.
    sub->callback( grenedalf_cli_callback(
        sub,
        {
            // Citation keys as needed
        },
        [ options ]() {
            run_cathedral_plot( *options );
        }
    ));
}

// =================================================================================================
//      Run
// =================================================================================================

void run_cathedral_plot( CathedralPlotOptions const& options )
{
    using namespace genesis::population;
    using namespace genesis::utils;

    // Check the out file target pattern, before we open any files.
    options.file_output.check_output_files_nonexistence(
        "cathedral-plot-*", std::vector<std::string>{ "bmp", "svg" }
    );

    // We only have json or csv files, never both. We just take either of them,
    // as the cathedral commands can deal with that.
    options.json_input.print();
    options.csv_input.print();
    auto const& json_paths = options.json_input.file_paths();
    auto const& csv_paths  = options.csv_input.file_paths();
    std::vector<std::string> file_paths;
    file_paths.insert(
        file_paths.end(), json_paths.begin(), json_paths.end()
    );
    file_paths.insert(
        file_paths.end(), csv_paths.begin(), csv_paths.end()
    );
    assert( file_paths.size() == json_paths.size() || file_paths.size() == csv_paths.size() );

    // Get the color map and heatmap params, which are the same for all plots.
    HeatmapParameters heatmap_params{ options.color_map.color_map() };
    heatmap_params.min_value = options.min_value.value;
    heatmap_params.max_value = options.max_value.value;
    heatmap_params.color_norm = get_enum_map_value( color_norm_map, options.color_norm.value );
    heatmap_params.normalization_range = get_enum_map_value(
        color_norm_range_map, options.color_norm_range.value
    );

    // Now we iterate each file and make a plot for it.
    for( auto const& file : file_paths ) {
        LOG_MSG << "Plotting " << file;
        auto const record = load_cathedral_plot_record_from_files( file );

        // Create the plot elements
        auto const img = make_cathedral_plot_heatmap( record, heatmap_params );
        auto const svg = make_cathedral_plot_svg( record, heatmap_params, img );

        // Write the plots to files
        auto const base_name = "cathedral-plot-" + record.plot_name + "-" + record.chromosome_name;
        BmpWriter().write( img, options.file_output.get_output_target( base_name, "bmp" ));
        svg.write( options.file_output.get_output_target( base_name, "svg" ));
    }
}
