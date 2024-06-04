# Description

Compute the data for an FST cathedral plot. This command does the main part of the computation:
reading the input sample data, computing per-position FST values along the chromosomes,
and computing the full matrix of values that make the plot (i.e., that are then turned into the
pixels of the plot). For each sample pair and each chromosome, the output of the command are two files:
a `csv` file containing the per-pixel value matrix for the plot, and a `json` file describing the
plot and its meta-data.

The matrix simply contains the FST values for each pixel of the plot. Each row corresponds to one
window size, and each cell along the columns is the FST value for a window centered on the position
along the genome that the column corresponds to.

Afterwards, use the [`cathedral-plot`](../wiki/Subcommand:-cathedral-plot) command to create plots
from these files, in a plain `bmp` format, as well as an `svg` image that also contains axis,
a legend, etc. Making a cathedral plot is split into these two commands: The computation (this
command here) can take a while, but one might want to iterate quickly over several color schemes
and other plot-specific settings. Having this as an independent command hence avoids recomputing
the value matrix each time.

Note that the command processes all given pairs of samples at the same time. We split the computation per chromosome, but for each chromosome, we need to keep data memory in the order of the number of variant
positions. Hence, for a large number of samples, this can result in high memory usage. If that is the
case, simply run this command individually on subsets of sample pairs.

Also, we here do not account for missing data in the windows. That is, the pi values involved in the computation of FST are not normalized per window, and instead simply their sum is used to compute the Nei/Hudson estimators for each window. Properly accounting for missing data would entail to keep information on that during the computation per pixel as well, which we have opted not to implement at this time. As cathedral plots are meant as a visual explorative tool anyway, this should be okay. It will not make a large difference unless there is extremely uneven coverage between the samples. If you need this though, please [let us know](https://github.com/lczech/grenedalf/issues).

# FST Cathedral Plots

![FST Cathedral Plot.](https://github.com/lczech/grenedalf/blob/master/doc/png/fst_cathedral.png?raw=true)

The above shows what we call a "cathedral plot" of FST on chromosome 1 between two pool sequencing
samples of *A. thaliana*, which reveals the differentiation structure between two populations across
different scales, with darker regions having a higher FST between the two populations.

The x-axis corresponds to chromosome 1 of *A. thaliana*, which is ≈30mio bases long. The y-axis position
(rows of pixels) determines the window size to be used for that row. From top to bottom, windows get
smaller, with an exponentially decaying size. The top row corresponds to windows of the full length of
the chromosome, and towards the bottom, windows get smaller until their size matches the resolution of
the image. Here, we plotted the image with a width of 1500 pixels, so that each pixel in the bottom row
represents ≈20k positions in the genome. For a given image width, this is the highest resolution that can
be displayed; at this point, each window is exactly as wide as the pixel it corresponds to, and windows
between pixels in the row are not overlapping with each other any more. Windows are centered around
their pixel, meaning that towards the left and right of the plot (and in particular towards the top of the
image, with larger window sizes), the actual number of positions within a window is truncated.

In other words, for each pixel, the plot involves accumulating values in a window centered on that pixel,
with the width of the window determined by the row in the image. Each pixel hence corresponds to
a distinct window width and position, and shows the FST value in that window. Our implementation
computes intermediate per-SNP values once; for FST for example, we compute π within, π between,
and π total for each SNP, as explained in our [equations document](https://github.com/lczech/pool-seq-pop-gen-stats/releases).
Per window, these values are then accumulated the resulting FST is computed.

The plot hence "zooms in" from larger to smaller windows, with increasing resolution along the genome
towards the bottom. This can help to visualize the scale of a statistic across different window sizes,
showing its broad and fine structure, and allows picking appropriate window sizes for further analyses.
Note that we chose exponential decay for the window size, as it offers a good balance between
revealing both broad structures and fine details. This results in the y-axis being log-scaled.

<!-- Other functions for the window are however also possible. Of course, similar types of plots can also be useful for other per-window statistics, such as diversity. -->
