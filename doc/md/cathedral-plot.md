# Description

Create cathedral plots. This is the second step after [`fst-cathedral`](../wiki/Subcommand:-fst-cathedral), and turns the matrix computed there into the actual plots. We split this into two commands
for efficiency, so that it is faster to iterate over different color schemes and other plotting settings.

The command takes either the `csv` or the `json` files produced by `fst-cathedral` as input (and infers the respective other file, both have to be in the same directory, with the same base name).
It then colors the value matrix according to the provided color map settings, and stores the result as
a `bmp` bitmap picture. Furthermore, it creates an `svg` file that additionally contains axes, a
legend for the color map, and a title, and can be edited and refined later with any vector graphics
program.

![FST Cathedral Plot.](https://github.com/lczech/grenedalf/blob/master/doc/png/fst_cathedral.png?raw=true)

See [`fst-cathedral`](../wiki/Subcommand:-fst-cathedral) for details on this plot.

Note that the pool-sequencing corrected estimator of FST applies a correction term that can yield
FST values below zero. This is expected, and a consequence of correcting for the statistical noise.
For details, see our [equations document](https://github.com/lczech/pool-seq-pop-gen-stats/releases).
In order to not have these artifacts influence the plot, and to create consistency in the plots,
we recommend to clip negative values to zero, by providing `--min-value 0 --clip-under`. The first
option limits the scale to non-negative values, and the second option makes sure that the negative
values are clipped to be 0, instead of being highlighted in the `--under-color`.

Similarly, it might be beneficial to use `--max-value X --clip-over` with some reasonable maximum
value `X`, if multiple plots are created that need to be compared to each other. That way, all
plots will have the same scale, and hence have comparable color values.

Lastly, at the moment, we only have implemented cathedral plots for FST. They are however also possible
for any other window-based statistic, such as the diversity metrics. If this is something that you
are interested in, please open an [issue](https://github.com/lczech/grenedalf/issues) to tell us.

# Colors

For the options of this command, the single colors and the main gradient can be specified as follows.

## Single Colors

Single colors can be specified

 * by name, as one of the 140 [web colors](https://en.wikipedia.org/wiki/Web_colors), that is, the basic 16 html color names and the extended 124 X11 color names. This is case-independent and insensitive to white spaces.
 * by name, as one of the 954 [xckd colors](https://xkcd.com/color/rgb/), again case- and white-space-insensitive.
  * by hex code in the format `#RRGGBB` or `#RRGGBBAA` (with alpha, which might be useful when producing svg files), using hexadecimal coding for each of the red, green, and blue values, case insensitive. For example, use `#000000` for black and `#ffffff` for white. Note that `#` also happens to denote the start of a comment in command lines; hence, you probably need to put this in quotation marks.

A typical color specification might hence look like this: `--under-color "#ff00ff"` or `--mask-color orange`.

## Lists and Gradients of Colors

Gradients and lists of colors can be specified as

  * a comma-separated list of colors following the above specifications for single colors (this list can either be provided in a file with one color per line, or directly as a string on the command line), or
  * as one of the following named color lists/gradients:

    ![Color lists in grenedalf.](https://github.com/lczech/genesis/raw/master/doc/png/utils/color_lists.png)

Depending on context, not all of these lists might be well suited; it does for example not make much sense to use a (categorical) qualitative color list as a (continuous) gradient.

When specifying individual colors to build a custom gradient, the specified colors are evenly spaced out across the range of values, and then linearly interpolated to create the gradient. For example, a gradient from black to red to yellow could be specified as `--color-list "#000000,#ff0000,#ffff00"`.

Our internal interpolation between colors to create a gradient (currently) is done linearly in RGB color space - this does not always yield the best looking results. We hence recommend to construct a gradient with several (5 or more) intermediate colors using external tools that operate in LCH space (e.g., [this gradient generator](https://colordesigner.io/gradient-generator)), and then use these intermediate colors as input. This way, we only need to interpolate between nearby similar colors in RGB, which works/looks better than RGB interpolation between vastly different colors.
