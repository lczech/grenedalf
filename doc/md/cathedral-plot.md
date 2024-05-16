# Description

Create cathedral plots. This is the second step after [`fst-cathedral`](../wiki/Subcommand:-fst-cathedral), and turns the matrix computed there into the actual plots. We split this into two commands
for efficiency, so that it is faster to iterate over different color schemes and other plotting settings.

The command takes either the `csv` or the `json` files as input (and infers the respective other).
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
