# Description

Compute Pool-seq corrected estimators of nucleotide diversity. In short, we correct for the two sources of noise in pool sequencing: First, when selecting individuals from a population, we have sampling bias due to the finite number of individuals in the pool; second, when selecting reads in the sequencing, a second level of sampling bias is introduced due to the finite coverage per position.

We compute three common estimators that correct for these noises:

  * Theta Pi
  * Theta Watterson
  * Tajima's D

For details on the statistics and the derivation of the estimators, please refer to our [publication](https://doi.org/10.1093/bioinformatics/btae508) and in particular our [equations document](https://github.com/lczech/pool-seq-pop-gen-stats).

As described there, we recommend to be careful when numerically interpreting Tajima's D for Pool seq data, as the estimator is not well suited to work well with the noises introduced by the nested sampling process.

The estimators use the following filter settings internally as well as part of their corrections: `--filter-sample-min-count`, `--filter-sample-min-read-depth`, and `--filter-sample-max-read-depth`. See the equations document above for details on how these parameters influence the estimators. Note that the `total` filters are applied for filtering here as well (as far as provided), but have no influence on the estimators themselves. We offer them here for convenience only.

The `--filter-sample-min-count` has to be set to exactly 2 when computing Tajima's D. This is a consequence of the equations, which only work out when exactly filtering out singleton counts, according to Kofler et al. Furthermore, the correction terms for Theta Pi and Theta Watterson make use of the `--filter-sample-min-read-depth`, by computing a sum over values in that range. Due to the exact way that this correction works, it is recommended to use a minimum read depth of at least twice the minimum per-base count:

    --filter-sample-min-count 2
    --filter-sample-min-read-depth 4

Lastly, the values per window are computed using the window averaging as described [here](../wiki/Windowing#window-averaging-policy), in order to scale the results per base pair. See [Output](../wiki/Output) for a description of the columns of the produced table.
