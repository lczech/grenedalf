# Description

Compute Pool-seq corrected FST. In short, we correct for the two sources of noise in pool sequencing: First, when selecting individuals from a population, we have sampling bias due to the finite number of individuals in the pool; second, when selecting reads in the sequencing, a second level of sampling bias is introduced due to the finite coverage per position.

We offer several flavors of FST to correct for those:

  * Pool-seq corrected variants of FST following the definitions of Nei and Hudson, respectively.
  * A re-implementation of the "Kofler" estimator of PoPoolation2.
  * A re-implementation of the Karlsson estimator as described in PoPoolation2.

For details on the statistics and the derivation of the estimators, please refer to our [publication](https://arxiv.org/abs/2306.11622) and in particular our [equations document](https://github.com/lczech/pool-seq-pop-gen-stats).

We recommend to use one of our two estimators, as they are correctly account for the noises of Pool-seq.
In simplified form (without the Pool-seq correction terms), the two variants are:

    FST_Nei = 1 - pi within / pi total
    FST_Hudson = 1 - pi within / pi between

Our recommendation is to use the Hudson estimator with a cutoff for low frequencies (either via `--filter-total-snp-min-count` or `--filter-total-snp-min-frequency`), see for instance [10.1101/gr.154831.113](https://doi.org/10.1101/gr.154831.113) for a rationale.
<!-- > Bhatia, G. et al. "Estimating and interpreting FST: The impact of rare variants". **Genome Research**, 2013. https://doi.org/10.1101/gr.154831.113 -->

The values per window are computed using the window averaging as described [here](../wiki/Windowing#window-averaging-policy), in order to scale the results per base pair.

<!-- we expect NaN if poolsize is 1. -->

# Output files

The output of the command is a table with values for each pair of samples, where each row contains the FST values in a window (according to `--window-type`) along the genome. See [Output](../wiki/Output) for a description of the columns.

## Pi within / between / total

When using one of our unbiased estimators (follow the FST definition of either Nei or Hudson), additional tables can be produced by specifying `--write-pi-tables`. They are named `pi-within`, `pi-between`, `pi-total`, respectively. These tables contain the same columns as the main output table, but instead of the final FST value, they contain the intermediate values of the involved estimates of pi that are used in the Nei and Hudson estimators.

## Whole genome output

When `--window-type genome` is set, that is, when we compute a single whole-genome FST value per pair of samples, we produce additionally output tables for convenience. As in that case we only have a single value per pair of samples, we can use a more concise format that is typically easier to use for post-processing.

First, we produce `fst-list`, which contains the pair of sample names in the first two columns, and their FST in the third. This is the same as the main output table, but in long instead of wide format. As this table contains each sample pair in a row, it is often easier to work with than a wide column when having many sample pairs.

Second, if no `--comparand` or `--comparand-list` is set (in other words, if whole-genome all-to-all FST is being computed), we additionally output a file `fst-matrix` that contains all pairwise whole-genome FST values for the input samples in a matrix.
