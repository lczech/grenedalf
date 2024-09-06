[![Bioconda install](https://img.shields.io/conda/vn/bioconda/grenedalf)](https://anaconda.org/bioconda/grenedalf)
[![Downloads](https://img.shields.io/conda/dn/bioconda/grenedalf)](https://anaconda.org/bioconda/grenedalf)
[![Release](https://img.shields.io/github/v/release/lczech/grenedalf.svg)](https://github.com/lczech/grenedalf/releases)
[![CI](https://github.com/lczech/grenedalf/workflows/CI/badge.svg?branch=master)](https://github.com/lczech/grenedalf/actions)
[![Softwipe Score](https://img.shields.io/badge/softwipe-9.0/10.0-blue)](https://github.com/adrianzap/softwipe/wiki/Code-Quality-Benchmark)
![Language](https://img.shields.io/badge/language-C%2B%2B11-lightgrey.svg)
[![License](https://img.shields.io/badge/license-GPLv3-blue.svg)](http://www.gnu.org/licenses/gpl.html)
<!-- [![Platforms](https://img.shields.io/conda/pn/bioconda/grenedalf)](https://anaconda.org/bioconda/grenedalf) -->
<!-- [![DOI](https://img.shields.io/badge/doi-10.1093%2Fbioinformatics%2Fbtaa070-blue)](https://doi.org/10.1093/bioinformatics/btaa070) -->
<!-- [![Build Status](https://travis-ci.org/lczech/grenedalf.svg?branch=master)](https://travis-ci.org/lczech/grenedalf) -->

![grenedalf](/doc/logo/grenedalf.png?raw=true "grenedalf")


Features
-------------------

grenedalf is a collection of commands for working with pool sequencing population genetic data.
Its main focus are statistical analyses such as Tajima's D and Fst, as well as convenience functions around pool sequencing data.
The statistics follow the approaches of [PoPoolation](https://sourceforge.net/projects/popoolation/)
and [PoPoolation2](https://sourceforge.net/projects/popoolation2/),
as well as [poolfstat](https://cran.r-project.org/web/packages/poolfstat/index.html)
and [npstat](https://github.com/lucaferretti/npstat). However, compared to those, grenedalf
is significantly more scalable, more user friendly, and offers more settings and input file formats.


Setup
-------------------

We recommend to first try the pre-compiled binaries of the latest [Release](https://github.com/lczech/grenedalf/releases), by downloading the binary for your system from the "Assets" list below the release. Alternatively, grenedalf can be installed via [bioconda](https://anaconda.org/bioconda/grenedalf). If that does not work, it can be build from source as described [here](https://github.com/lczech/grenedalf/wiki/Setup). If you want the best performance, we recommend to build from source anyway.


Usage and Documentation
-------------------

grenedalf is used via its command line interface, with commands for each task.
The commands have the general structure:
<!-- grenedalf <module> <subcommand> <options> -->

    grenedalf <command> <options>

See [**the Wiki pages**](https://github.com/lczech/grenedalf/wiki) for the list of all commands and their documentation.

<!-- # grenedalf
Genome Analyses of Differential Allele Frequencies -->


Support and Requests
-------------------

For **questions, bug reports, and feature requests**, please [open an issue](https://github.com/lczech/grenepipe/issues). Please do not send emails with questions or requests, as others might be having them as well, and so it is better to discuss them here where they can be found.

Feedback on functionality, interface, and features is highly appreciated, as we are always open to discussing new statistics and tools that might be relevant for the community. If you need other pool sequencing statistics than the ones already offered here, please reach out - happy to explore the potential for a collaboration there! We have a lot of code infrastructure in grenedalf that takes care of all the basics, and so it would be nice to implement further statistics!


Citation
-------------------

Please find the preprint here:

> grenedalf: population genetic statistics for the next generation of pool sequencing.<br />
> Lucas Czech, Jeffrey P. Spence, Moisés Expósito-Alonso.<br />
> Bioinformatics, 2024. doi:[10.1093/bioinformatics/btae508](https://doi.org/10.1093/bioinformatics/btae508)

Each command also prints out the relevant references for that command. Then, the command [`grenedalf citation`](https://github.com/lczech/grenedalf/wiki/Subcommand:-citation) can be used to obtain details on those references.
