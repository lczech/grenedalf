<!-- [![Conda install](https://img.shields.io/conda/vn/bioconda/grenedalf)](https://anaconda.org/bioconda/grenedalf) -->
<!-- [![Downloads](https://img.shields.io/conda/dn/bioconda/grenedalf)](https://anaconda.org/bioconda/grenedalf) -->
[![Release](https://img.shields.io/github/v/release/lczech/grenedalf.svg)](https://github.com/lczech/grenedalf/releases)
[![CI](https://github.com/lczech/grenedalf/workflows/CI/badge.svg?branch=master)](https://github.com/lczech/grenedalf/actions)
[![Softwipe Score](https://img.shields.io/badge/softwipe-9.0/10.0-blue)](https://github.com/adrianzap/softwipe/wiki/Code-Quality-Benchmark)
[![License](https://img.shields.io/badge/license-GPLv3-blue.svg)](http://www.gnu.org/licenses/gpl.html)
![Language](https://img.shields.io/badge/language-C%2B%2B11-lightgrey.svg)
<!-- [![Platforms](https://img.shields.io/conda/pn/bioconda/grenedalf)](https://anaconda.org/bioconda/grenedalf) -->
<!-- [![DOI](https://img.shields.io/badge/doi-10.1093%2Fbioinformatics%2Fbtaa070-blue)](https://doi.org/10.1093/bioinformatics/btaa070) -->
<!-- [![Build Status](https://travis-ci.org/lczech/grenedalf.svg?branch=master)](https://travis-ci.org/lczech/grenedalf) -->

![grenedalf](/doc/logo/grenedalf.png?raw=true "grenedalf")

Features
-------------------

grenedalf is a collection of commands for working with population genetic data,
in particular from pool sequencing.
Its main focus are statistical analyses such as Tajima's D and Fst.
The statistics follow the approaches of [PoPoolation](https://sourceforge.net/projects/popoolation/)
and [PoPoolation2](https://sourceforge.net/projects/popoolation2/),
as well as [poolfstat](https://cran.r-project.org/web/packages/poolfstat/index.html)
and [npstat](https://github.com/lucaferretti/npstat). However, compared to those, grenedalf
is significantly more scalable, more user friendly, and offers more settings and input file formats.

**Remark:** grenedalf is quite young and in activate development. The interface will change frequently
in the near future. Everything is considered to be in beta status.
[Feedback](https://github.com/lczech/grenedalf/issues) on functionality, interface, and features
is highly appreciated!

Setup
-------------------

We recommend to first try the pre-compiled binaries of the latest [Release](https://github.com/lczech/grenedalf/releases), by downloading the binary for your system from the "Assets" list below the release. If that does not work, grenedalf can be build from source as follows. If that does not work, grenedalf can be build from source as described [here](https://github.com/lczech/grenedalf/wiki/Build).

Usage and Documentation
-------------------

grenedalf is used via its command line interface, with commands for each task.
The commands have the general structure:
<!-- grenedalf <module> <subcommand> <options> -->

    grenedalf <command> <options>

Use the `--help` flag of grenedalf or of each command for usage information.
See [**the Wiki pages**](https://github.com/lczech/grenedalf/wiki)
for the full list of all subcommands and their documentation.

<!-- # grenedalf
Genome Analyses of Differential Allele Frequencies -->

Citation
-------------------

Please find the preprint here:

> grenedalf: population genetic statistics for the next generation of pool sequencing.<br />
> Lucas Czech, Jeffrey P. Spence, Moisés Expósito-Alonso.<br />
> arXiv, 2023. [arXiv:2306.11622](https://arxiv.org/abs/2306.11622)

Each command also prints out the relevant references for that command. Then, the command [`grenedalf citation`](https://github.com/lczech/grenedalf/wiki/Subcommand:-citation) can be used to obtain details on those references.
