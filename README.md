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

Build system requirements:

  * [Make](https://www.gnu.org/software/make/), [autotools](https://en.wikipedia.org/wiki/GNU_Autotools)
    (i.e., `autoconf` and `automake`), and [CMake](https://cmake.org/) 3.1 or higher.
  * A fairly up-to-date C++11 compiler, e.g.,
    [clang++ 5](http://clang.llvm.org/) or [GCC 7](https://gcc.gnu.org/), or higher<br />
    (older versions may work as well, but we do not test them any more).
  * There are some further compression library requirements, such as `libz`, `libbz2`, `liblzma`,
    and `libdeflate`, needed for compressed file support by [htslib](https://www.htslib.org/),
    which we internally use for `sam`/`bam`/`cram` and `VCF` reading.
    Usually, they should already be available on Unix/Linux/MacOS systems though.
    We include them below anyway, just in case.

To install these, run:

  * On MacOS with [homebrew](https://brew.sh/):
    ~~~.sh
    brew install make cmake autoconf automake libdeflate coreutils
    ~~~
  * On Ubuntu and other Debian systems:
    ~~~.sh
    apt-get install build-essential cmake autoconf automake zlib1g-dev libbz2-dev liblzma-dev libtinfo5
    ~~~

This works in our [CI tests](https://github.com/lczech/grenedalf/actions);
if it does not work for you, please open an [issue](https://github.com/lczech/grenedalf/issues).

Then, to build grenedalf on your machine, run:

~~~.sh
git clone --recursive https://github.com/lczech/grenedalf.git
cd grenedalf
make
~~~

You can also use the green "Code" button above or
[click here](https://github.com/lczech/grenedalf/archive/master.zip) to download the source
as a zip. Unpack, and call `make` in the main directory to build everything.

After building, the executable is stored in the `bin` directory, and used as below.

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
