![grenedalf](/doc/logo/grenedalf.png?raw=true "grenedalf")

Features
-------------------

grenedalf is a collection of commands for working with population genetic data,
in particular from pool sequencing.
Its main focus are statistical analyses such as Tajima's D and Fst, following the approaches of
[PoPoolation](https://sourceforge.net/projects/popoolation/) and
[PoPoolation2](https://sourceforge.net/projects/popoolation2/).

Setup
-------------------

To run grenedalf on your machine, simply get it, and build it:

~~~.sh
git clone --recursive https://github.com/lczech/grenedalf.git
cd grenedalf
make
~~~

You can also use the green "Code" button above or
[click here](https://github.com/lczech/grenedalf/archive/master.zip) to download the source as a zip
archive. Unpack, and call `make` in the main directory to build everything.

Requirements:

 *  [Make](https://www.gnu.org/software/make/) and [CMake](https://cmake.org/) 2.8.7 or higher.
 *  A fairly up-to-date C++11 compiler, e.g.,
    [clang++ 3.6](http://clang.llvm.org/) or [GCC 4.9](https://gcc.gnu.org/), or higher.

After building, the executable is stored in the `bin` directory, and used as follows.

Usage and Documentation
-------------------

grenedalf is used via its command line interface, with commands for each task.
The commands have the general structure:
<!-- grenedalf <module> <subcommand> <options> -->

    grenedalf <command> <options>

See [**the Wiki pages**](https://github.com/lczech/grenedalf/wiki)
for the full list of all subcommands and their documentation.

<!-- # grenedalf
Genome Analyses of Differential Allele Frequencies -->
