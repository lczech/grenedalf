grenedalf is a collection of commands for working with population genetic data.

Many commands in grenedalf are re-implementations of commands of the
[PoPoolation](https://sourceforge.net/projects/popoolation/) and
[PoPoolation2](https://sourceforge.net/projects/popoolation2/) tools.
However, being written in C++, grenedalf is much faster and needs less memory.

# General usage

  * [Setup](../wiki/Setup)
  * [Overview](../wiki/Overview)
  * [Input](../wiki/Input)
  * [Filtering](../wiki/Filtering)
  * [Windowing](../wiki/Windowing)
  * [Output](../wiki/Output)
  * [Example](../wiki/Example)

# Command Line Interface

grenedalf is used via its command line interface, with commands for each task.
The commands have the general structure:
<!-- grenedalf <module> <subcommand> <options> -->

    grenedalf <command> [options]

Use `grenedalf --help` to get an overview of all commands, and use `grenedalf <command> --help` for the list of options of each command. The pages below also list these options, along with additional information on each command.
