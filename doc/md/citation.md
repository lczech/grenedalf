# Description

The auxiliary command prints references that need to be cited when using grenedalf. By default, it simply prints the grenedalf reference itself.

It can however also be used to print other references used by commands in grenedalf. Each reference has a citation key (similar to how BibTeX uses them). If such `keys` are listed after the command, only the specified references are printed:

    grenedalf citation Czech2021-grenedalf

The `--list` flag prints a list of all available citation keys.
Furthermore, the `--all` flag prints out all available references that are used in grenedalf commands.
Using the `--format` option, the output format can be chosen as needed.
