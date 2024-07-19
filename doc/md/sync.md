# Description

Create a `sync` file from the given inputs. The file format is described in [Input: sync](../wiki/Input#sync).

This format is the fastest to parse in grenedalf, so we recommend to convert any complex input formats such as `sam`/`bam` files to `sync` once, if you plan to run multiple analyses on these files. For that purpose, it is best to not apply any numerical filters, as they can always be provided later on for the actual analysis, so that it's easier to adjust them to the needs then.

Alternatively, the created `sync` files can serve as a simple format that can be parsed externally with not too much effort, and hence is a simple way to obtain base counts for all nucleotides in `ACGT` for the input data.

See the [Settings](../wiki/Subcommand:-sync#settings) for options specific to this command.
