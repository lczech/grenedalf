#!/bin/bash

# Remove the generated data files
rm -r mpileup
rm -r sync
rm -r poolsizes

# Remove the run files
rm -r grenedalf/diversity
rm -r grenedalf/fst
rm -r grenedalf/logs
rm -r popoolation/diversity
rm -r popoolation/fst
rm -r popoolation/logs

# Remove the test result files
# rm -r test_results.csv
rm -r figures_*

# Also remove unnecessary py stuff
rm -r __pycache__
