#!/usr/bin/env python3

import sys, os
import pandas as pd
import tqdm

os.makedirs("mpileup", exist_ok=True)
os.makedirs("sync", exist_ok=True)

chromosome = "chr1"
ref_allele = 'A'
derived_allele = 'C'

# ------------------------------------------------------------
#     generate pileup
# ------------------------------------------------------------

def generate_mpileup(file_id, coverage, derived_count, window_size):
    output_file = f"mpileup/{file_id}.pileup"
    # print(f"Generating mpileup file: {output_file}")

    with open(output_file, 'w') as f:
        for pos in range(1, window_size + 1):
            if pos == 1:
                cov1 = coverage - derived_count
                cov2 = derived_count
            else:
                cov1 = coverage
                cov2 = 0
            bases = ref_allele * cov1 + derived_allele * cov2
            qualities = 'I' * coverage
            f.write(f"{chromosome}\t{pos}\t{ref_allele}\t{coverage}\t{bases}\t{qualities}\n")

def generate_mpileup_sample(coverage, derived_count, window_size):
    file_id = f"mpileup-cov_{coverage}-der_{derived_count}-win_{window_size}"
    generate_mpileup(file_id, coverage, derived_count, window_size)

# ------------------------------------------------------------
#     generate sync
# ------------------------------------------------------------

def write_sync_sample( filehandle, coverage, derived_count ):
    a_cnt = coverage - derived_count
    t_cnt = derived_count
    filehandle.write(f"\t{a_cnt}:{t_cnt}:0:0:0:0")

def generate_sync(file_id, coverage_1, coverage_2, derived_count_1, derived_count_2, window_size):
    output_file = f"sync/{file_id}.sync"
    # print(f"Generating sync file: {output_file}")

    with open(output_file, 'w') as f:
        for pos in range(1, window_size + 1):
            f.write(f"{chromosome}\t{pos}\t{ref_allele}")
            if pos == 1:
                write_sync_sample( f, coverage_1, derived_count_1 )
                write_sync_sample( f, coverage_2, derived_count_2 )
            else:
                write_sync_sample( f, coverage_1, 0 )
                write_sync_sample( f, coverage_2, 0 )
            f.write("\n")

def generate_sync_sample(coverage_1, coverage_2, derived_count_1, derived_count_2, window_size):
    file_id = f"sync-cov1_{coverage_1}-cov2_{coverage_2}-der1_{derived_count_1}-der2_{derived_count_2}-win_{window_size}"
    generate_sync(file_id, coverage_1, coverage_2, derived_count_1, derived_count_2, window_size)

# ------------------------------------------------------------
#     read table
# ------------------------------------------------------------

def create_from_df(df):
    # Iterate over each row
    for _, row in df.iterrows():
        generate_mpileup_sample(
            int(row['coverage_1']),
            int(row['derived_count_1']),
            int(row['window_size'])
        )
        generate_mpileup_sample(
            int(row['coverage_2']),
            int(row['derived_count_2']),
            int(row['window_size'])
        )
        generate_sync_sample(
            int(row['coverage_1']),
            int(row['coverage_2']),
            int(row['derived_count_1']),
            int(row['derived_count_2']),
            int(row['window_size'])
        )

def create_from_param_space(param_space):
    print("Creating files")
    for params in tqdm.tqdm(param_space):
        n1, n2, c1, c2, k1, k2, w = params

        generate_mpileup_sample(
            int(c1),
            int(k1),
            int(w)
        )
        generate_mpileup_sample(
            int(c2),
            int(k2),
            int(w)
        )
        generate_sync_sample(
            int(c1),
            int(c2),
            int(k1),
            int(k2),
            int(w)
        )

def main():
    # Read the tab-delimited table using pandas
    infile = "independent_check_statistics.tsv"
    df = pd.read_csv(infile, delimiter='\t')
    # print(df)
    create_from_df(df)

if __name__ == '__main__':
    main()
