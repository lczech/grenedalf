#!/usr/bin/env python3

# Import packages
import subprocess
import sys, os
import glob
from multiprocessing import Pool
from functools import partial
import datetime
import pandas as pd
import numpy as np
import tqdm
import warnings

# Import the other scripts in this directory
import independent_check
import create_files

# ------------------------------------------------------------
#     Parameter Space
# ------------------------------------------------------------

# File to write to
outtable = "test_results.csv"

# We here set up the param space, following the notation of independent_check.py

# For two samples (1, 2):
# n: pool sizes
# c: coverages
# f: derived frequencies
# w: window sizes
params_n1 = [10, 100]
params_n2 = [5, 50]
params_c1 = [10, 100]
params_c2 = [5, 50]
params_f1 = [0.1, 0.25, 0.5, 0.75, 0.9]
params_f2 = [0.2, 0.4, 0.6, 0.8]
params_w  = [1, 10, 100]

param_space = []
for n1 in params_n1:
    for n2 in params_n2:
        for c1 in params_c1:
            for c2 in params_c2:
                for f1 in params_f1:
                    for f2 in params_f2:
                        for w in params_w:
                            # k: derived allele count
                            k1 = round(f1*c1)
                            k2 = round(f2*c2)
                            param_space.append((n1, n2, c1, c2, k1, k2, w))

# ------------------------------------------------------------
#     Run Scipts
# ------------------------------------------------------------

# Run function, similar to the common/benchmarks.py function,
# but without all the printing and measuring, which we do not need here.
def run_script( script, args ):
    # Run a test script.
    # The function expects the script to be run, and a dict with key value pairs for the args.
    # We run bash scripts containing the tests, and they parse these args again.

    # Prepare the command to be run.
    # We cannot use normal `time` on Ubuntu, as it is a wrapper that does not have the -v option.
    # Need to run the underlying program instead.
    script=os.path.realpath(script)
    script_path = os.path.dirname(os.path.realpath(script))
    script_name = os.path.basename(os.path.realpath(script))
    script_args = " ".join([ key + "=" + str(val) for key, val in args.items() ])
    command = "/usr/bin/time -v " + script + " " + script_args + " 2>&1"

    # Start the process, capturing everything, and terminating on failure
    try:
        result = subprocess.run(
            command, cwd=script_path,
            capture_output=True, shell=True, text=True
        )
    except Exception as e:
        print("FAILED")
        print(e)
    if result.returncode != 0:
        print("FAILED")
        print(result.stdout)
        print(result.stderr)

# ------------------------------------------------------------
#     Independent Implementation
# ------------------------------------------------------------

def run_independent_implementation( params, df ):
    """
    params: the set of params to use, for simplicity as a tuple
    df: the dataframe to which to write results to,
        expected to have all columns already, and writing to the last row
        (assuemd to be present and empty at the cells that we write)
    """

    # Set up the parameters as needed
    # n: pool size
    # c: coverage
    # k: derived allele count
    # w: window size
    n1, n2, c1, c2, k1, k2, w = params
    c1s = [c1] * w
    c2s = [c2] * w
    k1s = [k1] + [0] * (w-1)
    k2s = [k2] + [0] * (w-1)
    row_index = len(df) - 1

    # Go through all statistics, compute them, and set them in the df columns of the last row
    df.loc[row_index, 'indep.theta_pi_1']   = independent_check.theta_pi(n1, c1s, k1s, 2)
    df.loc[row_index, 'indep.theta_pi_2']   = independent_check.theta_pi(n2, c2s, k2s, 2)
    df.loc[row_index, 'indep.theta_w_1']    = independent_check.theta_w(n1, c1s, k1s, 2)
    df.loc[row_index, 'indep.theta_w_2']    = independent_check.theta_w(n2, c2s, k2s, 2)
    df.loc[row_index, 'indep.achaz_var_1']  = independent_check.achaz_var_d(n1, c1s, k1s, 2)
    df.loc[row_index, 'indep.achaz_var_2']  = independent_check.achaz_var_d(n2, c2s, k2s, 2)
    with warnings.catch_warnings():
        # Some of the computations divide by zero, which spams our output...
        warnings.simplefilter("ignore")
        df.loc[row_index, 'indep.tajimas_d_1']  = independent_check.tajimas_d(n1, c1s, k1s, 2)
        df.loc[row_index, 'indep.tajimas_d_2']  = independent_check.tajimas_d(n2, c2s, k2s, 2)
    df.loc[row_index, 'indep.fst_nei']      = independent_check.fst_nei(n1, c1s, k1s, n2, c2s, k2s)
    df.loc[row_index, 'indep.fst_hudson']   = independent_check.fst_hudson(n1, c1s, k1s, n2, c2s, k2s)
    df.loc[row_index, 'indep.fst_karlsson'] = independent_check.fst_karlsson(n1, c1s, k1s, n2, c2s, k2s)
    df.loc[row_index, 'indep.fst_kofler']   = independent_check.fst_kofler(n1, c1s, k1s, n2, c2s, k2s)

# ------------------------------------------------------------
#     grenedalf
# ------------------------------------------------------------

def get_grenedalf_result(file, column):
    res = pd.read_csv(file, delimiter=',')
    return float(res[column].iloc[0])

def run_grenedalf_diversity_sample( n, c, k, w ):
    file_id = f"mpileup-cov_{c}-der_{k}-win_{w}"
    out_id = file_id + f"-pool_{n}"

    args = {}
    args["fileid"] = f"../mpileup/{file_id}.pileup"
    args["poolsize"] = n
    args["windowsize"] = w
    args["outid"] = out_id

    run_script("grenedalf/diversity.sh", args)
    return out_id

def run_grenedalf_diversity( n, c, k, w, df, sample ):
    # n: pool size
    # c: coverage
    # k: derived allele count
    # w: window size
    row_index = len(df) - 1

    # Run for all measures once. no need to repeat
    out_id = run_grenedalf_diversity_sample( n, c, k, w )
    file_id = f"mpileup-cov_{c}-der_{k}-win_{w}"

    # Names in the table and in grenedalf are different
    measure_name_map = {
        "theta_pi":  file_id + ".1.theta_pi",
        "theta_w":   file_id + ".1.theta_watterson",
        "tajimas_d": file_id + ".1.tajimas_d"
    }

    for measure in [ "theta_pi", "theta_w", "tajimas_d" ]:
        df_col = 'grenedalf.' + measure + '_' + str(sample)
        tab_col = measure_name_map[measure]
        val = get_grenedalf_result(f"grenedalf/diversity/diversity-{out_id}.csv", tab_col)
        df.loc[row_index, df_col] = val

def run_grenedalf_fst_sample( n1, n2, c1, c2, k1, k2, w, method ):
    file_id = f"sync-cov1_{c1}-cov2_{c2}-der1_{k1}-der2_{k2}-win_{w}"
    out_id = file_id + f"-pool1_{n1}-pool2_{n2}-method_{method}"

    # We need to make a pool size file, and specify the two sample names
    os.makedirs("poolsizes", exist_ok=True)
    with open(f"poolsizes/{out_id}.pool", "w") as f:
        f.write( f"{file_id}.1\t{n1}\n")
        f.write( f"{file_id}.2\t{n2}\n")

    args = {}
    args["fileid"] = f"../sync/{file_id}.sync"
    # args["poolsizes"] = '"' + str(n1) + ',' + str(n2) + '"'
    args["poolsizes"] = f"../poolsizes/{out_id}.pool"
    args["windowsize"] = w
    args["method"] = method
    args["outid"] = out_id

    run_script("grenedalf/fst.sh", args)
    return out_id

def run_grenedalf_fst( n1, n2, c1, c2, k1, k2, w, df ):
    # n: pool size
    # c: coverage
    # k: derived allele count
    # w: window size
    row_index = len(df) - 1
    file_id = f"sync-cov1_{c1}-cov2_{c2}-der1_{k1}-der2_{k2}-win_{w}"

    # Names in the table and in grenedalf are different
    method_name_map = {
        "unbiased-nei": "fst_nei",
        "unbiased-hudson": "fst_hudson",
        "karlsson": "fst_karlsson",
        "kofler": "fst_kofler"
    }

    # Need to run for each method independently
    for method in ["unbiased-nei","unbiased-hudson","kofler","karlsson"]:
        df_col = 'grenedalf.' + method_name_map[method]
        out_id = run_grenedalf_fst_sample( n1, n2, c1, c2, k1, k2, w, method )
        tab_col = f"{file_id}.1:{file_id}.2.fst"
        df.loc[row_index, df_col] = get_grenedalf_result(f"grenedalf/fst/fst-{out_id}.csv", tab_col)

def run_grenedalf( params, df ):
    # n: pool size
    # c: coverage
    # k: derived allele count
    # w: window size
    n1, n2, c1, c2, k1, k2, w = params
    row_index = len(df) - 1

    # Go through all statistics, compute them, and set them in the df columns of the last row
    run_grenedalf_diversity( n1, c1, k1, w, df, 1 )
    run_grenedalf_diversity( n2, c2, k2, w, df, 2 )
    run_grenedalf_fst( n1, n2, c1, c2, k1, k2, w, df )

# ------------------------------------------------------------
#     PoPoolation
# ------------------------------------------------------------

def get_popoolation_result(file, trim_fst_pair=False):
    with open(file) as f:
        line = f.readline()
    try:
        # PoPoolation stores the resulting value in the last column.
        # Go get it. If it's FST, we also need to remove the sample names.
        last = line.split('\t')[-1]
        if trim_fst_pair:
            last = last.split('=')[-1]
        val = float(last)
    except Exception as e:
        val = float("nan")
    return val



def run_popoolation_diversity_sample( n, c, k, w, measure ):
    file_id = f"mpileup-cov_{c}-der_{k}-win_{w}"
    out_id = file_id + f"-pool_{n}-measure_{measure}"

    args = {}
    args["fileid"] = f"../mpileup/{file_id}.pileup"
    args["poolsize"] = n
    args["windowsize"] = w
    args["measure"] = measure
    args["outid"] = out_id

    run_script("popoolation/diversity.sh", args)
    return out_id

def run_popoolation_diversity( n, c, k, w, df, sample ):
    # n: pool size
    # c: coverage
    # k: derived allele count
    # w: window size
    row_index = len(df) - 1

    # Names in the table and in popoolation are different
    measure_name_map = {
        "theta_pi":  "pi",
        "theta_w":   "theta",
        "tajimas_d": "d"
    }

    # Need to run for each measure independently
    for measure in [ "theta_pi", "theta_w", "tajimas_d" ]:
        df_col = 'popoolation.' + measure + '_' + str(sample)
        out_id = run_popoolation_diversity_sample( n, c, k, w, measure_name_map[measure] )
        val = get_popoolation_result(f"popoolation/diversity/{out_id}")
        df.loc[row_index, df_col] = val

def run_popoolation_fst_sample( n1, n2, c1, c2, k1, k2, w, method ):
    file_id = f"sync-cov1_{c1}-cov2_{c2}-der1_{k1}-der2_{k2}-win_{w}"
    out_id = file_id + f"-pool1_{n1}-pool2_{n2}-method_{method}"

    args = {}
    args["fileid"] = f"../sync/{file_id}.sync"
    args["poolsizes"] = '"' + str(n1) + ':' + str(n2) + '"'
    args["windowsize"] = w
    args["method"] = method
    args["outid"] = out_id

    run_script("popoolation/fst.sh", args)
    return out_id

def run_popoolation_fst( n1, n2, c1, c2, k1, k2, w, df ):
    # n: pool size
    # c: coverage
    # k: derived allele count
    # w: window size
    row_index = len(df) - 1

    # Names in the table and in popoolation are different
    method_name_map = {
        "karlsson": "fst_karlsson",
        "kofler": "fst_kofler"
    }

    # Need to run for each method independently
    for method in ["kofler","karlsson"]:
        df_col = 'popoolation.' + method_name_map[method]
        out_id = run_popoolation_fst_sample( n1, n2, c1, c2, k1, k2, w, method )
        df.loc[row_index, df_col] = get_popoolation_result(f"popoolation/fst/{out_id}.fst", trim_fst_pair=True)

def run_popoolation( params, df ):
    # n: pool size
    # c: coverage
    # k: derived allele count
    # w: window size
    n1, n2, c1, c2, k1, k2, w = params
    row_index = len(df) - 1

    # Go through all statistics, compute them, and set them in the df columns of the last row
    run_popoolation_diversity( n1, c1, k1, w, df, 1 )
    run_popoolation_diversity( n2, c2, k2, w, df, 2 )
    run_popoolation_fst( n1, n2, c1, c2, k1, k2, w, df )

# ------------------------------------------------------------
#     Main
# ------------------------------------------------------------

# Run everything but the bagel
def main():
    # Create the files that we need to run grenedalf and popoolation
    create_files.create_from_param_space( param_space )

    # Columns of the resulting dataframe that we are filling
    fixed_cols = [
        'pool_size_1', 'pool_size_2', 'coverage_1', 'coverage_2',
        'derived_count_1', 'derived_count_2', 'window_size'
    ]
    indep_cols = [
        'indep.theta_pi_1', 'indep.theta_pi_2', 'indep.theta_w_1', 'indep.theta_w_2',
        'indep.achaz_var_1', 'indep.achaz_var_2', 'indep.tajimas_d_1', 'indep.tajimas_d_2',
        'indep.fst_nei', 'indep.fst_hudson', 'indep.fst_karlsson', 'indep.fst_kofler'
    ]
    grenedalf_cols = [
        'grenedalf.theta_pi_1', 'grenedalf.theta_pi_2',
        'grenedalf.theta_w_1', 'grenedalf.theta_w_2',
        'grenedalf.tajimas_d_1', 'grenedalf.tajimas_d_2',
        'grenedalf.fst_nei', 'grenedalf.fst_hudson',
        'grenedalf.fst_karlsson', 'grenedalf.fst_kofler'
    ]
    popoolation_cols = [
        'popoolation.theta_pi_1', 'popoolation.theta_pi_2',
        'popoolation.theta_w_1', 'popoolation.theta_w_2',
        'popoolation.tajimas_d_1', 'popoolation.tajimas_d_2',
        'popoolation.fst_karlsson', 'popoolation.fst_kofler'
    ]

    # Create the dataframe for all results
    df = pd.DataFrame(columns=[fixed_cols + indep_cols + grenedalf_cols + popoolation_cols])

    # Run all tests of the whole param space
    print("Analyzing files")
    for params in tqdm.tqdm(param_space):
        # Set up the parameters as needed
        n1, n2, c1, c2, k1, k2, w = params
        # print("At", params)

        # Add an empty row to the df, and get its index
        df.loc[len(df)] = [np.nan] * len(df.columns)
        row_index = len(df) - 1

        # Add the param values to the cells of that row
        df.loc[row_index, 'pool_size_1']     = n1
        df.loc[row_index, 'pool_size_2']     = n2
        df.loc[row_index, 'coverage_1']      = c1
        df.loc[row_index, 'coverage_2']      = c2
        df.loc[row_index, 'derived_count_1'] = k1
        df.loc[row_index, 'derived_count_2'] = k2
        df.loc[row_index, 'window_size']     = w

        # Run all implementations
        run_independent_implementation( params, df )
        run_grenedalf( params, df )
        # run_popoolation( params, df )

    # Store the whole dataframe
    df.to_csv(outtable, index=False, sep='\t')

if __name__ == '__main__':
    main()
