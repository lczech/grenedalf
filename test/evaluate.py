#!/usr/bin/env python3

# Matplotlib shenennigans...
# https://stackoverflow.com/a/71511579/4184258
import matplotlib
matplotlib.use('Agg')
# matplotlib.use('TkAgg')

# libraries
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import numpy as np
import pandas as pd
import seaborn as sns
import os, sys
import scipy as sp
from matplotlib.ticker import FuncFormatter

# ------------------------------------------------------------
#     Settings
# ------------------------------------------------------------

# global names for our ins and outs
out_dir_png = "figures_png"
out_dir_pdf = "figures_pdf"
out_dir_svg = "figures_svg"

infile = "test_results.csv"
df = pd.read_csv(infile, sep='\t')

# ------------------------------------------------------------
#     Plotting
# ------------------------------------------------------------

plotnum = 1
def plot_corr(
    title, xlabel, ylabel, x, y, with_kernel = True
    # markers='o', dashes=(0, ""), palette='k', **kwargs
):
    global plotnum
    global out_dir_png
    global out_dir_pdf
    global out_dir_svg

    # Get proper name and print it
    testcase = title + '-' + xlabel + '-' + ylabel
    testcase = testcase.lower().replace( " ", "-" )
    testcase = ''.join(c for c in testcase if c.isalnum() or c == '-')
    testcase = testcase.replace("--", "-")
    # print("Plot: " + testcase)

    # print(x)
    # print(y)
    if len(x) != len(y):
        print("Warning: x has " + str(len(x)) + " rows, y has " + str(len(y)) + " rows")

    # Remake data as one frame, because of reasons
    # (being: we need to drop na, but in both at the same time)
    df = pd.DataFrame(data={'x': x, 'y': y})
    df = df[df['x'].notna()]
    df = df[df['y'].notna()]
    # print (df)

    # Create dirs for results
    if not os.path.exists(out_dir_png):
        os.makedirs(out_dir_png)
    if not os.path.exists(out_dir_pdf):
        os.makedirs(out_dir_pdf)
    if not os.path.exists(out_dir_svg):
        os.makedirs(out_dir_svg)

    # Compute the pearson correlation between both, and the mean squared error
    pearson_r, pearson_p = sp.stats.pearsonr(x=df['x'], y=df['y'])
    # mse = ((df['x'] - df['y'])**2).mean()
    mae = (np.absolute(df['x'] - df['y'])).mean()

    # Let's make a fancy density plot to be able to better see the dot density
    # Also, we subset for the kernel computation, as it takes waaaay to long otherwise...
    if with_kernel:
        max_df_size = 25000
        if len(df.index) > max_df_size:
            print("  subsetting")
            df = df.sample(n=max_df_size)
        values = np.vstack([df['x'], df['y']])
        colormap = plt.cm.viridis_r
        print("  kernel")
        kernel = sp.stats.gaussian_kde(values)(values)

        # Make the plot
        print("  plot", len(values[0]))
        plt.figure( plotnum, figsize=(8,8) )
        ax = sns.scatterplot(
            x=df["x"], y=df["y"],
            c=kernel,
            cmap=colormap,
            linewidth=0
        )
        # sns.kdeplot(
        #      x=df["x"], y=df["y"], fill=True, cmap=colormap
        # )
    else:
        # Make the plot
        plt.figure( plotnum, figsize=(8,8) )
        ax = sns.scatterplot(
            x=df["x"], y=df["y"],
            linewidth=0
        )

    # Make proper
    ax.set_box_aspect(1)
    ax.set_aspect('equal')

    # Draw a line of x=y
    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()
    lims = [max(x0, y0), min(x1, y1)]
    ax.plot( lims, lims, 'lightgray', zorder=-1)
    # ax.plot( lims, lims, 'k', alpha=0.2 )
    # plt.plot([0, 1], [0, 1])

    # annotate the pearson correlation coefficient text to 2 decimal places
    plt.text(.05, .95, 'r={:.2f}'.format(pearson_r), transform=ax.transAxes)
    # plt.text(.05, .92, 's={:.6f}'.format(mse), transform=ax.transAxes)
    plt.text(.05, .92, 'MAE={:.4f}'.format(mae), transform=ax.transAxes)

    # Nameing the plot and the axis.
    ax.set_title( title )
    ax.set_xlabel( xlabel )
    ax.set_ylabel( ylabel )

    # Save to files
    plt.savefig( out_dir_png + "/" + testcase + ".png", format='png', bbox_inches="tight" )
    # plt.savefig( out_dir_pdf + "/" + testcase + ".pdf", format='pdf', bbox_inches="tight" )
    # plt.savefig( out_dir_svg + "/" + testcase + ".svg", format='svg', bbox_inches="tight" )
    plt.close( plotnum )
    plotnum += 1

# ------------------------------------------------------------
#     Evaluation
# ------------------------------------------------------------

# We run an evaluation that tests correlation and mean absolute error.
# If either of them is out of the range of what we consider good,
# we set the success flag to 1 (false), indicating failure in the return value of the script.
success = 0
def eval_corr( title, xlabel, ylabel, x, y ):
    # print(title, xlabel, ylabel)

    # Prepare the data frame. Similar to the plotting above.
    df = pd.DataFrame(data={'x': x, 'y': y})
    df = df[df['x'].notna()]
    df = df[df['y'].notna()]

    # Compute the pearson correlation between both, and the mean squared error.
    # Again, same as above in the plotting. But well, we can live with that code duplication.
    pearson_r, pearson_p = sp.stats.pearsonr(x=df['x'], y=df['y'])
    mae = (np.absolute(df['x'] - df['y'])).mean()

    # If we fail, we just set the flag. We do not abort yet,
    # as we want all tests to run, so that we get all results.
    # Makes it easier for debugging.
    if pearson_r < 0.99 or mae > 0.001:
        print("FAIL", title, xlabel, ylabel)
        print("pearson_r == " + str(pearson_r))
        print("mae == " + str(mae))
        success = 1

    # Lastly also create a plot for evaluation purposes/
    plot_corr( title, xlabel, ylabel, x, y, with_kernel = False )

# ------------------------------------------------------------
#     grenedalf vs independent
# ------------------------------------------------------------

# Theta Pi
eval_corr( "theta_pi_1", "grenedalf", "independent", df["grenedalf.theta_pi_1"], df["indep.theta_pi_1"] )
eval_corr( "theta_pi_2", "grenedalf", "independent", df["grenedalf.theta_pi_2"], df["indep.theta_pi_2"] )

# Theta W
eval_corr( "theta_w_1", "grenedalf", "independent", df["grenedalf.theta_w_1"], df["indep.theta_w_1"] )
eval_corr( "theta_w_2", "grenedalf", "independent", df["grenedalf.theta_w_2"], df["indep.theta_w_2"] )

# Tajima's D
eval_corr( "tajimas_d_1", "grenedalf", "independent", df["grenedalf.tajimas_d_1"], df["indep.tajimas_d_1"] )
eval_corr( "tajimas_d_2", "grenedalf", "independent", df["grenedalf.tajimas_d_2"], df["indep.tajimas_d_2"] )

# Fst
eval_corr( "fst_nei", "grenedalf", "independent", df["grenedalf.fst_nei"], df["indep.fst_nei"] )
eval_corr( "fst_hudson", "grenedalf", "independent", df["grenedalf.fst_hudson"], df["indep.fst_hudson"] )
eval_corr( "fst_karlsson", "grenedalf", "independent", df["grenedalf.fst_karlsson"], df["indep.fst_karlsson"] )
eval_corr( "fst_kofler", "grenedalf", "independent", df["grenedalf.fst_kofler"], df["indep.fst_kofler"] )

# ------------------------------------------------------------
#     grenedalf vs PoPoolation
# ------------------------------------------------------------

# By default, we do not run popoolation, to save time, and because that's not our job here.
# Hence, we check here if the columns are filled. If not, we are done here.
if len(df["popoolation.theta_pi_1"].value_counts()) == 0:
    sys.exit( success )

# Furthermore, in the tests below we do not run the evaluation function,
# as we do not want to fail for PoPoolation reasons, and expect slight differences
# in the numerical values. In particular Tajima's D is known to be bugged in PoPoolation.
# So instead, we just plot the results here.

# Theta Pi
plot_corr( "theta_pi_1", "grenedalf", "popoolation", df["grenedalf.theta_pi_1"], df["popoolation.theta_pi_1"], with_kernel = False )
plot_corr( "theta_pi_2", "grenedalf", "popoolation", df["grenedalf.theta_pi_2"], df["popoolation.theta_pi_2"], with_kernel = False )

# Theta W
plot_corr( "theta_w_1", "grenedalf", "popoolation", df["grenedalf.theta_w_1"], df["popoolation.theta_w_1"], with_kernel = False )
plot_corr( "theta_w_2", "grenedalf", "popoolation", df["grenedalf.theta_w_2"], df["popoolation.theta_w_2"], with_kernel = False )

# Tajima's D
plot_corr( "tajimas_d_1", "grenedalf", "popoolation", df["grenedalf.tajimas_d_1"], df["popoolation.tajimas_d_1"], with_kernel = False )
plot_corr( "tajimas_d_2", "grenedalf", "popoolation", df["grenedalf.tajimas_d_2"], df["popoolation.tajimas_d_2"], with_kernel = False )

# Fst
plot_corr( "fst_karlsson", "grenedalf", "popoolation", df["grenedalf.fst_karlsson"], df["popoolation.fst_karlsson"], with_kernel = False )
plot_corr( "fst_kofler", "grenedalf", "popoolation", df["grenedalf.fst_kofler"], df["popoolation.fst_kofler"], with_kernel = False )

# ------------------------------------------------------------
#     independent vs PoPoolation
# ------------------------------------------------------------

# Theta Pi
plot_corr( "theta_pi_1", "independent", "popoolation", df["indep.theta_pi_1"], df["popoolation.theta_pi_1"], with_kernel = False )
plot_corr( "theta_pi_2", "independent", "popoolation", df["indep.theta_pi_2"], df["popoolation.theta_pi_2"], with_kernel = False )

# Theta W
plot_corr( "theta_w_1", "independent", "popoolation", df["indep.theta_w_1"], df["popoolation.theta_w_1"], with_kernel = False )
plot_corr( "theta_w_2", "independent", "popoolation", df["indep.theta_w_2"], df["popoolation.theta_w_2"], with_kernel = False )

# Tajima's D
plot_corr( "tajimas_d_1", "independent", "popoolation", df["indep.tajimas_d_1"], df["popoolation.tajimas_d_1"], with_kernel = False )
plot_corr( "tajimas_d_2", "independent", "popoolation", df["indep.tajimas_d_2"], df["popoolation.tajimas_d_2"], with_kernel = False )

# Fst
plot_corr( "fst_karlsson", "independent", "popoolation", df["indep.fst_karlsson"], df["popoolation.fst_karlsson"], with_kernel = False )
plot_corr( "fst_kofler", "independent", "popoolation", df["indep.fst_kofler"], df["popoolation.fst_kofler"], with_kernel = False )

# Final report of the success of this evaluation
sys.exit( success )
