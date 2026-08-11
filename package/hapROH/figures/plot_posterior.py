"""
Main Inference Class for HMM. Wrapper for Inerence of Posterior.
@ Author: Harald Ringbauer, 2019, All rights reserved
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colorbar as clb
from matplotlib import gridspec
import os as os
import sys as sys
import pandas as pd

import logging
import warnings

logger = logging.getLogger(__name__)

NUMPY_FILES = ["hap",           # shape (2, nb_snp), int, for compatibility with hapsburg only
               "readcounts",    # shape (nb_snp, 2), int
               "gt_count",      # shape (nb_snp), int
               "pos",           # shape (nb_snp), int
               "map",           # shape (nb_snp), float
               "posterior0"     # shape (nb_snp), int
               ]

ROH_FILES = ["roh", "roh_gt"]

def load_data(folder:str):
    """Load and return the Data from one Data Folder"""

    data = dict()
    # load numpy arrays
    for file in NUMPY_FILES:
        file_path = os.path.join(folder, file+".csv")
        if os.path.isfile(file_path):
            file_data = np.genfromtxt(file_path, dtype="float", delimiter=",")
            data[file] = file_data
    # load ROH dataframes
    for file in ROH_FILES:
        file_path = os.path.join(folder, file+".csv")
        if os.path.isfile(file_path):
            df_roh = pd.read_csv(file_path, delimiter=",")
            data[file] = df_roh

    logger.debug(f"Loaded following data: {data.keys()}")
    return data

def plot_posterior_cm(folder:str, savepath:None|str=None,
                      plot_calls=True, min_cm:float=1,
                      plot_hets:bool=True, plot_post:bool=True,
                      x_lim:None|tuple=None, min_reads=1,
                      figsize=(14,4), title:None|str=None, post_c="maroon", het_c="blue") -> plt.Figure:
    """
    Args:
    folder: Path to input data
    savepath: Path were to save resulting figure
    plot_calls: Whether to plot ROH Calls
    min_cm: Minimum length [centimorgans] for called ROH
    plot_hets: Whether to plot Heterozygote Markers
    plot_post: Whether to plot posterior probability
    x_lim: What area to zoom in (CentiMorgan)

    min_reads: How many reads of each REF and ALT to be considered heterozygous (only used if plot_hets and data contains readcounts but not gt_count)

    empirical: If true, do not load and plot latent states
    cm_lim: What Area to Zoom In (CentiMorgan)
    m: How many reads for ref and alt
    yticks: Where to place the Y ticks
    groundtruth: Whether to plot Ground Truth (saved as csv). 
    Only used in simulated data with known ROH
    plot: Whether to show the plot in python"""

    fs = 14  
    roh_lw = 6   # Linewidth for ROH

    data = load_data(folder)

    pos_x = 100*data["map"]

    fig = plt.figure(figsize=figsize)
    ax = plt.subplot()

    ### Plot the roh
    if plot_calls:
        if not "roh" in data:
            warnings.warn("No roh found in data. Argument plot_calls is ignored")
        else:
            df_roh = data["roh"]
            df_roh = df_roh[100*df_roh["lengthM"] >= min_cm]
            ax.hlines(y=np.full(len(df_roh), 1.2), xmin=100*df_roh["StartM"], xmax=100*df_roh["EndM"],
                    colors="blue", linewidth=roh_lw)

    # TODO: plot groundtruth ?

    ### Plot the posterior probability
    if plot_post:
        if not "posterior0" in data:
            warnings.warn("No posterior found in data. Argument plot_post is ignored")
        else:
            posterior0 = data["posterior0"]
            ax.plot(pos_x, posterior0, linewidth=2, color=post_c, zorder=1)
            ax.set_ylabel("Post. probability", fontsize=fs, color=post_c)

    ### Plot the heterozygotes
    if plot_hets:
        ylabel = f"Heterozygote (no/yes)"
        if "gt_count" in data:
            het = data["gt_count"] == 1
        elif "hap" in data:
            het = data["hap"][0] != data["hap"][1]
        elif "readcounts" in data:
            ylabel = r"Both REF and ALT reads $\geq$ {min_reads}"
            het = (data["readcounts"][:, 0] >= min_reads) & (data["readcounts"][:, 1] >= min_reads)
        else:
            warnings.warn("No genotype found in data. Argument plot_hets is ignored")
            het = None
        if het is not None:
            ax.plot(pos_x, (het * 1.1 - 0.05), "o", ms=1, alpha=0.3, zorder=0, color=het_c)
            ax2 = ax.twinx()
            ax2.set_ylim(ax.get_ylim())
            ax2.set_yticks(np.array([1,0]) * 1.1 - 0.05)
            ax2.set_yticklabels([])
            ax2.set_ylabel(ylabel)

    ### Customise the plot
    ax.set_xlabel("Genetic position (cM)", fontsize=fs)

    if x_lim is not None:
        ax.set_xlim(x_lim)

    if title is not None:
        plt.title(title, fontsize=fs)
        
    if savepath is not None:
        plt.savefig(savepath, bbox_inches='tight', pad_inches=0, dpi=300)
        print(f"Saved figure to: {savepath}")

    return fig