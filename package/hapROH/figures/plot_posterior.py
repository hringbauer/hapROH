import logging
import os
import warnings
from typing import Literal

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.figure import Figure

logger = logging.getLogger(__name__)


# Potential files created by hapROH
NUMPY_FILES = [
    "hap",  # shape (2, nb_snp), int
    "readcounts",  # shape (nb_snp, 2), int
    "pos",  # shape (nb_snp), int
    "map",  # shape (nb_snp), float
    "posterior0",  # shape (nb_snp), int
]

ROH_FILES = ["roh", "roh_gt"]


def load_data(folder: str):
    """Load and return the Data from one Data Folder"""
    data = {}
    # load numpy arrays
    for file in NUMPY_FILES:
        file_path = os.path.join(folder, file + ".csv")
        if os.path.isfile(file_path):
            file_data = np.genfromtxt(file_path, dtype="float", delimiter=",")
            data[file] = file_data
    # load ROH dataframes
    for file in ROH_FILES:
        file_path = os.path.join(folder, file + ".csv")
        if os.path.isfile(file_path):
            df_roh = pd.read_csv(file_path, delimiter=",")
            data[file] = df_roh

            required = {"ch", "StartM", "EndM", "lengthM"}
            # + "StartBP", "EndBP" when using unit = "BP" BUT columns not present in older version (StartPosGRCh37 instead)
            assert required.issubset(df_roh.columns), (
                f"Missing columns: {required - set(df_roh.columns)} in file {file_path}"
            )

    logger.debug(f"Loaded following data: {data.keys()}")
    return data


def plot_posterior(
    folder: str,
    iid: str,
    chrom: str,
    savepath: None | str = None,
    x_lim: None | tuple = None,
    unit: Literal["BP", "M", "cM"] = "cM",
    plot_calls=True,
    min_cm: float = 1,
    plot_hets: bool = True,
    min_reads=1,
    plot_post: bool = True,
    figsize=(14, 4),
    title: None | str = None,
    post_c="maroon",
    het_c="blue",
) -> Figure:
    """
    Args:
    folder: Path to hapROH output folder
    savepath: Path where to save resulting figure
    x_lim: What area to zoom in
    unit: One of 'BP' [base pairs], 'M' [Morgans] or 'cM' [centimorgans]
    plot_calls: Whether to plot ROH calls
    min_cm: Minimum length [centimorgans] for called ROH (only used if plot_calls)
    plot_hets: Whether to plot Heterozygote Markers
    min_reads: How many reads of each REF and ALT to be considered heterozygous
            (only used if plot_hets and data contains AD but not GT)
    plot_post: Whether to plot posterior probability
    """

    ### Plot settings
    fs = 14  # fontsize
    roh_lw = 6  # ROH line width

    folder_data = os.path.join(folder, iid, "chr" + str(chrom), "")
    data = load_data(folder_data)

    match unit:
        case "BP":
            pos_x = data["pos"]
            xlabel = "Physical position (BP)"
            roh_start, roh_end, roh_coeff = "StartBP", "EndBP", 1
        case "M":
            pos_x = data["map"]
            xlabel = "Genetic position (M)"
            roh_start, roh_end, roh_coeff = "StartM", "EndM", 1
        case "cM":
            pos_x = 100 * data["map"]
            xlabel = "Genetic position (cM)"
            roh_start, roh_end, roh_coeff = "StartM", "EndM", 100
        case _:
            raise ValueError(
                f"Unknown unit {unit}. Should be one of 'BP', 'M' or 'cM'."
            )

    fig = plt.figure(figsize=figsize)
    ax = plt.subplot()

    ### Plot the roh
    if plot_calls:
        if not "roh" in data:
            warnings.warn("No roh found in data. Argument plot_calls is ignored")
        else:
            df_roh = data["roh"]
            df_roh = df_roh[df_roh["lengthM"] >= 0.01 * min_cm]
            ax.hlines(
                y=np.full(len(df_roh), 1.2),
                xmin=roh_coeff * df_roh[roh_start],
                xmax=roh_coeff * df_roh[roh_end],
                colors=het_c,
                linewidth=roh_lw,
            )

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
        ylabel = "Heterozygote (no/yes)"
        if "hap" in data:
            het = data["hap"][0] != data["hap"][1]
        elif "readcounts" in data:
            ylabel = rf"Both REF and ALT reads $\geq$ {min_reads}"
            het = (data["readcounts"][:, 0] >= min_reads) & (
                data["readcounts"][:, 1] >= min_reads
            )
        else:
            warnings.warn("No genotype found in data. Argument plot_hets is ignored")
            het = None
        if het is not None:
            ax.plot(
                pos_x, (het * 1.1 - 0.05), "o", ms=1, alpha=0.3, zorder=0, color=het_c
            )
            ax2 = ax.twinx()
            ax2.set_ylim(ax.get_ylim())
            ax2.set_yticks(np.array([1, 0]) * 1.1 - 0.05)
            ax2.set_yticklabels([])
            ax2.set_ylabel(ylabel, color=het_c)

    ### Customise the plot
    ax.set_xlabel(xlabel, fontsize=fs)

    if x_lim is not None:
        ax.set_xlim(x_lim)

    if title is None:
        title = f"Individual {iid} — Chromosome {chrom}"
    plt.title(title, fontsize=fs)

    if savepath is not None:
        plt.savefig(savepath, bbox_inches="tight", pad_inches=0, dpi=600)
        print(f"Saved figure to: {savepath}")

    return fig
