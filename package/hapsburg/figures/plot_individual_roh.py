"""
Plotting function to plot karyotype plot of ROH from hapsburg output
@ Author: Harald Ringbauer, 2019, All rights reserved
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import gridspec
import matplotlib.patheffects as pe
import os as os
import sys as sys
from scipy.stats import gaussian_kde

### For Arial Font
from matplotlib import rcParams
rcParams['font.family'] = 'sans-serif'   # Set the defaul
rcParams['font.sans-serif'] = ['Arial']  # Make sure to have the font installed (it is on cluster for Harald)

from hapsburg.PackagesSupport.pp_individual_roh_csvs import merge_called_blocks, post_process_roh_df

def load_roh(iid, ch, path_folder = "./Empirical/1240k/", prefix_out = "e01/"):
    """Load the ROH Dataframe for Individual iid and 
    Chromosome ch"""
    path = path_folder +  iid + "/chr" + str(ch) + "/" + prefix_out+ "/roh.csv"
    roh_df = pd.read_csv(path)    
    return roh_df

def post_process_roh_df1(df, min_cm=4, snp_cm=100, output=False):
    """Post Process ROH Dataframe.
    min_cm: Minimum Length in CentiMorgan
    snp_cm: How many SNPs per CentiMorgan"""
    densities = df["length"] / (df["lengthM"] * 100)
    densities_ok = (densities > snp_cm)
    
    df["SNP_Dens"] = densities
    
    # Filter for SNP Density:
    df = df[densities_ok]
    
    # Filter for Length:
    length_okay = (df["lengthM"] * 100) > min_cm
    
    if output==True:
        print(f"Min SNPs per cM> {snp_cm}: {np.sum(densities_ok)}/{len(densities_ok)}")
        print(f"> {min_cm} cM: {np.sum(length_okay)}/{len(length_okay)}")
    
    df = df[length_okay]
    return df


def chrom_length(ch, ch_lengths=[], output=False):
    """Get and return length of Chromosome ch
    If ch_lengths not given, use default (from Eigenstrat map).
    Atm only do autosomes!"""
    if len(ch_lengths)==0:
        ch_lengths = [2.86273, 2.688325, 2.232573, 2.145423,
                      2.040858, 1.920325, 1.871526, 1.680022,
                      1.66301, 1.809153, 1.582171, 1.746799,
                      1.257046, 1.202023, 1.41346, 1.340263,
                      1.284738, 1.177099, 1.077316, 1.082134,
                      0.627865, 0.740762]
        return ch_lengths[ch-1]

def load_individual_roh(iid, min_cm=4, gap=0.0, snp_cm=50, path_folder = "./Empirical/1240k/", 
                        prefix_out = "", folder="./Empirical/1240k/", output=False):
    """Load ROH of one Individual"""
    df_rohs = []
    for i in range(1,23):
        df_roh = load_roh(iid=iid, ch=i, path_folder = folder, prefix_out = prefix_out)
        df_roh = merge_called_blocks(df_roh, max_gap=gap/100, output=output) # Merge Blocks
        df_roh = post_process_roh_df(df_roh, min_cm=min_cm, snp_cm=snp_cm, output=output) # Only use Blocks with high enough SNP density
        df_rohs.append(df_roh)
    return df_rohs

def load_bad_areas(path="./Data/1000Genomes/Markers/1240k/snp_density.csv", min_snps=50):
    """Load areas of low SNP density, and return list of Dataframes 
    (one for each chromosome)"""
    df_lows = []
    
    for i in range(1,23):
        df_t = pd.read_csv(path, sep="\t")
        df_t = df_t[df_t["chr"]==i]
        df_t = df_t[df_t["counts"]<min_snps]
        df_lows.append(df_t)
    return df_lows

def plot_chromosome(ax, l, x_pos, lw=24, df_roh = [], df_low = []):
    """Plot a Chromosome of length l with centromer ctr on ax 
    at x_pos"""
    ln, = ax.plot([x_pos, x_pos], [-0.05,l+0.05], lw = lw, color="lightgray",
                      solid_capstyle = 'round', zorder=0,
                  path_effects=[pe.Stroke(linewidth=lw+3, foreground='k'), pe.Normal()])
    
    ### Plot the ROH List if given
    if len(df_roh) > 0:
        starts, ends = df_roh["StartM"].values, df_roh["EndM"].values
        
        for i in range(len(df_roh)):
            ax.plot([x_pos, x_pos], [starts[i], ends[i]], lw=lw, color="maroon", 
                    zorder=1, alpha=1.0, solid_capstyle="butt")
    
    ### Plot shadows of bad areas
    if len(df_low)>0:
        starts, ends = df_low["StartM"].values, df_low["EndM"].values
        
        for i in range(len(df_low)):
            ax.plot([x_pos, x_pos], [starts[i], ends[i]], lw=lw, color="k", 
                    zorder=2, alpha=0.8, solid_capstyle="butt")
        
            
def plot_roh_individual(iid="MA89", fs=12, figsize=(8,8), savepath="", min_cm=4, snp_cm=50, gap=0.0,
                        folder="./Empirical/1240k/MarcusAncs/", prefix_out="", plot_bad=True,
                        title="True", output=False):
    """Plot ROH in one ancient Individual.
    gap: What Gap to Merge [in cM!]
    prefix_out: If there is a folder before indivual data, e.g. e01/"""
    
    ### Load the Data (could be done seperately)
    df_rohs = load_individual_roh(iid, min_cm=min_cm, snp_cm=snp_cm, gap=gap,
                                  folder=folder, prefix_out=prefix_out, output=output)
    if plot_bad:
        df_lows = load_bad_areas(min_snps=snp_cm)  # Load low density areas of the genome
    else:
        df_lows=[[] for _ in range(23)]
    
    plt.figure(figsize=figsize)
    gs = gridspec.GridSpec(2, 1, height_ratios=[3, 2])
    gs.update(hspace=0.1) # set the spacing between axes. 

    ax = plt.subplot(gs[0]) # The left subplot
    ax1 = plt.subplot(gs[1]) # The middle one
        
    ### Plot the First 11 Chromosomes
    for i in range(1,12):
        l = chrom_length(i)
        plot_chromosome(ax, l, x_pos=(i), df_roh=df_rohs[i-1], df_low=df_lows[i-1])  
    
    ### Plot the Second 11 Chromosomes
    for i in range(12,23):
        l = chrom_length(i)
        plot_chromosome(ax1, l, x_pos=(i - 11), df_roh=df_rohs[i-1], df_low=df_lows[i-1])

    ### Set the Plot Limits
    ax.set_xlim([0.3, 11.5])
    ax1.set_xlim([0.3, 11.5])

    ax.set_ylim([-0.3,3.3])
    ax1.set_ylim([-0.3, 2.05])

    ### Turn off the Y Axis
    for ax0 in [ax, ax1]:
        ax0.spines['right'].set_visible(False)
        ax0.spines['top'].set_visible(False)
        ax0.spines['bottom'].set_visible(False)
        ax0.yaxis.set_tick_params(labelsize=fs)

    rang = np.arange(1,12)
    ax.set_xticks(rang)
    ax.set_xticklabels(rang, fontsize=fs)

    ax1.set_xticks(rang)
    ax1.set_xticklabels(rang+11, fontsize=fs)
    ax1.set_xticklabels(rang+11, fontsize=fs)

    if len(savepath)>0:
        plt.savefig(savepath, bbox_inches = 'tight', pad_inches = 0, dpi=300)
        print(f"Saved figure to {savepath}")
    
    if title:
        ax.set_title(iid, fontsize=fs*2)
    plt.show()
    
    
##########################################################
### Plot Length Distributions

def expected_block_pdf(x, chr_l, m):
    """Gives back the pdfs for Blocks of Length l [Morgan]
    on a Chromosome of Length [Morgan].
    m: Nr of Meiosis.
    Return PDF (per Morgan)"""
    pdf0 = (chr_l - x) * m**2 * np.exp(-x * m)
    pdf1 = 2 * np.exp(- x * m) * m 
    return (pdf0 + pdf1) * (x < chr_l)  # If x < chr_l return 0

def expected_block_pdf_chromosomes(x, chr_lgts, m):
    """Calculate the PDF of ROH blocks of length x [Morgan]
    for m Recombination events
    x: Can be Array
    chr_lgts: Array of all chromosome lengths [in Morgan]
    Return PDF (per Morgan)"""
    pdfs = [expected_block_pdf(x, chr_l, m) for chr_l in chr_lgts]
    pdf_full = np.sum(pdfs, axis=0)
    return pdf_full

def coal_prob(m, comm_anc=1):
    """Calculate Coalescence Probability.
    m: Nr of Meiosis
    comm_anc: How many common ancestors"""
    c_prob = comm_anc * (1 / 2) ** m
    return c_prob

def exp_blocks_full_individual(x, m, comm_anc=1):
    """Calculates the Expected Nr of ROH Blocks per Morgan for full Individual
    x: Array of Block Lenths
    m: Nr of Meisois
    comm_anc: Nr of Ancestry Loops"""
    
    chr_lgts = [chrom_length(ch) for ch in range(1,23)]
    
    pdf_full = expected_block_pdf_chromosomes(x, chr_lgts, m)
    c_prob = coal_prob(m, comm_anc)
    exp_blocks = pdf_full * c_prob
    return exp_blocks


def plot_pde_individual(iid="MA89", figsize=(8,6),
                        min_cm=4, snp_cm=50, bw_cm=4, kde_plot=False, 
                        plotlim=[4,100], savepath="", 
                        folder="./Empirical/1240k/", prefix_out="e01/", output=False, gap=0.0, lw_curve=3,
                        comm_ancs=[4,4,4,2], ms=[6,5,4,3], labels=["First Cousins", "Aunt/Nephew", "Full Siblings", "Parent/Offpsring"],
                        cs=["red", "green", "orange", "gray"], title="", leg_loc="upper right"):
    """Plot Histograms/PDEs of ROH Distribution for one Individual (iid)
    If list of k iids given, load all ROH and plot all of them, and the expected value for k individuals.
    bw_cm: Length of one Bin (in cM)
    comm_ancs: How many common ancestors to plot [list]
    ms: How many meiosis to plot [list]
    labels: what labels do they have [list]
    cs: what colors to plot [list]"""

    ### Load And Prepare Data
    if isinstance(iid, str):
        df_rohs = load_individual_roh(iid=iid, min_cm=min_cm, snp_cm=snp_cm, 
                                      folder=folder, prefix_out=prefix_out, output=output, gap=gap)
        df_roh = pd.concat(df_rohs)
        k=1
    else: 
        df_iids = [pd.concat(load_individual_roh(iid=i, min_cm=min_cm, snp_cm=snp_cm, 
                                      folder=folder, prefix_out=prefix_out, output=output, gap=gap)) for i in iid]
        df_roh = pd.concat(df_iids)
        k = len(iid)
        print(f"Loaded ROH of {k} individuals. Plotting their sum...")

    bins = np.arange(plotlim[0], plotlim[1], bw_cm)
    bin_mean = (bins[1:] + bins[:-1]) / 2.0  # Mean of each bin
    
    ### Do the Plot
    fs = 16

    plt.figure(figsize=figsize)
    plt.hist(df_roh["lengthM"]*100, bins=bins, ec="k", fc="dodgerblue", label="Observed ROH")
    
    # Plot the Empirical Averages
    for i in range(len(labels)):
        block_pdf = exp_blocks_full_individual(bin_mean/100, m=ms[i], comm_anc=comm_ancs[i]) * k
        plt.plot(bin_mean, bw_cm * block_pdf/100, c=cs[i], label=labels[i], lw=lw_curve) # Plot Density Per cM (adjusted for bin width)
    
    ### DO KDE Plot
    if kde_plot==True:
        kde = gaussian_kde(df_roh["lengthM"]*100) # KDE per cM
        plt.plot(bin_mean, bw_cm * kde(bin_mean) * len(df_roh), "k--", label="KDE of observed ROH", lw=lw_curve)

    plt.xlabel("ROH Length [cm]", fontsize=fs)
    plt.ylabel(f"Number per {bw_cm} cM Bin", fontsize=fs)
    plt.title(title, fontsize=fs)
    leg = plt.legend(loc = leg_loc, fontsize=fs)
    leg.set_title("Parents being...", prop = {'size':fs})
    
    plt.tick_params(axis='both', which='major', labelsize=fs)
    
    if len(savepath)>0:
        plt.savefig(savepath, bbox_inches = 'tight', pad_inches = 0, dpi=300)
        print(f"Saved to {savepath}.")
    plt.show()