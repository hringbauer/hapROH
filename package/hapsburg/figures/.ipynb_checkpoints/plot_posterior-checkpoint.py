"""
Main Inference Class for HMM. Wrapper for Inerence of Posterior.
@ Author: Harald Ringbauer, 2019, All rights reserved
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colorbar as clb
from matplotlib import gridspec
import os as os
import sys as sys
import pandas as pd

# For nice Arial Fonts
from matplotlib import rcParams
rcParams['font.family'] = 'sans-serif'   # Set the default
rcParams['font.sans-serif'] = ['Arial']  # Make sure to have the font installed (it is on cluster for Harald)


def load_data(folder="../Simulated/Example0/", empirical=False, fullpost=False, viterbi=False, readcounts=True):
    """Load and return the Data from one Data Folder"""
    
    ob_stat = np.loadtxt(folder + "hap.csv", dtype="int", delimiter=",")
    
    if viterbi:
        viterbi_path = np.loadtxt(folder + "viterbi_path.csv", dtype="int", delimiter=",")
    else:
        viterbi_path = []
    
    if empirical==False:
        lats = np.loadtxt(folder + "lat.csv", dtype="int", delimiter=",")
        ref_states = np.loadtxt(folder + "refs.csv", dtype="int", delimiter=",")
        roh_df = 0
        read_counts = 0
        gmap = 0
        posterior = np.loadtxt(folder + "posterior.csv", dtype="float", delimiter=",")
        
    elif empirical:
        lats = 0
        roh_df  = pd.read_csv(folder + "roh.csv", delimiter=",")
        if readcounts==True:
            read_counts = np.loadtxt(folder + "readcounts.csv", delimiter=",")
        else:
            read_counts= 0
            
        gmap = np.loadtxt(folder + "map.csv", dtype="float", delimiter=",")
        
        if fullpost==True:
            posterior = np.loadtxt(folder + "posterior.csv", dtype="float", delimiter=",")  # For Plot of all Posteriors   
        else:
            posterior = np.loadtxt(folder + "posterior0.csv", dtype="float", delimiter=",")
            
        ref_states = np.zeros(np.shape(posterior)) # Just some Filler for the Moment
        
    else:
        raise RuntimeError("WTF DUDE")

    print(f"Successfully loaded Data from: {folder}")
    return (ref_states, ob_stat, lats, viterbi_path, 
            posterior, roh_df, read_counts, gmap)


def process_read_counts(read_counts, m=1):
    """Return Readcount that have at least m reads for both Ref and Alt"""
    hets = (read_counts[0,:]>=m) & (read_counts[1,:]>=m)
    return hets


def plot_posterior_cm(folder = "../Simulated/Test20r/", savepath="", empirical=True, 
                      plot=True, cm_lim=[], m=1, groundtruth=False, readcount=False, 
                     plot_hets=True, plot_calls=True, plot_post=True,
                      figsize=(14,4), title="", post_c="maroon", het_c="blue"):
    """Plot Viterbi Path of Haplotype copying.
    save: Whether to save the results.
    empirical: If true, do not load and plot latent states
    cm_lim: What Area to Zoom In (CentiMorgan)
    m: How many reads for ref and alt
    groundtruth: Whether to plot Ground Truth (saved as csv)
    plot_hets: Whether to plot Heterozygote Markers
    plot_calls: Whether to plot Calls
    plot_post: Whether to plot posterior"""
    
    cmap = cm.get_cmap("viridis")
    norm = plt.Normalize(-8, 0)

    fs = 14  
    lw = 6   # Linewidth for ROH
    
    _, ob_stat,_,_,posterior, roh_df, read_counts, gmap = load_data(folder, empirical, readcounts=readcount)
    assert(len(gmap)==np.shape(ob_stat)[1])
        
    ###########################
    ###########################
    # Do the Plot
    
    plt.figure(figsize=figsize)
    gs = gridspec.GridSpec(2, 1, height_ratios=[4, 8])
    gs.update(hspace=0.3) # set the spacing between axes. 

    #ax = plt.subplot(gs[0]) # The upper subplot
    ax1 = plt.subplot(gs[1]) # The lower subplot
    
    ### Depict Readcounts or GT Heteroyzgotes
    if readcount==False:
        het = (ob_stat[0,:] != ob_stat[1,:])
    elif readcount==True:
        het = process_read_counts(read_counts, m=m)
    else:
        raise RuntimeError("Invalid Mode!")
        
    #ax.plot(gmap*100, het, "bo", ms=2, alpha=0.3)
    #ax.set_ylabel("Het", fontsize=fs)
    #ax.set_title("Hetorzygosity", fontsize=fs)
    if plot_post:
        ax1.plot(gmap*100, np.exp(posterior), linewidth=2, color=post_c, label="State 0", zorder=1)
        ax1.set_ylabel("Post. Probability", fontsize=fs, color=post_c)
    
    if plot_calls:
        plt.hlines(y=[1.2]*len(roh_df), xmin=100 * roh_df["StartM"], xmax= 100 * roh_df["EndM"], 
                        colors="blue", linewidth=lw)
    
    ax1.set_xlabel("SNP", fontsize=fs)
    ax1.set_xlabel("CentiMorgan")
    
    if len(cm_lim)==2:
        ax1.set_xlim(cm_lim)
        
    ### Load & Plot Groundtruth (if given):
    if groundtruth:
        path = folder + "roh_gt.csv"
        dft = pd.read_csv(path, sep="\t")
        ### Plot them
        plt.hlines(y=[1.12]*len(dft), xmin=100 * dft["ROH_Begin"], xmax=100 * dft["ROH_End"], 
                    colors="green", linewidth=lw)
    
    ### Plot Heterozygotes, if plot_hets
    if plot_hets:
        ax1.plot(gmap*100, (het * 1.1 - 0.05), "o", ms=1, alpha=0.3, zorder=0, color=het_c)
        ax2 = ax1.twinx()
        ax2.set_ylim(ax1.get_ylim())
        ax2.set_yticks(np.array([1,0]) * 1.1 - 0.05)
        ax2.set_yticklabels([])
        ax1.set_yticklabels([])
        ax2.set_ylabel(f"$\geq/\geq$ {m} Ref/Alt Reads", fontsize=fs*0.7, color=het_c)
        
    if len(title)>0:
        plt.title(title, fontsize=fs)
        
    if len(savepath)>0:
        plt.savefig(savepath, bbox_inches = 'tight', pad_inches = 0, dpi=300)
        print(f"Saved figure to: {savepath}")
        #plt.savefig(folder + "posterior_cm.png", bbox_inches = 'tight', pad_inches = 0, dpi=300)
    
    if plot==True:
        plt.show()



