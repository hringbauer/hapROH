"""
Classes and functions to plot the simulated and processed Mosaic ROH Data
Python Functions to keep functions and analysis seperate
@ Author: Harald Ringbauer, 2019, All rights reserved
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import gridspec

def plot_power(bl_lens, df_call_vec1, powers, df_fp=[], fs = 12, fs_l=12, 
               xlim=(0,12.5), figsize=(10,6), n=100, ylim_pow=[0.5, 1.05],
               lw_power=1.2, color_fp="red", ec="silver", pw_yticks=[0.5,0.75,1.0], s=100,
               alpha=0.8, savepath="", title=""):
    """ bl_lens: Array of Block Lengths
        df_call_vec1: Array of Called Blocks
        powers: Array of Power to call Blocks
    Load, and plot power and power curves
    s: Size of ticks on y axis
    pw_ticks: Where to set the ticks for the power
    lw_power: Line width of power plot""" 
    assert(len(bl_lens)==len(df_call_vec1)==len(powers)) # Sanity Check
    
    bins = np.linspace(0, 15, 75)  # Bins of 0.1 cM
    
    ### Set Colors
    cmap = cm.get_cmap("viridis_r")
    colors = [cmap(x) for x in np.linspace(0,1, len(bl_lens))]
    
    ####### Do the actual Plot
    plt.figure(figsize=figsize)

    gs = gridspec.GridSpec(2, 1, height_ratios=[1, 4])
    gs.update(hspace=0.04) # set the spacing between axes. 

    ax = plt.subplot(gs[1]) # The lower subplot
    ax1 = plt.subplot(gs[0]) # The upper subplot
    
    ### Set Font Sizes
    ax.tick_params(axis='both', labelsize=fs)
    ax1.tick_params(axis='both', labelsize=fs)

    # Plot All Histograms        
    for i in range(len(bl_lens)):
        l = bl_lens[i]
        ax.hist(df_call_vec1[i]["CalledLength"], bins = bins, color = colors[i], alpha=alpha, 
                label= str(l) + " cM", ec=ec)
        ax.axvline(l, color = "gray", linewidth=2)
    
    if len(df_fp)>0:
        ax.hist(df_fp["lengthM"]*100, bins = bins, color = color_fp, 
                alpha=1.0, label= "FP", ec=ec)
        
    ax.set_xlim(xlim)
    ax.set_xlabel("Inferred ROH Length [cM]", fontsize = fs)
    ax.set_ylabel("# Inferred \nROH blocks", fontsize = fs)
    legend = ax.legend(loc = "upper right", fontsize = fs_l, title="Simulated ROH")
    legend.get_title().set_fontsize(fs_l)
    
    ### Plot the upper panel
    ax1.set_ylabel("Called at \n80% overlap", fontsize = fs)
    ax1.set_xticks([])
    ax1.yaxis.set_ticks(pw_yticks)
    for y in pw_yticks:
        ax1.axhline(y, linestyle='--', color='gray', lw=0.3, zorder=0)
    ax1.scatter(bl_lens, powers, c=colors, s=s, zorder=2, ec=ec)
    ax1.plot(bl_lens, powers, "gray", zorder=1, lw=lw_power)
    ax1.set_xlim(xlim)
    ax1.set_ylim(ylim_pow)

    plt.title(title, fontsize=fs)
    
    if len(savepath) >0:
        plt.savefig(savepath)
    plt.show()
    
    
def plot_fp_distribution(df_call_fp, fs = 14, fs_l = 12, bins = np.linspace(1, 3, 26), xlim = (0, 13), figsize = (10, 4),
                        title = "100 Mosaic Tuscany Samples - Reference: Rest EUR 1000G"):
    """ Plot the Distribution of false positive ROH calls"""
    ### The Actual Plot
    plt.figure(figsize=figsize)

    ax = plt.gca()
    ax.hist(df_call_fp["lengthM"] * 100, bins = bins, color = "green", alpha=0.9, ec="silver")
    ax.set_xlabel("Inferred block Length [cm]", fontsize=fs)
    ax.set_ylabel("Count", fontsize=fs)
    roh1cm = np.sum(df_call_fp["lengthM"]>0.01)
    ax.text(x=0.6, y=0.8, s=f"Tot. #ROH > 1cM: {roh1cm}", transform=ax.transAxes, fontsize=fs)  
    
    plt.title(title, fontsize=fs)
    plt.show()
    
    
    

