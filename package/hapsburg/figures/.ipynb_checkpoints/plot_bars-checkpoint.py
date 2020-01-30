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

from hapsburg.PackagesSupport.roh_expectations import Expected_Roh

def create_cousins_roh(degrees=[1,2,3], 
                       bins=[[0.04,0.08],[0.08,0.12],[0.12,0.2],[0.2,3.0]], bin_n=10000):
    """Create ROH sharing in list of bins (list of [begin,end]) 
    for Cousins of degree degrees [list]
    return sharing [len(degrees), len(bins)]"""
    e_roh = Expected_Roh()
    c_roh = np.zeros((len(degrees),len(bins))) # Container for results Cousins
    for i,c in enumerate(degrees):
        for j,b in enumerate(bins):
            m = c*2 + 4
            c_roh[i,j] = e_roh.exp_roh_len_in_bin_rel(l=b, m=m, comm_anc=4, bins=10000)
    return c_roh

#bins = [[0.04,0.08],[0.08,0.12],[0.12,3.00]]  # The bins I want to plot (eventually maybe do 12,16 as welll)
#cousins = [1, 2, 3]  # Which Cousins to Plot

def create_Ne_roh(Ns=[400, 800, 1600, 3200, 6400], 
                  bins=[[0.04,0.08],[0.08,0.12],[0.12,0.2],[0.2,3.0]], bin_n=10000):
    """Create ROH sharing in list of bins (list of [begin,end]) 
    for panmictic population sizes
    Ns: List of population sizes
    bins: Length Bins (in Morgan) to calculate expectations from
    return sharing [len(degrees), len(bins)]"""
    e_roh = Expected_Roh()
    n_roh = np.zeros((len(Ns),len(bins))) # Container for results Cousins
    for i,N in enumerate(Ns):
        for j,b in enumerate(bins):
            n_roh[i,j] = e_roh.exp_roh_len_in_bin_N(b, N=N, bins=bin_n)
    return n_roh

def std_Ne_roh(Ns=[400, 800, 1600, 3200, 6400], 
                  bins=[[0.04,0.08],[0.08,0.12],[0.12,0.2],[0.2,3.0]], bin_n=10000):
    """Create ROH sharing in list of bins (list of [begin,end]) 
    for panmictic population sizes
    Ns: List of population sizes
    bins: Length Bins (in Morgan) to calculate expectations from
    return sharing [len(degrees), len(bins)]"""
    e_roh = Expected_Roh()
    var_roh = np.zeros((len(Ns),len(bins))) # Container for results Cousins
    for i,N in enumerate(Ns):
        for j,b in enumerate(bins):
            var_roh[i,j] = e_roh.var_roh_len_in_bin_N(b, N=N, bins=bin_n)
    return np.sqrt(var_roh)  # Return Standard Deviation

######################################################################################

def plot_bar_ax(ax, fracs, bins=[], c=["#313695", "#abd9e9", "#fee090", "#d7191c"], x_ticks = [], 
                ec = "silver", fs_l=10, fs = 10, fs_x=8, barWidth=0.95, ylim = [0,220], stds = [], title="", 
                yticks=False, legend=False, r_title=0, hlines=[], ylabel="Sum Inferred ROH>4cM [cM]"):
    """Plot bars of ROH on Axis.
    ax: Where to Plot on
    fracs: Fractions of Expected ROH: [n Inds, n Bins]
    c: Which colors to plot
    bins: List of Bins (needed for legend - plotted if len()>0)
    yticks: Whether to plot Y tick Labels
    legend: Whether to plot Legend
    fs_x: Fontsize on the x Labels
    r_title: Rotation of the title
    hlines: List where to plot hlines"""
    y = fracs ### Fractions (along 1st axis length classes)
    x = np.arange(len(y))

    for i in range(len(y[0,:])): # From last to first (For Legend)
        b = np.sum(y[:,:i], axis=1)
        ax.bar(x, y[:,i], bottom=b, color=c[i], edgecolor=ec, width=barWidth, label=f"{bins[i,0]}-{bins[i,1]} cM", alpha=0.9)
        if len(stds)>0 and i>0: # Plot some standard deviations.
            ax1.errorbar(r, b, yerr=stds[:,i], fmt='none', linewidth=2, color="k")
    
    if len(hlines)>0:
        for y in hlines:
            ax.axhline(y=y, zorder=0, linestyle="--", color="gray", lw=0.5)     
    
    if legend:
        ax.legend(fontsize=fs_l, loc="upper right", title="Sum ROH in")
    ax.set_ylabel(ylabel, fontsize=fs)
    ax.set_ylim(ylim)
    ax.set_xlim(x[0] - 0.7*barWidth, x[-1] + 0.7*barWidth)
    if len(x_ticks)>0:
        ax.set_xticks(x)
        ax.set_xticklabels(x_ticks, fontsize=fs_x, rotation=270)
    else:
        ax.set_xticks([])
    if not yticks:
        ax.set_yticklabels([])
        ax.set_ylabel("")
    if len(title)>0:
        ax.set_title(title, fontsize=fs, rotation=r_title)
        
        
def plot_panel_row(plot_dfs, wspace=0.05, hspace=0.01, figsize=(24,3.5), savepath="", x_ticks=True,
                   c=["#313695", "#abd9e9", "#fee090", "#d7191c"], ylim = [0,250], r_title = 90, 
                   leg_pos = -2, show=True, title_col="clst", legends=True, hlines=[],
                   cols = ['sum_roh>4', 'sum_roh>8', 'sum_roh>12', 'sum_roh>20'],
                   bins = [[0.04, 0.08], [0.08, 0.12], [0.12, 0.2], [0.2, 3.0]],
                   degrees=[1, 2, 3], Ns=[400, 800, 1600, 3200, 6400]):
    """Plot row of ROH bin plots from plot_dfs (each df one panel)
    leg_pos = Where to plot legend (if outside range no legend plot)
    r_title: How much to rotate the title
    gs: Gridspec: If given plot on there.
    legends: Whether to plot the two legends
    hlines: Where to plot horizontal lines
    cols: List of Column Names for plot_dfs (assumes > in increasing order)
    bins: list of length bins to plot [[a1,a2],...[z1,z2]]
    Ns: What population sizes to plot in barplot [list]
    degrees: What degrees of Cousins to plot. [list]"""
    bins_cM=(np.array(bins)*100).astype("int")
    n_plots = len(plot_dfs)
    width_ratios = [len(df) for df in plot_dfs]
    if legends:     
        n_plots+= 2
        width_ratios+=[len(degrees), len(Ns)]
        
    fig = plt.figure(figsize=figsize)
    gs = gridspec.GridSpec(1, n_plots, width_ratios=width_ratios, figure=fig)
    
    gs.update(wspace=wspace, hspace=hspace) # set the spacing between axes

    for i,df in enumerate(plot_dfs):   
        if i==0:
            ytick=True
        else:
            ytick=False

        if i == (len(plot_dfs) + leg_pos):
            legend=True
        else:
            legend=False
    
        ax = plt.subplot(gs[i])    # Extract the Sub Plot to Plot onto
        obs_roh = df[cols].values
        
        ### Calculate the value in the Bins
        for i in range(len(cols)-1):
            obs_roh[:,i] = obs_roh[:,i] - obs_roh[:,i+1]
        
        if x_ticks:
            x_ticks0 = df["iid"].values
        else: 
            x_ticks0 = []
        plot_bar_ax(ax, obs_roh, bins_cM, yticks=ytick, legend=legend, r_title=r_title, c=c,
                    x_ticks = x_ticks0, title=df[title_col].values[0], ylim=ylim, hlines=hlines)

    #### Plot the legend bars
    if legends:
        ### 0) Calculate Expectations
        c_roh = create_cousins_roh(degrees = degrees, bins = bins)
        sum_roh = create_Ne_roh(Ns=Ns, bins = bins)
        
        ### 1) Small Pops
        ax = plt.subplot(gs[-2])    # The left subplot (Timeline)
        plot_bar_ax(ax, c_roh*100, bins_cM, yticks=False, legend=False, ylim=ylim, c=c, hlines=hlines,
                   x_ticks = ["1st C.", "2nd C.", "3rd C."], title="Recent Loops", r_title=r_title)
        ### 2) Cousins
        ax = plt.subplot(gs[-1])
        ticks_N = [f"2N={i}" for i in Ns]
        plot_bar_ax(ax, sum_roh*100, bins_cM, yticks=False, legend=False, ylim=ylim, c=c,
                    hlines=hlines, x_ticks = ticks_N, title="Small Pop. Size", r_title=r_title)

    if len(savepath)>0:
        plt.savefig(savepath, bbox_inches = 'tight', pad_inches = 0, dpi=300)
        print(f"Saved figure to {savepath}")
        
    if show:
        plt.show()
    return
    
    
def plot_legend_only(figsize=(7,6), wspace=0.05, hspace=0.01, savepath="",
                     fs=12, fs_l=10, fs_x=10, c=["#313695", "#abd9e9", "#fee090", "#d7191c"],
                     bins = [[0.04, 0.08], [0.08, 0.12], [0.12, 0.2], [0.2,3.0]],
                     degrees=[1, 2, 3], Ns=[400, 800, 1600, 3200, 6400]
                     ):
    """Plot Inbreeding from recent Cousins as well as small pop size.
    bins: list of length bins to plot [[a1,a2],...[z1,z2]]
    Ns: What population sizes to plot in barplot [list]
    degrees: What degrees of Cousins to plot. [list]"""
    width_ratios = [len(degrees), len(Ns)]
    bins_cM=(np.array(bins)*100).astype("int")
    fig = plt.figure(figsize=figsize)
    gs = gridspec.GridSpec(1, 2, width_ratios=width_ratios, figure=fig)
    ax_cousin = plt.subplot(gs[0])    # The left subplot (Timeline)
    ax_Ne = plt.subplot(gs[1])
    gs.update(wspace=wspace, hspace=hspace) # set the spacing between axes

    ### Calcualte Expectations Cousins:
    c_roh = create_cousins_roh(degrees = degrees, bins = bins)
    
    ### Calculate Expectations Ne:
    sum_roh = create_Ne_roh(Ns=Ns, bins = bins)
    
    plot_bar_ax(ax_cousin, c_roh*100, bins_cM, yticks=True, legend=False, 
                fs=fs, fs_l=fs_l, fs_x=fs_x, c=c, ylabel="Expected Sum ROH>4cM [cM]",
                x_ticks = ["1st C.", "2nd C.", "3rd C."], title="Recent Loops")

    ticks_N = [f"2N={i}" for i in Ns]
    plot_bar_ax(ax_Ne, sum_roh*100, bins_cM, yticks=False, legend=True, c=c, 
                fs=fs, fs_l=fs_l, fs_x=fs_x, x_ticks = ticks_N, 
                title="Small Pop. Size")
            
    if len(savepath)>0:
        plt.savefig(savepath, bbox_inches = 'tight', pad_inches = 0, dpi=300)
        print(f"Saved figure to {savepath}") 
    plt.show()
    
    
def prepare_dfs_plot(df, cms=[4,8,12], col_group="clst", split_modern=True,
                     mod_group="pop", sortcol=0):
    """Prepare the Dataframe which to plot
    Return split up (and sorted) list of df, return list of column names
    df: Master Dataframe
    cms: Minimum Length of Chromosomes in bins
    sortcol: By which > column to sort within group [int]
    if split_modern, split age==0 samples by mod_group
    """
    
    if split_modern:
        mod_idx = df["age"]==0  # Pull out modern idx
        df.loc[mod_idx, col_group]=df.loc[mod_idx, mod_group]
        
    plot_dfs = [dft for _, dft in df.groupby(col_group)]
    ### Sort by age
    idx = np.argsort([-df["age"].values[0] for df in plot_dfs])
    plot_dfs = [plot_dfs[i] for i in idx] ## Sort
    
    ### Split up blocks and sort by lowest
    cols = [f"sum_roh>{cm}" for cm in cms]
    for df in plot_dfs:
        df.sort_values(by=cols[sortcol], inplace=True, ascending=False)
    return plot_dfs, cols