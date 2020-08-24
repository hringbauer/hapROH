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
                       bins=[[0.04,0.08],
                             [0.08,0.12],
                             [0.12,0.2],
                             [0.2,3.0]], bin_n=10000):
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

def plot_bar_ax(ax, y, bins=[], c=["#313695", "#abd9e9", "#fee090", "#d7191c"], x_ticks = [], 
                ec = "silver", fs_l=10, fs_y = 10, fs_x=8, fs_t=10, alpha=1.0,
                barWidth=0.95, ylim = [0,220], stds = [], 
                title="", ha_title="left", bold_title=False,
                yticks=False, ylabel="Sum Inferred ROH>4cM [cM]",
                legend=False, r_title=0, hlines=[]):
    """Plot bars of ROH on Axis.
    ax: Where to Plot on
    y: Array of ROH to plot: [n Inds, k Legnth Bins]
    c: Which colors to plot
    bins: List of Bins (needed for legend - plotted if len()>0)
    yticks: Whether to plot Y tick Labels
    legend: Whether to plot Legend
    fs_x, fs_y: Fontsize on the x and yLabels
    r_title: Rotation of the title
    hlines: List where to plot hlines"""
    x = np.arange(len(y))

    for i in range(len(y[0,:])): # From last to first (For Legend)
        b = np.sum(y[:,:i], axis=1)
        ax.bar(x, y[:,i], bottom=b, color=c[i], edgecolor=ec, width=barWidth, 
                   label=f"{bins[i,0]}-{bins[i,1]} cM", alpha=alpha, zorder=4)
        if len(stds)>0 and i>0: # Plot some standard deviations.
            ax1.errorbar(r, b, yerr=stds[:,i], fmt='none', linewidth=2, color="k")
    
    if len(hlines)>0:
        for y in hlines:
            ax.axhline(y=y, zorder=0, linestyle="--", color="gray", lw=0.5)     
    
    if legend:
        ax.legend(fontsize=fs_l, loc="upper right", title="Sum ROH in")
    ax.set_ylabel(ylabel, fontsize=fs_y)
    ax.set_ylim(ylim)
    ax.set_xlim(x[0] - 0.7*barWidth, x[-1] + 0.7*barWidth)
    
    ### Do the xtricks
    ax.set_xticks(x)
    if len(x_ticks)>0:
        ax.set_xticklabels(x_ticks, fontsize=fs_x, rotation=270)
    else:
        ax.set_xticklabels([])
    if not yticks:
        ax.axes.yaxis.set_ticks([])
        #ax.set_yticklabels([])
    fw = 'normal'
    if bold_title:
        fw = "bold"
        
    if len(title)>0:
        ax.set_title(title, fontsize=fs_t, rotation=r_title, 
                     horizontalalignment=ha_title, fontweight=fw)
        
def plot_close_kin(ax, y, y_plot=100,cutoffs=[50,100], c="r", ec="k", lw=1,
                   ss=[10,14], m_cs=["v", "s"], bin_idx=-1):
    """Plot Symbols for close kin at height y_symbol
    ss: Sizes of markers [List]
    m_cs: marker_colors [List]
    cutoffs: cutoffs [List]
    """
    if len(y)==0:  # Move out if nothing to plot
        return
    x = np.arange(len(y))
    y0 = y[:,bin_idx]
    
    cutoffs_t = cutoffs + [4000]
    
    for i in range(len(cutoffs)):
        idx = (y0>cutoffs_t[i]) & (y0<cutoffs_t[i+1])
        if np.sum(idx)==0:
            continue
        ax.scatter(x[idx], [y_plot]*np.sum(idx), c=c, ec=ec, lw=lw, 
                   marker=m_cs[i], s=ss[i], zorder=5)
        
def plot_panel_row(plot_dfs, wspace=0.05, hspace=0.01, figsize=(24,3.5), savepath="", 
                   x_labels=[], ylabel="Sum Inferred ROH>4cM [cM]", c=["#313695", "#abd9e9", "#fee090", "#d7191c"], 
                   ylim = [0,250], r_title = 90, bolds=[],
                   fs_l=10, fs_y = 10, fs_x=8, fs_t=10, ha_title="left", hspace_leg=1,
                   leg_pos = -2, show=True, title_col="clst", titles=[], hlines=[],
                   cols = ['sum_roh>4', 'sum_roh>8', 'sum_roh>12', 'sum_roh>20'],
                   bins = [[0.04, 0.08], [0.08, 0.12], [0.12, 0.2], [0.2, 3.0]],
                   degrees=[1, 2, 3], Ns=[400, 800, 1600, 3200, 6400],
                   ticks_c=["1st C.", "2nd C.", "3rd C."],
                   cutoffs=[], ec="k", lw=1, ss=[40, 50], 
                   m_cs=["v", "s"], sym_ofst=-10, bin_idx=-1):
    """Plot row of ROH bin plots from plot_dfs (each df one panel)
    leg_pos = Where to plot legend (if outside range no legend plot)
    r_title: How much to rotate the title
    bolds: which titles to bold [list of bools]
    ha_title: Horizontal alignment of the titles
    fs_l, fs_y, fs_x, fs_t: Fontsize of legend, y and x labels and titles
    hspace_leg: Horizontol space between data plots and legend
    gs: Gridspec: If given plot on there
    legends: Whether to plot the two legends
    titles: If given, list of titles
    title_col: If not title use this column of the dataframe
    hlines: Where to plot horizontal lines
    cols: List of Column Names for plot_dfs (assumes > in increasing order)
    bins: list of length bins to plot [[a1,a2],...[z1,z2]]
    Ns: What population sizes to plot in barplot [list]
    degrees: What degrees of Cousins to plot. [list]
    ticks_c: Tick Labels for Cousin Legend [list]
    x_labels: If FALSE, don't plot any xlabels, if empty list default [list]"""
    bins_cM=(np.array(bins)*100).astype("int")
    n_plots0 = len(plot_dfs) # The original plots
    n_plots = len(plot_dfs) + 1 # Add 1 for the empty plot
    width_ratios = [len(df) for df in plot_dfs]
    width_ratios += [hspace_leg]

    ### Whether to make space for legends:
    if len(degrees)>0: 
        n_plots += 1
        width_ratios+=[len(degrees)]
    if len(Ns)>0: 	    
        n_plots += 1
        width_ratios+=[len(Ns)]
        
    fig = plt.figure(figsize=figsize)
    gs = gridspec.GridSpec(1, n_plots, width_ratios=width_ratios, figure=fig)
    gs.update(wspace=wspace, hspace=hspace) # set the spacing between axes

    for i,df in enumerate(plot_dfs):   
        if i == (len(plot_dfs) + leg_pos):
            legend=True
        else:
            legend=False
    
        ax = plt.subplot(gs[i])    # Extract the Sub Plot to Plot onto
        obs_roh = df[cols].values
        
        ### Calculate the value in the Bins
        for j in range(len(cols)-1):
            obs_roh[:,j] = obs_roh[:,j] - obs_roh[:,j+1]
        
        if x_labels == False:
            x_ticks0 = []
        else:
            if len(x_labels)>0:
                x_ticks0 = x_labels[i]
            else:
                x_ticks0 = df["iid"].values
        
        if len(titles)>0:
            title=titles[i]
        else:
            title=df[title_col].values[0]
        if len(bolds)>0:
            bold=bolds[i]
        else:
            bold=False
            
        if i==0:
            ylabel1 = ylabel
            yticks=True
        else:
            ylabel1 = ""
            yticks = False
                  
        plot_bar_ax(ax, obs_roh, bins_cM, legend=legend, r_title=r_title, c=c,
                    x_ticks = x_ticks0, title=title, ha_title=ha_title, bold_title=bold,
                    ylim=ylim, hlines=hlines, ylabel=ylabel1, yticks=yticks,
                    fs_l=fs_l, fs_y = fs_y, fs_x=fs_x, fs_t=fs_t)
        
        if len(cutoffs)>0:
            c_marker = np.array(c)[bin_idx]
            plot_close_kin(ax, obs_roh, y_plot=ylim[1]+sym_ofst, cutoffs=cutoffs, 
                           c=c_marker, ec=ec, lw=lw, 
                           ss=ss, m_cs=m_cs, bin_idx=bin_idx)
    
    ### Make the placeholder axis invisible
    ax_none = plt.subplot(gs[n_plots0])
    ax_none.set_visible(False)

    ###### Plot the legend bars
    ### 1) Cousins
    if len(degrees)>0:
        c_roh = create_cousins_roh(degrees = degrees, bins = bins)
        ax_c = plt.subplot(gs[n_plots0+1])    # Get the right axis
        plot_bar_ax(ax_c, c_roh*100, bins_cM, legend=False, ylim=ylim, c=c, 
                    hlines=hlines,x_ticks = ticks_c, 
                    yticks=False, ylabel="",
                    fs_l=fs_l, fs_y = fs_y, fs_x=fs_x, fs_t=fs_t,
                    title="Recent Loops", r_title=r_title)
    
    ### 2) Small Pops    
    if len(Ns)>0:
        ns_roh = create_Ne_roh(Ns=Ns, bins = bins)
        pos_leg_N = n_plots0 + 1 + (len(degrees)>0)  # Don't forget Python Indexing
        ax_N = plt.subplot(gs[pos_leg_N])
        ticks_N = [f"2N={i}" for i in Ns]
        plot_bar_ax(ax_N, ns_roh*100, bins_cM, legend=False, ylim=ylim, c=c,
                    yticks=False, ylabel="",
                    fs_l=fs_l, fs_y = fs_y, fs_x=fs_x, fs_t=fs_t, 
                    hlines=hlines, x_ticks = ticks_N, title="Small Pop. Size", r_title=r_title)

    if len(savepath)>0:
        plt.savefig(savepath, bbox_inches = 'tight', pad_inches = 0, dpi=300)
        print(f"Saved figure to {savepath}")
        
    if show:
        plt.show()
    return
    
    
def plot_legend_only(figsize=(7,6), wspace=0.05, hspace=0.01, savepath="", hlines=[],
                     fs_l=10, fs_y = 10, fs_x=8, fs_t=10, c=["#313695", "#abd9e9", "#fee090", "#d7191c"],
                     bins = [[0.04, 0.08], [0.08, 0.12], [0.12, 0.2], [0.2,3.0]], ha_title="center",
                     degrees=[1, 2, 3], Ns=[400, 800, 1600, 3200, 6400], ylim = [0,250],
                     x_ticks_c = ["1st C.", "2nd C.", "3rd C."], title_c = "Recent Loops",
                     title_N="Small Pop. Size", y_label="Expected Sum ROH>4cM [cM]"
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
    c_roh = create_cousins_roh(degrees = degrees, bins = bins) # in Morgan
    
    ### Calculate Expectations Ne:
    sum_roh = create_Ne_roh(Ns=Ns, bins = bins) # in Morgan
    
    plot_bar_ax(ax_cousin, c_roh*100, bins_cM, legend=False, 
                fs_l = fs_l, fs_y = fs_y, fs_x = fs_x, fs_t = fs_t, c = c,
                ylim = ylim, ylabel = y_label, yticks=True, x_ticks = x_ticks_c,
                hlines=hlines, title=title_c, ha_title=ha_title,
                )

    ticks_N = [f"2N={i}" for i in Ns]
    plot_bar_ax(ax_Ne, sum_roh*100, bins_cM, legend=True, 
                fs_l=fs_l, fs_y = fs_y, fs_x = fs_x, fs_t = fs_t, c = c, 
                x_ticks = ticks_N, ylim = ylim, 
                ylabel="", yticks=False,
                hlines=hlines, title=title_N, ha_title=ha_title)
            
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
        if np.sum(mod_idx)>0:
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

def prep_dfs_plot_exact(df, pops=[], col_group = "pop", 
                        mod_only=False, exact=True, cm_sort=4):
    """Pull all matching populations and return list of dataframes"""
    if mod_only:
        df = df[df["age"]==0]
    if exact:
        df_plots = [df[df[col_group].values==p] for p in pops]
    else:
        df_plots = [df[df[col_group].str.contains(p)] for p in pops]
    if cm_sort>0:
        df_plots = [df.sort_values(by=f"sum_roh>{cm_sort}", ascending=False) for df in df_plots]
    lgths = [len(df) for df in df_plots]
    return df_plots, lgths

def prep_xlabels(plot_dfs, col_age="age"):
    """Prepare vector of x labels
    corresponding to ages of individuals"""
    xlabels = []
    for df in plot_dfs:
        lbls = [str(int(x)) + " BP" for x in df[col_age].values]
        xlabels.append(lbls)
    return xlabels
