"""
Plotting function to plot IBD X Estimates from Ne
@ Author: Harald Ringbauer, 2020, All rights reserved
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import itertools as it
from scipy.special import gammaincinv
from hapsburg.PackagesSupport.roh_expectations import Expected_Roh


def new_column(df, df_meta, col="New Clade", col_new="", match_col="iid"):
    """Maps Entries from meta dataframe onto the IBD dataframe.
    Return modified dataframe"""
    if len(col_new)==0:
        col_new=col
    dct = pd.Series(df_meta[col].values, index=df_meta[match_col]).to_dict()
    
    for i in range(1,3):    
        df[col_new + str(i)] =df["iid" + str(i)].map(dct)
        df[col_new + str(i)] =df["iid" + str(i)].map(dct)
    return df

def give_sub_df(df, pop1="La Caleta", pop2="La Caleta", col="clst", 
                output=True, exact=False):
    """Return sub dataframe where pair across pop1 and pop2"""
    if exact:
        idx1 = (df[col+"1"]==pop1) & (df[col + "2"]==pop2)
        idx2 = (df[col+"1"]==pop2) & (df[col + "2"]==pop1)
    else:
        idx1 = (df[col+"1"].str.contains(pop1) & df[col + "2"].str.contains(pop2))
        idx2 = (df[col+"1"].str.contains(pop2) & df[col + "2"].str.contains(pop1))
        
    idx = idx1 | idx2
    if output:
        print(f"Found: {np.sum(idx)}\n")
    return df[idx]

def give_stats_cm_bin(df, cms = [4,6,8,10,12], col_b="n_roh>", output=True):
    """print some stats about IBD dataframe df
    col_b: Start of IBD/ROH column. If one of cms 0 use open
    upper bound"""
    counts = np.zeros(len(cms)-1)
    for i in range(len(cms)-1):
        cm1, cm2 = cms[i], cms[i+1]
        col1, col2 = f"{col_b}{cm1}", f"{col_b}{cm2}"
        if cm2>0:
            cts = df[col1].values-df[col2].values
        else:
            cts = df[col1].values
        tot_cts = np.sum(cts)
        counts[i] = tot_cts
        frac = np.mean(cts)
        if output:
            idx = df["sum_roh>4"]==0
            print(f"No ROH found: {np.sum(idx)}/{len(idx)}")
            print(f"#ROH in {cm1}-{cm2}: {tot_cts}/{len(cts)}: {frac*100:.2f}%")
    return counts, len(cts)


#################################################################################
#### Plotting Functions

def get_IBD_stats(df, pop1="", pop2="", col="clade", col_b="n_roh>",
                  cms=[4,6,8,10,12], output=False, a=0.05):
    """Get IBD fraction statistics.
    a: Signficance level
    Return fractions, confidence intervalls as well as 
    number of pairsise comparisons
    col_b: Start of IBD/ROH column"""
    if len(pop1)>0 or len(pop2)>0:
        df1 = give_sub_df(df, pop1=pop1, pop2=pop2, col=col, output=output)
    else:
        df1 = df
    
    counts, n = give_stats_cm_bin(df1, output=output, cms=cms, col_b=col_b)
    fracs = counts/n
    cis = get_ci_counts(counts, n, a=a)
    #stds = np.sqrt(counts)/n
    return fracs, cis, n

def get_IBD_stats_pops(df, pops1=[], pops2=[], col="clade", col_b="n_roh>", 
                       cms=[4,6,8,10,12], output=False, a=0.05):
    """Get IBD fraction statistics for list of pop pairs.
    Returns lists of fractions, confidence intervalls as well as 
    number of pairsise comparisons.
    a: Significance Level
    col_b: Start of IBD/ROH column"""
    assert(len(pops1)==len(pops2)) # Sanity Check
    fracss, ciss, ns = [], [], []
    
    for i in range(len(pops1)):
        fracs, cis, n = get_IBD_stats(df, pop1=pops1[i], pop2=pops2[i], col=col,
                                      col_b=col_b, cms=cms, output=output, a=a)
        fracss.append(fracs)
        ciss.append(cis)
        ns.append(n)
    return fracss, ciss, ns

def get_ci_counts(counts, n, a=0.05, minc=1e-4):
    """Get Confidence Intervalls from counts and
    trials. Return list of CIS (lentght 2 each)
    counts: Array of Counts
    n: Total number of trials
    a: Signficance level"""
    cis = []
    for k in counts:
        if k==0:  # In case of no Count
            c0=minc/n
        else:
            c0 = gammaincinv(k, 0.5 * a)     # Lowe Bound
        c1 = gammaincinv(k + 1, 1 - 0.5 * a) # Upper bound
        cis.append([c0/n,c1/n])
    cis = np.array(cis) # Transform to numpy
    return cis

#########################################
#### Function to calculate epxected ne

def create_Ne_roh_nr(Ns=[400, 800, 1600, 3200, 6400],  
                  bins=[[0.04,0.08],[0.08,0.12],[0.12,0.2],[0.2,3.0]], 
                  bin_n=10000, chr_lgts=[]):
    """Create ROH sharing in list of bins (list of [begin,end]) 
    for panmictic population sizes
    Ns: List of population sizes
    bins: Length Bins (in Morgan) to calculate expectations from
    return sharing [len(degrees), len(bins)]"""
    n_roh = np.zeros((len(Ns),len(bins))) # Container for results Cousins
    e_roh = Expected_Roh()
    for i,N in enumerate(Ns):
        for j,b in enumerate(bins):
            x_arr = np.linspace(b[0], b[1], bin_n)
            bw = x_arr[1] - x_arr[0]
            y = e_roh.roh_pdf_allchr_N(x_arr + bw/2, chr_lgts=chr_lgts, N=N) * bw
            n_roh[i,j] = np.sum(y)
            #n_roh[i,j] = e_roh.exp_roh_len_in_bin_N(b, N=N, bins=bin_n, chr_lgts=chr_lgts)
    return n_roh

def create_ne_dct(cms, ns, chr_lgts=[1.8085 * 2/3,]):
    """Create Dictionary of Ne"""
    bins = [cms[i:i+2]/100 for i in range(len(cms)-1)]  # Prepare bin vector
    roh_n = create_Ne_roh_nr(Ns=ns, bins=bins, bin_n=10000,
                          chr_lgts=chr_lgts)
    dct = {}
    for i in range(len(ns)):
        dct[ns[i]] = roh_n[i]
    return dct

########################################################
### Plotting Functions

def plot_IBD_fracs(ax, cms1, fracs, cis=[], n=0, c="maroon",
                   xlim=[4,12], ylim=[1e-3, 1e0], fs=(6,6),
                   lw=2, ms=6, yscale="log", show=False,
                   hlines=[], lw_h=0.8, label=""):
    """Plot IBD fractions in cM bins. Plots on axis.
    ax: Axis to plot on.
    cms1: Vector of cM mean bins.
    fracs: Fractions to plot on y
    cis: Confidence Intervals [n elements of length 2]"""
    ax.plot(cms1, fracs, "o-", lw=lw, c=c, label=label)
    if len(cis)>0:
        #ax.errorbar(cms1, fracs, yerr=sts, fmt='o', c=c, ms=ms)
        ax.vlines(cms1, ymin=cis[:,0], ymax=cis[:,1], color=c)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_yscale(yscale)
    
    for h in hlines:
        ax.axhline(h, linestyle="--", c="gray", lw=lw_h)
        
    if n>0:
        #ax.axhline(1/n, linestyle="-", c="k", lw=0.5)
        ax.fill_between(cms1,0,1/n, color="lightgray")
    
    if show:
        plt.show()
        
def plot_ibd_curves(figsize=(6,6), labels=[], cms=[], 
                    fracs=[], cis=[], fs_big=16, fs_l=12, 
                    fs_t=12, label_idx=-1,
                    xlim=[4,12], ylim=[0.5e-2, 1.1e0], colors=[],
                    ylabel='# IBD per pair of male X', 
                    xlabel='centimorgan [2 cM bins]',
                    savepath="", dpi=400, leg_loc="upper right",
                    hlines=[1e-3,1e-2,1e-1,1e0], offset=0,
                    dct_ne={}, frameon=False):
    """Plot IBD X curves on one plot.
    offset: Offset along x axis of each y value.
    label_idx: If>=0 Plot Ne labels on curves on according bin"""
    cms1 = (cms[1:] + cms[:-1])/2
    
    fig = plt.figure(figsize=figsize)
    ax=plt.gca()

    for i in range(len(labels)):
        cms2 = cms1 + offset*i - ((len(labels)-1)/2*offset)
        plot_IBD_fracs(ax, cms2, fracs[i], cis[i], n=0,  
                       xlim=xlim, ylim=ylim, show=False,
                       label=labels[i],c=colors[i], hlines=hlines)

    ax.legend(loc=leg_loc, fontsize=fs_l, frameon=frameon)
    if len(xlabel)>0:
        fig.text(0.5, 0.02, xlabel, ha='center', fontsize=fs_big)
    if len(ylabel)>0:
        fig.text(0.03, 0.5, ylabel, va='center', rotation='vertical', fontsize=fs_big)
    
    if len(dct_ne)>0:
        plot_ne_dct(ax, dct_ne, bins=cms, 
                    label_idx=label_idx, fs=fs_t)

    if len(savepath)>0:
        plt.savefig(savepath, bbox_inches = 'tight', pad_inches = 0, dpi=dpi)
        print(f"Saved to {savepath}")
    plt.show()
    
def plot_ne_dct(ax, dct_ne={}, bins=[], label_idx=-1, c="lightgray", fs=12):
    """Plot ne: ROH sharing onto axis"""
    cms1 = (bins[1:]+bins[:-1])/2 # cm Means
    for i,n in enumerate(dct_ne):
        roh_nr = dct_ne[n]
        label = str(int(n/2))
        if i==0:
            label = f"$N_e=$"+label
        
        assert(len(cms1)==len(roh_nr))
        ax.plot(cms1, roh_nr, label=label, 
                zorder=0, c=c)
        if label_idx>=0:
            ax.text(cms1[label_idx], roh_nr[label_idx], label, 
                    ha='left', va='bottom', fontsize=fs, 
                    zorder=0, color="silver")