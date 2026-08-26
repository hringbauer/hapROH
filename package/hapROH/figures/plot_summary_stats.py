import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def plot_summary_stats(
    df_stats: pd.DataFrame,
    L=(8, 12, 16, 20),
    L_colors=("#313695", "#abd9e9", "#fee090", "#d7191c"),
    legend: bool = True,
    x_ticks: None | str = None,
    y_ticks: bool = False,
    ax=None,
):
    """Plot the distribution of roh summary stats as a bar plot.
    Args:
        df_stats can be obtained by running create_stats(df_roh, L)
        L is the list of thresholds [in cM] considered (ie amount/length of roh >= x with x in L)
        legend: boolean, wether to plot legend or not
        x_ticks: if str, column of dataframe to use as x_ticks
        y_ticks: boolean wether to tick y axis or not"""

    if ax is None:
        fig, ax = plt.subplots()
    else:
        fig = ax.get_figure()

    # Prepare data
    df_stats = df_stats.sort_values(f"sum_ROH>{L[0]}", ascending=False)
    data = df_stats[[f"sum_ROH>{n}" for n in L]].values
    x = np.arange(len(df_stats))
    bottom = np.zeros((len(df_stats), len(L)))
    for i in range(1, len(L)):
        bottom[:, i] = data[:, 0] - data[:, i]

    # Make plot for each bin
    for i in range(len(L)):
        ax.bar(
            x,
            data[:, i],
            bottom=bottom[:, i],
            width=0.8,
            color=L_colors[i],
            edgecolor="black",
            label=f"{L[i]}-{L[i + 1]} cM" if i < len(L) - 1 else f">{L[i]} cM",
        )

    if legend:
        ax.legend(title="Sum of ROH in")
    if x_ticks:
        ax.set_xticks(x)
        ax.set_xticklabels(df_stats[x_ticks], rotation=270)
    else:
        ax.tick_params(labelbottom=False)
    ax.tick_params(axis="x", which="both", bottom=False, top=False)

    if y_ticks:
        ax.set_ylabel("Cumulative length of ROH [cM]")
    else:
        ax.tick_params(left=False, labelleft=False)
    ax.set_xlim(-1, len(df_stats))

    return fig, ax


def plot_summary_stats_panel(
    df_stats: list[pd.DataFrame],
    L=(8, 12, 16, 20),
    L_colors=("#313695", "#abd9e9", "#fee090", "#d7191c"),
    titles: None | list[str] = None,
    figsize=None,
    x_ticks: None | str = None,
):
    """Plot the distribution of multiple roh summary stats into different panels.
    Args:
        df_stats: list of dataframes containing the stats to plot. Each dataframe can be obtained by running create_stats(df_roh, L) on a different subset of data (eg different subspecies
        titles: list of titles for each panel
        For other arguments, refer to plot_summary_stats()"""
    if figsize is None:
        figsize = (6 * len(df_stats), 6)
    fig, axes = plt.subplots(
        1,
        len(df_stats),
        figsize=figsize,
        sharey=True,
        width_ratios=[len(df) for df in df_stats],
    )
    if titles is not None and len(df_stats) != len(titles):
        raise ValueError(
            f"Mismatch between the number of dataframes ({len(df_stats)}) and of titles ({len(titles)})"
        )
    for i, df in enumerate(df_stats):
        plot_summary_stats(
            df,
            L=L,
            L_colors=L_colors,
            legend=(i == len(df_stats) - 1),
            x_ticks=x_ticks,
            y_ticks=(i == 0),
            ax=axes[i],
        )
        title = titles[i] if titles is not None else ""
        axes[i].set_title(title)
    return fig, axes
