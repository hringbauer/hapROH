import logging

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

def call_roh(posterior0:np.ndarray, df_snp:pd.DataFrame, cutoff_post:float, snps_extend:int=0) -> pd.DataFrame:
    """Call ROH from the posterior probability.

    Args:
        posterior0: 1D np.ndarray of size nb_snp with the posterior probability of NOT being in ROH state
        df_snp: pandas.DataFrame with columns "pos" (BP) and "map" (Morgans)
        cutoff_post: minimal cutoff for the posterior probability to call ROH
        snps_extend: Number of snps to add to ROH blocks on both sides

    Returns:
    A dataframe with the ROH segments containing following columns:
        Start, End, length: SNP index
        StartBP, EndBP, lengthBP: position in basepairs
        StartM, EndM, lengthM: position in Morgans
    """

    roh = (1-posterior0) > cutoff_post
    roh = np.concatenate([[0], roh, [0]]).astype('int8')    # padding
    diff = np.diff(roh.astype('int8'))
    starts = np.flatnonzero(diff == 1)
    ends = np.flatnonzero(diff == -1) - 1

    starts = np.maximum(starts-snps_extend, 0)
    ends = np.minimum(ends+snps_extend, len(df_snp)-1)

    values = {"Start": starts, "End": ends, "length": ends-starts+1,
              "StartBP": df_snp["pos"].iloc[starts].values, "EndBP": df_snp["pos"].iloc[ends].values,
              "StartM": df_snp["map"].iloc[starts].values, "EndM": df_snp["map"].iloc[ends].values,
    }
    values["lengthBP"] = values["EndBP"] - values["StartBP"]
    values["lengthM"] = values["EndM"] - values["StartM"]

    return pd.DataFrame(values)

def merge_blocks(df_roh:pd.DataFrame, max_gap=0.005, min_len1=0.04, min_len2=0.02) -> pd.DataFrame:
    """Merge close adjacent ROH into longer ROH.

    Two adjacent ROH are merged if the gap between them (in Morgans) is at most
    max_gap, the longer of the two is at least min_len1, and the shorter of the
    two is at least min_len2.

    Args:
        df_roh: pandas.DataFrame of ROH segments, as returned by call_roh, with
            columns "Start", "End", "StartBP", "EndBP", "StartM", "EndM"
        max_gap: maximum gap (in Morgans) between two adjacent ROH for them to be merged
        min_len1: minimum length (in Morgans) required for the longer of the two adjacent ROH
        min_len2: minimum length (in Morgans) required for the shorter of the two adjacent ROH

    Returns:
    A dataframe with the merged ROH segments, with the same columns as df_roh:
        Start, End, length: SNP index
        StartBP, EndBP, lengthBP: position in basepairs
        StartM, EndM, lengthM: position in Morgans
    """
    if min_len1 < min_len2:
        min_len1, min_len2 = min_len2, min_len1

    if len(df_roh) == 0:
        return df_roh.copy()

    df_roh = df_roh[df_roh["lengthM"] >= min_len2]

    gap = df_roh["StartM"] - df_roh["EndM"].shift()
    len_1 = np.maximum(df_roh["lengthM"], df_roh["lengthM"].shift())
    len_2 = np.minimum(df_roh["lengthM"], df_roh["lengthM"].shift())

    merge = ((gap <= max_gap) & (len_1 >= min_len1) & (len_2 >= min_len2))[1:]
    group_id = np.concatenate([[0], np.cumsum(~merge)])

    out = df_roh.groupby(group_id).agg(
        Start=("Start", "first"),     End=("End", "last"),
        StartBP=("StartBP", "first"), EndBP=("EndBP", "last"),
        StartM=("StartM", "first"),   EndM=("EndM", "last"),
    ).reset_index(drop=True)

    out["length"]  = out["End"] - out["Start"] + 1
    out["lengthBP"] = out["EndBP"] - out["StartBP"]
    out["lengthM"]  = out["EndM"] - out["StartM"]

    logger.debug(f"Merged {len(df_roh)-len(out)} blocks")

    return out

def postproces(posterior0:np.ndarray, df_snp:pd.DataFrame,
               cutoff_post:float=0.999, snps_extend:int=0,
               max_gap=0.005, min_len1=0.04, min_len2=0.02):
    """Call ROH from the posterior probability and merge close adjacent ROH into longer ROH.

    Args:
        posterior0: 1D np.ndarray of size nb_snp with the posterior probability of NOT being in ROH state
        df_snp: pandas.DataFrame with columns "pos" (BP) and "map" (Morgans)
        cutoff_post: minimal cutoff for the posterior probability to call ROH
        snps_extend: Number of snps to add to ROH blocks on both sides
        max_gap: maximum gap (in Morgans) between two adjacent ROH for them to be merged
        min_len1: minimum length (in Morgans) required for the longer of the two adjacent ROH
        min_len2: minimum length (in Morgans) required for the shorter of the two adjacent ROH

    Returns:
    A dataframe with the merged ROH segments containing following columns:
        Start, End, length: SNP index
        StartBP, EndBP, lengthBP: position in basepairs
        StartM, EndM, lengthM: position in Morgans
    """
    df_roh = call_roh(posterior0, df_snp, cutoff_post, snps_extend)
    df_roh_merged = merge_blocks(df_roh, max_gap, min_len1, min_len2)
    return df_roh_merged

if __name__ == "__main__":
    segments = [(0,1), (1.2, 2), (2.01, 2.011)]
    segments = []

    df_roh = pd.DataFrame(segments, columns=["Start", "End"])
    df_roh["length"] = df_roh["End"] - df_roh["Start"]
    for column in ["Start", "End", "length"]:
        df_roh[f"{column}BP"] = 1 * df_roh[column]
        df_roh[f"{column}M"] = 1 * df_roh[column]
