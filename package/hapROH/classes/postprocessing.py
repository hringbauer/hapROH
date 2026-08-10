import numpy as np
import pandas as pd

def call_roh(posterior0:np.ndarray, df_snp:pd.DataFrame, cutoff_post:float) -> pd.DataFrame:
    """Call ROH from the posterior probability.
    Returns a dataframe with columns"""

    roh = (1-posterior0) > cutoff_post
    roh = np.concatenate([[0], roh, [0]]).astype('int8')    # padding
    diff = np.diff(roh.astype('int8'))
    starts = np.flatnonzero(diff == 1)
    ends = np.flatnonzero(diff == -1) - 1

    values = {"Start": starts, "End": ends, "length": ends-starts+1,
              "StartBP": df_snp["pos"].iloc[starts].values, "EndBP": df_snp["pos"].iloc[ends].values,
              "StartM": df_snp["map"].iloc[starts].values, "EndM": df_snp["map"].iloc[ends].values,
    }
    values["lengthBP"] = values["EndBP"] - values["StartBP"]
    values["lengthM"] = values["EndM"] - values["StartM"]

    return pd.DataFrame(values)

def merge_called_blocks(df_roh:pd.DataFrame) -> pd.DataFrame:
    """Merge close adjacent ROH into longer ROH"""
    raise NotImplementedError()

def postproces(posterior0:np.ndarray, df_snp:pd.DataFrame, cutoff_post:float=0.999):
    """Call roh"""
    df_roh = call_roh(posterior0, df_snp, cutoff_post)
    # df_roh_merged = merge_called_blocks(df_roh)
    return df_roh