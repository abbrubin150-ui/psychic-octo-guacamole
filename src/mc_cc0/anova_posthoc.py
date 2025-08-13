from __future__ import annotations
import pandas as pd
from statsmodels.stats.multicomp import pairwise_tukeyhsd

def tukey_hsd(df: pd.DataFrame, value_col: str, group_col: str, alpha: float=0.05):
    res = pairwise_tukeyhsd(endog=df[value_col], groups=df[group_col], alpha=alpha)
    return res.summary()
