from __future__ import annotations
import numpy as np
import pandas as pd
from typing import Literal, Tuple
from statsmodels.stats.multitest import multipletests

Method = Literal["bonferroni","holm","holm-sidak","sidak","hochberg","hommel","fdr_bh","fdr_by"]

def apply_correction(pvals: pd.Series | np.ndarray, alpha: float, method: Method="holm") -> Tuple[np.ndarray, np.ndarray]:
    p = np.asarray(pvals, dtype=float)
    reject, p_corr, _, _ = multipletests(p, alpha=alpha, method=method)
    return reject, p_corr

def bh_threshold(pvals: np.ndarray, q: float) -> float:
    p = np.sort(np.asarray(pvals, dtype=float))
    m = p.size
    crit = (np.arange(1, m+1) / m) * q
    mask = p <= crit
    if not mask.any():
        return 0.0
    k = np.max(np.where(mask)[0]) + 1
    return crit[k-1]
