from __future__ import annotations

import numpy as np
import pandas as pd
from typing import Literal, Tuple
from statsmodels.stats.multitest import multipletests

Method = Literal[
    "bonferroni",
    "holm",
    "holm-sidak",
    "sidak",
    "hochberg",
    "hommel",
    "fdr_bh",
    "fdr_by",
]


def apply_correction(
    pvals: pd.Series | np.ndarray, alpha: float, method: Method = "holm"
) -> Tuple[np.ndarray, np.ndarray]:
    p = np.asarray(pvals, dtype=float)
    reject, p_corr, _, _ = multipletests(p, alpha=alpha, method=method)
    return reject, p_corr


def ihw_weights(
    pvals: pd.Series | np.ndarray,
    covariate: pd.Series | np.ndarray,
    n_bins: int = 5,
    eps: float = 1e-12,
) -> np.ndarray:
    """Compute simple IHW-style weights by binning covariates.

    This is a lightweight approximation of the Independent Hypothesis
    Weighting (IHW) procedure. Hypotheses are split into quantile bins of
    the covariate, a weight inversely proportional to the mean p-value in
    each bin is assigned, and weights are scaled to sum to the number of
    hypotheses.
    """

    p = np.asarray(pvals, dtype=float)
    x = np.asarray(covariate, dtype=float)
    df = pd.DataFrame({"p": p, "x": x})
    df["bin"] = pd.qcut(df["x"], q=n_bins, duplicates="drop")
    mean_p = df.groupby("bin")["p"].mean() + eps
    df["weight"] = df["bin"].map(1.0 / mean_p).astype(float)
    # scale weights to sum to m
    df["weight"] *= len(df) / df["weight"].sum()
    return df["weight"].to_numpy()


def apply_ihw(
    pvals: pd.Series | np.ndarray,
    covariate: pd.Series | np.ndarray,
    alpha: float = 0.05,
    n_bins: int = 5,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Apply BY correction with IHW-derived weights.

    Parameters
    ----------
    pvals : array-like
        Raw p-values for each hypothesis.
    covariate : array-like
        Covariate used for weighting.
    alpha : float
        Target FDR level for BY correction.
    n_bins : int
        Number of quantile bins for the heuristic IHW weighting.

    Returns
    -------
    reject : ndarray
        Boolean array indicating rejected hypotheses.
    p_corr : ndarray
        BY-adjusted p-values after weighting.
    weights : ndarray
        Weights assigned to each hypothesis.
    """

    p = np.asarray(pvals, dtype=float)
    w = ihw_weights(p, covariate, n_bins=n_bins)
    p_weighted = p / w
    reject, p_adj, _, _ = multipletests(p_weighted, alpha=alpha, method="fdr_by")
    p_corr = np.minimum(p_adj * w, 1.0)
    return reject, p_corr, w


def bh_threshold(pvals: np.ndarray, q: float) -> float:
    p = np.sort(np.asarray(pvals, dtype=float))
    m = p.size
    crit = (np.arange(1, m + 1) / m) * q
    mask = p <= crit
    if not mask.any():
        return 0.0
    k = np.max(np.where(mask)[0]) + 1
    return crit[k - 1]
