from __future__ import annotations
from dataclasses import dataclass

@dataclass
class Context:
    purpose: str  # 'confirmatory' or 'exploratory'
    m: int
    independence: str  # 'independent','prds','arbitrary','unknown'
    online: bool
    anova: str  # 'none','pairwise','contrasts'
    cost_fp: str  # 'high','medium','low'

def recommend(context: Context) -> str:
    if context.anova == "pairwise":
        return "tukey_hsd"
    if context.anova == "contrasts":
        return "scheffe"
    if context.online:
        return "online_fdr"
    if context.purpose == "confirmatory":
        if context.independence == "independent":
            return "hochberg"
        if context.independence in ("unknown","arbitrary","prds"):
            return "holm"
        return "holm"
    # exploratory
    if context.independence in ("independent","prds"):
        return "fdr_bh"
    return "fdr_by"
