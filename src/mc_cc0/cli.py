from __future__ import annotations
import argparse, sys, json
import numpy as np
import pandas as pd

# ננסה לייבא גם apply_ihw (אם קיים בסניף שלך). אם לא – נבנה fallback פנימי בהמשך.
from mc_cc0.corrections import apply_correction, bh_threshold
try:
    from mc_cc0.corrections import apply_ihw  # אופציונלי: בסניף codex
    _HAS_APPLY_IHW = True
except Exception:
    _HAS_APPLY_IHW = False

from mc_cc0.anova_posthoc import tukey_hsd
from mc_cc0.decision_flow import Context, recommend
from mc_cc0.ssm_guardrail import check_guardrail

METHODS = ["bonferroni","holm","holm-sidak","sidak","hochberg","hommel","fdr_bh","fdr_by"]

def _fallback_apply_ihw(pvals: pd.Series, covariate: pd.Series, alpha: float = 0.05,
                        n_bins: int = 5, higher_is_better: bool = True):
    """
    Fallback מינימלי ל-BY+IHW אם apply_ihw לא קיים במודול corrections.
    מבצע:
      1) חלוקה ל-bins לפי קוונטים של הקו-ווריאנט
      2) משקולות לשכבות לפי ממוצע cov (גבוה=חזק אם higher_is_better=True; אחרת הפוך)
      3) נרמול sum(weights)=m
      4) p_weighted = p / w
      5) BY על p_weighted
    """
    from statsmodels.stats.multitest import multipletests
    p = np.asarray(pvals, dtype=float)
    cov = np.asarray(covariate, dtype=float)
    m = p.size
    # קוונטים
    qs = np.linspace(0, 1, n_bins + 1)
    bins = np.unique(np.quantile(cov, qs))
    if bins.size < 2:
        # דגנרטיבי – כל הקו-ווריאנטים זהים; שקול משקל אחיד
        weights = np.ones_like(p)
    else:
        strata = np.digitize(cov, bins[1:-1], right=True)
        # ממוצע cov לכל שכבה
        df = pd.DataFrame({"cov": cov, "stratum": strata})
        means = df.groupby("stratum")["cov"].mean().rename("mean_cov").reset_index()
        if higher_is_better:
            means["raw_w"] = (means["mean_cov"] - means["mean_cov"].min()) + 1e-8
        else:
            means["raw_w"] = (means["mean_cov"].max() - means["mean_cov"]) + 1e-8
        means["raw_w"] = means["raw_w"].clip(lower=1e-8)

        # גודל שכבה
        counts = df["stratum"].value_counts().sort_index()
        # סך משקל "לפני נרמול" = sum_j raw_w[j] * size_j
        total_raw = float((means["raw_w"] * counts.values).sum())
        scale = m / total_raw if total_raw > 0 else 1.0
        means["w_per_hypothesis"] = means["raw_w"] * scale
        w_map = dict(zip(means["stratum"], means["w_per_hypothesis"]))
        weights = np.array([w_map.get(s, 1.0) for s in strata], dtype=float)

    # בטיחות
    weights = np.clip(weights, 1e-8, None)
    weights *= (m / weights.sum())  # נרמול סופי

    p_weighted = p / weights
    reject, p_adj, _, _ = multipletests(p_weighted, alpha=alpha, method="fdr_by")
    return reject, p_adj, weights

def main():
    p = argparse.ArgumentParser(prog="mcc", description="Multiple-Comparisons CC0 Toolkit")
    sub = p.add_subparsers(dest="cmd", required=True)

    # correct
    c = sub.add_parser("correct", help="Apply multiple-comparisons correction")
    c.add_argument("--csv", required=True)
    c.add_argument("--col", required=True, help="Column name containing raw p-values")
    c.add_argument("--alpha", type=float, default=0.05)
    c.add_argument("--method", choices=METHODS, default="holm")
    c.add_argument("--out", default="results_corrected.csv")
    c.add_argument("--n", type=int, default=1000, help="Sample size for guardrail check")
    # דגלי תרשים BH
    c.add_argument("--bh-plot", action="store_true", help="Produce BH step-up plot (when using fdr_bh)")
    c.add_argument("--plot-out", default=None, help="Output path for BH plot image (e.g., bh_plot.png)")

    # tukey
    t = sub.add_parser("tukey", help="Tukey HSD post-hoc after ANOVA")
    t.add_argument("--csv", required=True)
    t.add_argument("--group", required=True)
    t.add_argument("--value", required=True)
    t.add_argument("--alpha", type=float, default=0.05)
    t.add_argument("--n", type=int, default=1000)

    # decide
    d = sub.add_parser("decide", help="Recommend a method from context")
    d.add_argument("--purpose", choices=["confirmatory","exploratory"], required=True)
    d.add_argument("--m", type=int, required=True)
    d.add_argument("--independence", choices=["independent","prds","arbitrary","unknown"], default="unknown")
    d.add_argument("--online", action="store_true")
    d.add_argument("--anova", choices=["none","pairwise","contrasts"], default="none")
    d.add_argument("--cost-fp", choices=["high","medium","low"], default="high")

    # ihw-example
    ihw = sub.add_parser("ihw-example", help="Illustrative IHW + BY correction")
    ihw.add_argument("--csv", required=True)
    ihw.add_argument("--p", required=True, help="Column with raw p-values")
    ihw.add_argument("--covariate", required=True, help="Covariate column for weighting")
    ihw.add_argument("--alpha", type=float, default=0.05)
    ihw.add_argument("--bins", type=int, default=5, help="Number of strata (quantile bins) for IHW")
    ihw.add_argument("--higher-is-better", action="store_true",
                     help="If covariate higher implies higher power; default: False (lower is better)")
    ihw.add_argument("--out", default="results_ihw.csv")

    args = p.parse_args()

    if args.cmd == "correct":
        guard = check_guardrail(args.n)
        df = pd.read_csv(args.csv)
        pvals = df[args.col]
        reject, p_corr = apply_correction(pvals, args.alpha, method=args.method)
        out = df.copy()
        out["reject"] = reject
        out["p_adj"] = p_corr
        out.to_csv(args.out, index=False)

        result = {"guardrail": guard, "out": args.out, "method": args.method, "alpha": args.alpha}

        # תרשים BH – רק כאשר method=fdr_bh והמשתמש ביקש
        if args.bh_plot and args.method == "fdr_bh":
            import matplotlib.pyplot as plt
            p_sorted = np.sort(pvals.to_numpy(dtype=float))
            m = p_sorted.size
            k = np.arange(1, m + 1)
            crit = (k / m) * args.alpha

            # cutoff לפי BH
            mask = p_sorted <= crit
            k_cut = (np.max(np.where(mask)[0]) + 1) if mask.any() else 0
            y_cut = (k_cut / m) * args.alpha if k_cut > 0 else None

            fig, ax = plt.subplots(figsize=(7, 5))
            ax.scatter(k, p_sorted, label="sorted p-values")
            ax.plot(k, crit, linestyle="--", label=f"BH line (q={args.alpha})")
            if k_cut > 0:
                ax.scatter([k_cut], [p_sorted[k_cut - 1]], s=60, zorder=5, label=f"cutoff k={k_cut}")
                ax.axvline(k_cut, color="red", linestyle="--", linewidth=1)
                ax.axhline(y_cut, color="red", linestyle="--", linewidth=1)
            ax.set_xlabel("rank k")
            ax.set_ylabel("p-value")
            ax.set_title("Benjamini–Hochberg step-up")
            ax.legend()
            ax.grid(True, linestyle=":", alpha=0.6)

            plot_path = args.plot_out or (args.out.rsplit(".", 1)[0] + "_bh_plot.png")
            fig.savefig(plot_path, dpi=150, bbox_inches="tight")
            result["bh_plot"] = plot_path

        print(json.dumps(result, ensure_ascii=False))
        return

    if args.cmd == "tukey":
        guard = check_guardrail(args.n)
        df = pd.read_csv(args.csv)
        print(json.dumps(guard, ensure_ascii=False))
        print(tukey_hsd(df, args.value, args.group, alpha=args.alpha))
        return

    if args.cmd == "decide":
        ctx = Context(args.purpose, args.m, args.independence, args.online, args.anova, args.cost_fp)
        print(recommend(ctx))
        return

    if args.cmd == "ihw-example":
        df = pd.read_csv(args.csv)
        pvals = df[args.p]
        covariate = df[args.covariate]

        if _HAS_APPLY_IHW:
            # שימוש במימוש שקיים אצלך במודול corrections (אם יש)
            reject, p_corr, weights = apply_ihw(pvals, covariate, alpha=args.alpha,
                                                n_bins=args.bins, higher_is_better=args.higher_is_better)
        else:
            # fallback פנימי
            reject, p_corr, weights = _fallback_apply_ihw(
                pvals, covariate, alpha=args.alpha,
                n_bins=args.bins, higher_is_better=args.higher_is_better
            )

        out = df.copy()
        out["weight"] = weights
        out["p_adj"] = p_corr
        out["reject"] = reject
        out.to_csv(args.out, index=False)
        print(json.dumps({"out": args.out, "alpha": args.alpha, "bins": args.bins,
                          "higher_is_better": args.higher_is_better}, ensure_ascii=False))
        return

if __name__ == "__main__":
    main()
