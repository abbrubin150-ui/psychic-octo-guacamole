from __future__ import annotations
import argparse, sys, json
import numpy as np
import pandas as pd
from mc_cc0.corrections import apply_correction, bh_threshold
from mc_cc0.anova_posthoc import tukey_hsd
from mc_cc0.decision_flow import Context, recommend
from mc_cc0.ssm_guardrail import check_guardrail

METHODS = ["bonferroni","holm","holm-sidak","sidak","hochberg","hommel","fdr_bh","fdr_by"]

def main():
    p = argparse.ArgumentParser(prog="mcc", description="Multiple-Comparisons CC0 Toolkit")
    sub = p.add_subparsers(dest="cmd", required=True)

    c = sub.add_parser("correct", help="Apply multiple-comparisons correction")
    c.add_argument("--csv", required=True)
    c.add_argument("--col", required=True, help="Column name containing raw p-values")
    c.add_argument("--alpha", type=float, default=0.05)
    c.add_argument("--method", choices=METHODS, default="holm")
    c.add_argument("--out", default="results_corrected.csv")
    c.add_argument("--n", type=int, default=1000, help="Sample size for guardrail check")
    c.add_argument("--plot", action="store_true", help="Generate BH diagnostic plot")
    
    t = sub.add_parser("tukey", help="Tukey HSD post-hoc after ANOVA")
    t.add_argument("--csv", required=True)
    t.add_argument("--group", required=True)
    t.add_argument("--value", required=True)
    t.add_argument("--alpha", type=float, default=0.05)
    t.add_argument("--n", type=int, default=1000)

    d = sub.add_parser("decide", help="Recommend a method from context")
    d.add_argument("--purpose", choices=["confirmatory","exploratory"], required=True)
    d.add_argument("--m", type=int, required=True)
    d.add_argument("--independence", choices=["independent","prds","arbitrary","unknown"], default="unknown")
    d.add_argument("--online", action="store_true")
    d.add_argument("--anova", choices=["none","pairwise","contrasts"], default="none")
    d.add_argument("--cost-fp", choices=["high","medium","low"], default="high")

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
        plot_path = None
        if args.plot:
            import matplotlib.pyplot as plt
            bh = bh_threshold(pvals.to_numpy(), args.alpha)
            p_sorted = np.sort(pvals.to_numpy())
            m = p_sorted.size
            k = np.arange(1, m + 1)
            crit = (k / m) * args.alpha
            fig, ax = plt.subplots()
            ax.scatter(k, p_sorted, label="p-values")
            ax.plot(k, crit, label="k/m·q")
            mask = p_sorted <= crit
            if mask.any():
                k_cut = np.max(np.where(mask)[0]) + 1
                ax.scatter(k_cut, p_sorted[k_cut-1], color="red", zorder=5, label="cutoff")
                ax.axvline(k_cut, color="red", linestyle="--")
                ax.axhline(bh, color="red", linestyle="--")
            ax.set_xlabel("k")
            ax.set_ylabel("p-value")
            ax.legend()
            plot_path = args.out.rsplit('.', 1)[0] + '.png'
            fig.savefig(plot_path, dpi=150, bbox_inches="tight")
            result["plot"] = plot_path
        print(json.dumps(result, ensure_ascii=False))
        if plot_path:
            print(f"![BH plot]({plot_path})")
        return

    if args.cmd == "tukey":
        guard = check_guardrail(args.n)
        df = pd.read_csv(args.csv)
        print(guard)
        print(tukey_hsd(df, args.value, args.group, alpha=args.alpha))
        return

    if args.cmd == "decide":
        ctx = Context(args.purpose, args.m, args.independence, args.online, args.anova, args.cost_fp)
        print(recommend(ctx))
        return

if __name__ == "__main__":
    main()
