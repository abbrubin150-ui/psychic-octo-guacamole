from __future__ import annotations
import argparse, sys, json
import pandas as pd
from mc_cc0.corrections import apply_correction, apply_ihw
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

    ihw = sub.add_parser("ihw-example", help="Illustrative IHW + BY correction")
    ihw.add_argument("--csv", required=True)
    ihw.add_argument("--p", required=True, help="Column with raw p-values")
    ihw.add_argument("--covariate", required=True, help="Covariate column for weighting")
    ihw.add_argument("--alpha", type=float, default=0.05)
    ihw.add_argument("--out", default="results_ihw.csv")

    args = p.parse_args()

    if args.cmd == "correct":
        guard = check_guardrail(args.n)
        df = pd.read_csv(args.csv)
        reject, p_corr = apply_correction(df[args.col], args.alpha, method=args.method)
        out = df.copy()
        out["reject"] = reject
        out["p_adj"] = p_corr
        out.to_csv(args.out, index=False)
        print(json.dumps({"guardrail": guard, "out": args.out, "method": args.method, "alpha": args.alpha}, ensure_ascii=False))
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

    if args.cmd == "ihw-example":
        df = pd.read_csv(args.csv)
        reject, p_corr, weights = apply_ihw(df[args.p], df[args.covariate], alpha=args.alpha)
        out = df.copy()
        out["weight"] = weights
        out["p_adj"] = p_corr
        out["reject"] = reject
        out.to_csv(args.out, index=False)
        print(json.dumps({"out": args.out, "alpha": args.alpha}, ensure_ascii=False))
        return

if __name__ == "__main__":
    main()
