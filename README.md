
Multiple-Comparisons CC0 Toolkit (FWER/FDR) — with SSM Guardrail
Project Purpose:
A lightweight CC0 package for multiple-comparisons control (FWER/FDR), including Python CLI code, a static web app (HTML+JS), and ready-to-use code snippets for reports.
Built on the principles of OpenDecision-SSM: transparent decision-making, controlled conservatism, and a clear separation between confirmatory and exploratory inference.

What’s Included
Python CLI: Corrections via Bonferroni / Holm / Hochberg / Hommel / BH / BY; Tukey HSD; automatic method selection via decision flow.

Static Web App (web/): Browser-only calculator for BH / Holm / BY (no server required).

SSM Guardrail: Default N = 1000, Var_eff = 1 + 2/(N-3) and gamma > 1; variance is never reported as exactly 1.000.

Board One-Pager (export/board.md): Concise decision sheet for executive teams.

Examples (examples/): Sample CSV files for quick runs.

License: CC0 1.0 Universal.

Decision Flow (Summary)
FWER if false positive cost is high → Holm (default); under strong independence → Hochberg; complex dependence → Westfall–Young / Hommel.

FDR for exploratory research → BH (independence/PRDS), BY (arbitrary dependence), optionally IHW / Storey when covariates are available.

ANOVA: All pairwise comparisons → Tukey HSD; general contrasts → Scheffé.

Online/Sequential: Online FDR (SAFFRON / LORD) or alpha-spending.

Simultaneous Confidence Intervals: FWER-only (Bonferroni / Holm / Tukey / Scheffé).

See export/board.md for a detailed flowchart and examples.
