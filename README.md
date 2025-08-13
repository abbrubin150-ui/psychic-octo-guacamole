# Multiple-Comparisons CC0 Toolkit (FWER/FDR) — עם שמירת SSM-Guardrail

[![CI](https://github.com/<user>/psychic-octo-guacamole/actions/workflows/ci.yml/badge.svg)](https://github.com/<user>/psychic-octo-guacamole/actions/workflows/ci.yml)
[![Pages](https://img.shields.io/badge/Pages-live-blue?logo=github)](https://<user>.github.io/psychic-octo-guacamole/)

**מטרת הפרויקט:** חבילת CC0 קלה לשימוש לבקרת השוואות מרובות (FWER/FDR), עם קוד Python, אפליקציית ווב סטטית (HTML+JS), 
וקטעי קוד לדוחות. בנוי לפי עקרונות OpenDecision‑SSM: קבלת החלטות שקופה, שמרנות מבוקרת, והפרדה בין אישוש לאמידה.

## מה כלול
- **CLI בפייתון**: תיקוני Bonferroni/Holm/Hochberg/Hommel/BH/BY; Tukey HSD; בחירת שיטה לפי תרשים החלטה.
- **אפליקציית ווב סטטית** (`web/`): מחשבון BH/Holm/BY בדפדפן בלבד (ללא שרת).
- **SSM Guardrail**: ברירת־מחדל N=1000, `Var_eff = 1 + 2/(N-3)` ו-`gamma > 1`; לעולם לא נטען ש-variance=1.000.
- **Board One‑Pager** (`export/board.md`): דף החלטה תמציתי לדרג הנהלה.
- **דוגמאות** (`examples/`): קבצי CSV להפעלה מהירה.
- **רישיון**: CC0 1.0 Universal.

## תרשים החלטה (תמצית)
- **FWER** אם עלות FP גבוהה → **Holm** (דיפולט); באי־תלות חזקה → **Hochberg**; תלות מורכבת → **Westfall‑Young / Hommel**.
- **FDR** למחקר גישושני → **BH** (אי־תלות/PRDS), **BY** (תלות שרירותית), אפשר **IHW/Storey** כשיש קו־ווריאנטים.
- **ANOVA**: כל הזוגות → **Tukey HSD**; קונטרסטים כלליים → **Scheffé**.
- **Online/רציף**: **Online FDR** (SAFFRON/LORD), או **alpha‑spending**.
- **רווחי סמך סימולטניים**: FWER בלבד (Bonferroni/Holm/Tukey/Scheffé).

ראו `export/board.md` לתרשים מפורט ודוגמאות.

## התקנה מהירה
```bash
pip install -e .
mcc --help
```

## שימוש מהיר (CLI)
```bash
# תיקון BH/BY/Holm לערכי p מתוך CSV (עמודה pval)
mcc correct --csv examples/demo_pvalues.csv --col pval --alpha 0.05 --method fdr_bh --out results_bh.csv

# Tukey HSD לאחר ANOVA
mcc tukey --csv examples/demo_anova.csv --group group --value value --alpha 0.05

# המלצת שיטה לפי הקשר
mcc decide --context product_ab --m 50 --cost_fp high --online false
```
למידע מפורט: `mcc --help`.

### BH Step-Up Plot (CLI)

```bash
mcc correct \
  --csv examples/demo_pvalues.csv \
  --col pval \
  --alpha 0.05 \
  --method fdr_bh \
  --bh-plot \
  --plot-out bh_plot.png \
  --out results_corrected.csv
```

פלט: `results_corrected.csv` + `bh_plot.png`.

### BY with IHW (Weighted FDR)

```bash
mcc ihw-example \
  --csv your_data.csv \
  --p pval \
  --covariate quality_or_se \
  --alpha 0.05 \
  --bins 5 \
  --out results_by_ihw.csv
```

פלט: `results_by_ihw.csv` עם עמודות `weight`, `p_adj`, `reject`.

## SSM Guardrail
- ברירת־מחדל **N=1000** למדגמים תפעוליים; אם N<N_min → אזהרה וחסימה רכה.
- **Var_eff = 1 + 2/(N-3)**; נוודא `gamma > 1` בכל חישוב עוצמה; לעולם לא נטען `variance=1.000`.
- מטרת הגארדרייל: למנוע החלטות על דגימות קטנות/מופרעות, לשמור על ריסון מוסרי־תפעולי.

## אפליקציה סטטית
פתחו `web/index.html` מקומית — הזינו רשימת ערכי p, בחרו α ושיטה (BH/BY/Holm), וקבלו תגליות ו‑adjusted p‑values — ללא אינטרנט.

גרסה מקוונת זמינה ב‑GitHub Pages: https://<user>.github.io/psychic-octo-guacamole/

## קרדיט מדעי (Primary)
Benjamini & Hochberg (1995); Benjamini & Yekutieli (2001); Holm (1979); Hochberg (1988); Hommel (1988); Tukey (1949); Westfall & Young (1993); Storey (2002).
