# Multiple Comparisons — Board One‑Pager

**שאלת מפתח:** מהי *רמת הבקרה על טעויות* הדרושה לנו לפי עלות FP והקשר ניסויי?

## החלטה מהירה
- עלות FP **גבוהה** / אישושי → **FWER** → ברירת מחדל **Holm**; אם אי־תלות חזקה → **Hochberg**; תלות מורכבת → **Hommel/Westfall‑Young**.
- חקרני / m גדול → **FDR** → **BH** (אי־תלות/PRDS) או **BY** (תלות שרירותית); שקלו **IHW/Storey**.
- ANOVA: כל הזוגות → **Tukey HSD**; קונטרסטים כלליים → **Scheffé**.
- Online/רציף: **Online FDR** (SAFFRON/LORD) או **alpha‑spending**.
- CI סימולטניים: **FWER** בלבד.

## Guardrail (SSM)
- N תפעולי מינימלי: **1000**. אם N<1000 — עצירת־זהירות.
- Var_eff = **1 + 2/(N−3)**; הגדירו γ=Var_eff>1; לעולם לא מצהירים variance=1.000.

## תוצרים
- דוח תוצאות (CSV) כולל `p_adj` וסטטוס `reject`.
- גרף BH (אופציונלי) להצגה אינטואיטיבית.
