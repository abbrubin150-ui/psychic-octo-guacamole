from __future__ import annotations

def var_eff(n: int) -> float:
    if n <= 3:
        return float("inf")
    return 1.0 + 2.0/(n-3)

def check_guardrail(n: int, n_min: int=1000) -> dict:
    ve = var_eff(n)
    gamma = max(1.0000001, ve)  # enforce gamma>1 and never claim variance=1.000
    ok = n >= n_min
    return {
        "ok": ok,
        "n": n,
        "n_min": n_min,
        "var_eff": ve,
        "gamma": gamma,
        "message": "OK" if ok else f"N too small for operational guardrail (need >= {n_min})"
    }
