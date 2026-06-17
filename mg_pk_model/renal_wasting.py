import numpy as np


def simulate_renal_wasting_patient(condition, route, dose_mg, duration_h, weight_kg):
    """
    Two-compartment magnesium PK simulation with renal wasting assumptions.
    Concentrations are returned in mmol/L over a 48-hour period.
    """
    if condition == "Healthy":
        c_base = 0.85
        fe_factor = 1.0
    elif condition == "Gitelman Syndrome":
        c_base = 0.55
        fe_factor = 2.5
    elif condition == "Bartter Syndrome Type 4":
        c_base = 0.45
        fe_factor = 7.0
    elif condition == "EAST/SeSAME Syndrome":
        c_base = 0.50
        fe_factor = 5.0
    elif condition == "CNNM2-Related":
        c_base = 0.40
        fe_factor = 6.0
    else:
        c_base = 0.85
        fe_factor = 1.0

    v1 = 0.20 * weight_kg
    v2 = 0.30 * weight_kg
    q = 0.12 * weight_kg

    cl_base = 0.0115 * weight_kg
    cl = cl_base * fe_factor

    k10 = cl / v1
    k12 = q / v1
    k21 = q / v2

    sum_k = k10 + k12 + k21
    prod_k = k21 * k10
    discriminant = np.sqrt(max(0, sum_k**2 - 4 * prod_k))

    alpha = (sum_k + discriminant) / 2.0
    beta = (sum_k - discriminant) / 2.0

    a_coeff = (alpha - k21) / (alpha - beta) if (alpha - beta) != 0 else 1.0
    b_coeff = 1.0 - a_coeff

    dose_mmol = dose_mg / 24.305
    t = np.linspace(0, 48, 481)
    c_plasma = np.zeros_like(t)

    if route == "IV Bolus":
        c_dynamic = (dose_mmol / v1) * (
            a_coeff * np.exp(-alpha * t)
            + b_coeff * np.exp(-beta * t)
        )
        c_plasma = c_base + c_dynamic

    elif route == "IV Infusion":
        r0 = dose_mmol / duration_h

        q_alpha = (a_coeff / alpha) * (1.0 - np.exp(-alpha * duration_h))
        q_beta = (b_coeff / beta) * (1.0 - np.exp(-beta * duration_h))

        mask_during = t <= duration_h
        mask_after = t > duration_h

        c_plasma[mask_during] = c_base + (r0 / v1) * (
            (a_coeff / alpha) * (1.0 - np.exp(-alpha * t[mask_during]))
            + (b_coeff / beta) * (1.0 - np.exp(-beta * t[mask_during]))
        )
        c_plasma[mask_after] = c_base + (r0 / v1) * (
            q_alpha * np.exp(-alpha * (t[mask_after] - duration_h))
            + q_beta * np.exp(-beta * (t[mask_after] - duration_h))
        )

    elif route in ["Oral Citrate", "Oral Oxide"]:
        if route == "Oral Citrate":
            bioavailability = 0.31
            ka = 0.90
        else:
            bioavailability = 0.055
            ka = 0.35

        if abs(ka - alpha) < 1e-5:
            ka += 1e-5
        if abs(ka - beta) < 1e-5:
            ka += 1e-5

        p_coeff = (k21 - alpha) / ((ka - alpha) * (beta - alpha))
        q_coeff = (k21 - beta) / ((ka - beta) * (alpha - beta))
        r_coeff = -(p_coeff + q_coeff)

        c_dynamic = ((bioavailability * dose_mmol * ka) / v1) * (
            p_coeff * np.exp(-alpha * t)
            + q_coeff * np.exp(-beta * t)
            + r_coeff * np.exp(-ka * t)
        )
        c_plasma = c_base + c_dynamic

    else:
        c_plasma = np.full_like(t, c_base)

    return t, c_plasma
