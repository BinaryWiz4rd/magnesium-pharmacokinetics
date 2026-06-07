from __future__ import annotations

import math
import numpy as np

MW_MG   = 24.305
C_LOW   = 0.85
C_HIGH  = 1.10

def _lognormal(rng: np.random.Generator, mean: float, cv: float, size: int) -> np.ndarray:
    if cv == 0.0:
        return np.full(size, mean)
    sigma2 = math.log(1.0 + cv ** 2)
    mu     = math.log(mean) - 0.5 * sigma2
    return np.exp(rng.normal(mu, math.sqrt(sigma2), size))


def _normal_clipped(
    rng: np.random.Generator, mean: float, cv: float, size: int,
    lo: float = 0.0, hi: float = float("inf"),
) -> np.ndarray:
    return np.clip(rng.normal(mean, mean * cv, size), lo, hi)

def _sim_1c_iv_bolus(t, c_base, d_iv_mmol, vd, k_el):
    dc0 = d_iv_mmol / vd
    return c_base + dc0 * np.exp(-k_el * t)


def _sim_1c_iv_infusion(t, c_base, d_iv_mmol, vd, k_el, t_inf):
    r   = d_iv_mmol / t_inf
    css = r / (k_el * vd)
    c_during = c_base + css * (1.0 - np.exp(-k_el * t))
    c_end    = c_base + css * (1.0 - math.exp(-k_el * t_inf))
    c_after  = c_base + (c_end - c_base) * np.exp(-k_el * (t - t_inf))
    return np.where(t <= t_inf, c_during, c_after)


def _sim_1c_oral(t, c_base, d_oral_mmol, vd, k_el, f, ka):
    _ka = ka
    if abs(_ka - k_el) < 1e-6:
        _ka = k_el + 1e-5
    coef = (f * d_oral_mmol * _ka) / (vd * (_ka - k_el))
    return c_base + coef * (np.exp(-k_el * t) - np.exp(-_ka * t))


def _build_2c_params(k_el_total, v1, v2, q):
    k10 = k_el_total * (v1 + v2) / v1
    k12 = q / v1
    k21 = q / v2
    s   = k10 + k12 + k21
    disc = max(s ** 2 - 4.0 * k21 * k10, 0.0)
    alpha = (s + math.sqrt(disc)) / 2.0
    beta  = (s - math.sqrt(disc)) / 2.0
    denom = alpha - beta if abs(alpha - beta) > 1e-12 else 1e-12
    A = (alpha - k21) / denom
    B = (k21  - beta) / denom
    return alpha, beta, A, B, k21


def _sim_2c_iv_bolus(t, c_base, d_iv_mmol, v1, alpha, beta, A, B):
    scale = d_iv_mmol / v1
    return c_base + scale * (A * np.exp(-alpha * t) + B * np.exp(-beta * t))


def _sim_2c_iv_infusion(t, c_base, d_iv_mmol, v1, alpha, beta, A, B, t_inf):
    R0    = d_iv_mmol / t_inf
    scale = R0 / v1
    qa = (A / alpha) * (1.0 - math.exp(-alpha * t_inf))
    qb = (B / beta)  * (1.0 - math.exp(-beta  * t_inf))
    c_during = c_base + scale * (
        (A / alpha) * (1.0 - np.exp(-alpha * t)) +
        (B / beta)  * (1.0 - np.exp(-beta  * t))
    )
    c_after = c_base + scale * (
        qa * np.exp(-alpha * (t - t_inf)) +
        qb * np.exp(-beta  * (t - t_inf))
    )
    return np.where(t <= t_inf, c_during, c_after)


def _sim_2c_oral(t, c_base, d_oral_mmol, v1, alpha, beta, A_coef, B_coef, k21, f, ka):
    _ka = ka
    if abs(_ka - alpha) < 1e-6:
        _ka = alpha + 1e-5
    if abs(_ka - beta) < 1e-6:
        _ka = beta  + 1e-5
    P  = (k21 - alpha) / ((_ka - alpha) * (beta  - alpha))
    Qc = (k21 - beta)  / ((_ka - beta)  * (alpha - beta))
    R  = -(P + Qc)
    scale = f * d_oral_mmol * _ka / v1
    return c_base + scale * (
        P  * np.exp(-alpha * t) +
        Qc * np.exp(-beta  * t) +
        R  * np.exp(-_ka   * t)
    )

def run_monte_carlo(
    *,
    n_sims: int               = 1000,
    route: str                = "iv_infusion",
    t: np.ndarray,

    bw_mean: float            = 70.0,# kg
    t_half_mean: float        = 30.0,# h  (terminal half-life)
    vd_kg_mean: float         = 0.50,# L/kg  (1-comp Vd or total for 2-comp)
    dose_iv_mg: float         = 500.0,# mg elementary Mg
    dose_oral_mg: float       = 300.0,# mg elementary Mg
    t_inf: float              = 4.0,# h
    f_cit: float              = 0.31,
    ka_cit: float             = 0.90,# h⁻¹
    f_ox: float               = 0.055,
    ka_ox: float              = 0.35,# h⁻¹
    c_base: float             = 0.85,# mmol/L  baseline Mg

    # inter-individual variability (coefficient of variation, unitless)
    cv_bw: float              = 0.15,
    cv_t_half: float          = 0.35,
    cv_vd: float              = 0.20,
    cv_f: float               = 0.25,
    cv_ka: float              = 0.30,

    # two-compartment specific central values
    vc_kg_mean: float         = 0.20,    # L/kg  central compartment
    vp_kg_mean: float         = 0.30,    # L/kg  peripheral compartment
    q_kg_mean: float          = 0.12,    # L/h/kg  intercompartmental clearance
    cv_vc: float              = 0.20,
    cv_vp: float              = 0.20,
    cv_q: float               = 0.30,

    percentiles: tuple[int, ...] = (5, 50, 95),
    seed: int                 = 42,
) -> dict:
    rng = np.random.default_rng(seed)

    bw_s      = _normal_clipped(rng, bw_mean,    cv_bw,    n_sims, lo=40.0, hi=120.0)
    t_half_s  = _lognormal(rng,      t_half_mean, cv_t_half, n_sims)
    vd_kg_s   = _lognormal(rng,      vd_kg_mean,  cv_vd,    n_sims)

    two_comp = route.startswith("2c_")

    if two_comp:
        vc_kg_s = _lognormal(rng, vc_kg_mean, cv_vc, n_sims)
        vp_kg_s = _lognormal(rng, vp_kg_mean, cv_vp, n_sims)
        q_kg_s  = _lognormal(rng, q_kg_mean,  cv_q,  n_sims)

    oral_route = route in ("oral_cit", "oral_ox", "2c_oral_cit")
    if oral_route:
        f_mean  = f_cit  if route in ("oral_cit",  "2c_oral_cit") else f_ox
        ka_mean = ka_cit if route in ("oral_cit",  "2c_oral_cit") else ka_ox
        f_s  = np.clip(_lognormal(rng, f_mean,  cv_f,  n_sims), 0.01, 0.99)
        ka_s = _lognormal(rng, ka_mean, cv_ka, n_sims)

    curves = np.empty((n_sims, len(t)))

    for i in range(n_sims):
        bw    = bw_s[i]
        k_el  = math.log(2.0) / t_half_s[i]
        vd    = vd_kg_s[i] * bw
        d_iv  = dose_iv_mg   / MW_MG
        d_ora = dose_oral_mg / MW_MG

        if route == "iv_bolus":
            curves[i] = _sim_1c_iv_bolus(t, c_base, d_iv, vd, k_el)

        elif route == "iv_infusion":
            curves[i] = _sim_1c_iv_infusion(t, c_base, d_iv, vd, k_el, t_inf)

        elif route == "oral_cit":
            curves[i] = _sim_1c_oral(t, c_base, d_ora, vd, k_el, f_s[i], ka_s[i])

        elif route == "oral_ox":
            curves[i] = _sim_1c_oral(t, c_base, d_ora, vd, k_el, f_s[i], ka_s[i])

        elif two_comp:
            v1    = vc_kg_s[i] * bw
            v2    = vp_kg_s[i] * bw
            q_val = q_kg_s[i]  * bw
            alpha, beta, A, B, k21 = _build_2c_params(k_el, v1, v2, q_val)

            if route == "2c_iv_bolus":
                curves[i] = _sim_2c_iv_bolus(t, c_base, d_iv, v1, alpha, beta, A, B)

            elif route == "2c_iv_infusion":
                curves[i] = _sim_2c_iv_infusion(
                    t, c_base, d_iv, v1, alpha, beta, A, B, t_inf
                )

            elif route == "2c_oral_cit":
                curves[i] = _sim_2c_oral(
                    t, c_base, d_ora, v1, alpha, beta, A, B, k21,
                    f_s[i], ka_s[i],
                )
        else:
            raise ValueError(f"unknown route: {route!r}")

    result = {
        "t":          t,
        "route":      route,
        "n_sims":     n_sims,
        "curves":     curves,
        "percentiles": {p: np.percentile(curves, p, axis=0) for p in percentiles},
    }
    return result