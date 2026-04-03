import numpy as np
import math
from core.derived import D_IV, VD, K_EL, C_BASE, T_INF


def cp_iv_bolus(t: np.ndarray) -> np.ndarray:
    """
    One-compartment IV bolus model
    """
    dC0 = D_IV / VD
    return C_BASE + dC0 * np.exp(-K_EL * t)


def cp_iv_infusion(t: np.ndarray) -> np.ndarray:
    """
    One-compartment IV infusion model
    """
    R   = D_IV / T_INF
    Css = R / (K_EL * VD)

    c_during = C_BASE + Css * (1 - np.exp(-K_EL * t))
    C_end    = float(C_BASE + Css * (1 - math.exp(-K_EL * T_INF)))
    c_after  = C_BASE + (C_end - C_BASE) * np.exp(-K_EL * (t - T_INF))

    return np.where(t <= T_INF, c_during, c_after)