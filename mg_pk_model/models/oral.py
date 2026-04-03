import numpy as np
from core.derived import D_ORAL, VD, K_EL, C_BASE


def cp_oral(t: np.ndarray, F: float, ka: float) -> np.ndarray:
    """
    One-compartment oral absorption model
    """
    coef = (F * D_ORAL * ka) / (VD * (ka - K_EL))
    return C_BASE + coef * (np.exp(-K_EL * t) - np.exp(-ka * t))