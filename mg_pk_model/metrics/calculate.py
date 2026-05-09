import numpy as np

def calculate_pk_metrics(t: np.ndarray, cp: np.ndarray) -> dict:
    """
    Calculates Cmax, Tmax, and AUC for a given concentration-time profile.
    """
    tmax_idx = np.argmax(cp)
    cmax = cp[tmax_idx]
    tmax = t[tmax_idx]

    # Area Under the Curve (AUC) using the trapezoidal rule
    auc = np.trapz(cp, t)

    return {
        "Cmax": cmax,
        "Tmax": tmax,
        "AUC": auc
    }

def time_to_therapeutic(t: np.ndarray, cp: np.ndarray, thresh_low: float) -> float:
    """
    finds the first time point where concentration meets the therapeutic threshold.
    """
    above_thresh = np.where(cp >= thresh_low)[0]
    if len(above_thresh) > 0:
        return t[above_thresh[0]]
    return None